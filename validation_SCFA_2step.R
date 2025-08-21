####### Validation of Sequential CFA on WJP Rule of Law Index ###

# Load required libraries
library(readr)
library(dplyr)
library(janitor)
library(lavaan)

setwd('C:/Users/zejpo/OneDrive/Documents/Work/SCFA Paper')

# Read data
wjp_data <- read_csv("wjp_data_csv_edited.csv")


# Clean variable names to snake_case, but keep composite factor columns as-is
composite_cols <- grep("^Factor \\d+", names(wjp_data), value = TRUE)
other_cols <- setdiff(names(wjp_data), composite_cols)
wjp_data <- wjp_data %>%
  rename_with(~ janitor::make_clean_names(.x), .cols = other_cols)

# Get the cleaned names for all columns
cleaned_names <- names(wjp_data)

# Identify factor indicators (excluding composite columns) using cleaned names
factor_vars <- list(
  constraints_on_government_powers = cleaned_names[grepl("^x1_", cleaned_names)],
  absence_of_corruption = cleaned_names[grepl("^x2_", cleaned_names)],
  open_government = cleaned_names[grepl("^x3_", cleaned_names)],
  fundamental_rights = cleaned_names[grepl("^x4_", cleaned_names)],
  order_and_security = cleaned_names[grepl("^x5_", cleaned_names)],
  regulatory_enforcement = cleaned_names[grepl("^x6_", cleaned_names)],
  civil_justice = cleaned_names[grepl("^x7_", cleaned_names)],
  criminal_justice = cleaned_names[grepl("^x8_", cleaned_names)]
)

# Remove any composite columns from factor_vars (shouldn't match, but for safety)
factor_vars <- lapply(factor_vars, function(vars) vars[!vars %in% janitor::make_clean_names(composite_cols)])

# Get unique years
years <- unique(wjp_data$year)


# Step 1: Run CFA for each factor and extract factor scores for each year
all_factor_scores <- list()
first_stage_fits <- list()
# Store first-stage omegas for summary
first_stage_omegas <- list()
for (yr in years) {
  cat("\n\nYear:", yr, "\n")
  data_year <- wjp_data %>% filter(year == yr)
  year_scores_list <- list()
  n_rows <- nrow(data_year)
  fitmat <- matrix(NA, nrow = length(factor_vars), ncol = 4)
  rownames(fitmat) <- names(factor_vars)
  colnames(fitmat) <- c("CFI", "TLI", "RMSEA", "SRMR")
  omegas <- rep(NA, length(factor_vars))
  names(omegas) <- names(factor_vars)
  for (i in seq_along(factor_vars)) {
    fct <- names(factor_vars)[i]
    indicators <- factor_vars[[fct]]
    if (length(indicators) > 1) {
      model_str <- paste0(fct, " =~ ", paste(indicators, collapse = " + "))
      fit <- tryCatch({
        cfa(model_str, data = data_year, std.lv = TRUE, missing = 'fiml')
      }, error = function(e) NULL)
      # Print fit indices and standardized loadings for first-stage CFA
      if (!is.null(fit)) {
        fitmeas <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
        fitmat[i,] <- as.numeric(fitmeas)
        std_loadings <- tryCatch({
          inspect(fit, 'std')$lambda[,1]
        }, error = function(e) NULL)
        error_vars <- tryCatch({
          diag(inspect(fit, 'std')$theta)
        }, error = function(e) NULL)
        # Calculate omega for this factor
        if (!is.null(std_loadings) && !is.null(error_vars)) {
          omega <- sum(std_loadings^2, na.rm = TRUE) / (sum(std_loadings^2, na.rm = TRUE) + sum(error_vars, na.rm = TRUE))
          omegas[i] <- omega
        }
        cat("\nFirst-stage CFA for ", fct, " in year ", yr, ":\n", sep = "")
        print(fitmeas)
        cat("Standardized loadings:\n")
        print(std_loadings)
        if (!is.na(omegas[i])) cat("Omega:", round(omegas[i], 3), "\n")
        fscores <- tryCatch({
          lavPredict(fit)
        }, error = function(e) NULL)
        if (!is.null(fscores)) {
          year_scores_list[[fct]] <- as.numeric(fscores)
        } else {
          year_scores_list[[fct]] <- rep(NA, n_rows)
        }
      } else {
        cat("First-stage CFA failed for ", fct, " in year ", yr, "\n", sep = "")
        year_scores_list[[fct]] <- rep(NA, n_rows)
      }
    } else {
      year_scores_list[[fct]] <- rep(NA, n_rows)
    }
  }
  year_scores <- as.data.frame(year_scores_list)
  all_factor_scores[[as.character(yr)]] <- year_scores
  first_stage_fits[[as.character(yr)]] <- fitmat
  first_stage_omegas[[as.character(yr)]] <- omegas
}

# Step 2: For each year, run CFA on the 8 factor scores
for (yr in years) {
  cat("\n\nSecond-step CFA for year:", yr, "\n")
  year_scores <- all_factor_scores[[as.character(yr)]]
  # Remove columns with all NA (in case any factor failed)
  year_scores <- year_scores[, colSums(is.na(year_scores)) < nrow(year_scores), drop = FALSE]
  # Build one-factor CFA model string
  if (ncol(year_scores) > 1) {
    model_str <- paste0("overall =~ ", paste(colnames(year_scores), collapse = " + "))
    fit2 <- tryCatch({
      cfa(model_str, data = year_scores, std.lv = TRUE, missing = 'fiml')
    }, error = function(e) NULL)
    if (!is.null(fit2)) {
      print(summary(fit2, fit.measures = TRUE, standardized = TRUE))
      # Extract overall factor scores
      fscores2 <- tryCatch({
        lavPredict(fit2)
      }, error = function(e) NULL)
      if (!is.null(fscores2)) {
        assign(paste0("overall_factor_scores_", yr), fscores2, envir = .GlobalEnv)
        cat("\nSummary of overall factor scores for year", yr, ":\n")
        print(summary(as.numeric(fscores2)))
      }
      # Calculate omega for overall factor
      std_loadings <- tryCatch({
        inspect(fit2, 'std')$lambda[,1]
      }, error = function(e) NULL)
      error_vars <- tryCatch({
        diag(inspect(fit2, 'std')$theta)
      }, error = function(e) NULL)
      if (!is.null(std_loadings) && !is.null(error_vars)) {
        omega <- sum(std_loadings^2, na.rm = TRUE) / (sum(std_loadings^2, na.rm = TRUE) + sum(error_vars, na.rm = TRUE))
        assign(paste0("overall_omega_", yr), omega, envir = .GlobalEnv)
        cat("Reliability (omega) for overall factor in year", yr, ":", round(omega, 3), "\n")
      }
    } else {
      cat("Second-step CFA failed for year", yr, "\n")
    }
  } else {
    cat("Not enough factors for second-step CFA in year", yr, "\n")
  }
}

# Print all overall factor scores and omegas at the end
cat("\n\nSummary of all overall factor scores and omegas by year:\n")
for (yr in years) {
  fs_name <- paste0("overall_factor_scores_", yr)
  om_name <- paste0("overall_omega_", yr)
  if (exists(fs_name, envir = .GlobalEnv)) {
    cat("\nYear:", yr, "\n")
    print(summary(get(fs_name, envir = .GlobalEnv)))
    if (exists(om_name, envir = .GlobalEnv)) {
      cat("Omega:", round(get(om_name, envir = .GlobalEnv), 3), "\n")
    }
  }
}

# Print summary of first-stage CFA fit indices by year
cat("\n\nSummary of first-stage CFA fit indices by year (mean across factors):\n")
for (yr in years) {
  fitmat <- first_stage_fits[[as.character(yr)]]
  omegas <- first_stage_omegas[[as.character(yr)]]
  # Remove rows with all NA (factors not estimated)
  fitmat <- fitmat[rowSums(is.na(fitmat)) < ncol(fitmat), , drop = FALSE]
  omegas <- omegas[!is.na(omegas)]
  if (nrow(fitmat) > 0) {
    means <- colMeans(fitmat, na.rm = TRUE)
    avg_omega <- if (length(omegas) > 0) round(mean(omegas, na.rm = TRUE), 3) else NA
    cat("\nYear:", yr, "\n")
    print(round(means, 3))
    cat("Average omega:", avg_omega, "\n")
  }
}


  
# First-stage average omega increases over time (from ~0.66 to ~0.75), indicating improving reliability of the subindices.
#CFI and SRMR are consistently good, while TLI is moderate and RMSEA is high (as discussed, likely due to few indicators per factor).
#By reporting both fit indices and average omega, you provide a transparent, comprehensive assessment of the measurement quality at the first stage.
#Interpretation:
  
#Your subindices are generally reliable (omega > 0.7 is acceptable, >0.75 is good).
#The second-stage CFA is built on a solid foundation, and your reporting allows readers to judge the quality of each step.
#Conclusion:
#You can confidently state that your overall factor is based on subindices that are, on average, acceptably reliable and reasonably well-fitting. This strengthens the validity of your sequential CFA approach and your results.


# --------- Compile First-Stage CFA Results ---------
first_stage_results <- data.frame()

for (yr in years) {
  fitmat <- first_stage_fits[[as.character(yr)]]
  omegas <- first_stage_omegas[[as.character(yr)]]
  
  if (!is.null(fitmat)) {
    for (fct in rownames(fitmat)) {
      if (!all(is.na(fitmat[fct,]))) {
        result_row <- data.frame(
          year = yr,
          factor = fct,
          CFI = fitmat[fct, "CFI"],
          TLI = fitmat[fct, "TLI"],
          RMSEA = fitmat[fct, "RMSEA"],
          SRMR = fitmat[fct, "SRMR"],
          Omega = omegas[fct]
        )
        first_stage_results <- bind_rows(first_stage_results, result_row)
      }
    }
  }
}

# --------- Compile Second-Stage CFA Results ---------
second_stage_results <- data.frame()

for (yr in years) {
  fs_name <- paste0("overall_factor_scores_", yr)
  om_name <- paste0("overall_omega_", yr)
  
  fit2_obj <- tryCatch({
    cfa(
      paste0("overall =~ ", paste(names(all_factor_scores[[as.character(yr)]]), collapse = " + ")),
      data = all_factor_scores[[as.character(yr)]],
      std.lv = TRUE, missing = "fiml"
    )
  }, error = function(e) NULL)
  
  if (!is.null(fit2_obj)) {
    fitmeas2 <- fitMeasures(fit2_obj, c("cfi", "tli", "rmsea", "srmr"))
    omega2 <- if (exists(om_name, envir = .GlobalEnv)) get(om_name, envir = .GlobalEnv) else NA
    
    second_stage_results <- bind_rows(second_stage_results, data.frame(
      year = yr,
      Overall_CFI = fitmeas2["cfi"],
      Overall_TLI = fitmeas2["tli"],
      Overall_RMSEA = fitmeas2["rmsea"],
      Overall_SRMR = fitmeas2["srmr"],
      Overall_Omega = omega2
    ))
  }
}

# --------- View or Export ---------
print("First-stage results:")
print(first_stage_results)

print("Second-stage results:")
print(second_stage_results)

# Optionally save to CSV
# write.csv(first_stage_results, "first_stage_results.csv", row.names = FALSE)
# write.csv(second_stage_results, "second_stage_results.csv", row.names = FALSE)

# Function to generate LaTeX table for all years (edit data as needed)

