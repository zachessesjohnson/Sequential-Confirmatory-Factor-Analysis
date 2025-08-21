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

# Loop through each year and factor, run CFA
for (yr in years) {

  cat("\n\nYear:", yr, "\n")
  data_year <- wjp_data %>% filter(year == yr)

  # Build hierarchical CFA model string
  # First-order factors
  first_order <- sapply(names(factor_vars), function(fct) {
    inds <- factor_vars[[fct]]
    if (length(inds) > 1) paste0(fct, " =~ ", paste(inds, collapse = " + ")) else NULL
  })
  # Second-order factor
  second_order <- paste0("overall =~ ", paste(names(factor_vars), collapse = " + "))
  model_str <- paste(c(first_order, second_order), collapse = "\n")

  cat("\nHierarchical CFA model for year:", yr, "\n")
  fit <- tryCatch({
    cfa(model_str, data = data_year, std.lv = TRUE, missing = 'fiml')
  }, error = function(e) {
    cat("Error for hierarchical CFA in year", yr, ":", e$message, "\n")
    NULL
  })
  if (!is.null(fit)) {
    print(summary(fit, fit.measures = TRUE, standardized = TRUE))
    # Extract overall factor scores
    fscores <- tryCatch({
      lavPredict(fit)[, "overall"]
    }, error = function(e) NULL)
    if (!is.null(fscores)) {
      assign(paste0("factor_scores_", yr), fscores, envir = .GlobalEnv)
      cat("\nSummary of overall factor scores for year", yr, ":\n")
      print(summary(as.numeric(fscores)))
    }
    # Calculate omega for overall factor (using standardized loadings from second-order)
    std_loadings <- tryCatch({
      inspect(fit, 'std')$lambda[, "overall"]
    }, error = function(e) NULL)
    error_vars <- tryCatch({
      diag(inspect(fit, 'std')$theta)
    }, error = function(e) NULL)
    if (!is.null(std_loadings) && !is.null(error_vars)) {
      omega <- sum(std_loadings^2, na.rm = TRUE) / (sum(std_loadings^2, na.rm = TRUE) + sum(error_vars, na.rm = TRUE))
      assign(paste0("omega_", yr), omega, envir = .GlobalEnv)
      cat("Reliability (omega) for overall factor in year", yr, ":", round(omega, 3), "\n")
    }
  } else {
    cat("Hierarchical CFA failed for year", yr, "\n")
  }
}

# Print all factor scores and omegas at the end
cat("\n\nSummary of all overall factor scores and omegas by year:\n")
for (yr in years) {
  fs_name <- paste0("factor_scores_", yr)
  om_name <- paste0("omega_", yr)
  if (exists(fs_name, envir = .GlobalEnv)) {
    cat("\nYear:", yr, "\n")
    print(summary(get(fs_name, envir = .GlobalEnv)))
    if (exists(om_name, envir = .GlobalEnv)) {
      cat("Omega:", round(get(om_name, envir = .GlobalEnv), 3), "\n")
    }
  }
}