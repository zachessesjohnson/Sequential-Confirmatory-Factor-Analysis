#### STRICT Bayes CFA (per-year, second-order: 44 -> 8 -> overall) + small-sample stress test ####
library(readr); library(dplyr); library(janitor); library(blavaan)


# Parallel chains (if rstan installed) and compiled model cache
if (requireNamespace("rstan", quietly = TRUE)) {
  options(mc.cores = parallel::detectCores(logical = TRUE))
  rstan::rstan_options(auto_write = TRUE)
}

set.seed(123)

# Blunt-force: Remove any masked 'character' in globalenv
if ("character" %in% ls(envir = .GlobalEnv)) {
  cat("Removing masked 'character' from globalenv()\n")
  rm(character, envir = .GlobalEnv)
}

setwd('C:/Users/zejpo/OneDrive/Documents/Work/SCFA Paper')
df <- readr::read_csv("wjp_data_csv_edited.csv")

# Keep composite columns as-is; clean non-composite names like your ML script

composite_cols <- grep("^Factor \\d+", names(df), value = TRUE)
df <- df %>%
  dplyr::rename_with(~ janitor::make_clean_names(.x), .cols = setdiff(names(df), composite_cols))

stopifnot(all(c("country","year") %in% names(df)))

# Identify the 44 item names by pillar prefix (no dropping)
cleaned <- names(df)
factor_vars <- list(
  constraints_on_government_powers = cleaned[grepl("^x1_", cleaned)],
  absence_of_corruption            = cleaned[grepl("^x2_", cleaned)],
  open_government                  = cleaned[grepl("^x3_", cleaned)],
  fundamental_rights               = cleaned[grepl("^x4_", cleaned)],
  order_and_security               = cleaned[grepl("^x5_", cleaned)],
  regulatory_enforcement           = cleaned[grepl("^x6_", cleaned)],
  civil_justice                    = cleaned[grepl("^x7_", cleaned)],
  criminal_justice                 = cleaned[grepl("^x8_", cleaned)]
)
items_all <- unique(unlist(factor_vars, use.names = FALSE))

# Build the exact hierarchical lavaan model string (items -> 8 pillars -> overall)
build_model <- function(factor_vars) {
  first_order <- vapply(names(factor_vars), function(fct) {
    inds <- factor_vars[[fct]]
    paste0(fct, " =~ ", paste(inds, collapse = " + "))
  }, base::character(1))
  second_order <- paste0("overall =~ ", paste(names(factor_vars), collapse = " + "))
  paste(c(first_order, second_order), collapse = "\n")
}

# --- Version-robust bcfa wrapper: filter unsupported args + target fallback ---
bcfa_safe <- function(args) {
  # Aggressively detach all non-essential packages except blavaan, dplyr, janitor, readr, base, stats, utils, methods
  base_pkgs <- c("package:blavaan","package:dplyr","package:janitor","package:readr","package:base","package:stats","package:utils","package:methods")
  attached <- search()
  to_detach <- setdiff(grep("^package:", attached, value=TRUE), base_pkgs)
  if (length(to_detach) > 0) {
    cat("\n--- DEBUG: Detaching non-essential packages ---\n"); print(to_detach)
    for (pkg in to_detach) try(detach(pkg, character.only=TRUE, unload=TRUE), silent=TRUE)
  }
  # Check for closures in all arguments
  cat("\n--- DEBUG: Argument types and values before bcfa call ---\n")
  for (nm in names(args)) {
    cat(sprintf("  %s: %s\n", nm, paste(class(args[[nm]]), collapse=",")))
    if (is.function(args[[nm]])) {
      cat(sprintf("\n*** ERROR: Argument '%s' is a function/closure!\n", nm))
      print(args[[nm]])
      stop(sprintf("Argument '%s' is a function/closure, not allowed!", nm))
    }
    if (is.list(args[[nm]]) && any(vapply(args[[nm]], is.function, logical(1)))) {
      cat(sprintf("\n*** ERROR: Argument '%s' contains a function/closure!\n", nm))
      print(args[[nm]])
      stop(sprintf("Argument '%s' contains a function/closure, not allowed!", nm))
    }
  }
  # Check model type explicitly
  if (!is.character(args[["model"]])) stop("model argument is not character!")
  bcfa_formals <- names(formals(blavaan::bcfa))
  cat("\n--- DEBUG: blavaan::bcfa formals ---\n"); print(bcfa_formals)
  # Remove any object named 'sample' from global environment to avoid masking
  if ("sample" %in% ls(envir = .GlobalEnv)) {
    rm(sample, envir = .GlobalEnv)
    cat("\n--- DEBUG: Removed 'sample' from globalenv() to avoid masking.\n")
  }
  # Rename data_ to data for the call, then check type
  if ("data_" %in% names(args)) {
    # Coerce to data.frame to ensure compatibility with blavaan::bcfa
    args[["data"]] <- as.data.frame(args[["data_"]])
    args[["data_"]] <- NULL
    # Now check data type after renaming
    if (!inherits(args[["data"]], c("data.frame", "tbl", "tbl_df"))) stop("data argument is not a data.frame!")
  }
  # Rename n_sample to sample for the call
  if ("n_sample" %in% names(args)) {
    args[["sample"]] <- args[["n_sample"]]
    args[["n_sample"]] <- NULL
  }
  # Extract model and data as positional arguments, remove model/data/std.lv from named args
  model_arg <- args[["model"]]; args[["model"]] <- NULL
  data_arg <- args[["data"]]; args[["data"]] <- NULL
  if ("std.lv" %in% names(args)) args[["std.lv"]] <- NULL
  # Check for unsupported arguments (now only named args)
  unsupported <- setdiff(names(args), bcfa_formals)
  if (length(unsupported) > 0) {
    stop(paste0("Unsupported argument(s) for blavaan::bcfa: ", paste(unsupported, collapse=", ")))
  }
  # Print the value and class of sample argument if present
  if ("sample" %in% names(args)) {
    cat("\n--- DEBUG: sample argument value ---\n"); print(args$sample); cat("class:", class(args$sample), "\n")
  }
  # Compose argument list: positional (model, data), then named (EXPLICITLY named)
  arglist <- c(list(model = model_arg, data = data_arg), args, list(target = "stanclassic"))
  cat("\n--- DEBUG: arglist to do.call ---\n")
  str(arglist)
  for (i in seq_along(arglist)) {
    if (is.function(arglist[[i]])) {
      cat(sprintf("*** ERROR: arglist[[%d]] is a function! Name: %s\n", i, names(arglist)[i]))
      print(arglist[[i]])
    }
  }
  try1 <- try(do.call(blavaan::bcfa, arglist), silent = TRUE)
  if (!inherits(try1, "try-error")) return(try1)
  # Fallback to default target
  arglist2 <- c(list(model = model_arg, data = data_arg), args)
  cat("\n--- DEBUG: arglist2 to do.call ---\n")
  str(arglist2)
  for (i in seq_along(arglist2)) {
    if (is.function(arglist2[[i]])) {
      cat(sprintf("*** ERROR: arglist2[[%d]] is a function! Name: %s\n", i, names(arglist2)[i]))
      print(arglist2[[i]])
    }
  }
  try2 <- try(do.call(blavaan::bcfa, arglist2), silent = TRUE)
  if (!inherits(try2, "try-error")) return(try2)
  stop(attr(try2, "condition")$message)
}

# Per-year strict BCFA with optional sub-sampling of countries (to test small-N)
fit_bcfa_year <- function(dat_year, n_countries = Inf, chains=3, burnin=750, n_sample=750) {
  # listwise complete cases on the 44 items to keep BCFA call simple & fair
  dat_cc <- dat_year %>% dplyr::filter(complete.cases(dplyr::across(all_of(items_all))))
  # sub-sample countries (fixed-size per year) to emulate small N
  if (is.finite(n_countries)) {
    ctry <- unique(dat_cc$country)
    if (length(ctry) == 0) return(list(success=FALSE, error="no complete-case data"))
    take <- min(n_countries, length(ctry))
    keep_c <- sample(ctry, size = take, replace = FALSE)
    dat_cc <- dat_cc %>% dplyr::filter(country %in% keep_c)
  }
  # Keep only the 44 item columns for CFA (avoid closure errors from other columns)
  dat_cc <- dat_cc %>% dplyr::select(dplyr::all_of(items_all))
  if (nrow(dat_cc) < 10) return(list(success=FALSE, error="too few rows after complete-case filter"))
  
  model_str <- build_model(factor_vars)
  t0 <- Sys.time()
  # Debugging output
  cat("\n--- DEBUG: Model string ---\n"); print(model_str)
  cat("\n--- DEBUG: Model string class ---\n"); print(class(model_str))
  cat("\n--- DEBUG: Data columns ---\n"); print(names(dat_cc))
  cat("\n--- DEBUG: Data class ---\n"); print(class(dat_cc))
  cat("\n--- DEBUG: First few rows ---\n"); print(utils::head(dat_cc))
  # Print global environment objects to check for masking
  cat("\n--- DEBUG: ls() in globalenv() ---\n"); print(ls(envir = .GlobalEnv))
  cat("\n--- DEBUG: objects() in globalenv() ---\n"); print(objects(envir = .GlobalEnv))
  # Print argument types for bcfa_safe
  bcfa_args <- list(
    model    = model_str,
    data     = dat_cc,
    n.chains = chains,
    burnin   = burnin,
    sample   = n_sample
    # No std.lv, no dp, no dpriors
  )
  cat("\n--- DEBUG: bcfa_safe argument types ---\n")
  for (nm in names(bcfa_args)) {
    cat(sprintf("  %s: %s\n", nm, paste(class(bcfa_args[[nm]]), collapse=",")))
  }
  fit <- tryCatch({
    bcfa_safe(bcfa_args)
  }, error = function(e) e)
  t1 <- Sys.time()
  
  # If it blew up before sampling
  if (inherits(fit, "error")) {
    return(list(success=FALSE, error=fit$message, n=nrow(dat_cc),
                elapsed=as.numeric(difftime(t1,t0,units="secs"))))
  }
  
  # If Stan reports no samples or pathologies, record that
  sfit <- try(blavaan::blavInspect(fit, "stanfit"), silent = TRUE)
  rhat_max <- NA_real_; divs <- NA_integer_
  if (!inherits(sfit, "try-error") && requireNamespace("rstan", quietly = TRUE)) {
    smry <- try(rstan::summary(sfit)$summary, silent = TRUE)
    if (!inherits(smry, "try-error")) {
      rhat_max <- suppressWarnings(max(smry[,"Rhat"], na.rm=TRUE))
    }
    sp <- try(rstan::get_sampler_params(sfit, inc_warmup=FALSE), silent=TRUE)
    if (!inherits(sp, "try-error")) {
      divs <- sum(vapply(sp, function(m) sum(m[,"divergent__"]), integer(1)))
    }
  }
  
  ppp  <- suppressWarnings(try(fitMeasures(fit, "ppp"),  silent = TRUE))
  waic <- suppressWarnings(try(fitMeasures(fit, "waic"), silent = TRUE))
  
  list(
    success = TRUE,
    n = nrow(dat_cc),
    ppp = if (!inherits(ppp, "try-error")) as.numeric(ppp) else NA_real_,
    waic = if (!inherits(waic,"try-error")) as.numeric(waic) else NA_real_,
    rhat_max = rhat_max,
    divergences = divs,
    elapsed = as.numeric(difftime(t1,t0,units="secs")),
    fit = fit # keep it so you can inspect/print later if it did run
  )
}


# ---------- Intensive environment cleanup and robust stress test ----------
# Remove any objects from globalenv that could mask bcfa arguments (but NOT required functions)
cleanup_names <- c("sample", "data", "model", "fit", "bcfa_args")
for (nm in cleanup_names) {
  if (nm %in% ls(envir = .GlobalEnv)) {
    rm(list=nm, envir = .GlobalEnv)
    cat(sprintf("\n--- DEBUG: Removed '%s' from globalenv() to avoid masking.\n", nm))
  }
}

years <- sort(unique(df$year))
n_grid <- c(10, 15, 20, 30, 40, 60, 80, Inf)  # Inf = use all available complete-case countries

results <- list()
for (yr in years) {
  dat_y <- df %>% filter(year == yr)
  for (n in n_grid) {
    cat(sprintf("\n[BCFA] Year %s | N target = %s\n", yr, ifelse(is.finite(n), n, "ALL")))
    # Run fit_bcfa_year in a local environment to avoid global masking
    out <- local({
      # Defensive: check for closures in arguments before calling bcfa_safe
      res <- fit_bcfa_year(dat_y, n_countries = n, chains=3, burnin=750, n_sample=750)
      # If bcfa_safe failed, print all argument types
      if (!isTRUE(res$success)) {
        cat("\n--- DEBUG: fit_bcfa_year returned failure.\n")
      }
      res
    })
    results[[paste0(yr,"__",n)]] <- c(list(year=as.character(yr), N_target=ifelse(is.finite(n), n, NA)), out[setdiff(names(out), "fit")])
    if (!out$success) cat("  -> FAILED:", out$error, "\n") else {
      cat("  -> OK | n=", out$n,
          "| PPP=", round(out$ppp,3),
          "| rhat_max=", ifelse(is.na(out$rhat_max),"NA",round(out$rhat_max,3)),
          "| divs=", ifelse(is.na(out$divergences),"NA",out$divergences),
          "| time=", round(out$elapsed,1), "s\n")
    }
  }
}

bcfa_summary <- dplyr::bind_rows(lapply(results, as.data.frame))
print(bcfa_summary)

# Optional: record package versions for reproducibility
cat("Versions — blavaan:", as.character(packageVersion("blavaan")),
    "| lavaan:", as.character(packageVersion("lavaan")),
    "| rstan:",  if (requireNamespace("rstan", quietly=TRUE)) as.character(packageVersion("rstan")) else "not installed",
    "\n")

#### Export + visualize BCFA stress-test results ################################
# Columns used: year, N_target, success, n, ppp, rhat_max, divergences, elapsed, error
library(dplyr); library(stringr); library(knitr); library(kableExtra); library(ggplot2)

# ---- Fix/normalize bcfa_summary columns so downstream code won't crash ----
required_cols <- c("year","N_target","success","n","ppp","rhat_max","divergences","elapsed","error")
missing_cols  <- setdiff(required_cols, names(bcfa_summary))
for (nm in missing_cols) {
  bcfa_summary[[nm]] <- if (nm == "error") NA_character_ else NA_real_
}
bcfa_summary <- bcfa_summary %>%
  dplyr::mutate(
    year        = as.character(.data$year),
    N_target    = suppressWarnings(as.numeric(.data$N_target)),
    success     = as.logical(.data$success),
    n           = suppressWarnings(as.integer(.data$n)),
    ppp         = suppressWarnings(as.numeric(.data$ppp)),
    rhat_max    = suppressWarnings(as.numeric(.data$rhat_max)),
    divergences = suppressWarnings(as.integer(.data$divergences)),
    elapsed     = suppressWarnings(as.numeric(.data$elapsed)),
    error       = as.character(.data$error)
  )

# 0) Clean and augment table data ----------------------------------------------
bcfa_tbl <- bcfa_summary %>%
  dplyr::mutate(
    N_label = ifelse(is.na(N_target), "ALL", as.character(N_target)),
    Status  = ifelse(success, "OK", "FAIL"),
    ppp     = ifelse(is.na(ppp), NA, round(ppp, 3)),
    rhat_max = ifelse(is.na(rhat_max), NA, round(rhat_max, 3)),
    divergences = ifelse(is.na(divergences), NA, divergences),
    elapsed = round(elapsed, 1),
    error_short = ifelse(success, "", stringr::str_trunc(as.character(error), 90))
  ) %>%
  dplyr::select(year, N_label, n, Status, ppp, rhat_max, divergences, elapsed, error_short) %>%
  dplyr::arrange(year, factor(N_label, levels=c("10","15","20","30","40","60","80","ALL")))

# 1) LaTeX table with red FAIL / green OK --------------------------------------
latex_tab <- bcfa_tbl %>%
  kableExtra::kbl(format = "latex", booktabs = TRUE, longtable = TRUE,
                  caption = "Bayesian CFA (traditional second-order) success by year and sample size.",
                  col.names = c("Year","N target","Obs used","Status","PPP","Max R-hat","Divergences","Time (s)","Error (truncated)")) %>%
  kableExtra::kable_styling(latex_options = c("hold_position","striped")) %>%
  kableExtra::row_spec(which(bcfa_tbl$Status == "FAIL"), background = "#ffe6e6") %>%
  kableExtra::row_spec(which(bcfa_tbl$Status == "OK"),   background = "#e8f5e9")

cat(latex_tab, file = "bcfa_summary_table.tex")

# 2) CSV export for replication -------------------------------------------------
write.csv(bcfa_tbl, "bcfa_summary_table.csv", row.names = FALSE)

# 3) Success-rate figure vs N (per year) ---------------------------------------
rate_df <- bcfa_summary %>%
  dplyr::mutate(N_label = ifelse(is.na(N_target), "ALL", as.character(N_target))) %>%
  dplyr::group_by(year, N_label) %>%
  dplyr::summarise(pass_rate = mean(success, na.rm = TRUE), .groups = "drop")

p <- ggplot2::ggplot(rate_df,
                     ggplot2::aes(x = factor(N_label, levels=c("10","15","20","30","40","60","80","ALL")),
                                  y = pass_rate, group = year)) +
  ggplot2::geom_line(alpha = 0.7) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(x = "Target countries per year (ALL = full sample)",
                y = "BCFA success rate",
                title = "Traditional Bayes CFA success vs. per-year sample size",
                subtitle = "44 items → 8 subfactors → 1 second-order factor; HMC via Stan") +
  ggplot2::ylim(0,1) +
  ggplot2::theme_minimal(base_size = 12)
ggplot2::ggsave("bcfa_success_rate_by_N.png", p, width = 8, height = 5, dpi = 300)

# 4) One-line summary for the text ---------------------------------------------
overall <- rate_df %>% dplyr::group_by(N_label) %>% dplyr::summarise(pass_rate = mean(pass_rate), .groups="drop")
print(overall)

