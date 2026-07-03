# =============================================================================
#  tests/test_scfa_propagation.R
#
#  Self-contained test suite for R/scfa_propagation.R.
#  Run from the repository root with:
#    Rscript tests/test_scfa_propagation.R
#
#  Requires: lavaan (>= 0.6)
# =============================================================================

library(lavaan)

# Locate R/scfa_propagation.R relative to this test file or the working directory
.this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
if (!is.null(.this_file)) {
  source(file.path(dirname(dirname(.this_file)), "R", "scfa_propagation.R"))
} else {
  source("R/scfa_propagation.R")
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
pass <- function(msg) cat("PASS:", msg, "\n")
fail <- function(msg) stop("FAIL: ", msg, call. = FALSE)

expect_true <- function(cond, msg) if (!isTRUE(cond)) fail(msg) else pass(msg)
expect_error <- function(expr, msg) {
  ok <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  if (!ok) fail(msg) else pass(msg)
}

# ---------------------------------------------------------------------------
# Shared fixture: two-factor simple-structure model
#   f1: loadings ~0.7, residuals ~0.51  -> lower reliability
#   f2: loadings ~0.9, residuals ~0.19  -> higher reliability
# ---------------------------------------------------------------------------
set.seed(42)

pop_model <- "
  f1 =~ 0.7*x1 + 0.7*x2 + 0.7*x3
  f2 =~ 0.9*x4 + 0.9*x5 + 0.9*x6
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f1 ~~ 0*f2
  x1 ~~ 0.51*x1
  x2 ~~ 0.51*x2
  x3 ~~ 0.51*x3
  x4 ~~ 0.19*x4
  x5 ~~ 0.19*x5
  x6 ~~ 0.19*x6
"
dat <- simulateData(pop_model, sample.nobs = 500, model.type = "sem")

fit1 <- cfa("f1 =~ x1 + x2 + x3
             f2 =~ x4 + x5 + x6", data = dat)

# Extract internals for analytic checks
lambda_hat <- lavInspect(fit1, "est")$lambda
theta_hat  <- diag(lavInspect(fit1, "est")$theta)
phi_hat    <- diag(lavInspect(fit1, "est")$psi)

# ---------------------------------------------------------------------------
cat("\n--- scfa_factor_information ---\n")
# ---------------------------------------------------------------------------
info <- scfa_factor_information(fit1)

expect_true(is.numeric(info) && length(info) == 2,
            "returns numeric vector of length 2")

expect_true(all(info > 0),
            "all information values are positive")

# Analytic check: I_k = sum_j lambda_jk^2 / theta_j
info_manual <- colSums(lambda_hat^2 / theta_hat)
expect_true(max(abs(info - info_manual)) < 1e-12,
            "matches manual colSums(lambda^2 / theta) formula")

expect_true(info["f2"] > info["f1"],
            "f2 (stronger loadings, lower residuals) has higher information than f1")

# ---------------------------------------------------------------------------
cat("\n--- scfa_propagation_variance ---\n")
# ---------------------------------------------------------------------------
psi <- scfa_propagation_variance(fit1)

expect_true(all(psi > 0),
            "all propagation variances are positive")

expect_true(max(abs(psi - 1 / info)) < 1e-14,
            "psi_nu_k = 1 / I_k identity holds exactly")

expect_true(psi["f1"] > psi["f2"],
            "f1 (less informative) has higher propagation variance")

# ---------------------------------------------------------------------------
cat("\n--- scfa_factor_reliability ---\n")
# ---------------------------------------------------------------------------
rho <- scfa_factor_reliability(fit1)

expect_true(all(rho > 0) && all(rho <= 1),
            "all reliabilities are in (0, 1]")

# Analytic check using actual estimated phi (lavaan uses marker-variable scaling)
rho_manual <- phi_hat * info / (phi_hat * info + 1)
expect_true(max(abs(rho - rho_manual)) < 1e-10,
            "matches analytic formula phi*I / (phi*I + 1)")

expect_true(rho["f2"] > rho["f1"],
            "f2 is more reliable than f1")

# ---------------------------------------------------------------------------
cat("\n--- scfa_propagation_diagnostics ---\n")
# ---------------------------------------------------------------------------
diag_tbl <- scfa_propagation_diagnostics(fit1, threshold = 0.90)

expect_true(is.data.frame(diag_tbl),
            "returns a data.frame")

expect_true(nrow(diag_tbl) == 2,
            "one row per factor (2 factors)")

expect_true(
  all(c("factor", "n_indicators", "I_k", "psi_nu", "phi_k", "rho_k", "flag")
      %in% names(diag_tbl)),
  "all expected columns present"
)

expect_true(all(diag_tbl$n_indicators == 3),
            "3 indicators per factor correctly counted")

expect_true(all(diag_tbl$I_k   == info[diag_tbl$factor]),  "I_k column matches scfa_factor_information()")
expect_true(all(diag_tbl$psi_nu == psi[diag_tbl$factor]),  "psi_nu column matches scfa_propagation_variance()")
expect_true(all(diag_tbl$rho_k  == rho[diag_tbl$factor]),  "rho_k column matches scfa_factor_reliability()")

# With threshold=0.90, f1 (rho~0.77) should be flagged; f2 (rho~0.92) should not
expect_true( diag_tbl$flag[diag_tbl$factor == "f1"],
             "f1 flagged when threshold = 0.90 and rho_f1 < 0.90")
expect_true(!diag_tbl$flag[diag_tbl$factor == "f2"],
             "f2 not flagged when rho_f2 > 0.90")

# Default threshold = 0.70: neither factor should be flagged
diag_default <- scfa_propagation_diagnostics(fit1)
expect_true(!any(diag_default$flag),
            "no flags with default threshold = 0.70 when both rho > 0.70")

# Invalid threshold rejected
expect_error(scfa_propagation_diagnostics(fit1, threshold = 1.5),
             "threshold > 1 correctly rejected")
expect_error(scfa_propagation_diagnostics(fit1, threshold = 0),
             "threshold = 0 correctly rejected")
expect_error(scfa_propagation_diagnostics(fit1, threshold = -0.1),
             "negative threshold correctly rejected")

# ---------------------------------------------------------------------------
cat("\n--- scfa_correct_loadings ---\n")
# ---------------------------------------------------------------------------
scores_reg   <- as.data.frame(lavPredict(fit1, method = "regression"))
fit2         <- cfa("g =~ f1 + f2", data = scores_reg)
corrected    <- suppressMessages(scfa_correct_loadings(fit2, fit1))
lambda2_raw  <- lavInspect(fit2, "est")$lambda

expect_true(is.matrix(corrected) && all(dim(corrected) == dim(lambda2_raw)),
            "corrected matrix has same dimensions as raw loading matrix")

nz <- lambda2_raw != 0
expect_true(all(corrected[nz] >= lambda2_raw[nz]),
            "corrected loadings >= raw loadings (attenuation reversed)")

# Zeros stay zero (cross-loadings fixed to zero are untouched)
expect_true(all(corrected[!nz] == 0),
            "zero entries remain zero after correction")

# Analytic check: each row k corrected = raw[k,] / rho_k
rho_vals <- rho[rownames(lambda2_raw)]
for (k in seq_len(nrow(lambda2_raw))) {
  cols_nz <- lambda2_raw[k, ] != 0
  expected <- lambda2_raw[k, cols_nz] / rho_vals[k]
  expect_true(max(abs(corrected[k, cols_nz] - expected)) < 1e-12,
              paste0("row ", k, " (", rownames(lambda2_raw)[k],
                     ") divided by rho_", k, " correctly"))
}

# Non-lavaan input rejected
expect_error(scfa_factor_information(list()),
             "non-lavaan object rejected by scfa_factor_information")
expect_error(scfa_propagation_variance(data.frame()),
             "non-lavaan object rejected by scfa_propagation_variance")
expect_error(scfa_factor_reliability("not lavaan"),
             "non-lavaan object rejected by scfa_factor_reliability")

# Mismatched row count rejected (Stage-2 model has wrong number of indicators)
fit_one <- cfa("f1 =~ x1 + x2 + x3", data = dat)
scores_one <- as.data.frame(lavPredict(fit_one, method = "regression"))
fit2_one   <- cfa("g =~ f1", data = scores_one)
expect_error(suppressMessages(scfa_correct_loadings(fit2_one, fit1)),
             "mismatched indicator count rejected by scfa_correct_loadings")

# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("All tests passed.\n")
cat("========================================\n")
