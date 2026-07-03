# =============================================================================
#  tests/test_scfa_end_to_end.R
#
#  End-to-end smoke tests that mirror the README usage examples and verify
#  the full sequential CFA workflow from source() through Stage-2 correction.
#
#  Run from the repository root with:
#    Rscript tests/test_scfa_end_to_end.R
#
#  Requires: lavaan (>= 0.6)
# =============================================================================

library(lavaan)
source("R/scfa_propagation.R")

pass <- function(msg) cat("PASS:", msg, "\n")
fail <- function(msg) stop("FAIL: ", msg, call. = FALSE)

# ---------------------------------------------------------------------------
# Fixture A: three-factor simple-structure model (3 indicators each)
#   f1: loadings 0.8, residuals 0.36  (high reliability,   rho ~ 0.83)
#   f2: loadings 0.6, residuals 0.64  (moderate reliability, rho ~ 0.64)
#   f3: loadings 0.4, residuals 0.84  (low reliability,    rho ~ 0.38)
# std.lv = TRUE fixes factor variances to 1 for stable identification.
# ---------------------------------------------------------------------------
set.seed(42)

pop_a <- "
  f1 =~ 0.8*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 0.6*y4 + 0.6*y5 + 0.6*y6
  f3 =~ 0.4*y7 + 0.4*y8 + 0.4*y9
  f1 ~~ 1*f1; f2 ~~ 1*f2; f3 ~~ 1*f3
  f1 ~~ 0*f2; f1 ~~ 0*f3; f2 ~~ 0*f3
  y1 ~~ 0.36*y1; y2 ~~ 0.36*y2; y3 ~~ 0.36*y3
  y4 ~~ 0.64*y4; y5 ~~ 0.64*y5; y6 ~~ 0.64*y6
  y7 ~~ 0.84*y7; y8 ~~ 0.84*y8; y9 ~~ 0.84*y9
"
dat_a  <- simulateData(pop_a, sample.nobs = 500, model.type = "sem")
fit1_a <- cfa("f1=~y1+y2+y3\nf2=~y4+y5+y6\nf3=~y7+y8+y9",
              data = dat_a, std.lv = TRUE)

if (!fit1_a@optim$converged) fail("Fixture A Stage-1 model did not converge")
pass("Fixture A Stage-1 model converges")

# ---------------------------------------------------------------------------
# Fixture B: README two-factor workflow (subfactor1, subfactor2)
# ---------------------------------------------------------------------------
set.seed(7)
pop_b <- "
  subfactor1 =~ 0.7*item1 + 0.7*item2 + 0.7*item3
  subfactor2 =~ 0.8*item4 + 0.8*item5 + 0.8*item6
  subfactor1 ~~ 1*subfactor1; subfactor2 ~~ 1*subfactor2
  subfactor1 ~~ 0*subfactor2
  item1 ~~ 0.51*item1; item2 ~~ 0.51*item2; item3 ~~ 0.51*item3
  item4 ~~ 0.36*item4; item5 ~~ 0.36*item5; item6 ~~ 0.36*item6
"
lower_data <- simulateData(pop_b, sample.nobs = 400, model.type = "sem")
fit1_b <- cfa("subfactor1 =~ item1+item2+item3
               subfactor2 =~ item4+item5+item6",
              data = lower_data, std.lv = TRUE)

if (!fit1_b@optim$converged) fail("Fixture B Stage-1 model did not converge")
pass("Fixture B Stage-1 model converges")

# ---------------------------------------------------------------------------
cat("\n--- README: source() exposes all public functions ---\n")
# ---------------------------------------------------------------------------
env_fresh <- new.env()
tryCatch(
  source("R/scfa_propagation.R", local = env_fresh),
  error = function(e) fail(paste("source() failed:", conditionMessage(e)))
)
fns <- c("scfa_factor_information", "scfa_propagation_variance",
         "scfa_factor_reliability", "scfa_propagation_diagnostics",
         "scfa_correct_loadings")
missing_fns <- fns[!vapply(fns, exists, logical(1), envir = env_fresh, inherits = FALSE)]
if (length(missing_fns) > 0)
  fail(paste("Missing after source():", paste(missing_fns, collapse = ", ")))
pass("source('R/scfa_propagation.R') exposes all 5 public functions")

# ---------------------------------------------------------------------------
cat("\n--- Three-factor ordering and magnitudes ---\n")
# ---------------------------------------------------------------------------
info_a <- scfa_factor_information(fit1_a)
psi_a  <- scfa_propagation_variance(fit1_a)
rho_a  <- scfa_factor_reliability(fit1_a)

# DGP ordering: f1 (best) > f2 > f3 (worst)
if (!(info_a["f1"] > info_a["f2"] && info_a["f2"] > info_a["f3"]))
  fail("Information ordering f1 > f2 > f3 not preserved")
pass("Information ordering matches DGP quality (f1 > f2 > f3)")

if (!(psi_a["f1"] < psi_a["f2"] && psi_a["f2"] < psi_a["f3"]))
  fail("Propagation variance ordering f1 < f2 < f3 not preserved")
pass("Propagation variance ordering matches DGP quality (f1 < f2 < f3)")

if (!(rho_a["f1"] > rho_a["f2"] && rho_a["f2"] > rho_a["f3"]))
  fail("Reliability ordering f1 > f2 > f3 not preserved")
pass("Reliability ordering matches DGP quality (f1 > f2 > f3)")

# With phi = 1 (std.lv=TRUE), rho = I/(I+1) exactly
rho_check_a <- info_a / (info_a + 1)
if (max(abs(rho_a - rho_check_a)) > 1e-10)
  fail("rho != I/(I+1) when phi = 1")
pass("rho = I/(I+1) exactly when factor variances = 1")

# Reasonable ranges given DGP: rho in (0, 1), rho_f1 > 0.70
if (rho_a["f1"] <= 0.70) fail("f1 reliability should be > 0.70 for strong DGP")
if (rho_a["f3"] >= 0.70) fail("f3 reliability should be < 0.70 for weak DGP")
pass("f1 reliability > 0.70 and f3 reliability < 0.70 as expected from DGP")

# ---------------------------------------------------------------------------
cat("\n--- scfa_propagation_diagnostics: indicator counts and flagging ---\n")
# ---------------------------------------------------------------------------
diag_a <- scfa_propagation_diagnostics(fit1_a, threshold = 0.80)

if (!all(diag_a$n_indicators == 3))
  fail("n_indicators should be 3 for each factor")
pass("n_indicators = 3 for all 3 factors")

# f3 (rho ~0.38) must be flagged at threshold = 0.80
if (!diag_a$flag[diag_a$factor == "f3"])
  fail("f3 should be flagged at threshold = 0.80")
pass("f3 correctly flagged as weak link at threshold = 0.80")

# f1 (rho ~0.83) should NOT be flagged at default threshold = 0.70
diag_default_a <- scfa_propagation_diagnostics(fit1_a)
if (diag_default_a$flag[diag_default_a$factor == "f1"])
  fail("f1 should NOT be flagged at default threshold = 0.70")
pass("f1 not flagged at default threshold = 0.70")

# f3 always flagged at default 0.70
if (!diag_default_a$flag[diag_default_a$factor == "f3"])
  fail("f3 should be flagged at default threshold = 0.70")
pass("f3 flagged at default threshold = 0.70")

# ---------------------------------------------------------------------------
cat("\n--- Regression score shrinkage: var(reg) < var(bartlett) ---\n")
# ---------------------------------------------------------------------------
# Bartlett scores are unbiased for the factor; regression scores are
# shrunk toward zero.  Their variance ratio equals rho^2 (in expectation).
scores_bart_a <- as.data.frame(lavPredict(fit1_a, method = "bartlett"))
scores_reg_a  <- as.data.frame(lavPredict(fit1_a, method = "regression"))

var_bart <- sapply(scores_bart_a, var)
var_reg  <- sapply(scores_reg_a,  var)

if (!all(var_reg < var_bart))
  fail("Regression score variances should all be smaller than Bartlett variances")
pass("var(regression scores) < var(Bartlett scores) for all factors")

# var_reg / var_bart should be close to rho^2
var_ratio  <- var_reg / var_bart
rho_sq     <- rho_a^2
# Allow 5% relative tolerance (finite-sample approximation)
if (any(abs(var_ratio - rho_sq) / rho_sq > 0.05))
  fail(paste("var ratio deviates from rho^2 by more than 5%:",
             paste(round(abs(var_ratio - rho_sq) / rho_sq, 3), collapse = ", ")))
pass("var(reg) / var(bartlett) is within 5% of rho^2 for all factors")

# ---------------------------------------------------------------------------
cat("\n--- scfa_correct_loadings: formula-level verification ---\n")
# ---------------------------------------------------------------------------
# Create a simple regression-based Stage-2 fit with Fixture B (2 factors)
scores_reg_b  <- as.data.frame(lavPredict(fit1_b, method = "regression"))
fit2_reg_b    <- suppressWarnings(
  cfa("g =~ subfactor1 + subfactor2", data = scores_reg_b, std.lv = TRUE)
)

rho_b       <- scfa_factor_reliability(fit1_b)
lambda2_raw <- lavInspect(fit2_reg_b, "est")$lambda
corrected   <- suppressMessages(scfa_correct_loadings(fit2_reg_b, fit1_b))

# Verify formula: corrected[k,] = raw[k,] / rho_b[k]  (row-wise)
for (k in seq_len(nrow(lambda2_raw))) {
  factor_name <- rownames(lambda2_raw)[k]
  cols_nz     <- lambda2_raw[k, ] != 0
  expected    <- lambda2_raw[k, cols_nz] / rho_b[factor_name]
  got         <- corrected[k, cols_nz]
  if (max(abs(got - expected)) > 1e-12)
    fail(paste("Row", k, "(", factor_name, "): corrected != raw/rho"))
  pass(paste("Row", k, "(", factor_name, "): corrected = raw / rho correctly"))
}

# Zeros remain zero
if (!all(corrected[lambda2_raw == 0] == 0))
  fail("Zero loadings changed after correction")
pass("Zero loadings remain zero after correction")

# Corrected loadings are always >= raw in absolute value (inflation, not shrinkage)
nz <- lambda2_raw != 0
if (!all(abs(corrected[nz]) >= abs(lambda2_raw[nz]) - 1e-12))
  fail("Corrected absolute loadings should be >= raw (attenuation reversed)")
pass("Absolute corrected loadings >= absolute raw loadings")

# ---------------------------------------------------------------------------
cat("\n--- Full README two-factor workflow ---\n")
# ---------------------------------------------------------------------------
diag_b <- scfa_propagation_diagnostics(fit1_b)
if (!is.data.frame(diag_b) || nrow(diag_b) != 2)
  fail("Diagnostic table wrong shape")
pass("scfa_propagation_diagnostics() returns 2-row table for 2-factor model")

# subfactor2 (loadings 0.8) should have higher reliability than subfactor1 (0.7)
rho_b_vec <- scfa_factor_reliability(fit1_b)
if (rho_b_vec["subfactor2"] <= rho_b_vec["subfactor1"])
  fail("subfactor2 should have higher reliability than subfactor1")
pass("subfactor2 more reliable than subfactor1 in Fixture B (as expected)")

# Stage 2: use Bartlett scores and extract final index scores
scores_stage2   <- as.data.frame(lavPredict(fit1_b))   # default = Bartlett
model_stage2    <- "higher_factor =~ subfactor1 + subfactor2"
fit_stage2      <- suppressWarnings(cfa(model_stage2, data = scores_stage2, std.lv = TRUE))
index_scores    <- lavPredict(fit_stage2)
if (nrow(index_scores) != 400)
  fail("Index scores should have 400 rows")
pass("Final index scores produced for all 400 observations")

# ---------------------------------------------------------------------------
cat("\n--- Error handling ---\n")
# ---------------------------------------------------------------------------
err_msg <- tryCatch(scfa_factor_information(NULL),
                    error = function(e) conditionMessage(e))
if (!grepl("lavaan", err_msg, ignore.case = TRUE))
  fail("Error for NULL should mention 'lavaan'")
pass("Error message for non-lavaan input mentions 'lavaan'")

# scfa_correct_loadings emits a usage note
msgs <- character(0)
withCallingHandlers(
  scfa_correct_loadings(fit2_reg_b, fit1_b),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
if (length(msgs) == 0) fail("scfa_correct_loadings should emit a usage note")
if (!grepl("regression", msgs[1], ignore.case = TRUE))
  fail("Usage note should mention 'regression'")
pass("scfa_correct_loadings emits usage note mentioning 'regression'")

# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("All end-to-end tests passed.\n")
cat("========================================\n")
