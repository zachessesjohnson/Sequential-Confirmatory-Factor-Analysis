# =============================================================================
#  tests/testthat/test-scfa-propagation.R
#
#  Unit tests for the error-propagation diagnostic functions in
#  R/scfa_propagation.R.
# =============================================================================

test_that("scfa_factor_information returns correct values", {
  fit1       <- .fixture_a$fit1
  lambda_hat <- lavaan::lavInspect(fit1, "est")$lambda
  theta_hat  <- diag(lavaan::lavInspect(fit1, "est")$theta)

  info <- scfa_factor_information(fit1)

  expect_true(is.numeric(info))
  expect_length(info, 2)
  expect_true(all(info > 0))

  info_manual <- colSums(lambda_hat^2 / theta_hat)
  expect_equal(info, info_manual, tolerance = 1e-12)

  expect_gt(info["f2"], info["f1"])
})

test_that("scfa_factor_information rejects non-lavaan objects", {
  expect_error(scfa_factor_information(list()),   "lavaan")
  expect_error(scfa_factor_information(NULL),     "lavaan")
  expect_error(scfa_factor_information("string"), "lavaan")
})

# ---------------------------------------------------------------------------

test_that("scfa_propagation_variance equals 1 / I_k", {
  fit1 <- .fixture_a$fit1
  info <- scfa_factor_information(fit1)
  psi  <- scfa_propagation_variance(fit1)

  expect_true(all(psi > 0))
  expect_equal(psi, 1 / info, tolerance = 1e-14)
  expect_gt(psi["f1"], psi["f2"])
})

# ---------------------------------------------------------------------------

test_that("scfa_factor_reliability is in (0, 1] and matches formula", {
  fit1     <- .fixture_a$fit1
  phi_hat  <- diag(lavaan::lavInspect(fit1, "est")$psi)
  info     <- scfa_factor_information(fit1)
  rho      <- scfa_factor_reliability(fit1)

  expect_true(all(rho > 0))
  expect_true(all(rho <= 1))

  rho_manual <- phi_hat * info / (phi_hat * info + 1)
  expect_equal(rho, rho_manual, tolerance = 1e-10)
  expect_gt(rho["f2"], rho["f1"])
})

test_that("scfa_factor_reliability = I/(I+1) when phi = 1", {
  fit1_b <- .fixture_b$fit1
  info   <- scfa_factor_information(fit1_b)
  rho    <- scfa_factor_reliability(fit1_b)

  rho_check <- info / (info + 1)
  expect_equal(rho, rho_check, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------

test_that("scfa_propagation_diagnostics returns correct structure", {
  fit1     <- .fixture_a$fit1
  info     <- scfa_factor_information(fit1)
  psi      <- scfa_propagation_variance(fit1)
  rho      <- scfa_factor_reliability(fit1)
  diag_tbl <- scfa_propagation_diagnostics(fit1, threshold = 0.90)

  expect_s3_class(diag_tbl, "scfa_diagnostics")
  expect_s3_class(diag_tbl, "data.frame")
  expect_equal(nrow(diag_tbl), 2)
  expect_named(diag_tbl, c("factor", "n_indicators", "I_k", "psi_nu",
                            "phi_k", "rho_k", "flag"))
  expect_true(all(diag_tbl$n_indicators == 3))
  expect_equal(diag_tbl$I_k,    unname(info[diag_tbl$factor]), tolerance = 1e-14)
  expect_equal(diag_tbl$psi_nu, unname(psi[diag_tbl$factor]),  tolerance = 1e-14)
  expect_equal(diag_tbl$rho_k,  unname(rho[diag_tbl$factor]),  tolerance = 1e-14)

  expect_true( diag_tbl$flag[diag_tbl$factor == "f1"])   # rho_f1 < 0.90
  expect_false(diag_tbl$flag[diag_tbl$factor == "f2"])   # rho_f2 > 0.90
})

test_that("scfa_propagation_diagnostics threshold attribute stored correctly", {
  fit1 <- .fixture_a$fit1
  d    <- scfa_propagation_diagnostics(fit1, threshold = 0.80)
  expect_equal(attr(d, "threshold"), 0.80)
})

test_that("scfa_propagation_diagnostics default threshold = 0.70 flags nothing", {
  fit1    <- .fixture_a$fit1
  d_def   <- scfa_propagation_diagnostics(fit1)
  expect_false(any(d_def$flag))
})

test_that("scfa_propagation_diagnostics rejects invalid thresholds", {
  fit1 <- .fixture_a$fit1
  expect_error(scfa_propagation_diagnostics(fit1, threshold = 1.5))
  expect_error(scfa_propagation_diagnostics(fit1, threshold = 0))
  expect_error(scfa_propagation_diagnostics(fit1, threshold = -0.1))
})

# ---------------------------------------------------------------------------

test_that("scfa_correct_loadings reverses attenuation row-wise", {
  fit1         <- .fixture_a$fit1
  rho          <- scfa_factor_reliability(fit1)
  scores_reg   <- as.data.frame(lavaan::lavPredict(fit1, method = "regression"))
  fit2         <- lavaan::cfa("g =~ f1 + f2", data = scores_reg)
  corrected    <- suppressMessages(scfa_correct_loadings(fit2, fit1))
  lambda2_raw  <- lavaan::lavInspect(fit2, "est")$lambda

  expect_true(is.matrix(corrected))
  expect_equal(dim(corrected), dim(lambda2_raw))

  nz <- lambda2_raw != 0
  expect_true(all(corrected[nz] >= lambda2_raw[nz]))
  expect_true(all(corrected[!nz] == 0))

  rho_vals <- rho[rownames(lambda2_raw)]
  for (k in seq_len(nrow(lambda2_raw))) {
    cols_nz  <- lambda2_raw[k, ] != 0
    expected <- lambda2_raw[k, cols_nz] / rho_vals[k]
    expect_equal(corrected[k, cols_nz], expected, tolerance = 1e-12)
  }
})

test_that("scfa_correct_loadings emits usage message mentioning 'regression'", {
  fit1       <- .fixture_a$fit1
  scores_reg <- as.data.frame(lavaan::lavPredict(fit1, method = "regression"))
  fit2       <- lavaan::cfa("g =~ f1 + f2", data = scores_reg)
  expect_message(scfa_correct_loadings(fit2, fit1), regexp = "regression",
                 ignore.case = TRUE)
})

test_that("scfa_correct_loadings rejects mismatched indicator count", {
  fit1      <- .fixture_a$fit1
  fit_one   <- lavaan::cfa("f1 =~ x1 + x2 + x3", data = .fixture_a$dat)
  scores    <- as.data.frame(lavaan::lavPredict(fit_one, method = "regression"))
  fit2_one  <- lavaan::cfa("g =~ f1", data = scores)
  expect_error(suppressMessages(scfa_correct_loadings(fit2_one, fit1)))
})

test_that("scfa_correct_loadings rejects non-lavaan objects", {
  expect_error(scfa_correct_loadings(list(), list()))
})
