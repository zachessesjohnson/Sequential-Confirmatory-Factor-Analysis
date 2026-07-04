# =============================================================================
#  tests/testthat/test-scfa-correct-residuals.R
#
#  Unit tests for scfa_correct_residuals().
# =============================================================================

test_that("scfa_correct_residuals returns a named numeric vector", {
  fit1         <- .fixture_c$fit1
  scores_bart  <- as.data.frame(lavaan::lavPredict(fit1, method = "bartlett"))
  fit2         <- suppressWarnings(
    lavaan::cfa("g =~ subfactor1 + subfactor2", data = scores_bart, std.lv = TRUE)
  )
  adj <- suppressMessages(scfa_correct_residuals(fit2, fit1))

  expect_true(is.numeric(adj))
  expect_length(adj, 2)
  expect_true(all(!is.na(adj)))
})

test_that("scfa_correct_residuals subtracts 1/I_k from each residual", {
  fit1        <- .fixture_c$fit1
  scores_bart <- as.data.frame(lavaan::lavPredict(fit1, method = "bartlett"))
  fit2        <- suppressWarnings(
    lavaan::cfa("g =~ subfactor1 + subfactor2", data = scores_bart, std.lv = TRUE)
  )
  theta2  <- diag(lavaan::lavInspect(fit2, "est")$theta)
  psi_nu  <- scfa_propagation_variance(fit1)
  adj     <- suppressMessages(scfa_correct_residuals(fit2, fit1))

  # Align by name
  psi_nu_ord <- psi_nu[names(theta2)]
  expected   <- pmax(theta2 - psi_nu_ord, 0)
  expect_equal(adj, expected, tolerance = 1e-14)
})

test_that("scfa_correct_residuals adjusted values are non-negative", {
  fit1        <- .fixture_c$fit1
  scores_bart <- as.data.frame(lavaan::lavPredict(fit1, method = "bartlett"))
  fit2        <- suppressWarnings(
    lavaan::cfa("g =~ subfactor1 + subfactor2", data = scores_bart, std.lv = TRUE)
  )
  adj <- suppressMessages(scfa_correct_residuals(fit2, fit1))
  expect_true(all(adj >= 0))
})

test_that("scfa_correct_residuals emits message mentioning 'Bartlett'", {
  fit1        <- .fixture_c$fit1
  scores_bart <- as.data.frame(lavaan::lavPredict(fit1, method = "bartlett"))
  fit2        <- suppressWarnings(
    lavaan::cfa("g =~ subfactor1 + subfactor2", data = scores_bart, std.lv = TRUE)
  )
  expect_message(scfa_correct_residuals(fit2, fit1), regexp = "Bartlett",
                 ignore.case = TRUE)
})

test_that("scfa_correct_residuals rejects mismatched stage sizes", {
  fit1      <- .fixture_a$fit1   # 2 factors
  fit_one   <- lavaan::cfa("f1 =~ x1 + x2 + x3", data = .fixture_a$dat)
  scores    <- as.data.frame(lavaan::lavPredict(fit_one, method = "bartlett"))
  fit2_one  <- lavaan::cfa("g =~ f1", data = scores)
  # fit2_one has 1 indicator but fit1 has 2 factors
  expect_error(suppressMessages(scfa_correct_residuals(fit2_one, fit1)))
})

test_that("scfa_correct_residuals rejects non-lavaan objects", {
  expect_error(scfa_correct_residuals(NULL, NULL))
})
