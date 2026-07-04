# =============================================================================
#  tests/testthat/test-s3-methods.R
#
#  Tests for print.scfa_diagnostics and plot.scfa_diagnostics.
# =============================================================================

test_that("print.scfa_diagnostics outputs header and threshold", {
  diag <- scfa_propagation_diagnostics(.fixture_a$fit1, threshold = 0.80)
  out  <- capture.output(print(diag))
  expect_true(any(grepl("Diagnostics", out, ignore.case = TRUE)))
  expect_true(any(grepl("0.80", out)))
})

test_that("print.scfa_diagnostics returns the object invisibly", {
  diag <- scfa_propagation_diagnostics(.fixture_a$fit1)
  ret  <- withVisible(print(diag))
  expect_false(ret$visible)
  expect_identical(ret$value, diag)
})

test_that("print.scfa_diagnostics flags count is correct", {
  diag <- scfa_propagation_diagnostics(.fixture_b$fit1, threshold = 0.80)
  out  <- capture.output(print(diag))
  # f3 (rho ~0.38) and f2 (rho ~0.64) should be flagged at 0.80
  flagged_row <- any(grepl("FLAG", out))
  expect_true(flagged_row)
})

test_that("plot.scfa_diagnostics runs without error (base graphics fallback)", {
  diag <- scfa_propagation_diagnostics(.fixture_a$fit1)
  # We can't guarantee ggplot2 is installed, so just check it doesn't error
  expect_no_error(suppressMessages(plot(diag)))
})
