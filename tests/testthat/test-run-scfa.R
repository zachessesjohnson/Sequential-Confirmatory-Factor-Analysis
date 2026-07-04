# =============================================================================
#  tests/testthat/test-run-scfa.R
#
#  Integration tests for the run_scfa() orchestration function.
# =============================================================================

# ---------------------------------------------------------------------------
# Two-stage Bartlett workflow
# ---------------------------------------------------------------------------
test_that("run_scfa two-stage Bartlett workflow produces valid result", {
  result <- run_scfa(
    stage_models = list(
      "subfactor1 =~ item1 + item2 + item3\nsubfactor2 =~ item4 + item5 + item6",
      "higher_factor =~ subfactor1 + subfactor2"
    ),
    data      = .fixture_c$dat,
    method    = "bartlett",
    threshold = 0.70,
    cfa_args  = list(std.lv = TRUE)
  )

  expect_s3_class(result, "scfa_result")

  # Fits
  expect_length(result$fits, 2)
  expect_true(all(vapply(result$fits, inherits, logical(1), what = "lavaan")))
  expect_true(result$fits$stage_1@optim$converged)

  # Diagnostics
  expect_length(result$diagnostics, 2)
  expect_s3_class(result$diagnostics$stage_1, "scfa_diagnostics")
  expect_equal(nrow(result$diagnostics$stage_1), 2)

  # Scores passed between stages
  expect_s3_class(result$scores$stage_1, "data.frame")
  expect_equal(nrow(result$scores$stage_1), 400)

  # Index scores
  expect_equal(nrow(result$index_scores), 400)

  # Correction (Bartlett -> scfa_correct_residuals)
  expect_true(!is.null(result$correction))
  expect_true(is.numeric(result$correction))
  expect_true(all(result$correction >= 0))

  # Propagation chain
  expect_length(result$propagation_chain, 2)
})

# ---------------------------------------------------------------------------
# Two-stage regression workflow
# ---------------------------------------------------------------------------
test_that("run_scfa two-stage regression workflow applies loading correction", {
  result <- run_scfa(
    stage_models = list(
      "subfactor1 =~ item1 + item2 + item3\nsubfactor2 =~ item4 + item5 + item6",
      "higher_factor =~ subfactor1 + subfactor2"
    ),
    data      = .fixture_c$dat,
    method    = "regression",
    threshold = 0.70,
    correct   = TRUE,
    cfa_args  = list(std.lv = TRUE)
  )

  expect_s3_class(result, "scfa_result")
  expect_true(is.matrix(result$correction))
})

# ---------------------------------------------------------------------------
# correct = FALSE
# ---------------------------------------------------------------------------
test_that("run_scfa with correct = FALSE returns NULL correction", {
  result <- run_scfa(
    stage_models = list(
      "subfactor1 =~ item1 + item2 + item3\nsubfactor2 =~ item4 + item5 + item6",
      "higher_factor =~ subfactor1 + subfactor2"
    ),
    data    = .fixture_c$dat,
    method  = "bartlett",
    correct = FALSE,
    cfa_args = list(std.lv = TRUE)
  )
  expect_null(result$correction)
})

# ---------------------------------------------------------------------------
# Single-stage workflow
# ---------------------------------------------------------------------------
test_that("run_scfa single-stage workflow works", {
  result <- run_scfa(
    stage_models = list(
      "subfactor1 =~ item1 + item2 + item3\nsubfactor2 =~ item4 + item5 + item6"
    ),
    data    = .fixture_c$dat,
    method  = "bartlett",
    correct = TRUE,
    cfa_args = list(std.lv = TRUE)
  )
  expect_length(result$fits, 1)
  expect_null(result$correction)   # correct = TRUE but only 1 stage
  expect_equal(nrow(result$index_scores), 400)
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
test_that("run_scfa rejects empty stage_models", {
  expect_error(run_scfa(list(), data = .fixture_c$dat))
})

test_that("run_scfa rejects invalid threshold", {
  expect_error(run_scfa(
    list("subfactor1 =~ item1 + item2 + item3"),
    data      = .fixture_c$dat,
    threshold = 1.5
  ))
})

test_that("run_scfa rejects blank model string", {
  expect_error(run_scfa(list("   "), data = .fixture_c$dat))
})

# ---------------------------------------------------------------------------
# print.scfa_result
# ---------------------------------------------------------------------------
test_that("print.scfa_result outputs without error", {
  result <- run_scfa(
    stage_models = list(
      "subfactor1 =~ item1 + item2 + item3\nsubfactor2 =~ item4 + item5 + item6",
      "higher_factor =~ subfactor1 + subfactor2"
    ),
    data    = .fixture_c$dat,
    method  = "bartlett",
    cfa_args = list(std.lv = TRUE)
  )
  expect_output(print(result), "Sequential CFA Result")
  expect_output(print(result), "bartlett")
})
