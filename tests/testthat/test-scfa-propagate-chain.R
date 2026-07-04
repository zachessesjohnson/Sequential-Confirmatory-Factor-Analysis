# =============================================================================
#  tests/testthat/test-scfa-propagate-chain.R
#
#  Unit tests for scfa_propagate_chain().
# =============================================================================

test_that("scfa_propagate_chain with single fit equals scfa_propagation_variance", {
  fit1  <- .fixture_a$fit1
  chain <- scfa_propagate_chain(list(fit1))

  expect_length(chain, 1)
  expect_equal(chain[[1]], scfa_propagation_variance(fit1), tolerance = 1e-14)
})

test_that("scfa_propagate_chain accumulates across two stages", {
  fit1         <- .fixture_c$fit1
  scores_bart  <- as.data.frame(lavaan::lavPredict(fit1, method = "bartlett"))
  fit2         <- suppressWarnings(
    lavaan::cfa("g =~ subfactor1 + subfactor2", data = scores_bart, std.lv = TRUE)
  )

  chain <- scfa_propagate_chain(list(fit1, fit2))

  expect_length(chain, 2)

  psi1 <- scfa_propagation_variance(fit1)
  psi2 <- scfa_propagation_variance(fit2)

  # Stage-1 element equals psi1
  expect_equal(chain[[1]], psi1, tolerance = 1e-14)

  # Stage-2 element = psi1 accumulated + psi2
  # chain[[2]] should equal sum(psi1) + psi2 (single factor in stage 2)
  # Because fit2 has 1 factor and fit1 has 2, the chain adds positionally
  # (stage-2 has 1 factor so only min_len = 1 entry is accumulated)
  expect_length(chain[[2]], 1)
  expect_true(chain[[2]] > 0)
})

test_that("scfa_propagate_chain names stages correctly", {
  fit1  <- .fixture_a$fit1
  chain <- scfa_propagate_chain(list(my_stage = fit1))
  expect_named(chain, "my_stage")
})

test_that("scfa_propagate_chain rejects non-list or empty input", {
  expect_error(scfa_propagate_chain(list()))
  expect_error(scfa_propagate_chain("not a list"))
})

test_that("scfa_propagate_chain rejects non-lavaan elements", {
  fit1 <- .fixture_a$fit1
  expect_error(scfa_propagate_chain(list(fit1, "not_lavaan")))
})

test_that("cumulative variance is always >= per-stage variance", {
  fit1         <- .fixture_c$fit1
  scores_bart  <- as.data.frame(lavaan::lavPredict(fit1, method = "bartlett"))
  fit2         <- suppressWarnings(
    lavaan::cfa("g =~ subfactor1 + subfactor2", data = scores_bart, std.lv = TRUE)
  )
  chain <- scfa_propagate_chain(list(fit1, fit2))
  psi2  <- scfa_propagation_variance(fit2)
  # chain[[2]] >= psi2 (cumulative >= marginal)
  expect_true(all(chain[[2]] >= psi2 - 1e-14))
})
