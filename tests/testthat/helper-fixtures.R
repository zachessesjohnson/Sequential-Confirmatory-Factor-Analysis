# =============================================================================
#  tests/testthat/helper-fixtures.R
#
#  Shared lavaan fixtures used across the SCFA test suite.
#  Loaded automatically by testthat before any test file.
# =============================================================================

library(lavaan)

# ---------------------------------------------------------------------------
# Fixture A – two-factor simple-structure model
#   f1: loadings ~0.7, residuals ~0.51  -> lower reliability
#   f2: loadings ~0.9, residuals ~0.19  -> higher reliability
# ---------------------------------------------------------------------------
.fixture_a <- local({
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
  dat  <- simulateData(pop_model, sample.nobs = 500, model.type = "sem")
  fit1 <- cfa("f1 =~ x1 + x2 + x3\nf2 =~ x4 + x5 + x6", data = dat)
  list(dat = dat, fit1 = fit1)
})

# ---------------------------------------------------------------------------
# Fixture B – three-factor simple-structure model (std.lv = TRUE)
#   f1: loadings 0.8, rho ~0.83   (high)
#   f2: loadings 0.6, rho ~0.64   (moderate)
#   f3: loadings 0.4, rho ~0.38   (low)
# ---------------------------------------------------------------------------
.fixture_b <- local({
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
  list(dat = dat_a, fit1 = fit1_a)
})

# ---------------------------------------------------------------------------
# Fixture C – README two-factor workflow (subfactor1, subfactor2)
# ---------------------------------------------------------------------------
.fixture_c <- local({
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
  fit1_b <- cfa("subfactor1 =~ item1+item2+item3\nsubfactor2 =~ item4+item5+item6",
                data = lower_data, std.lv = TRUE)
  list(dat = lower_data, fit1 = fit1_b)
})
