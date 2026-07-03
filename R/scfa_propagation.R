# =============================================================================
#  scfa_propagation.R
#  Error-propagation diagnostics for Sequential Confirmatory Factor Analysis
#
#  These functions implement the closed-form quantities derived in Section 4
#  ("Quantifying the Propagated Error") of the SCFA paper.  They all operate
#  on a lavaan CFA object fitted to the Stage-1 measurement model.
#
#  Dependencies: lavaan (>= 0.6)
# =============================================================================


# -----------------------------------------------------------------------------
#  Internal helper: extract simple-structure loadings, residual variances, and
#  factor variances from a lavaan fit object.
# -----------------------------------------------------------------------------
.scfa_extract_params <- function(fit) {
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.")
  }

  lambda <- lavaan::lavInspect(fit, "est")$lambda   # p x k loading matrix
  theta  <- diag(lavaan::lavInspect(fit, "est")$theta) # length-p residual variances
  phi    <- diag(lavaan::lavInspect(fit, "est")$psi)   # length-k factor variances

  list(lambda = lambda, theta = theta, phi = phi)
}


# -----------------------------------------------------------------------------
#' Compute per-factor Fisher information from a Stage-1 lavaan fit
#'
#' Under simple structure, the information that the indicators of factor \eqn{k}
#' carry about that factor is
#' \deqn{I_k = \sum_{j \in \mathcal{J}_k} \frac{\lambda_{jk}^2}{\theta_j},}
#' where \eqn{\mathcal{J}_k} is the set of indicators whose primary loading is
#' on factor \eqn{k} (i.e. all non-zero entries in column \eqn{k} of
#' \eqn{\boldsymbol{\Lambda}_1}).
#'
#' @param fit A fitted \code{lavaan} CFA object for the Stage-1 model.
#'
#' @return A named numeric vector of length equal to the number of first-order
#'   factors, giving \eqn{\hat{I}_k} for each factor.
#'
#' @examples
#' \dontrun{
#' library(lavaan)
#' model <- '
#'   f1 =~ x1 + x2 + x3
#'   f2 =~ x4 + x5 + x6
#' '
#' fit <- cfa(model, data = mydata)
#' scfa_factor_information(fit)
#' }
#'
#' @export
scfa_factor_information <- function(fit) {
  p <- .scfa_extract_params(fit)
  lambda <- p$lambda
  theta  <- p$theta

  # I_k = sum_j  lambda_{jk}^2 / theta_j  (over all j with lambda_{jk} != 0)
  # Using matrix operations: colSums(lambda^2 / theta), where theta is recycled
  # row-wise (each row j divided by theta_j).
  info <- colSums(lambda^2 / theta)
  info
}


# -----------------------------------------------------------------------------
#' Compute per-factor propagated estimation-error variance
#'
#' The variance of the factor-score estimation error (Bartlett scores) for
#' factor \eqn{k} is the reciprocal of its Fisher information:
#' \deqn{\hat\psi_{\nu_1,k} = \frac{1}{\hat{I}_k}.}
#'
#' This is the amount of noise that propagates into the Stage-2 residual
#' variance when Bartlett factor scores are used as Stage-2 inputs.
#'
#' @param fit A fitted \code{lavaan} CFA object for the Stage-1 model.
#'
#' @return A named numeric vector of length equal to the number of first-order
#'   factors, giving \eqn{\hat\psi_{\nu_1,k}} for each factor.
#'
#' @examples
#' \dontrun{
#' scfa_propagation_variance(fit_stage1)
#' }
#'
#' @export
scfa_propagation_variance <- function(fit) {
  info <- scfa_factor_information(fit)
  1 / info
}


# -----------------------------------------------------------------------------
#' Compute per-factor score reliability (squared factor-score determinacy)
#'
#' The reliability of the \eqn{k}th factor score is
#' \deqn{\hat\rho_k = \frac{\hat\phi_k \hat{I}_k}{\hat\phi_k \hat{I}_k + 1} \in (0,1],}
#' where \eqn{\hat\phi_k} is the estimated variance of factor \eqn{k} and
#' \eqn{\hat{I}_k} is its Fisher information (see \code{\link{scfa_factor_information}}).
#'
#' \eqn{\hat\rho_k = 1} means the factor score is measured without error;
#' \eqn{\hat\rho_k \approx 0} means it is almost pure noise and propagation
#' will severely distort Stage-2 estimates.
#'
#' @param fit A fitted \code{lavaan} CFA object for the Stage-1 model.
#'
#' @return A named numeric vector of length equal to the number of first-order
#'   factors, giving \eqn{\hat\rho_k} for each factor.
#'
#' @examples
#' \dontrun{
#' scfa_factor_reliability(fit_stage1)
#' }
#'
#' @export
scfa_factor_reliability <- function(fit) {
  p    <- .scfa_extract_params(fit)
  info <- scfa_factor_information(fit)
  phi  <- p$phi

  rho <- (phi * info) / (phi * info + 1)
  rho
}


# -----------------------------------------------------------------------------
#' Summarise all error-propagation diagnostics for a Stage-1 model
#'
#' Returns a data frame with one row per first-order factor, containing the
#' plug-in estimates of Fisher information \eqn{\hat{I}_k}, propagated error
#' variance \eqn{\hat\psi_{\nu_1,k}}, factor variance \eqn{\hat\phi_k}, and
#' score reliability \eqn{\hat\rho_k}.  This diagnostic can be computed
#' entirely from the Stage-1 output, before Stage 2 is ever run.
#'
#' A rule of thumb: factors with \eqn{\hat\rho_k < 0.70} are the weak links
#' in the propagation chain and warrant careful attention.
#'
#' @param fit A fitted \code{lavaan} CFA object for the Stage-1 model.
#' @param threshold Numeric scalar in (0, 1).  Factors with
#'   \eqn{\hat\rho_k} below this value are flagged in a \code{flag} column.
#'   Default is \code{0.70}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{factor}{Factor name.}
#'     \item{n_indicators}{Number of indicators loading on the factor.}
#'     \item{I_k}{Fisher information \eqn{\hat{I}_k}.}
#'     \item{psi_nu}{Propagated error variance \eqn{\hat\psi_{\nu_1,k} = 1/\hat{I}_k}.}
#'     \item{phi_k}{Estimated factor variance \eqn{\hat\phi_k}.}
#'     \item{rho_k}{Score reliability \eqn{\hat\rho_k}.}
#'     \item{flag}{Logical: \code{TRUE} if \eqn{\hat\rho_k < \code{threshold}}.}
#'   }
#'
#' @examples
#' \dontrun{
#' diag_tbl <- scfa_propagation_diagnostics(fit_stage1)
#' print(diag_tbl)
#' }
#'
#' @export
scfa_propagation_diagnostics <- function(fit, threshold = 0.70) {
  if (!is.numeric(threshold) || length(threshold) != 1 ||
      threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a single number strictly between 0 and 1.")
  }

  p      <- .scfa_extract_params(fit)
  lambda <- p$lambda
  info   <- scfa_factor_information(fit)
  phi    <- p$phi
  rho    <- scfa_factor_reliability(fit)
  psi_nu <- 1 / info

  # Count non-zero (primary) loadings per factor
  n_ind <- colSums(lambda != 0)

  data.frame(
    factor      = names(info),
    n_indicators = as.integer(n_ind),
    I_k         = unname(info),
    psi_nu      = unname(psi_nu),
    phi_k       = unname(phi),
    rho_k       = unname(rho),
    flag        = unname(rho < threshold),
    stringsAsFactors = FALSE,
    row.names   = NULL
  )
}


# -----------------------------------------------------------------------------
#' Correct Stage-2 loadings for attenuation due to regression-score shrinkage
#'
#' When regression factor scores are used as Stage-2 inputs, the estimated
#' loadings \eqn{\tilde\lambda_{2,k}} are attenuated by the factor-score
#' reliability:
#' \deqn{\tilde\lambda_{2,k} \approx \rho_k \lambda_{2,k,0}.}
#' This function reverses the attenuation:
#' \deqn{\hat\lambda_{2,k}^{\text{corr}} = \frac{\hat{\tilde\lambda}_{2,k}}{\hat\rho_k},}
#' producing estimates that match the Bartlett-score (unbiased) loadings.
#'
#' @param fit2 A fitted \code{lavaan} CFA object for the Stage-2 model.
#'   Its latent variables must correspond (in order) to the first-order
#'   factors of \code{fit1}.
#' @param fit1 A fitted \code{lavaan} CFA object for the Stage-1 model.
#'   Used to compute \eqn{\hat\rho_k} via \code{\link{scfa_factor_reliability}}.
#'
#' @return A matrix of the same dimensions as the Stage-2 loading matrix,
#'   with each column \eqn{k} divided by \eqn{\hat\rho_k}.  Entries that
#'   were exactly zero in the original matrix remain zero (cross-loadings
#'   fixed to zero are not corrected).
#'
#' @details
#' The correction is valid only when regression factor scores were used at
#' Stage 1.  Under Bartlett scores the loadings are already asymptotically
#' unbiased and no correction is needed.  The function issues a message
#' reminding the analyst of this distinction.
#'
#' @examples
#' \dontrun{
#' corrected_lambda <- scfa_correct_loadings(fit_stage2, fit_stage1)
#' }
#'
#' @export
scfa_correct_loadings <- function(fit2, fit1) {
  if (!inherits(fit2, "lavaan") || !inherits(fit1, "lavaan")) {
    stop("Both 'fit1' and 'fit2' must be fitted lavaan objects.")
  }

  message(
    "Note: scfa_correct_loadings() assumes regression factor scores were used\n",
    "at Stage 1.  Under Bartlett scores the loadings are already unbiased;\n",
    "applying this correction in that case would over-correct."
  )

  lambda2 <- lavaan::lavInspect(fit2, "est")$lambda
  rho     <- scfa_factor_reliability(fit1)

  k2 <- ncol(lambda2)
  k1 <- length(rho)
  if (k2 != k1) {
    stop(
      "The Stage-2 model has ", k2, " latent variable(s) but the Stage-1 model ",
      "has ", k1, " factor(s).  They must match."
    )
  }

  # Divide each column k by rho_k, but only for non-zero entries.
  corrected <- lambda2
  for (k in seq_len(k2)) {
    nz            <- lambda2[, k] != 0
    corrected[nz, k] <- lambda2[nz, k] / rho[k]
  }

  corrected
}
