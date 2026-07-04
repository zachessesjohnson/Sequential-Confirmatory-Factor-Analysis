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

  est    <- lavaan::lavInspect(fit, "est")
  lambda <- est$lambda                  # p x k loading matrix
  theta  <- diag(est$theta)            # length-p residual variances
  phi    <- diag(est$psi)              # length-k factor variances

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
#' @importFrom lavaan lavInspect
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

  data_out <- data.frame(
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
  attr(data_out, "threshold") <- threshold
  class(data_out) <- c("scfa_diagnostics", "data.frame")
  data_out
}


# -----------------------------------------------------------------------------
#' Correct Stage-2 loadings for attenuation due to regression-score shrinkage
#'
#' When regression factor scores are used as Stage-2 inputs, the estimated
#' loadings \eqn{\tilde\lambda_{2,k}} are attenuated by the factor-score
#' reliability:
#' \deqn{\tilde\lambda_{2,k} \approx \rho_k \lambda_{2,k,0}.}
#' Here \eqn{k} indexes the Stage-1 factors, each of which appears as an
#' *observed indicator* (a row) in the Stage-2 loading matrix.  This function
#' reverses the attenuation row-wise:
#' \deqn{\hat\lambda_{2,k}^{\text{corr}} = \frac{\hat{\tilde\lambda}_{2,k}}{\hat\rho_k},}
#' producing estimates that match the Bartlett-score (unbiased) loadings.
#'
#' @param fit2 A fitted \code{lavaan} CFA object for the Stage-2 model.
#'   Its observed variables must be the Stage-1 factor scores, one per
#'   first-order factor of \code{fit1} (in the same order).
#' @param fit1 A fitted \code{lavaan} CFA object for the Stage-1 model.
#'   Used to compute \eqn{\hat\rho_k} via \code{\link{scfa_factor_reliability}}.
#'
#' @return A matrix of the same dimensions as the Stage-2 loading matrix,
#'   with each row \eqn{k} divided by \eqn{\hat\rho_k}.  Entries that
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

  # lambda2 is (p2 x k2): rows are Stage-1 factor scores used as Stage-2
  # observed indicators; columns are Stage-2 latent factors.
  # The attenuation lambda_2k ~ rho_k * lambda_2k,0 is indexed by the
  # Stage-1 factor k, i.e. by *row* of lambda2.
  p2 <- nrow(lambda2)
  k1 <- length(rho)
  if (p2 != k1) {
    stop(
      "The Stage-2 loading matrix has ", p2, " row(s) (observed indicator(s)) ",
      "but the Stage-1 model has ", k1, " factor(s).  They must match: each ",
      "Stage-1 factor score should appear as one observed variable in Stage 2."
    )
  }

  # Divide each row k by rho_k, but only for non-zero entries.
  corrected <- lambda2
  for (k in seq_len(p2)) {
    nz              <- lambda2[k, ] != 0
    corrected[k, nz] <- lambda2[k, nz] / rho[k]
  }

  corrected
}


# -----------------------------------------------------------------------------
#' Correct Stage-2 residual variances for Bartlett-score propagation error
#'
#' When **Bartlett** factor scores are used as Stage-2 inputs, the Stage-2
#' loading estimates are asymptotically unbiased, but each Stage-2 residual
#' variance is inflated by the propagated estimation-error variance of the
#' corresponding Stage-1 factor:
#' \deqn{\hat\theta_{2,k}^{\text{adj}} = \hat\theta_{2,k} - \hat\psi_{\nu_1,k},}
#' where \eqn{\hat\psi_{\nu_1,k} = 1/\hat{I}_k} is the propagation variance
#' returned by \code{\link{scfa_propagation_variance}}.
#'
#' The adjustment is applied only to non-zero diagonal entries of
#' \eqn{\boldsymbol{\Theta}_2}.  Cross-indicator residual covariances (off-diagonal
#' entries) are left unchanged.
#'
#' @param fit2 A fitted \code{lavaan} CFA object for the Stage-2 model.
#'   Its observed variables must be the Stage-1 Bartlett factor scores, one per
#'   first-order factor of \code{fit1} (in the same order).
#' @param fit1 A fitted \code{lavaan} CFA object for the Stage-1 model.
#'   Used to compute \eqn{\hat\psi_{\nu_1,k}} via
#'   \code{\link{scfa_propagation_variance}}.
#'
#' @return A numeric vector of length equal to the number of Stage-1 factors,
#'   giving the adjusted Stage-2 residual variances
#'   \eqn{\hat\theta_{2,k}^{\text{adj}}}.  Entries that are already zero or
#'   would become negative after adjustment are clamped to zero, with a warning.
#'
#' @details
#' This correction is intended for **Bartlett** factor scores only.  Under
#' regression scores the loadings are already attenuated; use
#' \code{\link{scfa_correct_loadings}} instead (and note that the residual
#' inflation is then absorbed into the attenuation bias).  The function emits
#' a message reminding the analyst of this distinction.
#'
#' @examples
#' \dontrun{
#' scores_bart <- as.data.frame(lavPredict(fit_stage1, method = "bartlett"))
#' fit_stage2  <- cfa(model_stage2, data = scores_bart)
#' adj_theta   <- scfa_correct_residuals(fit_stage2, fit_stage1)
#' }
#'
#' @importFrom lavaan lavInspect
#' @export
scfa_correct_residuals <- function(fit2, fit1) {
  if (!inherits(fit2, "lavaan") || !inherits(fit1, "lavaan")) {
    stop("Both 'fit1' and 'fit2' must be fitted lavaan objects.")
  }

  message(
    "Note: scfa_correct_residuals() assumes Bartlett factor scores were used\n",
    "at Stage 1.  Under regression scores the loading attenuation absorbs\n",
    "the residual inflation; use scfa_correct_loadings() instead."
  )

  theta2  <- diag(lavaan::lavInspect(fit2, "est")$theta)  # Stage-2 residual variances
  psi_nu  <- scfa_propagation_variance(fit1)               # 1/I_k per Stage-1 factor

  p2 <- length(theta2)
  k1 <- length(psi_nu)
  if (p2 != k1) {
    stop(
      "The Stage-2 residual vector has ", p2, " element(s) but the Stage-1 ",
      "model has ", k1, " factor(s).  They must match: each Stage-1 factor ",
      "score should appear as one observed variable in Stage 2."
    )
  }

  # Align psi_nu to theta2 order if names are available
  if (!is.null(names(theta2)) && !is.null(names(psi_nu))) {
    common <- intersect(names(theta2), names(psi_nu))
    if (length(common) == p2) {
      psi_nu <- psi_nu[names(theta2)]
    }
  }

  adj <- theta2 - psi_nu

  # Clamp negative adjustments to zero with a warning
  neg <- which(adj < 0)
  if (length(neg) > 0) {
    warning(
      "Adjusted residual variance(s) for factor(s) [",
      paste(names(adj)[neg], collapse = ", "),
      "] would be negative after subtracting psi_nu; clamped to 0. ",
      "This may indicate the Stage-1 model is poorly identified or that ",
      "regression scores (not Bartlett) were used."
    )
    adj[neg] <- 0
  }

  adj
}


# -----------------------------------------------------------------------------
#' Accumulate propagation variances across a multi-stage sequential CFA chain
#'
#' In a hierarchy with more than two stages, estimation error accumulates at
#' each stage.  Given a list of fitted \code{lavaan} objects
#' \eqn{(\text{fit}_1, \text{fit}_2, \ldots, \text{fit}_S)}, this function
#' sums the propagated error variances stage by stage:
#' \deqn{\Psi_{\nu,k}^{(s)} = \sum_{t=1}^{s} \psi_{\nu,k}^{(t)},}
#' where \eqn{\psi_{\nu,k}^{(t)} = 1/I_k^{(t)}} is the propagation variance
#' contributed by stage \eqn{t}.
#'
#' This is the key diagnostic for three- or higher-level hierarchies (e.g.
#' item \eqn{\to} sub-factor \eqn{\to} factor \eqn{\to} index), where a
#' practitioner needs to know how much total noise has accumulated by the time
#' Stage-\eqn{S} inputs are formed.
#'
#' @param fits A named or unnamed list of fitted \code{lavaan} CFA objects, one
#'   per stage, ordered from lowest (Stage 1) to highest (Stage S) level.
#'   Each element must be a \code{lavaan} object.
#'
#' @return A list with one element per stage, each being a named numeric vector
#'   of cumulative propagation variances for the factors estimated at that
#'   stage.  The final element therefore contains the total accumulated
#'   propagation variance entering the final stage's inputs.
#'
#' @examples
#' \dontrun{
#' chain <- scfa_propagate_chain(list(fit_stage1, fit_stage2))
#' # Cumulative propagation variance after Stage 1:
#' chain[[1]]
#' # Cumulative propagation variance after Stage 2:
#' chain[[2]]
#' }
#'
#' @importFrom lavaan lavInspect
#' @export
scfa_propagate_chain <- function(fits) {
  if (!is.list(fits) || length(fits) < 1) {
    stop("'fits' must be a non-empty list of lavaan objects.")
  }
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "lavaan")) {
      stop("Element ", i, " of 'fits' is not a lavaan object.")
    }
  }

  stage_names <- names(fits)
  if (is.null(stage_names)) {
    stage_names <- paste0("stage_", seq_along(fits))
  }

  result <- vector("list", length(fits))
  names(result) <- stage_names

  cumulative <- NULL
  for (s in seq_along(fits)) {
    psi_s <- scfa_propagation_variance(fits[[s]])
    if (is.null(cumulative)) {
      cumulative <- psi_s
    } else {
      # Attempt to align by name; fall back to positional alignment
      if (!is.null(names(cumulative)) && !is.null(names(psi_s)) &&
          length(intersect(names(cumulative), names(psi_s))) == length(psi_s)) {
        cumulative <- cumulative[names(psi_s)] + psi_s
      } else {
        if (length(psi_s) != length(cumulative)) {
          warning(
            "Stage ", s, " has ", length(psi_s), " factor(s) but the ",
            "cumulative vector has ", length(cumulative), " element(s); ",
            "using positional alignment."
          )
          min_len     <- min(length(psi_s), length(cumulative))
          cumulative  <- cumulative[seq_len(min_len)] + psi_s[seq_len(min_len)]
        } else {
          cumulative <- cumulative + psi_s
        }
      }
    }
    result[[s]] <- cumulative
  }

  result
}


# -----------------------------------------------------------------------------
#' Print method for \code{scfa_diagnostics} objects
#'
#' Prints the diagnostic table returned by
#' \code{\link{scfa_propagation_diagnostics}}, with flagged rows clearly
#' marked.
#'
#' @param x An object of class \code{scfa_diagnostics}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.scfa_diagnostics <- function(x, ...) {
  thr <- attr(x, "threshold")
  cat("Sequential CFA – Error-Propagation Diagnostics\n")
  if (!is.null(thr)) {
    cat(sprintf("Reliability threshold: %.2f  (flagged if rho_k < threshold)\n", thr))
  }
  cat("\n")

  # Build a display copy with a ✓/✗ column replacing the logical flag
  disp            <- x
  disp$flag       <- ifelse(x$flag, "FLAG", "ok")
  disp$I_k        <- round(disp$I_k,   3)
  disp$psi_nu     <- round(disp$psi_nu, 4)
  disp$phi_k      <- round(disp$phi_k,  3)
  disp$rho_k      <- round(disp$rho_k,  3)

  print(as.data.frame(disp), row.names = FALSE)

  n_flagged <- sum(x$flag)
  if (n_flagged > 0) {
    cat(sprintf(
      "\n%d factor(s) flagged as weak links (rho_k < %.2f).\n",
      n_flagged, thr
    ))
  } else {
    cat(sprintf("\nAll factors meet the reliability threshold (%.2f).\n", thr))
  }

  invisible(x)
}


# -----------------------------------------------------------------------------
#' Plot method for \code{scfa_diagnostics} objects
#'
#' Produces a bar chart of per-factor score reliability (\eqn{\hat\rho_k})
#' with a horizontal line at the diagnostic threshold.  Bars below the
#' threshold are coloured red to draw attention to weak links.
#'
#' @param x An object of class \code{scfa_diagnostics}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot2} object (invisibly), or base-graphics output if
#'   \pkg{ggplot2} is not available.
#'
#' @export
plot.scfa_diagnostics <- function(x, ...) {
  thr <- attr(x, "threshold")
  if (is.null(thr)) thr <- 0.70

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(
      data = x,
      mapping = ggplot2::aes(
        x    = factor(.data$factor, levels = .data$factor),
        y    = .data$rho_k,
        fill = .data$flag
      )
    ) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_hline(yintercept = thr, linetype = "dashed", colour = "black") +
      ggplot2::annotate(
        "text",
        x     = 0.5, y = thr + 0.02,
        label = paste0("threshold = ", thr),
        hjust = 0, size = 3
      ) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
        labels = c("FALSE" = "OK", "TRUE" = "Flagged"),
        name   = NULL
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::labs(
        title    = "Sequential CFA – Factor Score Reliability",
        subtitle = sprintf("Dashed line = reliability threshold (%.2f)", thr),
        x        = "Factor",
        y        = expression(hat(rho)[k])
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "none")

    print(p)
    invisible(p)
  } else {
    # Fall back to base graphics
    bar_cols <- ifelse(x$flag, "firebrick", "steelblue")
    bp <- barplot(
      x$rho_k,
      names.arg = x$factor,
      col       = bar_cols,
      ylim      = c(0, 1),
      ylab      = expression(hat(rho)[k]),
      xlab      = "Factor",
      main      = "Sequential CFA – Factor Score Reliability"
    )
    abline(h = thr, lty = 2)
    invisible(bp)
  }
}
