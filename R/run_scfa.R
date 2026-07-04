# =============================================================================
#  run_scfa.R
#  High-level orchestration for the full Sequential CFA workflow
# =============================================================================


#' Run the full Sequential Confirmatory Factor Analysis workflow
#'
#' Orchestrates the complete Sequential CFA pipeline: fits each stage of a
#' hierarchical CFA model in turn, extracts factor scores to use as inputs for
#' the next stage, computes error-propagation diagnostics at every stage, and
#' returns a structured list of results.
#'
#' @details
#' **Workflow for a two-stage hierarchy:**
#' 1. Fit `stage_models[[1]]` to `data` → extract factor scores.
#' 2. Run propagation diagnostics on the Stage-1 fit.
#' 3. Fit `stage_models[[2]]` to the Stage-1 factor scores.
#' 4. (If `correct` is `TRUE` and `method = "regression"`) Correct Stage-2
#'    loadings for attenuation bias via \code{\link{scfa_correct_loadings}}.
#' 5. (If `correct` is `TRUE` and `method = "bartlett"`) Adjust Stage-2
#'    residual variances for propagated error via
#'    \code{\link{scfa_correct_residuals}}.
#' 6. Return all fitted models, scores, diagnostics, and (optionally)
#'    corrections.
#'
#' For hierarchies with three or more stages, steps 1–3 are repeated for each
#' stage and the cumulative propagation chain is computed via
#' \code{\link{scfa_propagate_chain}}.
#'
#' @param stage_models A non-empty character vector or list of character
#'   strings, each a valid **lavaan** model syntax string.  Element `s`
#'   specifies the CFA model for stage `s`.  The first element is fitted to
#'   `data`; subsequent elements are fitted to the factor scores from the
#'   previous stage.
#' @param data A \code{data.frame} (or coercible object) containing the
#'   observed indicators for the first stage.
#' @param method Character scalar: the factor-score extraction method to use at
#'   every stage.  Passed to \code{lavaan::lavPredict()}.  One of
#'   `"bartlett"` (default, asymptotically unbiased) or `"regression"`
#'   (shrunk toward zero; use \code{correct = TRUE} to de-attenuate loadings).
#' @param threshold Numeric scalar in (0, 1).  The reliability threshold passed
#'   to \code{\link{scfa_propagation_diagnostics}} at each stage.  Default
#'   `0.70`.
#' @param correct Logical.  If `TRUE` (default), apply the appropriate
#'   post-hoc correction at the final stage: \code{\link{scfa_correct_loadings}}
#'   for regression scores or \code{\link{scfa_correct_residuals}} for Bartlett
#'   scores.
#' @param cfa_args A named list of additional arguments forwarded to
#'   \code{lavaan::cfa()} at every stage (e.g. `list(std.lv = TRUE)`).
#'   Default `list()`.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{`fits`}{A named list of fitted \code{lavaan} objects, one per
#'       stage (`stage_1`, `stage_2`, …).}
#'     \item{`scores`}{A named list of \code{data.frame}s of factor scores, one
#'       per stage.  `scores[[s]]` are the scores extracted from `fits[[s]]`
#'       and used as inputs to stage `s + 1`.}
#'     \item{`diagnostics`}{A named list of \code{scfa_diagnostics} tables
#'       (see \code{\link{scfa_propagation_diagnostics}}), one per stage.}
#'     \item{`propagation_chain`}{The output of
#'       \code{\link{scfa_propagate_chain}} applied to `fits`: a list of
#'       cumulative propagation variances, one element per stage.}
#'     \item{`index_scores`}{A matrix or \code{data.frame} of final-stage
#'       factor scores (the output of the last stage's
#'       \code{lavaan::lavPredict()}).  These are the composite index scores.}
#'     \item{`correction`}{If `correct = TRUE`, the result of the appropriate
#'       correction function applied to the final two stages.  `NULL` if
#'       `correct = FALSE` or if there is only one stage.}
#'     \item{`method`}{The `method` argument as supplied.}
#'     \item{`threshold`}{The `threshold` argument as supplied.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(lavaan)
#'
#' # Simulate a two-stage hierarchy
#' set.seed(42)
#' pop <- '
#'   subfactor1 =~ 0.7*x1 + 0.7*x2 + 0.7*x3
#'   subfactor2 =~ 0.8*x4 + 0.8*x5 + 0.8*x6
#'   subfactor1 ~~ 1*subfactor1; subfactor2 ~~ 1*subfactor2
#'   subfactor1 ~~ 0*subfactor2
#' '
#' dat <- lavaan::simulateData(pop, sample.nobs = 300)
#'
#' result <- run_scfa(
#'   stage_models = list(
#'     "subfactor1 =~ x1 + x2 + x3
#'      subfactor2 =~ x4 + x5 + x6",
#'     "index =~ subfactor1 + subfactor2"
#'   ),
#'   data      = dat,
#'   method    = "bartlett",
#'   threshold = 0.70
#' )
#'
#' # Propagation diagnostics for Stage 1
#' print(result$diagnostics$stage_1)
#'
#' # Final composite index scores
#' head(result$index_scores)
#' }
#'
#' @importFrom lavaan cfa lavPredict
#' @export
run_scfa <- function(stage_models,
                     data,
                     method    = c("bartlett", "regression"),
                     threshold = 0.70,
                     correct   = TRUE,
                     cfa_args  = list()) {

  # ---- Input validation ------------------------------------------------------
  method <- match.arg(method)

  if (!is.list(stage_models)) {
    stage_models <- as.list(stage_models)
  }
  if (length(stage_models) < 1) {
    stop("'stage_models' must contain at least one model string.")
  }
  for (i in seq_along(stage_models)) {
    if (!is.character(stage_models[[i]]) || nchar(trimws(stage_models[[i]])) == 0) {
      stop("'stage_models[[", i, "]]' must be a non-empty character string.")
    }
  }

  if (!is.data.frame(data) && !is.matrix(data)) {
    data <- as.data.frame(data)
  }
  if (!is.numeric(threshold) || length(threshold) != 1 ||
      threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a single number strictly between 0 and 1.")
  }

  n_stages    <- length(stage_models)
  stage_labels <- paste0("stage_", seq_len(n_stages))

  fits        <- vector("list", n_stages)
  scores_list <- vector("list", n_stages)
  diag_list   <- vector("list", n_stages)
  names(fits)        <- stage_labels
  names(scores_list) <- stage_labels
  names(diag_list)   <- stage_labels

  current_data <- data

  # ---- Fit each stage --------------------------------------------------------
  for (s in seq_len(n_stages)) {
    fit_s <- do.call(
      lavaan::cfa,
      c(list(model = stage_models[[s]], data = current_data), cfa_args)
    )

    if (!fit_s@optim$converged) {
      warning(
        "Stage ", s, " CFA model did not converge.  Results may be unreliable."
      )
    }

    fits[[s]] <- fit_s

    diag_list[[s]] <- scfa_propagation_diagnostics(fit_s, threshold = threshold)

    if (s < n_stages) {
      scores_s      <- as.data.frame(lavaan::lavPredict(fit_s, method = method))
      scores_list[[s]] <- scores_s
      current_data  <- scores_s
    }
  }

  # ---- Final-stage index scores ----------------------------------------------
  index_scores <- lavaan::lavPredict(fits[[n_stages]], method = method)

  # ---- Correction (optional) -------------------------------------------------
  correction <- NULL
  if (correct && n_stages >= 2) {
    fit_last     <- fits[[n_stages]]
    fit_secondlast <- fits[[n_stages - 1]]

    correction <- if (method == "regression") {
      suppressMessages(scfa_correct_loadings(fit_last, fit_secondlast))
    } else {
      suppressMessages(scfa_correct_residuals(fit_last, fit_secondlast))
    }
  }

  # ---- Propagation chain -----------------------------------------------------
  prop_chain <- scfa_propagate_chain(fits)

  # ---- Return ----------------------------------------------------------------
  structure(
    list(
      fits              = fits,
      scores            = scores_list,
      diagnostics       = diag_list,
      propagation_chain = prop_chain,
      index_scores      = index_scores,
      correction        = correction,
      method            = method,
      threshold         = threshold
    ),
    class = "scfa_result"
  )
}


# -----------------------------------------------------------------------------
#' Print method for \code{scfa_result} objects
#'
#' Prints a concise summary of a Sequential CFA result returned by
#' \code{\link{run_scfa}}.
#'
#' @param x An object of class \code{scfa_result}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.scfa_result <- function(x, ...) {
  n_stages <- length(x$fits)
  cat("Sequential CFA Result\n")
  cat("=====================\n")
  cat(sprintf("Stages      : %d\n", n_stages))
  cat(sprintf("Score method: %s\n", x$method))
  cat(sprintf("Threshold   : %.2f\n", x$threshold))
  cat("\n")

  for (s in seq_len(n_stages)) {
    lbl  <- names(x$fits)[s]
    diag <- x$diagnostics[[s]]
    converged <- x$fits[[s]]@optim$converged
    n_flagged <- sum(diag$flag)

    cat(sprintf(
      "  [%s]  converged: %-5s  factors: %d  flagged: %d\n",
      lbl,
      as.character(converged),
      nrow(diag),
      n_flagged
    ))
  }

  cat("\n")
  n_obs <- nrow(x$index_scores)
  cat(sprintf("Index scores: %d observation(s)  x  %d factor(s)\n",
              n_obs, ncol(x$index_scores)))

  if (!is.null(x$correction)) {
    cat(sprintf("Correction  : applied (%s scores)\n", x$method))
  }

  invisible(x)
}
