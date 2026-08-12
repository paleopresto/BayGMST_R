#' Extract a tidy reconstruction from a fitted BayGMST model
#'
#' Builds a tidy `data.frame` combining the observed instrumental
#' temperatures with the posterior reconstruction of the pre-instrumental
#' years, equivalent to the `gmst_reconstruction_data.csv` produced as a
#' side effect by the original \code{BayGMST_v1.0.R} script -- but returned
#' as an object rather than written to disk. Save it yourself
#' (\code{write.csv(reconstruct(fit), "path/to/file.csv")}) if you want a
#' file.
#'
#' @param fit A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param probs_inner,probs_outer Numeric vectors of length 2, the lower/
#'   upper quantile levels for the "inner" and "outer" credible bands.
#'   Defaults `c(0.16, 0.84)` (68%) and `c(0.025, 0.975)` (95%), matching
#'   the original script.
#'
#' @return A `data.frame` with one row per year in
#'   `fit$years`, and columns:
#'   \item{year}{Calendar year.}
#'   \item{type}{`"instrumental"` or `"reconstruction"`.}
#'   \item{T_obs}{Observed instrumental temperature (`NA` for reconstructed
#'     years).}
#'   \item{T_mean}{Posterior mean (instrumental years: posterior predictive
#'     mean of the fitted value; reconstructed years: posterior mean of the
#'     latent temperature).}
#'   \item{T_lo_inner, T_hi_inner}{Inner credible band (`probs_inner`).}
#'   \item{T_lo_outer, T_hi_outer}{Outer credible band (`probs_outer`).}
#'
#' @seealso [plot_reconstruction()]
#' @export
reconstruct <- function(fit,
                         probs_inner = c(0.16, 0.84),
                         probs_outer = c(0.025, 0.975)) {
  stopifnot(inherits(fit, "baygmst_fit"))

  years   <- fit$years
  idx_obs <- fit$idx_obs
  idx_mis <- fit$idx_mis

  ins_summary <- fit$fit$summary(
    variables = "y_ins_fitted",
    ~ stats::setNames(
      as.list(posterior::quantile2(.x, probs = c(probs_inner, probs_outer))),
      c("lo_inner", "hi_inner", "lo_outer", "hi_outer")
    ),
    "mean"
  )
  ins <- data.frame(
    year       = years[idx_obs],
    type       = "instrumental",
    T_obs      = fit$data$y_obs,
    T_mean     = ins_summary$mean,
    T_lo_inner = ins_summary$lo_inner,
    T_hi_inner = ins_summary$hi_inner,
    T_lo_outer = ins_summary$lo_outer,
    T_hi_outer = ins_summary$hi_outer
  )

  if (length(idx_mis) > 0) {
    draws     <- fit$fit$draws("y_mis")
    idx_names <- paste0("y_mis[", seq_along(idx_mis), "]")
    mat       <- posterior::as_draws_matrix(draws)[, idx_names, drop = FALSE]

    recon <- data.frame(
      year       = years[idx_mis],
      type       = "reconstruction",
      T_obs      = NA_real_,
      T_mean     = apply(mat, 2, mean),
      T_lo_inner = apply(mat, 2, stats::quantile, probs_inner[1]),
      T_hi_inner = apply(mat, 2, stats::quantile, probs_inner[2]),
      T_lo_outer = apply(mat, 2, stats::quantile, probs_outer[1]),
      T_hi_outer = apply(mat, 2, stats::quantile, probs_outer[2])
    )
    out <- rbind(recon, ins)
  } else {
    out <- ins
  }

  out[order(out$year), , drop = FALSE]
}
