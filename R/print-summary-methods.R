#' Print a fitted BayGMST model
#'
#' @param x A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param ... Unused; present for S3 consistency.
#' @return `x`, invisibly.
#' @export
print.baygmst_fit <- function(x, ...) {
  cat("<baygmst_fit>\n")
  cat(sprintf(
    "  Years:          %d-%d (%d total, %d instrumental, %d reconstructed)\n",
    min(x$years), max(x$years), length(x$years),
    length(x$idx_obs), length(x$idx_mis)
  ))
  cat(sprintf(
    "  Chains/draws:   %d chains, %d post-warmup draws each\n",
    x$fit$num_chains(), x$fit$metadata()$iter_sampling
  ))
  cat("  Call:           ")
  print(x$call)
  invisible(x)
}

#' Summarize a fitted BayGMST model
#'
#' Reports posterior summaries for the structural parameters (matching the
#' variable set used in the original \code{BayGMST_v1.0.R} script) together
#' with simple instrumental-period fit performance statistics (RMSE, MAE,
#' bias, correlation, and R-squared, both raw and linearly detrended).
#'
#' @param object A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param ... Unused; present for S3 consistency.
#'
#' @return An object of class `"summary.baygmst_fit"`, a list with elements
#'   `posterior` (a data.frame of parameter posterior summaries) and
#'   `performance` (a one-row data.frame of instrumental-period fit
#'   statistics).
#' @export
summary.baygmst_fit <- function(object, ...) {
  posterior_summary <- object$fit$summary(variables = c(
    "alpha0", "alpha1", "phi_R", "phi_T",
    "beta0", "betaG", "betaS", "betaV",
    "sigma_y", "sigma_z"
  ))

  recon <- reconstruct(object)
  ins   <- recon[recon$type == "instrumental", , drop = FALSE]
  err   <- ins$T_mean - ins$T_obs

  dt <- stats::na.omit(ins[, c("year", "T_obs", "T_mean")])
  r2_detrended <- if (nrow(dt) > 2) {
    res_obs  <- stats::residuals(stats::lm(T_obs  ~ year, data = dt))
    res_mean <- stats::residuals(stats::lm(T_mean ~ year, data = dt))
    stats::cor(res_obs, res_mean)^2
  } else {
    NA_real_
  }

  performance <- data.frame(
    RMSE            = sqrt(mean(err^2, na.rm = TRUE)),
    MAE             = mean(abs(err), na.rm = TRUE),
    Bias            = mean(err, na.rm = TRUE),
    Correlation     = stats::cor(ins$T_obs, ins$T_mean, use = "complete.obs"),
    R2              = stats::cor(ins$T_obs, ins$T_mean, use = "complete.obs")^2,
    R2_detrended    = r2_detrended
  )

  structure(
    list(posterior = posterior_summary, performance = performance),
    class = "summary.baygmst_fit"
  )
}

#' @export
print.summary.baygmst_fit <- function(x, ...) {
  cat("Posterior parameter summary:\n")
  print(as.data.frame(x$posterior))
  cat("\nInstrumental-period fit performance:\n")
  print(x$performance, row.names = FALSE)
  invisible(x)
}
