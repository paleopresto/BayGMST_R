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
#'     mean of the proxy-only fitted value, `y_ins_fitted`, which is
#'     estimated from the proxies of the whole instrumental period, not from
#'     the instrumental data; reconstructed years: posterior mean of the
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

#' Extract posterior ensemble members of the reconstruction
#'
#' Returns individual posterior draws of the temperature series, one column
#' per ensemble member, for analyses that need more than the summary bands of
#' [reconstruct()] (e.g. propagating uncertainty into trends, extremes, or
#' comparisons with other reconstructions). Members are evenly thinned from
#' all post-warmup draws, so each column is one coherent draw of the whole
#' series.
#'
#' For the reconstructed (pre-instrumental) years each member is a draw of the
#' latent temperature `y_mis`. For the instrumental years, `instrumental`
#' chooses what goes in:
#' \describe{
#'   \item{`"proxy_only"`}{(default) draws of `y_ins_fitted`, the estimate of
#'     the instrumental years from the proxies alone (the quantity whose bands
#'     [reconstruct()] reports for those years). Use this to compare the
#'     proxy-based estimate with the instrumental record.}
#'   \item{`"observed"`}{the instrumental values themselves, identical in
#'     every member. Because the model treats instrumental temperature as
#'     error-free, this is the posterior of the temperature in those years;
#'     use it for a continuous GMST ensemble spanning the whole period.}
#' }
#'
#' @param fit A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param n_members Number of ensemble members to return. It cannot exceed
#'   the number of post-warmup draws in the fit (`chains * iter_sampling` in
#'   [fit_baygmst()]; 6000 with the defaults), so for a larger ensemble, e.g.
#'   10000 members, fit with more draws (e.g. `chains = 4, iter_sampling =
#'   2500`). `NULL` returns all draws.
#' @param instrumental What to use for the instrumental years:
#'   `"proxy_only"` (default) or `"observed"`. See Details.
#'
#' @return A `data.frame` with one row per year in `fit$years` and columns
#'   \item{year}{Calendar year.}
#'   \item{type}{`"reconstruction"` or `"instrumental"`.}
#'   \item{T_obs}{Observed instrumental temperature (`NA` for reconstructed
#'     years).}
#'   \item{m0001, m0002, ...}{One column per ensemble member.}
#'   The draw indices used are stored in `attr(x, "draws")`. Save it yourself,
#'   e.g. `write.csv(x, gzfile("gmst_ensemble.csv.gz"), row.names = FALSE)`.
#'
#' @seealso [reconstruct()] for summary bands.
#' @export
reconstruction_ensemble <- function(fit,
                                    n_members = 1000,
                                    instrumental = c("proxy_only",
                                                     "observed")) {
  stopifnot(inherits(fit, "baygmst_fit"))
  instrumental <- match.arg(instrumental)

  idx_obs <- fit$idx_obs
  idx_mis <- fit$idx_mis
  mis_cols <- paste0("y_mis[", seq_along(idx_mis), "]")
  ins_cols <- paste0("y_ins_fitted[", seq_along(idx_obs), "]")
  vars <- c(if (length(idx_mis) > 0) "y_mis",
            if (instrumental == "proxy_only") "y_ins_fitted")

  dm <- if (length(vars) > 0) {
    posterior::as_draws_matrix(fit$fit$draws(vars))
  }
  n_draws <- if (!is.null(dm)) {
    posterior::ndraws(dm)
  } else {
    fit$fit$num_chains() * fit$fit$metadata()$iter_sampling
  }
  if (is.null(n_members)) n_members <- n_draws
  if (length(n_members) != 1 || is.na(n_members) || n_members < 1 ||
      n_members != round(n_members)) {
    stop("n_members must be a single positive whole number or NULL.",
         call. = FALSE)
  }
  if (n_members > n_draws) {
    stop(sprintf(paste0(
      "n_members (%d) exceeds the %d post-warmup draws in this fit. ",
      "Refit with more draws (chains * iter_sampling >= %d in fit_baygmst()), ",
      "or use n_members = NULL for all draws."),
      n_members, n_draws, n_members), call. = FALSE)
  }
  keep <- unique(round(seq(1, n_draws, length.out = n_members)))

  ens <- matrix(NA_real_, nrow = length(fit$years), ncol = length(keep))
  if (length(idx_mis) > 0) {
    ens[idx_mis, ] <- t(unclass(dm[keep, mis_cols, drop = FALSE]))
  }
  if (instrumental == "proxy_only") {
    ens[idx_obs, ] <- t(unclass(dm[keep, ins_cols, drop = FALSE]))
  } else {
    ens[idx_obs, ] <- fit$data$y_obs
  }
  colnames(ens) <- sprintf("m%04d", seq_along(keep))

  T_obs <- rep(NA_real_, length(fit$years))
  T_obs[idx_obs] <- fit$data$y_obs
  out <- data.frame(
    year  = fit$years,
    type  = ifelse(seq_along(fit$years) %in% idx_obs, "instrumental",
                   "reconstruction"),
    T_obs = T_obs,
    ens,
    check.names = FALSE
  )
  attr(out, "draws") <- keep
  out
}
