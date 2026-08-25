#' Plot a BayGMST reconstruction
#'
#' Ported from the `p_ts` plotting block of the original
#' \code{BayGMST_v1.0.R} script. Unlike the original script, this returns a
#' `ggplot` object instead of calling \code{ggsave()} as a side effect --
#' save it yourself with \code{ggplot2::ggsave()} if you want a file.
#'
#' @param fit A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param title Plot title. Default `"GMST Reconstruction"`.
#'
#' @return A `ggplot` object.
#' @export
plot_reconstruction <- function(fit, title = "GMST Reconstruction") {
  stopifnot(inherits(fit, "baygmst_fit"))
  recon <- reconstruct(fit)

  ggplot2::ggplot(recon) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = .data$year, ymin = .data$T_lo_outer, ymax = .data$T_hi_outer),
      fill = "cyan3", alpha = 0.35
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = .data$year, ymin = .data$T_lo_inner, ymax = .data$T_hi_inner),
      fill = "cyan3", alpha = 0.75
    ) +
    ggplot2::geom_line(
      data = recon[recon$type == "instrumental", ],
      ggplot2::aes(x = .data$year, y = .data$T_obs),
      color = "orange", linewidth = 0.4, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$year, y = .data$T_mean),
      color = "darkorchid4", linewidth = 0.55
    ) +
    ggplot2::labs(x = "Year CE", y = "GMST Anomaly (\u00b0C)", title = title) +
    ggplot2::theme_light(base_size = 11)
}

#' Trace plots for a fitted BayGMST model
#'
#' Ported from the trace-plot block of the original \code{BayGMST_v1.0.R}
#' script.
#'
#' @param fit A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param parameters Character vector of Stan parameter names to plot.
#'   Default `c("alpha1", "betaG", "betaV", "betaS", "phi_R", "phi_T")`.
#'
#' @return A `ggplot` object, faceted by parameter.
#' @export
plot_trace <- function(fit,
                        parameters = c("alpha1", "betaG", "betaV", "betaS", "phi_R", "phi_T")) {
  stopifnot(inherits(fit, "baygmst_fit"))
  # plain data.frame: subsetting a posterior 'draws_df' emits a
  # "Dropping 'draws_df' class" warning downstream
  draws_df <- as.data.frame(fit$fit$draws(variables = parameters, format = "df"))

  trace_df <- tidyr::pivot_longer(
    draws_df,
    cols = dplyr::all_of(parameters),
    names_to = "parameter", values_to = "value"
  )

  ggplot2::ggplot(
    trace_df,
    ggplot2::aes(x = .data$.iteration, y = .data$value,
                 group = .data$.chain, color = factor(.data$.chain))
  ) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.3) +
    ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Iteration", y = "Draw value", color = "Chain", title = "Trace plots") +
    ggplot2::theme_bw()
}

#' Posterior density plots for a fitted BayGMST model
#'
#' Ported from the posterior-histogram block of the original
#' \code{BayGMST_v1.0.R} script.
#'
#' @inheritParams plot_trace
#'
#' @return A `ggplot` object overlaying the posterior densities of
#'   `parameters`.
#' @export
plot_posterior_densities <- function(fit,
                                      parameters = c("alpha1", "betaG", "betaV", "betaS", "phi_R", "phi_T")) {
  stopifnot(inherits(fit, "baygmst_fit"))
  draws_df <- as.data.frame(fit$fit$draws(variables = parameters, format = "df"))

  df_hist <- tidyr::pivot_longer(
    draws_df[, parameters, drop = FALSE],
    cols = dplyr::everything(),
    names_to = "parameter", values_to = "value"
  )

  ggplot2::ggplot(
    df_hist,
    ggplot2::aes(x = .data$value, fill = .data$parameter, color = .data$parameter)
  ) +
    ggplot2::geom_density(linewidth = 0.65, alpha = 0.35) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::labs(x = "Posterior dist.", y = "Post. density", fill = NULL, color = NULL) +
    ggplot2::theme_minimal(base_size = 10)
}

#' Plot k-fold cross-validation reconstructions
#'
#' Ported from the `p_ts_cv` plotting block of the original \code{cv_v0.1.R}
#' script.
#'
#' @param cv A `"baygmst_cv"` object, as returned by [cv_baygmst()].
#'
#' @return A `ggplot` object showing each fold's held-out reconstruction
#'   against its true instrumental value, faceted by fold.
#' @export
plot_cv <- function(cv) {
  stopifnot(inherits(cv, "baygmst_cv"))
  folds <- cv$folds

  ggplot2::ggplot(folds, ggplot2::aes(x = .data$year)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$T_lo_outer, ymax = .data$T_hi_outer),
      fill = "gray30", alpha = 0.25
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$T_lo_inner, ymax = .data$T_hi_inner),
      fill = "gray30", alpha = 0.5
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$T_true), color = "black", linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = .data$T_mean), color = "firebrick", linewidth = 0.6) +
    ggplot2::facet_wrap(~fold, scales = "free_x") +
    ggplot2::labs(
      x = "Year", y = "GMST Anomaly (\u00b0C)",
      title = "GMST Reconstructions via k-fold cross-validation",
      subtitle = sprintf("Median R2 = %.2f", stats::median(cv$r2))
    ) +
    ggplot2::theme_light(base_size = 10)
}
