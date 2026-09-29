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

#' Plot a fitted BayGMST model as a single summary figure
#'
#' Reproduces the combined figure of the original \code{BayGMST_v1.0.R}
#' script (and the PReSto manuscript): the reconstruction from
#' [plot_reconstruction()] on top, and below it three posterior-density
#' panels grouping the proxy coefficient (`alpha1`), the forcing
#' sensitivities (`betaG`, `betaV`, `betaS`), and the autoregressive
#' coefficients (`phi_R`, `phi_T`). The individual pieces remain available
#' through [plot_reconstruction()] and [plot_posterior_densities()].
#'
#' Requires the 'patchwork' package.
#'
#' @param x A `"baygmst_fit"` object, as returned by [fit_baygmst()].
#' @param title Figure title. Default
#'   `"GMST Reconstruction using a Reduced Proxy"`.
#' @param subtitle Figure subtitle. By default, states the instrumental
#'   period and the AR(1) model structure; pass `NULL` to omit it.
#' @param ... Unused; for compatibility with the [plot()] generic.
#'
#' @return A 'patchwork' object (which is also a `ggplot`); save it with
#'   [ggplot2::ggsave()] if a file is wanted.
#' @seealso [plot_reconstruction()], [plot_posterior_densities()],
#'   [plot_trace()]
#' @export
plot.baygmst_fit <- function(x,
                             title = "GMST Reconstruction using a Reduced Proxy",
                             subtitle = NULL,
                             ...) {
  stopifnot(inherits(x, "baygmst_fit"))
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("plot() on a baygmst_fit requires the 'patchwork' package. ",
         "Install it, or use plot_reconstruction() and ",
         "plot_posterior_densities() separately.", call. = FALSE)
  }
  if (missing(subtitle)) {
    obs_years <- range(x$years[x$idx_obs])
    subtitle <- sprintf("Instrumental period: (%s, %s);  AR(1) in T and R equations",
                        obs_years[1], obs_years[2])
  }

  groups <- list(
    alpha = "alpha1",
    beta  = c("betaG", "betaV", "betaS"),
    phi   = c("phi_R", "phi_T")
  )
  draws_df <- as.data.frame(
    x$fit$draws(variables = unlist(groups, use.names = FALSE), format = "df")
  )
  df_hist <- tidyr::pivot_longer(
    draws_df[, unlist(groups, use.names = FALSE), drop = FALSE],
    cols = dplyr::everything(),
    names_to = "parameter", values_to = "value"
  )

  p_ts <- plot_reconstruction(x, title = NULL)
  p_alpha <- posterior_density_panel(
    df_hist, groups$alpha,
    xlab = expression("Signed RP-T Coefficient (" * degree * C^{-1} * ")"),
    ylab = "Post. Density"
  )
  p_beta <- posterior_density_panel(
    df_hist, groups$beta,
    xlab = expression("Forcing Sensitivity (" * degree * C ~ m^2 ~ W^{-1} * ")")
  )
  p_phi <- posterior_density_panel(
    df_hist, groups$phi,
    xlab = expression("Autoregressive Components" * phantom(""^{-1})),
    xlim = c(-0.1, 1)
  )

  p_hist <- patchwork::wrap_plots(p_alpha, p_beta, p_phi, nrow = 1)
  patchwork::wrap_plots(p_ts, p_hist, ncol = 1, heights = c(5, 1)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
    )
}

# One posterior-density panel of plot.baygmst_fit(), styled as in the
# original BayGMST_v1.0.R script (colors, linetypes, plotmath labels).
posterior_density_panel <- function(df_hist, parameters, xlab, ylab = "",
                                    xlim = NULL) {
  cols <- c(alpha1 = "#5899E2", betaG = "#07D664", betaV = "#1CCAD8",
            betaS = "#7B287D", phi_R = "#625834", phi_T = "#FA9500")
  labs <- c(alpha1 = expression(alpha[T]), betaG = expression(beta[G]),
            betaV = expression(beta[V]), betaS = expression(beta[S]),
            phi_R = expression(phi[R]), phi_T = expression(phi[T]))
  ltys <- c(alpha1 = "solid", betaG = "solid", betaV = "42", betaS = "11",
            phi_R = "solid", phi_T = "42")

  dat <- df_hist[df_hist$parameter %in% parameters, , drop = FALSE]
  if (is.null(xlim)) xlim <- range(dat$value, na.rm = TRUE)

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = .data$value, fill = .data$parameter,
                 color = .data$parameter, linetype = .data$parameter)
  ) +
    ggplot2::geom_density(linewidth = 0.65) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_fill_manual(values = ggplot2::alpha(cols[parameters], 0.35),
                               breaks = parameters, labels = labs[parameters]) +
    ggplot2::scale_color_manual(values = cols[parameters],
                                breaks = parameters, labels = labs[parameters]) +
    ggplot2::scale_linetype_manual(values = ltys[parameters],
                                   breaks = parameters, labels = labs[parameters]) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = xlab, y = ylab, fill = NULL, color = NULL, linetype = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.35, "lines"),
      legend.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank()
    )
}
