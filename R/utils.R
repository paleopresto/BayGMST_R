# Internal utilities. Not exported.

#' Stop early, with installation instructions, if CmdStan is unavailable
#'
#' Called by [fit_baygmst()] and [cv_baygmst()] before touching the Stan
#' model, so users get an actionable message rather than a low-level
#' compilation error.
#' @noRd
check_cmdstan <- function(fun) {
  if (!instantiate::stan_cmdstan_exists()) {
    stop(
      fun, "() requires a working CmdStan installation. Install one with ",
      "cmdstanr::install_cmdstan(), then check ",
      "instantiate::stan_cmdstan_exists(). See ?BayGMST for details.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
