# Internal utilities. Not exported.

# Installation steps shown by the errors below. BayGMST compiles its Stan
# models when it is installed, and only if cmdstanr and CmdStan are present
# at that moment, so installing them afterwards means reinstalling BayGMST.
cmdstan_install_steps <- function() {
  paste0(
    "To enable model fitting:\n",
    "  1. install.packages(\"cmdstanr\", repos = c(",
    "\"https://stan-dev.r-universe.dev\", getOption(\"repos\")))\n",
    "  2. cmdstanr::install_cmdstan()\n",
    "  3. Reinstall BayGMST, so its Stan models are compiled.\n",
    "See ?BayGMST for details."
  )
}

#' Stop early, with installation instructions, if CmdStan is unavailable
#'
#' Called by [fit_baygmst()] and [cv_baygmst()] before touching the Stan
#' model, so users get an actionable message rather than a low-level
#' compilation error.
#' @noRd
check_cmdstan <- function(fun) {
  if (!instantiate::stan_cmdstan_exists()) {
    stop(
      fun, "() requires the 'cmdstanr' package and a working CmdStan ",
      "installation.\n", cmdstan_install_steps(),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Load a packaged Stan model, checking that it was compiled at install time
#'
#' `instantiate` compiles the models when BayGMST is installed, and silently
#' skips compilation if cmdstanr or CmdStan is missing then. Without this
#' check, a user who installs CmdStan afterwards would pass check_cmdstan()
#' and then hit cmdstanr's "Model not compiled" error from `$sample()`.
#' @noRd
baygmst_model <- function(name, fun) {
  check_cmdstan(fun)
  # A cmdstanr CmdStanModel object; cmdstanr is in Suggests (it is not on
  # CRAN), and check_cmdstan() has confirmed it is installed.
  model <- instantiate::stan_package_model(name = name, package = "BayGMST")
  if (!file.exists(model$exe_file())) {
    stop(
      "BayGMST's Stan model '", name, "' was not compiled when BayGMST was ",
      "installed, most likely because cmdstanr or CmdStan was not ",
      "available at the time.\n", cmdstan_install_steps(),
      call. = FALSE
    )
  }
  model
}
