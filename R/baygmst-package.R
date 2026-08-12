#' BayGMST: Bayesian Reconstruction of Global Mean Surface Temperature
#'
#' Fits a Bayesian hierarchical AR(1) state-space model relating global mean
#' surface temperature (GMST) to a reduced paleoclimate proxy network,
#' instrumental temperature observations, and radiative forcings
#' (greenhouse gases, volcanic aerosols, and solar irradiance). See
#' \code{vignette("baygmst-intro", package = "BayGMST")} for a worked
#' end-to-end example.
#'
#' @section Model fitting requires CmdStan:
#' [fit_baygmst()] and [cv_baygmst()] compile and sample a 'Stan' model via
#' the \pkg{instantiate} and \pkg{cmdstanr} packages, which in turn require a
#' working CmdStan installation on the user's machine. Use
#' \code{instantiate::stan_cmdstan_exists()} to check whether one is
#' available, and see \code{vignette("cmdstanr", package = "cmdstanr")} for
#' setup instructions. [reduce_proxies()] and [transform_forcings()] do not
#' require CmdStan.
#'
#' @section Provenance:
#' This package was restructured from a collection of top-level analysis
#' scripts; see \code{NEWS.md} for exactly what changed (fixed bugs, dropped
#' dead code, unresolved issues flagged but left unchanged) and
#' \code{system.file("legacy-scripts", package = "BayGMST")} for the
#' original, unmodified scripts.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom instantiate stan_package_model
#' @importFrom ggplot2 .data
## usethis namespace: end
NULL
