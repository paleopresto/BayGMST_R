#' BayGMST: Bayesian Reconstruction of Global Mean Surface Temperature
#'
#' Fits a Bayesian hierarchical AR(1) state-space model relating global mean
#' surface temperature (GMST) to a reduced paleoclimate proxy network,
#' instrumental temperature observations, and radiative forcings
#' (greenhouse gases, volcanic aerosols, and solar irradiance). See
#' \code{vignette("baygmst-intro", package = "BayGMST")} for a worked
#' end-to-end example, and
#' \code{vignette("baygmst-model", package = "BayGMST")} for the full
#' statistical model description, which builds on the reduced-proxy
#' Bayesian framework of Barboza et al. (2014, 2019) and Wang (2020).
#'
#' @section Model fitting requires CmdStan:
#' [fit_baygmst()] and [cv_baygmst()] compile and sample a 'Stan' model via
#' the \pkg{instantiate} and \pkg{cmdstanr} packages, which in turn require a
#' working CmdStan installation on the user's machine. \pkg{cmdstanr} is not
#' on CRAN; install it from \url{https://stan-dev.r-universe.dev}, then run
#' \code{cmdstanr::install_cmdstan()}. BayGMST compiles its Stan models when
#' it is installed, and only if \pkg{cmdstanr} and CmdStan are available at
#' that moment, so install them first, or reinstall BayGMST after adding
#' them. Use \code{instantiate::stan_cmdstan_exists()} to check whether
#' CmdStan is available. [reduce_proxies()] and [transform_forcings()] do
#' not require CmdStan.
#'
#' @section Provenance:
#' This package was restructured from a collection of top-level analysis
#' scripts; see \code{NEWS.md} for exactly what changed (fixed bugs, dropped
#' dead code, and modeling decisions reviewed by the authors) and
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
