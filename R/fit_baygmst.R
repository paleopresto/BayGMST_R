#' Fit the BayGMST Bayesian hierarchical AR(1) reconstruction model
#'
#' Prepares the data list expected by the \code{baygmst} 'Stan' model and
#' samples from it with \pkg{cmdstanr}, via \pkg{instantiate}. Ported from
#' the data-preparation and \code{mod$sample()} block of the original
#' \code{BayGMST_v1.0.R} script (see
#' \code{system.file("legacy-scripts", package = "BayGMST")}), restructured
#' to take its inputs as explicit arguments (rather than reading
#' \code{config.yml} and hardcoded relative paths) and to return an object
#' rather than leaving variables in the global environment.
#'
#' @param proxy The representative/reduced proxy series. Either a bare
#'   numeric vector with one value per year in `years` (aligned by
#'   position, no missing values), or any object accepted by
#'   [as_baygmst_proxy()] -- the output of [reduce_proxies()], a
#'   `data.frame` with year and value columns, or a composite object from
#'   another package (e.g. a 'compositeR' `paleoComposite`) -- in which
#'   case it is aligned to `years` by calendar year, with an informative
#'   error if any requested year is not covered.
#' @param instrumental_T Numeric vector, instrumental temperature, one value
#'   per year in `years`. `NA` for years without instrumental coverage (the
#'   pre-instrumental years the model reconstructs).
#' @param forcing_G,forcing_V,forcing_S Numeric vectors of *already
#'   transformed* greenhouse-gas, volcanic, and solar forcing -- see
#'   [transform_forcings()] -- one value per year in `years`. No missing
#'   values are allowed.
#' @param years Integer vector of calendar years, defining the row order for
#'   all of the above.
#' @param chains,parallel_chains,iter_warmup,iter_sampling,seed Passed to
#'   \code{cmdstanr}'s \code{$sample()} method. `iter_sampling` must be at
#'   least 1000 (the original script enforced the same minimum).
#' @param ... Additional arguments passed on to \code{$sample()}.
#'
#' @return An object of class `"baygmst_fit"`, a list with elements:
#'   \item{fit}{The underlying `CmdStanMCMC` object (see
#'     \code{cmdstanr::CmdStanMCMC}).}
#'   \item{years}{The `years` argument, echoed back.}
#'   \item{idx_obs, idx_mis}{Integer indices (into `years`) of the
#'     instrumental and pre-instrumental years, respectively.}
#'   \item{data}{The data list passed to Stan.}
#'   \item{cmdstan_version}{The CmdStan version used, for provenance.}
#'   \item{call}{The matched call.}
#'   Use [reconstruct()] to extract a tidy reconstruction `data.frame`,
#'   [plot_reconstruction()] / [plot_trace()] / [plot_posterior_densities()]
#'   to visualize it, and \code{summary()} for posterior parameter summaries.
#'
#' @section CmdStan required:
#' This function requires a working CmdStan installation (via \pkg{cmdstanr}
#' / \pkg{instantiate}). Check availability first with
#' \code{instantiate::stan_cmdstan_exists()}.
#'
#' @examples
#' \donttest{
#' if (instantiate::stan_cmdstan_exists()) {
#'   set.seed(1)
#'   years <- 1:200
#'   forcing <- transform_forcings(
#'     G = 280 + cumsum(rgamma(200, 0.05, 1)),
#'     V = abs(rnorm(200, 0, 0.05)),
#'     S = 1361 + rnorm(200, 0, 0.3)
#'   )
#'   true_T <- cumsum(rnorm(200, 0, 0.05))
#'   instrumental_T <- c(rep(NA_real_, 150), true_T[151:200])
#'   proxy <- true_T + rnorm(200, 0, 0.1)
#'
#'   fit <- fit_baygmst(
#'     proxy = proxy,
#'     instrumental_T = instrumental_T,
#'     forcing_G = forcing$G,
#'     forcing_V = forcing$V,
#'     forcing_S = forcing$S,
#'     years = years,
#'     chains = 1,
#'     iter_warmup = 200,
#'     iter_sampling = 1000
#'   )
#'   summary(fit)
#' }
#' }
#'
#' @export
fit_baygmst <- function(proxy,
                         instrumental_T,
                         forcing_G,
                         forcing_V,
                         forcing_S,
                         years,
                         chains = 4,
                         parallel_chains = 1,
                         iter_warmup = 500,
                         iter_sampling = 1500,
                         seed = NULL,
                         ...) {
  NT <- length(years)
  proxy <- resolve_proxy(proxy, years)
  check_model_inputs(proxy, instrumental_T, forcing_G, forcing_V, forcing_S,
                     NT)
  if (iter_sampling < 1000) {
    stop("iter_sampling must be at least 1000.", call. = FALSE)
  }

  idx_obs <- which(!is.na(instrumental_T))
  idx_mis <- which(is.na(instrumental_T))
  if (length(idx_obs) == 0) {
    stop("instrumental_T has no observed (non-NA) values.", call. = FALSE)
  }
  check_cmdstan("fit_baygmst")

  data_list <- list(
    NT      = NT,
    NT_obs  = length(idx_obs),
    NT_mis  = length(idx_mis),
    idx_obs = as.integer(idx_obs),
    idx_mis = as.integer(idx_mis),
    G       = as.vector(forcing_G),
    S       = as.vector(forcing_S),
    V       = as.vector(forcing_V),
    y_obs   = as.vector(instrumental_T[idx_obs]),
    z       = as.vector(proxy)
  )

  # `model` is a cmdstanr CmdStanModel object; its $sample() method is
  # provided by cmdstanr, not instantiate, which is why cmdstanr must stay a
  # real Imports (not just Suggests) even though it's never called via
  # `cmdstanr::` here -- see the instantiate package's own packaging guidance.
  model <- instantiate::stan_package_model(name = "baygmst", package = "BayGMST")
  fit <- model$sample(
    data             = data_list,
    chains           = chains,
    parallel_chains  = parallel_chains,
    iter_warmup      = iter_warmup,
    iter_sampling    = iter_sampling,
    seed             = seed,
    ...
  )

  structure(
    list(
      fit             = fit,
      years           = years,
      idx_obs         = idx_obs,
      idx_mis         = idx_mis,
      data            = data_list,
      cmdstan_version = as.character(cmdstanr::cmdstan_version()),
      call            = match.call()
    ),
    class = "baygmst_fit"
  )
}
