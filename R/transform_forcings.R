#' Transform raw radiative forcings for use with \code{fit_baygmst()}
#'
#' Applies the forcing transformations used throughout the BayGMST pipeline:
#' volcanic forcing is mapped to a negative saturating form, greenhouse-gas
#' (CO2) forcing is log-transformed relative to a reference concentration,
#' and solar forcing is centered on its own mean. Ported from the
#' normalization block near the top of the original \code{BayGMST_v1.0.R}
#' script (see \code{system.file("legacy-scripts", package = "BayGMST")}).
#'
#' @param G Numeric vector of raw greenhouse-gas forcing (e.g. atmospheric
#'   CO2 concentration in ppm).
#' @param V Numeric vector of raw volcanic forcing (e.g. aerosol optical
#'   depth). Same length as `G`.
#' @param S Numeric vector of raw solar forcing (e.g. total solar
#'   irradiance). Same length as `G`.
#' @param co2_c0 Reference/baseline CO2 concentration used in the log
#'   transform. Default `280` (a common preindustrial reference, ppm).
#' @param co2_coef Multiplier applied to the log CO2 ratio. Default `5.35`.
#' @param vol_coef Multiplier applied to the volcanic saturating transform.
#'   Default `25`.
#'
#' @details
#' The transforms applied are:
#' \deqn{V' = -|\mathrm{vol\_coef}| \times (1 - e^{-V})}
#' \deqn{G' = \mathrm{co2\_coef} \times \log(G / \mathrm{co2\_c0})}
#' \deqn{S' = S - \bar{S}}
#'
#' @return A list with three numeric vectors, `G`, `V`, and `S` -- the
#'   transformed forcings, each the same length as the corresponding input.
#'
#' @examples
#' fc <- transform_forcings(
#'   G = c(280, 300, 340, 400),
#'   V = c(0, 0.05, 0.02, 0.10),
#'   S = c(1360, 1361, 1360.5, 1361.2)
#' )
#' fc$G
#'
#' @export
transform_forcings <- function(G, V, S,
                                co2_c0 = 280,
                                co2_coef = 5.35,
                                vol_coef = 25) {
  stopifnot(
    is.numeric(G), is.numeric(V), is.numeric(S),
    length(G) == length(V), length(V) == length(S)
  )
  list(
    G = co2_coef * log(G / co2_c0),
    V = -abs(vol_coef) * (1 - exp(-V)),
    S = S - mean(S, na.rm = TRUE)
  )
}
