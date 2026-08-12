# Internal utilities. Not exported -- see man/ for the documented, exported
# API.

#' Validate a partition-year argument
#'
#' @param x A candidate year value.
#' @param name The argument name, used in error messages.
#' @return Invisibly `NULL`. Called for its side effect of raising an
#'   informative error if `x` is not a single, finite, positive integer.
#' @noRd
check_year <- function(x, name) {
  if (is.null(x)) {
    stop(sprintf("%s is missing.", name), call. = FALSE)
  }
  if (
    !is.numeric(x) ||
    length(x) != 1 ||
    is.na(x) ||
    !is.finite(x) ||
    x <= 0 ||
    x != as.integer(x)
  ) {
    stop(sprintf("%s must be a single positive integer.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
