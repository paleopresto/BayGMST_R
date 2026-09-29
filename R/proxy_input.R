#' Coerce a reduced-proxy series into BayGMST's input format
#'
#' Normalizes the many shapes a reduced/composited proxy series can arrive
#' in -- the output of [reduce_proxies()], a plain `data.frame`, a bare
#' vector, a time-by-ensemble matrix, or a composite object from another
#' package (such as a `paleoComposite` from the in-development 'compositeR'
#' package) -- into a single, year-aware `"baygmst_proxy"` object that
#' [fit_baygmst()] and [cv_baygmst()] accept directly and align by calendar
#' year (rather than by position).
#'
#' @param x The object to coerce. Supported inputs:
#'   \describe{
#'     \item{`baygmst_rp`}{The output of [reduce_proxies()]; uses its
#'       `$composite` series.}
#'     \item{`data.frame`}{Must contain a year column (named `year` or
#'       `Year`) and exactly one other numeric column (e.g. the
#'       `example_proxy.csv` shipped in `inst/extdata`, with columns
#'       `Year, RP1`).}
#'     \item{`numeric` vector}{One value per element of `years`, which must
#'       also be supplied.}
#'     \item{`matrix` (or `data.frame` of ensemble columns)}{Rows are time
#'       steps (matching `years`), columns are ensemble members; collapsed
#'       to a single series with `collapse`.}
#'     \item{composite `list`}{Any list with elements `$ages` and
#'       `$composite` (rows = time steps of `$ages`, columns = ensemble
#'       members), the shape produced by compositing packages in the LiPD
#'       ecosystem (e.g. `compositeR::compositeEnsembles2()`). Ages are
#'       assumed to be in years BP (before present, i.e. before 1950 CE)
#'       unless `age_units = "CE"`. The ensemble is collapsed with
#'       `collapse`.}
#'   }
#' @param years Numeric vector of calendar years (CE). Required when `x`
#'   carries no time axis of its own (bare vectors and matrices); ignored
#'   otherwise.
#' @param age_units For composite lists only: `"BP"` (default; ages are
#'   years before 1950 CE and are converted via `year = 1950 - age`) or
#'   `"CE"` (ages are already calendar years).
#' @param collapse How to collapse an ensemble (matrix or composite list)
#'   to a single series: `"median"` (default) or `"mean"`, taken across
#'   ensemble members within each time step, ignoring `NA`s. Time steps
#'   whose collapsed value is not finite (e.g. bins no record covers) are
#'   dropped with a message.
#' @param ... Passed between methods.
#'
#' @details
#' # Scale expectations
#' The BayGMST 'Stan' model regresses the proxy series on temperature with
#' an estimated intercept and slope, so the proxy does not need to be
#' pre-calibrated to temperature units. The (weakly informative) priors do,
#' however, assume the series is on roughly unit scale: a standardized
#' (z-score) composite or a temperature-anomaly-scale series is
#' appropriate; a raw-unit series (per mil, mm, ...) with large numerical
#' magnitude should be standardized first.
#'
#' # Uncertainty
#' In this version the ensemble dimension is collapsed before fitting, so
#' compositing/ensemble uncertainty is not propagated into the posterior.
#' A future release may add an ensemble-aware fitting pathway; the
#' `"baygmst_proxy"` object records `n_ensemble` so existing code will not
#' need to change shape when it does.
#'
#' @return An object of class `"baygmst_proxy"`: a `data.frame` with
#'   columns `year` (ascending) and `proxy`, plus attributes `source` (a
#'   short label for where the series came from) and `n_ensemble` (the
#'   number of ensemble members collapsed; `1L` for deterministic inputs).
#'
#' @examples
#' # from a bare vector
#' p <- as_baygmst_proxy(rnorm(100), years = 1901:2000)
#'
#' # from a data.frame like inst/extdata/example_proxy.csv
#' df <- data.frame(Year = 1901:2000, RP1 = rnorm(100))
#' p <- as_baygmst_proxy(df)
#'
#' # from a composite ensemble with ages in yr BP (compositeR-style)
#' comp <- list(
#'   ages = seq(950, 0, by = -10),               # 1000-1950 CE
#'   composite = matrix(rnorm(96 * 20), 96, 20)  # 20 ensemble members
#' )
#' p <- as_baygmst_proxy(comp)
#' head(p)
#'
#' @seealso [fit_baygmst()], [reduce_proxies()]
#' @export
as_baygmst_proxy <- function(x, ...) {
  UseMethod("as_baygmst_proxy")
}

#' @rdname as_baygmst_proxy
#' @export
as_baygmst_proxy.baygmst_proxy <- function(x, ...) {
  x
}

#' @rdname as_baygmst_proxy
#' @export
as_baygmst_proxy.baygmst_rp <- function(x, ...) {
  new_baygmst_proxy(
    year = x$composite$year,
    proxy = x$composite$RP1,
    source = sprintf("reduce_proxies(method = \"%s\")", x$method),
    n_ensemble = 1L
  )
}

#' @rdname as_baygmst_proxy
#' @export
as_baygmst_proxy.numeric <- function(x, years, ...) {
  if (missing(years)) {
    stop(
      "`years` is required when coercing a bare numeric vector: ",
      "as_baygmst_proxy(x, years = ...).",
      call. = FALSE
    )
  }
  if (length(years) != length(x)) {
    stop(
      sprintf(
        "length(years) (%d) must match length(x) (%d).",
        length(years), length(x)
      ),
      call. = FALSE
    )
  }
  new_baygmst_proxy(year = years, proxy = x, source = "numeric vector",
                    n_ensemble = 1L)
}

#' @rdname as_baygmst_proxy
#' @export
as_baygmst_proxy.matrix <- function(x, years,
                                    collapse = c("median", "mean"), ...) {
  collapse <- match.arg(collapse)
  if (missing(years)) {
    stop(
      "`years` is required when coercing an ensemble matrix: ",
      "as_baygmst_proxy(x, years = ...).",
      call. = FALSE
    )
  }
  if (length(years) != nrow(x)) {
    stop(
      sprintf(
        "length(years) (%d) must match nrow(x) (%d) -- rows are time ",
        length(years), nrow(x)
      ),
      "steps, columns are ensemble members.",
      call. = FALSE
    )
  }
  new_baygmst_proxy(
    year = years,
    proxy = collapse_ensemble(x, collapse),
    source = sprintf("ensemble matrix (%d members, %s)", ncol(x), collapse),
    n_ensemble = ncol(x)
  )
}

#' @rdname as_baygmst_proxy
#' @export
as_baygmst_proxy.data.frame <- function(x, ...) {
  year_col <- intersect(c("year", "Year"), names(x))
  if (length(year_col) == 0) {
    stop(
      "data.frame input must have a `year` (or `Year`) column; found: ",
      paste(names(x), collapse = ", "), ".",
      call. = FALSE
    )
  }
  year_col <- year_col[1]
  value_cols <- setdiff(names(x), year_col)
  if (length(value_cols) != 1) {
    stop(
      "data.frame input must have exactly one value column besides `",
      year_col, "`; found ", length(value_cols), " (",
      paste(value_cols, collapse = ", "), "). For an ensemble, pass a ",
      "matrix (plus `years`) or a composite list instead.",
      call. = FALSE
    )
  }
  new_baygmst_proxy(
    year = x[[year_col]],
    proxy = x[[value_cols]],
    source = sprintf("data.frame column `%s`", value_cols),
    n_ensemble = 1L
  )
}

#' @rdname as_baygmst_proxy
#' @export
as_baygmst_proxy.list <- function(x, age_units = c("BP", "CE"),
                                  collapse = c("median", "mean"), ...) {
  age_units <- match.arg(age_units)
  collapse <- match.arg(collapse)
  if (is.data.frame(x$composite) &&
      all(c("year", "RP1") %in% names(x$composite))) {
    # an (unclassed) reduce_proxies()-shaped list
    return(new_baygmst_proxy(
      year = x$composite$year,
      proxy = x$composite$RP1,
      source = "reduce_proxies()-style list",
      n_ensemble = 1L
    ))
  }
  if (is.null(x$ages) || is.null(x$composite)) {
    stop(
      "list input must be a composite object with `$ages` and `$composite` ",
      "elements (rows of `$composite` = time steps of `$ages`, columns = ",
      "ensemble members), e.g. the output of ",
      "compositeR::compositeEnsembles2().",
      call. = FALSE
    )
  }
  comp <- as.matrix(as.data.frame(x$composite))
  if (nrow(comp) != length(x$ages)) {
    stop(
      sprintf(
        "nrow(x$composite) (%d) must match length(x$ages) (%d).",
        nrow(comp), length(x$ages)
      ),
      call. = FALSE
    )
  }
  year <- if (age_units == "BP") 1950 - x$ages else x$ages
  new_baygmst_proxy(
    year = year,
    proxy = collapse_ensemble(comp, collapse),
    source = sprintf(
      "composite (%d members, %s; ages %s)", ncol(comp), collapse, age_units
    ),
    n_ensemble = ncol(comp)
  )
}

#' Internal constructor for `"baygmst_proxy"` objects
#'
#' Sorts by year, drops non-finite proxy values (with a message), and
#' validates that years are unique.
#' @noRd
new_baygmst_proxy <- function(year, proxy, source, n_ensemble) {
  year <- as.numeric(year)
  proxy <- as.numeric(proxy)
  keep <- is.finite(proxy) & is.finite(year)
  if (!all(keep)) {
    message(
      sum(!keep), " time step(s) with no finite proxy value dropped ",
      "(e.g. bins outside the composite's coverage)."
    )
  }
  year <- year[keep]
  proxy <- proxy[keep]
  if (length(year) == 0) {
    stop("No finite proxy values remain after dropping missing time steps.",
         call. = FALSE)
  }
  if (anyDuplicated(year)) {
    stop("Proxy input has duplicated years; each year must appear once.",
         call. = FALSE)
  }
  ord <- order(year)
  structure(
    data.frame(year = year[ord], proxy = proxy[ord]),
    source = source,
    n_ensemble = as.integer(n_ensemble),
    class = c("baygmst_proxy", "data.frame")
  )
}

#' @noRd
collapse_ensemble <- function(mat, collapse) {
  fun <- if (collapse == "median") stats::median else mean
  apply(as.matrix(mat), 1, fun, na.rm = TRUE)
}

#' @export
print.baygmst_proxy <- function(x, ...) {
  cat("<baygmst_proxy>\n")
  cat(sprintf("  Source:  %s\n", attr(x, "source")))
  cat(sprintf("  Years:   %s to %s CE (n = %d)\n",
              format(min(x$year)), format(max(x$year)), nrow(x)))
  if (isTRUE(attr(x, "n_ensemble") > 1)) {
    cat(sprintf("  Ensemble: collapsed from %d members\n",
                attr(x, "n_ensemble")))
  }
  cat(sprintf("  Range:   [%.3g, %.3g]\n", min(x$proxy), max(x$proxy)))
  invisible(x)
}

#' Validate the aligned model input vectors
#'
#' Shared by [fit_baygmst()] and [cv_baygmst()], after `proxy` has been
#' resolved to a bare numeric vector. Raises informative errors (rather
#' than `stopifnot()`'s terse condition text).
#' @noRd
check_model_inputs <- function(proxy, instrumental_T, forcing_G, forcing_V,
                               forcing_S, NT) {
  lens <- c(
    proxy = length(proxy),
    instrumental_T = length(instrumental_T),
    forcing_G = length(forcing_G),
    forcing_V = length(forcing_V),
    forcing_S = length(forcing_S)
  )
  bad_len <- lens != NT
  if (any(bad_len)) {
    stop(
      "All input series must have one value per year in `years` (length ",
      NT, "), but: ",
      paste(sprintf("%s has length %d", names(lens)[bad_len],
                    lens[bad_len]), collapse = "; "),
      ".",
      call. = FALSE
    )
  }
  nas <- c(
    proxy = anyNA(proxy),
    forcing_G = anyNA(forcing_G),
    forcing_V = anyNA(forcing_V),
    forcing_S = anyNA(forcing_S)
  )
  if (any(nas)) {
    stop(
      "Missing values are not allowed in ",
      paste(names(nas)[nas], collapse = ", "),
      " (only `instrumental_T` may contain NA, marking the years to ",
      "reconstruct).",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Resolve a proxy argument against the model's year vector
#'
#' Shared by [fit_baygmst()] and [cv_baygmst()]. A bare numeric vector is
#' used positionally (the original interface); anything else is coerced
#' with [as_baygmst_proxy()] and aligned to `years` by calendar year,
#' erroring informatively if any requested year is not covered.
#' @noRd
resolve_proxy <- function(proxy, years) {
  if (is.numeric(proxy) && is.null(dim(proxy))) {
    return(proxy)
  }
  proxy <- as_baygmst_proxy(proxy)
  idx <- match(years, proxy$year)
  if (anyNA(idx)) {
    missing_years <- years[is.na(idx)]
    stop(
      sprintf(
        "The proxy series (%s) does not cover %d of the requested years ",
        attr(proxy, "source"), length(missing_years)
      ),
      sprintf(
        "(e.g. %s). Restrict `years` to %s-%s or supply a proxy covering ",
        paste(utils::head(missing_years, 3), collapse = ", "),
        format(min(proxy$year)), format(max(proxy$year))
      ),
      "the full span.",
      call. = FALSE
    )
  }
  proxy$proxy[idx]
}
