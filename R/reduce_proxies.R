#' Reduce a multi-proxy network to a single representative proxy series
#'
#' Collapses a matrix of paleoclimate proxy records into one (or a small set
#' of) "representative proxy" series suitable for [fit_baygmst()], following
#' the segmented-composite approach of Barboza et al. (2019) as ported from
#' \code{utils/PAGES2k_reducedProxy_UNSC.R} (see
#' \code{system.file("legacy-scripts", package = "BayGMST")} for the
#' original). The proxy network is broken into overlapping time segments of
#' `chunk` years; within each segment, only proxies that already have data
#' are used, and are regressed against instrumental temperature during
#' `calib_years` using the chosen `method`. The final composite series
#' averages across all segments that cover a given year.
#'
#' @param proxy_matrix Numeric matrix or data.frame of proxy values, one row
#'   per year (in the same order as `years`), one column per proxy record.
#'   Do not include a year column here -- pass it separately via `years`.
#' @param years Integer vector of calendar years, one per row of
#'   `proxy_matrix`.
#' @param temp_calib Numeric vector of instrumental temperatures during the
#'   calibration period, one value per year in `calib_years`.
#' @param calib_years Integer vector of calendar years used for calibration
#'   (must be a subset of `years`, and `length(calib_years) ==
#'   length(temp_calib)`).
#' @param method One of `"PCR"` (principal component regression, the
#'   default), `"LASSO"`, `"sPCR"` (supervised PCR), `"SPLS"` (sparse partial
#'   least squares), or `"SIR"`. **`"SIR"` is experimental** -- see Details.
#' @param chunk Segment length, in years. Default `250`. Every segment (each
#'   of which extends from its own start year through to the *end* of the
#'   record, not just `chunk` years) must fully contain the calibration
#'   window, i.e. `chunk` must be small enough, relative to `length(years)`,
#'   that even the last, shortest segment still spans all of `calib_years` --
#'   otherwise `reduce_proxies()` errors rather than silently comparing
#'   mismatched vectors (this mirrors an explicit check present in the
#'   original script).
#' @param target_adj_r2 For `method = "PCR"`: the calibration-period adjusted
#'   R-squared at which to stop adding principal components (the smallest
#'   number of PCs reaching this threshold is used; if none reach it, the
#'   number maximizing adjusted R-squared is used). Default `0.70`.
#' @param max_na_frac Proxies missing more than this fraction of
#'   observations within a segment (evaluated both over the full segment and
#'   over its calibration window) are dropped from that segment. Default
#'   `0.05`.
#' @param n_components For `method = "sPCR"` only: number of supervised
#'   principal components to use for prediction. A single value (recycled
#'   across segments) or a vector with one value per segment. Default `3`.
#' @param spls_eta For `method = "SPLS"` only: grid of `eta` sparsity
#'   values passed to \code{spls::cv.spls()}. Default `seq(0.1, 0.9, 0.1)`.
#'
#' @details
#' # Method support
#' `"PCR"`, `"LASSO"`, and `"SPLS"` are fully adaptive (all tuning parameters
#' are chosen by calibration-period cross-validation within each segment) and
#' are considered supported. `"sPCR"` is supported but requires the caller to
#' choose `n_components` (it is not cross-validated in the underlying
#' \pkg{superpc} workflow, matching the original script).
#'
#' `"SIR"` is a simplified, **unvalidated** reimplementation using a single
#' sliced-inverse-regression direction from \pkg{dr}, with no per-segment
#' variable pre-selection. It intentionally does *not* reproduce the
#' original script's SIR branch, which called an undocumented
#' `edrSelec()` step from the non-CRAN package 'edrGraphicalTools' that the
#' original author flagged as never fully working ("NEED TO FIGURE OUT HOW TO
#' INSTALL edrGraphicalTools!!"). `config.yml` in the original repository
#' likewise documented SIR as "still under construction." Treat
#' `method = "SIR"` results with corresponding caution.
#'
#' # Required packages
#' `"LASSO"` requires \pkg{glmnet}; `"sPCR"` requires \pkg{superpc}; `"SPLS"`
#' requires \pkg{spls}; `"SIR"` requires \pkg{dr}. These are `Suggests`, not
#' hard dependencies of BayGMST, and are only required for the method you
#' actually call.
#'
#' # Reproducibility
#' `"LASSO"`, `"sPCR"`, and `"SPLS"` use randomized cross-validation
#' internally. Call \code{set.seed()} before calling this function if you
#' need reproducible output; `reduce_proxies()` itself never calls
#' `set.seed()`, so it never perturbs the caller's random number stream as a
#' side effect.
#'
#' @return An object of class `"baygmst_rp"`, a list with components:
#'   \item{segments}{A data.frame with columns `year`, `RP1`, ..., `RPns`
#'     (one column per time segment).}
#'   \item{composite}{A data.frame with columns `year` and `RP1`, the
#'     across-segment average -- the single representative proxy series to
#'     pass to [fit_baygmst()] (which accepts the whole object directly,
#'     via [as_baygmst_proxy()]).}
#'   \item{method}{The method used, echoed back.}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' years <- 1:500
#' proxy_matrix <- cbind(
#'   p1 = c(rep(NA_real_, 50), cumsum(rnorm(450))),
#'   p2 = cumsum(rnorm(500))
#' )
#' calib_years <- 401:500
#' temp_calib <- cumsum(rnorm(100)) / 10
#' rp <- reduce_proxies(
#'   proxy_matrix, years, temp_calib, calib_years,
#'   method = "PCR", chunk = 250
#' )
#' head(rp$composite)
#' }
#'
#' @export
reduce_proxies <- function(proxy_matrix,
                            years,
                            temp_calib,
                            calib_years,
                            method = c("PCR", "LASSO", "sPCR", "SPLS", "SIR"),
                            chunk = 250,
                            target_adj_r2 = 0.70,
                            max_na_frac = 0.05,
                            n_components = 3,
                            spls_eta = seq(0.1, 0.9, 0.1)) {
  method <- match.arg(method)

  proxy_matrix <- as.matrix(proxy_matrix)
  storage.mode(proxy_matrix) <- "double"
  stopifnot(
    is.numeric(proxy_matrix),
    length(years) == nrow(proxy_matrix),
    length(temp_calib) == length(calib_years),
    length(calib_years) > 0,
    all(calib_years %in% years)
  )

  if (method == "SIR") {
    if (!requireNamespace("dr", quietly = TRUE)) {
      stop("method = \"SIR\" requires the 'dr' package.", call. = FALSE)
    }
    warning(
      "method = \"SIR\" is an experimental, unvalidated reimplementation ",
      "that does not reproduce the original (incomplete) SIR branch. ",
      "See ?reduce_proxies Details before relying on its output.",
      call. = FALSE
    )
  }

  ny <- length(years)
  ns <- ceiling(ny / chunk)
  np <- ncol(proxy_matrix)

  # earliest available row (not calendar year) for each proxy column
  year_min <- vapply(seq_len(np), function(i) {
    finite <- which(is.finite(proxy_matrix[, i]))
    if (length(finite) == 0) Inf else min(finite)
  }, numeric(1))

  proxy_scaled <- scale(proxy_matrix)
  proxy_filled <- apply(proxy_scaled, 2, function(x) {
    seen_first <- cumsum(!is.na(x)) > 0
    x[is.na(x) & seen_first] <- 0
    x
  })
  rownames(proxy_filled) <- rownames(proxy_scaled)
  colnames(proxy_filled) <- colnames(proxy_scaled)

  is_calib <- years %in% calib_years
  if (sum(is_calib) != length(temp_calib)) {
    stop(
      "The number of years %in% calib_years (", sum(is_calib),
      ") does not match length(temp_calib) (", length(temp_calib), ").",
      call. = FALSE
    )
  }

  RP <- matrix(NA_real_, nrow = ny, ncol = ns)

  for (k in seq_len(ns)) {
    t_min       <- 1 + (k - 1) * chunk
    seg_proxies <- which(year_min <= t_min)
    time_span   <- which(seq_len(ny) >= t_min)

    if (length(seg_proxies) == 0) {
      warning(sprintf(
        "Segment %d (starting at row %d): no proxies available yet; skipping.",
        k, t_min
      ), call. = FALSE)
      next
    }

    proxy_seg <- proxy_filled[time_span, seg_proxies, drop = FALSE]
    keep_full <- colMeans(is.na(proxy_seg)) < max_na_frac
    proxy_seg <- proxy_seg[, keep_full, drop = FALSE]

    calib_idx_seg <- which(is_calib[time_span])
    proxy_calib   <- proxy_seg[calib_idx_seg, , drop = FALSE]
    keep_calib    <- colMeans(is.na(proxy_calib)) < max_na_frac
    proxy_seg     <- proxy_seg[, keep_calib, drop = FALSE]

    proxy_seg[!is.finite(proxy_seg)] <- 0

    if (ncol(proxy_seg) == 0) {
      warning(sprintf(
        "Segment %d: no proxies survived missingness filtering; skipping.", k
      ), call. = FALSE)
      next
    }

    if (length(calib_idx_seg) != length(temp_calib)) {
      stop(sprintf(
        paste(
          "Segment %d's time span does not fully contain the calibration",
          "window: found %d calibration years in this segment, but",
          "temp_calib/calib_years has %d. Every segment must fully contain",
          "the calibration window (this mirrors an explicit check in the",
          "original PAGES2k_reducedProxy_UNSC.R script) -- increase `chunk`",
          "so that even the last, shortest segment still spans the full",
          "calibration period. The simplest fix is a single segment:",
          "`chunk = length(years)`."
        ),
        k, length(calib_idx_seg), length(temp_calib)
      ), call. = FALSE)
    }

    n_comp_k <- if (length(n_components) == 1) n_components else n_components[k]

    RP[time_span, k] <- .reduce_segment(
      proxy_seg     = proxy_seg,
      calib_idx     = calib_idx_seg,
      temp_calib    = temp_calib,
      method        = method,
      target_adj_r2 = target_adj_r2,
      n_components  = n_comp_k,
      spls_eta      = spls_eta
    )
  }

  colnames(RP) <- paste0("RP", seq_len(ns))
  segments  <- data.frame(year = years, RP, check.names = FALSE)
  composite <- data.frame(year = years, RP1 = apply(RP, 1, mean, na.rm = TRUE))

  structure(
    list(segments = segments, composite = composite, method = method),
    class = "baygmst_rp"
  )
}

#' @export
print.baygmst_rp <- function(x, ...) {
  cat("<baygmst_rp>\n")
  cat(sprintf("  Method:   %s\n", x$method))
  cat(sprintf("  Years:    %s to %s (n = %d)\n",
              format(min(x$composite$year)), format(max(x$composite$year)),
              nrow(x$composite)))
  cat(sprintf("  Segments: %d\n", ncol(x$segments) - 1L))
  cat("  Pass to fit_baygmst() directly, or use x$composite$RP1.\n")
  invisible(x)
}

#' @param proxy_seg Numeric matrix (rows = years in this segment, columns =
#'   proxies available in this segment).
#' @param calib_idx Integer indices, into the rows of `proxy_seg`, of the
#'   calibration years.
#' @param temp_calib,method,target_adj_r2,n_components,spls_eta As in
#'   [reduce_proxies()].
#' @return Numeric vector, length `nrow(proxy_seg)`: the fitted/predicted
#'   reduced-proxy series for this segment.
#' @noRd
.reduce_segment <- function(proxy_seg, calib_idx, temp_calib, method,
                             target_adj_r2, n_components, spls_eta) {
  one <- rep(1, nrow(proxy_seg))

  if (method == "PCR") {
    pca <- stats::prcomp(proxy_seg, center = FALSE, scale. = FALSE)
    npcs_max <- min(ncol(pca$x), length(temp_calib) - 2)
    if (npcs_max < 1) {
      stop(
        "Not enough calibration observations to fit even one principal ",
        "component in this segment.", call. = FALSE
      )
    }
    adj_r2 <- vapply(seq_len(npcs_max), function(j) {
      pc_calib <- pca$x[calib_idx, seq_len(j), drop = FALSE]
      summary(stats::lm(temp_calib ~ pc_calib))$adj.r.squared
    }, numeric(1))
    hits <- which(adj_r2 >= target_adj_r2)
    npcs <- if (length(hits) > 0) hits[1] else which.max(adj_r2)
    pc_calib <- pca$x[calib_idx, seq_len(npcs), drop = FALSE]
    pc_all   <- pca$x[, seq_len(npcs), drop = FALSE]
    fit <- stats::lm(temp_calib ~ pc_calib)
    return(as.numeric(cbind(one, pc_all) %*% as.numeric(stats::coef(fit))))
  }

  if (method == "LASSO") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("method = \"LASSO\" requires the 'glmnet' package.", call. = FALSE)
    }
    x_calib  <- proxy_seg[calib_idx, , drop = FALSE]
    cv_lasso <- glmnet::cv.glmnet(x = x_calib, y = temp_calib, alpha = 1)
    return(as.numeric(stats::predict(cv_lasso, newx = proxy_seg, s = "lambda.min")))
  }

  if (method == "sPCR") {
    if (!requireNamespace("superpc", quietly = TRUE)) {
      stop("method = \"sPCR\" requires the 'superpc' package.", call. = FALSE)
    }
    datapc    <- list(x = t(proxy_seg[calib_idx, , drop = FALSE]), y = temp_calib)
    train_obj <- superpc::superpc.train(datapc, type = "regression")
    cv_obj    <- superpc::superpc.cv(train_obj, datapc)
    threshold <- cv_obj$thresholds[which.max(cv_obj$scor[1, ])]
    fit_obj <- superpc::superpc.predict(
      train_obj, datapc, list(x = t(proxy_seg)),
      threshold = threshold, n.components = n_components,
      prediction.type = "continuous"
    )
    pred_calib <- data.frame(temp_calib, fit_obj$v.pred[calib_idx, , drop = FALSE])
    fit_lm <- stats::lm(temp_calib ~ ., data = pred_calib)
    return(as.numeric(stats::predict(fit_lm, newdata = data.frame(fit_obj$v.pred))))
  }

  if (method == "SPLS") {
    if (!requireNamespace("spls", quietly = TRUE)) {
      stop("method = \"SPLS\" requires the 'spls' package.", call. = FALSE)
    }
    x_calib <- proxy_seg[calib_idx, , drop = FALSE]
    k_max   <- min(ncol(x_calib), floor(0.9 * nrow(x_calib) - 1))
    if (k_max < 1) {
      stop(
        "Not enough calibration observations/proxies for method = \"SPLS\" ",
        "in this segment.", call. = FALSE
      )
    }
    cv_spls  <- spls::cv.spls(x = x_calib, y = temp_calib, K = seq_len(k_max), eta = spls_eta)
    fit_spls <- spls::spls(x = x_calib, y = temp_calib, K = cv_spls$K.opt, eta = cv_spls$eta.opt)
    return(as.numeric(stats::predict(fit_spls, newx = proxy_seg)))
  }

  if (method == "SIR") {
    x_calib     <- proxy_seg[calib_idx, , drop = FALSE]
    sir_fit     <- dr::dr(temp_calib ~ x_calib, nslices = 8)
    direction   <- sir_fit$evectors[, 1, drop = FALSE]
    x_dir_calib <- as.matrix(x_calib %*% direction)
    fit_lm      <- stats::lm(temp_calib ~ x_dir_calib)
    x_dir_all   <- as.matrix(proxy_seg %*% direction)
    return(as.numeric(cbind(one, x_dir_all) %*% as.numeric(stats::coef(fit_lm))))
  }

  stop("Unreachable: unknown method.", call. = FALSE) # nocov
}
