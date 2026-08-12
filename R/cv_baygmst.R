#' K-fold cross-validate the BayGMST model over the instrumental period
#'
#' Partitions the instrumental (non-`NA`) years of `instrumental_T` into
#' `nfold` contiguous, roughly equal blocks; for each block in turn, refits
#' [fit_baygmst()]'s underlying model with that block's temperatures hidden
#' (treated as missing/to-be-reconstructed), and reports how well the model
#' recovers them. Ported from the fold loop in the original \code{cv_v0.1.R}
#' script (see \code{system.file("legacy-scripts", package = "BayGMST")}).
#'
#' Two blocks present in the original script were **not** ported: a
#' `Forcings.projections`/RCP-scenario block that was built but never
#' referenced again (dead code, and it hardcoded a path to a specific
#' contributor's machine), and a large commented-out "OLD" duplicate of the
#' fold loop at the end of the file. See `NEWS.md`.
#'
#' @inheritParams fit_baygmst
#' @param nfold Number of folds. Default `3`, matching the original script.
#'
#' @return An object of class `"baygmst_cv"`, a list with elements:
#'   \item{folds}{A single combined data.frame (one row per held-out year,
#'     across all folds) with columns `year`, `fold`, `T_true`, `T_mean`,
#'     `T_lo_inner`, `T_hi_inner`, `T_lo_outer`, `T_hi_outer`.}
#'   \item{r2}{Named numeric vector, the posterior mean of the Stan-computed
#'     `r2_cv` generated quantity for each fold (fully Bayesian; this differs
#'     slightly from the plug-in point-estimate R-squared computed in the
#'     original R script -- see Details in the package vignette).}
#'   \item{mse}{As `r2`, for `mse_cv`.}
#'   \item{years, call}{Echoed back.}
#'
#' @section CmdStan required:
#' See [fit_baygmst()].
#'
#' @seealso [plot_cv()]
#' @export
cv_baygmst <- function(proxy,
                        instrumental_T,
                        forcing_G,
                        forcing_V,
                        forcing_S,
                        years,
                        nfold = 3,
                        chains = 4,
                        parallel_chains = 2,
                        iter_warmup = 500,
                        iter_sampling = 1500,
                        seed = NULL,
                        ...) {
  NT <- length(years)
  stopifnot(
    length(proxy) == NT,
    length(instrumental_T) == NT,
    length(forcing_G) == NT,
    length(forcing_V) == NT,
    length(forcing_S) == NT,
    !anyNA(proxy), !anyNA(forcing_G), !anyNA(forcing_V), !anyNA(forcing_S),
    nfold >= 2
  )
  if (iter_sampling < 1000) {
    stop("iter_sampling must be at least 1000.", call. = FALSE)
  }

  idx_full_obs <- which(!is.na(instrumental_T))
  if (length(idx_full_obs) < nfold) {
    stop("Fewer instrumental observations than folds.", call. = FALSE)
  }
  nobs_fold <- floor(length(idx_full_obs) / nfold)

  model <- instantiate::stan_package_model(name = "baygmst_cv", package = "BayGMST")

  fold_frames <- vector("list", nfold)
  r2  <- stats::setNames(numeric(nfold), paste0("fold", seq_len(nfold)))
  mse <- stats::setNames(numeric(nfold), paste0("fold", seq_len(nfold)))

  for (i in seq_len(nfold)) {
    fold_start <- (i - 1) * nobs_fold + 1
    fold_end   <- if (i == nfold) length(idx_full_obs) else i * nobs_fold
    idx_cv     <- idx_full_obs[fold_start:fold_end]

    y_fold <- instrumental_T
    y_cv_true <- y_fold[idx_cv]
    y_fold[idx_cv] <- NA

    idx_obs <- which(!is.na(y_fold))
    idx_mis <- which(is.na(y_fold))

    data_list <- list(
      NT      = NT,
      NT_obs  = length(idx_obs),
      NT_mis  = length(idx_mis),
      NT_cv   = length(idx_cv),
      idx_obs = as.integer(idx_obs),
      idx_mis = as.integer(idx_mis),
      idx_cv  = as.integer(idx_cv),
      G       = as.vector(forcing_G),
      S       = as.vector(forcing_S),
      V       = as.vector(forcing_V),
      y_obs   = as.vector(y_fold[idx_obs]),
      y_cv_true = as.vector(y_cv_true),
      z       = as.vector(proxy)
    )

    fit_i <- model$sample(
      data            = data_list,
      chains          = chains,
      parallel_chains = parallel_chains,
      iter_warmup     = iter_warmup,
      iter_sampling   = iter_sampling,
      seed            = seed,
      ...
    )

    r2[i]  <- fit_i$summary(variables = "r2_cv")$mean
    mse[i] <- fit_i$summary(variables = "mse_cv")$mean

    draws     <- fit_i$draws("y_mis")
    mat_names <- paste0("y_mis[", match(idx_cv, idx_mis), "]")
    mat       <- posterior::as_draws_matrix(draws)[, mat_names, drop = FALSE]

    fold_frames[[i]] <- data.frame(
      year       = years[idx_cv],
      fold       = i,
      T_true     = y_cv_true,
      T_mean     = apply(mat, 2, mean),
      T_lo_inner = apply(mat, 2, stats::quantile, 0.16),
      T_hi_inner = apply(mat, 2, stats::quantile, 0.84),
      T_lo_outer = apply(mat, 2, stats::quantile, 0.025),
      T_hi_outer = apply(mat, 2, stats::quantile, 0.975)
    )
  }

  structure(
    list(
      folds = do.call(rbind, fold_frames),
      r2    = r2,
      mse   = mse,
      years = years,
      call  = match.call()
    ),
    class = "baygmst_cv"
  )
}

#' @export
print.baygmst_cv <- function(x, ...) {
  cat("<baygmst_cv>\n")
  cat(sprintf("  Folds:  %d\n", length(x$r2)))
  cat(sprintf("  R2 (posterior mean, per fold): %s\n",
              paste(sprintf("%.2f", x$r2), collapse = ", ")))
  cat(sprintf("  Median R2: %.2f\n", stats::median(x$r2)))
  invisible(x)
}
