# Changelog

## BayGMST 0.1.0

Initial packaging of the BayGMST model. Previously this repository was a
collection of top-level analysis scripts (`R_scripts/`, `utils/`) driven
by `config.yml`; this release converts that logic into a documented,
installable R package with an explicit function API. The original
scripts are preserved, unchanged, under `inst/legacy-scripts/` for
provenance.

### New

- [`transform_forcings()`](https://paleopresto.github.io/BayGMST_R/reference/transform_forcings.md)
  – ported from the normalization block in `BayGMST_v1.0.R`.
- [`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
  – ported from `utils/PAGES2k_reducedProxy_UNSC.R`. Supports `"PCR"`,
  `"LASSO"`, `"sPCR"`, and `"SPLS"`. `"SIR"` is included but marked
  experimental (see “Known limitations” below).
- [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
  – ported from the data-prep and `cmdstanr` sampling block in
  `BayGMST_v1.0.R`. Returns a `baygmst_fit` object instead of relying on
  variables left in the global environment.
- [`reconstruct()`](https://paleopresto.github.io/BayGMST_R/reference/reconstruct.md)
  – extracts a tidy reconstruction `data.frame` from a `baygmst_fit`
  (equivalent to the `gmst_reconstruction_data.csv` produced by the
  original script), without writing to disk.
- [`cv_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/cv_baygmst.md)
  – ported from the k-fold loop in `cv_v0.1.R`.
- [`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md),
  [`plot_trace()`](https://paleopresto.github.io/BayGMST_R/reference/plot_trace.md),
  [`plot_posterior_densities()`](https://paleopresto.github.io/BayGMST_R/reference/plot_posterior_densities.md),
  [`plot_cv()`](https://paleopresto.github.io/BayGMST_R/reference/plot_cv.md)
  – return `ggplot` objects instead of calling `ggsave()` as a side
  effect.
- [`print.baygmst_fit()`](https://paleopresto.github.io/BayGMST_R/reference/print.baygmst_fit.md),
  [`summary.baygmst_fit()`](https://paleopresto.github.io/BayGMST_R/reference/summary.baygmst_fit.md),
  `print.baygmst_cv()` S3 methods.
- [`plot.baygmst_fit()`](https://paleopresto.github.io/BayGMST_R/reference/plot.baygmst_fit.md):
  `plot(fit)` reproduces the combined figure of the original script and
  the PReSto manuscript (reconstruction on top; grouped posterior
  densities of `alpha1`, the forcing sensitivities, and the AR(1)
  coefficients below). Requires ‘patchwork’ (Suggests).
- [`as_baygmst_proxy()`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md)
  – normalizes the shapes a reduced/composited proxy series can arrive
  in (the output of
  [`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md),
  a `data.frame`, a bare vector, a time-by-ensemble matrix, or a
  composite object from another package such as a ‘compositeR’
  `paleoComposite`, with ages in yr BP converted to calendar years and
  ensembles collapsed to a median series) into a year-aware
  `baygmst_proxy` object.
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
  and
  [`cv_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/cv_baygmst.md)
  accept these objects (or anything coercible) directly and align them
  to `years` by calendar year; the original bare-vector, positional
  interface is unchanged. Ensemble inputs are collapsed before fitting
  in this release (compositing uncertainty is not yet propagated); the
  object records `n_ensemble` for a future ensemble-aware pathway.
- [`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
  output now carries class `baygmst_rp` (same list structure as before)
  with a [`print()`](https://rdrr.io/r/base/print.html) method.
- [`summary()`](https://rdrr.io/r/base/summary.html) on a fitted model
  now reports convergence diagnostics (`max_rhat`, `min_ess_bulk`) and
  its print method warns prominently when `max_rhat > 1.05` (surfaced by
  a Holocene-scale stress test in which a poorly identified fit
  previously printed without complaint).
- [`plot_trace()`](https://paleopresto.github.io/BayGMST_R/reference/plot_trace.md)
  and
  [`plot_posterior_densities()`](https://paleopresto.github.io/BayGMST_R/reference/plot_posterior_densities.md)
  no longer emit “Dropping ‘draws_df’ class” warnings.
- New vignette `baygmst-model` describing the statistical model (data
  and process levels, priors, forcing provenance, and the mapping to the
  Stan code), ported from the hand-typeset draft reference manual’s
  “Statistical Model” section. It explicitly credits the reduced-proxy
  Bayesian framework of Barboza et al. (2014, 2019) and Wang (2020).

### Fixed relative to the original scripts

- `load_proxies()`’s
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) choices omitted
  `"sPCR"`, even though `config.yml` and the bundled
  `data/barboza_rps/RP_new_All_sPCR.csv` both support it – selecting
  `sPCR` via the original script would error.
  `reduce_proxies(method = "sPCR")` now works.
- Removed the hardcoded, machine-specific path
  `/Users/tylerbagwell/Documents/GitHub/BayGMST_R/data/forcings_with_prediction_HanWang.csv`
  from the cross-validation logic. That block
  (`Forcings.projections`/`df_prj`) was also dead code in the original
  script – built but never referenced again – so it was not ported.
- Removed the `rm(list = ls())` call at the top of
  `utils/PAGES2k_reducedProxy_UNSC.R`. A package function must never
  clear its caller’s global environment.
- Package functions no longer write files (CSV/PNG) by default. Callers
  decide whether and where to save outputs, per CRAN policy on functions
  writing to the filesystem.

### Modeling decisions and known limitations

The following were reviewed and signed off by the package authors before
the first CRAN release:

- **`y[NT_mis]`/`z[NT_mis]` indexing, kept as is.** In
  `src/stan/baygmst.stan`, the observed-period fitted-value block uses
  the count `NT_mis` as a vector index at `t == 1`. Because missing
  (pre-instrumental) years always occupy positions `1:NT_mis`, this is
  numerically equivalent to `y[idx_mis[NT_mis]]`/`z[idx_mis[NT_mis]]`
  and is not a live bug. It would be wrong only if missing and observed
  years were interleaved, and it affects only `y_ins_fitted[1]`, a
  posterior-predictive diagnostic, not the parameter posteriors or the
  reconstruction `y_mis`.
- **Instrumental temperature is treated as error-free**, as in the
  models of Barboza et al. (2014) and Wang (2020). See
  [`vignette("baygmst-model")`](https://paleopresto.github.io/BayGMST_R/articles/baygmst-model.md).
- **Cross-validation R^2** is the fully Bayesian posterior mean of
  Stan’s `r2_cv`, rather than the plug-in estimate of the original
  `cv_v0.1.R`. The two are close but not identical.
- **`reduce_proxies(method = "SIR")` is experimental and unvalidated**
  and warns on use. The original SIR implementation depended on legacy
  libraries that are no longer available, so it could not be used as a
  reference.
- **Volcanic forcing coefficient.** `vol_coef = 25` (W m^-2 per unit
  AOD) is attributed to Hansen et al. (2005); the value is corroborated
  via IPCC AR5.
- Several packages
  [`library()`](https://rdrr.io/r/base/library.html)-loaded by the
  original scripts (`car`, `fda`, `ggmap`, `maps`) had no evident
  corresponding usage in the code as written, so they were not carried
  into `Imports`/`Suggests`.
- `data/HadCRUT.5.1.0.0.analysis.anomalies.ensemble_mean.nc` (31 MB) was
  not read by any ported function and has been removed from the
  repository tree (it remains in git history), along with the compiled
  CmdStan binaries in `inst/legacy-scripts/` and the `outputs/` build
  artifacts. The legacy `config.yml` moved to
  `inst/legacy-scripts/config.yml`; no package function reads it.
- The hand-typeset draft reference manual (`reference-manual/`) was
  retired. Its model description lives in
  [`vignette("baygmst-model")`](https://paleopresto.github.io/BayGMST_R/articles/baygmst-model.md),
  and the reference manual is the one R builds from `man/`.
