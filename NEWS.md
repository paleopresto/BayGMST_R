# BayGMST 0.1.0

Initial packaging of the BayGMST model. Previously this repository was a
collection of top-level analysis scripts (`R_scripts/`, `utils/`) driven by
`config.yml`; this release converts that logic into a documented, installable
R package with an explicit function API. The original scripts are preserved,
unchanged, under `inst/legacy-scripts/` for provenance.

## New

* `transform_forcings()` -- ported from the normalization block in
  `BayGMST_v1.0.R`.
* `reduce_proxies()` -- ported from `utils/PAGES2k_reducedProxy_UNSC.R`.
  Supports `"PCR"`, `"LASSO"`, `"sPCR"`, and `"SPLS"`. `"SIR"` is included but
  marked experimental (see "Known limitations" below).
* `fit_baygmst()` -- ported from the data-prep and `cmdstanr` sampling block
  in `BayGMST_v1.0.R`. Returns a `baygmst_fit` object instead of relying on
  variables left in the global environment.
* `reconstruct()` -- extracts a tidy reconstruction `data.frame` from a
  `baygmst_fit` (equivalent to the `gmst_reconstruction_data.csv` produced by
  the original script), without writing to disk.
* `cv_baygmst()` -- ported from the k-fold loop in `cv_v0.1.R`.
* `plot_reconstruction()`, `plot_trace()`, `plot_posterior_densities()`,
  `plot_cv()` -- return `ggplot`/`patchwork` objects instead of calling
  `ggsave()` as a side effect.
* `print.baygmst_fit()`, `summary.baygmst_fit()`, `print.baygmst_cv()` S3
  methods.

## Fixed relative to the original scripts

* `load_proxies()`'s `match.arg()` choices omitted `"sPCR"`, even though
  `config.yml` and the bundled `data/barboza_rps/RP_new_All_sPCR.csv` both
  support it -- selecting `sPCR` via the original script would error.
  `reduce_proxies(method = "sPCR")` now works.
* Removed the hardcoded, machine-specific path
  `/Users/tylerbagwell/Documents/GitHub/BayGMST_R/data/forcings_with_prediction_HanWang.csv`
  from the cross-validation logic. That block (`Forcings.projections`/`df_prj`)
  was also dead code in the original script -- built but never referenced
  again -- so it was not ported.
* Removed the `rm(list = ls())` call at the top of
  `utils/PAGES2k_reducedProxy_UNSC.R`. A package function must never clear its
  caller's global environment.
* Package functions no longer write files (CSV/PNG) by default. Callers
  decide whether and where to save outputs, per CRAN policy on functions
  writing to the filesystem.

## Known limitations (flagged for Tyler/Julien, not silently resolved)

* **Unresolved possible bug carried forward unchanged:** in
  `src/stan/baygmst.stan` (originally `BayGMST_v1.0.stan`), the
  observed-period fitted-value block uses `z[NT_mis]` as an index into the
  proxy vector at `t == 1`. `NT_mis` is a *count* of missing years, not a time
  index, so this looks like it should likely be `z[idx_mis[NT_mis]]` (the
  proxy value at the last pre-instrumental year) instead. The original author
  flagged this with `// CHECK THIS LINE` and never resolved it. This
  restructuring does not change model math without statistical sign-off, so
  the line is unchanged -- see the reference manual for the exact location.
* `reduce_proxies(method = "SIR")` mirrors the original script's SIR branch,
  which `config.yml` itself already documented as "still under construction."
  It is *not* recommended for use until validated.
* Several packages `library()`-loaded by the original scripts (`car`, `fda`,
  `ggmap`, `maps`) had no evident corresponding usage in the code as written,
  so they were not carried into `Imports`/`Suggests`. Re-add them if a hidden
  use is found.
* `data/HadCRUT.5.1.0.0.analysis.anomalies.ensemble_mean.nc` (31 MB) is not
  read by any ported function and is excluded from the installable package.
  Decide whether to keep it in the repo at all.
