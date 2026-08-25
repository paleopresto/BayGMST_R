# Legacy scripts (superseded)

These are the original, **unmodified** analysis scripts this repository
consisted of before the `BayGMST` package restructuring (branch `CRAN`,
2026-08). They are kept here for provenance only, are not part of the
package's public API, are not sourced by any exported function, and are not
run by `R CMD check`.

| Original file | Superseded by |
|---|---|
| `R_scripts/BayGMST_v1.0.R` | `fit_baygmst()`, `reconstruct()`, `plot_reconstruction()`, `plot_trace()`, `plot_posterior_densities()`, `transform_forcings()` |
| `R_scripts/cv_v0.1.R` | `cv_baygmst()`, `plot_cv()` |
| `R_scripts/calc_historic_annual_global_volcanismAOD.R` | Not ported -- reads NetCDF files from a specific contributor's local `~/Downloads`, not general-purpose. Kept for reference only. |
| `utils/PAGES2k_reducedProxy_UNSC.R` | `reduce_proxies()` |
| `utils/PAGES2k_reducedProxy_UNSC_OLD.R` | Already superseded by the file above *before* this restructuring; kept only because it was already in the repository. |
| `utils/PAGES2k_datagrabber.py` | Not ported (Python, outside an R package's scope). Kept for reference only. |
| `BayGMST_v1.0.stan` | `src/stan/baygmst.stan` (content-identical port, see `NEWS.md`) |
| `BayGMST_v1.0_5fcv.stan` | `src/stan/baygmst_cv.stan` (content-identical port, see `NEWS.md`) |
| `BayGMST_v1.0`, `BayGMST_v1.0_5fcv` | Not ported -- these were **compiled CmdStan binaries** (~2 MB each) that had been committed by accident; removed from the repository tree in 2026-08 (they remain in git history). |
| `config.yml` | The configuration file the original scripts read (moved here from the repository root). No package function reads it; its keys map onto function arguments now (e.g. `rp_method` -> `reduce_proxies(method = )`, `stan_params` -> `fit_baygmst()` sampling arguments). Note it contains a machine-specific `cmdstan_path` from the original author's machine. Excluded from the built package. |

See `NEWS.md` at the package root for the specific bugs fixed and dead code
dropped when porting these into the package's `R/` functions.
