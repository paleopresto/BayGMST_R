# BayGMST

R package implementing the BayGMST model of Bagwell et al. (in prep.): a
Bayesian hierarchical AR(1) state-space model reconstructing global mean
surface temperature (GMST) from a reduced paleoclimate proxy network,
instrumental observations, and radiative forcings.

**Status:** restructured (branch `CRAN`) from a collection of top-level
analysis scripts into an installable R package aimed at eventual CRAN
submission. See `NEWS.md` for exactly what changed in that restructuring, and
`CRAN-READINESS.md` for the remaining steps before an actual submission
(this restructuring was done in an environment with no R installed, so
several files -- `NAMESPACE`, `man/*.Rd` -- are hand-authored placeholders
that must be regenerated and verified with real R tooling; see that file).

## Installation

Fitting models requires [cmdstanr](https://mc-stan.org/cmdstanr/), which is
not on CRAN, and a CmdStan installation. Install both **before** BayGMST:
BayGMST compiles its Stan models when it is installed, and skips that step if
they are missing (proxy reduction and forcing transforms still work). If you
add them later, reinstall BayGMST.

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

Then install BayGMST itself:

```r
install.packages("BayGMST")                      # once on CRAN
# remotes::install_github("paleopresto/BayGMST_R") # development version
```

## Usage

```r
library(BayGMST)

rp <- reduce_proxies(proxy_matrix, years, temp_calib, calib_years, method = "PCR")
forcing <- transform_forcings(G = co2, V = volcanic, S = solar)

fit <- fit_baygmst(
  proxy = rp,
  instrumental_T = instrumental_T,
  forcing_G = forcing$G, forcing_V = forcing$V, forcing_S = forcing$S,
  years = years
)

reconstruct(fit)
plot(fit)   # reconstruction + posterior densities, as in the PReSto manuscript
summary(fit)
```

The proxy input does not have to come from `reduce_proxies()`:
`as_baygmst_proxy()` accepts bare vectors, `data.frame`s, time-by-ensemble
matrices, and composite objects from other packages (such as
[compositeR](https://github.com/nickmckay/compositeR), the intended primary
pathway once it is published), converting ages in yr BP and collapsing
ensembles as needed. See the "Using composites from other packages" section
of the vignette.

See `vignette("baygmst-intro", package = "BayGMST")` for a full worked
example using small bundled example datasets.

Note: `inst/legacy-scripts/config.yml` is the configuration file of the
original pre-package scripts, kept only for provenance; no package function
reads it.

## Model

The underlying model builds on the reduced-proxy Bayesian framework of
Barboza et al. (2014, 2019) and Wang (2020); see
`vignette("baygmst-model", package = "BayGMST")` for the full statistical
description (equations, priors, forcing data provenance, and the mapping
to the Stan code). For per-function documentation, see the package help
(`?BayGMST`); `inst/legacy-scripts/` holds the original, pre-packaging
scripts this was built from.
