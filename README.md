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

BayGMST depends on [cmdstanr](https://mc-stan.org/cmdstanr/), which is not on
CRAN and must be installed from Stan's own R-universe repository:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

Then install BayGMST itself (once published; for now, from source on this
branch):

```r
# install.packages("remotes")
remotes::install_local(".", dependencies = TRUE, build_vignettes = TRUE)
```

## Usage

```r
library(BayGMST)

rp <- reduce_proxies(proxy_matrix, years, temp_calib, calib_years, method = "PCR")
forcing <- transform_forcings(G = co2, V = volcanic, S = solar)

fit <- fit_baygmst(
  proxy = rp$composite$RP1,
  instrumental_T = instrumental_T,
  forcing_G = forcing$G, forcing_V = forcing$V, forcing_S = forcing$S,
  years = years
)

reconstruct(fit)
plot_reconstruction(fit)
summary(fit)
```

See `vignette("baygmst-intro", package = "BayGMST")` for a full worked
example using small bundled example datasets.

## Model

See the package-level help (`?BayGMST`) and the `reference-manual/`
directory for a full description of each function. The underlying model
follows Barboza et al. (2014, 2019) and Wang (2020); see
`inst/legacy-scripts/` for the original, pre-packaging scripts this was
built from.
