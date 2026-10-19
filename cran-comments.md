## Resubmission

This is a resubmission. In response to the review:

* Added the method references to the Description field of DESCRIPTION, in
  the form authors (year) <doi:...>.
* The legacy scripts under inst/legacy-scripts/ (kept for provenance) no
  longer call rm(list = ls()), and the one that changes par() now saves it
  with oldpar <- par(...) and restores it with par(oldpar). No examples,
  vignettes, or tests change options(), par(), or the working directory.

We also fixed a bug found since the first submission: reduce_proxies()
screened records for missing data after zero-filling their gaps, so the
screen had no effect. It now screens first, as the original method did.

## Test environments

- local macOS (ARM64), R 4.5.2, with CmdStan 2.39.0
- win-builder, R-devel (without cmdstanr: 1 NOTE, described below)
- R-hub: Linux, macOS (ARM64), and Windows, R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

* 'cmdstanr' is in Suggests and is not on CRAN. It is available from
  the Stan r-universe repository, declared in `Additional_repositories`,
  following the packaging pattern of the 'instantiate' package (on CRAN)
  and of 'brms'. Model fitting requires 'cmdstanr' and a CmdStan
  installation; the package checks for both via
  `instantiate::stan_cmdstan_exists()` and stops with installation
  instructions if either is missing. All examples, tests, and vignette
  chunks that fit a Stan model are skipped when CmdStan is unavailable,
  so the package installs and checks cleanly without it.

* Possibly misspelled words in DESCRIPTION (GMST, forcings, irradiance,
  paleoclimate, pre, radiative) are standard climate-science terms;
  GMST is spelled out in the Description.
