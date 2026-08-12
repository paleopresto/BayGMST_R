# CRAN readiness checklist

This package was restructured on the `CRAN` branch in a sandbox with **no R
installation available** -- every `.R`/`.stan`/`DESCRIPTION` file was
hand-authored and statically cross-checked (grep, not `R CMD check`), but
nothing here has actually been run. Treat this checklist as the bridge
between "carefully hand-authored" and "actually verified."

## 1. Regenerate documentation for real

`NAMESPACE` and every file under `man/` are hand-written placeholders that
mirror what `roxygen2` *should* produce from the `@export`/`@param`/etc. tags
in `R/`, but they were never actually run through roxygen2.

```r
install.packages(c("devtools", "roxygen2"))
devtools::document()
```

Then diff what changed. If `devtools::document()` produces something
different from what's already in `man/`/`NAMESPACE`, **trust the generated
output** -- it's authoritative; the hand-authored versions were only ever a
best-effort approximation.

## 2. Set up the Stan build scaffolding

Stan model files live in `src/stan/baygmst.stan` and
`src/stan/baygmst_cv.stan`, following the
[instantiate](https://wlandau.github.io/instantiate/) package's packaging
pattern. The five scaffold files instantiate needs
(`cleanup`, `cleanup.win`, `src/Makevars`, `src/Makevars.win`,
`src/install.libs.R`) were **deliberately not hand-written** -- they carry
exact, license-attributed shell/Makevars syntax that's easy to get subtly
wrong without being able to test it. Generate them for real:

```r
install.packages("instantiate")
instantiate::stan_package_configure()
```

Run this from the package root, then inspect what it wrote before
committing it.

## 3. Install and smoke-test

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()

devtools::install(build_vignettes = TRUE)
library(BayGMST)
instantiate::stan_cmdstan_exists() # should be TRUE now
```

**Do not use `devtools::load_all()` or `devtools::test()`** -- instantiate's
own documentation states `pkgload::load_all()` is incompatible with it. Use
real installs (`devtools::install()` / `R CMD INSTALL`) plus `R CMD check`'s
own test run instead.

## 4. Build and check

```sh
R CMD build .
R CMD check --as-cran BayGMST_0.1.0.tar.gz
```

Every example/test/vignette chunk that actually samples a model is already
guarded with `instantiate::stan_cmdstan_exists()`, so `R CMD check` should
not hard-fail in an environment without CmdStan (e.g. CRAN's own check
machines) -- but this has not been verified by actually running it.

## 5. Fill in placeholders that need your judgment, not mine

- **`DESCRIPTION` `Description:` field** -- drafted from the README and this
  package's own code. Read it and edit the wording to your satisfaction;
  it's public CRAN-facing text.
- **`RoxygenNote` in `DESCRIPTION`** -- set to a guessed placeholder
  (`7.3.2`); `devtools::document()` (step 1) will set the real value
  automatically.
- **`Authors@R` / `Maintainer`** -- currently Tyler Bagwell (`cre`) and
  Julien Emile-Geay (`aut`), per this session's discussion. Add/remove
  contributors as appropriate.
- **`URL`/`BugReports`** -- currently point at
  `github.com/julieneg/BayGMST_R`; update if the package moves to a
  different repository (e.g. under Tyler's own account) before submission.

## 6. Known unresolved issues (see `NEWS.md` for full detail)

- **Possible bug, left unchanged:** `src/stan/baygmst.stan`'s
  `y_ins_fitted` block indexes the proxy vector with `z[NT_mis]` (a *count*)
  at `t == 1`, where `z[idx_mis[NT_mis]]` looks more likely correct. Flagged
  by the original author, never resolved, not silently changed here --
  needs your statistical sign-off either way.
- **`reduce_proxies(method = "SIR")`** is experimental/unvalidated; see its
  documentation (`?reduce_proxies`) before using it for anything real.
- **`data/HadCRUT.5.1.0.0.analysis.anomalies.ensemble_mean.nc`** (31 MB, in
  the repo root) isn't read by any packaged function and is excluded from
  the built package via `.Rbuildignore`. Decide whether to keep it in the
  repo at all, or drop it (e.g. with `git lfs` or a data-download script if
  it's needed for some other purpose not captured in this restructuring).

## 7. Only after all of the above

- `devtools::build_manual()` (or `R CMD Rd2pdf .`) to produce the real,
  authoritative reference manual PDF, superseding
  `reference-manual/BayGMST-reference-manual.tex` (a hand-typeset preview
  built without access to R -- see the banner at the top of that document).
- `devtools::spell_check()`, `devtools::check_win_devel()` /
  `rhub::rhub_check()` before an actual CRAN submission.
