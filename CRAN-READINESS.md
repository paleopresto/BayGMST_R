# CRAN readiness checklist

**Status (2026-09-29):** the package passes `R CMD check --as-cran`
locally (macOS ARM64, R 4.5.2, CmdStan 2.39.0) with no ERRORs; all
tests, examples (including `--run-donttest`), and vignettes (including
real sampling) succeed. The remaining WARNINGs/NOTEs are documented in
`cran-comments.md`. Tyler Bagwell has confirmed the package reproduces
the results of the original scripts.

## Resolved by the authors (August 2026)

All six open items from the draft reference manual were reviewed by
Tyler Bagwell, Frederi Viens, and Julien Emile-Geay; the decisions are
recorded in `NEWS.md` under “Modeling decisions and known limitations”:

1.  `y[NT_mis]`/`z[NT_mis]` indexing shorthand: kept as is.
2.  Error-free instrumental temperature: kept, documented in the model
    vignette.
3.  Fully Bayesian cross-validation R^2: kept.
4.  `SIR` unvalidated, warns on use: kept.
5.  `vol_coef = 25` attributed to Hansen et al. (2005): kept.
6.  `DESCRIPTION` wording: replaced with the text Tyler and Frederi
    supplied.

Also resolved: `plot(fit)` now draws the combined reconstruction +
posterior-density figure used in the PReSto manuscript, and the
hand-typeset draft reference manual was retired (its model section lives
in
[`vignette("baygmst-model")`](https://paleopresto.github.io/BayGMST_R/articles/baygmst-model.md);
R builds the real manual from `man/`).

## Remaining before submission

**Make `paleopresto/BayGMST_R` public.** The `URL`/`BugReports` fields
404 while it is private, which fails CRAN’s URL check. Going public also
triggers the pkgdown deploy (`.github/workflows/pkgdown.yaml`).

Merge the `CRAN` branch into `main` (or decide which branch is
canonical) so the public repo matches the submitted package.

Run `devtools::spell_check()`, then the remote checks
`devtools::check_win_devel()` and `rhub::rhub_check()`. Win-builder
emails the maintainer (Nick McKay).

Re-run `R CMD check --as-cran` with the remote incoming checks on (after
the repo is public) and update `cran-comments.md` with the
win-builder/R-hub environments.

Nick McKay, as maintainer, submits via
<https://cran.r-project.org/submit.html> and confirms the email CRAN
sends to <nick@nau.edu>.

## Notes for working on the package

- Do not use `devtools::load_all()` or `devtools::test()`; `instantiate`
  is incompatible with
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html).
  Install for real (`R CMD INSTALL .`) and use `R CMD check`’s own test
  run.
- Every example, test, and vignette chunk that samples a Stan model is
  guarded with
  [`instantiate::stan_cmdstan_exists()`](https://wlandau.github.io/instantiate/reference/stan_cmdstan_exists.html),
  so the package checks cleanly on machines without CmdStan.
