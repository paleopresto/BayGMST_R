# Transform raw radiative forcings for use with `fit_baygmst()`

Applies the forcing transformations used throughout the BayGMST
pipeline: volcanic forcing is mapped to a negative saturating form,
greenhouse-gas (CO2) forcing is log-transformed relative to a reference
concentration, and solar forcing is centered on its own mean. Ported
from the normalization block near the top of the original
`BayGMST_v1.0.R` script (see
`system.file("legacy-scripts", package = "BayGMST")`).

## Usage

``` r
transform_forcings(G, V, S, co2_c0 = 280, co2_coef = 5.35, vol_coef = 25)
```

## Arguments

- G:

  Numeric vector of raw greenhouse-gas forcing (e.g. atmospheric CO2
  concentration in ppm).

- V:

  Numeric vector of raw volcanic forcing (e.g. aerosol optical depth).
  Same length as `G`.

- S:

  Numeric vector of raw solar forcing (e.g. total solar irradiance).
  Same length as `G`.

- co2_c0:

  Reference/baseline CO2 concentration used in the log transform.
  Default `280` (a common preindustrial reference, ppm).

- co2_coef:

  Multiplier applied to the log CO2 ratio. Default `5.35`.

- vol_coef:

  Multiplier applied to the volcanic saturating transform. Default `25`.

## Value

A list with three numeric vectors, `G`, `V`, and `S` – the transformed
forcings, each the same length as the corresponding input.

## Details

The transforms applied are: \$\$V' = -\|\mathrm{vol\\coef}\| \times (1 -
e^{-V})\$\$ \$\$G' = \mathrm{co2\\coef} \times \log(G /
\mathrm{co2\\c0})\$\$ \$\$S' = S - \bar{S}\$\$

## Examples

``` r
fc <- transform_forcings(
  G = c(280, 300, 340, 400),
  V = c(0, 0.05, 0.02, 0.10),
  S = c(1360, 1361, 1360.5, 1361.2)
)
fc$G
#> [1] 0.0000000 0.3691119 1.0387347 1.9082110
```
