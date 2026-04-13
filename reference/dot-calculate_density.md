# Estimate weighted kernel density

Estimate weighted kernel density

## Usage

``` r
.calculate_density(w, x, method, adjust = 1, map = TRUE)
```

## Arguments

- w:

  Vector with weights for each observation

- x:

  Matrix with dimensions where to calculate the density from. Only the
  first two dimensions will be used

- method:

  Kernel density estimation method:

  - `ks`: Computes density using the `kde` function from the `ks`
    package.

  - `wkde`: Computes density using a modified version of the `kde2d`
    function from the `MASS` package to allow weights. Bandwidth
    selection from the `ks` package is used instead.

- adjust:

  Numeric value to adjust to bandwidth. Default: 1. Not available for
  `ks` method

- map:

  Whether to map densities to individual observations

## Value

If `map` is `TRUE`, a vector with corresponding densities for each
observation is returned. Otherwise, a list with the density estimates
from the selected method is returned.

## Author

Jose Alquicira-Hernandez (modified)
