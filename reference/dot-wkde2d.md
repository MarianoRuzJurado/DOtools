# Weighted 2D kernel density estimation

Weighted 2D kernel density estimation

## Usage

``` r
.wkde2d(x, y, w, h, adjust = 1, n = 100, lims = c(range(x), range(y)))
```

## Arguments

- x:

  Dimension 1

- y:

  Dimension 2

- w:

  Weight variable

- h:

  vector of bandwidths for x and y directions. Defaults to normal
  reference bandwidth (ks::hpi). A scalar value will be taken to apply
  to both directions.

- adjust:

  Bandwidth adjustment

- n:

  Number of grid points in each direction. Can be scalar or a length-2
  integer vector.

- lims:

  The limits of the rectangle covered by the grid as c(xl, xu, yl, yu).

## Value

A list of three components.

- `x, y` The x and y coordinates of the grid points, vectors of length
  n.

- `z` A matrix of the weighted estimated density: rows correspond to the
  value of x, columns to the value of y.

## Author

Jose Alquicira-Hernandez (modified)
