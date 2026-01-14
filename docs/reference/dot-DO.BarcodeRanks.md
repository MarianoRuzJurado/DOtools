# DO.BarcodeRanks

Given a raw count matrix (e.g. from a CellRanger HDF5 file), estimate
the number of expected cells and droplets using the knee and inflection
points from barcodeRanks.

## Usage

``` r
.DO.BarcodeRanks(SCE_obj)
```

## Arguments

- SCE_obj:

  A Single cell experiment object.

## Value

A named numeric vector: \`c(xpc_cells = ..., total_cells = ...)

## Author

Mariano Ruz Jurado
