# DO.PyEnv

Sets up or connects to a conda Python environment for use with DOtools.
If no environment path is provided, it will create one at
`~/.venv/DOtools` and install required Python packages: `scvi-tools`,
`celltypist`, and `scanpro`.

## Usage

``` r
DO.PyEnv(conda_path = NULL)
```

## Arguments

- conda_path:

  character string specifying the path to an existing or new conda
  environment.

## Value

None

## Examples

``` r
# Creates DOtools environment at ~/.venv/DOtools if it doesn't exist
DO.PyEnv()
#> 2026-01-21 15:41:18 - Using existing conda environment at: /home/mariano/.venv/DOtools
#> 2026-01-21 15:41:18 - Python packages ready for DOtools!

# Use an existing conda environment at a custom location
# DO.PyEnv(conda_path = "~/miniconda3/envs/my_dotools_env")
```
