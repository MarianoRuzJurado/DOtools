# Plot density estimates

Plot density estimates

## Usage

``` r
.plot_density_(
  z,
  dens,
  feature,
  cell_embeddings,
  dim_names,
  shape = 16,
  size = 1,
  text_size = 14,
  density_quantile_threshold = 0.5,
  n_bands = 256,
  legend.position = "right",
  legend_title = "Density",
  pal = "Reds",
  raster = FALSE,
  ...
)
```

## Arguments

- z:

  Vector with density values for each cells

- dens:

  density grid retrieved from KDE or wkde2d

- feature:

  Name of the feature being plotted

- cell_embeddings:

  Matrix with cell embeddings

- dim_names:

  Names of the dimensions from the cell embeddings

- shape:

  Geom shape

- size:

  Geom size

- text_size:

  text_size given from parent function

- legend.position:

  legend.position from parent function

- legend_title:

  String used as legend title

- pal:

  String specifying the color palette to use, can be one of hcl.pals

- raster:

  Rasterise plot

- ...:

  Further scale arguments passed to scale_color_viridis_c

## Value

A ggplot object

## Author

Jose Alquicira-Hernandez (modified)
