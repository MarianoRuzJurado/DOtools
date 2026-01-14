# DO Heatmap of the mean expression of genes across a groups

Wrapper around heatmap.py, which generates a heatmap of showing the
average nUMI for a set of genes in different groups. Addiional an
argumnt can be made to show foldchanges between two conditions.
Differential gene expression analysis between the different groups can
be performed.

## Usage

``` r
DO.Heatmap(
  sce_object,
  features,
  assay_normalized = "RNA",
  group_by = "seurat_clusters",
  groups_order = NULL,
  value_plot = "expr",
  group_fc = "condition",
  group_fc_ident_1 = NULL,
  group_fc_ident_2 = NULL,
  clip_value = FALSE,
  max_fc = 5,
  z_score = NULL,
  path = NULL,
  filename = "Heatmap.svg",
  swap_axes = TRUE,
  cmap = "Reds",
  title = NULL,
  title_fontprop = NULL,
  clustering_method = "complete",
  clustering_metric = "euclidean",
  cluster_x_axis = FALSE,
  cluster_y_axis = FALSE,
  axs = NULL,
  figsize = c(5, 6),
  linewidth = 0.1,
  ticks_fontdict = NULL,
  xticks_rotation = NULL,
  yticks_rotation = NULL,
  vmin = 0,
  vcenter = NULL,
  vmax = NULL,
  legend_title = "LogMean(nUMI)\nin group",
  add_stats = TRUE,
  df_pvals = NULL,
  stats_x_size = NULL,
  square_x_size = NULL,
  test = "wilcox",
  pval_cutoff = 0.05,
  log2fc_cutoff = 0,
  only_pos = TRUE,
  square = TRUE,
  showP = TRUE,
  logcounts = TRUE
)
```

## Arguments

- sce_object:

  SCE object or Seurat with meta.data

- features:

  gene names or continuous value in meta data

- assay_normalized:

  Assay with raw counts

- group_by:

  meta data column name with categorical values

- groups_order:

  order for the categories in the group_by

- value_plot:

  plotted values correspond to expression values or foldchanges

- group_fc:

  if foldchanges specified than the groups must be specified that will
  be compared

- group_fc_ident_1:

  Defines the first group in the test

- group_fc_ident_2:

  Defines the second group in the test

- clip_value:

  Clips the colourscale to the 99th percentile, useful if one gene is
  driving the colourscale

- max_fc:

  Clips super high foldchanges to this value, so changes can still be
  appreciated

- z_score:

  apply z-score transformation, "group" or "var"

- path:

  path to save the plot

- filename:

  name of the file

- swap_axes:

  whether to swap the axes or not

- cmap:

  color map

- title:

  title for the main plot

- title_fontprop:

  font properties for the title (e.g., 'weight' and 'size')

- clustering_method:

  clustering method to use when hierarchically clustering the x and
  y-axis

- clustering_metric:

  metric to use when hierarchically clustering the x- and y-axis

- cluster_x_axis:

  hierarchically clustering the x-axis

- cluster_y_axis:

  hierarchically clustering the y-axis

- axs:

  matplotlib axis

- figsize:

  figure size

- linewidth:

  line width for the border of cells

- ticks_fontdict:

  font properties for the x and y ticks (e.g., 'weight' and 'size')

- xticks_rotation:

  rotation of the x-ticks

- yticks_rotation:

  rotations of the y-ticks

- vmin:

  minimum value

- vcenter:

  center value

- vmax:

  maximum value

- legend_title:

  title for the color bar

- add_stats:

  add statistical annotation, will add a square with an '\*' in the
  center if the expression is significantly different in a group with
  respect to the others

- df_pvals:

  dataframe with the p-values, should be gene x group or group x gene in
  case of swap_axes is False

- stats_x_size:

  size of the asterisk

- square_x_size:

  size and thickness of the square percentual, vector

- test:

  test to use for test for significance

- pval_cutoff:

  cutoff for the p-value

- log2fc_cutoff:

  minimum cutoff for the log2FC

- only_pos:

  if set to TRUE, only use positive genes in the condition

- square:

  whether to make the cell square or not

- showP:

  if set to false return a dictionary with the axis

- logcounts:

  whether the input is logcounts or not

## Value

Depending on `showP`, returns the plot if set to `TRUE` or a dictionary
with the axes.

## Author

Mariano Ruz Jurado & David Rodriguez Morales

## Examples

``` r
sce_data <-
  readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

DO.Heatmap(
sce_object = sce_data,
assay_normalized = "RNA",
group_by="seurat_clusters",
features = rownames(sce_data)[1:10],
z_score = NULL,
path = NULL,
filename = "Heatmap.svg",
swap_axes = TRUE,
cmap = "Reds",
title = NULL,
title_fontprop = NULL,
clustering_method = "complete",
clustering_metric = "euclidean",
cluster_x_axis = FALSE,
cluster_y_axis = FALSE,
axs = NULL,
figsize = c(5, 6),
linewidth = 0.1,
ticks_fontdict = NULL,
xticks_rotation = 45,
yticks_rotation = NULL,
vmin = 0.0,
vcenter = NULL,
vmax = NULL,
legend_title = "LogMean(nUMI)\nin group",
add_stats = TRUE,
df_pvals = NULL,
stats_x_size = NULL,
square_x_size = NULL,
test = "wilcox",
pval_cutoff = 0.05,
log2fc_cutoff = 0,
only_pos = TRUE,
square = TRUE,
showP = FALSE,
logcounts = TRUE
)
#> Calculating cluster 1
#> Warning: No features pass logfc.threshold threshold; returning empty data.frame
#> Calculating cluster 2
#> Calculating cluster 3
#> Calculating cluster 4
#> Calculating cluster 5
#> Calculating cluster 6
#> Calculating cluster 7
#> Warning: No features pass logfc.threshold threshold; returning empty data.frame
#> Calculating cluster 8

```
