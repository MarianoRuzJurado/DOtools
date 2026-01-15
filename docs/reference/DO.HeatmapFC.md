# DO Heatmap of the mean expression of genes across a groups

Wrapper around heatmap_foldchange, which generates a heatmap of showing
the foldchange for a set of gene expressions between specified groups.
Differential gene expression analysis between the different groups can
be performed.

## Usage

``` r
DO.HeatmapFC(
  sce_object,
  features,
  reference = NULL,
  assay_normalized = "RNA",
  group_by = "seurat_clusters",
  condition_key = "condition",
  groups_order = NULL,
  conditions_order = NULL,
  layer = NULL,
  figsize = c(5, 6),
  ax = NULL,
  swap_axes = TRUE,
  title = NULL,
  title_fontproperties = list(size = NULL, weight = NULL),
  palette = "RdBu_r",
  palette_conditions = "tab10",
  ticks_fontproperties = list(size = NULL, weight = NULL),
  xticks_rotation = NULL,
  yticks_rotation = NULL,
  vmin = NULL,
  vcenter = NULL,
  vmax = NULL,
  colorbar_legend_title = "Log2FC",
  groups_legend_title = "Comparison",
  group_legend_ncols = 1,
  path = NULL,
  filename = "Heatmap.svg",
  showP = TRUE,
  add_stats = TRUE,
  test = c("wilcox"),
  correction_method = c("bonferroni"),
  df_pvals = NULL,
  stats_x_size = NULL,
  square_x_size = NULL,
  pval_cutoff = 0.05,
  log2fc_cutoff = 0,
  linewidth = 0.1,
  color_axis_ratio = 0.15
)
```

## Arguments

- sce_object:

  A SingleCellExperiment or Seurat object containing expression data and
  metadata.

- features:

  Character vector of gene names or metadata column names to be
  visualized.

- reference:

  Reference condition used for fold-change calculation.

- assay_normalized:

  Name of the assay containing normalized expression values (default:
  "RNA").

- group_by:

  Metadata column defining the primary grouping variable (e.g.
  clusters).

- condition_key:

  Metadata column defining the condition or comparison variable.

- groups_order:

  Optional character vector specifying the order of groups in
  `group_by`.

- conditions_order:

  Optional character vector specifying the order of conditions.

- layer:

  Optional layer name to extract expression values from.

- figsize:

  Numeric vector of length two specifying figure width and height.

- ax:

  Optional matplotlib axis object (for Python backend usage).

- swap_axes:

  Logical; whether to swap x- and y-axes.

- title:

  Optional title for the heatmap.

- title_fontproperties:

  Named list specifying font properties for the title (e.g. size,
  weight).

- palette:

  Color palette used for the heatmap.

- palette_conditions:

  Color palette used for condition annotations.

- ticks_fontproperties:

  Named list specifying font properties for axis tick labels.

- xticks_rotation:

  Rotation angle for x-axis tick labels.

- yticks_rotation:

  Rotation angle for y-axis tick labels.

- vmin:

  Minimum value for the color scale.

- vcenter:

  Center value for the color scale.

- vmax:

  Maximum value for the color scale.

- colorbar_legend_title:

  Title for the color bar.

- groups_legend_title:

  Title for the group legend.

- group_legend_ncols:

  Number of columns in the group legend.

- path:

  Optional path to save the output figure.

- filename:

  Name of the output file.

- showP:

  Logical; whether to display the plot.

- add_stats:

  Logical; whether to add statistical annotations.

- test:

  Statistical test to use (currently "wilcox").

- correction_method:

  Multiple-testing correction method (currently "bonferroni").

- df_pvals:

  Optional data frame containing precomputed p-values (groups x genes or
  genes x groups depending on axis orientation).

- stats_x_size:

  Size of the statistical annotation symbol.

- square_x_size:

  Size of the square annotation.

- pval_cutoff:

  P-value significance threshold.

- log2fc_cutoff:

  Minimum absolute log2 fold-change cutoff.

- linewidth:

  Line width of heatmap cell borders.

- color_axis_ratio:

  Relative size of the color bar axis.

## Value

Depending on `showP`, returns the plot if set to `TRUE` or a dictionary
with the axes.

## Author

Mariano Ruz Jurado & David Rodriguez Morales

## Examples

``` r
sce_data <-
  readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

DO.HeatmapFC(
sce_object = sce_data,
features = c("HES4", "ISG15", "TNFRSF18", "TNFRSF4", "MMP23B"),
reference = NULL,
assay_normalized = "RNA",
group_by = "seurat_clusters",
condition_key = "condition",
groups_order = c("1", "2", "3", "4", "5", "6", "7", "8"),
conditions_order = NULL,
layer = NULL,

# Figure parameters
figsize = c(5, 6),
ax = NULL,
swap_axes = TRUE,
title = NULL,
title_fontproperties = list(size = NULL, weight = NULL),
palette = "RdBu_r",
palette_conditions = "tab10",
ticks_fontproperties = list(size = NULL, weight = NULL),
xticks_rotation = 45,
yticks_rotation = NULL,
vmin = NULL,
vcenter = NULL,
vmax = NULL,
colorbar_legend_title = "Log2FC",
groups_legend_title = "Comparison",
group_legend_ncols = 1,

# IO
path = NULL,
filename = "Heatmap.svg",
show = FALSE,

#   # Statistics
add_stats = TRUE,
test = c("wilcox"),
correction_method = c("bonferroni"),
df_pvals = NULL,
stats_x_size = NULL,
square_x_size = NULL,
pval_cutoff = 0.05,
log2fc_cutoff = 0.0,

# Fx specific
linewidth = 0.1,
color_axis_ratio = 0.15
)
#> 2026-01-14 16:45:34 - reference is set to NULL, using: healthy
```
