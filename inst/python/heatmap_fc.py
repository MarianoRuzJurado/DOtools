import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from collections.abc import Iterable
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.colors
from matplotlib.colors import ListedColormap
from matplotlib.cm import ScalarMappable
from scipy.stats import zscore
from typing import Literal, Union, Dict
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import scipy as sp
import warnings
from numba import njit
import platform

FAST_ARRAY_UTILS = False # Keep an eye on this, just in case the calculation time is very slow

#py_version = platform.python_version()  # Check the Python Version
#if int(py_version.split(".")[1]) > 10:
#    try:
#        import fast_array_utils
#        FAST_ARRAY_UTILS = True
#    except ImportError as e:
#        logger.warn("Python > 3.10 but fast_array_utils not available, consider installing it to speed up"
#                    "analysis")

warnings.filterwarnings("ignore")

def check_missing(adata: ad.AnnData, features: str | list = None, groups: str | list = None,
                  variables: str | list = None) -> None:
    """Check for missing features or columns in the observations from an AnnData Object.

    :param adata: AnnData Object.
    :param features: features to check for.
    :param groups: column names in the observations to check for.
    :param variables: column names in the variables to check for.
    :return: Returns None. Will raise an assertion if any feature or column name is missing.
    """

    if features:
        features = iterase_input(features)
        missing = [g for g in features if g not in adata.var_names]

        # features could be in .obs
        missing_x2 = []
        if len(missing) > 0:
            missing_x2 = [g for g in features if g not in adata.obs.columns]

        if len(missing_x2) > len(missing):
            assert len(missing) == 0, f"{missing} missing in the AnnData Object"
        else:
            assert len(missing_x2) == 0, f"{missing_x2} missing in the AnnData Object"

    if groups:
        groups = iterase_input(groups)
        missing = [g for g in groups if g not in adata.obs.columns]
        assert len(missing) == 0, f"{missing} missing in the AnnData Object"
    if variables:
        variables = iterase_input(variables)
        missing = [g for g in variables if g not in adata.var.columns]
        assert len(missing) == 0, f"{missing} missing in the AnnData Object"

    return None

def iterase_input(data: str | float | int | Iterable | None) -> list:
    """Convert input to list.

    :param data: string or iterable (list, tuple, index, etc.)
    :return: Returns a list.
    """
    if data is None:
        return []
    elif isinstance(data, str):
        return [data]
    elif isinstance(data, float):
        return [data]
    elif isinstance(data, int):
        return [data]
    elif isinstance(data, list):
        return data
    elif isinstance(data, Iterable):
        return list(data)
    else:
        raise Exception("Input is not a string or iterable object")

def convert_path(path: Path | str) -> Path:
    """Convert to Path format if string is provided.

    :param path: string or Path variable.
    :return: Path
    """
    if not isinstance(path, Path):
        return Path(path)
    else:
        return path

def sanitize_anndata(adata: ad.AnnData) -> None:
    """Transform string metadata to categorical.

    :param adata: AnnData
    :return None
    """
    adata._sanitize()
    return None

def make_grid_spec(
    ax_or_figsize,
    *,
    nrows: int,
    ncols: int,
    wspace: float = None,
    hspace: float = None,
    width_ratios: float | list = None,
    height_ratios: float | list = None,
):
    # Taken from Scanpy
    kw = {"wspace": wspace, "hspace": hspace, "width_ratios": width_ratios, "height_ratios": height_ratios}

    if isinstance(ax_or_figsize, tuple):
        fig = plt.figure(figsize=ax_or_figsize)
        return fig, gridspec.GridSpec(nrows, ncols, **kw)
    else:
        ax = ax_or_figsize
        ax.axis("off")
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])
        return ax.figure, ax.get_subplotspec().subgridspec(nrows, ncols, **kw)


def check_colornorm(vmin=None, vmax=None, vcenter=None, norm=None):
    from matplotlib.colors import Normalize

    try:
        from matplotlib.colors import TwoSlopeNorm as DivNorm
    except ImportError:
        # matplotlib<3.2
        from matplotlib.colors import DivergingNorm as DivNorm

    if norm is not None:
        if (vmin is not None) or (vmax is not None) or (vcenter is not None):
            raise ValueError("Passing both norm and vmin/vmax/vcenter is not allowed.")
    else:
        if vcenter is not None:
            norm = DivNorm(vmin=vmin, vmax=vmax, vcenter=vcenter)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)

    return norm


def square_color(rgba: list) -> str:
    """Determine if the background is dark or clear and return black or white.

    :param rgba: list with rgba values
    :return: black or white
    """
    r, g, b = rgba[:3]  # ignore alpha
    # Convert from 0 to 1 float to 0–255 int
    r, g, b = [int(c * 255) for c in (r, g, b)]
    # Use brightness heuristic
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return "black" if brightness > 128 else "white"


def small_squares(ax: plt.Axes, pos: list, color: list, size: float = 1, linewidth: float = 0.8,
                  zorder: int = 20) -> None:
    """Add small squares.

    :param ax: matplotlib axis
    :param pos: list of positions
    :param color: list of colors
    :param size:  size of the square
    :param linewidth: linewith of the square
    :param zorder: location of the square
    :return: None
    """
    for idx, xy in enumerate(pos):
        x, y = xy
        margin = (1 - size) / 2
        rect = patches.Rectangle(
            (y + margin, x + margin),
            size,
            size,
            linewidth=linewidth,
            edgecolor=color[idx],
            facecolor="none",
            zorder=zorder,
        )
        if zorder == 0:
            rect.set_alpha(0)  # Hide square if they should be in the back, for the dotplot
        ax.add_patch(rect)
    return None

def _expm1_anndata(adata: ad.AnnData) -> None:
    """Apply expm1 transformation for the X data.

    :param adata: annotated data matrix
    :return: None, changes are inplace
    """
    import scipy as sp

    if sp.sparse.issparse(adata.X):
        adata.X = adata.X.copy()
        adata.X.data = np.expm1(adata.X.data)
    else:
        adata.X = np.expm1(adata.X)


def expr(
    adata: ad.AnnData,
    features: str | list,
    groups: str | list | None = None,
    out_format: Literal["long", "wide"] = "long",
    layer: str | None = None,
) -> pd.DataFrame:
    """Extract the expression of features.

    This function extract the expression from an AnnData object and returns a DataFrame. If layer
    is not specified the expression in `X` will be extracted. Additionally, metadata from `obs` can be added
    to the dataframe.

    :param adata: Annotated data matrix.
    :param groups: Metadata column in `obs` to include in the DataFrame.
    :param features: Name of the features in `var_names` to extract the expression of.
    :param out_format: Format of the dataframe. The `wide` format will generate a DataFrame with shape n_obs x n_vars,
                      while the `long` format will generate an unpivot version.
    :param layer: Layer in the AnnData object to extract the expression from. If set to `None` the expression in
                  `X` will be used.

    Returns
    -------
    Returns a `DataFrame`. If `out_format` is set to `wide`, the index will be the cell barcodes and the column names
    will be set to the gene names. If `groups` are specified, extra columns will be present. If `out_format` is set to
    `long`, the following fields are included:

    `genes`
        Contains the gene names.
    `expr`
        Contains the expression values extracted.

    Example
    -------
    >>> import dotools_py as do
    >>> adata = do.dt.example_10x_processed()
    >>> df = do.get.expr(adata, "CD4", "annotation")
    >>> df.head(5)
      annotation genes  expr
    0    B_cells   CD4   0.0
    1         NK   CD4   0.0
    2    T_cells   CD4   0.0
    3    T_cells   CD4   0.0
    4    T_cells   CD4   0.0
    >>> df = do.get.expr(adata, "CD4", "annotation", out_format="wide")
    >>> df.head(5)
                                   CD4 annotation
    CAAAGAATCAGATTGC-1-batch2  0.0    B_cells
    AGCTTCCCAGTCAACT-1-batch1  0.0         NK
    GAGAGGTTCCCTCTAG-1-batch1  0.0    T_cells
    CTAACTTCAGATCATC-1-batch1  0.0    T_cells
    CATGGTACAAACGGCA-1-batch1  0.0    T_cells

    """
    sanitize_anndata(adata)
    features = iterase_input(features)
    groups = iterase_input(groups)
    assert len(features) != 0, "No features provided"
    assert out_format == "wide" or out_format == "long", f'{out_format} not recognize, try "long" or "wide"'

    # Remove features not present and warn
    check_missing(adata, features=features, groups=groups)

    # Set-up configuration
    adata = adata[:, features]  # Retain only the specified features
    if layer is not None:
        adata.X = adata.layers[layer].copy()  # Select the specified layer

    # Extract expression
    table_expr = adata.to_df().copy()

    # Add Metadata
    if groups is not None:
        table_expr[groups] = adata.obs[groups]
    if out_format == "long":
        table_expr = pd.melt(table_expr, id_vars=groups, var_name="genes", value_name="expr")
    free_memory()
    return table_expr


def mean_expr(
    adata: ad.AnnData,
    group_by: str | list,
    features: list | str | None = None,
    out_format: Literal["long", "wide"] = "long",
    layer: str | None = None,
    logcounts: bool = True,
) -> pd.DataFrame:
    """Calculate the average expression in an AnnData objects for features.

    This function calculates the average expression of a set of features grouping by one
    or several categories. Assume log-normalized counts. If logcounts is set to True, the
    log10 transformation is undone for the mean expression calculation. The reported mean
    expression is log-transformed.

    :param adata: Annotated data matrix.
    :param group_by: Metadata columns in `obs` to group by.
    :param features: List of features in `var_name` to use. If not set, it will be calculated over all the genes.
    :param out_format: Format of the Dataframe returned. This can be wide or long format.
    :param layer: Layer of the AnnData to use. If not set use `X`.
    :param logcounts: Set to `True` if the input is in log space.

    Returns
    -------
    Returns a `DataFrame` with the mean expression in log1p transformation. If `out_format` is set to `wide`, the index
    will be set to the gene names and the column names will be set to the groups. If `out_format` is set to `long`,
    the following fields are included:

    `gene`
        Contains the gene names.
    `groupN`
        Contains the groups (For each metadata column a new column will be added).
    `expr`
        Contains the mean expression values after log1p transformation.

    Example
    -------
    >>> import dotools_py as do
    >>> adata = do.dt.example_10x_processed()
    >>> df = do.get.mean_expr(adata, "annotation")
    >>> df.head(5)
             gene   group0      expr
    0  ATP2A1-AS1  B_cells  0.000000
    1      STK17A  B_cells  1.453713
    2    C19orf18  B_cells  0.000000
    3        TPP2  B_cells  0.126846
    4       MFSD1  B_cells  0.053630
    >>> df = do.get.mean_expr(adata, "annotation", out_format="wide")
    >>> df.head(5)
        group0   B_cells  Monocytes        NK   T_cells       pDC
    gene
    A4GALT  0.222505   0.000000  0.000000  0.000000  0.000000
    AAK1    0.000000   0.364976  1.126293  1.143016  0.128019
    ABAT    0.182251   0.146378  0.047404  0.045826  0.158761
    ABCB4   0.062785   0.000000  0.000000  0.000000  0.000000
    ABCB9   0.000000   0.000000  0.027683  0.057814  0.000000

    """

    sanitize_anndata(adata)
    features = list(adata.var_names) if len(iterase_input(features)) == 0 else iterase_input(features)
    group_by = iterase_input(group_by)
    check_missing(adata, features=features, groups=group_by)
    assert out_format == "wide" or out_format == "long", f'{out_format} not recognize, try "long" or "wide"'

    # Set-up configuration
    adata = adata[:, features]
    if layer is not None:
        adata.X = adata.layers[layer].copy()

    data = adata.copy()

    if logcounts:
        _expm1_anndata(data)

    # Group data by the specified values
    group_obs = adata.obs.groupby(group_by, as_index=False)

    # Compute AverageExpression
    main_df = pd.DataFrame([])
    for group_name, df in group_obs:
        if FAST_ARRAY_UTILS:
            current_mean = fast_array_utils.stats.mean(data[df.index].X, axis=0)
            current_mean = np.log1p(current_mean) if logcounts else current_mean
            df_tmp = pd.DataFrame(current_mean, columns=["expr"])
        else:
            if logcounts:
                df_tmp = np.log1p(
                    pd.DataFrame(data[df.index].X.mean(axis=0).T, columns=["expr"])
                )  # Mean expr per gene in groupN
            else:
                df_tmp = pd.DataFrame(data[df.index].X.mean(axis=0).T, columns=["expr"])

        df_tmp["gene"] = adata[df.index].var_names  # Update with Gene names
        group_name = iterase_input(group_name)
        for idx, name in enumerate(group_name):
            # df_tmp["group" + str(idx)] = str(name).replace("-", "_")  # Update with metadata
            df_tmp[group_by[idx]] = name
        main_df = pd.concat([main_df, df_tmp], axis=0)
    main_df["expr"] = pd.to_numeric(main_df["expr"])  # Convert to numeric values

    # Move expr column to last position
    expr_col = main_df.pop("expr")
    main_df["expr"] = expr_col

    # Change to wide format
    if out_format == "wide":
        main_df = pd.pivot_table(main_df, index="gene", columns=group_by, values="expr")
        if len(group_by) > 1:
            main_df.columns = main_df.columns.map("_".join)
    free_memory()
    return main_df


@njit(parallel=True)
def _get_log2fc(group: np.ndarray, ref: np.ndarray, psc=1e-9):
    return np.log2((np.expm1(group) + psc) / (np.expm1(ref) + psc))

def get_log2fc(
    adata: ad.AnnData,
    group_by: str,
    reference: str,
    groups: str | list | None = None,
    features: str | list | None = None,
    layer: str | None = None,
) -> pd.DataFrame:
    """Calculate the log2foldchanges for a set of groups.

    :param adata: Annotated data matrix.
    :param group_by: Column in `obs` to group by.
    :param reference: Reference condition to use for the calculation.
    :param groups: Alternative condititons to use. If `None`, all the condititons will be used.
    :param features: Features to use for calculating the log2foldchanges. If set to `None` all features will be used.
    :param layer: Layer in the AnnData to use for the calculation.

    Returns
    -------
    Returns a DataFrame with the log2-foldchanges. One column will be added for each condition in `groups`

    Example
    -------
    >>> import dotools_py as do
    >>> adata = do.dt.example_10x_processed()
    >>> df = do.get.log2fc(adata, group_by="condition", reference="healthy")
    >>> df.head(5)
            genes  log2fc_disease
    0  ATP2A1-AS1       26.073313
    1      STK17A       -0.429677
    2    C19orf18        0.775196
    3        TPP2      -22.599501
    4       MFSD1       -1.669137

    """

    # Get the data
    features = iterase_input(features)
    features = list(adata.var_names) if len(iterase_input(features)) == 0 else iterase_input(features)
    groups = list(adata.obs[group_by].unique()) if len(iterase_input(groups)) == 0 else iterase_input(groups)
    if reference in groups:
        groups.remove(reference)

    df_mean = mean_expr(adata, group_by=group_by, features=features, out_format="wide", layer=layer)

    logfoldchanges = pd.DataFrame([], index=list(df_mean.index))
    for group in groups:
        # Speed up with numba
        foldchanges = _get_log2fc(group=df_mean[group].to_numpy(), ref=df_mean[reference].to_numpy())
        logfoldchanges["log2fc_" + group] = foldchanges
    logfoldchanges.reset_index(inplace=True)
    logfoldchanges.rename(columns={"index": "genes"}, inplace=True)
    return logfoldchanges


def free_memory() -> None:
    """Garbage collector.

    :return:
    """
    import ctypes
    import gc

    gc.collect()

    system = platform.system()

    if system == "Linux":
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    else:
        pass
    return None

def get_hex_colormaps(colormap: str) -> list:
    """Get a list with Hexa IDs for a colormap.

    :param colormap: colormap name.
    :return: list with Hexa IDs.

    Example
    -------
    >>> import dotools_py as do
    >>> hex_list = do.utility.get_hex_colormaps("Reds")
    >>> hex_list[:5]
    ['#fff5f0', '#fff4ef', '#fff4ee', '#fff3ed', '#fff2ec']

    """
    cmap = plt.get_cmap(colormap)
    return [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]

def heatmap_foldchange(
    # Data
    adata: ad.AnnData,
    group_by: str | list,
    features: str | list,
    condition_key: str,
    reference: str,
    groups_order: list = None,
    conditions_order: list = None,
    layer: str = None,

    # Figure parameters
    figsize: tuple = (5, 6),
    ax: plt.Axes | None = None,
    swap_axes: bool = True,
    title: str = None,
    title_fontproperties: Dict[Literal["size", "weight"], str | int] = None,
    palette: str = "RdBu_r",
    palette_conditions: str | dict = "tab10",
    ticks_fontproperties: Dict[Literal["size", "weight"], str | int] = None,
    xticks_rotation: int = None,
    yticks_rotation: int = None,
    vmin: float = None,
    vcenter: float = None,
    vmax: float = None,
    colorbar_legend_title: str = "Log2FC",
    groups_legend_title: str = "Comparison",
    group_legend_ncols: int = 1,

    # IO
    path: str | Path = None,
    filename: str = "Heatmap.svg",
    show: bool = True,

    # Statistics
    add_stats: bool = False,
    test: Literal["wilcoxon", "t-test"] = "wilcoxon",
    correction_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg",
    df_pvals: pd.DataFrame = None,
    stats_x_size: float = None,
    square_x_size: dict = None,
    pval_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.0,

    # Fx specific
    linewidth: float = 0.1,
    color_axis_ratio=0.15,
    **kargs,
) -> dict | None:
    """Heatmap of the log2-foldchange of genes across a groups between two conditions.

    Generate a heatmap of showing the log2-foldchange for a set of genes in different groups between different
    conditions. Differential gene expression analysis between the different conditions can be performed.

    :param adata: Annotated data matrix.
    :param group_by:  Name of a categorical column in `adata.obs` to groupby.
    :param features: A valid feature in `adata.var_names` or column in `adata.obs` with continuous values.
    :param condition_key: Name of a categorical column in `adata.obs` to compare to compute the log2foldchanges
    :param reference: Category in `adata.obs[condition_key]` to use as the reference to compute the log2foldchanges
    :param groups_order: Order for the categories in `adata.obs[group_by]`.
    :param conditions_order: Order for the categories in `adata.obs[condition_key]`
    :param path: Path to the folder to save the figure.
    :param filename: Name of file to use when saving the figure.
    :param layer: Name of the AnnData object layer that wants to be plotted. By default, `adata.X` is plotted.
                 If layer is set to a valid layer name, then the layer is plotted.
    :param swap_axes: Whether to swap the axes or not.
    :param palette: String denoting matplotlib colormap.
    :param palette_conditions: String denoting matplotlib colormap for the comparisons.
    :param title: Title for the figure.
    :param title_fontproperties: Dictionary which should contain 'size' and 'weight' to define the fontsize and
                                fontweight of the title of the figure.
    :param ax: Matplotlib axes to use for plotting. If not set, a new figure will be generated.
    :param figsize: Figure size, the format is (width, height).
    :param linewidth: Linewidth for the border of cells.
    :param ticks_fontproperties: Dictionary which should contain 'size' and 'weight' to define the fontsize and fontweight of the font of the x/y-axis.
    :param xticks_rotation: Rotation of the x-ticks.
    :param yticks_rotation: Rotations of the y-ticks.
    :param vmin: The value representing the lower limit of the color scale.
    :param vcenter: The value representing the center of the color scale.
    :param vmax: The value representing the upper limit of the color scale.
    :param colorbar_legend_title: Title for the colorbar.
    :param groups_legend_title: Title for the comparison legend.
    :param group_legend_ncols: Number of columns for the comparison legend.
    :param add_stats: Add statistical annotation. Will add a square with an '*' in the center if the expression is significantly different in a group with respect to the others.
    :param df_pvals: Dataframe with the pvals. Should be `gene x group` or `group x gene` in case of swap_axes is `False`.
    :param stats_x_size: Scaling factor to control the size of the asterisk.
    :param square_x_size: Size and thickness of the square.
    :param test: Name of the method to test for significance.
    :param correction_method: Correction method for multiple testing.
    :param pval_cutoff: Cutoff for the p-value.
    :param log2fc_cutoff: Minimum cutoff for the log2FC.
    :param show: If set to `False`, returns a dictionary with the matplotlib axes.
    :param color_axis_ratio: Ratio of the axis reserved for the colors denoting the comparisons.
    :param kargs: Additional arguments pass to `sns.heatmap <https://seaborn.pydata.org/generated/seaborn.heatmap.html>`_.
    :return: Depending on ``show``, returns the plot if set to `True` or a dictionary with the axes.

    Example
    -------

    .. plot::
        :context: close-figs

        import dotools_py as do
        adata = do.dt.example_10x_processed()
        do.pl.heatmap_foldchange(adata, 'annotation', ['CD4', 'CD79A'], "condition", "healthy", add_stats=True, figsize=(4, 6))

    """
    import scanpy as sc  # Lazy load
    from matplotlib.colorbar import Colorbar

    # Checks
    assert reference is not None, "Provide reference condition"
    assert condition_key is not None, "Provide a column to compute log-foldchanges and statistics"
    features = iterase_input(features)
    assert all(item in list(adata.var_names) for item in features), \
        "log-foldchanges can only be computed for features in adata.var_names"

    sanitize_anndata(adata)
    check_missing(adata, features=features)

    # Compute the log-foldchanges,
    # Format gene|class x group_by (Genes duplicated for each class)
    df = []
    for g in adata.obs[group_by].unique():
        adata_subset = adata[adata.obs[group_by] == g]
        df_tmp = get_log2fc(
            adata_subset, group_by=condition_key, features=features, layer=layer, reference=reference
        )
        df_tmp = df_tmp.melt(id_vars="genes", value_name="log2fc", var_name="class")
        df_tmp["class"] = df_tmp["class"].str.replace("log2fc_", "")
        df_tmp["group_by"] = g
        df.append(df_tmp)
    df = pd.concat(df)
    df = df.pivot(index=["genes", "class"], columns="group_by", values="log2fc")
    # from dotools_py.utility import get_hex_colormaps
    if isinstance(palette_conditions, str):
        if len(adata.obs[condition_key].cat.categories) > len(get_hex_colormaps(palette_conditions)):
            print(f"There are {len(adata.obs[condition_key].cat.categories)} conditions but the colormap has"
                        f"{get_hex_colormaps(palette_conditions)} shades.")
        class_dictionaries = dict(zip(adata.obs[condition_key].cat.categories, get_hex_colormaps(palette_conditions)))
    elif isinstance(palette_conditions, dict):
        if all(item in list(adata.obs[condition_key].unique()) for item in palette_conditions):
            class_dictionaries = palette_conditions
        else:
            raise ValueError("Missing conditions in palette_conditions")
    else:
        raise NotImplementedError("palette_conditions is not a str or dict")

    # Set the order
    if conditions_order is not None:
        new_gene_class_order = [(g, c) for c in conditions_order for g in features]
    else:
        new_gene_class_order = [(g, c) for c in np.unique(df.index.get_level_values("class")) for g in features]
        conditions_order = np.unique(df.index.get_level_values("class"))
    groups_order = groups_order if groups_order is not None else list(adata.obs[group_by].unique())

    if isinstance(df.index, pd.MultiIndex):
        df = df.reindex(index=new_gene_class_order, columns=groups_order)
    else:
        df = df.reindex(columns=new_gene_class_order, index=groups_order)

    # Layout
    if swap_axes:
        df = df.T

    # Compute Statistics
    annot_pvals = None
    if add_stats:
        if df_pvals is None:
            if all(item in list(adata.var_names) for item in features):
                table_filt = []
                for g in adata.obs[group_by].unique():
                    adata_subset = adata[adata.obs[group_by] == g]
                    try:
                        rank_genes_groups(
                            adata_subset, groupby=condition_key, method=test, reference=reference, tie_correct=True,
                            corr_method=correction_method
                        )
                    except ValueError as e:
                        print(f"Failed to compute stats for {g}: {e}")
                        continue
                    table = sc.get.rank_genes_groups_df(
                        adata_subset, group=None, pval_cutoff=pval_cutoff, log2fc_min=log2fc_cutoff
                    )
                    table[group_by] = g
                    if "group" not in table.columns:
                        group = list(adata_subset.obs[condition_key].unique())
                        group.remove(reference)
                        table["group"] = group[0]

                    table_filt.append(table[table["names"].isin(features)])
                table_filt = pd.concat(table_filt)
            else:
                raise Exception('Not Implemented, all features needs to be in adata.var_names')

            # Dataframe with gene|class x groups with the pvals
            df_pvals = pd.DataFrame([], index=df.index, columns=df.columns)
            for idx, row in table_filt.iterrows():
                if isinstance(df.index, pd.MultiIndex):
                    if row["group"] in list(df.index.get_level_values("class")):
                        df_pvals.loc[(row["names"], row["group"]), row[group_by]] = row["pvals_adj"]
                else:
                    if row["group"] in list(df.columns.get_level_values("class")):
                        df_pvals.loc[row[group_by], (row["names"], row["group"])] = row["pvals_adj"]

            df_pvals[df_pvals.isna()] = 1
        else:  # TODO Check with R version

            df_pvals = (
                df_pvals
                .pivot(
                    index=group_by,
                    columns=["genes", condition_key],
                    values="value"
                )
                .sort_index()
            )

            if list(df.index)[0] in list(df_pvals.index):
                pass
            else:
                df_pvals = df_pvals.T
        # Replace pvals < 0.05 with an X
        annot_pvals = df_pvals.applymap(lambda x: "*" if x < pval_cutoff else "")

    # -------------------------------------------- Arguments for the layout --------------------------------------------

    # # # # # # # #
    # Figure Layout
    # # # # # # # #
    width, height = figsize
    legends_width_spacer = 0.7 / width
    mainplot_width = width - (1.5 + 0)

    min_figure_height = max([0.35, height])
    cbar_legend_height = min_figure_height * 0.08
    sig_legend = min_figure_height * 0.15
    foldchange_legend = min_figure_height * 0.1
    spacer_height = min_figure_height * 0.27

    height_ratios = [
        height - sig_legend - cbar_legend_height - spacer_height - foldchange_legend,
        foldchange_legend,
        sig_legend,
        spacer_height,
        cbar_legend_height,
    ]

    # # # # # # # # #
    # Text Properties
    # # # # # # # # #  ticks_fontproperties = {} if ticks_fontproperties is None else ticks_fontproperties
    title_fontproperties = {} if title_fontproperties is None else title_fontproperties
    ticks_fontproperties = {} if ticks_fontproperties is None else ticks_fontproperties
    ticks_fontproperties = {
        "weight": ticks_fontproperties.get("weight", "bold"), "size": ticks_fontproperties.get("size", 13)
    }
    title_fontprop = {
        "weight": title_fontproperties.get("weight", "bold"), "size": title_fontproperties.get("size", 15)
    }

    # # # # # #
    # Colorbar
    # # # # # #
    vmin =  df.min().min()  if vmin is None else vmin
    vmax = df.max().max() if vmax is None else vmax
    vcenter = 0.0 if vcenter is None else vcenter


    colormap = plt.get_cmap(palette)
    normalize = check_colornorm(vmin=vmin, vmax=vmax, vcenter=vcenter)
    mappable = ScalarMappable(norm=normalize, cmap=colormap)

    #mean_flat = df.T.values.flatten()
    #color = colormap(normalize(mean_flat))
    #color = [square_color(c) for c in color]

    # # # # # # #
    # Statistics
    # # # # # # #
    square_x_size = {} if square_x_size is None else square_x_size
    square_x_size = {
        "width": square_x_size.get("weight", 1), "size": square_x_size.get("size", 0.8)
    }
    stats_x_size = min(width / df.shape[1], height / df.shape[1]) * 10 if stats_x_size is None else min(
        width / df.shape[1], height / df.shape[1]) * stats_x_size

    # # # # # # # #
    # Colors axis
    # # # # # # # #
    height_ratios_main, width_ratios_main, nrows_main, ncols_main = [height], [mainplot_width], 1, 1
    if "genes" in df.columns.names:
        nrows_main, ncols_main = 2, 1
        height_ratios_main = [
            color_axis_ratio,
            height - color_axis_ratio
        ]
        pos_groups, pos_main = 0, 1
    else:
        nrows_main, ncols_main = 1, 2
        width_ratios_main = [
            mainplot_width - color_axis_ratio,
            color_axis_ratio
        ]
        pos_groups, pos_main = 1, 0

    # # # #
    # Axis
    # # # #
    return_ax_dict = {}

    # -------------------------------------------- Generate the Figure -------------------------------------------------
    fig, gs = make_grid_spec(
        ax or (width, height), nrows=1, ncols=2, wspace=legends_width_spacer, width_ratios=[mainplot_width, 1.5]
    )

    # Create Main Axis
    main_ax = fig.add_subplot(gs[0])
    fig, main_ax_gs = make_grid_spec(
        main_ax, nrows=nrows_main, ncols=ncols_main, hspace=0.01, wspace=0.01, height_ratios=height_ratios_main,
        width_ratios=width_ratios_main
    )
    main_ax = fig.add_subplot(main_ax_gs[pos_main])

    # Create Legend Axis
    legend_ax = fig.add_subplot(gs[1])
    fig, legend_gs = make_grid_spec(legend_ax, nrows=len(height_ratios), ncols=1, height_ratios=height_ratios)
    color_legend_ax = fig.add_subplot(legend_gs[4])

    # Create Colors Axis
    groups_ax = fig.add_subplot(main_ax_gs[pos_groups])

    # Create Significance axis
    if add_stats:
        sig_ax = fig.add_subplot(legend_gs[3])

    # Main Plot
    # Correction to remove class from the dataframe (gene|class x groups)
    df_copy = df.copy()
    df_copy = df_copy.droplevel("class") if "genes" in df.index.names else df_copy.T.droplevel("class").T

    hm = sns.heatmap(
        data=df_copy, cmap=palette, ax=main_ax, linewidths=linewidth, cbar=False, annot=annot_pvals, fmt="s",
        square=False, annot_kws=
        {"color": "black", "size": stats_x_size, "ha": "center", "va": "center", "fontfamily": 'DejaVu Sans Mono'},
        vmax=vmax, vmin=vmin, center=vcenter,
        **kargs,
    )
    if isinstance(df.index, pd.MultiIndex):
        classes = df.reset_index()["class"].astype(
            pd.CategoricalDtype(categories=conditions_order, ordered=True)
        )
        tmp = pd.DataFrame(classes.cat.codes)
    else:
        classes = df.T.reset_index()["class"].astype(
            pd.CategoricalDtype(categories=conditions_order, ordered=True)
        )
        tmp = pd.DataFrame(classes.cat.codes).T

    unique_color_vector = [class_dictionaries[c] for c in classes.cat.categories]
    cmap = ListedColormap(unique_color_vector)
    sns.heatmap(
        tmp, cmap=cmap, ax=groups_ax, xticklabels=False, yticklabels=False, cbar=False
    )

    # Add colors group Legend
    groups_ax_legend = fig.add_subplot(legend_gs[1])
    handles = []

    for lab in conditions_order:
        txt = lab + " Vs " + reference
        c = class_dictionaries[lab]
        handles.append(
            mlines.Line2D(
                [0], [0], marker=".", color=c, lw=0, label=txt, markerfacecolor=c, markeredgecolor=None,
                markersize=18
            )
        )
    groups_ax_legend.legend(
        handles=handles, frameon=False, loc="center", ncols=group_legend_ncols, prop={"size": "small", "weight": "bold"}, title=groups_legend_title,
        title_fontproperties={"size": "small", "weight": "bold"}, borderaxespad=0.2, bbox_transform=groups_ax_legend.transAxes, bbox_to_anchor=(0.5, 0.5),

    )
    groups_ax_legend.axis("off")  # Hide axes for clean display
    return_ax_dict["color_group_ax"] = groups_ax
    return_ax_dict["legend_group_ax"] = groups_ax_legend


    # Add Colorbar Legend
    Colorbar(color_legend_ax, mappable=mappable, orientation="horizontal")
    color_legend_ax.xaxis.set_tick_params(labelsize="small")
    color_legend_ax.set_title(colorbar_legend_title, fontsize="small", fontweight="bold")
    return_ax_dict["legend_ax"] = color_legend_ax

    # Significance Legend
    if add_stats:
        x, y = 0, 0.5
        sig_ax.scatter(x, y, s=500, facecolors="none", edgecolors="black", marker="s")
        sig_ax.text(x, y, "*", fontsize=18, ha="center", va="center", color="black", fontfamily='DejaVu Sans Mono')
        sig_ax.text(x + 0.03, y, "FDR < 0.05", fontsize=12, va="center", fontweight="bold")
        sig_ax.set_xlim(x - 0.02, x + 0.1)
        sig_ax.set_title("Significance", fontsize="small", fontweight="bold")
        plt.gca().set_aspect("equal")
        sig_ax.axis("off")  # Hide axes for clean display
        return_ax_dict["signifiance_ax"] = sig_ax

    # Modify layout from main plot
    hm.spines[["top", "right", "bottom", "left"]].set_visible(True)
    hm.set_xlabel("")
    hm.set_ylabel("")

    rotation_props_x, rotation_props_y = {"rotation": None}, {"rotation": None}
    rotation_props_x = (
        {"rotation": xticks_rotation, "va": "top", "ha": "right"} if xticks_rotation is not None else rotation_props_x
    )
    rotation_props_y = (
        {"rotation": yticks_rotation, "va": "top", "ha": "right"} if yticks_rotation is not None else rotation_props_y
    )
    hm.set_xticklabels(hm.get_xticklabels(), fontdict={"weight": ticks_fontproperties["weight"], "size": ticks_fontproperties["size"]}, **rotation_props_x)
    hm.set_yticklabels(hm.get_yticklabels(), fontdict={"weight": ticks_fontproperties["weight"], "size": ticks_fontproperties["size"]}, **rotation_props_y)
    groups_ax.set_title(title, **title_fontprop)
    return_ax_dict["mainplot_ax"] = hm

    # Add Square around the Xs
    if add_stats:
        df_x = pd.DataFrame([], index=df.index, columns=df.columns)
        df_x[df_x.isna()] = "black"
        df_x = df.map(lambda x: square_color(colormap(normalize(x))))
        pos_rows, pos_cols = np.where(df_pvals < 0.05)
        pos = list(zip(pos_rows, pos_cols, strict=False))
        colors = [df_x.iloc[row, col] for row, col in pos]

        small_squares(
            hm,
            color=colors,
            pos=pos,
            size=square_x_size["size"],
            linewidth=square_x_size["width"],
        )

        # Now set colors manually on each annotation text base on the background
        for text, color in zip(hm.texts, df_x.values.flatten(), strict=False):
            text.set_color(color)

    if path is not None:
        plt.savefig(convert_path(path) / filename, bbox_inches="tight")
    if show:
        return plt.show()
    else:
        return return_ax_dict
