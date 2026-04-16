# Polished UMAP function using Dimplot or FeaturePlot function from Seurat
#' @author Mariano Ruz Jurado
#' @title DO.UMAP
#' @description Creates a polished UMAP plot using Seurat's DimPlot or
#' FeaturePlot functions. In addition a density plot can be made in a similar
#' way to nebulosa R package. It allows customization of colors, labels,
#' and other plot elements for better visualisation. The function handles both
#' cluster-based visualisations and gene-based visualisations in a UMAP plot.
#' Ideal for refining UMAP outputs with added flexibility and enhanced
#' presentation.
#' @param sce_object The seurat or SCE object
#' @param FeaturePlot Is it going to be a FeaturePlot?
#' @param DensityPlot Is it going to be a DensityPlot?
#' @param features features for Featureplot
#' @param reduction reduction to use
#' @param group.by grouping of plot in DImplot and defines in featureplot the
#' labels
#' @param dims Dimensions to plot, must be a two-length numeric vector
#' specifying x- and y-dimensions description
#' @param layer Layer to use for DensityPlot, default data
#' @param umap_colors what colors to use for UMAP, specify as vector
#' @param text_size Size of text
#' @param label label the clusters on the plot by group.by column
#' @param order Boolean determining whether to plot cells in order of
#' expression.
#' @param plot.title title for UMAP
#' @param legend.position specify legend position
#' @param method Kernel density estimation method, can be "ks" or "wkde"
#' @param ... Further arguments passed to DimPlot, FeaturePlot or DensityPlot
#' functions
#' @return Plot with Refined colors and axes
#'
#' @import Seurat
#' @import ggplot2
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.UMAP(
#'     sce_object = sce_data,
#'     group.by = "seurat_clusters"
#' )
#'
#' DO.UMAP(
#'     sce_object = sce_data,
#'     FeaturePlot = TRUE,
#'     features = c("BAG2", "CD74")
#' )
#'
#' DO.UMAP(
#'     sce_object = sce_data,
#'     DensityPlot = TRUE,
#'     features = c("CD74")
#' )
#'
#' @export
DO.UMAP <- function(
    sce_object,
    features = NULL,
    group.by = "seurat_clusters",
    FeaturePlot = FALSE,
    DensityPlot = FALSE,
    reduction = NULL,
    dims = c(1, 2),
    layer = NULL,
    umap_colors = NULL,
    text_size = 14,
    label = TRUE,
    order = TRUE,
    plot.title = TRUE,
    legend.position = "none",
    method = c("ks", "wkde"),
    ...
) {
    # support for single cell experiment objects
    if (methods::is(sce_object, "SingleCellExperiment")) {
        sce_object <- .suppressAllWarnings(as.Seurat(sce_object))
    }

    # Dimplot
    if (isFALSE(FeaturePlot) && isFALSE(DensityPlot)) {
        if (is.null(umap_colors)) {
            umap_colors <- rep(
                c(
                    "#1f77b4",
                    "#ff7f0e",
                    "#2ca02c",
                    "tomato2",
                    "#9467bd",
                    "chocolate3",
                    "#e377c2",
                    "#ffbb78",
                    "#bcbd22",
                    "#17becf",
                    "darkgoldenrod2",
                    "#aec7e8",
                    "#98df8a",
                    "#ff9896",
                    "#c5b0d5",
                    "#c49c94",
                    "#f7b6d2",
                    "#c7c7c7",
                    "#dbdb8d",
                    "#9edae5",
                    "sandybrown",
                    "moccasin",
                    "lightsteelblue",
                    "darkorchid",
                    "salmon2",
                    "forestgreen",
                    "bisque"
                ),
                5
            )
        }

        p <- DimPlot(
            sce_object,
            reduction = reduction,
            group.by = group.by,
            cols = umap_colors,
            dims = dims,
            ...
        ) +
            labs(x = "UMAP1", y = "UMAP2") +
            theme(
                plot.title = element_blank(),
                # text = element_text(face = "bold",size = 20),
                axis.title.x = element_text(
                    size = text_size,
                    family = "Helvetica"
                ),
                axis.title.y = element_text(
                    size = text_size,
                    family = "Helvetica"
                ),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = legend.position,
                legend.text = element_text(face = "bold")
            )

        if (label == TRUE) {
            p <- LabelClusters(p, id = group.by, fontface = "bold", box = FALSE)
        }
        return(p)
    }

    if (FeaturePlot & DensityPlot) {
        stop("Please define only FeaturePlot or DensityPlot as TRUE!")
    }

    # FeaturePlot
    if (FeaturePlot == TRUE) {
        if (is.null(features)) {
            stop("Please provide any gene names if using FeaturePlot=TRUE.")
        }

        if (is.null(umap_colors)) {
            umap_colors <- c("lightgrey", "red2")
        }

        Idents(sce_object) <- group.by
        p <- .suppressDeprecationWarnings(FeaturePlot(sce_object,
            reduction = reduction,
            features = features,
            cols = umap_colors,
            label = label,
            order = order,
            dims = dims,
            ...
        ) &
            labs(x = "UMAP1", y = "UMAP2") &
            theme(
                axis.title.x = element_text(size = 14, family = "Helvetica"),
                axis.title.y = element_text(size = 14, family = "Helvetica"),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = legend.position,
                legend.text = element_text(face = "bold")
            ))

        if (plot.title == FALSE) {
            p <- p & theme(plot.title = element_blank())
        }

        return(p)
    }

    # DensityPlot after Nebulosa package, minor adjustments
    if (DensityPlot == TRUE) {
        if (is.null(features)) {
            stop("Please provide any gene names if using DensityPlot=TRUE.")
        }

        # Validate existence of reduction
        if (!is.null(reduction)) {
            if (!reduction %in% Reductions(sce_object)) {
                stop("No reduction named '", reduction, "' found in sce_object")
            }
        }

        # Set up default reduction
        reduction_list <- Reductions(sce_object)
        if (!length(reduction_list)) {
            stop("No reduction has been computed!")
        }

        if (is.null(umap_colors)) {
            umap_colors <- "Reds"
        }

        if (is.null(reduction)) {
            # change to search for umap/tsne as default

            reduction <- grepv("umap|tsne", reduction_list, ignore.case = TRUE)

            if (!length(reduction)) {
                # afterwards take last
                reduction <- reduction_list[length(reduction_list)]
            }
        }

        cell_embeddings <- as.data.frame(Embeddings(sce_object[[reduction]]))

        # check if all dimensions where found
        i <- dims %in% seq_len(ncol(cell_embeddings))
        if (!all(i)) {
            missing_dims <- dims[which(!i)]
            stop(
                "Dimension(s) ",
                missing_dims,
                " not present in",
                reduction,
                "\n"
            )
        }

        cell_embeddings <- cell_embeddings[, dims]

        # Layer set to data
        if (is.null(layer)) {
            layer <- "data"
        }

        # Match kde method in ks
        method <- match.arg(method)

        # Extract metadata
        metadata <- sce_object[[]]

        # Determine type of feature and extract
        if (all(features %in% colnames(metadata))) {
            vars <- metadata[, features, drop = FALSE]
        } else {
            exp_data <- FetchData(sce_object, vars = features, layer = layer)

            # Extract data for input features
            i <- colnames(exp_data) %in% features

            # Test existence of feature in gene expression data
            j <- !features %in% colnames(exp_data)
            if (any(j)) {
                stop(
                    "'",
                    paste(features[j], collapse = ", "),
                    "' feature(s) not present in meta.data or expression data"
                )
            }
            vars <- exp_data[, i, drop = FALSE]
            vars <- vars[, features, drop = FALSE]
        }

        ## plotting
        dim_names <- colnames(cell_embeddings)
        # for multiple features
        if (ncol(vars) > 1) {
            # simple way that works
            plot_list <- list()
            for (feature in colnames(vars)) {
                vars_sub <- vars[, feature, drop = FALSE]
                res <- .calculate_density(
                    vars_sub[, 1],
                    cell_embeddings,
                    method
                )
                plot_list[[feature]] <- .plot_density_(
                    z = res$z,
                    dens = res$dens,
                    feature,
                    pal = umap_colors,
                    cell_embeddings,
                    dim_names,
                    text_size = text_size,
                    legend.position = legend.position,
                    ...
                )
            }

            return(plot_list)
        } else {
            res <- .calculate_density(vars[, 1], cell_embeddings, method)
            p <- .plot_density_(
                z = res$z,
                dens = res$dens,
                features,
                pal = umap_colors,
                cell_embeddings,
                dim_names,
                text_size = text_size,
                legend.position = legend.position,
                ...
            )

            return(p)
        }
    }
}


#' @title Estimate weighted kernel density
#' @author Jose Alquicira-Hernandez (modified)
#' @param w Vector with weights for each observation
#' @param x Matrix with dimensions where to calculate the density from. Only
#' the first two dimensions will be used
#' @param method Kernel density estimation method:
#' \itemize{
#' \item \code{ks}: Computes density using the \code{kde} function from the
#'  \code{ks} package.
#' \item \code{wkde}: Computes density using a modified version of the
#'  \code{kde2d} function from the \code{MASS}
#' package to allow weights. Bandwidth selection from the \code{ks} package
#'  is used instead.
#' }
#' @param adjust Numeric value to adjust to bandwidth. Default: 1. Not available
#'  for \code{ks} method
#' @param map Whether to map densities to individual observations
#' @return If \code{map} is \code{TRUE}, a vector with corresponding densities
#'  for each observation is returned. Otherwise,
#' a list with the density estimates from the selected method is returned.
#' @importFrom ks kde hpi
#' @keywords internal
.calculate_density <- function(w, x, method, adjust = 1, map = TRUE) {
    if (method == "ks") {
        dens <- kde(x[, c(1, 2)],
            w = w / sum(w) * length(w)
        )
    } else if (method == "wkde") {
        dens <- .wkde2d(
            x = x[, 1],
            y = x[, 2],
            w = w / sum(w) * length(w),
            adjust = adjust
        )
    }

    if (map) {
        if (method == "ks") {
            ix <- findInterval(x[, 1], dens$eval.points[[1]])
            iy <- findInterval(x[, 2], dens$eval.points[[2]])
            ii <- cbind(ix, iy)
            z <- dens$estimate[ii]
        } else if (method == "wkde") {
            ix <- findInterval(x[, 1], dens$x)
            iy <- findInterval(x[, 2], dens$y)
            ii <- cbind(ix, iy)
            z <- dens$z[ii]
        }
    } else {
        return(dens)
    }
    return(list(z = z, dens = dens))
}


#' @title Plot density estimates
#' @author Jose Alquicira-Hernandez (modified)
#' @param z Vector with density values for each cells
#' @param dens density grid retrieved from KDE or wkde2d
#' @param feature Name of the feature being plotted
#' @param cell_embeddings Matrix with cell embeddings
#' @param dim_names Names of the dimensions from the cell embeddings
#' @param shape Geom shape
#' @param size Geom size
#' @param text_size text_size given from parent function
#' @param legend.position legend.position from parent function
#' @param legend_title String used as legend title
#' @param pal String specifying the color palette to use, can be one of hcl.pals
#' @param raster Rasterise plot
#' @param ... Further scale arguments passed to scale_color_viridis_c
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab ggtitle labs
#' guide_legend theme element_text element_line element_rect element_blank
#' scale_color_viridis_c scale_color_gradientn
#' @importFrom ggrastr rasterise
#' @keywords internal
.plot_density_ <- function(
    z,
    dens,
    feature,
    cell_embeddings,
    dim_names,
    shape = 16,
    size = 1,
    text_size = 14,
    density_quantile_threshold = 0.5, # for now
    n_bands = 256, # matches the number of colours in palette
    legend.position = "right",
    legend_title = "Density",
    pal = "Reds",
    raster = FALSE,
    rev_colours_density = FALSE,
    ...
) {
    # build grid with density -> has no 0 values since gaussian KDE
    if (!is.null(dens$eval.points)) {
        grid_df <- expand.grid(
            x = dens$eval.points[[1]],
            y = dens$eval.points[[2]]
        )
        grid_df$z <- as.vector(dens$estimate)
    } else {
        grid_df <- expand.grid(
            x = dens$x,
            y = dens$y
        )
        grid_df$z <- as.vector(dens$z)
    }

    palette_colors <- rev(hcl.colors(256, palette = pal))

    if (rev_colours_density == TRUE) {
        palette_colors <- rev(palette_colors)
    }

    # Keep only densities above a quantile threshold
    # (removes low-density clutter)
    threshold <- quantile(grid_df$z, density_quantile_threshold, na.rm = TRUE)
    grid_df$z[grid_df$z < threshold] <- NA

    # breaks for filled contours
    z_range <- range(grid_df$z, na.rm = TRUE)
    if (diff(z_range) == 0) z_range <- c(0, 1)
    breaks <- seq(z_range[1], z_range[2], length.out = n_bands + 1)

    # Map breaks to colours
    fill_colors <- colorRampPalette(palette_colors)(length(breaks) - 1)

    p <- ggplot(data.frame(cell_embeddings, feature = z)) +
        aes(x = .data[[dim_names[1]]], y = .data[[dim_names[2]]]) +
        geom_point(
            aes(color = .data[["feature"]]),
            shape = shape,
            size = size,
            alpha = 1
        ) +


        # contour lines using the same colour scale
        # geom_contour(
        #     data = grid_df,
        #     aes(x = x, y = y, z = z, colour = after_stat(level)),
        #     inherit.aes = FALSE,
        #     linewidth = 0.25,
        #     alpha = 0.5,
        #     bins = 20
        # ) +

        geom_contour_filled(
            data = grid_df,
            aes(x = x, y = y, z = z, fill = after_stat(level)),
            inherit.aes = FALSE,
            breaks = breaks,
            alpha = 0.25, # low opacity so points remain visible
            show.legend = FALSE # avoid duplicate legend
        ) +

        # alternative
        # geom_raster(
        #     data = grid_df,
        #     aes(x = x, y = y, fill = z),
        #     inherit.aes = FALSE,
        #     # alpha = 0.3
        # )+


        # stat_density2d( # calculates density newly based on point density
        #     aes(fill = after_stat(level)),
        #     geom = "polygon",
        #     alpha = 0.3,
        #     # colour = "white",
        #     lineend = "round",
        #     h = c(0.5, 0.5)) +

        ggtitle(feature) +
        labs(
            title = feature,
            x = "UMAP1",
            y = "UMAP2",
            color = legend_title
        ) +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(linewidth = 0.5),
            strip.background = element_rect(color = "black", fill = "#ffe5cc"),
            axis.title.x = element_text(size = text_size, family = "Helvetica"),
            axis.title.y = element_text(size = text_size, family = "Helvetica"),
            plot.title = element_text(
                size = text_size,
                face = "bold",
                family = "Helvetica",
                hjust = 0.5
            ),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = legend.position,
            legend.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold")
        ) +
        scale_colour_gradientn(
            colours = palette_colors,
            breaks = c(min(z), max(z)),
            labels = c("Min", "Max")
        ) +
        scale_fill_manual(values = fill_colors)
    p
    # p <-
    # p + scale_color_continuous(palette = rev(hcl.colors(256, palette = pal)))
    # p <-
    # p + scale_fill_continuous(palette = rev(hcl.colors(256, palette = pal)))

    if (raster) {
        # Make sure R ggrastr package is installed
        rt <- system.file(package = "ggrastr")
        ifelse(nzchar(rt), "", stop(
            "Install ggrastr R package for rastor option in ggplot2 usage ",
        ))
        ggrastr::rasterise(p, dpi = 300)
    } else {
        return(p)
    }
}

#' @title Weighted 2D kernel density estimation
#' @author Jose Alquicira-Hernandez (modified)
#' @param x Dimension 1
#' @param y Dimension 2
#' @param w Weight variable
#' @param h vector of bandwidths for x and y directions.
#' Defaults to normal reference bandwidth (ks::hpi).
#' A scalar value will be taken to apply to both directions.
#' @param adjust Bandwidth adjustment
#' @param n Number of grid points in each direction. Can be scalar or a
#' length-2 integer vector.
#' @param lims The limits of the rectangle covered by the grid as
#' c(xl, xu, yl, yu).
#' @return A list of three components.
#' \itemize{
#' \item \code{x, y} The x and y coordinates of the grid points, vectors of
#' length n.
#' \item \code{z} A matrix of the weighted estimated density:
#' rows correspond to the value of x, columns to the value of y.
#' }
#' @importFrom Matrix Matrix
#' @importFrom stats dnorm
#' @importFrom methods is
#' @keywords internal
.wkde2d <- function(
    x, y, w, h, adjust = 1, n = 100,
    lims = c(range(x), range(y))
) {
    # Validate values and dimensions
    nx <- length(x)
    if (!all(all(nx == length(y)), all(nx == length(w)))) {
        stop("data vectors must be the same length")
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
        stop("missing or infinite values in the data are not allowed")
    }
    if (any(!is.finite(lims))) {
        stop("only finite values are allowed in 'lims'")
    }

    h <- c(
        ks::hpi(x),
        ks::hpi(y)
    )
    h <- h * adjust

    # Get grid
    gx <- seq.int(lims[1L], lims[2L], length.out = n)
    gy <- seq.int(lims[3L], lims[4L], length.out = n)

    # weight
    ax <- outer(gx, x, "-") / h[1L]
    ay <- outer(gy, y, "-") / h[2L]

    w <- Matrix::Matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE)

    z <- Matrix::tcrossprod(dnorm(ax) * w, dnorm(ay) * w) /
        (sum(w) * h[1L] * h[2L])

    dens <- list(x = gx, y = gy, z = z)
    dens
}
