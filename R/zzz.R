.onLoad <- function(libname, pkgname) {
    if (getRversion() >= "2.15.1") {
        utils::globalVariables(c(
            ".", "Mean", "UMI", "all_zeros", "avg.exp", "avg_log2FC",
            "avg_log2FC_PB_DESeq2", "cluster", "condition", "freq", "gene",
            "group", "group_count", "label", "nCount_RNA", "nFeature_RNA",
            "orig.ident", "p_val", "p_val_PB_DESeq2", "p_val_adj", "stars",
            "p_val_adj_PB_DESeq2", "pct.exp", "prob", "proportion",
            "run_scvi", "scDblFinder_class", "significance", "heatmap",
            "split_bar_gsea", "useMyCol", "value", "variable", "xaxis",
            "pivot_longer", "heatmap_foldchange", "p.adjust", "group1", "group2"
            , ".y.", "n1","n2","statistic", "p.adj", "p.adj.signif",
            "setNames"))
    }
    invisible()
}
