#perform SEM Graphs
#' @author Mariano Ruz Jurado
#' @title SEM Graph with wilcox test on single cell level
#' @description Perform SEM-based graphs with Wilcox test on single-cell level for Seurat objects.
#' Calculates mean expression values and SEM for the selected feature, and visualizes them.
#' Performs pairwise Wilcox tests comparing conditions, with optional custom control condition and clustering.
#' Optionally returns a summary data frame, statistical test results, and the generated plot.
#' @param sce_object combined SCE object or Seurat
#' @param Feature name of the feature/gene
#' @param ListTest List for which conditions wilcoxon test will be performed, if NULL always CTRL group against everything
#' @param returnValues return data frames needed for the plot, containing df.melt, df.melt.sum, df.melt.orig and wilcoxstats
#' @param ctrl.condition set your ctrl condition, relevant if running with empty comparison List
#' @param group.by select the seurat object slot where your conditions can be found, default conditon
#' @param bar_colours colour vector
#' @param plotPvalue plot the non adjusted p-value without correcting for multiple tests
#' @param stat_pos_mod Defines the distance to the graphs of the statistic
#' @param step_mod Defines the distance between each statistics bracket
#' @param x_label_rotation Rotation of x-labels
#' @param log1p_nUMI If nUMIs should be log1p transformed
#' @param y_limits set limits for y-axis
#' @param wilcox_test perform wilcox test. boolean default TRUE
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#' @importFrom SeuratObject as.Seurat
#'
#' @return a ggplot or a list with plot and data frame
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' ListTest <- list()
#' ListTest[[1]] <- c("healthy", "disease")
#'
#' DO.BarplotWilcox(
#'   sce_object = sce_data,
#'   Feature = "NKG7",
#'   ListTest = ListTest,
#'   ctrl.condition = "healthy",
#'   group.by="condition"
#' )
#'
#' @export
DO.BarplotWilcox <- function(sce_object,
                             Feature,
                             ListTest=NULL,
                             returnValues=FALSE,
                             ctrl.condition=NULL,
                             group.by = "condition",
                             wilcox_test=TRUE,
                             bar_colours=NULL,
                             stat_pos_mod = 1.15,
                             step_mod=0.2,
                             x_label_rotation=45,
                             plotPvalue=FALSE,
                             y_limits = NULL,
                             log1p_nUMI=TRUE){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    sce_object <- as.Seurat(sce_object)
  }

  if (!(Feature %in% rownames(sce_object)) && !(Feature %in% names(sce_object@meta.data))) {
    stop("Feature not found in SCE Object!")
  }

  if (wilcox_test == TRUE) {
    rstat <- system.file(package = "rstatix") # Make sure package is installed
    ifelse(nzchar(rstat), "", stop("Install rstatix R package for wilcox statistic!"))
  }

  #SEM function defintion
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided sce_object, aswell as original identifier of samples
  df<-data.frame(condition=setNames(sce_object[[group.by]][,group.by], rownames(sce_object[[group.by]]))
                 ,orig.ident = sce_object$orig.ident)
  #get expression values for genes from individual cells, add to df
  #if (SeuV5==FALSE) {
  #  for(i in Feature){
  #    df[,i] <- expm1(sce_object@assays$RNA@data[i,])
  #
  #  }
  #}
  #For Seuratv5 where everything is a layer now
  #if (SeuV5==TRUE) {
  #rlang::warn("\nSeuV5 set to TRUE, if working with Seuratv4 or below change SeuV5 to FALSE", .frequency = "once", .frequency_id = "v5Mean")
  if (Feature %in% rownames(sce_object)) {
    df[,Feature] <- expm1(FetchData(sce_object, vars = Feature))
  }else{
    df[,Feature] <- FetchData(sce_object, vars = Feature)
  }
  #}

  # stat.df$condition <- factor(stat.df$condition)
  # stat.df$variable <- factor(stat.df$variable)

  # stat.df.test <- data.frame(Mean = sce_object[["RNA"]]@data["Figf",],condition = sce_object[[group.by]][,group.by])

  #melt results
  df.melt <- melt(df)
  #group results and summarize, also add/use SEM
  df.melt.sum <- df.melt %>%
    dplyr::group_by(condition, variable) %>%
    dplyr::summarise(Mean = mean(value))
  #second dataframe containing mean values for individual samples
  df.melt.orig <- df.melt %>%
    dplyr::group_by(condition, variable, orig.ident) %>%
    dplyr::summarise(Mean = mean(value))

  if (Feature %in% names(sce_object@meta.data)) {
    plot.title <- paste("Mean", Feature, sep = " ")
  }  else if (log1p_nUMI==TRUE) {
    df.melt.sum$Mean <- log1p(df.melt.sum$Mean)
    df.melt.orig$Mean <- log1p(df.melt.orig$Mean)
    plot.title <- "Mean log(nUMI)"
  } else{
    plot.title <- "Mean nUMI"
  }


  #create comparison list for wilcox, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list
    # and define ListTest comparing every other condition with that ctrl condition
    .logger("ListTest empty, comparing every sample with each other")
    conditions <- unique(sce_object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) {
      ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","ctrl","WT","Wt","wt"),collapse ="|")
                                        ,conditions)[1]]
    }

    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in seq_along(conditions)) {
      cndtn <- as.character(conditions[i])
      if (cndtn != ctrl.condition) {
        ListTest[[length(ListTest) + 1]] <- c(ctrl.condition, cndtn)
      }
    }
  }
  #delete Null values, created by count index also reorder for better p-value depiction
  ListTest <- ListTest[!vapply(ListTest, is.null, logical(1))]
  indices <- vapply(ListTest, function(x) match(x[2], df.melt.sum$condition), integer(1))
  ListTest <- ListTest[order(indices)]

  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$condition %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(sprintf(
          "Removing Test %s vs %s since both values are 0",
          elements[1], elements[2]
        ))
        }
    }
    return(lst_filtered)
  }

  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df.melt.sum)


  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE) {
    stat.test <- df.melt %>%
      ungroup() %>%
      rstatix::wilcox_test(value ~ condition, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(sce_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
  }
  #add SEM calculated over sample means
  df.melt.sum$SEM <- NA
  for (condition in df.melt.sum$condition) {
    df.melt.orig.con <- df.melt.orig[df.melt.orig$condition %in% condition,] # condition wise
    for (gene in Feature) {
      df.melt.orig.con.gen <- df.melt.orig.con[df.melt.orig.con$variable %in% gene,] #gene wise
      df.melt.sum[df.melt.sum$condition %in% condition & df.melt.sum$variable %in% gene,]$SEM <- SEM(df.melt.orig.con.gen$Mean)
    }
  }

  if (is.null(bar_colours)) {
    bar_colours <- rep(c("#1f77b4","#ea7e1eff","royalblue4","tomato2","darkgoldenrod","palegreen4","maroon","thistle3"),5)#20 colours set for more change number
  }

  if (x_label_rotation == 45) {
    hjust <- 1
  } else{hjust <- 0.5}

  #create barplot with significance
  p<-ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = Mean, ymax = Mean+SEM), width = 0.2)+
    geom_point(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1, position = "jitter")+
    #ordering, control always first
    scale_x_discrete(limits=c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
    #t-test, always against control, using means from orig sample identifier
    #stat_compare_means(data=stat.df,aes(x=condition, y=Mean, group=variable),comparisons = ListTest,method = "wilcox",size=3, label.y = max(df.melt.sum$Mean)*1.4)+
    facet_wrap(~variable, ncol = 9, scales = "free") +
    scale_fill_manual(values =  bar_colours
                      , name = "Condition")+
    labs(title = "", y = plot.title) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",angle = x_label_rotation,hjust = hjust, size = 16),
          axis.text.y = element_text(color = "black", size = 16),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5),
          axis.line = element_line(color = "black"),
          strip.text.x = element_text(size = 16, color = "black"),
          legend.position = "none")
  if (!is.null(y_limits)) {
    p <- p + ylim(y_limits)
  }
  if (wilcox_test == TRUE) {

    #Adjustments when ylim is changed manually
    y_pos_test <- max(df.melt.orig$Mean)*stat_pos_mod
    if (!is.null(y_limits) && y_pos_test > max(y_limits)) {
      y_pos_test <- max(y_limits)* stat_pos_mod - 0.1 * diff(y_limits)
    }
    if (plotPvalue==TRUE) {
      p <- p + stat_pvalue_manual(stat.test, label = "p = {p}", y.position = y_pos_test, step.increase = step_mod)
    }
    else{
      p <- p + stat_pvalue_manual(stat.test, label = "p = {p.adj}", y.position = y_pos_test, step.increase = step_mod)
    }
  }


  if (returnValues==TRUE) {
    returnList <- list(p, df.melt, df.melt.orig, df.melt.sum, stat.test)
    names(returnList) <- c("plot","df.melt", "df.melt.orig", "df.melt.sum", "stat.test")
    return(returnList)
  }
  return(p)
}

