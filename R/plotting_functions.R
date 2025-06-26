#perform SEM Graphs
#' @author Mariano Ruz Jurado
#' @title SEM Graph with t test on cluster level
#' @description Perform SEM-based graphs with t-test on cluster level for Seurat objects.
#' Calculates mean expression values and SEM for selected features, and visualizes them.
#' Performs pairwise t-tests comparing conditions, with optional custom control condition and clustering.
#' Optionally returns a summary data frame.
#' @param Seu_object Combined Seu_object
#' @param Features Vector containing featurenames
#' @param ListTest List for which conditions t-test will be performed, if NULL always against provided CTRL
#' @param returnValues return df.melt.sum data frame containing means and SEM for the set group
#' @param ctrl.condition set your ctrl condition, relevant if running with empty comparison List
#' @param group.by select the seurat object slot where your conditions can be found, default conditon
#' @param bar_colours colour vector
#' @param stat_pos_mod Defines the distance to the graphs of the statistic
#' @param step_mod Defines the distance between each statistics bracket
#' @param x_label_rotation Rotation of x-labels
#' @param log1p_nUMI If nUMIs should be log1p transformed
#' @param y_limits set limits for y-axis
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#'
#' @return a ggplot or a dataframe
#'
#' @examples
#' \dontrun{
#'
#' (optional)
#' ListTest <- list()
#' ListTest[[1]] <- c("CTRL", "CONDITION")
#'
#' DO.BarplotClustert(
#'   Seu_object = Seurat,
#'   Features = "CDH5",
#'   ListTest = ListTest,
#'   ctrl.condition = "CTRL",
#'   group.by="condition"
#' )
#' }
#'
#' @export
DO.BarplotClustert <- function(Seu_object,
                               Feature,
                               ListTest=NULL,
                               returnValues=FALSE,
                               ctrl.condition=NULL,
                               group.by = "condition",
                               returnPlot=TRUE,
                               bar_colours=NULL,
                               stat_pos_mod = 1.15,
                               step_mod=0.2,
                               x_label_rotation=45,
                               y_limits = NULL,
                               log1p_nUMI=T){

  if (!(Feature %in% rownames(Seu_object)) && !(Feature %in% names(Seu_object@meta.data))) {
    stop("Feature not found in Seurat Object!")
  }

  #SEM function defintion
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided Seu_object, aswell as original identifier of samples
  df<-data.frame(condition=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]]))
                 ,orig.ident = Seu_object$orig.ident)
  #get expression values for genes from individual cells, add to df
  #if (SeuV5==F) {
  #  for(i in Feature){
  #    df[,i] <- expm1(Seu_object@assays$RNA@data[i,])
  #
  #  }
  #}
  #For Seuratv5 where everything is a layer now
  #if (SeuV5==T) {
  #rlang::warn("\nSeuV5 set to TRUE, if working with Seuratv4 or below change SeuV5 to FALSE", .frequency = "once", .frequency_id = "v5Mean")
  if (Feature %in% rownames(Seu_object)) {
    df[,Feature] <- expm1(FetchData(Seu_object, vars = Feature))
  }else{
    df[,Feature] <- FetchData(Seu_object, vars = Feature)
  }
  #}

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


  if (Feature %in% names(Seu_object@meta.data)) {
    plot.title <- paste("Mean", Feature, sep = " ")
  }  else if (log1p_nUMI==T) {
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
    cat("ListTest empty, comparing every sample with each other\n")
    conditions <- unique(Seu_object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) {
      ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","ctrl","WT","Wt","wt"),collapse ="|")
                                        ,conditions)[1]]
    }

    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in 1:length(conditions)) {
      cndtn <- as.character(conditions[i])
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }

  #delete Null values, created by count index also reorder for betetr p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  indices <- sapply(ListTest, function(x) match(x[2], df.melt.sum$condition))
  ListTest <- ListTest[order(indices)]

  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$condition %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
      }
    }
    return(lst_filtered)
  }

  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df.melt.sum)


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


  #Adjustments when ylim is changed manually
  y_pos_test <- max(df.melt.orig$Mean)*stat_pos_mod
  if (!is.null(y_limits) && y_pos_test > max(y_limits)) {
    y_pos_test <- max(y_limits)* stat_pos_mod - 0.1 * diff(y_limits)
  }


  #create barplot with significance
  p<-ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = Mean, ymax = Mean+SEM), width = 0.2)+
    geom_point(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1, position = "jitter")+
    #ordering, control always first
    scale_x_discrete(limits=c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
    #t-test, always against control, using means from orig sample identifier
    stat_compare_means(data=df.melt.orig, comparisons = ListTest, method = "t.test", size=3, y.position = y_pos_test, step.increase = step_mod)+
    facet_wrap(~variable, ncol = 9, scales = "free") +
    scale_fill_manual(values = bar_colours #20 colours set for more change number
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
  # print(p)
  if (!is.null(y_limits)) {
    p = p + ylim(y_limits)
  }

  if (returnValues==TRUE) {
    return(df.melt.sum)
  }
  if (returnPlot==TRUE) {
    return(p)
  }
}


#perform SEM Graphs
#' @author Mariano Ruz Jurado
#' @title SEM Graph with wilcox test on single cell level
#' @description Perform SEM-based graphs with Wilcox test on single-cell level for Seurat objects.
#' Calculates mean expression values and SEM for the selected feature, and visualizes them.
#' Performs pairwise Wilcox tests comparing conditions, with optional custom control condition and clustering.
#' Optionally returns a summary data frame, statistical test results, and the generated plot.
#' @param Seu_object combined Seurat object
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
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#'
#' @return a ggplot or a list with plot and data frame
#'
#' @examples
#' \dontrun{
#'
#' (optional)
#' ListTest <- list()
#' ListTest[[1]] <- c("CTRL", "CONDITION")
#'
#' DO.BarplotWilcox(
#'   Seu_object = Seurat,
#'   Feature = "CDH5",
#'   ListTest = ListTest,
#'   ctrl.condition = "CTRL",
#'   group.by="condition"
#' )
#' }
#'
#' @export
DO.BarplotWilcox <- function(Seu_object,
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
                             log1p_nUMI=T){

  if (!(Feature %in% rownames(Seu_object)) && !(Feature %in% names(Seu_object@meta.data))) {
    stop("Feature not found in Seurat Object!")
  }

  if (wilcox_test == T) {
    rstat <- system.file(package = "rstatix") # Make sure package is installed
    ifelse(nzchar(rstat), "", stop("Install rstatix R package for wilcox statistic!"))
  }

  #SEM function defintion
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided Seu_object, aswell as original identifier of samples
  df<-data.frame(condition=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]]))
                 ,orig.ident = Seu_object$orig.ident)
  #get expression values for genes from individual cells, add to df
  #if (SeuV5==F) {
  #  for(i in Feature){
  #    df[,i] <- expm1(Seu_object@assays$RNA@data[i,])
  #
  #  }
  #}
  #For Seuratv5 where everything is a layer now
  #if (SeuV5==T) {
  #rlang::warn("\nSeuV5 set to TRUE, if working with Seuratv4 or below change SeuV5 to FALSE", .frequency = "once", .frequency_id = "v5Mean")
  if (Feature %in% rownames(Seu_object)) {
    df[,Feature] <- expm1(FetchData(Seu_object, vars = Feature))
  }else{
    df[,Feature] <- FetchData(Seu_object, vars = Feature)
  }
  #}

  # stat.df$condition <- factor(stat.df$condition)
  # stat.df$variable <- factor(stat.df$variable)

  # stat.df.test <- data.frame(Mean = Seu_object[["RNA"]]@data["Figf",],condition = Seu_object[[group.by]][,group.by])

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

  if (Feature %in% names(Seu_object@meta.data)) {
    plot.title <- paste("Mean", Feature, sep = " ")
  }  else if (log1p_nUMI==T) {
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
    cat("ListTest empty, comparing every sample with each other\n")
    conditions <- unique(Seu_object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) {
      ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","ctrl","WT","Wt","wt"),collapse ="|")
                                        ,conditions)[1]]
    }

    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in 1:length(conditions)) {
      cndtn <- as.character(conditions[i])
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }

  #delete Null values, created by count index also reorder for betetr p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  indices <- sapply(ListTest, function(x) match(x[2], df.melt.sum$condition))
  ListTest <- ListTest[order(indices)]

  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$condition %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
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
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(Seu_object)))
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
    p = p + ylim(y_limits)
  }
  if (wilcox_test == TRUE) {

    #Adjustments when ylim is changed manually
    y_pos_test <- max(df.melt.orig$Mean)*stat_pos_mod
    if (!is.null(y_limits) && y_pos_test > max(y_limits)) {
      y_pos_test <- max(y_limits)* stat_pos_mod - 0.1 * diff(y_limits)
    }
    if (plotPvalue==TRUE) {
      p = p + stat_pvalue_manual(stat.test, label = "p = {p}", y.position = y_pos_test, step.increase = step_mod)
    }
    else{
      p = p + stat_pvalue_manual(stat.test, label = "p = {p.adj}", y.position = y_pos_test, step.increase = step_mod)
    }
  }


  if (returnValues==TRUE) {
    returnList <- list(p, df.melt, df.melt.orig, df.melt.sum, stat.test)
    names(returnList) <- c("plot","df.melt", "df.melt.orig", "df.melt.sum", "stat.test")
    return(returnList)
  }
  return(p)
}




#perform better Violins
#' @author Mariano Ruz Jurado
#' @title Violin Graph with wilcox test on single cell level
#' @description Creates a violin plot to compare gene expression across different conditions or groups within a Seurat object.
#' It incorporates Wilcoxon rank-sum tests to evaluate statistical differences between conditions.
#' The plot can be customized with options for data transformation, jitter display, and significance annotations.
#' The function also supports multiple conditions and allows for visualization of statistical results from wilcoxon test.
#' @param Seu_object combined Seurat object
#' @param Feature name of the feature
#' @param ListTest List for which conditions wilcox will be performed, if NULL always CTRL group against everything
#' @param returnValues return df.melt.sum data frame containing means and SEM for the set group
#' @param ctrl.condition set your ctrl condition, relevant if running with empty comparison List
#' @param group.by select the seurat Seu_object slot where your conditions can be found, default conditon
#' @param group.by.2 relevant for multiple group testing, e.g. for each cell type the test between each of them in two conditions provided
#' @param geom_jitter_args vector for dots visualisation in vlnplot: size, width, alpha value
#' @param vector_colours specify a minimum number of colours as you have entries in your condition, default 2
#' @param wilcox_test Bolean if TRUE a bonferoni wilcoxon test will be carried out between ctrl.condition and the rest
#' @param stat_pos_mod value for modifiyng statistics height
#' @param SeuV5 Seuratv5 object? (TRUE or FALSE)
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#'
#' @return a ggplot or a list used data frames
#'
#' @examples
#' \dontrun{
#'
#' (optional)
#' ListTest <- list()
#' ListTest[[1]] <- c("CTRL", "CONDITION")
#'
#' DO.VlnPlot(
#'   Seu_object = Seurat,
#'   SeuV5=T,
#'   Feature = "CDH5",
#'   ListTest = ListTest,
#'   ctrl.condition = "CTRL",
#'   group.by="condition"
#' )
#' }
#'
#' @export
DO.VlnPlot <- function(Seu_object,
                               SeuV5=T,
                               Feature,
                               ListTest=NULL,
                               returnValues=FALSE,
                               ctrl.condition=NULL,
                               group.by = "condition",
                               group.by.2 = NULL,
                               geom_jitter_args = c(0.20, 0.25, 0.25),
                               geom_jitter_args_group_by2 = c(0.1, 0.1, 1),
                               vector_colors = c("#1f77b4","#ea7e1eff","royalblue4","tomato2","darkgoldenrod","palegreen4","maroon","thistle3"),
                               wilcox_test = T,
                               stat_pos_mod = 1.15,
                               hjust.wilcox=0.8,
                               vjust.wilcox = 2.0,
                               size.wilcox=3.33,
                               step_mod=0,
                               hjust.wilcox.2=0.5,
                               vjust.wilcox.2=0,
                               sign_bar=0.8){

  if (!(Feature %in% rownames(Seu_object)) && !(Feature %in% names(Seu_object@meta.data))) {
    stop("Feature not found in Seurat Object!")
  }

  if (wilcox_test == T) {
    rstat <- system.file(package = "rstatix") # Make sure package is installed
    ifelse(nzchar(rstat), "", stop("Install rstatix R package for wilcox statistic!"))
  }

  if (is.null(ctrl.condition)) {
    stop("Please specify the ctrl condition as string!")
  }

  if (SeuV5==T) {
    rlang::warn("SeuV5 set to TRUE, if working with Seuratv4 or below change SeuV5 to FALSE", .frequency = "once", .frequency_id = "v5Mean")

    if (Feature %in% rownames(Seu_object)) {
      vln.df = data.frame(Feature = Seu_object[["RNA"]]$data[Feature,],
                          cluster = Seu_object[[group.by]])
    }else{
      vln.df = data.frame(Feature = FetchData(Seu_object, vars = Feature)[,1],
                          cluster = Seu_object[[group.by]])
    }



    df<-data.frame(group=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]]))
                   ,orig.ident = Seu_object$orig.ident)

    # add a second group for individual splitting and testing in the wilcoxon
    if (!is.null(group.by.2)) {
      vln.df[group.by.2] <- Seu_object[[group.by.2]]
      df[group.by.2] <- Seu_object[[group.by.2]]
    }

    #get expression values for genes from individual cells, add to df
    if (Feature %in% rownames(Seu_object)) {
      df[,Feature] <- Seu_object@assays$RNA$data[Feature,]
    }else{
      df[,Feature] <- FetchData(Seu_object, vars = Feature)[,1]
    }



  }

  if (SeuV5==F) {
    if (Feature %in% rownames(Seu_object)) {
      vln.df = data.frame(Feature = Seu_object[["RNA"]]@data[Feature,],
                          cluster = Seu_object[[group.by]])
    }else{
      vln.df = data.frame(Feature = FetchData(Seu_object, vars = Feature)[,1],
                          cluster = Seu_object[[group.by]])
    }

    df<-data.frame(group=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]]))
                   ,orig.ident = Seu_object$orig.ident)
    # add a second group for individual splitting and testing in the wilcoxon
    if (!is.null(group.by.2)) {
      vln.df[group.by.2] <- Seu_object[[group.by.2]]
      df[group.by.2] <- Seu_object[[group.by.2]]
    }

    #get expression values for genes from individual cells, add to df
    if (Feature %in% rownames(Seu_object)) {
      df[,Feature] <- Seu_object@assays$RNA@data[Feature,]
    }else{
      df[,Feature] <- FetchData(Seu_object, vars = Feature)[,1]
    }
  }


  df.melt <- melt(df)

  vln.df$group <- factor(vln.df[[group.by]]
                         ,levels = c(as.character(ctrl.condition),levels(factor(vln.df[[group.by]]))[!(levels(factor(vln.df[[group.by]])) %in% ctrl.condition)]))
  #create comparison list for wilcox, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with each other")
    group <- unique(Seu_object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) {
      ctrl.condition <- group[grep(pattern = paste(c("CTRL","Ctrl","ctrl","WT","Wt","wt"),collapse ="|")
                                   ,group)[1]]
    }


    #create ListTest
    ListTest <- list()
    for (i in 1:length(group)) {
      cndtn <- as.character(group[i])
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }

  #delete Null values, created by count index also reorder for better p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  if (!is.null(group.by.2)) {
    indices <- sapply(ListTest, function(x) match(x[2], vln.df[[group.by.2]]))
  } else{
    indices <- sapply(ListTest, function(x) match(x[2], vln.df[[group.by]]))
  }
  ListTest <- ListTest[order(indices)]

  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$group %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
      }
    }
    return(lst_filtered)
  }

  #group results and summarize
  if (is.null(group.by.2)) {
    df.melt.sum <- df.melt %>%
      dplyr::group_by(group, variable) %>%
      dplyr::summarise(Mean = mean(value))
  } else{
    df.melt.sum <- df.melt %>%
      dplyr::group_by(group, !!sym(group.by.2), variable) %>%
      dplyr::summarise(Mean = mean(value))
  }


  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df.melt.sum)

  # if (!is.null(group.by.2) && length(ListTest) > 1) {
  #   stop("The provided Seurat has more than two groups in group.by and you specified group.by.2, currently not supported (to crowded)!")
  # }


  # artificially set a group to 0 in all their expression values to create the error in the test
  # df.melt$value[df.melt$group == "Human-HFrEF"] <- 0
  # df.melt$value[df.melt$group == "Human-CTRLlv"] <- 0
  # df.melt$value[df.melt$group == "Human-HFrEF" & df.melt$cell_type == "Cardiomyocytes"] <- 0
  # df.melt$value[df.melt$group == "Human-CTRLlv" & df.melt$cell_type == "Cardiomyocytes"] <- 0

  #check before test if there are groups in the data which contain only 0 values and therefore let the test fail
  if (is.null(group.by.2)) {
    group_of_zero <- df.melt %>%
      dplyr::group_by(group) %>%
      summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      filter(all_zeros)

    if (nrow(group_of_zero) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df.melt <- df.melt %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(value = if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup()
    }
  } else{

    group_of_zero <- df.melt %>%
      dplyr::group_by(group, !!sym(group.by.2)) %>%
      summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      filter(all_zeros)

    #check now the result for multiple entries in group.by.2
    groupby2_check <- group_of_zero %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      summarise(group_count = n_distinct(group), .groups = "drop") %>%
      filter(group_count > 1)

    if (nrow(groupby2_check) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df.melt <- df.melt %>%
        dplyr::group_by(group, !!sym(group.by.2)) %>%
        dplyr::mutate(value = if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup()
    }
  }

  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE & is.null(group.by.2)) {
    stat.test <- df.melt %>%
      ungroup() %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(Seu_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p))

  }

  if (wilcox_test == TRUE & !is.null(group.by.2)) {
    stat.test <- df.melt %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(Seu_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p))

  }

  if (length(unique(vln.df[[group.by]])) >  length(vector_colors)) {
    stop(paste0("Only ", length(vector_colors)," colors provided, but ", length(unique(vln.df[[group.by]])), " needed!"))
  }

  #normal violin
  if(is.null(group.by.2)){
    p <- ggplot(vln.df, aes(x = group, y = Feature))+
      geom_violin(aes(fill = group), trim = T, scale = "width", )+
      geom_jitter(size = geom_jitter_args[1], width = geom_jitter_args[2], alpha = geom_jitter_args[3])+
      labs(title = Feature, y = "Expression Level")+
      xlab("")+
      ylab("")+
      theme_classic()+
      theme(plot.title = element_text(face = "bold", color = "black", hjust = 0.5, size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
            legend.position = "none")+
      scale_fill_manual(values = vector_colors)

    if (Feature %in% rownames(Seu_object)) {
      p_label="p = {p.adj}"
    } else{
      p_label="p = {p}"
    }

    if (wilcox_test == TRUE) {
      p = p + stat_pvalue_manual(stat.test, label = p_label, y.position = max(vln.df$Feature)*1.15, step.increase = 0.2)
    }
    return(p)
  }


  if(!is.null(group.by.2)){
    #plot
    p <- ggplot(vln.df, aes(x = !!sym(group.by.2), y = Feature, fill = !!sym(group.by)))+
      geom_violin(aes(fill = group), trim = T, scale = "width")+
      labs(title = Feature, y = "log(nUMI)")+
      xlab("")+
      theme_classic()+
      theme(plot.title = element_text(face = "bold", color = "black", hjus = 0.5, size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
            legend.position = "bottom",
            panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
            axis.line = element_line(colour = "black"),
            strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
            strip.text = element_text(colour = "black", size = 12),
      )+
      scale_fill_manual(values = vector_colors)

    p2 <- ggplot(vln.df, aes(x = !!sym(group.by.2), y = Feature, fill = factor(!!sym(group.by),levels = levels(vln.df$group))))+
      geom_boxplot(width=.1,color="grey", position = position_dodge(width = 0.9), outlier.shape = NA)+
      xlab("")+
      scale_fill_manual(values = rep("black",length(vector_colors)), name=group.by)+
      theme_classic()+
      theme(plot.title = element_text(face = "bold", color = "transparent", hjus = 0.5, size = 14),
            axis.title.y = element_text(face = "bold", color = "transparent", size = 14),
            axis.text.x = element_text(face = "bold", color = "transparent", angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(face = "bold", color = "transparent", hjust = 1, size = 14),
            legend.position = "bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            strip.background = element_rect(fill = "transparent", colour = NA),
            axis.ticks = element_blank(),
            legend.background = element_rect(fill = "transparent", color = NA),  # Transparent legend background
            legend.key = element_rect(fill = "transparent", color = NA),          # Transparent legend keys
            legend.title = element_text(face = "bold", color = "transparent"),          # Legend title styling
            legend.text = element_text(color = "transparent"),
            # strip.text = element_blank(),
      )

    if (wilcox_test == TRUE & !is.null(group.by.2)) {

      if (Feature %in% rownames(Seu_object)) {
        stat.test_plot <- stat.test %>%
          mutate(y.position = seq(from = max(Seu_object@assays$RNA$data[Feature,][!is.na(Seu_object@assays$RNA$data[Feature,])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
          mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
                 xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
                 xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)
      } else{
        stat.test_plot <- stat.test %>%
          mutate(y.position = seq(from = max(Seu_object[[Feature]][,1][!is.na(Seu_object[[Feature]][,1])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
          mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
                 xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
                 xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)
      }

      # dplyr::select(x.axis, y.position, p.adj)

      #if Feature not a gene than use the uncorrected p
      if (Feature %in% rownames(Seu_object)) {
        p_label="p = {p.adj}"
      } else{
        p_label="p = {p}"
      }

      p = p + stat_pvalue_manual(stat.test_plot,
                                 label = p_label,
                                 y.position = "y.position",
                                 # x="x",
                                 xmin = "xmin",
                                 xmax = "xmax",
                                 # xend="xend",
                                 # step.increase = 0.2,
                                 inherit.aes = FALSE,
                                 size = size.wilcox,
                                 angle= 0,
                                 hjust= hjust.wilcox.2,
                                 vjust = vjust.wilcox.2,
                                 tip.length = 0.02,
                                 bracket.size = sign_bar)

      p2 = p2 + stat_pvalue_manual(stat.test_plot,
                                   label = p_label,
                                   y.position = "y.position",
                                   # x="x",
                                   xmin = "xmin",
                                   xmax = "xmax",
                                   # xend="xend",
                                   # step.increase = 0.2,
                                   inherit.aes = FALSE,
                                   size = size.wilcox,
                                   angle= 0,
                                   hjust= hjust.wilcox.2,
                                   vjust = vjust.wilcox.2,
                                   tip.length = 0.02,
                                   bracket.size = sign_bar,
                                   color = "transparent")

    }
    plot_p <- cowplot::ggdraw() + cowplot::draw_plot(p) + cowplot::draw_plot(p2)
    print(plot_p)
  }



  if (returnValues==TRUE) {
    returnList <- list(vln.df, df.melt, stat.test)
    names(returnList) <- c("vln.df", "df.melt", "stat.test")
    return(returnList)
  }

  return(plot_p)
}



# Boxplot function for one or two given groups per gene, using a pseudo seurat approach
#' @author Mariano Ruz Jurado
#' @title Box Graph with wilcox test on single cell level
#' @description Creates a box plot using a pseudo-bulk approach and performs a Wilcoxon test on single-cell level.
#' Allows customization of outlier removal, statistical labels, and color schemes.
#' Supports comparison of conditions with optional second grouping.
#' Useful for visualizing gene expression and statistical differences.
#' @param Seu_object The seurat object
#' @param group.by group name to look for in meta data
#' @param group.by.2 second group name to look for in meta data
#' @param ctrl.condition select condition to compare to
#' @param outlier_removal Outlier calculation
#' @param vector_colors get the colours for the plot
#' @param wilcox_test If you want to have wilcoxon performed between ctrl.condition and given ones
#' @param stat_pos_mod modificator for where the p-value is plotted increase for higher
#' @param hjust.wilcox value for adjusting height of the text
#' @param vjust.wilcox value for vertical of text
#' @param size.wilcox value for size of text of statistical test
#' @param step_mod value for defining the space between one test and the next one
#' @param orderAxis vector for xaxis sorting, alphabetically by default
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import dplyr
#' @import magrittr
#' @import Seurat
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#'
#' (optional)
#' ListTest <- list()
#' ListTest[[1]] <- c("CTRL", "CONDITION")
#'
#' DO.BoxPlot(
#'   Seu_object = Seurat,
#'   Feature = "CDH5",
#'   sample.column="orig.ident",
#'   ListTest = ListTest,
#'   group.by="condition",
#'   ctrl.condition = "CTRL",
#' )
#' }
#'
#' @export
DO.BoxPlot <- function(Seu_object,
                               Feature,
                               sample.column = "orig.ident",
                               ListTest=NULL,
                               group.by = "condition",
                               group.by.2 = NULL,
                               ctrl.condition=NULL,
                               outlier_removal = T,
                               plot_sample=T,
                               vector_colors = c("#1f77b4","#ea7e1eff","royalblue4","tomato2","darkgoldenrod","palegreen4","maroon","thistle3"),
                               wilcox_test = T,
                               stat_pos_mod = 1.15,
                               step_mod = 0,
                               hjust.wilcox=0.5,
                               vjust.wilcox = 0.25,
                               size.wilcox=3.33,
                               hjust.wilcox.2=0.5,
                               vjust.wilcox.2=0,
                               sign_bar=0.8,
                               orderAxis=NULL
){


  #aggregate expression, pseudobulk to visualize the boxplot
  if (is.null(group.by.2)) {
    pseudo_Seu <- Seurat::AggregateExpression(Seu_object,
                                              assays = "RNA",
                                              return.seurat = T,
                                              group.by = c(group.by, sample.column),
                                              verbose = F)

    pseudo_Seu$celltype.con <- pseudo_Seu[[group.by]]

  } else{
    pseudo_Seu <- Seurat::AggregateExpression(Seu_object,
                                              assays = "RNA",
                                              return.seurat = T,
                                              group.by = c(group.by, group.by.2, sample.column),
                                              verbose = F)

    #cover the case of subsetted to only have one cell type
    if (length(unique(Seu_object@meta.data[[group.by.2]])) == 1) {
      pseudo_Seu@meta.data[[group.by.2]] <- unique(Seu_object@meta.data[[group.by.2]])
      pseudo_Seu@meta.data[[group.by]] <- gsub(paste0(".*(", paste(unique(Seu_object$condition), collapse = "|"), ").*"), "\\1", pseudo_Seu@meta.data[[group.by]])
    }

    pseudo_Seu$celltype.con <- paste(pseudo_Seu[[group.by]][,1], pseudo_Seu[[group.by.2]][,1], sep = "_")

  }


  if (Feature %in% rownames(Seu_object)) {
    df_Feature = data.frame(group=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]])),
                            Feature = Seu_object[["RNA"]]$data[Feature,],
                            cluster = Seu_object[[sample.column]][,1])
    df_Feature[,Feature] <- expm1(Seu_object@assays$RNA$data[Feature,])

  }else{
    df_Feature = data.frame(group=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]])),
                            Feature = FetchData(Seu_object, vars = Feature)[,1],
                            cluster = Seu_object[[sample.column]][,1])

    if (is.null(group.by.2)) {
      aggregated_meta <- Seu_object@meta.data %>%
        group_by(!!sym(group.by), !!sym(sample.column)) %>%
        summarise(Feature = mean(!!sym(Feature), na.rm = TRUE))

      aggregated_meta$comb <- paste(aggregated_meta[[group.by]],
                                    aggregated_meta[[sample.column]],
                                    sep = "_")

      #check if the matching agrees with the dash names if there were mutliple dashes in the names
      if(!any(is.na(aggregated_meta$Feature[match(pseudo_Seu[[sample.column]][,1], aggregated_meta$comb)]))){
        pseudo_Seu[[Feature]] <- aggregated_meta$Feature[match(pseudo_Seu[[sample.column]][,1], aggregated_meta$comb)]
      } else{
        unified_sep_column <- gsub("-", "_", pseudo_Seu[[sample.column]][,1])
        pseudo_Seu[[Feature]] <- aggregated_meta$Feature[match(unified_sep_column, aggregated_meta$comb)]
      }


    } else{
      #Compute mean for each orig.ident
      aggregated_meta <- Seu_object@meta.data %>%
        group_by(!!sym(group.by), !!sym(sample.column), !!sym(group.by.2)) %>%
        summarise(Feature = mean(!!sym(Feature), na.rm = TRUE))

      aggregated_meta$comb <- paste(aggregated_meta[[group.by]],
                                    aggregated_meta[[group.by.2]],
                                    aggregated_meta[[sample.column]],
                                    sep = "_")

      #cover the case of subsetted to only have one cell type
      if (!length(unique(Seu_object@meta.data[[group.by.2]])) == 1) {
        pseudo_Seu[[sample.column]][,1] <- gsub("_", "-", pseudo_Seu[[sample.column]][,1])
      } else{
        pseudo_Seu[[sample.column]][,1] <- paste0(gsub("_", "-", pseudo_Seu[[sample.column]][,1]), "-", unique(Seu_object@meta.data[[group.by.2]]))
        pseudo_Seu[[sample.column]][,1] <- gsub("_", "-", pseudo_Seu[[sample.column]][,1])
        rownames(pseudo_Seu@meta.data) <- pseudo_Seu[[sample.column]][,1]
      }
      aggregated_meta$comb <- gsub("_", "-", aggregated_meta$comb)
      pseudo_Seu[[Feature]] <- aggregated_meta$Feature[match(pseudo_Seu[[sample.column]][,1], aggregated_meta$comb)]

    }


  }

  # df_Feature<-data.frame(group=setNames(Seu_object[[group.by]][,group.by], rownames(Seu_object[[group.by]]))
  # ,orig.ident = Seu_object$orig.ident)
  # df_Feature[,Feature] <- expm1(Seu_object@assays$RNA$data[Feature,])

  #group results and summarize
  if (is.null(group.by.2)) {
    df_melt <- reshape2::melt(df_Feature) # melt in conditon since the second group might need to get added before the melt
    df_melt_sum <- df_melt %>%
      dplyr::group_by(group, variable) %>%
      dplyr::summarise(Mean = mean(value))
  } else{
    df_Feature[,{group.by.2}] <- setNames(Seu_object[[group.by.2]][,group.by.2], rownames(Seu_object[[group.by.2]]))
    df_melt <- reshape2::melt(df_Feature)
    df_melt_sum <- df_melt %>%
      dplyr::group_by(group, !!sym(group.by.2), variable) %>% #!!sym(), gets the actual variable name useable for dplyr functions
      dplyr::summarise(Mean = mean(value))
  }

  #create comparison list for wilcox, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with each other")
    group <- unique(Seu_object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) {
      ctrl.condition <- group[grep(pattern = paste(c("CTRL","Ctrl","WT","Wt","wt"),collapse ="|")
                                   ,group)[1]]
    }


    #create ListTest
    ListTest <- list()
    for (i in 1:length(group)) {
      cndtn <- as.character(group[i])
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }

  #delete Null values, created by count index also reorder for betetr p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  indices <- sapply(ListTest, function(x) match(x[2], df_melt_sum$group))
  ListTest <- ListTest[order(indices)]

  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$group %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
      }
    }
    return(lst_filtered)
  }

  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df_melt_sum)

  # if (!is.null(group.by.2) && length(ListTest) > 1) {
  #   stop("The provided Seurat has more than two groups in group.by and you specified group.by.2, currently not supported (to crowded)!")
  # }

  #check before test if there are groups in the data which contain only 0 values and therefore let the test fail
  if (is.null(group.by.2)) {
    group_of_zero <- df_melt %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      dplyr::filter(all_zeros)

    if (nrow(group_of_zero) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df_melt <- df_melt %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(value = dplyr::if_else(dplyr::row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup()
    }
  } else{

    group_of_zero <- df_melt %>%
      dplyr::group_by(group, !!sym(group.by.2)) %>%
      dplyr::summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      dplyr::filter(all_zeros)

    #check now the result for multiple entries in group.by.2
    groupby2_check <- group_of_zero %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      dplyr::summarise(group_count = dplyr::n_distinct(group), .groups = "drop") %>%
      dplyr::filter(group_count > 1)

    if (nrow(groupby2_check) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df_melt <- df_melt %>%
        dplyr::group_by(group, !!sym(group.by.2)) %>%
        dplyr::mutate(value = dplyr::if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        dplyr::ungroup()
    }
  }

  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE & is.null(group.by.2)) {
    stat.test <- df_melt %>%
      dplyr::ungroup() %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(Seu_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p))

  }

  #do statistix with rstatix + stats package add second group
  if (wilcox_test == TRUE & !is.null(group.by.2)) {
    stat.test <- df_melt %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(Seu_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p))

  }


  #pseudobulk boxplot
  theme_box <- function(){
    theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
        axis.line = element_line(colour = "black"),
        #facet_grid colors
        strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
        strip.text = element_text(colour = "black", size = 12),
        # legend.background = element_rect(colour = "grey", fill = "white"),
        # legend.box.background = element_rect(colour = "grey", size = 0.5),
      )
  }

  #TODO there might be a problem with the quantiles, recheck
  #SCpubr does not support outlier removal, therefore identify manually
  if (outlier_removal == T) {

    if (Feature %in% rownames(Seu_object)) {
      data_matrix <- pseudo_Seu@assays$RNA$data[Feature,]
    } else{
      data_matrix <- pseudo_Seu[[Feature]]
      data_matrix <- setNames(data_matrix[[Feature]], rownames(data_matrix))
    }

    for (grp2 in unique(pseudo_Seu[[group.by.2]][,1])) {
      for(grp in unique(pseudo_Seu[[group.by]][,1])){
        group_cells <- pseudo_Seu@meta.data[[group.by.2]] == grp2 & pseudo_Seu@meta.data[[group.by]] == grp
        subset_mat <- data_matrix[group_cells]

        Q1 <- quantile(subset_mat, 0.25)
        Q3 <- quantile(subset_mat, 0.75)
        IQR <- Q3 - Q1  # interquartile range calculation

        lower_bound <- Q1 - 1.5 * IQR #empirical rule derived from statistics. 1.5 as a default threshold
        upper_bound <- Q3 + 1.5 * IQR

        data_matrix_sub <- ifelse(subset_mat >= lower_bound & subset_mat <= upper_bound,
                                  subset_mat,
                                  NA)

        if (Feature %in% rownames(Seu_object)) {
          pseudo_Seu@assays$RNA$data[Feature, group_cells] <- data_matrix_sub
        } else{
          pseudo_Seu@meta.data[group_cells,Feature] <- data_matrix_sub
        }
      }
    }
  }

  if (!is.null(orderAxis)) {
    pseudo_Seu[[group.by]][,1] <- factor(pseudo_Seu[[group.by]][,1], levels = orderAxis)
  }


  if (is.null(group.by.2)) {
    p <- SCpubr::do_BoxPlot(sample = pseudo_Seu,
                            feature = Feature,
                            group.by = group.by,
                            order = F,
                            boxplot.width = 0.8,
                            legend.position = "right")

  } else {
    p <- SCpubr::do_BoxPlot(sample = pseudo_Seu,
                            feature = Feature,
                            group.by = group.by.2,
                            split.by = group.by,
                            boxplot.width = 0.8,
                            order = F)
  }

  if (plot_sample == T) {
    p <- p + geom_point(size=2, alpha=1, position = position_dodge(width = 0.8))
  }

  p <- p +
    scale_fill_manual(values = rep(vector_colors, 2))+ # 16 colours by default, 8 repeat after it
    theme_box()+
    theme(axis.text.x = element_text(color = "black",angle = 45,hjust = 1, size = 16, family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 16, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16,family = "Helvetica",face = "bold"),
          axis.title = element_text(size = 16, color = "black", family = "Helvetica"),
          plot.title = element_text(size = 16, hjust = 0.5, family = "Helvetica"),
          plot.subtitle = element_text(size = 16, hjust = 0, family = "Helvetica"),
          axis.line = element_line(color = "black"),
          strip.text.x = element_text(size = 16, color = "black", family = "Helvetica"),
          legend.text = element_text(size = 14, color = "black", family = "Helvetica"),
          legend.title = element_text(size = 14, color = "black", family = "Helvetica", face = "bold", hjust =0.5),
          legend.position = "bottom")

  #p
  # for only one group
  if (wilcox_test == TRUE & is.null(group.by.2)) {

    if (Feature %in% rownames(Seu_object)) {
      p_label="p = {p.adj}"
    } else{
      p_label="p = {p}"
    }


    if (Feature %in% rownames(Seu_object)) {
      stat.test_plot <- stat.test %>%
        mutate(y.position = seq(from= max(pseudo_Seu@assays$RNA$data[Feature,][!is.na(pseudo_Seu@assays$RNA$data[Feature,])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))
    } else{
      stat.test_plot <- stat.test %>%
        mutate(y.position = seq(from= max(pseudo_Seu@meta.data[,Feature][!is.na(pseudo_Seu@meta.data[,Feature])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))

    }


    p = p + stat_pvalue_manual(stat.test_plot,
                               label = "p = {p.adj}",
                               y.position = "y.position",
                               # step.increase = 0.2,
                               inherit.aes = FALSE,
                               size = size.wilcox,
                               angle= 0,
                               hjust= hjust.wilcox,
                               vjust = vjust.wilcox)
  }

  if (wilcox_test == TRUE & !is.null(group.by.2)) {


    if (Feature %in% rownames(Seu_object)) {
      p_label="p = {p.adj}"
    } else{
      p_label="p = {p}"
    }


    if (Feature %in% rownames(Seu_object)) {
      stat.test_plot <- stat.test %>%
        mutate(y.position = seq(from = max(pseudo_Seu@assays$RNA$data[Feature,][!is.na(pseudo_Seu@assays$RNA$data[Feature,])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
        mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
               xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
               xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)
    } else{
      stat.test_plot <- stat.test %>%
        mutate(y.position = seq(from = max(pseudo_Seu@meta.data[,Feature][!is.na(pseudo_Seu@meta.data[,Feature])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
        mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
               xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
               xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)
    }

    # dplyr::select(x.axis, y.position, p.adj)




    p = p + stat_pvalue_manual(stat.test_plot,
                               label = p_label,
                               y.position = "y.position",
                               # x="x",
                               xmin = "xmin",
                               xmax = "xmax",
                               # xend="xend",
                               # step.increase = 0.2,
                               inherit.aes = FALSE,
                               size = size.wilcox,
                               angle= 0,
                               hjust= hjust.wilcox.2,
                               vjust = vjust.wilcox.2,
                               tip.length = 0.02,
                               bracket.size = sign_bar)
  }
  return(p)
}


# Dotplot function for one or two given groups for multiple genes, using expression values
#' @author Mariano Ruz Jurado
#' @title DO Dot plot
#' @description This function generates a dot plot for multiple genes, comparing expression levels across one or two specified groups.
#' It supports both individual and pseudobulk expression calculations.
#' Highly variable customization options allow control over dot size, color scaling, annotations, and axis orientation.
#' The function integrates seamlessly with Seurat objects for single-cell RNA-seq analysis.
#' @param Seu_object The seurat Seu_object
#' @param group.by.x group name to plot on x-axis
#' @param group.by.y group name to look for in meta data
#' @param group.by.y2 second group name to look for in meta data
#' @param across.group.by.x calculate a pseudobulk expression approach for the x-axis categories
#' @param sort_x Vector sorting the xaxis
#' @param across.group.by.y calculate a pseudobulk expression approach for the y-axis categories, if group.by.y2 is specified Pseudobulk will be splitted by it
#' @param dot.size Vector of dot size
#' @param plot.margin = plot margins
#' @param midpoint midpoint in color gradient
#' @param Feature Genes or DF of interest, Data frame should have columns with gene and annotation information, e.g. output of FindAllMarkers
#' @param limits_colorscale Set manually colorscale limits
#' @param scale_gene IF True calculates the Z-score of the average expression per gene
#' @param hide_zero Removes dots for genes with 0 expression
#' @param annotation_x Adds annotation on top of x axis instead on y axis
#' @param point_stroke Defines the thickness of the black stroke on the dots
#' @param ... Further arguments passed to annoSegment function if annotation_x == T
#'
#' @import ggplot2
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#' @import ggtext
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.Dotplot(
#'   Seu_object = Seurat,
#'   Feature = c("CDH5","TTN","MALAT1"),
#'   group.by.x="condition"
#' )
#' }
#'
#'
#' @export
DO.Dotplot <- function(Seu_object,
                       Feature,
                       group.by.x = NULL,
                       group.by.y = NULL,
                       group.by.y2 = NULL,
                       across.group.by.x=F,
                       across.group.by.y=F,
                       sort_x=NULL,
                       dot.size = c(1,6),
                       plot.margin = c(1, 1, 1, 1),
                       midpoint = 0.5,
                       scale_gene=F,
                       returnValue = F,
                       log1p_nUMI=T,
                       hide_zero=T,
                       annotation_x=F,
                       annotation_x_position=0.25,
                       annotation_x_rev=F,
                       point_stroke=0.2,
                       limits_colorscale=NULL,
                       coord_flip=F,
                       ... ){

  if(!is.vector(Feature) && !is.data.frame(Feature)){
    stop("Feature is not a vector of strings or a data frame!")
  }

  if(across.group.by.x==T && across.group.by.y==T){
    stop("Both Pseudobulk options are set to true, please define just one!")
  }

  #type of Feature
  FeatureType <- mode(Feature)

  #check if Feature is a vector and annotation specified as true -> no cluster information provided for annotation
  if (is.vector(Feature) && annotation_x == T) {
    stop("Feature is a vector, but annotation_x is set to TRUE. If annotation on xaxis is wanted with specific cluster names you need to provide a dataframe with a column containing cluster names for the genes, like in Seurat::FindAllMarkers!")
  }

  #check the input if it is a data frame
  if (!is.vector(Feature)) {
    orig_DF <- Feature # save original df for annotation purposes
    orig_DF$cluster <- as.vector(orig_DF$cluster)
    cluster_name <- grep("cluster|group|annotation|cell", colnames(Feature), value = T) # Grep name of column, relevant for downstream assignments
    cluster <- unique(Feature[[grep("cluster|group|annotation|cell", colnames(Feature))]])
    Feature <- unique(Feature[[grep("gene|feature", colnames(Feature))]])

    if (is.null(cluster) || is.null(Feature)) {
      stop("Couldn't derive Cluster and Feature information from the provided Dataframe. \n Supported names for cluster: cluster|group|annotation|cell\n Supported names for Feature: gene|feature. \n Please make sure that your colnames in the provided Dataframe are supported.")
    }
  }

  # Create Feature expression data frame with grouping information
  geneExp <- expm1(FetchData(object = Seu_object, vars = Feature, layer = "data")) #

  #catch wrong handling of arguments
  if (is.null(group.by.x) && !is.null(group.by.y) && is.null(group.by.y2)) {
    stop("If you want to make a Marker Plot with just one provided group.by then please use group.by.x!")
  }

  geneExp$xaxis <- Seu_object@meta.data[[group.by.x]]

  if (!is.null(group.by.y) && is.null(group.by.y2)) {
    geneExp$id <- paste(Seu_object@meta.data[[group.by.y]], sep = "")
  } else if(!is.null(group.by.y) && !is.null(group.by.y2)){
    geneExp$id <- paste(Seu_object@meta.data[[group.by.y]], " (",
                        Seu_object@meta.data[[group.by.y2]], ")", sep = "")
  } else if(is.null(group.by.y) && is.null(group.by.y2)){
    geneExp$id <- paste(Seu_object@meta.data[[group.by.x]], sep = "")
  }

  # Include xaxis in the overall grouping
  data.plot <- lapply(X = unique(geneExp$id), FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, ]

    lapply(X = unique(data.use$xaxis), FUN = function(x_axis) {
      data.cell <- data.use[data.use$xaxis == x_axis, 1:(ncol(geneExp) - 2), drop = FALSE]
      avg.exp <- apply(X = data.cell, MARGIN = 2, FUN = function(x) {
        return(mean(x))
      })
      pct.exp <- apply(X = data.cell, MARGIN = 2, FUN = PercentAbove,
                       threshold = 0)

      res <- data.frame(id = ident, xaxis = x_axis, avg.exp = avg.exp, pct.exp = pct.exp * 100)
      res$gene <- rownames(res)
      return(res)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .) %>% data.frame()

  data.plot.res <- data.plot

  #add the cluster information to the plot if annotation_x is set to true and a dataframe was provided with cluster annotation
  #TODO Clean this part a bit up
  if (annotation_x==T && !is.null(cluster)) {
    data.plot.res <- purrr::map_df(seq_len(nrow(data.plot.res)), function(x){
      tmp <- data.plot.res[x,]
      tmp$celltype <- orig_DF[which(orig_DF[[grep("gene|feature", colnames(orig_DF))]] == tmp[[grep("gene|feature", colnames(tmp))]]), cluster_name][[1]]
      return(tmp)
    })

    # data.plot.res <- data.plot.res %>%
    #   dplyr::arrange(celltype)
  }

  #sort x by provided x-axis
  if (is.null(sort_x) && !is.null(group.by.y)) {
    data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = sort(unique(data.plot.res$xaxis)))
  } else if (!is.null(sort_x) && !is.null(group.by.y)){
    data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = sort_x)
  }

  #sort x by provided x-axis for just one group.by.x
  if (is.null(sort_x) && is.null(group.by.y)) {
    data.plot.res$gene <- factor(data.plot.res$gene, levels = sort(unique(data.plot.res$gene)))
  } else if (!is.null(sort_x) && is.null(group.by.y)){
    data.plot.res$gene <- factor(data.plot.res$gene, levels = sort_x)
  }

  #create grouping column for multiple grouping variables on the y-axis
  if (!is.null(group.by.y2)) {
    data.plot.res$group <- sapply(strsplit(as.character(data.plot.res$id),
                                           split = "\\(|\\)"), "[", 2)
  }

  if (hide_zero==T) {
    data.plot.res$pct.exp <- ifelse(data.plot.res$pct.exp == 0, NA, data.plot.res$pct.exp) # so fraction 0 is not displayed in plot
    data.plot.res <- data.plot.res[complete.cases(data.plot.res$pct.exp),]# remove empty lines
  }

  #create bulk expression for group.by.x
  if (across.group.by.x == T) {
    bulk_tmp <- data.plot.res %>%
      dplyr::group_by(id, gene) %>%
      summarise(avg.exp = mean(avg.exp),
                pct.exp = mean(pct.exp))
    bulk_tmp$xaxis <- "Pseudobulk"
    data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
    data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = c("Pseudobulk", setdiff(sort(unique(data.plot.res$xaxis)), "Pseudobulk")))
  }

  #create bulk expression for group.by.y, will divide by group.by.y2 if provided
  if (across.group.by.y == T) {
    if (is.null(group.by.y2)) {
      bulk_tmp <- data.plot.res %>% dplyr::group_by(xaxis, gene) %>%
        summarise(avg.exp = mean(avg.exp), pct.exp = mean(pct.exp))
      bulk_tmp$id <- "Pseudobulk"
      pseudo_levels <- unique(bulk_tmp$id)
      data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
      data.plot.res$id <- factor(data.plot.res$id, levels = c(pseudo_levels, setdiff(sort(unique(data.plot.res$id)), pseudo_levels)))
    } else{
      bulk_tmp <- data.plot.res %>% dplyr::group_by(xaxis, gene, group) %>%
        summarise(avg.exp = mean(avg.exp), pct.exp = mean(pct.exp))
      bulk_tmp$id <- paste0("Pseudobulk (", bulk_tmp$group, ")")
      data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
      pseudo_levels <- unique(bulk_tmp$id)
      data.plot.res$id <- factor(data.plot.res$id, levels = c(pseudo_levels, setdiff(sort(unique(data.plot.res$id)), pseudo_levels)))
    }
  }

  # get the scale pvalue for plotting
  if (log1p_nUMI==T) {
    data.plot.res$avg.exp.plot <- log1p(data.plot.res$avg.exp) # reapply the log transformation if wanted
  } else{
    data.plot.res$avg.exp.plot <- data.plot.res$avg.exp
  }

  #define how expression values are transformed
  if (scale_gene==T) {
    data.plot.res <- data.plot.res %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(z_avg_exp = (avg.exp - mean(avg.exp, na.rm=TRUE)) / sd(avg.exp, na.rm=TRUE)) %>%
      ungroup()
    exp.title = "Z-score expression \n per gene"
    fill.values = data.plot.res$z_avg_exp
    ###
  } else if(log1p_nUMI ==T){
    exp.title = "Mean log(nUMI) \n in group"
    fill.values = data.plot.res$avg.exp.plot
  } else{
    exp.title = "Mean nUMI \n in group"
    fill.values = data.plot.res$avg.exp.plot
  }

  #Define which columns to take for dotplot, it should be able to correctly capture one group.by.x
  if (identical(as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis)) && FeatureType=="list") { # go over input type
    # get rid of previous factoring to set new one, first alphabetical order on y
    data.plot.res$xaxis <- as.vector(data.plot.res$xaxis)
    data.plot.res$id <- factor(data.plot.res$id, levels = sort(unique(data.plot.res$id)))
    data.plot.res$gene <- factor(data.plot.res$gene, levels = orig_DF[order(orig_DF$cluster, decreasing = F),]$gene)

    if (annotation_x_rev == T) {
      data.plot.res$id <- factor(data.plot.res$id, levels = rev(sort(unique(data.plot.res$id))))
    }

    aes_var <- c("gene", "id")
    # there is a second case here for providing just a gene list which need to be adressed with the same aes_var
  } else if(identical(as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis))){

    # data.plot.res$id <- factor(data.plot.res$id, levels = sort(unique(data.plot.res$id)))
    # data.plot.res$gene <- factor(data.plot.res$gene, levels = sort(unique(data.plot.res$gene)))

    # #sort x by provided x-axis
    # if (is.null(sort_x)) {
    #   data.plot.res$id <- factor(data.plot.res$id, levels = sort(unique(data.plot.res$id)))
    # } else{
    #   data.plot.res$id <- factor(data.plot.res$id, levels = sort_x)
    # }

    aes_var <- c("gene", "id")

  } else{ # all other cases where group.by.y is specified
    data.plot.res$id <- factor(data.plot.res$id, levels = rev(sort(unique(data.plot.res$id))))
    aes_var <- c("xaxis", "id")

  }



  pmain <- ggplot2::ggplot(data.plot.res, ggplot2::aes(x = !!sym(aes_var[1]),y = !!sym(aes_var[2])))+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::xlab("") + ggplot2::ylab("")+ ggplot2::coord_fixed(clip = "off")+
    ggplot2::theme(plot.margin = ggplot2::margin(t = plot.margin[1],
                                                 r = plot.margin[2],
                                                 b = plot.margin[3],
                                                 l = plot.margin[4],
                                                 unit = "cm"),
                   axis.text = ggplot2::element_text(color = "black"),
                   legend.direction = "horizontal",
                   axis.text.x = element_text(color = "black",angle = 90,hjust = 1,vjust = 0.5, size = 14, family = "Helvetica"),
                   axis.text.y = element_text(color = "black", size = 14, family = "Helvetica"),
                   axis.title.x = element_text(color = "black", size = 14, family = "Helvetica"),
                   axis.title = element_text(size = 14, color = "black", family = "Helvetica"),
                   plot.title = element_text(size = 14, hjust = 0.5,face="bold", family = "Helvetica"),
                   plot.subtitle = element_text(size = 14, hjust = 0, family = "Helvetica"),
                   axis.line = element_line(color = "black"),
                   strip.text.x = element_text(size = 14, color = "black", family = "Helvetica", face = "bold"),
                   legend.text = element_text(size = 10, color = "black", family = "Helvetica"),
                   legend.title = element_text(size = 10, color = "black", family = "Helvetica", hjust =0),
                   legend.position = "right",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),)

  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = exp.title,
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.8,"cm"), # changes the width of the color legend
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 2),
                                  size = ggplot2::guide_legend(title = "Fraction of cells \n in group (%)",
                                                               title.position = "top", title.hjust = 0.5, label.position = "bottom",
                                                               override.aes = list(color = "black", fill = "grey50"),
                                                               keywidth = ggplot2::unit(0.5, "cm"), # changes the width of the precentage dots in legend
                                                               order = 1))

  dot.col = c("#fff5f0","#990000") # TODO change the scalefillgradient to +n in the else part
  gradient_colors <- c("#fff5f0", "#fcbba1", "#fc9272", "#fb6a4a", "#990000")
  # "#FFFFFF","#08519C","#BDD7E7" ,"#6BAED6", "#3182BD",
  if (length(dot.col) == 2) {
    breaks <- scales::breaks_extended(n=5)(range(fill.values))

    if (is.null(limits_colorscale)) {
      limits_colorscale <- c(min(range(fill.values))*.99,max(range(fill.values))*1.01)
    }

    if (max(breaks) > max(limits_colorscale)) {
      limits_colorscale[length(limits_colorscale)] <- breaks[length(breaks)]
    }

    pmain <- pmain + ggplot2::scale_fill_gradientn(colours = gradient_colors,
                                                   breaks = breaks,
                                                   #breaks = pretty(as.vector(quantile(fill.values)), n =10),
                                                   limits = limits_colorscale)
  }else{

    pmain <- pmain + ggplot2::scale_fill_gradient2(low = dot.col[1],
                                                   mid = dot.col[2],
                                                   high = dot.col[3],
                                                   midpoint = midpoint, name = "Gradient")
  }
  if (across.group.by.x == T) {

    pmain <- pmain +
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = pct.exp),shape = 21,stroke=point_stroke)+
      guides.layer +
      facet_grid(cols = vars(gene %>% factor(levels = Feature)), scales = "fixed")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(data.plot.res$pct.exp)*1.05,max(data.plot.res$pct.exp)*1.05))+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x=ggtext::element_markdown(color = "black",angle = 90,hjust = 1,vjust = 0.5, size = 14, family = "Helvetica"))+
      scale_x_discrete(labels = function(labels){
        labels <- ifelse(labels== "Pseudobulk", paste0("<b>", labels, "</b>"),labels)
        return(labels)
      })

  } else if(across.group.by.y == T){

    pmain <- pmain + ggplot2::geom_point(ggplot2::aes(fill = fill.values, size = pct.exp), shape = 21, stroke = point_stroke)+
      guides.layer + facet_grid(cols = vars(gene %>% factor(levels = Feature)), scales = "fixed") +
      ggplot2::scale_size(range = c(dot.size[1], dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n = 10)[seq(1, 10, by = 2)],
                                     limits = c(min(data.plot.res$pct.exp) * 1.05, max(data.plot.res$pct.exp) * 1.05))+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.y = ggtext::element_markdown(color = "black", angle = 0, hjust = 1, vjust = 0.5, size = 14, family = "Helvetica"))+
      scale_y_discrete(labels = function(labels) {
        labels <- ifelse(labels %in% pseudo_levels, paste0("<b>", labels, "</b>"), labels)
        return(labels)
      })

  } else if(identical(as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis))){

    pmain <- pmain +
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = pct.exp),shape = 21,stroke=point_stroke)+
      guides.layer +
      # facet_wrap(~facet_group, scales="free_x")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)])*.95,max(data.plot.res$pct.exp)*1.05))+
      theme(panel.spacing = unit(0, "lines"))

    if (annotation_x==T) {
      #use the annotation function from jjAnno package
      # jjA <- system.file(package = "jjAnno") # Make sure package is installed
      # ifelse(nzchar(jjA), "", stop("Install jjAnno CRAN package for annotation on xaxis:!"))

      # pmain <- pmain +
      #   ggplot2::theme(axis.text.y = element_blank(),
      #                  axis.ticks.y = element_blank())

      plot_max_y <- ggplot_build(pmain)
      plot_max_y <- plot_max_y$layout$panel_params[[1]]$y.range[2] + annotation_x_position

      pmain <- .annoSegment(
        object = pmain,
        annoPos = "top",
        aesGroup = T,
        aesGroName = "celltype",
        fontface = "bold",
        fontfamily = "Helvetica",
        pCol = rep("black", length(cluster)),
        textCol = rep("black", length(cluster)),
        addBranch = T,
        branDirection = -1,
        addText = T,
        yPosition = plot_max_y,
        # textSize = 14,
        # hjust = 0.5,
        # vjust = 0,
        # textRot = 0,
        # segWidth = 0.3,
        # lwd = 3
        ...
      )


    }

    if (coord_flip==T) {
      pmain <- pmain + ggplot2::coord_flip()
    }

    if (coord_flip==T && annotation_x==T) {
      warning("Annotation_x and coord_flip set on TRUE, might result in unwanted behaviour!")
    }

  } else{

    pmain <- pmain +
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = pct.exp),shape = 21,stroke=point_stroke)+
      guides.layer +
      facet_grid(cols = vars(gene %>% factor(levels = Feature)), scales = "fixed")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)])*.95,max(data.plot.res$pct.exp)*1.05))+
      theme(panel.spacing = unit(0, "lines"))
  }

  if(returnValue == T){
    return(data.plot.res)
  }
  return(pmain)

}

theme_box <- function(){
  theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
      panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
      axis.line = element_line(colour = "black"),
      #facet_grid colors
      strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
      strip.text = element_text(colour = "black", size = 12),
      # legend.background = element_rect(colour = "grey", fill = "white"),
      # legend.box.background = element_rect(colour = "grey", size = 0.5),
    )
}

#' @author Mariano Ruz Jurado
#' @title DO CellComposition
#' @description Computes and visualizes cell composition changes from a Seurat object using Scanpro.
#' Integrates R and Python via reticulate to convert Seurat to AnnData and run Scanpro.
#' Generates customizable plots with options for transformation, grouping, and bootstrapping.
#' Returns a ggplot object or a list containing plot data and visualization.
#' @param Seu_object The seurat object
#' @param assay_normalized Assay with raw counts
#' @param cluster_column Column in meta data which will be used to segment the bar plot
#' @param sample_column Column in meta data containing individual sample names
#' @param condition_column Column in meta data plotted on the xaxis
#' @param transform_method Method of transformation of proporties, default: "logit"
#' @param sort_x Vector sorting the xaxis
#' @param sort_fill Vector sorting the bar segments
#' @param sub_ident vector to subset the whole plot by
#' @param scanpro_plots Boolean, will create plots provided in scanpro package, default: FALSE
#' @param scanpro_group Defines the clusters showed in scanpro plots
#' @param outputFolder Scanpro plots will be saved in this directory, defaults to current working directory
#' @param return_df Boolean, makes a list with first entry beeaing the dataframe used to create the ggplot and second entry is the plot
#' @param bar_colors (named) vector of colors used for plotting, the names will be sued to factorize the segments
#' @param n_reps Number of replicates generated by scanpro
#' @param legend.pos.x adjusts the position of the legend in horizontal
#' @param legend.pos.y adjusts the position of the legend in vertical
#' @param cowplot_width Changes the width on the plotting device for plot
#' @param cowlegend_width Changes the width on the plotting device for legend, useful for adjusting the legend position in combination with legend.pos.x and y
#' @param ... Further arguments passed to scanpro plotting functions
#'
#' @return ggplot object or list
#'
#' @import reticulate
#' @import ggplot2
#' @import ggalluvial
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.CellComposition(
#'   Seu_object = Seurat,
#'   cluster_column="cell_type",
#'   condition_column="condition",
#'   scanpro_plots=T,
#'   bar_colors=c("#ff7f0e","tomato2", "darkgoldenrod2","#aec7e8"), # can also be a named vector
#'   n_reps=5)
#'
#'
#'
#' }
#'
#' @export
DO.CellComposition <- function(Seu_object,
                               assay_normalized="RNA",
                               cluster_column="seurat_clusters",
                               sample_column="orig.ident",
                               condition_column="condition",
                               transform_method="logit",
                               sort_x=NULL,
                               sub_ident=NULL,
                               sort_fill=NULL,
                               scanpro_plots=F,
                               scanpro_group=NULL,
                               outputFolder=NULL,
                               return_df=FALSE,
                               bar_colors=NULL,
                               n_reps=NULL,
                               legend.pos.x=0.48,
                               legend.pos.y=0,
                               cowplot_width=0.9,
                               cowlegend_width=0.9,
                               ...
){

  if(scanpro_plots==T && is.null(outputFolder)){
    warning(paste0("scanpro_plots will be saved in: ", getwd(), "\n"))
    outputFolder <- getwd()
  }

  if(!is.null(n_reps)){
    warning(paste0("Bootstrapping method activated with ", n_reps, " simulated replicates!\n."))
  }

  # Make sure R Zellkonverter package is installed
  zk <- system.file(package = "zellkonverter")
  ifelse(nzchar(zk), "", stop("Install zellkonverter R package for Seurat tranformation!"))

  # Make sure R reticulate package is installed
  rt <- system.file(package = "reticulate")
  ifelse(nzchar(rt), "", stop("Install reticulate R package for Python usage in R!"))

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not correctly configured. Run "usethis::edit_r_environ()" to specify your Python instance'))
  }

  #Make Anndata object
  DefaultAssay(Seu_object) <- assay_normalized
  Seu_object_sce <- Seurat::as.SingleCellExperiment(Seu_object, assay = assay_normalized)
  AnnData_counts <- zellkonverter::SCE2AnnData(Seu_object_sce)



  # Import scanpro
  sc <- import("scanpro.scanpro")
  merge_design_props_fct <- py_get_attr(sc$ScanproResult, "_merge_design_props")

  if(is.null(n_reps)){
    out = sc$scanpro(AnnData_counts,
                     clusters_col=cluster_column,
                     samples_col=sample_column,
                     conds_col=condition_column,
                     transform = transform_method)
  } else{
    out = sc$scanpro(AnnData_counts,
                     clusters_col=cluster_column,
                     conds_col=condition_column,
                     transform = transform_method,
                     n_reps=as.integer(n_reps))
  }

  if (scanpro_plots==T) {
    if (!dir.exists(paste0(outputFolder,"/scanpro"))) {
      dir.create(paste0(outputFolder,"/scanpro"))
    }
  }

  #Plot scanpro plots if requested
  if(scanpro_plots==T && is.null(scanpro_group)){
    out$plot(kind="boxplot", ...)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_complete_box_plot.svg"), dpi=300, bbox_inches="tight")
    plt$close()

    out$plot(kind="barplot", ...)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_complete_bar_plot.svg"), dpi=300, bbox_inches="tight")
    plt$close()

    #As stacked bar plot
    out$plot_samples(stacked=TRUE)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_complete_stacked_bar.svg"), dpi=300, bbox_inches="tight")
    plt$close()

    #As NON-stacked bar plot
    out$plot_samples(x=cluster_column)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_complete_notStacked_bar.svg"), dpi=300, bbox_inches="tight")
    plt$close()

  }

  if(scanpro_plots==T && !is.null(scanpro_group)){
    out$plot(kind="boxplot",clusters=scanpro_group, ...)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_",scanpro_group,"_box_plot.svg"), dpi=300, bbox_inches="tight")
    plt$close()

    out$plot(kind="barplot",clusters=scanpro_group, ...)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_",scanpro_group,"_bar_plot.svg"), dpi=300, bbox_inches="tight")
    plt$close()

    #As stacked bar plot
    out$plot_samples(stacked=TRUE)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_complete_stacked_bar.svg"), dpi=300, bbox_inches="tight")
    plt$close()

    #As NON-stacked bar plot
    out$plot_samples(x=cluster_column)
    plt$savefig(paste0(outputFolder, "/scanpro/CellComposition_complete_notStacked_bar.svg"), dpi=300, bbox_inches="tight")
    plt$close()
  }

  #assign results df
  results_df <- out$results
  results_counts <- out$counts
  results_counts <- tibble::rownames_to_column(results_counts, var="orig.ident")
  results_counts$condition <- result_vector <- apply(out$design, 1, function(x) names(x)[x == 1])
  results_counts <- reshape2::melt(results_counts)

  #per sample proportions
  prop_merge <- merge_design_props_fct(out)
  prop_merge_melt <- reshape2::melt(prop_merge)

  #add count numbers
  prop_merge_melt$counts <- results_counts$value

  summary_df <- prop_merge_melt %>%
    group_by(!!sym(condition_column), variable) %>%
    summarize(freq = sum(counts), .groups = 'drop')

  #convert Frequency to proportion
  prop_df <- summary_df %>%
    group_by(!!sym(condition_column)) %>%
    mutate(proportion = freq / sum(freq)) %>%
    ungroup()

  #assign p-values # keep an eye out if this causes bugs
  if (!is.null(results_df$adjusted_p_values)) {
    prop_df$p_val <- rep(results_df$adjusted_p_values, length(unique(prop_df[[condition_column]])))
  } else{
    prop_df$p_val <- rep(results_df$p_values, length(unique(prop_df[[condition_column]])))
  }

  prop_df$p_val <- ifelse(as.numeric(prop_df$p_val) > 0.05, ">0.05", "<0.05")

  prop_df$variable<- as.vector(prop_df$variable)

  # sum(prop_df[prop_df[[condition_column]] == colnames(out$design)[1],]$proportion)
  # sum(prop_df[prop_df[[condition_column]] == colnames(out$design)[2],]$proportion)

  if (is.null(bar_colors)) {
    bar_colors <- rep(c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
      "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7",
      "#dbdb8d","#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
    ),5)
  } else {

    #if bar_colors is a named vector it can be used for adjusting the fill colors
    if (!is.null(names(bar_colors))) {
      prop_df$variable <- factor(prop_df$variable, levels = names(bar_colors))
    }
  }

  # length(prop_df$variable)

  if (!is.null(sort_x)) {
    prop_df[[condition_column]] <- factor(prop_df[[condition_column]], levels = sort_x)
  }

  #sub by sub_ident if given
  if(!is.null(sub_ident)){
    prop_df <- subset(prop_df, subset = variable %in% sub_ident)
  }

  #sort fill
  if (!is.null(sort_fill)) {
    prop_df$variable <- factor(prop_df$variable, levels = sort_fill)
  }

  p1 <- ggplot(prop_df,aes(y = proportion, x = !!sym(condition_column))) +
    geom_col(aes(fill = variable), width = 0.5, color = "black") +
    scale_fill_manual(values = bar_colors)+
    scale_y_continuous(expand = c(0,0)) +
    xlab("")+
    theme_classic()+
    theme(plot.title = element_text(face = "bold", color = "black", hjus = 0.5, size = 14),
          axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Helvetica"),
          axis.text.x = element_text(face = "bold", color = "black", angle = 0, hjust = 0.5, size = 14, family = "Helvetica"),
          axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14, family = "Helvetica"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
          strip.text = element_text(colour = "black", size = 12))

  #make just the lines in a different plot
  p2 <- ggplot(prop_df, aes(y = proportion, x = !!sym(condition_column))) +
    geom_flow(aes(alluvium = variable, color = ifelse(p_val == "<0.05", "sig", "ns")),
              alpha= .9, lty = 2, fill = "white",
              curve_type = "linear", width = .5) +
    geom_col(aes(fill = variable), width = .5, color = NA, fill = "transparent") +
    #TODO add percentage numbers to plot
    # geom_text(aes(label = scales::percent(proportion),
    #               color=ifelse(p_val == "<0.05", "sig", "ns")),
    #           position = position_stack(vjust=0.5),
    #           size = 4) +
    # geom_text(aes(label= scales::percent(after_stat(y)),
    #               color=ifelse(p_val == "<0.05", "sig", "ns")),
    #           stat = "summary", fun="sum",
    #           position = position_fill(vjust = 0.5), size = 4)+
    scale_color_manual(values = c("sig" = "black", "ns" = "transparent")) +
    scale_fill_manual(values = bar_colors)+
    scale_y_continuous(expand = c(0,0)) +
    xlab("")+
    theme_classic()+
    theme(plot.title = element_text(face = "bold", color = "transparent", hjus = 0.5, size = 14),
          axis.title.y = element_text(face = "bold", color = "transparent", size = 14),
          axis.text.x = element_text(face = "bold", color = "transparent", angle = 0, hjust = 0.5, size = 14),
          axis.text.y = element_text(face = "bold", color = "transparent", hjust = 1, size = 14),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "transparent", colour = NA),
          axis.ticks = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.title = element_text(face = "bold", color = "transparent"),
          legend.text = element_text(color = "transparent"))


  #add legend
  p3 <- ggplot(prop_df,aes(y = proportion, x = !!sym(condition_column))) +
    geom_col(aes(fill = variable), width = .5, color="black") +
    scale_fill_manual(values = bar_colors)+
    scale_y_continuous(expand = c(0,0)) +
    xlab("")+
    theme_classic()+
    theme(plot.title = element_text(face = "bold", color = "transparent", hjus = 0.5, size = 14),
          axis.title.y = element_text(face = "bold", color = "transparent", size = 14),
          axis.text.x = element_text(face = "bold", color = "transparent", angle = 0, hjust = 0.5, size = 14),
          axis.text.y = element_text(face = "bold", color = "transparent", hjust = 1, size = 14),
          legend.position = "right",
          legend.title = element_text(size = 10, color = "black",hjust = 0.5,face = "bold", family = "Helvetica"),
          legend.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "transparent", colour = NA),
          axis.ticks = element_blank())

  guides.layer <- ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell type",
                                                               title.position = "top",
                                                               title.hjust = 0.5,
                                                               label.position = "right",
                                                               keywidth = ggplot2::unit(0.2, "cm"),
                                                               keyheight = ggplot2::unit(0.5, "cm"),
                                                               order = 1))

  p3 <- p3 + guides.layer
  legend_only <- cowplot::get_legend(p3)


  plot_p <- cowplot::ggdraw() +
    cowplot::draw_plot(p1, width = cowplot_width) +
    cowplot::draw_plot(p2, width = cowplot_width) +
    cowplot::draw_plot(legend_only, x = legend.pos.x, y = legend.pos.y, width = cowlegend_width) # TODO find a way to elongate the plot so that the legend automatically fits
  plot_p



  if (return_df == T) {
    return(list(prop_df, plot_p))
  } else{
    return(plot_p)
  }

}

#' @author Mariano Ruz Jurado
#' @title DO Bar plot for GSEA df result
#' @description This function generates a split barplot. This is a plot where the top 10 Go terms
#'are shown, sorted based on a column ('col_split'). Two conditions are shown at the same
#'time. One condition is shown in the positive axis, while the other in the negative one.
#'The condition to be shown as positive is set with 'pos_col'.
#'
#'The GO terms will be shown inside the bars, if the term is too long, using 'cutoff',
#'you can control the maximum number of characters per line.
#'
#' Pre-filter of the dataframe to contain significant Terms is recommended
#'
#' @param df_GSEA dataframe with the results of a gene set enrichment analysis
#' @param term_col column in the dataframe that contains the terms
#' @param col_split: column in the dataframe that will be used to sort and split the plot
#' @param cond_col: column in the dataframe that contains the condition information
#' @param pos_cond: condition that will be shown in the positive side of the plot
#' @param cutoff: maximum number of characters per line
#' @param log10_transform: if col_split contains values between 0 and 1, assume they are pvals and apply a -log10 transformation
#' @param figsize: figure size
#' @param topN: how many terms are shown
#' @param path: path to save the plot
#' @param filename: filename for the plot
#' @param spacing: space to add between bars and origin. It is a percentage value, indicating that the bars start at 5 % of the maximum X axis value.
#' @param txt_size: size of the go terms text
#' @param alpha_colors: alpha value for the colors of the bars
#' @param colors_pairs: colors for each condition (1st color --> negative axis; 2nd color --> positive axis)
#' @param title: title of the plot
#' @param showP: if False, the axis is return
#'
#'
#' @return: None or the axis
#'
#' @import reticulate
#'
#' @examples
#' \dontrun{
#' #df with GSEA results e.g. Metascape results
#' df_GSEA <- read.csv("path/to/GSEA_results.csv")
#'
#' #Please add a 'celltype' and 'condition' column
#' head(df_GSEA)
#' #Columns should include: "Term", "Combined.Score", "conditon", "celltype", etc.
#'
#' #Run SplitBarGSEA visualization
#' DO.SplitBarGSEA(df_GSEA = df_GSEA,
#'                 term_col = "Term",
#'                 col_split = "Combined.Score",
#'                 cond_col = "condition",          # condition column (e.g., experimental group)
#'                 pos_cond = "Treatment",
#'                 cutoff = 40,                 # max term length before truncation
#'                 log10_transform = TRUE,
#'                 figsize = c(12, 8),
#'                 topN = 10,                   # top N GO terms to show
#'                 color_pairs = c("sandybrown", "royalblue"),
#'                 alpha_colors = 0.3,
#'                 path = "results/plots",
#'                 spacing = 5,
#'                 txt_size = 12,
#'                 filename = "SplitBar.svg",
#'                 title = "Top 10 GO Terms in each Condition: ",
#'                 showP = FALSE,
#'                 celltype = "all")
#' }
#'
#' @export
DO.SplitBarGSEA <- function(df_GSEA,
                            term_col,
                            col_split,
                            cond_col,
                            pos_cond,
                            cutoff=40,
                            log10_transform=TRUE,
                            figsize=c(12,8),
                            topN=10,
                            colors_pairs=c("sandybrown","royalblue"),
                            alpha_colors=0.3,
                            path=NULL,
                            spacing=5,
                            txt_size=12,
                            filename="SplitBar.svg",
                            title="Top 10 GO Terms in each Condition: ",
                            showP=F,
                            celltype="all")
{

  if (!"celltype" %in% colnames(df_GSEA)) {
    stop("Provided data frame has no column named 'celltype'! Please add a column with cell type information.")
  }

  #source PATH to python script in install folder
  path_py <- system.file("python", "GSEA_split_bar.py", package = "DOtools")
  source_python(path_py)

  #Initialize matplot package
  plt <- import("matplotlib.pyplot")

  #Subset by given celltype names
  if (!celltype=="all") {
    df_GSEA <- df_GSEA[df_GSEA$celltype %in% celltype,]
  }

  #Run for each celltype the split_bar_gsea function
  for (celltype in unique(df_GSEA$celltype)) {
    df_GSEA_sub <- df_GSEA[df_GSEA$celltype == celltype,]
    df_GSEA_sub$Term <- sapply(strsplit(df_GSEA_sub$Term, " \\(GO"), function(x)x[1])

    #Conversion of the data.frame
    df_GSEA_sub_pd <- r_to_py(df_GSEA_sub)

    #run python function with given arguments
    plot <- split_bar_gsea(df = df_GSEA_sub_pd,
                           term_col = term_col,
                           col_split = col_split,
                           cond_col = cond_col,
                           pos_cond = pos_cond,
                           title = paste0(title, celltype),
                           cutoff = cutoff,
                           log10_transform = log10_transform,
                           figsize = figsize,
                           topN = topN,
                           colors_pairs = colors_pairs,
                           alpha_colors = alpha_colors,
                           path = path,
                           spacing = spacing,
                           txt_size = txt_size,
                           filename = filename,
                           showP = showP)

    if (showP == TRUE) {
      #x-title settings
      # plot$set_xlabel("Combined Score")
      plt$show()
    }

    if (!is.null(path)) {
    #Save under provided PATH
    plt$savefig(paste0(path, celltype, "_", filename), dpi=300, bbox_inches="tight")
    }

  }
}

#' @author Mariano Ruz Jurado
#' @title DO Correlation Plot for visualizing similarity between categories
#' @description#' Generates a correlation heatmap from expression data to visualize similarity across sample groups.
#' Allows customization of plot type, correlation method, and color scaling using the ggcorrplot2 and ggplot2 architectures.
#' Ideal for comparing transcriptional profiles between conditions or clusters.
#' @param Seu_obj Seurat Object
#' @param group_by Column to aggregate the expression over it, default "orig.ident"
#' @param assay Assay in object to use, default "RNA"
#' @param features What genes to include by default all, default "None"
#' @param method Correlation method, default "spearman"
#' @param plotdesign Plot design, default "circle"
#' @param plottype Show the full plot or only half of it, default "full"
#' @param auto_limits Automatically rescales the colour bar based on the values in the correlation matrix, default "TRUE"
#' @param outline.color the outline color of square or circle. Default value is "white".
#' @param colormap Defines the colormap used in the plot, default c("royalblue4", "royalblue2","firebrick","firebrick4")
#' @param lab_size Size to be used for the correlation coefficient labels. used when lab = TRUE.
#' @param lab logical value. If TRUE, add correlation coefficient on the plot.
#' @param lab_col color to be used for the correlation coefficient labels. used when lab = TRUE.
#' @param ... Additionally arguments passed to ggcorrplot function
#'
#'
#' @return ggplot2
#'
#' @import ggcorrplot
#' @import ggplot2
#' @import Seurat
#'
#' @examples
#' \dontrun{
#' DO.Correlation(
#'   Seu_obj = seurat_object,
#'   group_by = "orig.ident",
#'   assay = "RNA",
#'   features = NULL,
#'   method = "spearman",
#'   plotdesign = "square",
#'   plottype = "full",
#'   auto_limits = TRUE,
#'   outline.color = "white",
#'   colormap = c("royalblue4", "royalblue2", "firebrick", "firebrick4"),
#'   lab_size = 10,
#'   lab = TRUE,
#'   lab_col = "white"
#' )
#' }
#'
#' @export
DO.Correlation <- function(Seu_obj,
                           group_by="orig.ident",
                           assay="RNA",
                           features=NULL,
                           method="spearman",
                           plotdesign="square",
                           plottype="full",
                           auto_limits=TRUE,
                           outline.color="white",
                           colormap=c("royalblue4", "royalblue2","firebrick","firebrick4"),
                           lab_size=10,
                           lab=TRUE,
                           lab_col="white",
                           ...){

  #Aggregate Expression, creating Pseudobulk
  corr_frame <- AggregateExpression(Seu_obj,
                                    assays = assay,
                                    features = features,
                                    group.by = group_by)$RNA

  corr_frame <- as.data.frame(corr_frame)

  #Calculate correlation
  corr_df <- cor(corr_frame, method = method)


  # Auto color limits based on matrix
  if (auto_limits == TRUE) {
    min_val <- min(corr_df, na.rm = TRUE)*0.99
    max_val <- max(corr_df, na.rm = TRUE)*1.01
  } else {
    min_val <- -1
    max_val <- 1
  }

  pmain <- ggcorrplot(corr_df,
                      method = plotdesign,
                      type = plottype,
                      colors = colormap,
                      outline.color = outline.color,
                      lab_col = lab_col,
                      lab_size = lab_size,
                      lab = lab,
                      ...)+
    scale_fill_gradientn(colors = colormap,
                        limits = c(min_val, max_val),
                        name = paste0(tools::toTitleCase(method), "Correlation"))+
    ggplot2::theme(axis.text = ggplot2::element_text(color = "black"),
                   legend.direction = "horizontal",
                   axis.text.x = element_text(color = "black",angle = 0,hjust = 0.5, size = 14, family = "Helvetica"),
                   axis.text.y = element_text(color = "black",angle=90,,hjust = 0.5, size = 14, family = "Helvetica"),
                   axis.title = element_text(size = 14, color = "black", family = "Helvetica"),
                   plot.title = element_text(size = 14, hjust = 0.5, face="bold", family = "Helvetica"),
                   plot.subtitle = element_text(size = 14, hjust = 0, family = "Helvetica"),
                   strip.text.x = element_text(size = 14, color = "black", family = "Helvetica", face = "bold"),
                   legend.text = element_text(size = 10, color = "black", family = "Helvetica"),
                   legend.title = element_text(size = 10, color = "black", family = "Helvetica", hjust =0),
                   legend.position = "bottom")

  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = paste0(tools::toTitleCase(method), " \nCorrelation"),
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.8,"cm"), # changes the width of the color legend
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 1))
  pmain <- pmain + guides.layer

  return(pmain)

}

# AnnoSegment function conservation: the original author is no longer maintaining it on CRAN and it would be a shame to lose it
#' @author Mariano Ruz Jurado (edited from: Jun Zhang)
#' @title Annotation modifier for plots
#' @description Used for segment the plot for further annotations
#' @param object ggplot list. Default(NULL).
#' @param relSideDist The relative distance ratio to the y axis range. Default(0.1).
#' @param aesGroup Whether use your group column to add rect annotation. Default("FALSE").
#' @param aesGroName The mapping column name. Default(NULL).
#' @param annoPos The position for the annotation to be added. Default("top").
#' @param xPosition The x axis coordinate for the segment. Default(NULL).
#' @param yPosition The y axis coordinate for the segment. Default(NULL).
#' @param pCol The segment colors. Default(NULL).
#' @param segWidth The relative segment width. Default(1).
#' @param lty The segment line type. Default(NULL).
#' @param lwd The segment line width. Default(NULL).
#' @param alpha The segment color alpha. Default(NULL).
#' @param lineend The segment line end. Default("square").
#' @param annoManual Whether annotate by yourself by supplying with x and y coordinates. Default(FALSE).
#' @param mArrow Whether add segment arrow. Default(FALSE).
#' @param addBranch Whether add segment branch. Default(FALSE).
#' @param bArrow Whether add branch arrow. Default(FALSE).
#' @param branDirection The branch direction. Default(1).
#' @param branRelSegLen The branch relative length to the segment. Default(0.3).
#' @param addText Whether add text label on segment. Default(FALSE).
#' @param textCol The text colors. Default(NULL).
#' @param textSize The text size. Default(NULL).
#' @param fontfamily The text fontfamily. Default(NULL).
#' @param fontface The text fontface. Default(NULL).
#' @param textLabel The text textLabel. Default(NULL).
#' @param textRot The text angle. Default(NULL).
#' @param textHVjust The text distance from the segment. Default(0.2).
#' @param hjust The text hjust. Default(NULL).
#' @param vjust The text vjust. Default(NULL).
#' @param myFacetGrou Your facet group name to be added with annotation when object is a faceted object. Default(NULL).
#' @param aes_x = NULL You should supply the plot X mapping name when annotate a facetd plot. Default(NULL).
#' @param aes_y = NULL You should supply the plot Y mapping name when annotate a facetd plot. Default(NULL).
#' @keywords internal
.annoSegment <- function (object = NULL,
                         relSideDist = 0.1,
                         aesGroup = FALSE,
                         aesGroName = NULL,
                         annoPos = "top",
                         xPosition = NULL,
                         yPosition = NULL,
                         pCol = NULL,
                         segWidth = 1,
                         lty = NULL,
                         lwd = 10,
                         alpha = NULL,
                         lineend = "square",
                         annoManual = FALSE,
                         mArrow = NULL,
                         addBranch = FALSE,
                         bArrow = NULL,
                         branDirection = 1,
                         branRelSegLen = 0.3,
                         addText = FALSE,
                         textCol = NULL,
                         textSize = NULL,
                         fontfamily = NULL,
                         fontface = NULL,
                         textLabel = NULL,
                         textRot = 0,
                         textHVjust = 0.2,
                         hjust = NULL,
                         vjust = NULL,
                         myFacetGrou = NULL,
                         aes_x = NULL,
                         aes_y = NULL)
{
  facetName <- names(object$facet$params$facets)
  if (is.null(myFacetGrou) & !is.null(facetName)) {
    myFacetGrou <- unique(data[, facetName])[1]
  }
  else if (!is.null(myFacetGrou) & !is.null(facetName)) {
    myFacetGrou <- myFacetGrou
  }
  else {
  }
  data <- object$data
  if (is.null(facetName)) {
    aes_x <- ggiraphExtra::getMapping(object$mapping, "x")
    aes_y <- ggiraphExtra::getMapping(object$mapping, "y")
  }
  else {
    aes_x <- aes_x
    aes_y <- aes_y
  }
  data_x <- data[, c(aes_x)]
  data_y <- data[, c(aes_y)]
  if (annoManual == FALSE) {
    if (annoPos %in% c("top", "botomn")) {
      if (aesGroup == FALSE) {
        nPoints <- length(xPosition)
        xPos <- xPosition
        xmin <- xPos - segWidth/2
        xmax <- xPos + segWidth/2
      }
      else {
        groupInfo <- data %>% dplyr::select(.data[[aes_x]],
                                            .data[[aesGroName]]) %>% unique() %>% dplyr::select(.data[[aesGroName]]) %>%
          table() %>% data.frame()
        start <- c(1, groupInfo$Freq[1:(length(groupInfo$Freq) -
                                          1)]) %>% cumsum()
        end <- cumsum(groupInfo$Freq)
        xmin <- start - segWidth/2
        xmax <- end + segWidth/2
        nPoints <- length(start)
      }
      if (is.null(yPosition)) {
        if (is.numeric(data_y)) {
          if (annoPos == "top") {
            ymax <- max(data_y) + relSideDist * max(data_y)
            ymin <- ymax
          }
          else {
            ymin <- min(data_y) - relSideDist * max(data_y)
            ymax <- ymin
          }
        }
        else {
          if (annoPos == "top") {
            ymax <- length(unique(data_y)) + relSideDist *
              length(unique(data_y))
            ymin <- ymax
          }
          else {
            ymin <- -relSideDist * length(unique(data_y))
            ymax <- ymin
          }
        }
      }
      else {
        ymax <- yPosition[1]
        ymin <- yPosition[1]
      }
    }
    else if (annoPos %in% c("left", "right")) {
      if (aesGroup == FALSE) {
        nPoints <- length(yPosition)
        yPos <- yPosition
        ymin <- yPos - segWidth/2
        ymax <- yPos + segWidth/2
      }
      else {
        groupInfo <- data %>% dplyr::select(.data[[aes_y]],
                                            .data[[aesGroName]]) %>% unique() %>% dplyr::select(.data[[aesGroName]]) %>%
          table() %>% data.frame()
        start <- c(1, groupInfo$Freq[1:(length(groupInfo$Freq) -
                                          1)]) %>% cumsum()
        end <- cumsum(groupInfo$Freq)
        ymin <- start - segWidth/2
        ymax <- end + segWidth/2
        nPoints <- length(start)
      }
      if (is.null(xPosition)) {
        if (is.numeric(data_x)) {
          if (annoPos == "left") {
            xmin <- min(data_x) - relSideDist * max(data_x)
            xmax <- xmin
          }
          else {
            xmax <- max(data_x) + relSideDist * max(data_x)
            xmin <- xmax
          }
        }
        else {
          if (annoPos == "left") {
            xmin <- -relSideDist * length(unique(data_x))
            xmax <- xmin
          }
          else {
            xmax <- length(unique(data_x)) + relSideDist *
              length(unique(data_x))
            xmin <- xmax
          }
        }
      }
      else {
        xmin <- xPosition[1]
        xmax <- xPosition[1]
      }
    }
  }
  else {
    if (annoPos %in% c("top", "botomn")) {
      xmin <- xPosition[[1]] - segWidth/2
      xmax <- xPosition[[2]] + segWidth/2
      ymax <- yPosition[[1]]
      ymin <- yPosition[[1]]
    }
    else {
      xmin <- xPosition[[1]]
      xmax <- xPosition[[1]]
      ymin <- yPosition[[1]] - segWidth/2
      ymax <- yPosition[[2]] + segWidth/2
    }
    nPoints <- max(length(xmin), length(ymin))
  }
  annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf,
                                 ymin = -Inf, ymax = Inf, data) {
    ggplot2::layer(data = data, stat = StatIdentity, position = PositionIdentity,
                   geom = ggplot2::GeomCustomAnn, inherit.aes = TRUE,
                   params = list(grob = grob, xmin = xmin, xmax = xmax,
                                 ymin = ymin, ymax = ymax))
  }
  if (is.null(pCol)) {
    pCol <- useMyCol("stallion", n = nPoints)
  }
  else {
    pCol <- pCol
  }
  if (is.null(facetName)) {
    if (annoPos %in% c("top", "botomn")) {
      for (i in 1:nPoints) {
        object <- object + ggplot2::annotation_custom(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol[i],
                                                                                                fill = pCol[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                                alpha = alpha), arrow = mArrow), xmin = ggplot2::unit(xmin[i],
                                                                                                                                                      "native"), xmax = ggplot2::unit(xmax[i], "native"),
                                                      ymin = ggplot2::unit(ymin, "native"), ymax = ggplot2::unit(ymax,
                                                                                                                 "native"))
      }
    }
    else if (annoPos %in% c("left", "right")) {
      for (i in 1:nPoints) {
        object <- object + ggplot2::annotation_custom(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol[i],
                                                                                                fill = pCol[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                                alpha = alpha), arrow = mArrow), xmin = ggplot2::unit(xmin,
                                                                                                                                                      "native"), xmax = ggplot2::unit(xmax, "native"),
                                                      ymin = ggplot2::unit(ymin[i], "native"), ymax = ggplot2::unit(ymax[i],
                                                                                                                    "native"))
      }
    }
    else {
    }
  }
  else {
    facet_data <- data.frame(myFacetGrou)
    colnames(facet_data) <- facetName
    if (annoPos %in% c("top", "botomn")) {
      for (i in 1:nPoints) {
        object <- object + annotation_custom2(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol[i],
                                                                                        fill = pCol[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                        alpha = alpha), arrow = mArrow), data = facet_data,
                                              xmin = xmin[i], xmax = xmax[i], ymin = ymin,
                                              ymax = ymax)
      }
    }
    else if (annoPos %in% c("left", "right")) {
      for (i in 1:nPoints) {
        object <- object + annotation_custom2(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol[i],
                                                                                        fill = pCol[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                        alpha = alpha), arrow = mArrow), data = facet_data,
                                              xmin = xmin, xmax = xmax, ymin = ymin[i], ymax = ymax[i])
      }
    }
    else {
    }
  }
  if (addBranch == TRUE) {
    if (annoPos %in% c("top", "botomn")) {
      brXmin <- c(xmin, xmax)
      brXmax <- c(xmin, xmax)
      brYmin <- ymax + branRelSegLen * segWidth * branDirection
      brYmax <- ymax
    }
    else {
      brXmin <- xmax
      brXmax <- xmax + branRelSegLen * segWidth * branDirection
      brYmin <- c(ymin, ymax)
      brYmax <- c(ymin, ymax)
    }
    pCol2 <- rep(pCol, 2)
  }
  if (is.null(facetName)) {
    if (addBranch == TRUE & annoPos %in% c("top", "botomn")) {
      for (i in 1:(2 * nPoints)) {
        object <- object + ggplot2::annotation_custom(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol2[i],
                                                                                                fill = pCol2[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                                alpha = alpha), arrow = bArrow), xmin = ggplot2::unit(brXmin[i],
                                                                                                                                                      "native"), xmax = ggplot2::unit(brXmax[i],
                                                                                                                                                                                      "native"), ymin = ggplot2::unit(brYmin, "native"),
                                                      ymax = ggplot2::unit(brYmax, "native"))
      }
    }
    else if (addBranch == TRUE & annoPos %in% c("left", "right")) {
      for (i in 1:(2 * nPoints)) {
        object <- object + ggplot2::annotation_custom(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol2[i],
                                                                                                fill = pCol2[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                                alpha = alpha), arrow = bArrow), xmin = ggplot2::unit(brXmin,
                                                                                                                                                      "native"), xmax = ggplot2::unit(brXmax, "native"),
                                                      ymin = ggplot2::unit(brYmin[i], "native"),
                                                      ymax = ggplot2::unit(brYmax[i], "native"))
      }
    }
    else {
    }
  }
  else {
    if (addBranch == TRUE & annoPos %in% c("top", "botomn")) {
      for (i in 1:(2 * nPoints)) {
        object <- object + annotation_custom2(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol2[i],
                                                                                        fill = pCol2[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                        alpha = alpha), arrow = bArrow), data = facet_data,
                                              xmin = brXmin[i], xmax = brXmax[i], ymin = brYmin,
                                              ymax = brYmax)
      }
    }
    else if (addBranch == TRUE & annoPos %in% c("left", "right")) {
      for (i in 1:(2 * nPoints)) {
        object <- object + annotation_custom2(grob = grid::segmentsGrob(gp = grid::gpar(col = pCol2[i],
                                                                                        fill = pCol2[i], lty = lty, lwd = lwd, lineend = lineend,
                                                                                        alpha = alpha), arrow = bArrow), data = facet_data,
                                              xmin = brXmin, xmax = brXmax, ymin = brYmin[i],
                                              ymax = brYmax[i])
      }
    }
    else {
    }
  }
  if (is.null(textCol)) {
    textCol <- useMyCol("stallion", n = nPoints)
  }
  else {
    textCol <- textCol
  }
  if (aesGroup == FALSE) {
    textLabel <- textLabel
  }
  else {
    textLabel <- groupInfo[, 1]
  }
  if (is.null(facetName)) {
    if (addText == TRUE & annoPos %in% c("top", "botomn")) {
      for (i in 1:nPoints) {
        object <- object + ggplot2::annotation_custom(grob = grid::textGrob(gp = grid::gpar(col = textCol[i],
                                                                                            fontsize = textSize, fontfamily = fontfamily,
                                                                                            fontface = fontface), hjust = hjust, vjust = vjust,
                                                                            label = textLabel[i], check.overlap = T, just = "centre",
                                                                            rot = textRot), xmin = ggplot2::unit(xmin[i],
                                                                                                                 "native"), xmax = ggplot2::unit(xmax[i], "native"),
                                                      ymin = ggplot2::unit(ymin + textHVjust, "native"),
                                                      ymax = ggplot2::unit(ymax + textHVjust, "native"))
      }
    }
    else if (addText == TRUE & annoPos %in% c("left", "right")) {
      for (i in 1:nPoints) {
        object <- object + ggplot2::annotation_custom(grob = grid::textGrob(gp = grid::gpar(col = textCol[i],
                                                                                            fontsize = textSize, fontfamily = fontfamily,
                                                                                            fontface = fontface), hjust = hjust, vjust = vjust,
                                                                            label = textLabel[i], check.overlap = T, just = "centre",
                                                                            rot = textRot), xmin = ggplot2::unit(xmin +
                                                                                                                   textHVjust, "native"), xmax = ggplot2::unit(xmax +
                                                                                                                                                                 textHVjust, "native"), ymin = ggplot2::unit(ymin[i],
                                                                                                                                                                                                             "native"), ymax = ggplot2::unit(ymax[i], "native"))
      }
    }
    else {
    }
  }
  else {
    if (addText == TRUE & annoPos %in% c("top", "botomn")) {
      for (i in 1:nPoints) {
        object <- object + annotation_custom2(grob = grid::textGrob(gp = grid::gpar(col = textCol[i],
                                                                                    fontsize = textSize, fontfamily = fontfamily,
                                                                                    fontface = fontface), hjust = hjust, vjust = vjust,
                                                                    label = textLabel[i], check.overlap = T, just = "centre",
                                                                    rot = textRot), data = facet_data, xmin = xmin[i],
                                              xmax = xmax[i], ymin = ymin + textHVjust, ymax = ymax +
                                                textHVjust)
      }
    }
    else if (addText == TRUE & annoPos %in% c("left", "right")) {
      for (i in 1:nPoints) {
        object <- object + annotation_custom2(grob = grid::textGrob(gp = grid::gpar(col = textCol[i],
                                                                                    fontsize = textSize, fontfamily = fontfamily,
                                                                                    fontface = fontface), hjust = hjust, vjust = vjust,
                                                                    label = textLabel[i], check.overlap = T, just = "centre",
                                                                    rot = textRot), data = facet_data, xmin = xmin +
                                                textHVjust, xmax = xmin + textHVjust, ymin = ymin[i],
                                              ymax = ymax[i])
      }
    }
    else {
    }
  }
  print(object)
}
