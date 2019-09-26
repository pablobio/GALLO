#' Plot enrichment results for QTL enrichment analysis 
#' 
#' Takes the output from overlapping_amoung_groups function and creates a heatmap with the overlapping between groups
#' @param qtl_enrich The output from qtl_enrich function
#' @param x Id column to be used from the qtl_enrich output
#' @param pval P-value to be used in the plot. If "p_value" informed, a non-adjusted pvalue will be plot. If "padj" informed, the adjusted p-value from the qtl enrichment analysis will be plotted.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 theme_gray
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @export

QTLenrich_plot<-function(qtl_enrich,x,pval){
  pvalue<-qtl_enrich[,pval]
  ggplot(qtl_enrich,aes(x=(qtl_enrich[,"Number_QTLs"]/qtl_enrich[,"Average_exp"]), y=qtl_enrich[,x], size=qtl_enrich[,"Number_QTLs"], color=-log10(pvalue))) +
    geom_point(alpha=0.8) + scale_color_gradient2(low="white", high="red") +
    scale_size(range = c(.1, 24), name="Number of QTLs") + theme_gray(base_size = 14) +  theme(axis.title.x=element_text(size=15),axis.text.x=element_text(size=15),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=15)) + labs(x="Richness factor")
}