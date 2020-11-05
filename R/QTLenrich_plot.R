#' Plot enrichment results for QTL enrichment analysis
#'
#' Takes the output from qtl_enrich function and creates a bubble plot with enrichment results
#' @param qtl_enrich The output from qtl_enrich function
#' @param x Id column to be used from the qtl_enrich output
#' @param pval P-value to be used in the plot. The name informed to this argument must match the p-value column name in the enrichment table
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
#' @return A plot with the QTL enrichment results
#' @examples
#' \donttest{data(QTLmarkers)
#' data(gffQTLs)
#' out.qtls<-find_genes_qtls_around_markers(db_file=gffQTLs,
#' marker_file=QTLmarkers, method = "qtl",
#' marker = "snp", interval = 500000,
#' nThreads = 1)
#'
#' out.enrich<-qtl_enrich(qtl_db=gffQTLs,
#' qtl_file=out.qtls, qtl_type = "Name",
#' enrich_type = "genome", chr.subset = NULL, padj = "fdr",nThreads = 1)
#'
#' out.enrich.filtered<-out.enrich[which(out.enrich$adj.pval<0.05),]
#' QTLenrich_plot(out.enrich.filtered, x="QTL", pval="adj.pval")}
#' @export

QTLenrich_plot<-function(qtl_enrich,x,pval){
  Pvalue<-qtl_enrich[,pval]
  breaks.qtl<-qtl_enrich[,"N_QTLs"]
  breaks.qtl<-as.numeric(breaks.qtl)
  breaks.qtl<-unique(breaks.qtl)
  breaks.qtl<-breaks.qtl[order(breaks.qtl, decreasing = F)]
  ggplot(qtl_enrich,aes(x=(qtl_enrich[,"N_QTLs"]/qtl_enrich[,"N_QTLs_db"]), y=qtl_enrich[,x], size=qtl_enrich[,"N_QTLs"], color=-log10(Pvalue))) + geom_point(alpha=0.8) + scale_color_gradient2(low="white", high="red") + scale_size(limits = c(min(breaks.qtl), max(breaks.qtl)), breaks = breaks.qtl, name="Number of QTLs") + theme(axis.title.x=element_text(size=15),axis.text.x=element_text(size=15),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=15), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + labs(x="Richness factor") + labs(size="Number of QTLs",col="-log10(P-value)")
}
