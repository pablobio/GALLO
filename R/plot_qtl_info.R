#' Plot QTLs information from the find_genes_qtls_around_markers output
#'
#' Takes the output from find_genes_qtls_around_markers and create plots for the frequency of each QTL type and trait
#' @param qtl_file The output from find_genes_qtls_around_markers function
#' @param qtl_plot "qtl_type" or"qtl_name"
#' @param n Number of QTLs to be plotted when the qtl_name option is selected
#' @param qtl_class Class of QTLs to be plotted when the qtl_name option is selected
#' @param horiz The legend of the pie plot for the qtl_type should be plotted vertically or horizontally. The default is FALSE. Therefore, the legend is plotted vertically.
#' @param ... Arguments to be passed to/from other methods. For the default method these can include further arguments (such as axes, asp and main) and graphical parameters (see par) which are passed to plot.window(), title() and axis.
#' @return A plot with the requested information
#' @importFrom graphics pie
#' @importFrom graphics legend
#' @importFrom graphics barplot
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' data(QTLmarkers)
#' data(gffQTLs)
#'
#' out.qtls<-find_genes_qtls_around_markers(db_file=gffQTLs,
#' marker_file=QTLmarkers, method = "qtl",
#' marker = "snp", interval = 500000,
#' nThreads = 1)
#'
#' plot_qtl_info(out.qtls, qtl_plot = "qtl_type", cex=2)
#' @export
plot_qtl_info<-function(qtl_file,qtl_plot=c("qtl_type","qtl_name"), n="all",qtl_class=NULL,horiz=FALSE,...){
  #check method
  qtl_plot <- match.arg(qtl_plot)

  qtl_file=qtl_file
  if(qtl_plot=="qtl_type"){
    qtl_file<-qtl_file[order(qtl_file$QTL_type),]
    data.prop<-as.data.frame(table(qtl_file$QTL_type))
    data.prop$Var1<-gsub("_"," ",data.prop$Var1)
    labels.prop<-c(round(((data.prop[,2])/sum(data.prop[,2]))*100,2))
    labels.comp<-paste(labels.prop, "%", sep="")
    col.pallet<-RColorBrewer::brewer.pal(n=nrow(data.prop), name="Set3")

    pie(x=data.prop[,2], labels=labels.comp, col=col.pallet, ...)
    legend(x=-2.2,y=0.2, pch=15, col=col.pallet, legend=unique(data.prop$Var1),bty="n",horiz = horiz, xpd = TRUE,inset = c(0,0),y.intersp=1.2,xjust=0,yjust=0, ...)

  }

  if(qtl_plot=="qtl_name"){
    data.prop<-as.data.frame(table(qtl_file$Name))
    data.prop$percentage<-c(round(((data.prop[,2])/sum(data.prop[,2]))*100,2))
    data.prop<-data.prop[order(data.prop$percentage, decreasing=TRUE),]
    if(length(which(qtl_class %in% "all"))!=0){

      if(n=="all"){
        n_final<-nrow(data.prop)
      }else{
        n_final<-n
      }
      barplot(data.prop$percentage[seq_along(1:n_final)], names.arg=data.prop$Var1[seq_along(1:n_final)],horiz=TRUE, xlab="QTL Names (%)", las=1, xlim=c(0,round((max(data.prop$percentage[seq_along(1:n_final)])+0.5),0)), ...)

    }

    if(length(which(qtl_class %in% "all"))==0){
      data.prop<-as.data.frame(table(qtl_file$Name))
      qtl_file<-qtl_file[-duplicated(qtl_file[,c("QTL_type","Name")]),]
      data.prop$qtl_type<-qtl_file[match(as.character(data.prop[,1]),as.character(qtl_file$Name)),"QTL_type"]
      data.prop$percentage<-c(round(((data.prop[,2])/sum(data.prop[,2]))*100,2))
      data.prop<-data.prop[order(data.prop$percentage, decreasing=TRUE),]
      data.prop<-data.prop[which(data.prop$qtl_type %in% qtl_class),]

      if(n=="all"){
        n_final<-nrow(data.prop)
      }else{
        n_final<-n
      }

      barplot(data.prop$percentage[seq_along(1:n_final)], names.arg=data.prop$Var1[seq_along(1:n_final)],horiz=TRUE, xlab="QTL Names (%)", las=1, xlim=c(0,round((max(data.prop$percentage[seq_along(1:n_final)])+0.5),0)),...)

    }
  }
}
