#' Plot overlapping between data and grouping factors
#'
#' Takes the output from overlapping_among_groups function and creates a heatmap with the overlapping between groups
#' @param overlapping_matrix The object obtained in overlapping_amoung_groups function
#' @param nmatrix An interger from 1 to 3 indicating which matrix will be used to plot the overlapping, where: 1) A matrix with the number of overllaping data; 2) A matrix with the percentage of overlapping; 3) A matrix with the combination of the two previous one
#' @param ntext An interger from 1 to 3 indicating which matrix will be used as the text matrix for the heatmap, where: 1) A matrix with the number of overllaping data; 2) A matrix with the percentage of overlapping; 3) A matrix with the combination of the two previous one
#' @param group A vector with the size of groups. This vector will be plotted as row and column names in the heatmap
#' @param labelcex A numeric value indicating the size of the row and column labels
#' @return A heatmap with the overlapping between groups
#' @importFrom grDevices colorRampPalette
#' @examples
#' data(QTLmarkers)
#' data(gtfGenes)
#' genes.out <- find_genes_qtls_around_markers(db_file=gtfGenes,
#' marker_file=QTLmarkers,method="gene",
#' marker="snp",interval=100000, nThreads=NULL)
#' overlapping.out<-overlapping_among_groups(file=genes.out,x="Reference",y="gene_id")
#' plot_overlapping(overlapping.out,nmatrix=2,ntext=2,group=unique(genes.out$Reference))
#' @export
plot_overlapping<-function(overlapping_matrix,nmatrix,ntext,group,labelcex=1){
overlapping<-as.matrix(overlapping_matrix[[nmatrix]])
colnames(overlapping)<-group
rownames(overlapping)<-group

myPanel <- function(x, y, z, ...) {
lattice::panel.levelplot(x,y,z,...)
lattice::panel.text(x, y,  overlapping_matrix[[ntext]][cbind(x,y)])
}

x.scale <- list(cex=labelcex, alternating=1, col='black',rot=90)
y.scale <- list(cex=labelcex, alternating=1, col='black')

    if(nmatrix==1){
    my_palette <- colorRampPalette(c("white", "red"))(n = max(overlapping))
    lattice::levelplot(overlapping, panel=myPanel,xlab="",ylab="",col.regions=my_palette,at=seq(0,max(overlapping),1),scales=list(x=x.scale, y=y.scale,tck = c(1,0)))
    }else{
        my_palette <- colorRampPalette(c("white", "red"))(n = 1000)
        lattice::levelplot(overlapping, panel=myPanel,xlab="",ylab="",col.regions=my_palette,at=seq(0,1,0.01),scales=list(x=x.scale, y=y.scale,tck = c(1,0)))
    }
}
