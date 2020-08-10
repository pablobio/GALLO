#' Overlapping between grouping factors
#'
#' Takes a dataframe with a column of genes, QTLs (or other data) and a grouping column and create some matrices with the overlapping information
#' @param file A dataframe with the data and grouping factor
#' @param x The grouping factor to be compared
#' @param y The data to be compared among the levels of the grouping factor
#' @return A list with three matrices: 1) A matrix with the number of overlapping data; 2) A matrix with the percentage of overlapping; 3) A matrix with the combination of the two previous one
#' @examples
#' data(QTLmarkers)
#' data(gtfGenes)
#' genes.out <- find_genes_qtls_around_markers(db_file=gtfGenes,
#' marker_file=QTLmarkers,method="gene",
#' marker="snp",interval=100000, nThreads=NULL)
#'overlapping.out<-overlapping_among_groups(file=genes.out,x="Reference",
#'y="gene_id")
#' @export
overlapping_among_groups<-function(file,x,y){
out.matrix.N<-matrix(ncol=length(unique(file[,x])), nrow=length(unique(file[,x])), NA, dimnames = list(unique(file[,x]),unique(file[,x])))
out.matrix.perc<-matrix(ncol=length(unique(file[,x])), nrow=length(unique(file[,x])), NA, dimnames = list(unique(file[,x]),unique(file[,x])))
out.matrix.merged<-matrix(ncol=length(unique(file[,x])), nrow=length(unique(file[,x])), NA, dimnames = list(unique(file[,x]),unique(file[,x])))
trait<-unique(file[,x])
    for(i in seq_along(1:length(trait))){
        for(k in seq_along(1:length(trait))){
        tmp_perc<-(round(length(which(file[which(file[,x]==trait[i]),y] %in% file[which(file[,x]==trait[k]),y]))/length(file[which(file[,x]==trait[i]),y]),2))
        tmp_N<-(round(length(which(file[which(file[,x]==trait[i]),y] %in% file[which(file[,x]==trait[k]),y])),2))
        out.matrix.N[i,k]<-tmp_N
        out.matrix.perc[i,k]<-tmp_perc
        out.matrix.merged[i,k]<-paste(tmp_N," ","(",tmp_perc,")",sep="")
        }
    }
out.matrix.list<-list(N=out.matrix.N,percentage=out.matrix.perc,combined=out.matrix.merged)
return(out.matrix.list)
}
