#' Overlapping between grouping factors
#' 
#' Takes a dataframe with a column of genes, QTLs (or other data) and a grouping column and create some matrices with the ovelapping information 
#' @param file A dataframe with the data and grouping factor. The columns "data" and "group" are mandatories
#' @return A list with three matrices: 1) A matrix with the number of overllaping data; 2) A matrix with the percentage of overlapping; 3) A matrix with the combination of the two previous one
#' @export
overlapping_among_groups<-function(file){
  file<-file
  out.matrix.N<-matrix(ncol=length(unique(file$group)), nrow=length(unique(file$group)), NA)
  
  out.matrix.perc<-matrix(ncol=length(unique(file$group)), nrow=length(unique(file$group)), NA)
  
  out.matrix.merged<-matrix(ncol=length(unique(file$group)), nrow=length(unique(file$group)), NA)
  
  trait<-unique(file$group)
  
  for(i in 1:length(trait)){
    for(k in 1:length(trait)){
      tmp_perc<-(round(length(which(file[which(file$group==trait[i]),"data"] %in% file[which(file$group==trait[k]),"data"]))/length(file[which(file$group==trait[i]),"data"]),2))
      tmp_N<-(round(length(which(file[which(file$group==trait[i]),"data"] %in% file[which(file$group==trait[k]),"data"])),2))
      out.matrix.N[i,k]<-tmp_N
      out.matrix.perc[i,k]<-tmp_perc
      out.matrix.merged[i,k]<-paste(tmp_N," ","(",tmp_perc,")",sep="")
    }
    
  }
  out.matrix.list<-list(N=out.matrix.N,percentage=out.matrix.perc,combined=out.matrix.merged)
  return(out.matrix.list)
}
