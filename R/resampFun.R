#' Function to performn resampling from qtl output
#' 
#' Takes a list of candidate markers and search for genes a determined interval
#' @param sub.qtl temporary qtl file splitted by chromosome created by enrich_qtl function
#' @param i number of interactions to be used during the bootstrap analysis

resampFun<-function(sub.qtl,i){
  tmp.resamp<-sub.qtl[i,]
  n.tmp.qtl<-nrow(tmp.resamp[grep(tmp.class,tmp.resamp$QTL_type),])
}