#' Function to performn resampling from qtl output
#'\alias resampFun


resampFun<-function(sub.qtl,i){
  tmp.resamp<-sub.qtl[i,]
  n.tmp.qtl<-nrow(tmp.resamp[grep(tmp.class,tmp.resamp$QTL_type),])
}
