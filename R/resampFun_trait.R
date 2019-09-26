#' Function to performn resampling from qtl output
#'\alias resampFun.trait


resampFun.trait<-function(sub.qtl,i){
  tmp.resamp<-sub.qtl[i,]
  n.tmp.qtl<-nrow(tmp.resamp[grep(pattern=tmp.class,x=tmp.resamp$extra_info),])
}
