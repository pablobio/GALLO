#' Sub-function to search QTLs around candidate markers
#'
#' Takes a list of candidate markers and search for genes a determined interval
#' @param chr_list "Object with the chromosomes to be analyzed"
#' @param db_file Data frame with the information from .gff file
#' @param marker_file Data frame with the information from the candidate regions file
#'@param nThreads The number of threads to be used
#' @param int The interval in base pair
#' @name sub_qtl_markers
#' @importFrom dplyr do
#' @importFrom data.table setkey
#' @importFrom data.table key
#' @importFrom data.table as.data.table
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @keywords internal
#' @return A dataframe with the QTLs mapped within the specified intervals

sub_qtl_markers<-function(chr_list,db_file,marker_file,nThreads=NULL,int=0){
    foreach::foreach(i=seq_along(1:length(chr_list)),.combine="rbind")%dopar%{
        chr<-chr_list[i]
        tmp_qtl<-data.table::as.data.table(db_file[which(db_file$chr==chr),])
        tmp_markers<-data.table::as.data.table(marker_file[which(marker_file$CHR==chr),])
        tmp_markers$tmpBP1<-tmp_markers$BP-int
        tmp_markers$tmpBP2<-tmp_markers$BP+int
        data.table::setkey(tmp_markers, tmpBP1,tmpBP2)
        out<-data.table::foverlaps(tmp_qtl, tmp_markers, by.x = c("start_pos","end_pos"), by.y = data.table::key(tmp_markers),nomatch = 0)
        out[,-c("tmpBP1","tmpBP2")]
    }
}
