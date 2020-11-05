#' Sub-function to split comment column from QTL output
#'
#' Takes a list of candidate markers and search for genes a determined interval
#' @param output_qtls Output from QTL annotation
#' @return A dataframe with the extra_info column content, from the gff file, broken in several additional columns
#' @importFrom stringr str_split_fixed
#' @keywords internal

splitQTL_comment<-function(output.final){
    output_qtls<-as.data.frame(output.final)
    output_qtls$QTL_type<-as.character(output_qtls$QTL_type)
    tmp.type<-vapply(strsplit(as.character(output_qtls[grep("_Association",output_qtls$QTL_type),"QTL_type"]), "_Association"),'[',1,FUN.VALUE=character(1))
    if(length(tmp.type)!=0){
        output_qtls[grep("_Association",output_qtls$QTL_type),"Association_type"]<-"Association"
        output_qtls[grep("_Association",output_qtls$QTL_type),"QTL_type"]<-tmp.type
    }
    tmp.type<-vapply(strsplit(as.character(output_qtls[grep("_QTL",output_qtls$QTL_type),"QTL_type"]), "_QTL"), '[', 1,FUN.VALUE=character(1))
    if(length(tmp.type)!=0){
        output_qtls[grep("_QTL",output_qtls$QTL_type),"Association_type"]<-"QTL"
        output_qtls[grep("_QTL",output_qtls$QTL_type),"QTL_type"]<-tmp.type
    }
    tmp.type<-vapply(strsplit(as.character(output_qtls[grep("_Mendelian",output_qtls$QTL_type),"QTL_type"]), "_Mendelian"), '[', 1,FUN.VALUE=character(1))
    if(length(tmp.type)!=0){
        output_qtls[grep("_Mendelian",output_qtls$QTL_type),"Association_type"]<-"Mendelian"
        output_qtls[grep("_Mendelian",output_qtls$QTL_type),"QTL_type"]<-tmp.type
    }
    output_qtls$extra_info<-iconv(output_qtls$extra_info, "latin1", "ASCII", "")
    tmp.split<-stringr::str_split_fixed(output_qtls$extra_info, ";",n=Inf)
    output_qtls[,"QTL_ID"]<-stringr::str_split_fixed(tmp.split[,1], "=",n=2)[,2]
    output_qtls[,"trait_ID"]<-stringr::str_split_fixed(tmp.split[,6], "=",n=2)[,2]
    output_qtls[,"breed"]<-stringr::str_split_fixed(tmp.split[,8], "=",n=2)[,2]
    output_qtls[,"Name"]<-stringr::str_split_fixed(tmp.split[,2], "=",n=2)[,2]
    output_qtls[,"Abbrev"]<-stringr::str_split_fixed(tmp.split[,3], "=",n=2)[,2]
    output_qtls[,"Model"]<-stringr::str_split_fixed(tmp.split[,12], "=",n=2)[,2]
    output_qtls[,"Test_Base"]<-stringr::str_split_fixed(tmp.split[,13], "=",n=2)[,2]
    output_qtls[,"pubmed_id"]<-stringr::str_split_fixed(tmp.split[,4], "=",n=2)[,2]
    output_qtls[,"p_value"]<-stringr::str_split_fixed(tmp.split[,15], "=",n=2)[,2]
    output_qtls[,"bayes_value"]<-stringr::str_split_fixed(tmp.split[,14], "=",n=2)[,2]
    output_qtls[,"Flank_Markers"]<-stringr::str_split_fixed(tmp.split[,7], "=",n=2)[,2]
    return(output_qtls[,-which(colnames(output_qtls)%in%"extra_info")])
}
