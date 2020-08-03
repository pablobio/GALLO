#' Performs a QTL enrichment analysis based in a Bootstrap simulation for each QTL class using the QTL information across the whole genome
#'
#' Takes the output from find_genes_qtls_around_markers and run a QTL enrichment analysis
#' @param qtl_file The output from find_genes_qtls_around_markers function
#' @param qtl.file.types A vector with the observed QTL classes
#' @param qtl_type A string indicating with QTL enrichment will be performed: "QTL_type" or "Name"
#' @param table.qtl.class An frequency table for the number of each QTL in each chromosome
#' @param qtl_db The QTL annotation database
#' @param search_qtl The column to perform the QTL searching in counting from the QTL annotation database
#' @param nThreads Number of threads for parallel processing
#' @return A data frame with the p-value for th enrichment result
#' @name sub_qtlEnrich_geno
#' @importFrom stats p.adjust
#' @importFrom stats phyper
#' @keywords internal
#'
sub_qtlEnrich_geno<-function(qtl_file,qtl_type,qtl.file.types,table.qtl.class,padj,qtl_db,search_qtl,nThreads){
    MultiCores(nThreads)
    out.final<-foreach::foreach(k=qtl.file.types,.combine="rbind")%dopar%{
        tmp.qtl<-k
        tmp.qtl.file<-qtl_file[which(qtl_file[,qtl_type]==tmp.qtl),]
        n.qtls<-nrow(tmp.qtl.file[!duplicated(tmp.qtl.file[,c("Name","QTL_ID","trait_ID")]),])
        n.qtls.db<-nrow(qtl_db[grep(pattern=tmp.qtl,x=qtl_db[,search_qtl], fixed = TRUE),])
        Total_annotated_QTLs<-nrow(qtl_file[!duplicated(qtl_file[,c("Name","QTL_ID","trait_ID")]),])
        pvalue<-phyper(n.qtls-1,n.qtls.db,(nrow(qtl_db)-n.qtls.db),Total_annotated_QTLs, lower.tail = F)
        out.enrich.tmp<-data.frame(QTL=tmp.qtl, N_QTLs=n.qtls,N_QTLs_db=n.qtls.db,Total_annotated_QTLs=Total_annotated_QTLs,Total_QTLs_db=nrow(qtl_db),pvalue=pvalue)
    }
    out.final$adj.pval<-p.adjust(out.final$pvalue,method=padj,n=nrow(out.final))
    return(out.final)
}
