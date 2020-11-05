#' Performs a QTL enrichment analysis based in a Bootstrap simulation for each QTL class using the QTL information per chromosome
#'
#' Takes the output from find_genes_qtls_around_markers and run a QTL enrichment analysis
#' @param qtl_file The output from find_genes_qtls_around_markers function
#' @param qtl_type A string indicating with QTL enrichment will be performed: "QTL_type" or "Name"
#' @param qtl.file.types A vector with the observed QTL classes
#' @param table.qtl.class An frequency table for the number of each QTL in each chromosome
#' @param qtl_db The QTL annotation database
#' @param search_qtl The column to perform the QTL searching in counting from the QTL annotation database
#' @param nThreads Number of threads for parallel processing
#' @return A data frame with the p-value for th enrichment result
#' @name sub_qtlEnrich_chrom
#' @importFrom stats p.adjust
#' @importFrom stats phyper
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @keywords internal
#'
sub_qtlEnrich_chom<-function(qtl_file,qtl_type,qtl.file.types,table.qtl.class,padj,qtl_db,search_qtl,nThreads){
    nCores<-parallel::detectCores()
    if (!is.null(nThreads)){
        cl<-doParallel::registerDoParallel(nThreads)
    }else{
        nThreads<- nCores
        cl<-doParallel::registerDoParallel(nThreads)
    }
    out.final<-foreach::foreach(k=qtl.file.types,.combine="rbind")%dopar%{
        tmp.qtl<-k
        table.qtl.class.tmp<-table.qtl.class[which(table.qtl.class$Var1==tmp.qtl),]
        out.chr<-foreach::foreach(j=table.qtl.class.tmp$Var2,.combine="rbind")%dopar%{
            tmp.chr<-j
            sub.qtl<-qtl_db[which(qtl_db$chr==tmp.chr),]
            tmp.qtl.file<-qtl_file[which(qtl_file$CHR==tmp.chr),]
            Total_annotated_QTLs<-nrow(tmp.qtl.file[!duplicated(tmp.qtl.file[,c("Name","QTL_ID","trait_ID")]),])
            tmp.qtl.file<-tmp.qtl.file[which(tmp.qtl.file[,qtl_type]==tmp.qtl),]
            n.qtls<-nrow(tmp.qtl.file[!duplicated(tmp.qtl.file[,c("Name","QTL_ID","trait_ID")]),])
            n.qtls.db<-nrow(sub.qtl[grep(pattern=tmp.qtl,x=sub.qtl[,search_qtl], fixed = TRUE),])
            pvalue<-phyper(n.qtls-1,n.qtls.db,(nrow(sub.qtl)-n.qtls.db),Total_annotated_QTLs, lower.tail = F)
            out.enrich.tmp<-data.frame(QTL=tmp.qtl,CHR=tmp.chr, N_QTLs=n.qtls,N_QTLs_db=n.qtls.db,Total_annotated_QTLs=Total_annotated_QTLs,Total_QTLs_db=nrow(sub.qtl),pvalue=pvalue)
        }
    }
    out.final$adj.pval<-p.adjust(out.final$pvalue,method=padj,n=nrow(out.final))
    return(out.final)
    return(cl)
}
