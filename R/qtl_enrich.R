#' Performs a QTL enrichment analysis based in a Bootstrap simulation for each QTL class
#'
#' Takes the output from find_genes_qtls_around_markers and run a QTL enrichment analysis
#' @param qtl_db The object obtained using the import_gff_gtf() function
#' @param qtl_file The output from find_genes_qtls_around_markers function
#' @param qtl_type A character indicating which type of enrichment will be performed. QTL_type indicates that the enrichment processes will be performed for the QTL classes, while Name indicates that the enrichment analysis will be performed for each trait individually
#' @param enrich_type A character indicating if the enrichment analysis will be performed for all the chromosomes ("genome") or for a subset of chromosomes ("chromosome). If the "genome" option is selected, the results reported are the merge of all chromosomes
#' @param chr.subset If enrich_type equal "chromosome", it is possible to define a subset of chromosomes to be analyzed. The default is equal NULL. Therefore, all the chromosomes will be analyzed
#' @param nThreads The number of threads to be used.
#' @param padj The algorithm for multiple testing correction to be adopted ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
#' @param  verbose Logical value defining if messages should of not be printed during the analysis (default=TRUE)
#' @details The simple bias of investigation for some traits (such as milk production related traits in the QTL database for cattle) may result in a larger proportion of records in the database. Consequently, the simple investigation of the proportion of each QTL type might not be totally useful. In order to reduce the impact of this bias, a QTL enrichment analysis can be performed. The QTL enrichment analysis performed by GALLO package is based in a hypergeometric test using the number of annoatted QTLs within the candidate regions and the total number of the same QTL in the QTL database.
#' @return A data frame with the p-value for the enrichment result
#' @name qtl_enrich
#' @importFrom dynamicTreeCut printFlush
#' @examples
#' \donttest{data(QTLmarkers)
#' data(gffQTLs)
#' out.qtls<-find_genes_qtls_around_markers(db_file=gffQTLs,
#' marker_file=QTLmarkers, method = "qtl",
#' marker = "snp", interval = 500000, nThreads = NULL)
#' out.enrich<-qtl_enrich(qtl_db=gffQTLs, qtl_file=out.qtls,
#' qtl_type = "Name", enrich_type = "chromosome",
#' chr.subset = NULL, padj = "fdr",nThreads = NULL)}
#' @export
qtl_enrich<-function(qtl_db,qtl_file,qtl_type=c("QTL_type","Name"),enrich_type=c("genome","chromosome"),chr.subset=NULL,nThreads=NULL,padj=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"), verbose=TRUE){
nThreads<-nThreads
    if(is.null(chr.subset)){
        chr.subset<-unique(qtl_file$CHR)
    }
qtl_file<-qtl_file[which(qtl_file$CHR%in%chr.subset),]
    if(qtl_type=="QTL_type"){
        qtl.file.types<-unique(qtl_file$QTL_type)
        search_qtl<-"QTL_type"
    }else{
        qtl.file.types<-unique(qtl_file$Name)
        search_qtl<-"extra_info"
    }
    if(verbose==TRUE){
        cat("Staring QTL enrichment analysis for QTL class")
        cat("\n")
    }
    if(enrich_type=="genome"){
        table.qtl.class<-as.data.frame(table(qtl_file[,qtl_type]))
        qtl.file.types<-as.character(unique(table.qtl.class$Var1))
        out.enrich<-sub_qtlEnrich_geno(qtl_file,qtl_type,qtl.file.types,table.qtl.class,padj,qtl_db,search_qtl,nThreads)
    }
    if(enrich_type=="chromosome"){
        table.qtl.class<-as.data.frame(table(qtl_file[,qtl_type],qtl_file$CHR))
        table.qtl.class<-table.qtl.class[which(table.qtl.class$Freq!=0),]
        qtl.file.types<-as.character(unique(table.qtl.class$Var1))
        out.enrich<-sub_qtlEnrich_chom(qtl_file,qtl_type,qtl.file.types,table.qtl.class,padj,qtl_db,search_qtl,nThreads)
    }
out.enrich$QTL<-gsub("_"," ", out.enrich$QTL)
    if(qtl_type=="Name"){
    out.enrich$QTL_type<-qtl_file[match(out.enrich$QTL,qtl_file[,qtl_type]),"QTL_type"]
    out.enrich$QTL_type<-gsub("_"," ", out.enrich$QTL_type)
    }
    if(verbose==TRUE){
        cat("\n")
        message("End of QTL enrichment analysis")
    }
return(out.enrich)
}
