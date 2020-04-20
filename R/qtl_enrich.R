#' Performs a QTL enrichment analysis based in a Bootstrap simulation for each QTL class
#'
#' Takes the output from find_genes_qtls_around_markers and run a QTL enrichment analysis
#' @param qtl_db The .gff file that can be downloaded from Animal QTlLdb
#' @param qtl_file The output from find_genes_qtls_around_markers function
#' @param qtl_type A character indicating which type of enrichment will be performed. QTL_type indicates that the enrichment processes will be performed for the QTL classes, while trait indicates that the enrichment analysis will be performed for each trait individually.
#' @param enrich_type A character indicating if the enrichment analysis will be performed for all the chromosomes ("genome") or for a subset of chromosomes ("chromosome). If the "genome" option is selected, the results reported are the merge of all chromosomes.
#' @param chr.subset If enrich_type equal "chromosome", it is possible to define a subset of chromosomes to be analyzed. The default is equal NULL. Therefore, all the chromosomes will be analyzed.
#' @param n.it number of iterations for the bootstrap simulation.
#' @param nThreads The number of threads used.
#' @param padj The alogorithm for multiple testing correction to be adopted ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none").
#' @param parallel The type of parallel operation to be used.
#' @return A data frame with the p-value for th enrichment result
#' @importFrom utils read.delim
#' @importFrom stats p.adjust
#' @importFrom stats sd
#' @importFrom parallel detectCores
#' @importFrom dynamicTreeCut printFlush
#' @importFrom doParallel registerDoParallel
#' @importFrom boot boot
#' @examples
#' data(QTLwindows)
#'\donttest{qtl.out <- find_genes_qtls_around_markers(db_file="QTL_db.gff",
#'marker_file=QTLwindows,method="qtl",
#'marker="haplotypes",interval=100000)}
#'\donttest{qtl.enrich.out<-qtl_enrich(qtl_db="QTL_db.gff",
#'qtl_file=qtl.out,qtl_type="QTL_type",
#'enrich_type="genome",chr.subset=NULL,
#'n.it=1000,padj="fdr")}
#'\donttest{head(qtl.enrich.out)}
#' @export
qtl_enrich<-function(qtl_db,qtl_file,qtl_type=c("QTL_type","trait"),enrich_type=c("genome","chromosome"),chr.subset=NULL,n.it=NULL,nThreads=NULL,padj=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),parallel=c("no", "multicore", "snow")){
  nCores = detectCores()
  if (is.null(nThreads)) {
    if (nCores < 4)
      nThreads = nCores
    else nThreads = nCores - 1
    registerDoParallel(nThreads)
  }

  if (!is.numeric(nThreads) || nThreads < 2){
    stop("nThreads must be numeric and at least 2.")
  }
  if (nThreads > nCores){
    printFlush(paste("Warning in number of threads: Requested number of threads is higher than number\n","of available processors (or cores).","It is recommended that the number of threads is no more than number\n","of available processors.\n"))
    registerDoParallel(nThreads)
  }

  if (file.exists(qtl_db)){
    qtl=read.delim(qtl_db, header=F, comment.char="#",stringsAsFactors = F)
    qtl<-qtl[,c(1:5,9)]
    names(qtl)<-c("chr", "database","QTL_type", "start_pos", "end_pos","extra_info")
    qtl$chr<-gsub("Chr.", "", qtl$chr)
    qtl<-qtl[!is.na(qtl$chr),]
    qtl<-qtl[!is.na(qtl$start_pos),]
    qtl<-qtl[!is.na(qtl$end_pos),]
  }else{
    stop(paste("file", qtl_db, "doesn't exists"))
  }


  if(qtl_type=="QTL_type"){
    if(enrich_type=="genome"){
      message("Staring QTL enrichment analysis for QTL class")
      qtl.file.types<-unique(qtl_file$QTL_type)
      table.qtl<-as.data.frame(table(qtl_file$CHR))

      n.qtls<-NULL
      Average_exp<-NULL
      sd_exp<-NULL
      out.enrich<-foreach::foreach(k=1:length(qtl.file.types),.combine="rbind")%dopar%{
        tmp.class<<-qtl.file.types[k]
        table.qtl.class<-as.data.frame(table(qtl_file[which(qtl_file$QTL_type==tmp.class),"chr"]))
        n.qtls<-nrow(qtl_file[which(qtl_file$QTL_type==tmp.class),])

        out.chr<-foreach::foreach(j=1:length(table.qtl$Var1),.combine="cbind")%dopar%{
          sub.qtl<-qtl[which(qtl$chr==table.qtl$Var1[j]),]
          b <- boot(sub.qtl, resampFun, R=n.it, parallel=parallel,ncpus=nThreads,sim = "parametric",ran.gen=function(sub.qtl,p) sub.qtl[sample(1:nrow(sub.qtl), table.qtl$Freq[j],replace = T), ])
          b$t
        }
        out.int<-rowSums(out.chr)
        hyp.fert<-(sum(table.qtl.class$Freq))-out.int
        pvalue<-mean(hyp.fert>0)
        pvalue[pvalue==0]<-1
        Average_exp<-mean(out.int)
        sd_exp<-sd(out.int)
        data.frame(QTL_type=tmp.class,Number_QTLs=n.qtls,Average_exp=Average_exp,sd_exp=sd_exp,p_value=pvalue)
      }
      out.enrich$QTL_type<-gsub("_"," ", out.enrich$QTL_type)
      out.enrich$adjpval<-p.adjust(out.enrich$p_value,method=padj,n=length(out.enrich$p_value))
    }
    if(enrich_type=="chromosome"){
      message("Staring QTL enrichment analysis for QTL class")

      if(is.null(chr.subset)){
        chr.subset<-unique(qtl_file$CHR)
      }

      qtl_file<-qtl_file[which(qtl_file$CHR%in%chr.subset),]
      qtl.file.types<-unique(qtl_file$QTL_type)
      table.qtl<-as.data.frame(table(qtl_file$CHR))
      table.qtl$Var1<-as.character(table.qtl$Var1)

      n.qtls<-NULL
      Average_exp<-NULL
      sd_exp<-NULL
      out.enrich<-foreach::foreach(k=1:length(qtl.file.types),.combine="rbind")%dopar%{
        tmp.class<<-qtl.file.types[k]
        table.qtl.class<-as.data.frame(table(qtl_file[which(qtl_file$QTL_type==tmp.class),"chr"]))

        fake.table<-data.frame(Var1=chr.subset,Freq=rep(0,length(chr.subset)))
        fake.table[match(table.qtl.class$Var1,fake.table$Var1),"Freq"]<-table.qtl.class$Freq

        table.qtl.class<-fake.table
        table.qtl.class$Var1<-as.character(table.qtl.class$Var1)
        out.chr<-foreach::foreach(j=table.qtl$Var1,.combine="cbind")%dopar%{
          tmp.chr<-j
          sub.qtl<-qtl[which(qtl$chr==tmp.chr),]
          b <- boot(sub.qtl, resampFun, R=n.it, parallel=parallel,ncpus=nThreads,sim = "parametric",ran.gen=function(sub.qtl,p) sub.qtl[sample(1:nrow(sub.qtl), table.qtl[which(table.qtl$Var1==tmp.chr),"Freq"],replace = T), ])
          b$t
        }
        hyp.fert<-t(apply(out.chr,1,"-",table.qtl.class[match(table.qtl$Var1,table.qtl.class$Var1),"Freq"]))
        pvalue<-colSums(hyp.fert>0)/n.it
        pvalue[pvalue==0]<-1
        Average_exp<-colMeans(out.chr)
        sd_exp<-apply(out.chr, 2, sd)
        data.frame(QTL_type=rep(tmp.class,length(table.qtl.class$Var1)),Chr=table.qtl.class$Var1,Number_QTLs=table.qtl.class[match(table.qtl$Var1,table.qtl.class$Var1),"Freq"],Average_exp=Average_exp,sd_exp=sd_exp,p_value=pvalue)
      }
      out.enrich<-out.enrich[which(out.enrich$Number_QTLs!=0),]
      out.enrich$QTL_type<-gsub("_"," ", out.enrich$QTL_type)
      out.enrich$adjpval<-p.adjust(out.enrich$p_value,method=padj,n=length(out.enrich$p_value))
      qtl_file$coord_qtl<-paste(qtl_file$QTL_type,"_",qtl_file$CHR,sep="")
      out.enrich$coord_qtl<-paste(out.enrich$QTL_type,"_",out.enrich$Chr,sep="")
      out.enrich<-out.enrich[which(out.enrich$coord_qtl%in%qtl_file$coord_qtl),]
      out.enrich<-out.enrich[,-which(names(out.enrich)=="coord_qtl")]
    }
  }

  if(qtl_type=="trait"){
    if(enrich_type=="genome"){
      message("Staring QTL enrichment analysis for trait")
      qtl.file.types<-unique(qtl_file$Name)
      trait.ID<-unique(qtl_file$trait_ID)
      table.qtl<-as.data.frame(table(qtl_file$CHR))

      n.qtls<-NULL
      Average_exp<-NULL
      sd_exp<-NULL
      out.enrich<-foreach::foreach(k=1:length(qtl.file.types),.combine="rbind")%dopar%{
        tmp.class<<-qtl.file.types[k]
        table.qtl.class<-as.data.frame(table(qtl_file[which(qtl_file$Name==tmp.class),"chr"]))
        n.qtls<-nrow(qtl_file[which(qtl_file$Name==tmp.class),])

        out.chr<-foreach::foreach(j=1:length(table.qtl$Var1),.combine="cbind")%dopar%{
          sub.qtl<-qtl[which(qtl$chr==table.qtl$Var1[j]),]
          b <- boot(sub.qtl, resampFun.trait, R=n.it, parallel=parallel,ncpus=nThreads,sim = "parametric",ran.gen=function(sub.qtl,p) sub.qtl[sample(1:nrow(sub.qtl), table.qtl$Freq[j],replace = T), ])
          b$t
        }
        out.int<-rowSums(out.chr)
        hyp.fert<-(sum(table.qtl.class$Freq))-out.int
        pvalue<-mean(hyp.fert>0)
        pvalue[pvalue==0]<-1
        Average_exp<-mean(out.int)
        sd_exp<-sd(out.int)
        data.frame(QTL_type=tmp.class,Number_QTLs=n.qtls,Average_exp=Average_exp,sd_exp=sd_exp,p_value=pvalue)
      }
      out.enrich$QTL_type<-gsub("_"," ", out.enrich$QTL_type)
      out.enrich$adjpval<-p.adjust(out.enrich$p_value,method=padj,n=length(out.enrich$p_value))
    }
    if(enrich_type=="chromosome"){
      message("Staring QTL enrichment analysis for trait")

      if(is.null(chr.subset)){
        chr.subset<-unique(qtl_file$CHR)
      }
      qtl_file<-qtl_file[which(qtl_file$CHR%in%chr.subset),]
      qtl.file.types<-unique(qtl_file$Name)
      trait.ID<-unique(qtl_file$trait_ID)
      table.qtl<-as.data.frame(table(qtl_file$CHR))
      table.qtl$Var1<-as.character(table.qtl$Var1)

      n.qtls<-NULL
      Average_exp<-NULL
      sd_exp<-NULL
      out.enrich<-foreach::foreach(k=1:length(qtl.file.types),.combine="rbind")%dopar%{
        tmp.class<<-qtl.file.types[k]
        table.qtl.class<-as.data.frame(table(qtl_file[which(qtl_file$Name==tmp.class),"chr"]))

        fake.table<-data.frame(Var1=chr.subset,Freq=rep(0,length(chr.subset)))
        fake.table[match(table.qtl.class$Var1,fake.table$Var1),"Freq"]<-table.qtl.class$Freq
        table.qtl.class<-fake.table

        out.chr<-foreach::foreach(j=table.qtl$Var1,.combine="cbind")%dopar%{
          tmp.chr<-j
          sub.qtl<-qtl[which(qtl$chr==tmp.chr),]
          b <- boot(sub.qtl, resampFun.trait, R=n.it, parallel=parallel,ncpus=nThreads,sim = "parametric",ran.gen=function(sub.qtl,p) sub.qtl[sample(1:nrow(sub.qtl), table.qtl[which(table.qtl$Var1==tmp.chr),"Freq"],replace = T), ])
          b$t
        }
        hyp.fert<-t(apply(out.chr,1,"-",table.qtl.class[match(table.qtl$Var1,table.qtl.class$Var1),"Freq"]))
        pvalue<-colSums(hyp.fert>0)/n.it
        pvalue[pvalue==0]<-1
        Average_exp<-colMeans(out.chr)
        sd_exp<-apply(out.chr, 2, sd)
        data.frame(QTL_type=rep(tmp.class,length(unique(table.qtl.class$Var1))),Chr=unique(table.qtl.class$Var1),Number_QTLs=table.qtl.class[match(table.qtl$Var1,table.qtl.class$Var1),"Freq"],Average_exp=Average_exp,sd_exp=sd_exp,p_value=pvalue)
      }
      out.enrich<-out.enrich[which(out.enrich$Number_QTLs!=0),]
      out.enrich$QTL_type<-gsub("_"," ", out.enrich$QTL_type)
      out.enrich$adj.pval<-p.adjust(out.enrich$p_value,method=padj,n=length(out.enrich$p_value))
      qtl_file$coord_qtl<-paste(qtl_file$QTL_type,"_",qtl_file$CHR,sep="")
      out.enrich$coord_qtl<-paste(out.enrich$QTL_type,"_",out.enrich$Chr,sep="")
      out.enrich<-out.enrich[which(out.enrich$coord_qtl%in%qtl_file$coord_qtl),]
      out.enrich<-out.enrich[,-which(names(out.enrich)=="coord_qtl")]
    }
  }
  cat("\n")
  message("End of QTL enrichment analysis")
  cat("\n")
  out.enrich
}

