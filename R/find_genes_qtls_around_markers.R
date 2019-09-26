#' Search genes and QTLs around candidate regions
#'
#' Takes a list of candidate markers and or regions (haplotypes) and search for genes or QTLs in a determined interval
#' @param db_file file with the gene mapping or QTL information. For the gene mapping, you should use the .gtf file download from Ensembl data base. For the QTL search, you need to inform the .gff file that can be downloaded from Animal QTlLdb.
#' @param marker_file The file with the SNP or haplotype positions. Detail: For SNP files, you must have a column called “CHR” and a column called “BP” with the chromosome and base pair position, respectively. For the haplotype, you must have three columns: “CHR”, “BP1” and “BP2”. All the columns names are capitals.
#' @param method “gene” or “qtl”
#' @param marker "snp" or "haplotype"
#' @param interval The interval in base pair
#' @param nThreads Number of threads to be used
#' @return A dataframe with the genes or QTLs mapped within the specified intervals
#' @name find_genes_qtls_around_markers
#' @importFrom utils read.delim
#' @examples
#' data(QTLwindows)
#'\donttest{qtl.out <- find_genes_qtls_around_markers(db_file="QTL_db.gff",
#'marker_file=QTLwindows,method="qtl",
#'marker="haplotypes",interval=100000)}
#'\donttest{head(qtl.out)}
#' @export

find_genes_qtls_around_markers<-function(db_file,marker_file,method=c("gene","qtl"),marker=c("snp","haplotype"),interval=0,nThreads=NULL){
  options(stringsAsFactors = F)
  #Number of threads to be used
  interval=interval
  nThreads=nThreads
  #check method
  method <- match.arg(method)
  #Print method selected
  message(paste("You are using the method:", method, "with", marker))
  cat("\n")
  #Checking which kind of markers will be used
  if(marker=="snp"){
    #Creating gene data frame
    if (method=="gene"){
      #load file
      #Checking if gft file was imported
      if (!(file.exists(db_file))){
        stop(paste("file ", db_file, " doesn't exists"))
      }


      markers=marker_file

      #reading gtf file
      gtf.file<-rtracklayer::import(db_file)
      gtf.file<-as.data.frame(gtf.file)
      gene<-gtf.file[which(gtf.file$type=="gene"),]
      gene<-gene[,c("seqnames","start","end","width","strand","gene_id","gene_name","gene_biotype")]
      colnames(gene)[c(1,2,3)]<-c("chr","start_pos","end_pos")
      gene$chr<-as.character(gene$chr)
      output_genes<-NULL
      out_final<-NULL
      cat("\n")
      message(paste("Starting Gene searching using ", interval, " bp", " as interval", sep=""))

      chr_list<-unique(markers$CHR)
      #Sub-setting tables by chromosome and searching within intervals
      output.final<-sub_genes_markers(chr_list,gene,markers,nThreads=nThreads,int=interval)
    }else {#Running the QTL searching
      #load file
      #Checking if SNP file was imported
      if (file.exists(db_file)){
        qtl=read.delim(db_file, header=F, comment.char="#")
        qtl<-qtl[,c(1:5,9)]
        names(qtl)<-c("chr", "database","QTL_type", "start_pos", "end_pos","extra_info")
        qtl$chr<-gsub("Chr.", "", qtl$chr)
        qtl<-qtl[!is.na(qtl$chr),]
        qtl<-qtl[!is.na(qtl$start_pos),]
        qtl<-qtl[!is.na(qtl$end_pos),]
      }else{
        stop(paste("file", marker_file, "doesn't exists"))
      }

      #Checking if SNP or haplotype file was imported
      markers=marker_file

      message(paste("Starting QTL searching using ", interval, " bp", " as interval", sep=""))

      chr_list<-unique(markers$CHR)
      #Sub-setting tables by chromosome and searching within intervals
      output.final<-sub_qtl_markers(chr_list,qtl,markers,nThreads=nThreads,int=interval)

      #Splitting extra_info column
      cat("\n")
      message("Preparing output file for QTL annotation")
      cat("\n")
      output.final<-splitQTL_comment(output.final)
    }
  }else{#Runnuning the analysis for haplotypes

    #Creating gene data frame
    if (method=="gene"){
      #load file
      #Checking if gft file was imported

      if (!(file.exists(db_file))){
        stop(paste("file", db_file, "doesn't exists"))
      }

      #Checking if SNP or haplotype file was imported

      markers=marker_file


      #read gtf file
      gtf.file<-rtracklayer::import(db_file)
      gtf.file<-as.data.frame(gtf.file)
      gene<-gtf.file[which(gtf.file$type=="gene"),]
      gene<-gene[,c("seqnames","start","end","width","strand","gene_id","gene_name","gene_biotype")]
      colnames(gene)[c(1,2,3)]<-c("chr","start_pos","end_pos")
      gene$chr<-as.character(gene$chr)
      output_genes<-NULL
      cat("\n")
      message(paste("Starting Gene searching using ", interval, " bp", " as interval", sep=""))

      chr_list<-unique(markers$CHR)
      #Sub-setting tables by chromosome and searching within intervals
      output.final<-sub_genes_windows(chr_list,gene,markers,nThreads=nThreads,int=interval)
    }else{
      #Running the QTL searching
      #load file
      #Checking if SNP file was imported
      if (file.exists(db_file)){
        qtl=read.delim(db_file, header=F, comment.char="#")
        qtl<-qtl[,c(1:5,9)]
        names(qtl)<-c("chr", "database","QTL_type", "start_pos", "end_pos","extra_info")
        qtl$chr<-gsub("Chr.", "", qtl$chr)
        qtl<-qtl[!is.na(qtl$chr),]
        qtl<-qtl[!is.na(qtl$start_pos),]
        qtl<-qtl[!is.na(qtl$end_pos),]
      }else{
        stop(paste("file", marker_file, "doesn't exists"))
      }

      #Checking if SNP or haplotype file was imported
      markers=marker_file

      message(paste("Starting QTL searching using ", interval, " bp", " as interval", sep=""))

      chr_list<-unique(markers$CHR)
      #Sub-setting tables by chromosome and searching within intervals
      output.final<-sub_qtl_windows(chr_list,qtl,markers,nThreads=nThreads,int=interval)

      #Splitting extra_info column
      cat("\n")
      message("Preparing output file for QTL annotation")
      cat("\n")
      output.final<-splitQTL_comment(output.final)
    }
  }
  return(as.data.frame(output.final))
}
