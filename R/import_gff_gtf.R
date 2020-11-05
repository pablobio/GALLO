#' Import .gtf and .gff files to be used during gene and QTL annotation, respectively
#'
#' Takes a .gft or .gff file and import into a dataframe
#' @param db_file File with the gene mapping or QTL information. For gene mapping, a .gtf file  from Ensembl database must be used. For the QTL search, a .gff file from Animal QTlLdb must be used. Both files must use the same reference annotation used in the original study
#' @param file_type "gtf" or "gff"
#' @return A dataframe with the gtf or gtf content
#' @name import_gff_gtf
#' @importFrom utils read.delim
#' @examples
#' gffpath <- system.file("extdata", "example.gff", package="GALLO")
#'
#' qtl.inp <- import_gff_gtf(db_file=gffpath,file_type="gff")
#' @export

import_gff_gtf<-function(db_file,file_type){
    if (!(file.exists(db_file))){
    stop(paste("file ", db_file, " doesn't exists"))
    }else{
        if(file_type=="gtf"){
        gtf.file<-rtracklayer::import(db_file)
        gtf.file<-as.data.frame(gtf.file)
        gene<-gtf.file[which(gtf.file$type=="gene"),]
        gene<-gene[,c("seqnames","start","end","width","strand","gene_id","gene_name","gene_biotype")]
        colnames(gene)[c(1,2,3)]<-c("chr","start_pos","end_pos")
        gene$chr<-as.character(gene$chr)
        return(gene)
        }

        if(file_type=="gff"){
        qtl=read.delim(db_file, header=F, comment.char="#")
        qtl<-qtl[,c(seq_along(1:5),9)]
        names(qtl)<-c("chr", "database","QTL_type", "start_pos", "end_pos","extra_info")
        qtl$chr<-gsub("Chr.", "", qtl$chr)
        qtl<-qtl[!is.na(qtl$chr),]
        qtl<-qtl[!is.na(qtl$start_pos),]
        qtl<-qtl[!is.na(qtl$end_pos),]
        return(as.data.frame(qtl))
        }
    }
}
