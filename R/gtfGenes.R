#'A gtf example for gene annotation
#'
#'Data from the Ensembl comprasing the gene annoation for the bovine genome
#'
#'@docType data
#'
#'@usage data(gtfGenes)
#'
#'@format A data frame with 24616 rows and 8 variables:
#'\itemize{
#'  \item chr: Chromosome
#'  \item start_pos: Start position in the genome for each geme
#'  \item end_pos: End position in the genome for each gene
#'  \item width Gene length
#'  \item strand Strand which the gene is mapped (+ or -)
#'  \item gene_id Ensemble gene ID
#'  \item gene_name Gene symbol
#'  \item gene_biotype Gene biotype
#' }
#'
#'@keywords datasets
#'
#'@examples
#'data(gtfGenes)
#'data(QTLmarkers)
#'\donttest{gene.out <- find_genes_qtls_around_markers(db_file="gffQTLs",
#'marker_file=QTLmarkers,method="gene",
#'marker="snp",interval=100000)}
#'\donttest{head(gene.out)}
#' @keywords internal
"gtfGenes"
