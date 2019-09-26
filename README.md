# GALLO: Genomic Annotation in Livestock for positional candidate LOci

The Genomic Annotation in Livestock for positional candidate LOci (GALLO) is an R package designed to provide a straightforward environment for gene and QTL annotation, as well as data integration from multiple data sources. The QTL enrichment analyses can additionally be performed directly by GALLO using the output obtained from the QTL annotation step. In addition, GALLO also provide a set of functions for graphical visualization for the annotation, comparison, integration and QTL enrichment results. Consequently, GALLO is a useful package for the annotation, identification of hidden pattern across datasets, datamining of previous reported associations, as well as the efficient scrutinization of the genetic architecture of complex traits.

## Usage

- The package starts with the outputs from the most common high-throughput genomic association software such as PLINK, BLUPF90, DESeq2, etc. All that it is necessary in the input file is a column named "CHR", indicating the chromosome number and a column named "BP" with the chromosomal position in base pairs, when punctual positions are analyzed (i.e., SNPs). On the other hand, when chromosomal windows are evaluated, the column "CHR" is still mandatory, however, now it is necessary the presence of the columns "BP1" and "BP2". Additionally, the .gtf and .gff files must be provided for gene and QTL annotation respectivelly. 

*.gtf files can de found, for example, in ensembl FTP website: https://www.ensembl.org/info/data/ftp/index.html*
*.gtf files for QTL annotation can be found in Animal QTLdb: https://www.animalgenome.org/QTLdb/

- The outputs from each function within GALLO can be used in the downstream functions. For example, the output from find_genes_and_qtls_around_markers() can be used as input of qtl_enrich(), relationship_plot(), etc.

**To install the package, the following command line can be use in R:**
```
library(devtools)
install_github("pablobio/GALLO")
```

