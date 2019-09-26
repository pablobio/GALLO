# GALLO: Genomic Annotation in Livestock for positional candidate LOci

The Genomic Annotation in Livestock for positional candidate LOci (GALLO) is an R package designed to provide a straightforward environment for gene and QTL annotation, as well as data integration from multiple data sources. The QTL enrichment analyses can additionally be performed directly by GALLO using the output obtained from the QTL annotation step. In addition, GALLO also provide a set of functions for graphical visualization for the annotation, comparison, integration and QTL enrichment results. Consequently, GALLO is a useful package for the annotation, identification of hidden pattern across datasets, datamining of previous reported associations, as well as the efficient scrutinization of the genetic architecture of complex traits.

## Requirements

*Depends:* R (>= 3.5.0)

*Imports:* rtracklayer, data.table, doParallel, parallel, foreach, lattice, utils, graphics, dplyr, grDevices, boot, RColorBrewer, circlize, dynamicTreeCut, ggplot2, unbalhaar

License: MIT + file LICENSE

## Usage

- The package starts with the outputs from the most common high-throughput genomic association software such as PLINK, BLUPF90, DESeq2, etc. All that it is necessary in the input file is a column named "CHR", indicating the chromosome number and a column named "BP" with the chromosomal position in base pairs, when punctual positions are analyzed (i.e., SNPs). On the other hand, when chromosomal windows are evaluated, the column "CHR" is still mandatory, however, now it is necessary the presence of the columns "BP1" and "BP2". Additionally, the .gtf and .gff files must be provided for gene and QTL annotation respectivelly. 

*.gtf files for gene annotation can de found, for example, in ensembl FTP website: https://www.ensembl.org/info/data/ftp/index.html*

*.gtf files for QTL annotation can be found in Animal QTLdb: https://www.animalgenome.org/QTLdb/*

- The outputs from each function within GALLO can be used in the downstream functions. For example, the output from find_genes_and_qtls_around_markers() can be used as input of qtl_enrich(), relationship_plot(), etc.

**To install the package, the following command line can be use in R:**
```
#install.packages("devtools")
library(devtools)
install_github("pablobio/GALLO")
```

## Functions description

- find_genes_qtls_around_markers:	Search genes and QTLs around candidate regions
- overlapping_among_groups:	Overlapping between grouping factors
- plot_overlapping:	Plot overlapping between data and grouping factors
- plot_qtl_info:	Plot QTLs information from the find_genes_qtls_around_markers output
- QTLenrich_plot:	Plot enrichment results for QTL enrichment analysis
- qtl_enrich:	Performs a QTL enrichment analysis based in a Bootstrap simulation for each QTL class
- relationship_plot:	Plot relationship between data and grouping factors

**Annotation of QTLs overlapping genomic windows, a short example**
```
#Loading package
library(GALLO)

#Loading example dataset
data(QTLwindows)

#Performing QTL annotation (method="qtl") for genomic windows (marker="haplotype"), 
#using an interval of 100Kb upstream and downstream (interval=100000)

qtl.out <- find_genes_qtls_around_markers(db_file="QTL_db.gff",'marker_file=QTLwindows,
method="qtl",'marker="haplotypes",interval=100000)

head(qtl.out)
```

## Contact

For more informations, suggestions, dicussions, and bug reports, contact pfonseca@uoguelph.ca
