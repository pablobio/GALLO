# GALLO: Genomic Annotation in Livestock for positional candidate LOci

The Genomic Annotation in Livestock for positional candidate LOci (GALLO) is an R package designed to provide a straightforward environment for gene and QTL annotation, as well as data integration from multiple data sources. The QTL enrichment analyses can additionally be performed directly by GALLO using the output obtained from the QTL annotation step. In addition, GALLO also provide a set of functions for graphical visualization for the annotation, comparison, integration and QTL enrichment results. Consequently, GALLO is a useful package for the annotation, identification of hidden pattern across datasets, data mining of previous reported associations, as well as the efficient scrutinization of the genetic architecture of complex traits.

## Requirements

*Depends:* R (>= 3.5.0)

*Imports:* circlize, DT, data.table, doParallel, dplyr, dynamicTreeCut, ggplot2, graphics, grDevices, foreach, lattice , parallel, RColorBrewer, rtracklayer, stats, stringr, unbalhaar, utils

License: GPL-3

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

1. _import_gff_gtf():_ Takes a .gft or .gff file and import into a dataframe

2. _find_genes_qtls_around_markers:_ Takes a dataframe with candidate markers and/or regions (haplotypes, windows, CNVs, etc) and search for genes or QTLs in a specified interval

3. _overlapping_among_groups:_ Takes a dataframe with a column for genes, QTLs (or any other data) and a grouping column and create matrices with the ovelapping information

4. _plot_overlapping:_ Takes the output from overlapping_amoung_groups function and creates a heatmap with the overlapping between groups

5. _plot_qtl_info:_ Takes the output from find_genes_qtls_around_markers and create plots for the frequency of each QTL type and trait

6. _qtl_enrich:_ Takes the output from find_genes_qtls_around_markers and perform a QTL enrichment analysis

7. _QTLenrich_plot:_ Takes the output from _find_genes_qtls_around_markers function and creates a heatmap with the overlapping between groups

8. _relationship_plot:_ Takes the output from find_genes_qtls_around_markers function and creates a chord plot with the relationship between groups

**A tutorial for GALLO usage can be found at:**

https://rpubs.com/pablo_bio/GALLO_vignette

## Contact

For more informations, suggestions, discussions, and bug reports, contact pfonseca@uoguelph.ca
