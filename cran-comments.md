## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

This is a resubmission. In this version I have:

* Removed the cran-commments.md file from .Rbuildignore. Therefore, the file should be available. 

* Included the gff (inst/extdata/QTL_UMD_3.1.gff.txt) file needed to run the example of import_gff_gtf function. This file has around 40 Mb and can not be split for the example.

* Unwrapped all examples from \dontrun{} or \donttest{}. The qtl_enrich and relationship_plot functions are running in 6.1 and 6.3 seconds of elapsed time, respectively. Subset this examples in order to reduce the elapsed time will result in meaningless examples.

* Removed on.exit(par(oldpar)) from example plot_qtl_info.Rd and provided the suggested alterations.

* Removed on.exit(par(oldpar)) from example

In the previous review round I had:

* Rephrased the description section to avoid to start the section with "this package" or package name.

* Included the https link to the current Github repository where more information regarding the package is provided. The manuscript fully describing the methods and the package usage is currently waiting the CRAN approval to be submitted to GigaScience Journal (as suggested by the editor). 

* Replaced the \dontrun{} wrapping to \donttest{}, except for the function import_gff_gtf as the function depends of an external file which will not be provided during installation in order to reduce the package size.

* Included the option for the user choose include or not the messages during the analysis using if(verbose)cat(). Additionally, all the message() were replaced by cat()

* Make sure that options user`s par or working directory were changed during the code execution. 

* Make sure that all the examples, vigenettes, etc. are running sequentially in this version.
