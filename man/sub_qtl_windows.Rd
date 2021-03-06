% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sub_qtls_windows.R
\name{sub_qtl_windows}
\alias{sub_qtl_windows}
\title{Sub-function to search QTLs around candidate markers}
\usage{
sub_qtl_windows(chr_list, db_file, marker_file, nThreads = NULL, int = 0)
}
\arguments{
\item{chr_list}{"Object with the chromosomes to be analyzed"}

\item{db_file}{Data frame with the information from .gff file}

\item{marker_file}{Data frame with the information from the candidate regions file}

\item{nThreads}{The number of threads to be used}

\item{int}{The interval in base pair}
}
\value{
A dataframe with the QTLs mapped within the specified intervals
}
\description{
Takes a list of candidate markers and search for genes a determined interval
}
\keyword{internal}
