#' Sub-function to allow multithreading
#'
#' Takes the number of cores to be used and check if it is possible to allocate the work in the in the respective cores
#' @param nThreads The number of threads to be used
#' @param int The interval in base pair
#' @name MultiCores
#' @importFrom parallel detectCores
#' @importFrom dynamicTreeCut printFlush
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @keywords internal
#' @return Registers the parallel work in the required cores

MultiCores<-function(nThreads){
nCores = parallel::detectCores()
    if (!is.null(nThreads)){
    if (!is.numeric(nThreads) || nThreads < 2){
        stop("nThreads must be numeric and at least 2.")
    }
    if (nThreads > nCores){
        printFlush(paste("Warning in number of threads: Requested number of threads is higher than number\n","of available processors (or cores).","It is recommended that the number of threads is no more than number\n","of available processors.\n"))
        doParallel::registerDoParallel(nThreads)
    }
    if (nThreads < nCores){
        doParallel::registerDoParallel(nThreads)
    }
  }
}
