#' Sub-function to auto-stop the clusters created for parallel processes
#' @param cl The cluster created by the makePSOCKcluster function
#' @importFrom utils capture.output
#' @keywords internal
autoStopCluster <- function(cl) {
  stopifnot(inherits(cl, "cluster"))
  env <- new.env()
  env$cluster <- cl
  attr(cl, "gcMe") <- env
  reg.finalizer(env, function(e) {
    try(parallel::stopCluster(e$cluster), silent = TRUE)
  })
  cl
}
