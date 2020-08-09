#' @rdname logLik
#'
#' @description
#' This method is used to retrieve all the log likelihood ratio test results for all pairs tested for all genes.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param ... Additional arguments for the generic
#' @return A matrix of log likelihood test results for all the tested model comparisons
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' logLik(nascentInspObj10)
setMethod('logLik', 'INSPEcT_model', function(object, ...) {
	out <- logLik_internal(object)
	colnames(out) <- c("no-reg","s","p","d","sd","sp","pd","spd")
	return(out)
	})

#' @rdname logLik
setMethod('logLik', 'INSPEcT', function(object, ...) {
	checkINSPEcTObjectversion(object)
	logLik(object@model)
	})

logLik_internal <- function(object) {
	t(sapply(object@ratesSpecs, function(x) 
		tryCatch(sapply(x, '[[', 'logLik'),
						 error=function(e) c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA))
	))
}