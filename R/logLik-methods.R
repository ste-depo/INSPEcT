#' @rdname logLik
#'
#' @description
#' This method is used to retrieve all the log likelihood ratio test results for all pairs tested for all genes.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param ... Additional arguments for the generic
#' @return A matrix of log likelihood test results for all the tested model comparisons
#' @examples
#' data('nascentInspObj10', package='INSPEcT')
#' logLik(nascentInspObj10)
setMethod('logLik', 'INSPEcT_model', function(object, ...) {
	t(sapply(object@ratesSpecs, function(x) 
		tryCatch(sapply(x, '[[', 'logLik'),
			error=function(e) c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA))
		))
	})

#' @rdname logLik
setMethod('logLik', 'INSPEcT', function(object, ...) {
	logLik(object@model)
	})
