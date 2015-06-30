#' @rdname nCores
#'
#' @description
#' This method is a fast way to get and set the number of cores that the modeling part will 
#' use (method \code{\link{modelRates}}). The same task can be achieved using the more general
#' \code{\link{modelingParams}} method
#' @param object An object of class INSPEcT
#' @param ... Additional arguments for the generic
#' @return Number of cores to be used
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' nCores(mycerIds10)
#' ## use the more general modelingParams method
#' modelingParams(mycerIds10)$nCores
setMethod('nCores', 'INSPEcT', function(object, ...) { object@params$nCores })
#' @rdname nCores
#'
#' @param value A numeric value with the number of cores to be used
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' mycerIds10 <- removeModel(mycerIds10)
#' nCores(mycerIds10) <- detectCores()
#' ## use the more general modelingParams method
#' modelingParams(mycerIds10)$nCores <- detectCores()
setReplaceMethod('nCores', 
	signature(object='INSPEcT', value='numeric'), 
	function(object, value) { 
		object@params$nCores <- as.integer(value) 
		return(object)
		})
