#' @rdname getModel
#'
#' @description
#' A method to get or set the INSPEcT_model object within an INSPEcT object. This method
#' is particularly useful to get and set testing parameters of the INSPEcT_model object
#' within the INSPEcT object.
#' @param object An object of class INSPEcT
#' @return An object of class INSPEcT model
#' @seealso \code{\link{testingParams}}
#' @examples
#' data('nascentInspObj10', package='INSPEcT')
#' getModel(nascentInspObj10)
setMethod('getModel', 'INSPEcT', function(object) {
	return(object@model)
	})
#' @rdname getModel
#'
#' @param value An object of class INSPEcT model
setReplaceMethod('getModel', 'INSPEcT', function(object, value) {
	object@model <- value
	object
	})
