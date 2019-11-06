#' @rdname getModel
#'
#' @description
#' A method to get or set the INSPEcT_model object within an INSPEcT object. This method
#' is particularly useful to get and set testing parameters of the INSPEcT_model object
#' within the INSPEcT object.
#' @param object An object of class INSPEcT
#' @return An object of class INSPEcT model
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' getModel(nascentInspObj10)
setMethod('getModel', 'INSPEcT', function(object) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	return(object@model)
	})
#' @rdname getModel
#'
#' @param value An object of class INSPEcT model
setReplaceMethod('getModel', 'INSPEcT', function(object, value) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	object@model <- value
	object
	})
