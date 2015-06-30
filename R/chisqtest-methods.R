c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA)

#' @rdname chisqtest
#'
#' @description
#' This method is used to retrieve all the chi-squared test results for all models tested for all genes.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param ... Additional arguments for the generic
#' @return A matrix of chi-squared test results for all the tested models
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' chisqtest(mycerIds10)
setMethod('chisqtest', 'INSPEcT_model', function(object, ...) {
	exp(t(sapply(object@ratesSpecs, function(x) 
		tryCatch(sapply(x, '[[', 'test'), 
			error=function(e) c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA))
		)))
	})
#' @rdname chisqtest
setMethod('chisqtest', 'INSPEcT', function(object, ...) {
	chisqtest(object@model, ...)
	})

#' @rdname chisqmodel
#'
#' @description
#' This method is used to retrieve the chi-squared test results for the models
#' that have been selected to better represent the behavior of each gene.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param ... Additional arguments for the generic
#' @return A vector of chi-squared test results
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' chisqmodel(mycerIds10)
setMethod('chisqmodel', 'INSPEcT_model', function(object, ...) {
	# chisqtest <- exp(t(sapply(object@ratesSpecs, function(x) sapply(x, '[[', 'test'))))
	chisqtest <- chisqtest(object)
	gc <- geneClass(object)
	colid <- sapply(gc, function(x) which(colnames(chisqtest)==x))
	chisqtest[cbind(1:nrow(chisqtest),colid)]
	})
#' @rdname chisqmodel
setMethod('chisqmodel', 'INSPEcT', function(object, ...) {
	chisqmodel(object@model)
	})

#' @name AIC-INSPEcT-method
#' @title Akaike information criterion calculated for the models evaluated by INSPEcT
#' @description
#' This method is used to retrieve AIC values for all models tested for all genes.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param ... Additional arguments for the generic
#' @param k Additional parameter for the generic
#' @return A matrix of AIC values
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' AIC(mycerIds10)
NULL

#' @rdname AIC-INSPEcT-method
setMethod('AIC', 'INSPEcT_model', function(object, ..., k=2) {
	t(sapply(object@ratesSpecs, function(x) 
		tryCatch(sapply(x, '[[', 'AIC'),
			error=function(e) c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA))
		))
	})
#' @rdname AIC-INSPEcT-method
setMethod('AIC', 'INSPEcT', function(object, ..., k=2) {
	AIC(object@model, ..., k)
	})
