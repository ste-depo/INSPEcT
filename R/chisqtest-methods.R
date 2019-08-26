#' @rdname chisqtest
#'
#' @description
#' This method is used to retrieve all the chi-squared test results for all models tested for all genes.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param ... Additional arguments for the generic
#' @return A matrix of chi-squared test results for all the tested models
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' chisqtest(nascentInspObj10)
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
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' chisqmodel(nascentInspObj10)
setMethod('chisqmodel', 'INSPEcT_model', function(object, ...) {
	# chisqtest <- exp(t(sapply(object@ratesSpecs, function(x) sapply(x, '[[', 'test'))))
	if(is.null(gc)&is.null(tpts))
	{
		gc <- geneClass(object)
		NoNascent <- TRUE
	}else{
		NoNascent <- FALSE
	}
	
	chisqtest <- chisqtest(object)	
	
	if(NoNascent)
	{
		colid <- sapply(gc, function(x) which(colnames(chisqtest)==x))
		chisqmodel <- chisqtest[cbind(1:nrow(chisqtest),colid)]
		names(chisqmodel) <- rownames(chisqtest)
		return(chisqmodel)
	}else{
		colid <- chisqtest[,"abc"]
		oldDF <- sapply(object@ratesSpecs,function(g)
		{
			dataTmp <- unlist(g[[8]])
			sum(unlist(dataTmp[grep(".df",names(dataTmp))]))
		})
		newDF <- ((oldDF/3)-1)*sapply(gc,function(c)length(strsplit(c,"")[[1]])) + 3

		oldDF <- sapply(oldDF,function(oldDF)max(0,3*length(tpts) - oldDF))
		newDF <- sapply(newDF,function(newDF)max(0,3*length(tpts) - newDF))

		pchisq(qchisq(colid,oldDF),newDF)
	}
	})
#' @rdname chisqmodel
setMethod('chisqmodel', 'INSPEcT', function(object, ...) {
	if(object@NoNascent)
	{
		chisqmodel(object@model)
	}else{
		chisqmodel(object@model, gc=geneClass(object), tpts=tpts(object))
	}
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
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' AIC(nascentInspObj10)
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
