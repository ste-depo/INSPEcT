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
	out <- chisqtest_internal(object)
	colnames(out) <- c("no-reg","s","p","d","sd","sp","pd","spd")
	return(out)
	})
#' @rdname chisqtest
setMethod('chisqtest', 'INSPEcT', function(object, ...) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	chisqtest(object@model, ...)
})

chisqtest_internal <- function(object) {
	exp(t(sapply(object@ratesSpecs, function(x) 
		tryCatch(sapply(x, '[[', 'test'), 
						 error=function(e) c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA))
	)))
}

#' @rdname convergence
#'
#' @description
#' This method is used to retrieve the convergence of the models
#' that have been selected to better represent the behavior of each gene.
#' 0 is converged, 1 not converged, 10 degenerated
#' @param object An object of class INSPEcT or INSPEcT_model
#' @return A vector of numeric
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' convergence(nascentInspObj10)
setMethod('convergence', 'INSPEcT', function(object) {
	out <- t(sapply(object@model@ratesSpecs, function(gene) sapply(gene, function(model) model$convergence)))
	colnames(out) <- c("no-reg","s","p","d","sd","sp","pd","spd")
	return(out)
})

#' @rdname chisqmodel
#'
#' @description
#' This method is used to retrieve the chi-squared test results for the models
#' that have been selected to better represent the behavior of each gene.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param gc Additional arguments for the generic
#' @param tpts Additional arguments for the generic
#' @param ... Additional arguments for the generic
#' @return A vector of chi-squared test results
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' chisqmodel(nascentInspObj10)
setMethod('chisqmodel', 'INSPEcT', function(object, ...) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	if(object@NoNascent)
	{
		chisqmodelfun(object@model)
	}else{
		chisqmodelfun(object@model, gc=geneClassInternal(object), tpts=tpts(object))
	}
	})

chisqmodelfun <- function(object, gc=NULL, tpts=NULL, ...) 
{
	
	if(is.null(gc)&is.null(tpts))
	{
		gc <- geneClassInternal(object)
		NoNascent <- TRUE
	}else{
		NoNascent <- FALSE
	}
	
	chisqtest <- chisqtest_internal(object)	
	
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
}


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
	out <- AIC_internal(object)
	colnames(out) <- c("no-reg","s","p","d","sd","sp","pd","spd")
	return(out)
	})
#' @rdname AIC-INSPEcT-method
setMethod('AIC', 'INSPEcT', function(object, ..., k=2) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	AIC(object@model, ..., k)
	})

AIC_internal <- function(object) {
	t(sapply(object@ratesSpecs, function(x) 
		tryCatch(sapply(x, '[[', 'AIC'),
						 error=function(e) c("0"=NA,"a"=NA,"b"=NA,"c"=NA,"ab"=NA,"bc"=NA,"ac"=NA,"abc"=NA))
	))
}
