#' @rdname makeModelRates
#'
#' @description
#' This function is used to evaluate rates and concentrations after modeling of the 
#' rates has been run with \code{\link{modelRates}}. The modeled rates are in functional 
#' form and can be evaluated at any time points.
#' @param object An object of class INSPEcT_model
#' @param ... additional arguments
#'  tpts : A vector of time points where rates and concentrations have to be evaluated
#' @return An object of class ExpressionSet containing the modeled rates and concentrations

setMethod('makeModelRates', 'INSPEcT_model', function(object, ...) {

	args <- list(...)
	tpts <- args$tpts

	if( is.null(tpts) )
		stop('makeModelRates: missing "tpts" argument with no default.')

	models <- object@ratesSpecs
	
	if( names(models[[1]][[1]])[1] == 'alpha') {
		realSynthesis <- t(sapply(models,function(model)model[[1]][['alpha']]$fun$value(tpts,model[[1]][['alpha']]$par)))
		realProcessing <- t(sapply(models,function(model)model[[1]][['gamma']]$fun$value(tpts,model[[1]][['gamma']]$par)))
		realDegradation <- t(sapply(models,function(model)model[[1]][['beta']]$fun$value(tpts,model[[1]][['beta']]$par)))
	
		rownames(realSynthesis) <- rownames(realProcessing) <- rownames(realDegradation) <- seq_along(models)
		colnames(realSynthesis) <- colnames(realProcessing) <- colnames(realDegradation) <- tpts
	
		return(list('synthesis'=realSynthesis,'processing'=realProcessing,'degradation'=realDegradation))
	} else {
		stop('makeModelRates: not working for derivative models.')
	}
})

#' @rdname makeModelRates
#' @description
#' This method can be used to regenerate the rates assiciated to the modeling, in case
#' some testing parameters has changed.
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' viewModelRates(nascentInspObj10, 'degradation')
#' ## force every degradation rate to be accepted as variable (makeModelRates is called internally)
#' nascentInspObj10 <- calculateRatePvals(nascentInspObj10, p_variability = c(.05,.05,1))
#' viewModelRates(nascentInspObj10, 'degradation')

setMethod(f='makeModelRates', 'INSPEcT', definition=function(object, ...) {
	checkINSPEcTObjectversion(object)
	
	## get ratesSpec field
	ratesSpecs <- object@model@ratesSpecs
	allGenes <- names(ratesSpecs)

	tpts <- object@tpts	

	# choose the best model for each gene in ratesSpecs 
	if( any(sapply(ratesSpecs, length)!=1) )
	{
		if(object@NoNascent)
		{
			bestModels <- geneClassInternal(object)
		}else{
			bestModels <- rep('abc', length(allGenes))
			names(bestModels) <- allGenes
		}
	}
	ratesSpecs <- lapply(seq_along(ratesSpecs), function(i) ratesSpecs[[i]][bestModels[i]])
	names(ratesSpecs) <- allGenes
	
	if(object@params$estimateRatesWith=="der")
	{
		eiModelRates <- lapply(allGenes, function(i)
		{
			tryCatch(.makeModel_Derivative(tpts = tpts, hyp = ratesSpecs[[i]][[1]], geneBestModel = bestModels[i])
		   ,error=function(e) .makeEmptyModel(tpts))
		})
		names(eiModelRates) <- allGenes
	}else{
		eiModelRates <- lapply(allGenes, function(i)
		{
			tryCatch(.makeModel(tpts = tpts, hyp = ratesSpecs[[i]][[1]], nascent=FALSE)
		   ,error=function(e) .makeEmptyModel(tpts))
		})
		names(eiModelRates) <- allGenes
	}
	modelRates <- eiModelRates
	## make an objec of ExpressionSet class
	total <- t(sapply(modelRates, function(x) x$total))
	preMRNA <- t(sapply(modelRates, function(x) x$preMRNA))
	alpha <- t(sapply(modelRates, function(x) x$alpha))
	beta <- t(sapply(modelRates, function(x) x$beta))
	gamma <- t(sapply(modelRates, function(x) x$gamma))
	## put eGenes preMRNA and gamma (k2) values to NAs
	firstGuess_preMRNA <- ratesFirstGuess(object, 'preMRNA')
	eGenes <- rownames(firstGuess_preMRNA[!is.finite(firstGuess_preMRNA[,1]),])
	if( !is.null(eGenes) ) {
	  preMRNA[eGenes,] <- NA
	  gamma[eGenes,] <- NA
	}
	## finalize ExpressionSet
	exprData <- cbind(total
					, preMRNA
					, alpha
					, beta
					, gamma)
	nTpts <- length(tpts)
	pData <- data.frame(
		feature=c(
			rep('total',nTpts)
			, rep('preMRNA',nTpts)
			, rep('synthesis',nTpts)
			, rep('degradation',nTpts)
			, rep('processing',nTpts)
			)
		, time=rep(tpts, 5))
	colnames(exprData) <- paste(pData$feature, 
		signif(pData$time,2), sep='_')
	rownames(exprData) <- rownames(exprs(object@ratesFirstGuess))
	rownames(pData) <- colnames(exprData)
	phenoData <- new('AnnotatedDataFrame', data=pData)
	modelRates <- ExpressionSet(
		assayData=exprData
		, phenoData=phenoData
		)
	## update the object 
	object@modelRates <- modelRates
	return(object)
})