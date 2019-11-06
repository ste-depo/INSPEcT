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
#' ## force every degradation rate to be accepted as variable
#' modelSelection(getModel(nascentInspObj10))$thresholds$brown <- 
#'   c(synthesis=.01, degradation=1, processing=.01)
#' nascentInspObj10 <- makeModelRates(nascentInspObj10)
#' viewModelRates(nascentInspObj10, 'degradation')

setMethod(f='makeModelRates', 'INSPEcT', definition=function(object, ...) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	## get ratesSpec field
	ratesSpecs <- object@model@ratesSpecs

	# I split the genes in eiGenes and eGenes according to the message saved after the modeling
	allGenes <- names(ratesSpecs)
	eGenes <- allGenes[sapply(ratesSpecs,"[[","abc")["message",]=="eGene"]
	eiGenes <- setdiff(allGenes,eGenes)

	tpts <- object@tpts	

	## in case some elements of ratesSpecs are longer than one,
	# meaning that a unique choiche for a model has not been done yet,
	# choose one using "bestModel" method
	if( any(sapply(ratesSpecs, length)!=1) )
	{	
		if(object@NoNascent)
		{
			ratesSpecs <- .bestModel(object@model)@ratesSpecs			
		}else{
			ratesSpecs <- .bestModel(object@model,Nascent=TRUE)@ratesSpecs # Always VVV
		}
	}
	names(ratesSpecs) <- allGenes
	bestModels <- sapply(ratesSpecs,names)

	if(object@params$estimateRatesWith=="der")
	{
		eiModelRates <- lapply(eiGenes, function(i)
		{
			tryCatch(.makeModel_Derivative(tpts = tpts, hyp = ratesSpecs[[i]][[1]], geneBestModel = bestModels[i])
		   ,error=function(e) .makeEmptyModel(tpts))
		})
		names(eiModelRates) <- eiGenes

		if(length(eGenes)>0)
		{
			print("makeModelRates: simple still to be implemented!")
			# eModelRates <- lapply(eGenes, function(i)
			# {
			# 	tryCatch(.makeModel_Derivative_Simple(tpts = tpts, hyp = ratesSpecs[[i]][[1]], geneBestModel = bestModels[i])
			#    ,error=function(e) .makeEmptyModel(tpts))
			# })				
			# names(eModelRates) <- eGenes
			# eiModelRates <- c(eiModelRates,eModelRates)
		}
	}else{
		eiModelRates <- lapply(eiGenes, function(i)
		{
			tryCatch(.makeModel(tpts = tpts, hyp = ratesSpecs[[i]][[1]], nascent=FALSE)
		   ,error=function(e) .makeEmptyModel(tpts))
		})
		names(eiModelRates) <- eiGenes

		if(length(eGenes)>0)
		{
			print("makeModelRates: simple still to be implemented!")
			# eModelRates <- lapply(eGenes, function(i)
			# {
			# 	tryCatch(.makeModel(tpts = tpts, hyp = ratesSpecs[[i]][[1]], nascent = FALSE)
			#    ,error=function(e) .makeEmptyModel(tpts))
			# })
		
			# names(eModelRates) <- eGenes
			# eiModelRates <- c(eiModelRates,eModelRates)
		}
	}
	modelRates <- eiModelRates
	## make an objec of ExpressionSet class
	exprData <- cbind(t(sapply(modelRates, function(x) x$total))
					, t(sapply(modelRates, function(x) x$preMRNA))
					, t(sapply(modelRates, function(x) x$alpha))
					, t(sapply(modelRates, function(x) x$beta))
					, t(sapply(modelRates, function(x) x$gamma)))
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
	## update and return the object
	object@modelRates <- modelRates

	return(object)
})