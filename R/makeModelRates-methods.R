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
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' tpts <- mycerIds10@tpts
#' eSet <- makeModelRates(getModel(mycerIds10), tpts=tpts)
#' exprs(eSet)
setMethod('makeModelRates', 'INSPEcT_model', function(object, ...) {
	args <- list(...)
	tpts <- args$tpts
	if( is.null(tpts) )
		stop('makeModelRates: missing "tpts" argument with no default.')
	## get ratesSpec field
	ratesSpecs <- object@ratesSpecs
	## in case some elements of ratesSpecs are longer than one,
	# meaning that a unique choiche for a model has not been done yet,
	# choose one using "bestModel" method
	if( any(sapply(ratesSpecs, length)!=1) )
		ratesSpecs <- .bestModel(object)@ratesSpecs
	## solve the differential equation model for each gene
	nGenes <- length(ratesSpecs)
	log_shift <- .find_tt_par(tpts)
	modelRates <- lapply(1:nGenes, function(i) {
		tryCatch(
			.makeModel(tpts, ratesSpecs[[i]][[1]], log_shift, 
				.time_transf, deSolve::ode, .rxnrate)
			, error=function(e)
				tryCatch(
					.makeSimpleModel(tpts, ratesSpecs[[i]][[1]], log_shift, 
						.time_transf, deSolve::ode, .rxnrateSimple)
					, error=function(e) .makeEmptyModel(tpts)
					)
				)
		})
	## make an objec of ExpressionSet class
	exprData <- cbind(
		t(sapply(modelRates, function(x) x$total))
		, t(sapply(modelRates, function(x) x$preMRNA))
		, t(sapply(modelRates, function(x) x$alpha))
		, t(sapply(modelRates, function(x) x$beta))
		, t(sapply(modelRates, function(x) x$gamma))
		)
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
	rownames(exprData) <- names(object@ratesSpecs)
	rownames(pData) <- colnames(exprData)
	phenoData <- new('AnnotatedDataFrame', data=pData)
	modelRates <- ExpressionSet(
		assayData=exprData
		, phenoData=phenoData
		)
	## return the ExpressionSet object
	return(modelRates)
	})

# setGeneric('makeModelRates', function(object, ...) standardGeneric('makeModelRates'))
#' @rdname makeModelRates
#' @description
#' This method can be used to regenerate the rates assiciated to the modeling, in case
#' some testing parameters has changed.
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' viewModelRates(mycerIds10, 'degradation')
#' ## force every degradation rate to be accepted as variable
#' thresholds(getModel(mycerIds10))$brown <- c(synthesis=.01, degradation=1, processing=.01)
#' mycerIds10 <- makeModelRates(mycerIds10)
#' viewModelRates(mycerIds10, 'degradation')

setMethod(f='makeModelRates', 'INSPEcT', definition=function(object, ...) {
	## get ratesSpec field
	ratesSpecs <- object@model@ratesSpecs
	tpts <- object@tpts
	log_shift <- .find_tt_par(tpts)
	
	## in case some elements of ratesSpecs are longer than one,
	# meaning that a unique choiche for a model has not been done yet,
	# choose one using "bestModel" method
	if( any(sapply(ratesSpecs, length)!=1) )
		ratesSpecs <- .bestModel(object@model)@ratesSpecs
		
	bestModels <- sapply(ratesSpecs,names)

	nGenes <- length(ratesSpecs)

	if(object@params$NoNascent & object@params$estimateRatesWith == "int")
	{
		a <- log_shift
		c <- abs(min(.time_transf(tpts,a)))

		## solve the differential equation model for each gene
		modelRates <- lapply(1:nGenes, function(i) {
			tryCatch(
				.makeModel(tpts = tpts
						 , hyp = ratesSpecs[[i]][[1]]
						 , log_shift = log_shift
						 , .time_transf = .time_transf_NoNascent
						 , ode = deSolve::ode
						 , .rxnrate = .rxnrate
						 , c = c)
				, error=function(e) .makeEmptyModel(tpts)
			)
		})

	}else if(!object@params$NoNascent)
	{
		## solve the differential equation model for each gene
		modelRates <- lapply(1:nGenes, function(i) {
			tryCatch(
				.makeModel(tpts, ratesSpecs[[i]][[1]], log_shift, 
					.time_transf, deSolve::ode, .rxnrate)
				, error=function(e)
					tryCatch(
						.makeSimpleModel(tpts, ratesSpecs[[i]][[1]], log_shift, 
							.time_transf, deSolve::ode, .rxnrateSimple)
						, error=function(e) .makeEmptyModel(tpts)
						)
					)
		})
	}else{

		a <- log_shift
		c <- abs(min(.time_transf(tpts,a)))
		## solve the differential equation model for each gene
		modelRates <- lapply(1:nGenes, function(i) {
			if(any(bestModels[i]==c("b","c","bc")))
			{
				tryCatch(
					.makeModel(tpts = tpts
							 , hyp = ratesSpecs[[i]][[1]]
							 , log_shift = log_shift
							 , .time_transf = .time_transf_NoNascent
							 , ode = deSolve::ode
							 , .rxnrate = .rxnrate
							 , c = c)
					, error=function(e) .makeEmptyModel(tpts)
				)
			}else{
				tryCatch(
					.makeModel_Derivative(tpts = tpts
							 , hyp = ratesSpecs[[i]][[1]]
							 , log_shift = log_shift
							 , .time_transf = .time_transf_NoNascent
							 , ode = deSolve::ode
							 , .rxnrate = .rxnrate
							 , c = c
							 , geneBestModel = bestModels[i])
					, error=function(e) .makeEmptyModel(tpts)
				)
			}
		})
	}
	## make an objec of ExpressionSet class
	exprData <- cbind(
		t(sapply(modelRates, function(x) x$total))
		, t(sapply(modelRates, function(x) x$preMRNA))
		, t(sapply(modelRates, function(x) x$alpha))
		, t(sapply(modelRates, function(x) x$beta))
		, t(sapply(modelRates, function(x) x$gamma))
		)
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
