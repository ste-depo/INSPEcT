#' @rdname setConfidenceIntervals
#'
#' @description
#' This function is used to set the confidence intervals in the nascent RNA mode.
#' @param object An object of class INSPEcT_model
#' @param confidenceIntervals list of confidence intervals.
#' @return An object of class ExpressionSet containing the confidence intervals.

setMethod(f='setConfidenceIntervals', 'INSPEcT', definition=function(object, confidenceIntervals) {

	tpts <- tpts(object)

	confidenceIntervals <- lapply(confidenceIntervals,function(g){cbind(g[[1]][,c(1,3,4)],g[[2]][,c(1,3,4)],g[[3]][,c(1,3,4)])})

	## make an objec of ExpressionSet class
	exprData <- cbind(t(sapply(confidenceIntervals, function(x) x[,1]))
					, t(sapply(confidenceIntervals, function(x) x[,2]))
					, t(sapply(confidenceIntervals, function(x) x[,3]))
					, t(sapply(confidenceIntervals, function(x) x[,4]))
					, t(sapply(confidenceIntervals, function(x) x[,5]))
					, t(sapply(confidenceIntervals, function(x) x[,6]))
					, t(sapply(confidenceIntervals, function(x) x[,7]))
					, t(sapply(confidenceIntervals, function(x) x[,8]))
					, t(sapply(confidenceIntervals, function(x) x[,9])))

	nTpts <- length(tpts)
	pData <- data.frame(
		feature=c(rep('synthesis_left',nTpts)
				, rep('synthesis_right',nTpts)
				, rep('synthesis_constant',nTpts)
				, rep('processing_left',nTpts)
				, rep('processing_right',nTpts)
				, rep('processing_constant',nTpts)
				, rep('degradation_left',nTpts)
				, rep('degradation_right',nTpts)
				, rep('degradation_constant',nTpts)
				)
		, time=rep(tpts))

	colnames(exprData) <- paste(pData$feature,signif(pData$time,2), sep='_')
	rownames(exprData) <- rownames(exprs(object@ratesFirstGuess))
	rownames(pData) <- colnames(exprData)
	phenoData <- new('AnnotatedDataFrame', data=pData)
	confidenceIntervals <- ExpressionSet(
		assayData=exprData
		, phenoData=phenoData
		)
	## update and return the object
	object@confidenceIntervals <- confidenceIntervals
	return(object)
})