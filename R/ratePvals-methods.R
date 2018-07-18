#' @rdname ratePvals
#'
#' @description
#' This method is used to retrieve all the p-values combined with Brown's method that combines the results of the 
#' log likelihood ratio test results for all pairs tested for each rate and all genes. P-values will change 
#' according to the threshold set for the chi-squared test because it influences the model that will be taken into 
#' consideration to perform log likelihood ratio tests. To have a sense of the best parameter to choose, a sythetic
#' data-set can be built and tested (\code{\link{makeSimModel}}, \code{\link{makeSimDataset}})
#' In case 'aic' has been selected via \code{\link{modelSelection}} method, 
#' this method assigns the chi-squared test result of the model selected by AIC
#' to the respective variable rates
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param cTsh A numeric representing the threshold for the chi-squared test to consider a model as valid
#' @details ratePvlas retrieve a single p-value for each rate thanks to multiple log likelihood tests performed on 
#' nested models that has a chi-squared test below the selected threshold. 
#' Among the many p-values that log likelihood ratio test calculate, a single p-value is obtaied applying Brown's method 
#' for combining dependent p-values.
#' @return A matrix containing p-values calculated for each rate
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}
#' @examples
#' data('nascentInspObj10', package='INSPEcT')
#' ratePvals(nascentInspObj10)
#' # calculate agin the p-values with Brown with a different threshold 
#' # for considering a model valid for the log likelihood ratio test
#' ratePvals(nascentInspObj10, cTsh=.2)
#' # Set permaenently the chi-squared threshold at .2 for nascentInspObj10 object
#' thresholds(nascentInspObj10)$chisquare <- .2
setMethod('ratePvals', 'INSPEcT_model', function(object, cTsh=NULL) {
	## calculates the pval to be varying per rate per gene, 
	## according to the threshold set for the chisq masking step)

	if( object@params$modelSelection=='llr' ) {
		if( is.null(cTsh) )
			cTsh <- object@params$thresholds$chisquare
		## retrieve models
		ratesSpecs <- object@ratesSpecs
		## define which models have been tested
		availableTests <- names(ratesSpecs[[1]])
		## generic tests
		chisq_pvals <- chisqtest(object)
		logLik_vals <- logLik(object)
		#########
		# alpha #
		#########
		alphaCols <- object@params$llrtests$synthesis
		alphaLLRtestPvlas <- sapply(alphaCols, function(compare) {
			null <- compare[1]; alt <- compare[2]
			sapply(ratesSpecs, function(x) 
				tryCatch(.logLikRatioTest(x[[null]], x[[alt]]), error=function(e) NA)
				)
			})
		alphaChisqMask <- sapply(alphaCols, function(compare) {
			chisq_pvals[,compare[1]] <= cTsh | 
				chisq_pvals[,compare[2]] <= cTsh})

		if( !is.matrix(alphaLLRtestPvlas) ) {
			alphaLLRtestPvlas <- t(as.matrix(alphaLLRtestPvlas))
			alphaChisqMask <- t(as.matrix(alphaChisqMask))
			rownames(alphaLLRtestPvlas) <- rownames(alphaChisqMask) <-
				names(ratesSpecs)
		}
		colnames(alphaLLRtestPvlas) <- colnames(alphaChisqMask) <-
			sapply(alphaCols, paste, collapse='_VS_')
		########
		# beta #
		########
		betaCols <- object@params$llrtests$degradation
		## make the logLikRatio tests
		betaLLRtestPvlas <- sapply(betaCols, function(compare) {
			null <- compare[1]; alt <- compare[2]
			sapply(ratesSpecs, function(x) 
				tryCatch(.logLikRatioTest(x[[null]], x[[alt]]), error=function(e) NA)
				)
			})
		betaChisqMask <- sapply(betaCols, function(compare) {
			chisq_pvals[,compare[1]] <= cTsh | 
				chisq_pvals[,compare[2]] <= cTsh})
		if( !is.matrix(betaLLRtestPvlas) ) {
			betaLLRtestPvlas <- t(as.matrix(betaLLRtestPvlas))
			betaChisqMask <- t(as.matrix(betaChisqMask))
			rownames(betaLLRtestPvlas) <- rownames(betaChisqMask) <-
				names(ratesSpecs)
		}
		colnames(betaLLRtestPvlas) <- colnames(betaChisqMask) <-
			sapply(betaCols, paste, collapse='_VS_')
		#########
		# gamma #
		#########
		gammaCols <- object@params$llrtests$processing
		## make the logLikRatio tests
		gammaLLRtestPvlas <- sapply(gammaCols, function(compare) {
			null <- compare[1]; alt <- compare[2]
			sapply(ratesSpecs, function(x) 
				tryCatch(.logLikRatioTest(x[[null]], x[[alt]]), error=function(e) NA)
				)
			})
		gammaChisqMask <- sapply(gammaCols, function(compare) {
			chisq_pvals[,compare[1]] <= cTsh | 
				chisq_pvals[,compare[2]] <= cTsh})
		if( !is.matrix(gammaLLRtestPvlas) ) {
			gammaLLRtestPvlas <- t(as.matrix(gammaLLRtestPvlas))
			gammaChisqMask <- t(as.matrix(gammaChisqMask))
			rownames(gammaLLRtestPvlas) <- rownames(gammaChisqMask) <-
				names(ratesSpecs)
		}
		colnames(gammaLLRtestPvlas) <- colnames(gammaChisqMask) <-
			sapply(gammaCols, paste, collapse='_VS_')

		# remove comparisons which have null variance due to the impossibility to evaluate a correlation mandatory for the brown's test

		alphaChisqMask <- alphaChisqMask[,which(apply(alphaLLRtestPvlas,2,var,na.rm=T)!=0)]
		betaChisqMask <- betaChisqMask[,which(apply(betaLLRtestPvlas,2,var,na.rm=T)!=0)]
		gammaChisqMask <- gammaChisqMask[,which(apply(gammaLLRtestPvlas,2,var,na.rm=T)!=0)]

		alphaLLRtestPvlas <- alphaLLRtestPvlas[,which(apply(alphaLLRtestPvlas,2,var,na.rm=T)!=0)]
		betaLLRtestPvlas <- betaLLRtestPvlas[,which(apply(betaLLRtestPvlas,2,var,na.rm=T)!=0)]
		gammaLLRtestPvlas <- gammaLLRtestPvlas[,which(apply(gammaLLRtestPvlas,2,var,na.rm=T)!=0)]

		#If we have a matrix of pValues we proceed with the Brown's test, otherwise it is useless
		if(is.matrix(alphaLLRtestPvlas)){synthesisBP <- .brown_method_mask(alphaLLRtestPvlas, alphaChisqMask)}
		else{synthesisBP <- alphaLLRtestPvlas}
		if(is.matrix(betaLLRtestPvlas)){degradationBP <- .brown_method_mask(betaLLRtestPvlas, betaChisqMask)}
		else{degradationBP <- betaLLRtestPvlas}
		if(is.matrix(gammaLLRtestPvlas)){processingBP <- .brown_method_mask(gammaLLRtestPvlas, gammaChisqMask)}
		else{processingBP <- gammaLLRtestPvlas}

		if(object@params$padjG)
		{
			synthesisBP <- p.adjust(synthesisBP,method="BH",n=length(synthesisBP))
			degradationBP <- p.adjust(degradationBP,method="BH",n=length(degradationBP))
			processingBP <- p.adjust(processingBP,method="BH",n=length(processingBP))
		}
		ratePvals <- data.frame(
			synthesis=synthesisBP
			, degradation=degradationBP
			, processing=processingBP
			)		
	} else if( object@params$modelSelection=='aic' ) {
		## assign the chi-squared test result of the model selected by AIC
		## to the respective variable rates
		x2test <- chisqtest(object)
		aictest <- AIC(object)
		aictest[is.na(aictest)] <- Inf
		modelNames <- colnames(aictest)
		ix <- apply(aictest, 1, which.min)
		geneBestModel <- modelNames[ix]
		genePvals <- x2test[cbind(1:length(ix),ix)]
		synthesisPval <- rep(1, length(ix))
		ix <- grepl('a', geneBestModel)
		synthesisPval[ix] <- genePvals[ix]
		degradationPval <- rep(1, length(ix))
		ix <- grepl('b', geneBestModel)
		degradationPval[ix] <- genePvals[ix]
		processingPval <- rep(1, length(ix))
		ix <- grepl('c', geneBestModel)
		processingPval[ix] <- genePvals[ix]
		ratePvals <- data.frame(
			synthesis=synthesisPval
			, degradation=degradationPval
			, processing=processingPval
			)		
	}
	# return
	return(ratePvals)
	})

#' @rdname ratePvals
setMethod('ratePvals', 'INSPEcT', function(object, cTsh=NULL) {
	return(ratePvals(object@model, cTsh))
	})

