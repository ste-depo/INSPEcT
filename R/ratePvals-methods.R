#' @rdname ratePvals
#'
#' @description
#' This method is used to retrieve all the p-values relative to the variability of synthesis, processing and degradation rates.
#' @param object An object of class INSPEcT
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' ratePvals(nascentInspObj10)
setMethod('ratePvals', 'INSPEcT', function(object) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	if( !is.numeric(tpts(object)) ) {
		stop("Run 'compareSteady' method on this object to evaluate differential rate ")
	}
	object@ratePvals
})

#' @rdname calculateRatePvals
#' @description
#' This method is used to calculate all the p-values relative to the variability of synthesis, processing and degradation rates.
#' For object modeled with nascent RNA or when non-functional modeling was used, the variability is calculated using the 
#' confidence intervals. For objects modeled without nascent RNA, model selection is performed by comparing the likelihood of 
#' different (nested) models.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param modelSelection 'aic' compares nested models closest to the one with lowest AIC, 'llr' compares all nested models, 
#' 'hib' is a mix between the previous two. (default 'aic')
#' @param preferPValue a logical, if TRUE (default) limit the search for best models among the ones with succeded the goodness of fit test.
#' @param padj a logical, if TRUE (default) correct the p-values for multiple testing
#' @param p_goodness_of_fit a numeric, the threshold for the goodness-of-fit test (default = .1)
#' @param p_variability a numeric, a vector with the thresholds for the p-value of the variability test (one threshold for each rate, default = rep(.05, 3))
#' @param limitModelComplexity a logical that limits the complexity of the function used to describe dynamics to the length of the time-course (default = FALSE)
#' @details ratePvals retrieve a single p-value for each rate and gene associated to its variability (null hypothesis = the rate is not
#' changing between the conditions)
#' @return A matrix containing p-values calculated for each rate
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' # Set the chi-squared threshold at .2 for nascentInspObj10 object
#' nascentInspObj10 <- calculateRatePvals(nascentInspObj10, p_goodness_of_fit=.2)
setMethod('calculateRatePvals', 'INSPEcT', function(object
																										, modelSelection = c('aic','llr','hib')
																										, preferPValue = TRUE
																										, padj = TRUE
																										, p_goodness_of_fit = .1
																										, p_variability=rep(.05, 3)
																										, limitModelComplexity = FALSE) 
	{
	
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	# if( nrow(object@modelRates) == 0 )
	# 	stop('ratePvals: run modelRates or modelRatesNF before.')
	modelSelection <- modelSelection[1]
	if( !modelSelection %in% c('aic','llr','hib') )
		stop('calculateRatePvals: modelSelection argument must be either "aic", "llr" or "hib"')
	if( !is.logical(preferPValue) )
		stop('calculateRatePvals: preferPValue argument must be a logical')		
	if( !is.logical(padj) )
		stop('calculateRatePvals: padj argument must be a logical')
	if( !is.numeric(p_goodness_of_fit) )
		stop('calculateRatePvals: p_goodness_of_fit argument must be a numeric')
	if( !is.numeric(p_variability) )
		stop('calculateRatePvals: p_variability argument must be a numeric')
	if( length(p_variability) != 3 )
		stop('calculateRatePvals: p_variability argument must have length 3')
	if( !is.logical(limitModelComplexity) )
		stop('calculateRatePvals: limitModelComplexity argument must be a logical')		
	message('Calculating rate p-values...')

	## store parameters used for model selection
	object@model@params <- list(
		modelSelection =modelSelection
		, preferPValue = preferPValue
		, padj = padj
		, p_goodness_of_fit = p_goodness_of_fit
		, p_variability = p_variability
		, limitModelComplexity = limitModelComplexity)
	
	if( !object@NoNascent | object@NF )
	# in case of Nascent mode or Non-Functional (both Nascent and No-Nascent)
	# perform the estimation of variability using the confidence intervals
	{

		synthesis_left <- viewConfidenceIntervals(object,"synthesis_left")
		synthesis_center <- viewModelRates(object,"synthesis")
		synthesis_right <- viewConfidenceIntervals(object,"synthesis_right")

		processing_left <- viewConfidenceIntervals(object,"processing_left")
		processing_center <- viewModelRates(object,"processing")
		processing_right <- viewConfidenceIntervals(object,"processing_right")

		degradation_left <- viewConfidenceIntervals(object,"degradation_left")
		degradation_center <- viewModelRates(object,"degradation")
		degradation_right <- viewConfidenceIntervals(object,"degradation_right")
		
		fitResults_synthesis <- unlist(lapply(featureNames(object),function(g)
		{
			rate_conf_int <- cbind(synthesis_left[g,],synthesis_center[g,],synthesis_right[g,])
			k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
			if(!is.finite(k_start)) return(NaN) #return(list(par=NaN, value=NaN))
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int)
			pchisq(k_scores_out$value,length(tpts(object))-1,lower.tail=FALSE)
		}))

		fitResults_processing <- unlist(lapply(featureNames(object),function(g)
		{
			rate_conf_int <- cbind(processing_left[g,],processing_center[g,],processing_right[g,])
			k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
			if(!is.finite(k_start)) return(NaN) #return(list(par=NaN, value=NaN))
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int)
			pchisq(k_scores_out$value,length(tpts(object))-1,lower.tail=FALSE)
		}))

		fitResults_degradation <- unlist(lapply(featureNames(object),function(g)
		{
			rate_conf_int <- cbind(degradation_left[g,],degradation_center[g,],degradation_right[g,])
			k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
			if(!is.finite(k_start)) return(NaN) #return(list(par=NaN, value=NaN))
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int)
			pchisq(k_scores_out$value,length(tpts(object))-1,lower.tail=FALSE)
		}))
		
		if(modelSelection(object)$padj)
		{
			rate_pvals <- data.frame(
				"synthesis"=p.adjust(fitResults_synthesis, method="BH", n=length(fitResults_synthesis)),
				"processing"=p.adjust(fitResults_processing, method="BH", n=length(fitResults_processing)),
				"degradation"=p.adjust(fitResults_degradation, method="BH", n=length(fitResults_degradation)),
				row.names=featureNames(object))
		} else {
			rate_pvals <- data.frame(
				"synthesis"=fitResults_synthesis,
				"processing"=fitResults_processing,
				"degradation"=fitResults_degradation,
				row.names=featureNames(object))
		}
	} else { # NoNascent
		if( modelSelection(object)$limitModelComplexity ) {
			dfmax = length(tpts(object))
			rate_pvals = calculate_rates_pvalues(object@model, bTsh=p_variability, cTsh=p_goodness_of_fit, dfmax)
		} else {
			rate_pvals = calculate_rates_pvalues(object@model, bTsh=p_variability, cTsh=p_goodness_of_fit)
		}
	}
	object@ratePvals <- rate_pvals
	# in case of NoNascent update the modeled rates after the update of ratePvals
	if( object@NoNascent ) object <- makeModelRates(object)
	return(object)
	})

calculate_rates_pvalues <- function(object, bTsh, cTsh, dfmax=Inf) {
	## calculates the pval to be varying per rate per gene, 
	## according to the threshold set for the chisq masking step)
	# priors <- object@params$priors
	# if(is.null(priors)) 
		priors<-c("synthesis"=1,"processing"=1,"degradation"=1)
	
	## in case a rate has a threshold set to 0, avoid testing it
	rates_to_avoid <- names(bTsh)[bTsh == 0]
	rates_to_avoid <- c('synthesis'='a','degradation'='b','processing'='c')[rates_to_avoid]
	llrtests=list(
		synthesis=list(c('0','a')
			,c('b','ab')
			,c('c','ac')
			,c('bc','abc'))
		, degradation=list(c('0','b')
			,c('a','ab')
			,c('c','bc')
			,c('ac','abc'))
		, processing=list(c('0','c')
			,c('a','ac')
			,c('b','bc')
			,c('ab','abc'))
		)

	## retrieve models
	ratesSpecs <- object@ratesSpecs

	llr_temp_function <- function()
	{
			
		## generic tests
		chisq_pvals <- chisqtest_internal(object)
		logLik_vals <- logLik_internal(object)

		### New lines "copied" from AIC method - control repeated in logLikRatioTestInscpectModels
		chisq_pvals[!is.finite(chisq_pvals)] <- 1
		logLik_vals[!is.finite(logLik_vals)] <- -Inf

		#########
		# alpha #
		#########
		alphaCols <- llrtests$synthesis
		if( length(rates_to_avoid)>0 )
			alphaCols <- alphaCols[!sapply(alphaCols, function(x) any(grepl(rates_to_avoid, x)))]
		if( length(alphaCols)>0 ) {
			alphaLLRtestPvlas <- sapply(alphaCols, function(compare) {
				null <- compare[1]; alt <- compare[2]
				sapply(ratesSpecs, function(x) 
					tryCatch(logLikRatioTestInscpectModels(x[[null]], x[[alt]], dfmax = dfmax, constantProbability=priors["synthesis"]), error=function(e) NA)
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
			# alphaChisqMask <- alphaChisqMask[,which(apply(alphaLLRtestPvlas,2,var,na.rm=T)!=0)]
			# alphaLLRtestPvlas <- alphaLLRtestPvlas[,which(apply(alphaLLRtestPvlas,2,var,na.rm=T)!=0)]
		} else {
			alphaLLRtestPvlas <- rep(1, length(ratesSpecs))
		}
		########
		# beta #
		########
		betaCols <- llrtests$degradation
		if( length(rates_to_avoid)>0 )
			betaCols <- betaCols[!sapply(betaCols, function(x) any(grepl(rates_to_avoid, x)))]
		## make the logLikRatio tests
		if( length(betaCols)>0 ) {
			betaLLRtestPvlas <- sapply(betaCols, function(compare) {
				null <- compare[1]; alt <- compare[2]
				sapply(ratesSpecs, function(x) 
					tryCatch(logLikRatioTestInscpectModels(x[[null]], x[[alt]], dfmax = dfmax, constantProbability=priors["degradation"]), error=function(e) NA)
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

			# betaChisqMask <- betaChisqMask[,which(apply(betaLLRtestPvlas,2,var,na.rm=T)!=0)]
			# betaLLRtestPvlas <- betaLLRtestPvlas[,which(apply(betaLLRtestPvlas,2,var,na.rm=T)!=0)]
		} else {
			betaLLRtestPvlas <- rep(1, length(ratesSpecs))
		}
		#########
		# gamma #
		#########
		gammaCols <- llrtests$processing
		if( length(rates_to_avoid)>0 )
			gammaCols <- gammaCols[!sapply(gammaCols, function(x) any(grepl(rates_to_avoid, x)))]
		## make the logLikRatio tests
		if( length(gammaCols)>0 ) {
			gammaLLRtestPvlas <- sapply(gammaCols, function(compare) {
				null <- compare[1]; alt <- compare[2]
				sapply(ratesSpecs, function(x) 
					tryCatch(logLikRatioTestInscpectModels(x[[null]], x[[alt]], dfmax = dfmax, constantProbability=priors["processing"]), error=function(e) NA)
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
			# gammaChisqMask <- gammaChisqMask[,which(apply(gammaLLRtestPvlas,2,var,na.rm=T)!=0)]
			# gammaLLRtestPvlas <- gammaLLRtestPvlas[,which(apply(gammaLLRtestPvlas,2,var,na.rm=T)!=0)]
		} else {
			gammaLLRtestPvlas <- rep(1, length(ratesSpecs))
		} 
		#If we have a matrix of pValues we proceed with the Brown's test, otherwise it is useless
		if(is.matrix(alphaLLRtestPvlas)){synthesisBP <- brown_method_mask(alphaLLRtestPvlas, alphaChisqMask)}else{synthesisBP <- alphaLLRtestPvlas}
		if(is.matrix(betaLLRtestPvlas)){degradationBP <- brown_method_mask(betaLLRtestPvlas, betaChisqMask)}else{degradationBP <- betaLLRtestPvlas}
		if(is.matrix(gammaLLRtestPvlas)){processingBP <- brown_method_mask(gammaLLRtestPvlas, gammaChisqMask)}else{processingBP <- gammaLLRtestPvlas}
		
		if(modelSelection(object)$padj)
		{
			synthesisBP <- p.adjust(synthesisBP,method="BH",n=length(synthesisBP))
			degradationBP <- p.adjust(degradationBP,method="BH",n=length(degradationBP))
			processingBP <- p.adjust(processingBP,method="BH",n=length(processingBP))
		}

		ratePvals <- data.frame(
			synthesis=synthesisBP
			, processing=processingBP
			, degradation=degradationBP
			)

		ratePvals <- replace(ratePvals,is.na(ratePvals),1)
		### New controls to uniform ratePvals and geneClass
		if(object@params$preferPValue)
		{
			for(i in rownames(chisq_pvals))
			{
				acceptableModelsTemp <- which(chisq_pvals[i,] <= cTsh)
				if(length(acceptableModelsTemp)==0)
				{
					if(is.finite(sort(chisq_pvals[i,])[1]))
					{
						if(grepLogic("a",names(sort(chisq_pvals[i,])[1]))){ratePvals[i,"synthesis"] <- 0}else{ratePvals[i,"synthesis"] <- 1}
						if(grepLogic("c",names(sort(chisq_pvals[i,])[1]))){ratePvals[i,"processing"] <- 0}else{ratePvals[i,"processing"] <- 1}
						if(grepLogic("b",names(sort(chisq_pvals[i,])[1]))){ratePvals[i,"degradation"] <- 0}else{ratePvals[i,"degradation"] <- 1}
					}else{
						ratePvals[i,"synthesis"] <- 1
						ratePvals[i,"processing"] <- 1
						ratePvals[i,"degradation"] <- 1
					}
				}
			
				if(length(acceptableModelsTemp)==1)
				{
					if(grepLogic("a",colnames(chisq_pvals)[acceptableModelsTemp])){ratePvals[i,"synthesis"] <- 0}else{ratePvals[i,"synthesis"] <- 1}
					if(grepLogic("c",colnames(chisq_pvals)[acceptableModelsTemp])){ratePvals[i,"processing"] <- 0}else{ratePvals[i,"processing"] <- 1}
					if(grepLogic("b",colnames(chisq_pvals)[acceptableModelsTemp])){ratePvals[i,"degradation"] <- 0}else{ratePvals[i,"degradation"] <- 1}
				}
			}
		}
		ratePvals
	}
	aic_temp_function <- function()
	{
		## assign the chi-squared test result of the model selected by AIC
		## to the respective variable rates
		aictest <- AIC_internal(object)
		aictest[is.na(aictest)] <- Inf
		if( length(rates_to_avoid)>0 )
			aictest = aictest[,grep(rates_to_avoid, colnames(aictest), invert=TRUE)]
		gene_class <- colnames(aictest)[apply(aictest, 1, which.min)]
		ratePvals <- data.frame(t(sapply(seq_along(ratesSpecs), function(i) {
			switch(gene_class[i],
				'0' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['0']] ,ratesSpecs[[i]][['a']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['0']] ,ratesSpecs[[i]][['c']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['0']] ,ratesSpecs[[i]][['b']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'a' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['0']] ,ratesSpecs[[i]][['a']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['a']] ,ratesSpecs[[i]][['ac']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['a']] ,ratesSpecs[[i]][['ab']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'c' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['c']] ,ratesSpecs[[i]][['ac']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['0']] ,ratesSpecs[[i]][['c']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['c']] ,ratesSpecs[[i]][['bc']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'b' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['b']] ,ratesSpecs[[i]][['ab']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['b']] ,ratesSpecs[[i]][['bc']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['0']] ,ratesSpecs[[i]][['b']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'ac' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['c']] ,ratesSpecs[[i]][['ac']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['a']] ,ratesSpecs[[i]][['ac']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['ac']],ratesSpecs[[i]][['abc']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'ab' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['b']] ,ratesSpecs[[i]][['ab']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['ab']],ratesSpecs[[i]][['abc']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['a']] ,ratesSpecs[[i]][['ab']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'bc' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['bc']],ratesSpecs[[i]][['abc']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['b']], ratesSpecs[[i]][['bc']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['c']], ratesSpecs[[i]][['bc']], dfmax = dfmax, constantProbability=priors["degradation"])
					),
				'abc' = c(
					synthesis=logLikRatioTestInscpectModels(ratesSpecs[[i]][['bc']],ratesSpecs[[i]][['abc']], dfmax = dfmax, constantProbability=priors["synthesis"]),
					processing=logLikRatioTestInscpectModels(ratesSpecs[[i]][['ab']],ratesSpecs[[i]][['abc']], dfmax = dfmax, constantProbability=priors["processing"]),
					degradation=logLikRatioTestInscpectModels(ratesSpecs[[i]][['ac']],ratesSpecs[[i]][['abc']], dfmax = dfmax, constantProbability=priors["degradation"])
					)
				)
			})))
		rownames(ratePvals) <- rownames(aictest)
		if( length(rates_to_avoid)>0 ) {
			ratePvals[,c('a'='synthesis','b'='degradation','c'='processing')[rates_to_avoid]] <- 1
		}
 		if(modelSelection(object)$padj) {
			ratePvals <- data.frame(apply(ratePvals, 2, p.adjust, method="BH", n=nrow(ratePvals)))
			if(ncol(ratePvals)!=3){ratePvals <- t(ratePvals); rownames(ratePvals) <- rownames(aictest)}
		}
		return(ratePvals)
	}
	if( object@params$modelSelection=='llr' ) {
		ratePvals <- llr_temp_function()
	} else if( object@params$modelSelection=='aic' ) {
		ratePvals <- aic_temp_function()
	} else if( object@params$modelSelection=='hib' ) {
		
		llrMatrix <- llr_temp_function()
		aicMatrix <- aic_temp_function()
		ratePvals <- cbind(synthesis=aicMatrix[,"synthesis"],processing=llrMatrix[,"processing"],degradation=llrMatrix[,"degradation"])
		rownames(ratePvals) <- rownames(llrMatrix)
		ratePvals <- as.data.frame(ratePvals)
	}
	# return
	return(ratePvals)

}

brown_method <- function(y, na.rm=FALSE) {
	## initialize answer as a vector of NA
	ans <- rep(NA, nrow(y))
	names(ans) <- rownames(y)
	## count per each row (gene) the number of tests that
	## could actually being performed
	k <- apply(y,1,function(x) length(na.omit(x)))
	## if it's zero for every gene return the NA vector
	if( all(k==0) ) return(ans)
	## otherwise calculate the brown's test for the available genes
	## dividing them in cathegories, according to the availbale tests
	y_logical <- matrix(FALSE, nrow(y), ncol(y))
	y_logical[!is.na(y)] <- TRUE
	groups <- unique(as.data.frame(y_logical))
	if( nrow(groups)==1 )
		belongs_to <- rep(1, nrow(y_logical))
	else
		belongs_to <- apply(apply(y_logical,1,function(x) 
			apply(groups, 1, function(y) all(x==y))),2,which)
	fisherSum <- function(x) 
		apply(x,1,function(z) -2*sum(log(z), na.rm=na.rm))
	tmp <- c(-2.59,-2.382,-2.17,-1.946,-1.709,-1.458,-1.194,
		-.916,-.625,-.320,0,.334,.681,1.044,1.421,1.812,
		2.219,2.641,3.079,3.531,4)
	for( i in 1:nrow(groups) ) {
		acceptedTests <- as.logical(groups[i,])
		if( any(acceptedTests) ) {
			ix <- belongs_to==i
			y_group <- y[,acceptedTests, drop=FALSE]
			fX2 <- fisherSum(y_group[ix, , drop=FALSE])
			k_group <- length(which(acceptedTests))
			xout_matrix <- suppressWarnings(stats::cor(y_group, use='complete.obs'))
			xout_vector <- xout_matrix[which(as.vector(lower.tri(xout_matrix)))]
			s2X2 <- 4 * k_group + 2 * sum(approx(seq(-1,1,.1),tmp,xout=xout_vector)$y)
			f <- 2 * (2 * k_group)^2 /s2X2
			correction <- s2X2/(2*2*k_group)
			ans[ix] <- pchisq(fX2/correction, f, lower.tail=FALSE)
		}
	}
	return(ans)
}

brown_method_mask <- function(y, mask) {
	# if there is only one gene, no assumption about
	# the correlation of the tests (unless it is set and fixed)
	# can be done, therefore fisher is done
	if( nrow(y)==1 ) {
		y_filtered <- y[mask]
		fX2 <- -2*sum(log(y_filtered))
		f <- 2*length(y_filtered)
		out <- pchisq(fX2, f, lower.tail=FALSE)
		names(out) <- rownames(y)
	} else {
		y[c(!mask)] <- NA
		out <- brown_method(y, na.rm=TRUE)
	}
	return(out)
}

logLikRatioTestInscpectModels <- function(null, alt, dfmax=Inf, constantProbability=1)
{
	if(is.finite(null$logLik)){nullLL <- null$logLik}else{nullLL <- (-Inf)}
	if(is.finite(alt$logLik)){altnullLL <- alt$logLik}else{altnullLL <- (-Inf)}
	
	D <- - 2*nullLL + 2*altnullLL

	df_null <- min(dfmax,sum(sapply(null[names(null)%in%c("alpha","beta","gamma","total","mature","premature")],"[[","df")))
	df_alt <- min(dfmax,sum(sapply(alt[names(alt)%in%c("alpha","beta","gamma","total","mature","premature")],"[[","df")))

	if(constantProbability == 1){return(pchisq(D, df_alt-df_null, lower.tail=FALSE))
	}else{
		errorFunction <- function(x,constantProbability,df_alt,df_null)
		{
			sum(sqrt((qchisq(c(0.0001,0.001,0.01,0.05,0.10,0.25),df_alt-df_null,lower.tail=FALSE)*constantProbability-qchisq(c(0.0001,0.001,0.01,0.05,0.10,0.25),x,lower.tail=FALSE))^2))
		}
		deltaDF <- optimize(errorFunction,lower=0,upper=10^3,constantProbability=constantProbability,df_alt=df_alt,df_null=df_null)$minimum

		pchisq(D, deltaDF, lower.tail=FALSE)
	}
}


