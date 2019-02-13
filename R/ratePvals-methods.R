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
#' @param bTsh A numeric representing the p-value threshold for considering a rate as variable. P-values are calculated through \code{\link{ratePvals}}
#' @param cTsh A numeric representing the threshold for the chi-squared test to consider a model as valid
#' @details ratePvlas retrieve a single p-value for each rate thanks to multiple log likelihood tests performed on 
#' nested models that has a chi-squared test below the selected threshold. 
#' Among the many p-values that log likelihood ratio test calculate, a single p-value is obtaied applying Brown's method 
#' for combining dependent p-values.
#' @return A matrix containing p-values calculated for each rate
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' ratePvals(nascentInspObj10)
#' # calculate agin the p-values with Brown with a different threshold 
#' # for considering a model valid for the log likelihood ratio test
#' ratePvals(nascentInspObj10, cTsh=.2)
#' # Set permaenently the chi-squared threshold at .2 for nascentInspObj10 object
#' modelSelection(nascentInspObj10)$thresholds$chisquare <- .2
setMethod('ratePvals', 'INSPEcT_model', function(object, bTsh=NULL, cTsh=NULL) {
	## calculates the pval to be varying per rate per gene, 
	## according to the threshold set for the chisq masking step)

	## in case a rate has a threshold set to 0, avoid testing it
	if( is.null(bTsh) )
		bTsh <- object@params$thresholds$brown
	rates_to_aviod <- names(bTsh)[bTsh == 0]
	rates_to_aviod <- c('synthesis'='a','degradation'='b','processing'='c')[rates_to_aviod]
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
		browser()
		alphaCols <- llrtests$synthesis
		if( length(rates_to_aviod)>0 )
			alphaCols <- alphaCols[!sapply(alphaCols, function(x) any(grepl(rates_to_aviod, x)))]
		if( length(alphaCols)>0 ) {
			alphaLLRtestPvlas <- sapply(alphaCols, function(compare) {
				null <- compare[1]; alt <- compare[2]
				sapply(ratesSpecs, function(x) 
					tryCatch(logLikRatioTestInscpectModels(x[[null]], x[[alt]]), error=function(e) NA)
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
			alphaChisqMask <- alphaChisqMask[,which(apply(alphaLLRtestPvlas,2,var,na.rm=T)!=0)]
			alphaLLRtestPvlas <- alphaLLRtestPvlas[,which(apply(alphaLLRtestPvlas,2,var,na.rm=T)!=0)]
		} else {
			alphaLLRtestPvlas <- rep(1, length(ratesSpecs))
		}
		########
		# beta #
		########
		betaCols <- llrtests$degradation
		if( length(rates_to_aviod)>0 )
			betaCols <- betaCols[!sapply(betaCols, function(x) any(grepl(rates_to_aviod, x)))]
		## make the logLikRatio tests
		if( length(betaCols)>0 ) {
			betaLLRtestPvlas <- sapply(betaCols, function(compare) {
				null <- compare[1]; alt <- compare[2]
				sapply(ratesSpecs, function(x) 
					tryCatch(logLikRatioTestInscpectModels(x[[null]], x[[alt]]), error=function(e) NA)
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

			betaChisqMask <- betaChisqMask[,which(apply(betaLLRtestPvlas,2,var,na.rm=T)!=0)]
			betaLLRtestPvlas <- betaLLRtestPvlas[,which(apply(betaLLRtestPvlas,2,var,na.rm=T)!=0)]
		} else {
			betaLLRtestPvlas <- rep(1, length(ratesSpecs))
		}
		#########
		# gamma #
		#########
		gammaCols <- llrtests$processing
		if( length(rates_to_aviod)>0 )
			gammaCols <- gammaCols[!sapply(gammaCols, function(x) any(grepl(rates_to_aviod, x)))]
		## make the logLikRatio tests
		if( length(gammaCols)>0 ) {
			gammaLLRtestPvlas <- sapply(gammaCols, function(compare) {
				null <- compare[1]; alt <- compare[2]
				sapply(ratesSpecs, function(x) 
					tryCatch(logLikRatioTestInscpectModels(x[[null]], x[[alt]]), error=function(e) NA)
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
			gammaChisqMask <- gammaChisqMask[,which(apply(gammaLLRtestPvlas,2,var,na.rm=T)!=0)]
			gammaLLRtestPvlas <- gammaLLRtestPvlas[,which(apply(gammaLLRtestPvlas,2,var,na.rm=T)!=0)]
		} else {
			gammaLLRtestPvlas <- rep(1, length(ratesSpecs))
		} 
		#If we have a matrix of pValues we proceed with the Brown's test, otherwise it is useless
		if(is.matrix(alphaLLRtestPvlas)){synthesisBP <- brown_method_mask(alphaLLRtestPvlas, alphaChisqMask)}
		else{synthesisBP <- alphaLLRtestPvlas}
		if(is.matrix(betaLLRtestPvlas)){degradationBP <- brown_method_mask(betaLLRtestPvlas, betaChisqMask)}
		else{degradationBP <- betaLLRtestPvlas}
		if(is.matrix(gammaLLRtestPvlas)){processingBP <- brown_method_mask(gammaLLRtestPvlas, gammaChisqMask)}
		else{processingBP <- gammaLLRtestPvlas}
		if(object@params$padj)
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
	} else if( object@params$modelSelection=='aic' ) {
		## assign the chi-squared test result of the model selected by AIC
		## to the respective variable rates
		diz <- c('KKK', 'VKK', 'KVK', 'KKV', 'VVK', 'VKV', 'KVV', 'VVV') 
		ord = c('0','a','b','c','ac','ab','bc','abc')
		aictest <- AIC(object)[,ord]
		log_liks <- logLik(object)[,ord]
		colnames(aictest) = colnames(log_liks) = diz
		aictest[is.na(aictest)] <- Inf
		if( length(rates_to_aviod)>0 )
			aictest = aictest[,grep(rates_to_aviod, ord, invert=TRUE)]
		gene_class <- diz[apply(aictest, 1, which.min)]
		ratePvals <- pvals_from_aic(gene_class, log_liks)
		colnames(ratePvals) <- c('synthesis','processing','degradation')
		rownames(ratePvals) <- rownames(aictest)
		if( length(rates_to_aviod)>0 ) {
			ratePvals[,c('a'='synthesis','b'='degradation','c'='processing')[rates_to_aviod]] <- 1
		}
 		if(object@params$padj) {
			ratePvals <- apply(ratePvals, 2, p.adjust, method="BH")
		}
	}
	# return
	return(ratePvals)
	})

pvals_from_aic <- function(gene_class, log_liks) {

	rates_pvals <- t(sapply(seq_along(gene_class), function(i) {
		switch(gene_class[i],
			'KKK' = c(
				k1=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'VKK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KKV'], 1, lower.tail=FALSE)
				),
			'VKK' = c(
				k1=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'VKK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE)
				),
			'KVK' = c(
				k1=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE)
				),
			'KKV' = c(
				k1=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KKV'], 1, lower.tail=FALSE)
				),
			'VVK' = c(
				k1=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VVK'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE)
				),
			'VKV' = c(
				k1=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE)
				),
			'KVV' = c(
				k1=pchisq(- 2*log_liks[i,'KVV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE)
				),
			'VVV' = c(
				k1=pchisq(- 2*log_liks[i,'KVV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VVK'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE)
				)
			)
		}))
	return(rates_pvals)

}


#' @rdname ratePvals
setMethod('ratePvals', 'INSPEcT', function(object, cTsh=NULL) {
	return(ratePvals(object@model, cTsh))
	})

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

logLikRatioTestInscpectModels <- function(null, alt, dfmax=Inf)
{
	D <- - 2*null$logLik + 2*alt$logLik
	df_null <- min(dfmax, null$alpha$df) + min(dfmax, null$beta$df) +
		ifelse(!is.null(min(dfmax, null$gamma$df)), min(dfmax, null$gamma$df) , 0)
	df_alt <- min(dfmax, alt$alpha$df) + min(dfmax, alt$beta$df) +
		ifelse(!is.null(min(dfmax, alt$gamma$df)), min(dfmax, alt$gamma$df) , 0)
	#chisq.test.inspect(D,  df_alt - df_null)
	pchisq(D, df_alt-df_null, lower.tail=FALSE)	
}
