#' @rdname plotGene
#'
#' @description
#' A method to see the shapes of the estimated synthesis, degradation and processing rates, pre-RNA and total RNA 
#' concentrations (solid thin lines) their variances (dashed lines) and the modeled rates and concentrations 
#' (ticker solid line) of a single gene. 
#' @param object An object of class INSPEcT
#' @param ix Eiher a rowname or a row number to select one single gene
#' @param fix.yaxis A logical, indicating whether the limits for y-axis of degradation and processing rates should be fixed
#' relative to their distributions
#' @param priors A logical, if true the priors of the rates are plotted
#' @return A list containing total RNA levels and their confidence interval (levels plus and minus
#' one standard deviation), pre-RNA lelevs and their confidence intervals, synthsis rates and 
#' their confidence intervals, degradation rates and processing rates of the selected gene.
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' plotGene(nascentInspObj10, 1)
setMethod('plotGene', 'INSPEcT', function(object, ix, fix.yaxis=FALSE, priors=TRUE) {

	ix <- ix[1]
	tpts <- object@tpts

	oneGene <- object[ix]

	ratesFirstGuessTotalTmp <- ratesFirstGuess(oneGene, 'total')
	ratesFirstGuessPreTmp <- ratesFirstGuess(oneGene, 'preMRNA')

	ratesFirstGuessTotalVarTmp <- ratesFirstGuessVar(oneGene, 'total')
	ratesFirstGuessPreVarTmp <- ratesFirstGuessVar(oneGene, 'preMRNA')

	ratesFirstGuessSynthesisTmp <- ratesFirstGuess(oneGene, 'synthesis')
	ratesFirstGuessProcessingTmp <- ratesFirstGuess(oneGene, 'processing')
	ratesFirstGuessDegradationTmp <- ratesFirstGuess(oneGene, 'degradation')

	ratesFirstGuessSynthesisVarTmp <- ratesFirstGuessVar(oneGene, 'synthesis')

	if(!priors)
	{
		ratesFirstGuessProcessingTmp[is.finite(ratesFirstGuessProcessingTmp)] <- NaN
		ratesFirstGuessDegradationTmp[is.finite(ratesFirstGuessDegradationTmp)] <- NaN

		if(object@NoNascent)
		{
			ratesFirstGuessSynthesisTmp[is.finite(ratesFirstGuessSynthesisTmp)] <- NaN
			ratesFirstGuessSynthesisVarTmp[is.finite(ratesFirstGuessSynthesisVarTmp)] <- NaN			
		}
	}

	if(!object@NoNascent)
	{
		alpha_left <- viewConfidenceIntervals(object,"synthesis_left")[ix,]
		alpha_right <- viewConfidenceIntervals(object,"synthesis_right")[ix,]

		gamma_left <- viewConfidenceIntervals(object,"processing_left")[ix,]
		gamma_right <- viewConfidenceIntervals(object,"processing_right")[ix,]

		beta_left <- viewConfidenceIntervals(object,"degradation_left")[ix,]
		beta_right <- viewConfidenceIntervals(object,"degradation_right")[ix,]
	}
	if( length(object@model@ratesSpecs) > 0 ) {

			total <- t(rbind(
				ratesFirstGuessTotalTmp
				, ratesFirstGuessTotalTmp + 
					sqrt(ratesFirstGuessTotalVarTmp)
				, ratesFirstGuessTotalTmp - 
					sqrt(ratesFirstGuessTotalVarTmp)
				, viewModelRates(oneGene, 'total')
				))
			preMRNA <- t(rbind(
				ratesFirstGuessPreTmp
				, ratesFirstGuessPreTmp + 
					sqrt(ratesFirstGuessPreVarTmp)
				, ratesFirstGuessPreTmp - 
					sqrt(ratesFirstGuessPreVarTmp)
				, viewModelRates(oneGene, 'preMRNA')
				))
			if(!object@NoNascent)
			{
				alpha <- t(rbind(ratesFirstGuessSynthesisTmp
						, alpha_left#ratesFirstGuessSynthesisTmp + sqrt(ratesFirstGuessSynthesisVarTmp)
						, alpha_right#ratesFirstGuessSynthesisTmp - sqrt(ratesFirstGuessSynthesisVarTmp)
						, viewModelRates(oneGene, 'synthesis')
				))

				gamma <- t(rbind(ratesFirstGuessProcessingTmp
						, gamma_left#ratesFirstGuessSynthesisTmp + sqrt(ratesFirstGuessSynthesisVarTmp)
						, gamma_right#ratesFirstGuessSynthesisTmp - sqrt(ratesFirstGuessSynthesisVarTmp)
						, viewModelRates(oneGene, 'processing')
				))

				beta <- t(rbind(ratesFirstGuessDegradationTmp
						, beta_left#ratesFirstGuessSynthesisTmp + sqrt(ratesFirstGuessSynthesisVarTmp)
						, beta_right#ratesFirstGuessSynthesisTmp - sqrt(ratesFirstGuessSynthesisVarTmp)
						, viewModelRates(oneGene, 'degradation')
				))

			}else{
				alpha <- t(rbind(
					ratesFirstGuessSynthesisTmp
					, viewModelRates(oneGene, 'synthesis')
					))			
				beta <- t(rbind(
					ratesFirstGuessDegradationTmp
					, viewModelRates(oneGene, 'degradation')
					))
				gamma <- t(rbind(
					ratesFirstGuessProcessingTmp
					, viewModelRates(oneGene, 'processing')
					))
			}

		if( fix.yaxis ) {

			degradationYlim <- quantile(ratesFirstGuess(object, 'degradation'), probs=c(.02, .98), na.rm=TRUE)
			processingYlim <- quantile(ratesFirstGuess(object, 'processing'), probs=c(.02, .98), na.rm=TRUE)

			x <- seq_along(tpts)

			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3), col=1, main='total RNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)

			functionTmp <- function()
			{
				matplot(x, preMRNA, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3), col=2, main='pre-RNA', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)
			}
			tryCatch(functionTmp(),error=function(e)NaN)

			if(object@NoNascent)
			{
				matplot(x, alpha, type='l', lty=c(1,1), lwd=c(1,3)
					, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')		
				axis(1, at=x, labels=signif(tpts, 2), las=3)

				matplot(x, gamma, type='l', lty=c(1,1), lwd=c(1,3)
					, col=5, main='processing', xaxt='n', xlab='time', ylab='')		
				axis(1, at=x, labels=signif(tpts, 2), las=3)

				matplot(x, beta, type='l', lty=c(1,1), lwd=c(1,3)
					, col=4, main='degradation', xaxt='n', xlab='time', ylab='')		
				axis(1, at=x, labels=signif(tpts, 2), las=3)
			}else{
				matplot(x, alpha, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
					, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')		
				axis(1, at=x, labels=signif(tpts, 2), las=3)

				functionTmp <- function()
				{
					matplot(x, gamma, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3), col=5
						, main='processing', xaxt='n', xlab='time', ylab='', ylim=processingYlim)
					axis(1, at=x, labels=signif(tpts, 2), las=3)		
				}
				tryCatch(functionTmp(),error=function(e)NaN)

				matplot(x, beta, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3), col=4
					, main='degradation', xaxt='n', xlab='time', ylab='', ylim=degradationYlim)
				axis(1, at=x, labels=signif(tpts, 2), las=3)
			}
		} else {

			x <- seq_along(tpts)

			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=1, main='total RNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)

			functionTmp <- function()
			{
				matplot(x, preMRNA, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
					, col=2, main='pre-RNA', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)			
			}
			tryCatch(functionTmp(),error=function(e)NaN)

			if(object@NoNascent)
			{
				matplot(x, alpha, type='l', lty=c(1,1), lwd=c(1,3)
					, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)

				matplot(x, gamma, type='l', lty=c(1,1), lwd=c(1,3)
					, col=5, main='processing', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)

				matplot(x, beta, type='l', lty=c(1,1), lwd=c(1,3)
					, col=4, main='degradation', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)

			}else{
				matplot(x, alpha, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
					, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')		
				axis(1, at=x, labels=signif(tpts, 2), las=3)

				functionTmp <- function()
				{
					matplot(x, gamma, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3), col=5
						, main='processing', xaxt='n', xlab='time', ylab='')
					axis(1, at=x, labels=signif(tpts, 2), las=3)		
				}

				tryCatch(functionTmp(),error=function(e)NaN)

				matplot(x, beta, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3), col=4
					, main='degradation', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)
			}
		}

	} else {

		total <- t(rbind(
			ratesFirstGuessTotalTmp
			, ratesFirstGuessTotalTmp + 
				sqrt(ratesFirstGuessTotalVarTmp)
			, ratesFirstGuessTotalTmp - 
				sqrt(ratesFirstGuessTotalVarTmp)
			))
		preMRNA <- t(rbind(
			ratesFirstGuessPreTmp
			, ratesFirstGuessPreTmp + 
				sqrt(ratesFirstGuessPreVarTmp)
			, ratesFirstGuessPreTmp - 
				sqrt(ratesFirstGuessPreVarTmp)
			))
		if(!object@NoNascent)
		{
			alpha <- t(rbind(
				ratesFirstGuessSynthesisTmp
				, ratesFirstGuessSynthesisTmp + 
					sqrt(ratesFirstGuessSynthesisVarTmp)
				, ratesFirstGuessSynthesisTmp - 
					sqrt(ratesFirstGuessSynthesisVarTmp)
				))
		}else{
		alpha <- t(rbind(
			ratesFirstGuessSynthesisTmp
			))			
		}
		beta <- t(rbind(
			ratesFirstGuessDegradationTmp
			))
		gamma <- t(rbind(
			ratesFirstGuessProcessingTmp
			))

		if( fix.yaxis ) {

			degradationYlim <- quantile(ratesFirstGuess(object, 'degradation'), probs=c(.02, .98), na.rm=TRUE)
			processingYlim <- quantile(ratesFirstGuess(object, 'processing'), probs=c(.02, .98), na.rm=TRUE)
			
			x <- seq_along(tpts)

			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=1, main='total RNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			
			functionTmp() <- function()
			{
			matplot(x, preMRNA, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=2, main='pre-RNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)				
			}
			tryCatch(functionTmp(),error=function(e)NaN)

			matplot(x, alpha, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)

			functionTmp <- function()
			{
				matplot(x, gamma, type='l', lty=c(1), lwd=c(1), col=5
					, main='processing', xaxt='n', xlab='time', ylab='', ylim=processingYlim)
				axis(1, at=x, labels=signif(tpts, 2), las=3)		
			}
			tryCatch(functionTmp(),error=function(e)NaN)

			matplot(x, beta, type='l', lty=c(1), lwd=c(1), col=4
				, main='degradation', xaxt='n', xlab='time', ylab='', ylim=degradationYlim)
			axis(1, at=x, labels=signif(tpts, 2), las=3)
		} else {
			
			x <- seq_along(tpts)

			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=1, main='total RNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)

			functionTmp <- function()
			{
				matplot(x, preMRNA, type='l', lty=c(1,2,2), lwd=c(1,1,1)
					, col=2, main='pre-RNA', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)
			}
			tryCatch(functionTmp(),error=function(e)NaN)

			matplot(x, alpha, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)

			functionTmp <- function()
			{
				matplot(x, gamma, type='l', lty=c(1), lwd=c(1), col=5
					, main='processing', xaxt='n', xlab='time', ylab='')
				axis(1, at=x, labels=signif(tpts, 2), las=3)		
			}

			tryCatch(functionTmp(),error=function(e)NaN)

			matplot(x, beta, type='l', lty=c(1), lwd=c(1), col=4
				, main='degradation', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
		}

	}

	out <- list(total=total, premature=preMRNA, synthesis=alpha, processing=gamma, degradation=beta)

	})
