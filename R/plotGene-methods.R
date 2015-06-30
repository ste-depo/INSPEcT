#' @rdname plotGene
#'
#' @description
#' A method to see the shapes of the estimated synthesis, degradation and processing rates, pre-mRNA and total mRNA 
#' concentrations (solid thin lines) their variances (dashed lines) and the modeled rates and concentrations 
#' (ticker solid line) of a single gene. 
#' @param object An object of class INSPEcT
#' @param ix Eiher a rowname or a row number to select one single gene
#' @param fix.yaxis A logical, indicating whether the limits for y-axis of degradation and processing rates should be fixed
#' relative to their distributions
#' @return A list containing total mRNA levels and their confidence interval (levels plus and minus
#' one standard deviation), pre-mRNA lelevs and their confidence intervals, synthsis rates and 
#' their confidence intervals, degradation rates and processing rates of the selected gene.
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' plotGene(mycerIds10, 1)
setMethod('plotGene', 'INSPEcT', function(object, ix, fix.yaxis=FALSE) {

	ix <- ix[1]
	tpts <- object@tpts
	oneGene <- object[ix]

	if( length(object@model@ratesSpecs) > 0 ) {

		total <- t(rbind(
			ratesFirstGuess(oneGene, 'total')
			, ratesFirstGuess(oneGene, 'total') + 
				sqrt(ratesFirstGuessVar(oneGene, 'total'))
			, ratesFirstGuess(oneGene, 'total') - 
				sqrt(ratesFirstGuessVar(oneGene, 'total'))
			, viewModelRates(oneGene, 'total')
			))
		preMRNA <- t(rbind(
			ratesFirstGuess(oneGene, 'preMRNA')
			, ratesFirstGuess(oneGene, 'preMRNA') + 
				sqrt(ratesFirstGuessVar(oneGene, 'preMRNA'))
			, ratesFirstGuess(oneGene, 'preMRNA') - 
				sqrt(ratesFirstGuessVar(oneGene, 'preMRNA'))
			, viewModelRates(oneGene, 'preMRNA')
			))
		alpha <- t(rbind(
			ratesFirstGuess(oneGene, 'synthesis')
			, ratesFirstGuess(oneGene, 'synthesis') + 
				sqrt(ratesFirstGuessVar(oneGene, 'synthesis'))
			, ratesFirstGuess(oneGene, 'synthesis') - 
				sqrt(ratesFirstGuessVar(oneGene, 'synthesis'))
			, viewModelRates(oneGene, 'synthesis')
			))
		beta <- t(rbind(
			ratesFirstGuess(oneGene, 'degradation')
			, viewModelRates(oneGene, 'degradation')
			))
		gamma <- t(rbind(
			ratesFirstGuess(oneGene, 'processing')
			, viewModelRates(oneGene, 'processing')
			))

		if( fix.yaxis ) {

			degradationYlim <- quantile(ratesFirstGuess(object, 'degradation'), probs=c(.02, .98), na.rm=TRUE)
			processingYlim <- quantile(ratesFirstGuess(object, 'processing'), probs=c(.02, .98), na.rm=TRUE)

			log_shift <- .find_tt_par(tpts)
			x <- .time_transf(tpts, log_shift)
			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=1, main='total mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, preMRNA, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=2, main='pre-mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, alpha, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, beta, type='l', lty=c(1,1), lwd=c(1,3), col=4
				, main='degradation', xaxt='n', xlab='time', ylab='', ylim=degradationYlim)
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, gamma, type='l', lty=c(1,1), lwd=c(1,3), col=5
				, main='processing', xaxt='n', xlab='time', ylab='', ylim=processingYlim)
			axis(1, at=x, labels=signif(tpts, 2), las=3)		
		} else {
			log_shift <- .find_tt_par(tpts)
			x <- .time_transf(tpts, log_shift)
			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=1, main='total mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, preMRNA, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=2, main='pre-mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, alpha, type='l', lty=c(1,2,2,1), lwd=c(1,1,1,3)
				, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, beta, type='l', lty=c(1,1), lwd=c(1,3), col=4
				, main='degradation', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, gamma, type='l', lty=c(1,1), lwd=c(1,3), col=5
				, main='processing', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)		
		}


	} else {

		total <- t(rbind(
			ratesFirstGuess(oneGene, 'total')
			, ratesFirstGuess(oneGene, 'total') + 
				sqrt(ratesFirstGuessVar(oneGene, 'total'))
			, ratesFirstGuess(oneGene, 'total') - 
				sqrt(ratesFirstGuessVar(oneGene, 'total'))
			))
		preMRNA <- t(rbind(
			ratesFirstGuess(oneGene, 'preMRNA')
			, ratesFirstGuess(oneGene, 'preMRNA') + 
				sqrt(ratesFirstGuessVar(oneGene, 'preMRNA'))
			, ratesFirstGuess(oneGene, 'preMRNA') - 
				sqrt(ratesFirstGuessVar(oneGene, 'preMRNA'))
			))
		alpha <- t(rbind(
			ratesFirstGuess(oneGene, 'synthesis')
			, ratesFirstGuess(oneGene, 'synthesis') + 
				sqrt(ratesFirstGuessVar(oneGene, 'synthesis'))
			, ratesFirstGuess(oneGene, 'synthesis') - 
				sqrt(ratesFirstGuessVar(oneGene, 'synthesis'))
			))
		beta <- t(rbind(
			ratesFirstGuess(oneGene, 'degradation')
			))
		gamma <- t(rbind(
			ratesFirstGuess(oneGene, 'processing')
			))

		if( fix.yaxis ) {

			degradationYlim <- quantile(ratesFirstGuess(object, 'degradation'), probs=c(.02, .98), na.rm=TRUE)
			processingYlim <- quantile(ratesFirstGuess(object, 'processing'), probs=c(.02, .98), na.rm=TRUE)

			log_shift <- .find_tt_par(tpts)
			x <- .time_transf(tpts, log_shift)
			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=1, main='total mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, preMRNA, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=2, main='pre-mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, alpha, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, beta, type='l', lty=c(1), lwd=c(1), col=4
				, main='degradation', xaxt='n', xlab='time', ylab='', ylim=degradationYlim)
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, gamma, type='l', lty=c(1), lwd=c(1), col=5
				, main='processing', xaxt='n', xlab='time', ylab='', ylim=processingYlim)
			axis(1, at=x, labels=signif(tpts, 2), las=3)		
		} else {
			log_shift <- .find_tt_par(tpts)
			x <- .time_transf(tpts, log_shift)
			par(mfrow=c(1,5))
			matplot(x, total, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=1, main='total mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, preMRNA, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=2, main='pre-mRNA', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, alpha, type='l', lty=c(1,2,2), lwd=c(1,1,1)
				, col=3, main='synthesis', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, beta, type='l', lty=c(1), lwd=c(1), col=4
				, main='degradation', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)
			matplot(x, gamma, type='l', lty=c(1), lwd=c(1), col=5
				, main='processing', xaxt='n', xlab='time', ylab='')
			axis(1, at=x, labels=signif(tpts, 2), las=3)		
		}

	}

	out <- list(total, preMRNA, alpha, beta, gamma)

	})