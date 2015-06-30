#' @rdname rocThresholds
#'
#' @description
#' A method to visualize the performance in the classification of synthesis, degradation
#' and processing rates based on the comparison of the original simulated rates and the one
#' obtained by the function \code{\link{modelRates}}. For each rate, classification performance is measured 
#' in terms of sensitivity and specificity using a ROC curve analysis. False negatives (FN) represent cases 
#' where the rate is identified as constant while it was simulated as varying. False positives (FP) represent 
#' cases where INSPEcT identified a rate as varying while it was simulated as constant. On the contrary, 
#' true positives (TP) and negatives (TN) are cases of correct classification of varying and constant rates, respectively. 
#' Consequently, at increasing brown p-values different sensitivity and specificity can be achieved.
#' @param object An object of class INSPEcT_model, with true rates
#' @param object2 An object of class INSPEcT or INSPEcT_model, with modeled rates
#' @param cTsh A numeric representing the threshold for the chi-squared test to consider a model as valid
#' @param bTsh A numeric representing the threshold for the Brown's method to consider a rate as varying
#' @param xlim A numeric representing limits for the x-axis (default is c(1-e-5,1))
#' @return None
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}, \code{\link{rocCurve}}
#' @examples
#' data('simRates', package='INSPEcT')
#' data('simData3rep', package='INSPEcT')
#' rocThresholds(simRates, simData3rep)
#' # Increase the Brown threshold for all rates (be more relaxed)
#' thresholds(simData3rep)$brown <- c(alpha=.05, beta=.05, gamma=.05)
#' rocThresholds(simRates, simData3rep)
setMethod('rocThresholds', signature(object='INSPEcT_model', object2='INSPEcT_model'), 
	function(object, object2, cTsh=NULL, bTsh=NULL, xlim=c(1e-5,1)) {
	if( is.null(bTsh) ) {
		bTsh <- object2@params$thresholds$brown
	} else {
		if( !is.numeric(bTsh) )
			stop('rocThresholds: bTsh must be either numeric or NULL')
		if( length(bTsh) != 3 )
			stop('rocThresholds: bTsh must be a vector of length 3')
		if( any(bTsh < 0 | bTsh > 1) )
			stop('rocThresholds: every element of bTsh must be between 0 and 1')
	}
	if( is.null(cTsh) ) {
		cTsh <- object2@params$thresholds$chisquare
	} else {
		if( !is.numeric(cTsh) )
			stop('rocThresholds: cTsh must be either numeric or NULL')
		if( length(cTsh) != 1 )
			stop('rocThresholds: cTsh must be a vector of length 1')
		if( cTsh<0 | cTsh>1 )
			stop('rocThresholds: cTsh must be between 0 and 1')
	}
	if( !is.numeric(xlim) )
		stop('rocThresholds: xlim must be either numeric or NULL')
	if( length(xlim) != 2 )
		stop('rocThresholds: xlim must be a vector of length 2')
	if( any(xlim<=0 | xlim>1) )
		stop('rocThresholds: xlim must be between 0 and 1')
	## obtain the response
	allResponses <- geneClass(object)
	ratePvals <- ratePvals(object2, cTsh)
	rAlpha <- roc(response=grepl('a', allResponses)
		, predictor=ratePvals$synthesis)
	rBeta <- roc(response=grepl('b', allResponses)
		, predictor=ratePvals$degradation)
	rGamma <- roc(response=grepl('c', allResponses)
		, predictor=ratePvals$processing)
	par(mfrow=c(1,3))
	ix <- is.finite(rAlpha$thresholds)
	matplot(rAlpha$thresholds[ix], cbind(rAlpha$sensitivities[ix], rAlpha$specificities[ix])
		, type='l', lty=1, lwd=4, main='synthesis rate', log='x', xlab='Brown threshold'
		, ylab='', xlim=xlim)
	abline(v=bTsh[1], lwd=3, lty=2, col='blue')
	ix <- is.finite(rBeta$thresholds)
	matplot(rBeta$thresholds[ix], cbind(rBeta$sensitivities[ix], rBeta$specificities[ix])
		, type='l', lty=1, lwd=4, main='degradation rate', log='x', xlab='Brown threshold'
		, ylab='', xlim=xlim)
	abline(v=bTsh[2], lwd=3, lty=2, col='blue')
	ix <- is.finite(rGamma$thresholds)
	matplot(rGamma$thresholds[ix], cbind(rGamma$sensitivities[ix], rGamma$specificities[ix])
		, type='l', lty=1, lwd=4, main='processing rate', log='x', xlab='Brown threshold'
		, ylab='', xlim=xlim)
	abline(v=bTsh[3], lwd=3, lty=2, col='blue')
	legend('left', col=c('blue','black','red'), lty=c(2,1,1), lwd=c(2,4,4)
		, legend=c('selected Brown\nthreshold', 'sensitivities', 'specificities'))
	})

#' @rdname rocThresholds
setMethod('rocThresholds', signature(object='INSPEcT_model', object2='INSPEcT'),
	function(object, object2, cTsh=NULL, bTsh=NULL, xlim=c(1e-5,1)) {
		rocThresholds(object, object2@model, cTsh, bTsh, xlim)
	# if( is.null(bTsh) ) {
	# 	bTsh <- object2@model@params$thresholds$brown
	# } else {
	# 	if( !is.numeric(bTsh) )
	# 		stop('rocThresholds: bTsh must be either numeric or NULL')
	# 	if( length(bTsh) != 3 )
	# 		stop('rocThresholds: bTsh must be a vector of length 3')
	# 	if( any(bTsh < 0 | bTsh > 1) )
	# 		stop('rocThresholds: every element of bTsh must be between 0 and 1')
	# }
	# if( is.null(cTsh) ) {
	# 	cTsh <- object2@model@params$thresholds$chisquare
	# } else {
	# 	if( !is.numeric(cTsh) )
	# 		stop('rocThresholds: cTsh must be either numeric or NULL')
	# 	if( length(cTsh) != 1 )
	# 		stop('rocThresholds: cTsh must be a vector of length 1')
	# 	if( cTsh<0 | cTsh>1 )
	# 		stop('rocThresholds: cTsh must be between 0 and 1')
	# }
	# if( !is.numeric(xlim) )
	# 	stop('rocThresholds: xlim must be either numeric or NULL')
	# if( length(xlim) != 2 )
	# 	stop('rocThresholds: xlim must be a vector of length 2')
	# if( any(xlim<=0 | xlim>1) )
	# 	stop('rocThresholds: xlim must be between 0 and 1')
	# ## obtain the response
	# allResponses <- geneClass(object)
	# ratePvals <- ratePvals(object2@model, cTsh)
	# rAlpha <- roc(response=grepl('a', allResponses)
	# 	, predictor=ratePvals$synthesis)
	# rBeta <- roc(response=grepl('b', allResponses)
	# 	, predictor=ratePvals$degradation)
	# rGamma <- roc(response=grepl('c', allResponses)
	# 	, predictor=ratePvals$processing)
	# par(mfrow=c(1,3))
	# ix <- is.finite(rAlpha$thresholds)
	# matplot(rAlpha$thresholds[ix], cbind(rAlpha$sensitivities[ix], rAlpha$specificities[ix])
	# 	, type='l', lty=1, lwd=4, main='synthesis rate', log='x', xlab='Brown threshold'
	# 	, ylab='', xlim=xlim)
	# abline(v=bTsh[1], lwd=3, lty=2, col='blue')
	# ix <- is.finite(rBeta$thresholds)
	# matplot(rBeta$thresholds[ix], cbind(rBeta$sensitivities[ix], rBeta$specificities[ix])
	# 	, type='l', lty=1, lwd=4, main='degradation rate', log='x', xlab='Brown threshold'
	# 	, ylab='', xlim=xlim)
	# abline(v=bTsh[2], lwd=3, lty=2, col='blue')
	# ix <- is.finite(rGamma$thresholds)
	# matplot(rGamma$thresholds[ix], cbind(rGamma$sensitivities[ix], rGamma$specificities[ix])
	# 	, type='l', lty=1, lwd=4, main='processing rate', log='x', xlab='Brown threshold'
	# 	, ylab='', xlim=xlim)
	# abline(v=bTsh[3], lwd=3, lty=2, col='blue')
	# legend('left', col=c('blue','black','red'), lty=c(2,1,1), lwd=c(2,4,4)
	# 	, legend=c('selected Brown\nthreshold', 'sensitivities', 'specificities'))
	})
