#' @rdname rocCurve
#'
#' @description
#' A method to visualize the performance in the classification of synthesis, degradation
#' and processing rates based on the comparison of the original simulated rates and the one
#' obtained by the function \code{\link{modelRates}}. For each rate, classification performance is measured 
#' in terms of sensitivity and specificity using a ROC curve analysis. False negatives (FN) represent cases 
#' where the rate is identified as constant while it was simulated as varying. False positives (FP) represent 
#' cases where INSPEcT identified a rate as varying while it was simulated as constant. On the contrary, 
#' true positives (TP) and negatives (TN) are cases of correct classification of varying and constant rates, respectively. 
#' Consequently, sensitivity and specificity are computed using increasing thresholds for the brown p-values, 
#' and the ability of correctly classifying a rate is measured through the area under the curve (AUC) for each rate.
#' @param object An object of class INSPEcT_model, with true rates
#' @param object2 An object of class INSPEcT or INSPEcT_model, with modeled rates
#' @param cTsh A numeric representing the threshold for the chi-squared test to consider a model as valid;
#' if NULL the value is taken from the INSPEcT_model object
#' @param plot A logical indicating whether ROC curves should be plotted or not
#' @return A list of objects of class pROC with summary of each roc curve
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}, \code{\link{rocThresholds}}
#' @examples
#' data('simRates', package='INSPEcT')
#' data('simData3rep_4su', package='INSPEcT')
#' rocCurve(simRates[1:10], simData3rep_4su)
setMethod('rocCurve', signature(object='INSPEcT_model', object2='INSPEcT_model'), 
	function(object, object2, cTsh=NULL, plot=TRUE) {
	## obtain the response
	allResponses <- geneClass(object)
	ratePvals <- ratePvals(object2, cTsh)
	# AT LEAST ONE FINITE PVALUE FOR EACH CONTROL CLASS!!!
	rAlpha <- roc(response=grepl('a', allResponses)
		, predictor=ratePvals$synthesis)
	rBeta <- roc(response=grepl('b', allResponses)
		, predictor=ratePvals$degradation)
	rGamma <- roc(response=grepl('c', allResponses)
		, predictor=ratePvals$processing)
	if( plot ) {
		legendText <- paste(
			c('synthesis', 'degradation', 'processing')
			, ' - AUC='
			, signif(c(as.numeric(rAlpha$auc), as.numeric(rBeta$auc), as.numeric(rGamma$auc)), 3)
			, sep=''
			)
		plot.roc(rAlpha, col='red', lwd=4)
		plot.roc(rBeta, col='deepskyblue', lwd=4, add=TRUE)
		plot.roc(rGamma, col='navy', lwd=4, add=TRUE)
		legend('bottomright', legend=legendText
			, col=c('red', 'deepskyblue', 'navy'), lty=1, lwd=4)
	}
	## return the roc objects
	out <- list(synthesis=rAlpha, degradation=rBeta, processing=rGamma)
	})

#' @rdname rocCurve
setMethod('rocCurve', signature(object='INSPEcT_model', object2='INSPEcT'), 
	function(object, object2, cTsh=NULL, plot=TRUE) {
		out <- rocCurve(object, object2@model, cTsh, plot)
	# ## obtain the response
	# allResponses <- geneClass(object)
	# ratePvals <- ratePvals(object2@model, cTsh)
	# rAlpha <- roc(response=grepl('a', allResponses)
	# 	, predictor=ratePvals$synthesis)
	# rBeta <- roc(response=grepl('b', allResponses)
	# 	, predictor=ratePvals$degradation)
	# rGamma <- roc(response=grepl('c', allResponses)
	# 	, predictor=ratePvals$processing)
	# if( plot ) {
	# 	legendText <- paste(
	# 		c('synthesis', 'degradation', 'processing')
	# 		, ' - AUC='
	# 		, signif(c(as.numeric(rAlpha$auc), as.numeric(rBeta$auc), as.numeric(rGamma$auc)), 3)
	# 		, sep=''
	# 		)
	# 	plot.roc(rAlpha, col='red', lwd=4)
	# 	plot.roc(rBeta, col='deepskyblue', lwd=4, add=TRUE)
	# 	plot.roc(rGamma, col='navy', lwd=4, add=TRUE)
	# 	legend('bottomright', legend=legendText
	# 		, col=c('red', 'deepskyblue', 'navy'), lty=1, lwd=4)
	# }
	# ## return the roc objects
	# out <- list(synthesis=rAlpha, degradation=rBeta, processing=rGamma)
	})
