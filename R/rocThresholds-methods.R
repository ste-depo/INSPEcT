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
#' @param xlim A numeric representing limits for the x-axis (default is c(1-e-5,1))
#' @param plot A logical that indicates whether to plot or not. (default=TRUE)
#' @return The thresholds that maximize both sensitivity and specificity
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}, \code{\link{rocCurve}}
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   nascentInspObj <- readRDS(system.file(package='INSPEcT', 'nascentInspObj.rds'))
#'  
#'   simRates<-makeSimModel(nascentInspObj, 1000, seed=1)
#'    
#'   # newTpts<-simRates@params$tpts
#'   # nascentSim2replicates<-makeSimDataset(object=simRates
#'   #                                    ,tpts=newTpts
#'   #                                    ,nRep=3
#'   #                                    ,NoNascent=FALSE
#'   #                                    ,seed=1)
#'   # nascentSim2replicates<-modelRates(nascentSim2replicates[1:100]
#'   #                                ,seed=1)
#'   # (not evaluated to save computational time)
#'  
#'   data("nascentSim2replicates",package='INSPEcT')
#'  
#'   rocThresholds(simRates[1:100],nascentSim2replicates)
#' }
setMethod('rocThresholds', signature(object='INSPEcT_model', object2='INSPEcT'), function(object, object2, xlim=c(1e-5,1), plot=TRUE)
{
	if( !.hasSlot(object2, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	# reduce the set of genes to the one present in the modeled object
	if( length(object@ratesSpecs) !=  length(featureNames(object2)) )
		object <- object[as.numeric(featureNames(object2))]
	# manage arguments
	if( !is.numeric(xlim) )
		stop('rocThresholds: xlim must be either numeric or NULL')
	if( length(xlim) != 2 )
		stop('rocThresholds: xlim must be a vector of length 2')
	if( any(xlim<=0 | xlim>1) )
		stop('rocThresholds: xlim must be between 0 and 1')
	## obtain the response (true class)
	allResponses <- geneClass(object)
	## compare with p-values of rate variability
	ratePvals <- ratePvals(object2)
	rAlpha <- roc(response=grepl('a', allResponses)
		, predictor=ratePvals$synthesis)
	rBeta <- roc(response=grepl('b', allResponses)
		, predictor=ratePvals$degradation)
	rGamma <- roc(response=grepl('c', allResponses)
		, predictor=ratePvals$processing)

	outTmp <- list()

	ix <- is.finite(rAlpha$thresholds)
	ths <- rAlpha$thresholds[ix]
	outTmp$synthesis <- ths[which.min(sqrt((rAlpha$sensitivities[ix]-rAlpha$specificities[ix])^2))]

	ix <- is.finite(rGamma$thresholds)
	ths <- rGamma$thresholds[ix]
	outTmp$processing <- ths[which.min(sqrt((rGamma$sensitivities[ix]-rGamma$specificities[ix])^2))]

	ix <- is.finite(rBeta$thresholds)
	ths <- rBeta$thresholds[ix]
	outTmp$degradation <- ths[which.min(sqrt((rBeta$sensitivities[ix]-rBeta$specificities[ix])^2))]

	if(plot)
	{

		bTsh <- modelSelection(object2)$thresholds$brown
		
		oldMfrow <- par()$mfrow
		oldMar <- par()$mar
	
		par(mfrow=c(1,3))
		ix <- is.finite(rAlpha$thresholds)
		matplot(rAlpha$thresholds[ix], cbind(rAlpha$sensitivities[ix], rAlpha$specificities[ix])
			, type='l', lty=1, lwd=4, main='synthesis rate', log='x', xlab='Brown threshold'
			, ylab='', xlim=xlim)
		abline(v=bTsh[1], lwd=3, lty=2, col='blue')
		abline(v=outTmp$synthesis, lwd=3, lty=3, col="blue")

		ix <- is.finite(rBeta$thresholds)
		matplot(rBeta$thresholds[ix], cbind(rBeta$sensitivities[ix], rBeta$specificities[ix])
			, type='l', lty=1, lwd=4, main='degradation rate', log='x', xlab='Brown threshold'
			, ylab='', xlim=xlim)
		abline(v=bTsh[2], lwd=3, lty=2, col='blue')
		abline(v=outTmp$degradation, lwd=3, lty=3, col="blue")

		ix <- is.finite(rGamma$thresholds)
		matplot(rGamma$thresholds[ix], cbind(rGamma$sensitivities[ix], rGamma$specificities[ix])
			, type='l', lty=1, lwd=4, main='processing rate', log='x', xlab='Brown threshold'
			, ylab='', xlim=xlim)
		abline(v=bTsh[3], lwd=3, lty=2, col='blue')
		abline(v=outTmp$processing, lwd=3, lty=3, col="blue")

		legend('left', col=c('blue','blue','black','red'), lty=c(3,2,1,1), lwd=c(2,2,4,4)
			, legend=c('optimal Brown\nthreshold','selected Brown\nthreshold', 'sensitivities', 'specificities'))

		par(mfrow=oldMfrow, mar=oldMar)
		
		}

		return(unlist(outTmp))
})