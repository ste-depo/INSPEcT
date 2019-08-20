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
#' @param object2 An modeled object of class INSPEcT
#' @param cTsh A numeric representing the threshold for the chi-squared test to consider a model as valid;
#' if NULL the value is taken from the INSPEcT_model object
#' @param plot A logical indicating whether ROC curves should be plotted or not
#' @param comparative A logical indicating whether the cross-prediction should be visualized. When this mode is selected, 
#' the p-values assigned to the variability of one rate (e.g. synthesis) are tested against the variability the other rates 
#' (e.g. processing and degradation). Cross-prediction ROC curves are plotted with dashed lines.
#' @return A list of objects of class pROC with summary of each roc curve
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}, \code{\link{rocThresholds}}
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
#'   rocCurve(simRates[1:100],nascentSim2replicates)
#'   title("3rep. 11t.p. Total and nascent RNA", line=3)
#' }

setMethod('rocCurve', signature(object='INSPEcT_model', object2='INSPEcT'), 
	function(object, object2, cTsh=NULL, plot=TRUE, comparative=FALSE) {
	if(!comparative)
	{		
		## obtain the response
		allResponses <- geneClass(object)

		## in case the classification is based on model selction (i.e. functional NoNascent)
		## plot also the AIC information
		plotAIC <- object2@NoNascent & !object2@NF
		if( plotAIC ) {
			## Pure AIC selection
			AICsTmp <- AIC(object2)
			AICclass <- colnames(AICsTmp)[apply(AICsTmp,1,which.min)]

			aCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="a")
			bCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="b")
			cCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="c")			
		}

		ratePvals <- ratePvals(object2, cTsh)
		# AT LEAST ONE FINITE PVALUE FOR EACH CONTROL CLASS!!!
		responseAlpha <- as.numeric(grepl('a', allResponses))
		if( length(table(responseAlpha))>1 ) {
			rAlpha <- roc(response=responseAlpha, predictor=ratePvals$synthesis,direction=">")
		} else {
			rAlpha <- list(auc=NA)
		}
		responseBeta <- as.numeric(grepl('b', allResponses))
		if( length(table(responseBeta))>1 ) {
			rBeta <- roc(response=responseBeta, predictor=ratePvals$degradation,direction=">")
		} else {
			rBeta <- list(auc=NA)
		}
		responseGamma <- as.numeric(grepl('c', allResponses))
		if( length(table(responseGamma))>1 ) {
			rGamma <- roc(response=responseGamma, predictor=ratePvals$processing,direction=">")
		} else {
			rGamma <- list(auc=NA)
		}
		if( plot ) {
			legendText <- paste(
				c('synthesis', 'processing', 'degradation')
				, ' - AUC='
				, signif(c(as.numeric(rAlpha$auc), as.numeric(rGamma$auc), as.numeric(rBeta$auc)), 3)
				, sep=''
				)
			legendText <- c(legendText,'p = 0.05')
			if( plotAIC ) legendText <- c(legendText,'AIC')

			if( !is.na(rAlpha[[1]]) ) {
				plot.roc(rAlpha, col='red', lwd=4)
				alphaIdx <- which.min(abs(0.05-rAlpha[[4]]))
				points(rAlpha[[3]][alphaIdx],rAlpha[[2]][alphaIdx],pch=16,col='red',cex=2.5)
				if( plotAIC ) points(aCoordinatesAIC["specificity"],aCoordinatesAIC["sensitivity"],pch=17,col='red',cex=2.5)
			}

			if( !is.na(rGamma[[1]]) ) {
				plot.roc(rGamma, col='navy', lwd=4, add=!(is.na(rAlpha[[1]]) & is.na(rBeta[[1]])))
				gammaIdx <- which.min(abs(0.05-rGamma[[4]]))
				points(rGamma[[3]][gammaIdx],rGamma[[2]][gammaIdx],pch=16,col='navy',cex=2.5)
				if( plotAIC ) points(cCoordinatesAIC["specificity"],cCoordinatesAIC["sensitivity"],pch=17,col='navy',cex=2.5)
			}

			if( !is.na(rBeta[[1]]) ) {
				plot.roc(rBeta, col='deepskyblue', lwd=4, add=!is.na(rAlpha[[1]]))
				betaIdx <- which.min(abs(0.05-rBeta[[4]]))
				points(rBeta[[3]][betaIdx],rBeta[[2]][betaIdx],pch=16,col='deepskyblue',cex=2.5)
				if( plotAIC ) points(bCoordinatesAIC["specificity"],bCoordinatesAIC["sensitivity"],pch=17,col='deepskyblue',cex=2.5)
			}

			legend('bottomright', legend=legendText
				, col=c('red', 'navy', 'deepskyblue','grey70','grey70'), 
				lty=c(1,1,1,NA,NA), 
				pch=c(NA,NA,NA,16,17), lwd=4)
		}
		## return the roc objects
		return(list(synthesis=rAlpha, degradation=rBeta, processing=rGamma))

	} else { ########## comparative mode, show also cross-classifications

		## obtain the response
		allResponses <- geneClass(object)
		ratePvals <- ratePvals(object2)

		rAlphaAlpha <- tryCatch(roc(response=as.numeric(grepl('a', allResponses))
			, predictor=ratePvals[,"synthesis"],direction=">"), error=function(e) list(auc=NA))
		rGammaAlpha <- tryCatch(roc(response=as.numeric(grepl('c', allResponses))
			, predictor=ratePvals[,"synthesis"],direction=">"), error=function(e) list(auc=NA))
		rBetaAlpha <- tryCatch(roc(response=as.numeric(grepl('b', allResponses))
			, predictor=ratePvals[,"synthesis"],direction=">"), error=function(e) list(auc=NA))

		rAlphaGamma <- tryCatch(roc(response=as.numeric(grepl('a', allResponses))
			, predictor=ratePvals[,"processing"],direction=">"), error=function(e) list(auc=NA))
		rGammaGamma <- tryCatch(roc(response=as.numeric(grepl('c', allResponses))
			, predictor=ratePvals[,"processing"],direction=">"), error=function(e) list(auc=NA))
		rBetaGamma <- tryCatch(roc(response=as.numeric(grepl('b', allResponses))
			, predictor=ratePvals[,"processing"],direction=">"), error=function(e) list(auc=NA))

		rAlphaBeta <- tryCatch(roc(response=as.numeric(grepl('a', allResponses))
			, predictor=ratePvals[,"degradation"],direction=">"), error=function(e) list(auc=NA))
		rGammaBeta <- tryCatch(roc(response=as.numeric(grepl('c', allResponses))
			, predictor=ratePvals[,"degradation"],direction=">"), error=function(e) list(auc=NA))
		rBetaBeta <- tryCatch(roc(response=as.numeric(grepl('b', allResponses))
			, predictor=ratePvals[,"degradation"],direction=">"), error=function(e) list(auc=NA))

		if( plot ) {
	
			oldMfrow <- par()$mfrow
			par(mfrow=c(1,3))

			if( !is.na(rAlphaAlpha[[1]]) )
				plot.roc(rAlphaAlpha, col='red', lwd=4, add=FALSE, main="Synthesis rate classification")
			if( !is.na(rGammaAlpha[[1]]) ) 
				plot.roc(rGammaAlpha, col='navy', lwd=4, add=!is.na(rAlphaAlpha[[1]]), lty=3)
			if( !is.na(rBetaAlpha[[1]]) ) 
				plot.roc(rBetaAlpha, col='deepskyblue', lwd=4, add=!(is.na(rAlphaAlpha[[1]]) & is.na(rGammaAlpha[[1]])), lty=3)

			legendText <- paste(
				c('synthesis', 'processing', 'degradation')
				, ' - AUC='
				, signif(c(as.numeric(rAlphaAlpha$auc), as.numeric(rGammaAlpha$auc), as.numeric(rBetaAlpha$auc)), 3)
				, sep=''
				)

			legend('bottomright', legend=legendText
				, col=c('red', 'navy', 'deepskyblue'), lty=c(1,3,3), lwd=4)
	
			if( !is.na(rAlphaGamma[[1]]) )
				plot.roc(rAlphaGamma, col='red', lwd=4, lty=3, add=FALSE, main="Processing rate classification")
			if( !is.na(rGammaGamma[[1]]) )
				plot.roc(rGammaGamma, col='navy', lwd=4, add=!is.na(rAlphaGamma[[1]]))
			if( !is.na(rBetaGamma[[1]]) )
				plot.roc(rBetaGamma, col='deepskyblue', lwd=4, lty=3, add=!(is.na(rAlphaGamma[[1]]) & is.na(rGammaGamma[[1]])))
	
			legendText <- paste(
				c('synthesis', 'processing', 'degradation')
				, ' - AUC='
				, signif(c(as.numeric(rAlphaGamma$auc), as.numeric(rGammaGamma$auc), as.numeric(rBetaGamma$auc)), 3)
				, sep=''
				)
	
			legend('bottomright', legend=legendText
				, col=c('red', 'navy', 'deepskyblue'), lty=c(3,1,3), lwd=4)
	
			if( !is.na(rAlphaBeta[[1]]) )
				plot.roc(rAlphaBeta, col='red', lwd=4, lty=3, add=FALSE, main="Degradation rate classification")
			if( !is.na(rGammaBeta[[1]]) )
				plot.roc(rGammaBeta, col='navy', lwd=4, lty=3, add=!is.na(rAlphaBeta[[1]]))
			if( !is.na(rBetaBeta[[1]]) )
				plot.roc(rBetaBeta, col='deepskyblue', lwd=4, add=!(is.na(rAlphaBeta[[1]]) & is.na(rGammaBeta[[1]])))
	
			legendText <- paste(
				c('synthesis', 'processing', 'degradation')
				, ' - AUC='
				, signif(c(as.numeric(rAlphaBeta$auc), as.numeric(rGammaBeta$auc), as.numeric(rBetaBeta$auc)), 3)
				, sep=''
				)
	
			legend('bottomright', legend=legendText
				, col=c('red', 'navy', 'deepskyblue'), lty=c(3,3,1), lwd=4)
		
			# restore par settings
			par(mfrow=oldMfrow)
		}

		out <- list("ClassifyAlphaWithAlpha"=rAlphaAlpha
				   ,"ClassifyBetaWithAlpha"=rBetaAlpha
				   ,"ClassifyGammaWithAlpha"=rGammaAlpha
				   ,"ClassifyAlphaWithGamma"=rAlphaGamma
				   ,"ClassifyBetaWithGamma"=rBetaGamma
				   ,"ClassifyGammaWithGamma"=rGammaGamma
				   ,"ClassifyAlphaWithBeta"=rAlphaBeta
				   ,"ClassifyBetaWithBeta"=rBetaBeta
				   ,"ClassifyGammaWithBeta"=rGammaBeta)
		return(out)
	}
})

coordinatesAIC <- function(AICclass,allResponses,class)
{
	TP <- length(which(grepLogic(class,allResponses)&grepLogic(class,AICclass)))
	P <- length(grep(class,allResponses))

	TN <- length(which(!grepLogic(class,allResponses)&!grepLogic(class,AICclass)))
	N <- length(allResponses) - P

	specificityAIC <- TP/P
	sensitivityAIC <- TN/N

	c(specificity=specificityAIC,sensitivity=sensitivityAIC)
}

# rocCurve_confidenceIntervals<-function(object, object2, plot, comparative=FALSE)
# {
# 	if(!comparative)
# 	{
# 		## obtain the response
# 		allResponses <- geneClass(object)
# 		ratePvals <- ratePvals(object2)

# 		rAlpha <- roc(response=as.numeric(grepl('a', allResponses))
# 			, predictor=ratePvals[,"synthesis"],direction=">")
# 		rBeta <- roc(response=as.numeric(grepl('b', allResponses))
# 			, predictor=ratePvals[,"degradation"],direction=">")
# 		rGamma <- roc(response=as.numeric(grepl('c', allResponses))
# 			, predictor=ratePvals[,"processing"],direction=">")
# 		if( plot ) {
# 			legendText <- paste(
# 				c('synthesis', 'degradation', 'processing')
# 				, ' - AUC='
# 				, signif(c(as.numeric(rAlpha$auc), as.numeric(rBeta$auc), as.numeric(rGamma$auc)), 3)
# 				, sep=''
# 				)
# 			legendText <- c(legendText,paste(c('synthesis', 'degradation', 'processing'), ' - threshold 0.05', sep=''))
# 			plot.roc(rAlpha, col='red', lwd=4)
# 			alphaIdx <- which.min(abs(0.05-rAlpha[[4]]))
# 			points(rAlpha[[3]][alphaIdx],rAlpha[[2]][alphaIdx],pch=16,col='red',cex=2.5)

# 			plot.roc(rBeta, col='deepskyblue', lwd=4, add=TRUE)
# 			betaIdx <- which.min(abs(0.05-rBeta[[4]]))
# 			points(rBeta[[3]][betaIdx],rBeta[[2]][betaIdx],pch=16,col='deepskyblue',cex=2.5)

# 			plot.roc(rGamma, col='navy', lwd=4, add=TRUE)
# 			gammaIdx <- which.min(abs(0.05-rGamma[[4]]))
# 			points(rGamma[[3]][gammaIdx],rGamma[[2]][gammaIdx],pch=16,col='navy',cex=2.5)

# 			legend('bottomright', legend=legendText
# 				, col=rep(c('red', 'deepskyblue', 'navy'),2), lty=c(1,1,1,NA,NA,NA), pch=c(NA,NA,NA,16,16,16), lwd=4)
# 		}
# 		## return the roc objects
# 		out <- list(synthesis=rAlpha, degradation=rBeta, processing=rGamma)		
# 		return(out)
# 	}else{
# 		par(mfrow=c(1,3))

# 		## obtain the response
# 		allResponses <- geneClass(object)
# 		ratePvals <- ratePvals(object2)

# 		rAlphaAlpha <- roc(response=as.numeric(grepl('a', allResponses))
# 			, predictor=ratePvals[,"synthesis"],direction=">")
# 		rBetaAlpha <- roc(response=as.numeric(grepl('b', allResponses))
# 			, predictor=ratePvals[,"synthesis"],direction=">")
# 		rGammaAlpha <- roc(response=as.numeric(grepl('c', allResponses))
# 			, predictor=ratePvals[,"synthesis"],direction=">")

# 		rAlphaBeta <- roc(response=as.numeric(grepl('a', allResponses))
# 			, predictor=ratePvals[,"degradation"],direction=">")
# 		rBetaBeta <- roc(response=as.numeric(grepl('b', allResponses))
# 			, predictor=ratePvals[,"degradation"],direction=">")
# 		rGammaBeta <- roc(response=as.numeric(grepl('c', allResponses))
# 			, predictor=ratePvals[,"degradation"],direction=">")

# 		rAlphaGamma <- roc(response=as.numeric(grepl('a', allResponses))
# 			, predictor=ratePvals[,"processing"],direction=">")
# 		rBetaGamma <- roc(response=as.numeric(grepl('b', allResponses))
# 			, predictor=ratePvals[,"processing"],direction=">")
# 		rGammaGamma <- roc(response=as.numeric(grepl('c', allResponses))
# 			, predictor=ratePvals[,"processing"],direction=">")

# 		plot.roc(rAlphaAlpha, col='red', lwd=4, add=FALSE, main="Synthesis rate classification")
# 		plot.roc(rBetaAlpha, col='deepskyblue', lwd=4, add=TRUE, lty=2)
# 		plot.roc(rGammaAlpha, col='navy', lwd=4, add=TRUE, lty=2)

# 		legendText <- paste(
# 			c('synthesis', 'degradation', 'processing')
# 			, ' - AUC='
# 			, signif(c(as.numeric(rAlphaAlpha$auc), as.numeric(rBetaAlpha$auc), as.numeric(rGammaAlpha$auc)), 3)
# 			, sep=''
# 			)

# 		legend('bottomright', legend=legendText
# 			, col=c('red', 'deepskyblue', 'navy'), lty=c(1,2,2), lwd=4)

# 		plot.roc(rAlphaGamma, col='red', lwd=4, lty=2, add=FALSE, main="Processing rate classification")
# 		plot.roc(rBetaGamma, col='deepskyblue', lwd=4, lty=2, add=TRUE)
# 		plot.roc(rGammaGamma, col='navy', lwd=4, add=TRUE)

# 		legendText <- paste(
# 			c('synthesis', 'degradation', 'processing')
# 			, ' - AUC='
# 			, signif(c(as.numeric(rAlphaGamma$auc), as.numeric(rBetaGamma$auc), as.numeric(rGammaGamma$auc)), 3)
# 			, sep=''
# 			)

# 		legend('bottomright', legend=legendText
# 			, col=c('red', 'deepskyblue', 'navy'), lty=c(2,2,1), lwd=4)

# 		plot.roc(rAlphaBeta, col='red', lwd=4, lty=2, add=FALSE, main="Degradation rate classification")
# 		plot.roc(rBetaBeta, col='deepskyblue', lwd=4, add=TRUE)
# 		plot.roc(rGammaBeta, col='navy', lwd=4, lty=2, add=TRUE)

# 		legendText <- paste(
# 			c('synthesis', 'degradation', 'processing')
# 			, ' - AUC='
# 			, signif(c(as.numeric(rAlphaBeta$auc), as.numeric(rBetaBeta$auc), as.numeric(rGammaBeta$auc)), 3)
# 			, sep=''
# 			)

# 		legend('bottomright', legend=legendText
# 			, col=c('red', 'deepskyblue', 'navy'), lty=c(2,1,2), lwd=4)

# 		out <- list("ClassifyAlphaWithAlpha"=rAlphaAlpha
# 				   ,"ClassifyBetaWithAlpha"=rBetaAlpha
# 				   ,"ClassifyGammaWithAlpha"=rGammaAlpha
# 				   ,"ClassifyAlphaWithGamma"=rAlphaGamma
# 				   ,"ClassifyBetaWithGamma"=rBetaGamma
# 				   ,"ClassifyGammaWithGamma"=rGammaGamma
# 				   ,"ClassifyAlphaWithBeta"=rAlphaBeta
# 				   ,"ClassifyBetaWithBeta"=rBetaBeta
# 				   ,"ClassifyGammaWithBeta"=rGammaBeta)
# 		return(out)
# 	}
# }

# #' @rdname rocCurve
# setMethod('rocCurve', signature(object='INSPEcT_model', object2='INSPEcT'), function(object, object2, cTsh=NULL, plot=TRUE, comparative=FALSE)
# {
# 	# if( object2@NoNascent & !object2@NF )
# 	# {
# 	# 	out <- rocCurve(object, object2@model, cTsh, plot, comparative=comparative)
# 	# }else{
# 	# 	out <- rocCurve_confidenceIntervals(object, object2, plot, comparative=comparative)
# 	# }
# 	out <- rocCurve(object, object2@model, cTsh, plot, comparative=comparative)
# 	return(out)
# })

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
	# })
