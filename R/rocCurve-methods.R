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
	function(object, object2, plot=TRUE, comparative=FALSE) {
		if( !.hasSlot(object2, 'version') ) {
			stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
		}
		
		# reduce the set of genes to the one present in the modeled object
		if( length(object@ratesSpecs) !=  length(featureNames(object2)) )
			object <- object[as.numeric(featureNames(object2))]
		
		if(!comparative)
		{		
			## obtain the response (true class)
			allResponses <- geneClass(object)
			
			## in case the classification is based on model selction (i.e. functional NoNascent)
			## plot also the AIC information
			plotAIC <- (object2@NoNascent & !object2@NF)
			if( plotAIC ) {
				## Pure AIC selection
				AICsTmp <- AIC_internal(object2@model)
				AICclass <- colnames(AICsTmp)[apply(AICsTmp,1,which.min)]
				
				aCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="a")
				bCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="b")
				cCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="c")			
			}
			
			ratePvals <- ratePvals(object2)
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
			out <- list(synthesis=rAlpha, degradation=rBeta, processing=rGamma)
			
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
			
			## in case the classification is based on model selction (i.e. functional NoNascent)
			## plot also the AIC information
			plotAIC <- object2@NoNascent & !object2@NF
			if( plotAIC ) {
				## Pure AIC selection
				AICsTmp <- AIC(object2)
				AICclass <- colnames(AICsTmp)[apply(AICsTmp,1,which.min)]
				
				aaCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="a")
				gaCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="a",predictor="c")
				baCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="a",predictor="b")
				
				acCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="c",predictor="a")
				ccCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="c")			
				bcCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="c",predictor="b")
				
				abCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="b",predictor="a")
				cbCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="b",predictor="c")
				bbCoordinatesAIC <- coordinatesAIC(AICclass=AICclass,allResponses=allResponses,class="b")
				
			}
			
			if( plot ) {
				
				oldMfrow <- par()$mfrow
				par(mfrow=c(1,3))
				
				if( !is.na(rAlphaAlpha[[1]]) ) {
					plot.roc(rAlphaAlpha, col='red', lwd=4, add=FALSE, main="Synthesis rate classification")
					alphaAlphaIdx <- which.min(abs(0.05-rAlphaAlpha[[4]]))
					points(rAlphaAlpha[[3]][alphaAlphaIdx],rAlphaAlpha[[2]][alphaAlphaIdx],pch=16,col='red',cex=2.5)
					if( plotAIC ) points(aaCoordinatesAIC["specificity"],aaCoordinatesAIC["sensitivity"],pch=17,col='red',cex=2.5)
				}
				if( !is.na(rGammaAlpha[[1]]) ) {
					plot.roc(rGammaAlpha, col='navy', lwd=4, add=!is.na(rAlphaAlpha[[1]]), lty=3)
					gammaAlphaIdx <- which.min(abs(0.05-rGammaAlpha[[4]]))
					points(rGammaAlpha[[3]][gammaAlphaIdx],rGammaAlpha[[2]][gammaAlphaIdx],pch=16,col='navy',cex=2.5)
					if( plotAIC ) points(gaCoordinatesAIC["specificity"],gaCoordinatesAIC["sensitivity"],pch=17,col='navy',cex=2.5)
				}
				if( !is.na(rBetaAlpha[[1]]) ) {
					plot.roc(rBetaAlpha, col='deepskyblue', lwd=4, add=!(is.na(rAlphaAlpha[[1]]) & is.na(rGammaAlpha[[1]])), lty=3)
					betaAlphaIdx <- which.min(abs(0.05-rBetaAlpha[[4]]))
					points(rBetaAlpha[[3]][betaAlphaIdx],rBetaAlpha[[2]][betaAlphaIdx],pch=16,col='deepskyblue',cex=2.5)
					if( plotAIC ) points(baCoordinatesAIC["specificity"],baCoordinatesAIC["sensitivity"],pch=17,col='deepskyblue',cex=2.5)
				}
				
				legendText <- paste(
					c('synthesis', 'processing', 'degradation')
					, ' - AUC='
					, signif(c(as.numeric(rAlphaAlpha$auc), as.numeric(rGammaAlpha$auc), as.numeric(rBetaAlpha$auc)), 3)
					, sep=''
				)
				legendText <- c(legendText,'p = 0.05')
				if( plotAIC ) legendText <- c(legendText,'AIC')
				
				legend('bottomright', legend=legendText
							 , col=c('red', 'navy', 'deepskyblue','grey70','grey70'), 
							 lty=c(1,3,3,NA,NA), 
							 pch=c(NA,NA,NA,16,17), lwd=4)
				
				if( !is.na(rAlphaGamma[[1]]) ) {
					plot.roc(rAlphaGamma, col='red', lwd=4, lty=3, add=FALSE, main="Processing rate classification")
					alphaGammaIdx <- which.min(abs(0.05-rAlphaGamma[[4]]))
					points(rAlphaGamma[[3]][alphaGammaIdx],rAlphaGamma[[2]][alphaGammaIdx],pch=16,col='red',cex=2.5)
					if( plotAIC ) points(acCoordinatesAIC["specificity"],acCoordinatesAIC["sensitivity"],pch=17,col='red',cex=2.5)
				}
				if( !is.na(rGammaGamma[[1]]) ) {
					plot.roc(rGammaGamma, col='navy', lwd=4, add=!is.na(rAlphaGamma[[1]]))
					gammaGammaIdx <- which.min(abs(0.05-rGammaGamma[[4]]))
					points(rGammaGamma[[3]][gammaGammaIdx],rGammaGamma[[2]][gammaGammaIdx],pch=16,col='navy',cex=2.5)
					if( plotAIC ) points(ccCoordinatesAIC["specificity"],ccCoordinatesAIC["sensitivity"],pch=17,col='navy',cex=2.5)
				}
				if( !is.na(rBetaGamma[[1]]) ) {
					plot.roc(rBetaGamma, col='deepskyblue', lwd=4, lty=3, add=!(is.na(rAlphaGamma[[1]]) & is.na(rGammaGamma[[1]])))
					betaGammaIdx <- which.min(abs(0.05-rBetaGamma[[4]]))
					points(rBetaGamma[[3]][betaGammaIdx],rBetaGamma[[2]][betaGammaIdx],pch=16,col='deepskyblue',cex=2.5)
					if( plotAIC ) points(bcCoordinatesAIC["specificity"],bcCoordinatesAIC["sensitivity"],pch=17,col='deepskyblue',cex=2.5)
				}
				
				legendText <- paste(
					c('synthesis', 'processing', 'degradation')
					, ' - AUC='
					, signif(c(as.numeric(rAlphaGamma$auc), as.numeric(rGammaGamma$auc), as.numeric(rBetaGamma$auc)), 3)
					, sep=''
				)
				legendText <- c(legendText,'p = 0.05')
				if( plotAIC ) legendText <- c(legendText,'AIC')
				
				legend('bottomright', legend=legendText
							 , col=c('red', 'navy', 'deepskyblue','grey70','grey70'), 
							 lty=c(3,1,3,NA,NA), 
							 pch=c(NA,NA,NA,16,17), lwd=4)
				
				if( !is.na(rAlphaBeta[[1]]) ) {
					plot.roc(rAlphaBeta, col='red', lwd=4, lty=3, add=FALSE, main="Degradation rate classification")
					alphaBetaIdx <- which.min(abs(0.05-rAlphaBeta[[4]]))
					points(rAlphaBeta[[3]][alphaBetaIdx],rAlphaBeta[[2]][alphaBetaIdx],pch=16,col='red',cex=2.5)
					if( plotAIC ) points(abCoordinatesAIC["specificity"],abCoordinatesAIC["sensitivity"],pch=17,col='red',cex=2.5)
				}
				if( !is.na(rGammaBeta[[1]]) ) {
					plot.roc(rGammaBeta, col='navy', lwd=4, lty=3, add=!is.na(rAlphaBeta[[1]]))
					gammaBetaIdx <- which.min(abs(0.05-rGammaBeta[[4]]))
					points(rGammaBeta[[3]][gammaBetaIdx],rGammaBeta[[2]][gammaBetaIdx],pch=16,col='navy',cex=2.5)
					if( plotAIC ) points(cbCoordinatesAIC["specificity"],cbCoordinatesAIC["sensitivity"],pch=17,col='navy',cex=2.5)
				}
				if( !is.na(rBetaBeta[[1]]) ) {
					plot.roc(rBetaBeta, col='deepskyblue', lwd=4, add=!(is.na(rAlphaBeta[[1]]) & is.na(rGammaBeta[[1]])))
					betaBetaIdx <- which.min(abs(0.05-rBetaBeta[[4]]))
					points(rBetaBeta[[3]][betaBetaIdx],rBetaBeta[[2]][betaBetaIdx],pch=16,col='deepskyblue',cex=2.5)
					if( plotAIC ) points(bbCoordinatesAIC["specificity"],bbCoordinatesAIC["sensitivity"],pch=17,col='deepskyblue',cex=2.5)
				}
				
				legendText <- paste(
					c('synthesis', 'processing', 'degradation')
					, ' - AUC='
					, signif(c(as.numeric(rAlphaBeta$auc), as.numeric(rGammaBeta$auc), as.numeric(rBetaBeta$auc)), 3)
					, sep=''
				)
				legendText <- c(legendText,'p = 0.05')
				if( plotAIC ) legendText <- c(legendText,'AIC')
				
				legend('bottomright', legend=legendText
							 , col=c('red', 'navy', 'deepskyblue','grey70','grey70'), 
							 lty=c(3,3,1,NA,NA), 
							 pch=c(NA,NA,NA,16,17), lwd=4)
				
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
		}
	})

coordinatesAIC <- function(AICclass,allResponses,class,predictor=NULL)
{
	if(is.null(predictor)) predictor <- class
	TP <- length(which(grepLogic(predictor,allResponses)&grepLogic(class,AICclass)))
	P <- length(grep(predictor,allResponses))
	
	TN <- length(which(!grepLogic(predictor,allResponses)&!grepLogic(class,AICclass)))
	N <- length(allResponses) - P
	
	specificityAIC <- TP/P
	sensitivityAIC <- TN/N
	
	c(specificity=specificityAIC,sensitivity=sensitivityAIC)
}

