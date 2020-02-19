#' @rdname compareSteadyNoNascent
#'
#' @description
#' This function compare exons and introns expression level matrices, from two up to an arbitrary number of samples, in order
#' to identify genes which are oddly regluated, compared to an expected standard behaviour, from the post transcriptional point of view.
#' @param inspectIds An object of class INSPEcT.
#' @param expressionThreshold A parameter which sets how many log2 fold changes of distance from the median behaviour are imputable to noise.
#' @param log2FCThreshold A parameter which sets the log2 fold change distance from the median behaviour that is imputable to noise.
#' @param trivialAngle A numeric between 0 and 90 to define the standard behavior, if NULL (default) it is computed internally from the data.
#' @param returnNormScores A logical, if TRUE returned the deviations from the standard behavior normalized by the sd.
#' @param referenceCondition The label of the condition to use as reference, if NULL (default) the medians are used.
#' @examples
#' data('allcounts', package='INSPEcT')
#' data('featureWidths', package='INSPEcT')
#' data('libsizes', package='INSPEcT')
#' 
#' nascentCounts<-allcounts$nascent
#' matureCounts<-allcounts$mature
#' conditions<-letters[1:11]
#' expDes<-rep(conditions,3)
#' 
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'       allcounts=matureCounts
#'       ,libsize=totalLS
#'       ,exonsWidths=exWdths
#'       ,intronsWidths=intWdths
#'       ,experimentalDesign=expDes)
#'
#' matureInspObj <- newINSPEcT(tpts=conditions,matureExpressions=matExp_DESeq2)
#' 
#' regGenes<-compareSteadyNoNascent(inspectIds=matureInspObj
#' 								   ,expressionThreshold=0.25
#'								   ,log2FCThreshold=.5)
#' head(regGenes)
#' table(regGenes)
setMethod('compareSteadyNoNascent', 'INSPEcT', function(inspectIds,
	expressionThreshold=0.25, log2FCThreshold=2., trivialAngle=NULL, 
	returnNormScores=FALSE, referenceCondition=NULL, plot=FALSE)
{
	if( !.hasSlot(inspectIds, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	# Mature, premature and total rpkms
	premature <- ratesFirstGuess(inspectIds, 'preMRNA')
	total <- ratesFirstGuess(inspectIds, 'total')
	mature <- total - premature
	# Mature, premature and total variances
	prematureVar <- ratesFirstGuessVar(inspectIds, 'preMRNA')
	totalVar <- ratesFirstGuessVar(inspectIds, 'total')
	matureVar <- total + premature

	prematureMedian <- apply(premature,1,function(r)median(r,na.rm=T))
	matureMedian <- apply(mature,1,function(r)median(r,na.rm=T))
	
	if( is.null(trivialAngle) ) {
		suppressWarnings(standardCurveFit <- standardCurveFitFunction(p=prematureMedian
																	, m=matureMedian
																	, err=log2FCThreshold))
		message(paste0("Trivial angle: ",standardCurveFit))
	} else {
		standardCurveFit <- trivialAngle		
	}
	
	if( plot ) {
		## palette for 2D density plots (similar to heat colors)
		denscols <- c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100")
		denspalette <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))
		plot(log2(prematureMedian), log2(matureMedian), pch='.', xlab = 'log2 premature expression', ylab = 'log2 mature expression',
				 col = densCols(log2(prematureMedian), log2(matureMedian), colramp = denspalette))
		coef <- tan(standardCurveFit*pi/180)
		intercept <- median(log2(matureMedian), na.rm=TRUE) - coef * median(log2(prematureMedian), na.rm=TRUE)
		delta_intercept <- log2FCThreshold/cos(standardCurveFit*pi/180)
		abline(intercept, coef); abline(intercept - delta_intercept, coef, lty=2); abline(intercept + delta_intercept, coef, lty=2)
		abline(0,1,col=2)
		legend('topleft', legend = c(paste0("Trivial angle: ",standardCurveFit), 'log2 FC threshold', '', 'k2 = k3 limit'), text.col=c(1,1,NA,2), bty='n', lty=c(1,2,NA,1), col=c(1,1,NA,2))
	}

	perc <- function(x) round(length(which(x))/length(x)*100)
	perc_premature <- perc(premature<=expressionThreshold)
	if(perc_premature>0) message(paste0(perc_premature, '% of premature expressions below expressionThreshold set to NA'))
	premature[premature<=expressionThreshold] <- NA
	perc_mature <- perc(mature<=expressionThreshold)
	if(perc_mature>0) message(paste0(perc_mature, '% of mature expressions below expressionThreshold set to NA'))
	mature[mature<=expressionThreshold] <- NA

	if( is.null(referenceCondition) ) {
		suppressWarnings(log2maturemodel <- classificationFunction(p=premature,m=mature,
			alpha=standardCurveFit))
	} else {
		ref <- which(tpts(inspectIds)==referenceCondition)
		if( length(ref) != 1 ) stop('not existing referenceCondition')
		suppressWarnings(log2maturemodel <- classificationFunction(p=premature,m=mature,
			alpha=standardCurveFit, ref=ref))
	}

	if(returnNormScores) {
		maturemodel <- 2^log2maturemodel
		normScores <- (mature - maturemodel)/sqrt(matureVar)
		colnames(normScores) <- tpts(inspectIds)
		return(normScores)
	} else {
		scores <- log2(mature) - log2maturemodel
		colnames(scores) <- tpts(inspectIds)
		pi_angle <- standardCurveFit * pi/180
		threshold <- log2FCThreshold/cos(pi_angle)
		classification <- abs(scores)>threshold
		return(classification)
	}

	# return(list(scores=scores, norm=normScores, lognorm=lognormScores))
	# return(classificationTmp)
})