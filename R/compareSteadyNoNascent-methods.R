#' @rdname compareSteadyNoNascent
#'
#' @description
#' This function compare exons and introns expression level matrices, from two up to an arbitrary number of samples, in order
#' to identify genes which are oddly regluated, compared to an expected standard behaviour, from the post transcriptional point of view.
#' @param inspectIds An object of class INSPEcT.
#' @param expressionThreshold A parameter which sets how many log2 fold changes of distance from the median behaviour are imputable to noise.
#' @param log2FCThreshold A parameter which sets the log2 fold change distance from the median behaviour that is imputable to noise.
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
#'								   ,log2FCThreshold=2.)
#' head(regGenes)
#' table(regGenes)
setMethod('compareSteadyNoNascent', 'INSPEcT', function(inspectIds,
	expressionThreshold=0.25, log2FCThreshold=2.)
{
	# Mature, premature and total rpkms
	premature <- ratesFirstGuess(inspectIds, 'preMRNA')
	total <- ratesFirstGuess(inspectIds, 'total')
	mature <- total - premature

	prematureMedian <- apply(premature,1,function(r)median(r,na.rm=T))
	matureMedian <- apply(mature,1,function(r)median(r,na.rm=T))
	
	suppressWarnings(standardCurveFit <- standardCurveFitFunction(p=prematureMedian
																, m=matureMedian
																, err=log2FCThreshold))
	message(paste0("Trivial angle: ",standardCurveFit))

	premature[premature<=expressionThreshold] <- NA
	mature[mature<=expressionThreshold] <- NA

	suppressWarnings(classificationTmp <- classificationFunction(p=premature,m=mature,
		alpha=standardCurveFit,err=log2FCThreshold))

	return(classificationTmp)
})
