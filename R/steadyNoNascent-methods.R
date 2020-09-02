# FIX PM gene plot to include the confidence interval in the plot
# FIX class slot for refernece condition to include the NaN / median option

#' @rdname premature
#' 
#' @description
#' Extract premature RNA expressions
#' @param object An object of class INSPEcT_steadyNoNascent
#' @return A matrix containing premature RNA expressions
setMethod('premature', 'INSPEcT_steadyNoNascent', function(object)
{
	return(object@premature)
})
#' @rdname mature
#' 
#' @description
#' Extract mature RNA expressions
#' @param object An object of class INSPEcT_steadyNoNascent
#' @return A matrix containing mature RNA expressions
setMethod('mature', 'INSPEcT_steadyNoNascent', function(object)
{
	return(object@mature)
})

#' @rdname PTratio
#' 
#' @description
#' Extract the ratio between mature and premature RNAs
#' @param object An object of class INSPEcT_steadyNoNascent
#' @param infToNA A logical indicating whether infinite values (originating from zero 
#' valued premature expressions) should be set artificially to NA or not
#' @return A matrix containing the PTratios
setMethod('PTratio', 'INSPEcT_steadyNoNascent', function(object, infToNA=TRUE)
{
	premature <- object@premature
	mature <- object@mature
	ptratio <- mature/premature
	if( infToNA ) ptratio[!is.finite(ptratio)] <- NA
	return(ptratio)
})

#' @rdname prematureVar
#' 
#' @description
#' Extract premature RNA expressions variances
#' @param object An object of class INSPEcT_steadyNoNascent
#' @return A matrix containing premature RNA expressions variances
setMethod('prematureVar', 'INSPEcT_steadyNoNascent', function(object)
{
	return(object@prematureVar)
})
#' @rdname matureVar
#' 
#' @description
#' Extract mature RNA expressions variances
#' @param object An object of class INSPEcT_steadyNoNascent
#' @return A matrix containing mature RNA expressions variances
setMethod('matureVar', 'INSPEcT_steadyNoNascent', function(object)
{
	return(object@matureVar)
})

#' @rdname compareSteadyNoNascent
#'
#' @description
#' This function compare exons and introns expression level matrices, from two up to an arbitrary number of samples, in order
#' to identify genes which are oddly regluated, compared to an expected standard behaviour, from the post transcriptional point of view.
#' @param inspectIds An object of class INSPEcT_steadyNoNascent
#' @param expressionThreshold A parameter which sets how many log2 fold changes of distance from the median behaviour are imputable to noise.
#' @param log2FCThreshold A parameter which sets the log2 fold change distance from the median behaviour that is imputable to noise.
#' @param trivialAngle A numeric between 0 and 90 to define the standard behavior, if NaN (default) it is computed internally from the data.
#' @param returnNormScores A logical, if TRUE returned the deviations from the standard behavior normalized by the sd.
#' @param referenceCondition The label of the condition to use as reference, if NaN (default) the medians are used.
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
#' matureInspObj<-compareSteadyNoNascent(inspectIds=matureInspObj
#' 								 ,expressionThreshold=0.25
#'								 ,log2FCThreshold=.5)
#' regGenes <- PTreg(matureInspObj)
#' head(regGenes)
#' table(regGenes)
setMethod('compareSteadyNoNascent', 'INSPEcT_steadyNoNascent', function(inspectIds,
																																	 expressionThreshold=0.25, log2FCThreshold=2., trivialAngle=NaN, 
																																	 returnNormScores=FALSE, referenceCondition='median')
{
  checkINSPEcTObjectversion(inspectIds)
	# Mature, premature and total rpkms
	premature <- premature(inspectIds)
	mature <- mature(inspectIds)
	# Mature, premature and total variances
	prematureVar <- prematureVar(inspectIds)
	matureVar <- matureVar(inspectIds)

	premature[premature==0] <- NaN
	mature[mature==0] <- NaN
	
	prematureMedian <- apply(premature,1,function(r)median(r,na.rm=T))
	matureMedian <- apply(mature,1,function(r)median(r,na.rm=T))
	
	if( is.na(trivialAngle) ) {
		suppressWarnings(standardCurveFit <- standardCurveFitFunction(p=prematureMedian
																																	, m=matureMedian
																																	, err=log2FCThreshold))
		message(paste0("Trivial angle: ",standardCurveFit))
	} else {
		standardCurveFit <- trivialAngle		
	}
	
	perc <- function(x) round(length(which(x))/length(x)*100)
	perc_premature <- perc(premature<=expressionThreshold)
	if(perc_premature>0) message(paste0(perc_premature, '% of premature expressions below expressionThreshold set to NA'))
	premature[premature<=expressionThreshold] <- NA
	perc_mature <- perc(mature<=expressionThreshold)
	if(perc_mature>0) message(paste0(perc_mature, '% of mature expressions below expressionThreshold set to NA'))
	mature[mature<=expressionThreshold] <- NA
	
	if( referenceCondition == 'median' ) {
		suppressWarnings(log2maturemodel <- classificationFunction(p=premature,m=mature,
																															 alpha=standardCurveFit))
	} else {
		ref <- which(inspectIds@sampleNames==referenceCondition)
		if( length(ref) != 1 ) stop('not existing referenceCondition')
		suppressWarnings(log2maturemodel <- classificationFunction(p=premature,m=mature,alpha=standardCurveFit, ref=ref))
	}
	
	# update the objects
	inspectIds@trivialAngle <- standardCurveFit
	inspectIds@log2FCThreshold <- log2FCThreshold
	inspectIds@expressionThreshold <- expressionThreshold
	inspectIds@referenceCondition <- referenceCondition
	
	# identify pt regulated genes
	if(returnNormScores) {
		maturemodel <- 2^log2maturemodel
		normScores <- (mature - maturemodel)/sqrt(matureVar)
		colnames(normScores) <- inspectIds@sampleNames
		inspectIds@ptreg <- normScores
	} else {
		scores <- log2(mature) - log2maturemodel
		colnames(scores) <- inspectIds@sampleNames
		pi_angle <- standardCurveFit * pi/180
		threshold <- log2FCThreshold/cos(pi_angle)
		classification <- abs(scores)>abs(threshold)
		inspectIds@ptreg <- classification
	}
	return(inspectIds)
})

#' @rdname PTreg
#' 
#' @description
#' Extract the post-transcriptional regulation matrix
#' @param object An object of class INSPEcT_steadyNoNascent
#' @return A matrix containing the post-transcriptional regulated genes.
#' This matrix is generated by the method compareSteadyStateNoNascent. 
#' It generally report 1 for regulated genes in specific samples, 0 for 
#' non regulated genes and NA for genes that do not passed the 
#' expression threshold. In case the argument returnNormScores was set 
#' to TRUE, instead of discretes values, the deviations from the 
#' expected model normalized by the experimental standard deviation
#' is reported.
setMethod('PTreg', 'INSPEcT_steadyNoNascent', function(object)
{
	return(object@ptreg)
})

#' @rdname plotPMtrend
#' 
#' @description
#' Plot the null model estimated for the specific dataset. The null model
#' is the trend between premature and mature expression, which is usually
#' linear in the log-log scale and generally points to an increase in the
#' ratio between premature and mature RNA at increased levels of 
#' expression
#' @param inspectIds An object of class INSPEcT_steadyNoNascent
setMethod('plotPMtrend', 'INSPEcT_steadyNoNascent', function(inspectIds) 
{
	premature <- premature(inspectIds)
	mature <- mature(inspectIds)
	##
	prematureMedian <- apply(premature,1,function(r)median(r,na.rm=T))
	matureMedian <- apply(mature,1,function(r)median(r,na.rm=T))
	##
	log2FCThreshold <- inspectIds@log2FCThreshold
	if( length(log2FCThreshold) == 0 ) log2FCThreshold <- 2
	standardCurveFit <- inspectIds@trivialAngle
	if( length(standardCurveFit) == 0 )
		suppressWarnings(standardCurveFit <- standardCurveFitFunction(p=prematureMedian
																																	, m=matureMedian
																																	, err=log2FCThreshold))
	## palette for 2D density plots (similar to heat colors)
	denscols <- c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100")
	denspalette <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))
	##
	if( length(prematureMedian)<1000 ) pch = 1 else pch='.'
	plot(log2(prematureMedian), log2(matureMedian), pch=pch, xlab = 'log2 premature expression', ylab = 'log2 mature expression',
			 col = densCols(log2(prematureMedian), log2(matureMedian), colramp = denspalette))
	##
	coef <- tan(standardCurveFit*pi/180)
	intercept <- median(log2(matureMedian), na.rm=TRUE) - coef * median(log2(prematureMedian), na.rm=TRUE)
	delta_intercept <- log2FCThreshold/cos(standardCurveFit*pi/180)
	abline(intercept, coef)
	abline(intercept - delta_intercept, coef, lty=2)
	abline(intercept + delta_intercept, coef, lty=2)
	## 
	abline(0,1,col=2)
	legend('topleft', legend = c(paste0("Trivial angle: ",standardCurveFit), 'log2 FC threshold'), text.col=c(1,1), bty='n', lty=c(1,2), col=c(1,1))
	legend('bottomright', legend = c('k2 = k3 limit'), text.col=c(2), bty='n', lty=c(1), col=c(2))
})

#' @rdname plotPMgene
#' 
#' @description
#' Plot the premature and mature expressions of a specific gene in the 
#' different samples of the dataset along with the null model and the 
#' log2 fold change threshold. Individal observations that fall outside of the
#' dashed lines are considered post-transcriptional events.
#' @param object An object of class INSPEcT_steadyNoNascent
#' @param gene_id A numeric that indicated the index of the gene to be plotted
#' @param samples_colors The color code relative to the samples
setMethod('plotPMgene', 'INSPEcT_steadyNoNascent', function(object, gene_id, samples_colors=1)
{
	p <- log2(premature(object)[gene_id,])
	m <- log2(mature(object)[gene_id,])
	##
	referenceCondition <- object@referenceCondition
	standardCurveFit <- object@trivialAngle
	log2FCThreshold <- object@log2FCThreshold
	##
	if( referenceCondition == 'median' ) {
		m_reference <- median(m, na.rm=TRUE)
		p_reference <- median(p, na.rm=TRUE)
	} else {
		m_reference <- m[referenceCondition]
		p_reference <- p[referenceCondition]
	}
	coef <- tan(standardCurveFit*pi/180)
	intercept <- m_reference - coef * p_reference
	delta_intercept <- log2FCThreshold/cos(standardCurveFit*pi/180)
	##
	plot(p, m, xlab='log2 premature expression', ylab='log2 mauture expression', col=samples_colors,
			 main=paste0('gene: ', object@geneNames[gene_id]), pch=19)
	abline(intercept, coef)
	abline(intercept + delta_intercept, coef, lty=2)
	abline(intercept - delta_intercept, coef, lty=2)
})

#' @rdname INSPEcT_steadyNoNascent-class
#' @param x An object of class INSPEcT_steadyNoNascent
#' @param i A numeric, a vector of logicals indicating the rows to be extracted
#' @param j A numeric, a vector of logicals indicating the columns to be extracted
setMethod('[', 'INSPEcT_steadyNoNascent', function(x, i, j) {
	if( !missing(i) ) {
		if( is.logical(i) ) i <- which(i)
		# x@sampleNames <- x@sampleNames[]
		x@geneNames <- x@geneNames[i]
		x@premature <- x@premature[i,]
		x@mature <- x@mature[i,]
		x@prematureVar <- x@prematureVar[i,]
		x@matureVar <- x@matureVar[i,]
		# x@trivialAngle <- x@trivialAngle[]
		# x@log2FCThreshold <- x@log2FCThreshold[]
		# x@expressionThreshold <- x@expressionThreshold[]
		# x@referenceCondition <- x@referenceCondition[]
		if( nrow(x@ptreg)>0 ) x@ptreg <- x@ptreg[i,]
	}
	if( !missing(j) ) {
		if( is.logical(j) ) j <- which(j)
		x@sampleNames <- x@sampleNames[j]
		# x@geneNames <- x@geneNames[]
		x@premature <- x@premature[,j]
		x@mature <- x@mature[,j]
		x@prematureVar <- x@prematureVar[,j]
		x@matureVar <- x@matureVar[,j]
		# x@trivialAngle <- x@trivialAngle[]
		# x@log2FCThreshold <- x@log2FCThreshold[]
		# x@expressionThreshold <- x@expressionThreshold[]
		if( length(x@referenceCondition)>0 ) {
			if( x@referenceCondition %in% j ) {
				x@referenceCondition <- which( j == x@referenceCondition )
			} else {
				warning('referenceCondition is not among selected columns')
				x@referenceCondition <- numeric(0)
			}
		}
		if( ncol(x@ptreg)>0 ) x@ptreg <- x@ptreg[,j]
	}
	return(x)
})

#' @rdname INSPEcT_steadyNoNascent-class
#' @param object An object of class INSPEcT_steadyNoNascent
#' @return Method show for objects of class INSPEcT_steadyNoNascent
setMethod('show', 'INSPEcT_steadyNoNascent', function(object) {
	str(object)
})
