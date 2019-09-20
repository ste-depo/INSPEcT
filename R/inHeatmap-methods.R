#' @rdname inHeatmap
#'
#' @description
#' A method to see as an heatmap the logRatios of synthesis, degradation and processing rates and 
#' pre-RNA and total RNA concentration of a population of genes, either at the level of 
#' etimated or modeled rates.
#' @param object An object of class INSPEcT
#' @param type Eiher "pre-model" or "model" to switch between pre-modeled or modeled features
#' @param breaks A vector of breaks for the heatmap
#' @param palette A color generating function, output of colorRampPalette
#' @param plot_matureRNA A logical. If set to TRUE, matrue-RNA is displayed instead of 
#'   total-RNA (default: FALSE)
#' @param absoluteExpression A logical. If set to FALSE, the plot representing the 
#'   intensity of expression is omitted. (default=TRUE)
#' @param show_rowLabels A character that represent the label names that will be
#'   shown on the y-axis of the heatmap. If NULL featureNames(object) will be shown
#'   (default is NULL)
#' @param clustering A logical. If set to FALSE, it displys genes the order they are, 
#'   with no clustering (default: TRUE)
#' @param clustIdx A numeric. Indicates which of the features are used for the 
#'   clustering. 0=absoluteExpression; 1=total-RNA/mature-RNA; 2=preMRNA; 
#'   3=synthesis; 4=degradation; 5=processing (default=3:5, meaning that
#'   synthesis, degradation and processing are used for the clustering)
#' @return A list of matrices containing the logRatios for total RNA levels, pre-RNA levels,
#' synthesis rates, degradation rates and processing rates. Matrices are ordered according to
#' the clustering.
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' inHeatmap(nascentInspObj10, 'pre-model')
#' inHeatmap(nascentInspObj10, 'model')
setMethod('inHeatmap', 'INSPEcT', function(object, type='pre-model'
	, breaks=seq(-1,1,length.out=51)
	, palette=colorRampPalette(c("green", "black", "firebrick3"))
	, plot_matureRNA=FALSE, absoluteExpression=TRUE
	, show_rowLabels=TRUE, clustering=TRUE, clustIdx=3:5) 
{
	if( !is.character(type) )
		stop('inHeatmap: type must be character')
	if( !type %in% c('model', 'pre-model') )
		stop('inHeatmap: type must either "model" or "pre-model"')
	# if( !is.null(rowLabels) ) {
		if( !is.logical(show_rowLabels) )
			stop('inHeatmap: rowLabels must be logical')
		# if( length(rowLabels) != nGenes(object) )
		# 	stop('inHeatmap: rowLabels must be a vector of length equal to nGenes(object)')
	# }
	if( !is.logical(plot_matureRNA) )
		stop('inHeatmap: plot_matureRNA must be logical')
	if( !is.logical(absoluteExpression) )
		stop('inHeatmap: absoluteExpression must be logical')
	if( !is.numeric(breaks) )
		stop('inHeatmap: breaks must be a numeric')
	if( class(palette) != 'function' )
		stop('inHeatmap: palette must be a (color generating) function')
	if( !is.logical(clustering) )
		stop('inHeatmap: clustering must be logical')
	if( !is.numeric(clustIdx) )
		stop('inHeatmap: clustIdx must be numeric')
	if( !all(clustIdx %in% 0:5) )
		stop('inHeatmap: clustIdx must contain only values between 0 and 5')

	oldMfrow <- par()$mfrow
	oldMar <- par()$mar
	nBreaks <- length(breaks)

	if( type=='pre-model') {
		total <- ratesFirstGuess(object, 'total')
		preMRNA <- ratesFirstGuess(object, 'preMRNA')
		alpha <- ratesFirstGuess(object, 'synthesis')
		beta <- ratesFirstGuess(object, 'degradation')
		gamma <- ratesFirstGuess(object, 'processing')
	} else {
		if( length(object@model@ratesSpecs) == 0 )
			stop('inHeatmap: modeled rates have not been computed yet. Use modelRates method')
		total <- viewModelRates(object, 'total')
		preMRNA <- viewModelRates(object, 'preMRNA')
		alpha <- viewModelRates(object, 'synthesis')
		beta <- viewModelRates(object, 'degradation')
		gamma <- viewModelRates(object, 'processing')
	}

	if( plot_matureRNA ) {
		total <- total-preMRNA
		total[total<0] <- 0
	}

	total[total<=0] <- NA
	preMRNA[preMRNA<=0] <- NA
	alpha[alpha<=0] <- NA
	beta[beta<=0] <- NA
	gamma[gamma<=0] <- NA

	total_l2fc <- log2(total[,-1,drop=FALSE]) - log2(total[,1])
	preMRNA_l2fc <- log2(preMRNA[,-1,drop=FALSE]) - log2(preMRNA[,1])
	alpha_l2fc <- log2(alpha[,-1,drop=FALSE]) - log2(alpha[,1])
	beta_l2fc <- log2(beta[,-1,drop=FALSE]) - log2(beta[,1])
	gamma_l2fc <- log2(gamma[,-1,drop=FALSE]) - log2(gamma[,1])

	## expression
	eLevel <- log2(total[,1,drop=FALSE])
	eRange <- quantile(eLevel, na.rm=TRUE, probs=c(.05,.95))
	eLevel[eLevel > max(eRange)] <- max(eRange)
	eLevel[eLevel < min(eRange)] <- min(eRange)

	clustMat <- do.call('cbind', list(
		sapply(1:ncol(total_l2fc), function(x) eLevel)
		, total_l2fc
		, preMRNA_l2fc
		, alpha_l2fc
		, beta_l2fc
		, gamma_l2fc
		)[clustIdx+1])

	## remove all NAs rows	
	# ix <- !apply(clustMat, 1, function(x) any(!is.finite(x)))
	ix <- !apply(clustMat, 1, function(x) all(!is.finite(x)))
	total <- total[ix, , drop=FALSE] ## for expression
	clustMat <- clustMat[ix, , drop=FALSE]
	total_l2fc <- total_l2fc[ix, , drop=FALSE]
	preMRNA_l2fc <- preMRNA_l2fc[ix, , drop=FALSE]
	alpha_l2fc <- alpha_l2fc[ix, , drop=FALSE]
	beta_l2fc <- beta_l2fc[ix, , drop=FALSE]
	gamma_l2fc <- gamma_l2fc[ix, , drop=FALSE]
	eLevel <- eLevel[ix, , drop=FALSE]

	if( nrow(clustMat)>1 & clustering ) {
		dd <- dist(clustMat)
		# dd <- as.dist((1 - cor(t(clustMat), use='complete.obs'))/2)
		geneOrder <- hclust(dd)$order		
	} else {
		geneOrder <- seq(1, nrow(clustMat))
	}

	## set the break boundaries into the data
	# min
	total_l2fc[total_l2fc<min(breaks)] <- min(breaks)
	preMRNA_l2fc[preMRNA_l2fc<min(breaks)] <- min(breaks)
	alpha_l2fc[alpha_l2fc<min(breaks)] <- min(breaks)
	beta_l2fc[beta_l2fc<min(breaks)] <- min(breaks)
	gamma_l2fc[gamma_l2fc<min(breaks)] <- min(breaks)
	# max
	total_l2fc[total_l2fc>max(breaks)] <- max(breaks)
	preMRNA_l2fc[preMRNA_l2fc>max(breaks)] <- max(breaks)
	alpha_l2fc[alpha_l2fc>max(breaks)] <- max(breaks)
	beta_l2fc[beta_l2fc>max(breaks)] <- max(breaks)
	gamma_l2fc[gamma_l2fc>max(breaks)] <- max(breaks)

	## the plotting function image will display it in reverse order
	## (bottom-up)
	geneOrder <- rev(geneOrder)

	if( absoluteExpression ) {
		## plot
		l <- matrix(rep(c(2,3,3,4,4,0,5,5,6,6,7,7), 10), nrow=10, byrow=TRUE)
		l <- rbind(c(rep(0,10),1,1), l)
		layout( l )
	} else {
		## plot
		l <- matrix(rep(c(2,2,3,3,0,4,4,5,5,6,6), 10), nrow=10, byrow=TRUE)
		l <- rbind(c(rep(0,9),1,1), l)
		layout( l )
	}

	## legend
	par(mar=c(2,1,.5,1))
	##
	cols <- palette(nBreaks-1)
	image(as.matrix(breaks), col=cols, xaxt='n', yaxt='n')
	axis(1, at=c(0,.5,1), labels=signif(c(min(breaks),min(breaks)/2+max(breaks)/2,max(breaks)),2))

	## data
	par(mar=c(4,1,2,1))

	if( absoluteExpression ) {
		## expression
		eBreaks <- seq(min(eRange),max(eRange),length.out=nBreaks)
		eCols <- colorRampPalette(c("black", "firebrick3"))(nBreaks-1)
		image(t(eLevel[geneOrder,]), col=eCols, breaks=eBreaks
			, xaxt='n', yaxt='n', main='exprs')
	}

	##
	if( plot_matureRNA ) {totalMain <- 'mature RNA'} else {totalMain <- 'total RNA'}

	rowLabels <- featureNames(object)
	rowLabels <- rowLabels[geneOrder]
	# colLabels <- paste('t', signif(object@tpts[-1],2), sep='_')
	colLabels <- signif(object@tpts[-1],2)
	atX <- seq(0, 1, length.out=length(colLabels))
	atY <- seq(0, 1, length.out=length(rowLabels))

	image(t(total_l2fc[geneOrder,,drop=FALSE]), col=cols, breaks=breaks
		, xaxt='n', yaxt='n', main=totalMain, xlab='time')
	axis(1, at=atX, labels=colLabels, las=3)
	image(t(preMRNA_l2fc[geneOrder,,drop=FALSE]), col=cols, breaks=breaks
		, xaxt='n', yaxt='n', main='pre-RNA', xlab='time')
	axis(1, at=atX, labels=colLabels, las=3)
	image(t(alpha_l2fc[geneOrder,,drop=FALSE]), col=cols, breaks=breaks
		, xaxt='n', yaxt='n', main='synthesis', xlab='time')
	axis(1, at=atX, labels=colLabels, las=3)
	if( show_rowLabels ) axis(2, at=atY, labels=rowLabels, las=1, tick=FALSE, line=1)
	image(t(gamma_l2fc[geneOrder,,drop=FALSE]), col=cols, breaks=breaks
		, xaxt='n', yaxt='n', main='processing', xlab='time')
	axis(1, at=atX, labels=colLabels, las=3)
	image(t(beta_l2fc[geneOrder,,drop=FALSE]), col=cols, breaks=breaks
		, xaxt='n', yaxt='n', main='degradation', xlab='time')
	axis(1, at=atX, labels=colLabels, las=3)

	if( plot_matureRNA  ) {
		out <- list(mature_l2fc=total_l2fc[rev(geneOrder),,drop=FALSE]
			, preMRNA_l2fc=preMRNA_l2fc[rev(geneOrder),,drop=FALSE]
			, synthesis_l2fc=alpha_l2fc[rev(geneOrder),,drop=FALSE]
			, processing_l2fc=gamma_l2fc[rev(geneOrder),,drop=FALSE]
			, degradation_l2fc=beta_l2fc[rev(geneOrder),,drop=FALSE]
		)
	} else {
		out <- list(total_l2fc=total_l2fc[rev(geneOrder),,drop=FALSE]
			, preMRNA_l2fc=preMRNA_l2fc[rev(geneOrder),,drop=FALSE]
			, synthesis_l2fc=alpha_l2fc[rev(geneOrder),,drop=FALSE]
			, processing_l2fc=gamma_l2fc[rev(geneOrder),,drop=FALSE]
			, degradation_l2fc=beta_l2fc[rev(geneOrder),,drop=FALSE]
		)
	}
	par(mfrow=oldMfrow, mar=oldMar)
	out <- out

	})
