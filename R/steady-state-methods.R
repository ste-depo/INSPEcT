#' @rdname compareSteady
#'
#' @description
#' This method compares two object of class INSPEcT in order to identify differential usage
#' of synthesis, processing or degradation rates in two different steady-state conditions. 
#' The two INSPEcT objects must have been profiled with replicates in order to provide
#' a statistical significance to the differences between their rates.
#' @param inspectIds1 An object of calss INSPEcT
#' @param inspectIds2 A second object of calss INSPEcT
#' @return An object of class INSPEcT_diffsteady which contains both the absolute 
#' quantification of the rates as well as the comparison with the statistcal significance
#' associated for each gene and rate. (See \code{\link{INSPEcT_diffsteady-class}})
#' @examples
#' ## load data
#' data('simRates', package='INSPEcT')
#' data('simData3rep', package='INSPEcT')
#' ## generate a new data set with 3 replicate to make the comparison
#' newTpts <- c(0, 1/6)
#' simData3rep_2 <- makeSimDataset(simRates, newTpts, 3, seed=2)
#' ## compare
#' diffrates <- compareSteady(simData3rep, simData3rep_2)
#' ## visualize results
#' diffrates
setMethod('compareSteady', signature('INSPEcT','INSPEcT'), 
	function(inspectIds1, inspectIds2) 
{

	if( all(is.na(fData(inspectIds1@ratesFirstGuess)$nRepT0)) |
		all(fData(inspectIds1@ratesFirstGuess)$nRepT0<2) )
		stop('Object "inspectIds1" has less than two replicates at steady state.')

	if( all(is.na(fData(inspectIds2@ratesFirstGuess)$nRepT0)) |
		all(fData(inspectIds2@ratesFirstGuess)$nRepT0<2) )
		stop('Object "inspectIds2" has less than two replicates at steady state.')

	cGenes <- intersect(featureNames(inspectIds1), featureNames(inspectIds2))
	if( length(cGenes)==0 ) {
		stop('compareSteady: The two datasets do not share any gene.')
	} else {
		inspectIds1 <- inspectIds1[cGenes]
		inspectIds2 <- inspectIds2[cGenes]
	}
	
	## synthesis
	s1var <- fData(inspectIds1@ratesFirstGuess)$synthesis_t0
	s1mean <- ratesFirstGuess(inspectIds1, 'synthesis')[,1]
	s2var <- fData(inspectIds2@ratesFirstGuess)$synthesis_t0
	s2mean <- ratesFirstGuess(inspectIds2, 'synthesis')[,1]
	## degradation
	d1var <- fData(inspectIds1@ratesFirstGuess)$degradation_t0
	d1mean <- ratesFirstGuess(inspectIds1, 'degradation')[,1]
	d2var <- fData(inspectIds2@ratesFirstGuess)$degradation_t0
	d2mean <- ratesFirstGuess(inspectIds2, 'degradation')[,1]
	## processing
	p1var <- fData(inspectIds1@ratesFirstGuess)$processing_t0
	p1mean <- ratesFirstGuess(inspectIds1, 'processing')[,1]
	p2var <- fData(inspectIds2@ratesFirstGuess)$processing_t0
	p2mean <- ratesFirstGuess(inspectIds2, 'processing')[,1]

	# transform into log scale
	## synthesis
	s1logvar <- sapply(1:length(s1var), function(i) (1/(log(2)*s1mean[i]))^2*s1var[i])
	s1logmean <- log2(s1mean)
	s2logvar <- sapply(1:length(s2var), function(i) (1/(log(2)*s2mean[i]))^2*s2var[i])
	s2logmean <- log2(s2mean)
	## degradation
	d1logvar <- sapply(1:length(d1var), function(i) (1/(log(2)*d1mean[i]))^2*d1var[i])
	d1logmean <- log2(d1mean)
	d2logvar <- sapply(1:length(d2var), function(i) (1/(log(2)*d2mean[i]))^2*d2var[i])
	d2logmean <- log2(d2mean)
	## processing
	p1logvar <- sapply(1:length(p1var), function(i) (1/(log(2)*p1mean[i]))^2*p1var[i])
	p1logmean <- log2(p1mean)
	p2logvar <- sapply(1:length(p2var), function(i) (1/(log(2)*p2mean[i]))^2*p2var[i])
	p2logmean <- log2(p2mean)

	## sample size
	s1n <- d1n <- p1n <- fData(inspectIds1@ratesFirstGuess)$nRepT0[1]
	s2n <- d2n <- p2n <- fData(inspectIds2@ratesFirstGuess)$nRepT0[1]

	rate.t.test <- function(m1, m2, s1, s2, n1, n2) {
		sx1x2 <- sqrt( ( (n1-1)*s1^2 + (n2-1)*s2^2 ) / (n1+n2-2) )
		t <- (m1 - m2) / ( sx1x2 * sqrt(1/n1+1/n2) )
		df <- n1+n2-2
		## two tailed test
		pt(-abs(t), df)*2
	}

	s.test <- sapply(1:length(s1mean), function(i) 
		rate.t.test(s1logmean[i], s2logmean[i], sqrt(s1logvar[i]), sqrt(s2logvar[i]), s1n, s2n)
		)
	d.test <- sapply(1:length(d1mean), function(i) 
		rate.t.test(d1logmean[i], d2logmean[i], sqrt(d1logvar[i]), sqrt(d2logvar[i]), d1n, d2n)
		)
	p.test <- sapply(1:length(p1mean), function(i) 
		rate.t.test(p1logmean[i], p2logmean[i], sqrt(p1logvar[i]), sqrt(p2logvar[i]), p1n, p2n)
		)

	synthesis_res <- data.frame(
		condition1=s1logmean,
		variance1=s1logvar,
		condition2=s2logmean,
		variance2=s2logvar,
		samplesize1=s1n,
		samplesize2=s2n,
		log2mean=log2(sqrt(s1mean*s2mean)),
		log2fc=log2(s2mean/s1mean),
		pval=s.test,
		padj=p.adjust(s.test, method='BH')
		)

	degradation_res <- data.frame(
		condition1=d1mean,
		variance1=d1var,
		condition2=d2mean,
		variance2=d2var,
		samplesize1=d1n,
		samplesize2=d2n,
		log2mean=log2(sqrt(d1mean*d2mean)),
		log2fc=log2(d2mean/d1mean),
		pval=d.test,
		padj=p.adjust(d.test, method='BH')
		)

	processing_res <- data.frame(
		condition1=p1mean,
		variance1=p1var,
		condition2=p2mean,
		variance2=p2var,
		samplesize1=p1n,
		samplesize2=p2n,
		log2mean=log2(sqrt(p1mean*p2mean)),
		log2fc=log2(p2mean/p1mean),
		pval=p.test,
		padj=p.adjust(p.test, method='BH')
		)

	new_object <- new('INSPEcT_diffsteady')
	new_object@synthesis <- synthesis_res
	new_object@degradation <- degradation_res
	new_object@processing <- processing_res

	return(new_object)

})

#' @rdname INSPEcT_diffsteady-class
#' @examples
#' data('simData3rep', package='INSPEcT')
#' data('simData3rep_2', package='INSPEcT')
#' diffrates <- compareSteady(simData3rep, simData3rep_2)
#' head(synthesis(diffrates))
setMethod('synthesis', 'INSPEcT_diffsteady', function(object) slot(object, 'synthesis'))
#' @rdname INSPEcT_diffsteady-class
#' @examples
#' head(processing(diffrates))
setMethod('processing', 'INSPEcT_diffsteady', function(object) slot(object, 'processing'))
#' @rdname INSPEcT_diffsteady-class
#' @examples
#' head(degradation(diffrates))
setMethod('degradation', 'INSPEcT_diffsteady', function(object) slot(object, 'degradation'))

#' @name plotMA
#' @title MA-plot from base means and log fold changes
NULL

#' @rdname plotMA
#' @description Visualize the comparison between the rates calculated from two different INSPEcT objects
#' profiled in steady-state conditions.
#' @param object An object of calss INSPEcT_diffsteady
#' @param ... Additional parameters, see Details section
#' @details
#' Possible arguments to "plotMA":
#' \itemize{
#' \item "rate" - A character, which represent the rate to be visualized, either "synthesis", "processing" or "degradation". By default, "synthesis" is chosen.
#' \item "alpha" - A numeric, The confidence interval for significance (FDR), by default 0.1
#' \item "xlim" - A numeric vector of length 2, limits of x-axis, by default the range of the data.
#' \item "xlab" - A character, the label of x-axis, by default "log2 geometric mean"
#' \item "ylim" - A numeric vector of length 2, limits of y-axis, by default the range of the data.
#' \item "ylab" - A character, the label of y-axis, by default "log2 fold change"
#' \item "main" - A character, the title of the plot, by default the name of the visualized rate.
#' }
#' @seealso \url{http://en.wikipedia.org/wiki/MA_plot}
#' @examples
#' data('simData3rep', package='INSPEcT')
#' data('simData3rep_2', package='INSPEcT')
#' diffrates <- compareSteady(simData3rep, simData3rep_2)
#' plotMA(diffrates, alpha=.5)
setMethod('plotMA', 'INSPEcT_diffsteady', function(object, ...) {
	addargs <- list(...)

	## argument "rate"
	feasiblerates <- c('synthesis','processing','degradation')
	if( any(names(addargs) == 'rate') ) {
		rate <- addargs[['rate']]
		if( ! rate %in% feasiblerates )
			stop('plotMA: Unrecognized "rate" argument.')
	} else rate <- "synthesis"
	
	## argument "alpha"
	if( any(names(addargs) == 'alpha') ) {
		alpha <- addargs[['alpha']]
		if( !is.numeric(alpha) | alpha<0 | alpha>1 )
			stop('plotMA: "alpha" must be numeric and between 0 and 1.')
	} else alpha <- .1

	data <- slot(object, rate)
	x <- data$log2mean
	y <- data$log2fc
	signif_genes <- data$padj<alpha
	ix <- is.na(x) | is.na(y)
	if( any(ix) ) {
		x <- x[!ix]
		y <- y[!ix]
		signif_genes <- signif_genes[!ix]
	}

	## argument "xlim"
	if( any(names(addargs) == 'xlim') ) {
		xlim <- addargs[['xlim']]
		if( !(is.numeric(xlim) & length(xlim)==2) )
			stop('plotMA: "xlim" must be a numeric of length 2.')
		x[x<xlim[1]] <- xlim[1]
		x[x>xlim[2]] <- xlim[2]
	} else xlim <- range(x, na.rm=TRUE)

	## argument "ylim"
	if( any(names(addargs) == 'ylim') ) {
		ylim <- addargs[['ylim']]
		if( !(is.numeric(ylim) & length(ylim)==2) )
			stop('plotMA: "ylim" must be a numeric of length 2.')
		y[y<ylim[1]] <- ylim[1]
		y[y>ylim[2]] <- ylim[2]
	} else ylim <- range(y, na.rm=TRUE)

	## argument "xlab"
	if( any(names(addargs) == 'xlab') ) {
		xlab <- addargs[['xlab']]
		if( !is.character(xlab) )
			stop('plotMA: "xlab" must be a character.')	 		
	} else xlab <- 'log2 geometric mean'

	## argument "ylab"
	if( any(names(addargs) == 'ylab') ) {
		ylab <- addargs[['ylab']]
		if( !is.character(ylab) )
			stop('plotMA: "ylab" must be a character.')	 		
	} else ylab <- 'log2 fold change'

	## argument "main"
	if( any(names(addargs) == 'main') ) {
		main <- addargs[['main']]
		if( !is.character(main) )
			stop('plotMA: "main" must be a character.')	 		
	} else main <- rate

	smoothScatter(x, y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
	points(x[signif_genes], y[signif_genes], col='orange', pch=2)

	})

