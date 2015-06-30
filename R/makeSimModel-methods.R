#' @rdname makeSimModel
#'
#' @description
#' This method allow the creation of synthesis, degradation and processing rates for a certain number of genes.
#' The rates are created according to the distributions of the real data-set which is given as an input of the
#' method. Different proportions of constant varying rates can be set and a new vector of time points can be
#' provided. This method has to be used before the \code{\link{makeSimDataset}} method.
#' @param object An object of class INSPEcT
#' @param nGenes A numeric with the number of synthtic genes to be created
#' @param newTpts A numeric verctor with time points of the synthtic dataset, if NULL the time points of the real dataset will be used
#' @param probs A numeric vector wich describes the probability of a rate to be constant, shaped like a sigmoid or like an impulse model
#' @param na.rm A logical that set whether missing values in the real dataset should be removed
#' @param seed A numeric to obtain reproducible results
#' @details The method \code{\link{makeSimModel}} generates an object of class INSPEcT_model that stores the parametric functions to genrate clean rates of a time-course. To any of the rates also a noise variance is associate but not used yet. In a typical workflow the output of \code{\link{makeSimModel}} is the input of the method \code{\link{makeSimDataset}}, that build the noisy rates and concentrations, given a specified number of replicates.
#' @return An object of class INSPEcT_model with synthetic rates
#' @seealso \code{\link{makeSimDataset}}
#' @examples
#' data('rpkms', package='INSPEcT')
#' tpts <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
#' tL <- 1/6
#' mycerIds <- newINSPEcT(tpts, tL, rpkms$rpkms_4su_exons, rpkms$rpkms_total_exons, 
#'	rpkms$rpkms_4su_introns, rpkms$rpkms_total_introns)
#' ## generate a synthtic data-set of 10 genes based on the real data-set
#' simRates <- makeSimModel(mycerIds, 10)
#' simData <- makeSimDataset(simRates, tpts, 1)
#' ## measure sensitivity/sensibility of synthesis, degradation and processing
#' ## rates identification
#' data('simRates', package='INSPEcT')
#' data('simData3rep', package='INSPEcT')
#' rocCurve(simRates, simData3rep)
#' ## measure classification with a different threshold for the chi-suared 
#' ## test acceptance of models
#' rocCurve(simRates, simData3rep, cTsh=.2)
#' ## generate a synthtic data-set of 10 genes based on the real data-set
#' ## with more replicates and more time points
#' newTpts <- c(0, 1/6, 1/3, 1/2, 1, 1.5, 2, 4, 8, 12, 16, 24)
#' simRates <- makeSimModel(mycerIds, 10, newTpts=newTpts)
#' simData <- makeSimDataset(simRates, newTpts, 3)
setMethod('makeSimModel', 'INSPEcT', 
	function(object, nGenes, newTpts=NULL, 
		probs=c(constant=.5,sigmoid=.3,impulse=.2), na.rm=TRUE, seed=NULL) {

	tpts <- object@tpts
	concentrations <- list(
		total=ratesFirstGuess(object, 'total')
		, total_var=ratesFirstGuessVar(object, 'total')
		, preMRNA=ratesFirstGuess(object, 'preMRNA')
		, preMRNA_var=ratesFirstGuessVar(object, 'preMRNA')
		)
	rates <- list(
		alpha=ratesFirstGuess(object, 'synthesis')
		, alpha_var=ratesFirstGuessVar(object, 'synthesis')
		, beta=ratesFirstGuess(object, 'degradation')
		, gamma=ratesFirstGuess(object, 'processing')
		)
	#
	out <- .makeSimData(nGenes, tpts, concentrations, rates
		, newTpts=newTpts, probs=probs, na.rm=na.rm, seed=seed)
	# arrange simdataSpecs form .makeSimData
	simdataSpecs <- out$simdataSpecs
	simdataSpecs <- lapply(simdataSpecs, function(x) list(x))
	#
	newObject <- new('INSPEcT_model')
	newObject@ratesSpecs <- simdataSpecs
	newObject@params$sim$flag <- TRUE
	newObject@params$sim$foldchange <- out$simulatedFC
	newObject@params$sim$noiseVar <- out$noiseVar

	return(newObject)

	})
