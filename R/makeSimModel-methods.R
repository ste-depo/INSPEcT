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
#' data('simRates', package='INSPEcT')
#' 
#' tpts <- simRates@params$tpts
#' 
#' simData2rep_4su <- makeSimDataset(object=simRates,tpts=tpts,nRep=3,No4sU=FALSE,seed=1)
#' simData2rep_4su <- modelRates(simData2rep_4su[1:10], seed=1)
setMethod('makeSimModel', 'INSPEcT', function(object
											, nGenes
											, probs=c(constant=.5
													, sigmoid=.3
													, impulse=.2)
											, na.rm=TRUE
											, seed=NULL)
{

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
	suppressWarnings(
		out <- .makeSimData(nGenes
						  , tpts
						  , concentrations
						  , rates
						  , probs = probs
						  , na.rm = na.rm
						  , seed = seed)
	)

	# arrange simdataSpecs form .makeSimData
	simdataSpecs <- out$simdataSpecs
	simdataSpecs <- lapply(simdataSpecs, function(x) list(x))
	#
	newObject <- new('INSPEcT_model')
	newObject@ratesSpecs <- simdataSpecs
	newObject@params$sim$flag <- TRUE
	newObject@params$sim$foldchange <- out$simulatedFC
	newObject@params$sim$noiseVar <- out$noiseVar
	newObject@params$sim$noiseFunctions <- out$noiseFunctions
	newObject@params$tpts <- tpts

	if(length(out$simdataSpecs)>nGenes){return(newObject[1:nGenes])} #Return only the required number of genes or less
	return(newObject)

	})
