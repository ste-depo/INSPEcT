#' @rdname makeSimDataset
#'
#' @description
#' This method generates rates and concentrations where noise is added according to the desired number of 
#' replicates that the user set as an arguments from the INSPEcT_model object that has been created by the 
#' method of the class INSPEcT \code{\link{makeSimModel}}. Rates and concentrations can be generated at the 
#' time-points of interest. This method generates an INSPEcT object that can be modeled and the performance of
#' the modeling can be tested directly aginst the INSPEcT_model object created by \code{\link{makeSimModel}}.
#' @param object An object of class INSPEcT_model, usually the output of \code{\link{makeSimModel}}
#' @param tpts A numeric vector of time points where rates and concentrations have to be evaluated
#' @param nRep Number of replicates to simulate
#' @param seed A numeric to obtain reproducible results
#' @return An object of the class ExpressionSet containing rates and concentrations
#' @seealso \code{\link{makeSimModel}}
#' @examples
#' ## generate a synthtic data-set of 10 genes based on the real data-set
#' data('rpkms', package='INSPEcT')
#' tpts <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
#' tL <- 1/6
#' mycerIds <- newINSPEcT(tpts, tL, rpkms$foursu_exons, rpkms$total_exons, 
#' 	rpkms$foursu_introns, rpkms$total_introns)
#' simRates <- makeSimModel(mycerIds, 10)
#' simData <- makeSimDataset(simRates, tpts, 1)
#' ## load simulated datasets
#' data('simRates', package='INSPEcT')
#' data('simData3rep', package='INSPEcT')
#' ## measure sensitivity/sensibility of synthesis, degradation and processing
#' ## rates identification
#' dev.new()
#' rocCurve(simRates, simData3rep)
#' ## measure classification with a different threshold for the chi-squared 
#' ## test acceptance of models
#' rocCurve(simRates, simData3rep, cTsh=.2)
setMethod('makeSimDataset', 'INSPEcT_model', function(object, tpts, nRep, seed=NULL) {
	## create the clean concentrations and rates for each gene
	ratesSpecs <- object@ratesSpecs
	nGenes <- length(ratesSpecs)
	log_shift <- .find_tt_par(tpts)
	cleanRates <- lapply(1:nGenes, function(i) {
		tryCatch(
			.makeModel(tpts, ratesSpecs[[i]][[1]], log_shift, 
				.time_transf, deSolve::ode, .rxnrate)
			, error=function(e)
				.makeEmptyModel(tpts)
			)
		})
	## store total, preMRNA and alpha
	totalSim <- t(sapply(cleanRates, function(x) x$total))
	preMRNASim <- t(sapply(cleanRates, function(x) x$preMRNA))
	alphaSim <- t(sapply(cleanRates, function(x) x$alpha))
	## get noise variance form the object
	totalSim_noisevar <- object@params$sim$noiseVar$total
	preMRNASim_noisevar <- object@params$sim$noiseVar$pre
	alphaSim_noisevar <- object@params$sim$noiseVar$alpha
	## simulate the noise
	addNoise <- function(signal, noiseVar) {
		nConditions <- ncol(signal)
		noise <- t(sapply(sqrt(noiseVar), 
			function(sd) rnorm(nConditions,mean=0,sd=sd)))
		out <- signal + noise
		out[out < 0] = 0
		return(out)
	}
	if( !is.null(seed) ) set.seed(seed)
	totalSimReplicates <- do.call('cbind', 
		lapply(1:nRep, function(i) addNoise(totalSim,totalSim_noisevar)))
	rownames(totalSimReplicates) <- 1:nGenes
	preMRNASimReplicates <- do.call('cbind', 
		lapply(1:nRep, function(i) addNoise(preMRNASim,preMRNASim_noisevar)))
	rownames(preMRNASimReplicates) <- 1:nGenes
	alphaSimReplicates <- do.call('cbind', 
		lapply(1:nRep, function(i) addNoise(alphaSim,alphaSim_noisevar)))
	rownames(alphaSimReplicates) <- 1:nGenes
	tpts <- rep(tpts, nRep)
	## create the INSPEcT object
	newObject <- newINSPEcT(tpts=tpts, labeling_time=1, 
		, rpkms_4su_exons=alphaSimReplicates, rpkms_total_exons=totalSimReplicates
		, rpkms_total_introns=preMRNASimReplicates, simulatedData=TRUE)
	return(newObject)
	})