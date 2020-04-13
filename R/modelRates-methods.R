#' @rdname modelRates
#'
#' @description
#' This method model the synthesis, degradation and processing rates after their estimation by the constructor function
#' \code{\link{newINSPEcT}}. Estimated rates are not guaranteed to optimally describes provided input data yet. 
#' To this purpose, modeled rates can be generated and genes can be assigned to a transcriptional regulatory mechanism.
#' Modeled rates can be accessed via the method \code{\link{viewModelRates}} and gene classification according 
#' to the regulatory mechanism can be accessed by \code{\link{geneClass}}. The modeling options used for the 
#' modeling can be later accessed by the user via \code{\link{modelingParams}}. After modeling, model selection is run
#' by the method \code{\link{calculateRatePvals}} with default parameters.
#' @param object An object of class INSPEcT
#' @param estimateRatesWith Either "int" or "der". With "int" the degradation and processing
#'    rates are estimated integrating the system between one time point and the following. 
#'    With "der" degradation and processing rates are estimated using the derivative of total
#'    and pre mRNA. (default is "der")
#' @param useSigmoidFun A logical, whether to choose between sigmoid and impulse function 
#'    to fit rates and concentrations. In case not, always impulse function is used. 
#'    (default is TRUE)
#' @param nInit number of optimization to find the best functional representation of each rate (by default 10)
#' @param nIter number of max iteration during optimization (default is 300)
#' @param Dmin lower bondary for degradation rates in the NoNascent mode (default 1e-06)
#' @param Dmax upper bondary for degradation rates in the NoNascent mode (default 10)
#' @param seed A numeric, indicatindg the seed to be set for reproducible results. If NULL it is randomly selected (default NULL)
#' @param BPPARAM Parallelization parameters for bplapply. By default SerialParam()
#' @return An object of class INSPEcT with modeled rates
#' @seealso \code{\link{viewModelRates}}, \code{\link{calculateRatePvals}}, \code{\link{geneClass}}
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#' 	nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' 	## models removal
#' 	nascentInspObjThreeGenes <- removeModel(nascentInspObj10[1:3])
#' 	nascentInspObjThreeGenes <- modelRates(nascentInspObjThreeGenes, 
#' 	  seed=1, BPPARAM=SerialParam())
#' 	## view modeled synthesis rates
#' 	viewModelRates(nascentInspObjThreeGenes, 'synthesis')
#' 	## view gene classes
#' 	geneClass(nascentInspObjThreeGenes)
#' }
setMethod('modelRates', 'INSPEcT', function(object
																						, estimateRatesWith = c('der', 'int')
																						, useSigmoidFun = TRUE
																						, nInit = 10
																						, nIter = 300
																						, Dmin = 1e-06
																						, Dmax = 10
																						, seed = NULL
																						, BPPARAM = SerialParam())
{
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	if( length(object@model@ratesSpecs) > 0 )
		stop('Remove the model before running the model again. (See "?removeModel")')
	estimateRatesWith <- estimateRatesWith[1]
	if( !estimateRatesWith %in% c('int', 'der') )
		stop('modelRates: estimateRatesWith argument must be either "int" or "der"')
	if( !is.logical(useSigmoidFun) )
		stop('modelRates: useSigmoidFun argument must be a logical')		
	# if( !is.logical(testOnSmooth) )
	# 	stop('modelRates: testOnSmooth argument must be a logical')
	if( !is.numeric(nInit) )
		stop('modelRates: nInit argument must be a numeric')
	if( !is.numeric(nIter) )
		stop('modelRates: nIter argument must be a numeric')
	if( !is.numeric(Dmin) )
		stop('modelRates: Dmin argument must be a numeric')
	if( !is.numeric(Dmax) )
		stop('modelRates: Dmax argument must be a numeric')
	
	# if not set, select randomly a seed and store it in modelingParams
	if( is.null(seed) ) {
		seed <- sample(1:10000,1)
	}
	
	object@params <- list(
		estimateRatesWith = estimateRatesWith,
		useSigmoidFun = useSigmoidFun,
		# testOnSmooth = testOnSmooth,
		nInit = nInit,
		nIter = nIter,
		Dmin = Dmin,
		Dmax = Dmax,
		seed = seed
	)
	
	llConfidenceThreshold <- qchisq(0.95,1)
	initialPenalityRelevance <- 1
	derivativePenalityRelevance <- 10^-50

	NoNascent <- object@NoNascent

	if(NoNascent){message("No nascent RNA data mode.")}
	if(!NoNascent){message("Nascent RNA data mode.")}

	tpts <- tpts(object)

	concentrations <- list(total = ratesFirstGuess(object, 'total')
						 , total_var = ratesFirstGuessVar(object, 'total')
						 , preMRNA = ratesFirstGuess(object, 'preMRNA')
						 , preMRNA_var = ratesFirstGuessVar(object, 'preMRNA')
						 , mature = ratesFirstGuess(object, 'total') - ratesFirstGuess(object, 'preMRNA')
						 , mature_var = ratesFirstGuessVar(object, 'total') + ratesFirstGuessVar(object, 'preMRNA'))

	rates <- list(alpha=ratesFirstGuess(object, 'synthesis')
				, alpha_var=ratesFirstGuessVar(object, 'synthesis')
				, beta=ratesFirstGuess(object, 'degradation')
				, gamma=ratesFirstGuess(object, 'processing'))

	if(NoNascent){

		if( !is.numeric(tpts(object)) ) {
			stop("modelRates is not supported in steady-state analysis without nascent, run compareSteadyNoNascent instead.")
		}
	
		if( !is.null(seed) && !is.numeric(seed) )
			stop('Seed argument must be either NULL or numeric.')

		if(object@params$estimateRatesWith=="der")
		{
			ratesSpecs <- .inspect.engine_Derivative_NoNascent(tpts=tpts
															 , concentrations=concentrations
															 , rates=rates
															 , BPPARAM=BPPARAM
															 , na.rm=TRUE
															 , verbose=TRUE
															 # , testOnSmooth=testOnSmooth
															 , seed=seed
															 , nInit=nInit
															 , nIter=nIter
															 , computeDerivatives=TRUE
															 , useSigmoidFun=useSigmoidFun
															 , initialPenalityRelevance = initialPenalityRelevance
															 , derivativePenalityRelevance = derivativePenalityRelevance
															 , llConfidenceThreshold = NULL)
		}else if(object@params$estimateRatesWith=="int"){
			message("Integrative modeling")
			ratesSpecs <- .inspect.engine_Integrative_NoNascent(tpts=tpts
															  , concentrations=concentrations
															  , rates=rates
															  , BPPARAM=BPPARAM
															  , na.rm=TRUE
															  , verbose=TRUE
															  # , testOnSmooth=testOnSmooth
															  , seed=seed
															  , nInit=nInit
															  , nIter=nIter
															  , computeDerivatives=TRUE
															  , useSigmoidFun=useSigmoidFun
															  , initialPenalityRelevance = initialPenalityRelevance
															  , derivativePenalityRelevance = derivativePenalityRelevance
															  , llConfidenceThreshold = NULL)
		}else{stop('modelRates: the user must set the variable "estimateRatesWith" either equal to "int" or equal to "der" (default).')}
		object@NF <- FALSE
		object@model@ratesSpecs <- ratesSpecs
		object <- calculateRatePvals(object)
		# object <- makeModelRates(object) # called from calculateRatePvals
		return(object)

	}else{
		if( object@degDuringPulse ) stop('modelRates: degDuringPulse mode not implemented yet.')

		if( !is.null(seed) && !is.numeric(seed) )
			stop('Seed argument must be either NULL or numeric.')

	# Split between genes with and without intronic signal
		eiGenes <- rownames(concentrations$preMRNA)[is.finite(concentrations$preMRNA[,1])]
		eGenes <- NULL # rownames(concentrations$preMRNA[!is.finite(concentrations$preMRNA[,1]),])

		if(!is.null(eiGenes))
		{
			eiConcentrations <- lapply(concentrations,function(o){o[eiGenes,,drop=FALSE]})
			eiRates <- lapply(rates,function(o){o[eiGenes,,drop=FALSE]})

			if(object@params$estimateRatesWith=="der")
			{
				eiRatesSpecs <- .inspect.engine_Derivative_Nascent(tpts=tpts
															     , concentrations=eiConcentrations
															     , rates=eiRates
															     , BPPARAM=BPPARAM
															     , na.rm=TRUE
															     , verbose=TRUE
															     # , testOnSmooth=testOnSmooth
															     , seed=seed
															     , nInit=nInit
															     , nIter=nIter
															     , computeDerivatives=TRUE
															     , useSigmoidFun=useSigmoidFun
															     , initialPenalityRelevance=initialPenalityRelevance
															     , derivativePenalityRelevance=derivativePenalityRelevance
															     , llConfidenceThreshold=llConfidenceThreshold)
			}else if(object@params$estimateRatesWith=="int"){
				message("Integrative modeling")
				eiRatesSpecs <- .inspect.engine_Integrative_Nascent(tpts=tpts
																  , concentrations=eiConcentrations
																  , rates=eiRates
																  , BPPARAM=BPPARAM
																  , na.rm=TRUE
																  , verbose=TRUE
																  # , testOnSmooth=testOnSmooth
																  , seed=seed
																  , nInit=nInit
																  , nIter=nIter
																  , computeDerivatives=TRUE
																  , useSigmoidFun=useSigmoidFun
																  , initialPenalityRelevance=initialPenalityRelevance
																  , derivativePenalityRelevance=derivativePenalityRelevance
																  , llConfidenceThreshold=llConfidenceThreshold)
			}else{stop('modelRates: the user must set the variable "estimateRatesWith" either equal to "int" or equal to "der" (default).')}

			eiconfidenceIntervals <- eiRatesSpecs$confidenceIntervals
			eiRatesSpecs <- eiRatesSpecs$ratesSpecs

			eiGenes <- names(eiRatesSpecs)
		#	names(eiRatesSpecs) <- eiGenes
		}

		if(!is.null(eGenes))
		{
			message("Modeling for genes without introns to setup")
			# 			eConcentrations <- lapply(concentrations,function(o){o[eGenes,]})
			# 			eRates <- lapply(rates,function(o){o[eGenes,]})
			# 
			# 			eRatesSpecs <- .inspect.engine_Derivative_Nascent_Simple(tpts=tpts
			# 																   , concentrations=eConcentrations
			# 																   , rates=eRates
			# 																   , BPPARAM=BPPARAM
			# 																   , na.rm=TRUE
			# 																   , verbose=verbose
			# 																   , testOnSmooth=testOnSmooth
			# 																   , seed=seed
			# 																   , nInit=nInit
			# 																   , nIter=nIter
			# 																   , computeDerivatives=computeDerivatives
			# 																   , useSigmoidFun=useSigmoidFun)
			# 			names(eRatesSpecs) <- eGenes
			eRatesSpecs <- NULL
			econfidenceIntervals <- NULL
		}
		
		if(!is.null(eiGenes)&!is.null(eGenes)){
			ratesSpecs <- c(eiRatesSpecs,eRatesSpecs)
		}else if(!is.null(eiGenes)){
			ratesSpecs <- eiRatesSpecs
		}else if(!is.null(eGenes)){
			ratesSpecs <- eRatesSpecs
		}else{stop("No genes suitable for the analysis! ")}

		if(!is.null(eiGenes)&!is.null(eGenes)){
			confidenceIntervals <- c(eiconfidenceIntervals,econfidenceIntervals)
		}else if(!is.null(eiGenes)){
			confidenceIntervals <- eiconfidenceIntervals
		}else if(!is.null(eGenes)){
			confidenceIntervals <- econfidenceIntervals
		}else{stop("No genes suitable for the analysis! ")}
		
		## update and return the object
		object@NF <- FALSE
		object <- object[names(ratesSpecs)]
		object@model@ratesSpecs <- ratesSpecs
		object <- setConfidenceIntervals(object=object,confidenceIntervals=confidenceIntervals)
		object <- makeModelRates(object)
		object <- calculateRatePvals(object)
		return(object)
	}
})
