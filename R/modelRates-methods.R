#' @rdname modelRates
#'
#' @description
#' This method model the synthesis, degradation and processing rates after their estimation by the constructor function
#' \code{\link{newINSPEcT}}. Estimated rates are not guaranteed to optimally describes provided input data yet. 
#' To this purpose, modeled rates can be generated and genes can be assigned to a transcriptional regulatory mechanism.
#' Modeled rates can be accessed via the method \code{\link{viewModelRates}} and gene classification according 
#' to the regulatory mechanism can be accessed by \code{\link{geneClass}}. The modeling procedure can be set by the
#' user changing the parameters via \code{\link{modelingParams}}
#' @param object An object of class INSPEcT
#' @param seed A numeric, indicatindg the seed to be set for reproducible results
#' @param BPPARAM Parallelization parameters for bplapply. By default bpparam()
#' @param verbose Either NULL or logical. If logical indicates whether to output some text
#' during computation or not, if NULL  it takes the information from the object
#' (see \code{\link{modelingParams}}) (Default: NULL)
#' @return An object of class INSPEcT with modeled rates
#' @seealso \code{\link{viewModelRates}}, \code{\link{geneClass}}, \code{\link{modelingParams}}
#' @details
#' When modeling many genes, parallelization is strongly suggested to reduce computational time.
#' Since all genes run independently, the computational time is diveded by the number of
#' cores used/available.
#' However, when modeling more than 500 genes, it may happen that 
#' a single gene returns an error that escapes the try/catch controls of INSPEcT. 
#' With the parallel mode, the error will propagate on all genes that have been computed 
#' with the same processor (or core). To avoid this, the computation could be splitted in chunks
#' and the whole data set can be obtaied by combining the chunks (see Examples).
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#' 	nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' 	## models removal
#' 	nascentInspObjThreeGenes <- removeModel(nascentInspObj10[1:3])
#' 	nascentInspObjThreeGenes <- modelRates(nascentInspObjThreeGenes, seed=1, BPPARAM=SerialParam())
#' 	## view modeled synthesis rates
#' 	viewModelRates(nascentInspObjThreeGenes, 'synthesis')
#' 	## view gene classes
#' 	geneClass(nascentInspObjThreeGenes)
#' }
setMethod('modelRates', 'INSPEcT', function(object
										  , seed=NULL
										  , BPPARAM=bpparam()
										  , verbose=NULL)
{
	if( length(object@model@ratesSpecs) > 0 )
		stop('Remove the model before running the model again. (See "?removeModel")')

	llConfidenceThreshold <- object@model@params$logLikelihoodConfidenceThreshold
	if(is.null(llConfidenceThreshold)) llConfidenceThreshold <- 0.95
	llConfidenceThreshold <- qchisq(llConfidenceThreshold,1)

	initialPenalityRelevance <- object@model@params$initialPenalityRelevance
	if(is.null(initialPenalityRelevance)) initialPenalityRelevance <- 1

	derivativePenalityRelevance <- object@model@params$derivativePenalityRelevance
	if(is.null(derivativePenalityRelevance)) derivativePenalityRelevance <- 10^-50

	NoNascent <- object@NoNascent

	if(NoNascent){message("No nascent RNA data mode.")}
	if(!NoNascent){message("Nascent RNA data mode.")}

	nInit <- object@params$nInit
	nIter <- object@params$nIter
	Dmax <- object@params$Dmax
	Dmin <- object@params$Dmin
	na.rm <- object@params$na.rm
	verbose <- object@params$verbose
	testOnSmooth <- object@params$testOnSmooth
	limitModelComplexity <- object@model@params$limitModelComplexity
	sigmoid <- object@params$useSigmoidFun
	object@params$seed <- seed

	computeDerivatives <- object@model@params$computeDerivatives
	if(is.null(computeDerivatives)) computeDerivatives <- TRUE

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
		if( !is.null(verbose) && !is.logical(verbose) )
			stop('verbose argument must be either NULL or logical.')
	
								# Time transformation
						#		a <- find_tt_par(tpts)
						#		tptsOriginal <- tpts
						#		tptsLinear <- time_transf(tptsOriginal,a)
						#		c <- abs(min(tptsLinear))
						#		tptsLinear <- tptsLinear + abs(min(tptsLinear))

						#		if(object@params$estimateRatesWith=="int")
						#		{
						#			ratesSpecs <- .inspect.engine_Integrative_NoNascent(tptsOriginal = tptsOriginal
						#															,tptsLinear = tptsLinear
						#															,a = a
						#															,c = c
						#															,concentrations = concentrations
						#															,rates = rates
						#															,BPPARAM = BPPARAM
						#															,na.rm = na.rm
						#															,verbose = verbose
						#															,testOnSmooth = testOnSmooth
						#															,seed = seed
						#															,nInit = nInit
						#															,nIter = nIter
						#															,limitModelComplexity=limitModelComplexity
						#															,sigmoid=sigmoid)
						#		}else{
						#			ratesSpecs <- .inspect.engine_Derivative_NoNascent(tptsOriginal = tptsOriginal
						#														  ,tptsLinear = tptsLinear
						#														  ,a = a
						#														  ,c = c
						#														  ,concentrations = concentrations
						#														  ,rates = rates
						#														  ,BPPARAM = BPPARAM
						#														  ,na.rm = na.rm
						#														  ,verbose = verbose
						#														  ,testOnSmooth = testOnSmooth
						#														  ,seed = seed
						#														  ,nInit = nInit
						#														  ,nIter = nIter
						#														  ,limitModelComplexity=limitModelComplexity)

		if(object@params$estimateRatesWith=="der")
		{
			ratesSpecs <- .inspect.engine_Derivative_NoNascent(tpts=tpts
															 , concentrations=concentrations
															 , rates=rates
															 , BPPARAM=BPPARAM
															 , na.rm=na.rm
															 , verbose=verbose
															 , testOnSmooth=testOnSmooth
															 , seed=seed
															 , nInit=nInit
															 , nIter=nIter
															 , limitModelComplexity=limitModelComplexity
															 , computeDerivatives=computeDerivatives
															 , useSigmoidFun=sigmoid
															 , initialPenalityRelevance = 1
															 , derivativePenalityRelevance = 10^-50
															 , llConfidenceThreshold = NULL)
		}else if(object@params$estimateRatesWith=="int"){
			ratesSpecs <- .inspect.engine_Integrative_NoNascent(tpts=tpts
															  , concentrations=concentrations
															  , rates=rates
															  , BPPARAM=BPPARAM
															  , na.rm=na.rm
															  , verbose=verbose
															  , testOnSmooth=testOnSmooth
															  , seed=seed
															  , nInit=nInit
															  , nIter=nIter
															  , limitModelComplexity=limitModelComplexity
															  , computeDerivatives=computeDerivatives
															  , useSigmoidFun=sigmoid
															  , initialPenalityRelevance = 1
															  , derivativePenalityRelevance = 10^-50
															  , llConfidenceThreshold = NULL)
		}else{stop('modelRates: the user must set the variable "estimateRatesWith" either equal to "int" or equal to "der" (default).')}

		## update and return the object
		# names(ratesSpecs) <- featureNames(object)
		object@model@ratesSpecs <- ratesSpecs
		object <- makeModelRates(object)
		return(object)

	}else{

		if( object@degDuringPulse ) stop('modelRates: degDuringPulse mode not implemented yet.')

		if( !is.null(seed) && !is.numeric(seed) )
			stop('Seed argument must be either NULL or numeric.')
		if( !is.null(verbose) && !is.logical(verbose) )
			stop('verbose argument must be either NULL or logical.')

	# Split between genes with and without intronic signal
		eiGenes <- rownames(concentrations$preMRNA[is.finite(concentrations$preMRNA[,1]),])
		eGenes <- rownames(concentrations$preMRNA[!is.finite(concentrations$preMRNA[,1]),])

		if(!is.null(eiGenes))
		{
			eiConcentrations <- lapply(concentrations,function(o){o[eiGenes,]})
			eiRates <- lapply(rates,function(o){o[eiGenes,]})

			if(object@params$estimateRatesWith=="der")
			{
				eiRatesSpecs <- .inspect.engine_Derivative_Nascent(tpts=tpts
															     , concentrations=eiConcentrations
															     , rates=eiRates
															     , BPPARAM=BPPARAM
															     , na.rm=na.rm
															     , verbose=verbose
															     , testOnSmooth=testOnSmooth
															     , seed=seed
															     , nInit=nInit
															     , nIter=nIter
															     , limitModelComplexity=limitModelComplexity
															     , computeDerivatives=computeDerivatives
															     , useSigmoidFun=sigmoid
															     , initialPenalityRelevance=initialPenalityRelevance
															     , derivativePenalityRelevance=derivativePenalityRelevance
															     , llConfidenceThreshold=llConfidenceThreshold)
			}else if(object@params$estimateRatesWith=="int"){
				eiRatesSpecs <- .inspect.engine_Integrative_Nascent(tpts=tpts
																  , concentrations=eiConcentrations
																  , rates=eiRates
																  , BPPARAM=BPPARAM
																  , na.rm=na.rm
																  , verbose=verbose
																  , testOnSmooth=testOnSmooth
																  , seed=seed
																  , nInit=nInit
																  , nIter=nIter
																  , limitModelComplexity=limitModelComplexity
																  , computeDerivatives=computeDerivatives
																  , useSigmoidFun=sigmoid
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
			message("Simple methods still to setup")
			# 			eConcentrations <- lapply(concentrations,function(o){o[eGenes,]})
			# 			eRates <- lapply(rates,function(o){o[eGenes,]})
			# 
			# 			eRatesSpecs <- .inspect.engine_Derivative_Nascent_Simple(tpts=tpts
			# 																   , concentrations=eConcentrations
			# 																   , rates=eRates
			# 																   , BPPARAM=BPPARAM
			# 																   , na.rm=na.rm
			# 																   , verbose=verbose
			# 																   , testOnSmooth=testOnSmooth
			# 																   , seed=seed
			# 																   , nInit=nInit
			# 																   , nIter=nIter
			# 																   , limitModelComplexity=limitModelComplexity
			# 																   , computeDerivatives=computeDerivatives
			# 																   , useSigmoidFun=sigmoid)
			# 			names(eRatesSpecs) <- eGenes
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
		object <- object[names(ratesSpecs)]
		object@model@ratesSpecs <- ratesSpecs
		object <- makeModelRates(object)
		object <- setConfidenceIntervals(object=object,confidenceIntervals=confidenceIntervals)

		return(object)
	}
})
