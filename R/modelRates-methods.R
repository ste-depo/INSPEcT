#' @rdname modelRates
#'
#' @description
#' This method model the synthesis, degradation and processign rates after their estimation by the constructor function
#' \code{\link{newINSPEcT}}. Estimated rates are not guaranteed to optimally describes provided input data yet. 
#' To this purpose, modeled rates can be generated and genes can be assigned to a transcriptional regulatory mechanism.
#' Modeled rates can be accessed via the method \code{\link{viewModelRates}} and gene classification according 
#' to the regulatory mechanism can be accessed by \code{\link{geneClass}}. The modeling procedure can be set by the
#' user by modyging the parameters via \code{\link{modelingParams}}
#' @param object An object of class INSPEcT
#' @param seed A numeric, indicatindg the seed to be set for reproducible results
## #' @param nCores Either NULL or numeric. If numeric indicates the number of cores to be
## #' used for parallelization (if nCores=1 doesn't parallelize). If NULL takes the information
## #' from the object (see \code{\link{modelingParams}})
#' @param BPPARAM Parallelization parameters for bplapply. By default bpparam()
#' @param verbose Either NULL or logical. If logical indicates whether to output some text
#' during computation or not, if NULL  it takes the information from the object
#' (see \code{\link{modelingParams}}) (Default: NULL)
#' used for parallelization (if nCores=1 doesn't parallelize). If NULL takes the information
#' from the object (see \code{\link{modelingParams}}) (Default: NULL)
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
#' 
#' data('mycerIds10', package='INSPEcT')
#' ## models removal
#' mycerIdsThreeGenes <- removeModel(mycerIds10[1:3])
#' mycerIdsThreeGenes <- modelRates(mycerIdsThreeGenes, seed=1, BPPARAM=SerialParam())
#' ## view modeled synthesis rates
#' viewModelRates(mycerIdsThreeGenes, 'synthesis')
#' ## view gene classes
#' geneClass(mycerIdsThreeGenes)
#'
#' ## Divide a parallel computation into chunks
#' \dontrun{
#' nCores(mycerIds10) <- parallel::detectCores()
#' chunkSize <- 100
#' splitIdx <- ceiling(c(1:nGenes(mycerIds10))/chunkSize)
#' chunks <- lapply(split(mycerIds10, splitIdx), modelRates)
#' mycerIdsModeled <- do.call('combine', chunks)
#' }
setMethod('modelRates', 'INSPEcT', function(object
										  , seed=NULL
										  , BPPARAM=bpparam()
										  , verbose=NULL)
{

	No4sU <- object@params$No4sU

	if(No4sU){print("No nascent RNA data mode.")}
	if(!No4sU){print("Nascent RNA data mode.")}

	object@params$seed <- seed
	
	if(No4sU){
		nInit <- object@params$nInit
		nIter <- object@params$nIter
		cores <- object@params$cores
		Dmax <- object@params$Dmax
		Dmin <- object@params$Dmin
		na.rm <- object@params$na.rm
		verbose <- object@params$verbose
		testOnSmooth <- object@params$testOnSmooth
	
		if( length(object@model@ratesSpecs) > 0 )
			stop('Remove the model before running the model again. (See "?removeModel")')
		if( !is.null(seed) && !is.numeric(seed) )
			stop('Seed argument must be either NULL or numeric.')
		if( !is.null(verbose) && !is.logical(verbose) )
			stop('verbose argument must be either NULL or logical.')
	
		# Time transformation
		tpts <- object@tpts
		a <- .find_tt_par(tpts)
		tptsOriginal <- tpts
		tptsLinear <- .time_transf(tptsOriginal,a)
		c <- abs(min(tptsLinear))
		tptsLinear <- tptsLinear + abs(min(tptsLinear))
		concentrations <- list(total = ratesFirstGuess(object, 'total')
							 , total_var = ratesFirstGuessVar(object, 'total')
							 , preMRNA = ratesFirstGuess(object, 'preMRNA')
							 , preMRNA_var = ratesFirstGuessVar(object, 'preMRNA')
							 , mature = ratesFirstGuess(object, 'total') - ratesFirstGuess(object, 'preMRNA')
							 , mature_var = ratesFirstGuessVar(object, 'total') + ratesFirstGuessVar(object, 'preMRNA'))
	
		rates <- list(alpha=ratesFirstGuess(object, 'synthesis')
					, beta=ratesFirstGuess(object, 'degradation')
					, gamma=ratesFirstGuess(object, 'processing'))

		if(object@params$estimateRatesWith=="int")
		{

			ratesSpecs <- .inspect.engine_Integrative_No4sU(tptsOriginal = tptsOriginal
															,tptsLinear = tptsLinear
															,a = a
															,c = c
															,concentrations = concentrations
															,rates = rates
															,BPPARAM = BPPARAM
															,na.rm = na.rm
															,verbose = verbose
															,testOnSmooth = testOnSmooth
															,seed = seed
															,nInit = nInit
															,nIter = nIter)
		}else{

			ratesSpecs <- .inspect.engine_Derivative_No4sU(tptsOriginal = tptsOriginal
														  ,tptsLinear = tptsLinear
														  ,a = a
														  ,c = c
														  ,concentrations = concentrations
														  ,rates = rates
														  ,BPPARAM = BPPARAM
														  ,na.rm = na.rm
														  ,verbose = verbose
														  ,testOnSmooth = testOnSmooth
														  ,seed = seed
														  ,nInit = nInit
														  ,nIter = nIter)
		}
	
		names(ratesSpecs) <- featureNames(object@ratesFirstGuess)[seq_along(ratesSpecs)]
	
		## update and return the object
		object@model@ratesSpecs <- ratesSpecs
		object <- makeModelRates(object)

		return(object)
	}else{
		if( length(object@model@ratesSpecs) > 0 )
			stop('Remove the model before running the model again. (See "?removeModel")')
		if( !is.null(seed) && !is.numeric(seed) )
			stop('Seed argument must be either NULL or numeric.')
		# if( !is.null(nCores) && !is.numeric(nCores) )
		# 	stop('nCores argument must be either NULL or numeric.')
		if( !is.null(verbose) && !is.logical(verbose) )
			stop('verbose argument must be either NULL or logical.')
		tpts <- object@tpts
		log_shift <- .find_tt_par(tpts)
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
		ratesSpecs <- .inspect.engine(tpts, log_shift, concentrations, rates
			, nInit=object@params$nInit
			, nIter=object@params$nIter
			, na.rm=object@params$na.rm
			, BPPARAM=BPPARAM
			# , nCores=if(is.null(nCores)) object@params$nCores else nCores
			, verbose=if(is.null(verbose)) object@params$verbose else verbose
			, estimateRatesWith=object@params$estimateRatesWith
			, sigmoidDegradation=object@params$useSigmoidFun
			, sigmoidSynthesis=object@params$useSigmoidFun
			, sigmoidTotal=object@params$useSigmoidFun
			, sigmoidProcessing=object@params$useSigmoidFun
			, sigmoidPre=object@params$useSigmoidFun
			, testOnSmooth=object@params$testOnSmooth
			, seed=seed
			)
		names(ratesSpecs) <- featureNames(object@ratesFirstGuess)
		## update and return the object
		object@model@ratesSpecs <- ratesSpecs
		object <- makeModelRates(object)
		return(object)
	}
})
