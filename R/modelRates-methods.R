#' @rdname modelRates
#'
#' @description
#' This method model the synthesis, degradation and processing rates after their estimation by the constructor function
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

	# .getPriorRatesFromConcentrations <- function(tpts,concentrations,Dmax,Dmin,genesFilter=TRUE)
	# {
	# 	eiGenes <- rownames(concentrations$total)
	
	# 	mature <- concentrations$mature
	# 	premature <- concentrations$preMRNA
	
	# 	matureVariance <- concentrations$mature_var
	# 	prematureVariance <- concentrations$preMRNA_var
	
	# 	# Constant post transcriptional rates and fixed post transcriptional ratio 
	# 	k3Prior <- firstStep_NoNascent(tpts = tpts
	# 							  ,mature = mature
	# 							  ,premature = premature
	# 							  ,matureVariance = matureVariance
	# 							  ,Dmin = Dmin
	# 							  ,Dmax = Dmax)
	# 	rownames(k3Prior) <- eiGenes
	
	# 	# Constant post transcriptional rates and variable post transcriptiona ratio
	# 	fits <- t(mcsapply(1:nrow(mature), function(row)
	# 	{
	# 		unlist(
	# 			tryCatch(
	# 	    		optim(par = c(mature[row,1]/premature[row,1]*k3Prior[row,'k3'], k3Prior[row,'k3'])
	# 	    				   ,fn = secondStepError_NoNascent
	# 	        			   ,tpts = tpts
	# 	        			   ,premature = premature[row,]
	# 	        			   ,mature = mature[row,]
	# 	        			   ,matureVariance = matureVariance[row,])
	# 			,error=function(e)list(par = c(NaN,NaN), value = NaN, counts = NaN, convergence = NaN, message = "Optimization error."))[1:4])
	# 	},BPPARAM = BPPARAM))
		
	# 	fits[,3] <- pchisq(fits[,3], length(tpts)-3)
	# 	colnames(fits) <- c('k2','k3','p','counts','gradient','convergence')
	# 	rownames(fits) <- eiGenes
	
	# 	# Correction of negative priors with the median
	# 	if(genesFilter){
	# 		fits[fits[,'k2']<0,'k2'] <- NaN
	# 		fits[fits[,'k3']<0,'k3'] <- NaN
	
	# 		notFiniteRates <- !is.finite(fits[,'k2']) | !is.finite(fits[,'k3'])
	
	# 		fits[notFiniteRates,'k2'] <- median(fits[is.finite(fits[,'k2']),'k2'])
	# 		fits[notFiniteRates,'k3'] <- median(fits[is.finite(fits[,'k3']),'k3'])
	
	# 		fits[notFiniteRates,'p'] <- NaN
	# 	}
	
	# 	# Data formatting
	# 	constantModels <- list(models = fits
	# 						 , premature = premature
	# 						 , mature = mature
	# 						 , prematureVariance = prematureVariance
	# 						 , matureVariance = matureVariance)
	
	# 	ratesConstantPriors <- constantModels$models  
	
	# 	# alphaTC <- t(sapply(seq_along(ratesConstantPriors[,'k3']),function(g)
	# 	# {
	# 	# 	sapply(tpts,function(t){k1KKK_NoNascent(t,par = c(mean(mature[g,],na.rm = T),ratesConstantPriors[g,'k2'],ratesConstantPriors[g,'k3']))})
	# 	# }))
	
	# 	betaTC <- matrix(rep(ratesConstantPriors[,'k3'],length(tpts)),ncol=length(tpts))
	# 	gammaTC <- matrix(rep(ratesConstantPriors[,'k2'],length(tpts)),ncol=length(tpts))
	
	# 	prematureDer <- as.matrix(t(sapply(1:nrow(premature),function(i){
	# 		if(all(is.finite(premature[i,]))){
	# 			spfun <- splinefun(tpts, premature[i,])
	# 			return(spfun(tpts, deriv=1))
	# 		} else return(rep(NA,length(tpts)))
	# 	})))
	
	# 	alphaTC <- prematureDer + gammaTC * premature
	# 	alphaTC[alphaTC<0] <- NaN
	
	# 	return(list(alpha=alphaTC
	# 			  , beta=betaTC
	# 			  , gamma=gammaTC))
	
	# }

	NoNascent <- object@params$NoNascent

	if(NoNascent){print("No nascent RNA data mode.")}
	if(!NoNascent){print("Nascent RNA data mode.")}

	object@params$seed <- seed
	
	if(NoNascent){
		nInit <- object@params$nInit
		nIter <- object@params$nIter
		cores <- object@params$cores
		Dmax <- object@params$Dmax
		Dmin <- object@params$Dmin
		na.rm <- object@params$na.rm
		verbose <- object@params$verbose
		testOnSmooth <- object@params$testOnSmooth
	
 #		if( length(object@model@ratesSpecs) > 0 )
 #			stop('Remove the model before running the model again. (See "?removeModel")')
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
		rates <- list(
			alpha=ratesFirstGuess(object, 'synthesis')
			, beta=ratesFirstGuess(object, 'degradation')
			, gamma=ratesFirstGuess(object, 'processing')
			)

		# Code to reanalize old simualted data - TEMP
		
		#  nGenes <- nrow(concentrations$total)
		#  nSamples <- ncol(concentrations$total)
		#  
		#  concentrations$total_var <- matrix(rep(concentrations$total_var,nSamples),nrow=nGenes,ncol=nSamples,byrow=FALSE)
		#  concentrations$preMRNA_var <- matrix(rep(concentrations$preMRNA_var,nSamples),nrow=nGenes,ncol=nSamples,byrow=FALSE)
		#  concentrations$mature_var <- matrix(rep(concentrations$mature_var,nSamples),nrow=nGenes,ncol=nSamples,byrow=FALSE)
		#  
		#  rownames(concentrations$total_var) <- rownames(concentrations$preMRNA_var) <- rownames(concentrations$mature_var) <- rownames(concentrations$total)
		#  colnames(concentrations$total_var) <- colnames(concentrations$preMRNA_var) <- colnames(concentrations$mature_var) <- colnames(concentrations$total)

		# rates <- .getPriorRatesFromConcentrations(tpts=tpts,concentrations=concentrations,Dmax=Dmax,Dmin=Dmin,genesFilter=TRUE)

		# if(length(object@model@ratesSpecs)==0)
		# {
			if(object@params$estimateRatesWith=="int")
			{
				ratesSpecs <- .inspect.engine_Integrative_NoNascent(tptsOriginal = tptsOriginal
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
				ratesSpecs <- .inspect.engine_Derivative_NoNascent(tptsOriginal = tptsOriginal
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
		# }else{ratesSpecs <- object@model@ratesSpecs}

		## update and return the object
		object@model@ratesSpecs <- ratesSpecs
		
		object <- makeModelRates(object)

		return(object)
	}else{
		if( length(object@model@ratesSpecs) > 0 )
			stop('Remove the model before running the model again. (See "?removeModel")')
		if( !is.null(seed) && !is.numeric(seed) )
			stop('Seed argument must be either NULL or numeric.')
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
