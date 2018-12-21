#' @rdname makeSimDataset
#'
#' @description
#' This method generates rates and concentrations where noise is added according to the desired number of 
#' replicates that the user set as an arguments from the INSPEcT_model object that has been created by the 
#' method of the class INSPEcT \code{\link{makeSimModel}}. Rates and concentrations can be generated at the 
#' time-points of interest within the original time window. This method generates an INSPEcT object that can
#' be modeled and the performance of the modeling can be tested directly aginst the INSPEcT_model object
#' created by \code{\link{makeSimModel}}.
#' @param object An object of class INSPEcT_model, usually the output of \code{\link{makeSimModel}}
#' @param tpts A numeric vector of time points where rates and concentrations have to be evaluated
#' @param nRep Number of replicates to simulate
#' @param NoNascent A logical which, if true, makes the output of the method suitable for an analysis wothout Nascent
#' @param seed A numeric to obtain reproducible results
#' @return An object of the class ExpressionSet containing rates and concentrations
#' @seealso \code{\link{makeSimModel}}
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   nascentInspObj <- readRDS(system.file(package='INSPEcT', 'nascentInspObj.rds'))
#'   simRates<-makeSimModel(nascentInspObj, 1000, seed=1)
#'   tpts <- tpts(nascentInspObj)
#'   nascentInspObj_sim3 <- makeSimDataset(object=simRates,tpts=tpts,nRep=3,NoNascent=FALSE,seed=1)
#' }
setMethod('makeSimDataset', 'INSPEcT_model', function(object
													, tpts
													, nRep
													, NoNascent = FALSE
													, seed=NULL)
{

	if(tpts[[1]]!=object@params$tpts[[1]]){stop("makeSimDataset: new and old tpts starts must coincide.")}
	if(tpts[[length(tpts)]]!=object@params$tpts[[length(object@params$tpts)]]){stop("makeSimDataset: new and old tpts ends must coincide.")}

	ratesSpecs <- object@ratesSpecs
	nGenes <- length(ratesSpecs)

	## create the clean concentrations and rates for each gene
	ratesSpecs <- object@ratesSpecs
	nGenes <- length(ratesSpecs)
	log_shift <- find_tt_par(tpts)
	cleanRates <- lapply(1:nGenes, function(i) {
		tryCatch(
			.makeModel(tpts
					 , ratesSpecs[[i]][[1]]
					 , log_shift
					 , time_transf
					 , deSolve::ode
					 , .rxnrate)
			, error=function(e)
				.makeEmptyModel(tpts)
			)
		})
	
	## store total, preMRNA and alpha
	totalSim <- t(sapply(cleanRates, function(x) x$total))
	preMRNASim <- t(sapply(cleanRates, function(x) x$preMRNA))
	alphaSim <- t(sapply(cleanRates, function(x) x$alpha))

	#evaluate variance for the new time points
	totalFitVarianceLaw <- object@params$sim$noiseFunctions$total
	preFitVarianceLaw <- object@params$sim$noiseFunctions$pre
	alphaFitVarianceLaw <- object@params$sim$noiseFunctions$alpha

	totalSim_noisevar <- t(apply(totalSim,1,function(r)
	{
		totalFitVarianceLaw(r)
	}))
	preMRNASim_noisevar <- t(apply(preMRNASim,1,function(r)
	{
		preFitVarianceLaw(r)
	}))
	alphaSim_noisevar <- t(apply(alphaSim,1,function(r)
	{
		alphaFitVarianceLaw(r)
	}))

	addNoise <- function(signal, noiseVar){
		t(sapply(1:nrow(signal),function(r){
			sapply(1:ncol(signal),function(c)
			{
				rnorm(1,mean=signal[r,c],sd=sqrt(noiseVar[r,c]))
			})
		}))
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
	experimentalDesign <- rep(tpts, nRep)
	colnames(totalSimReplicates)<-colnames(preMRNASimReplicates)<-colnames(alphaSimReplicates) <- experimentalDesign

	nascentExpressions <- quantifyExpressionsFromTrAbundance(
											 trAbundaces = list(
											 	exonsAbundances = alphaSimReplicates
										     	, intronsAbundances = NULL
										     	)
											 , experimentalDesign = experimentalDesign
											 , varSamplingCondition = as.character(tpts[[1]])
											 , simulatedData = TRUE)

	matureExpressions <- quantifyExpressionsFromTrAbundance(
											 trAbundaces = list(
														 exonsAbundances = totalSimReplicates
										     			 , intronsAbundances = preMRNASimReplicates
										     			 )
								 			 , experimentalDesign = experimentalDesign
								 			 , varSamplingCondition = as.character(tpts[[1]]))

	if(!NoNascent)
	{
		## create the INSPEcT object
		newObject <- newINSPEcT(tpts = tpts
							  , labeling_time = 1
							  , nascentExpressions = nascentExpressions
							  , matureExpressions = matureExpressions
							  , simulatedData = TRUE)
	}else{
		## create the INSPEcT object
		newObject <- newINSPEcT(tpts = tpts
							  , labeling_time = NULL
							  , nascentExpressions = NULL
							  , matureExpressions = matureExpressions
							  , simulatedData = TRUE)
	}

	return(newObject)

})

