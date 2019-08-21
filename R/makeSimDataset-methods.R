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
#' @param object2 An object of class INSPEcT_model, usually the output of \code{\link{makeSimModel}}
#' @param tpts A numeric vector of time points where rates and concentrations have to be evaluated
#' @param nRep Number of replicates to simulate
#' @param NoNascent A logical which, if true, makes the output of the method suitable for an analysis wothout Nascent
#' @param seed A numeric to obtain reproducible results
#' @param a A numeric which represents the probability of contamination of the labeled sample due to the unlabled one
#' @param b A numeric which represents the probability of contamination of the unlabeled sample due to the labled one
#' @param tL A numeric which represents the labeling time for an ideal nascent RNA profiling, it is required for the contamination analysis
#' @return An object of the class ExpressionSet containing rates and concentrations
#' @seealso \code{\link{makeSimModel}}
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   nascentInspObj <- readRDS(system.file(package='INSPEcT', 'nascentInspObj.rds'))
#'   simRates<-makeSimModel(nascentInspObj, 1000, seed=1)
#'   tpts <- tpts(nascentInspObj)
#'   nascentSim2replicates <- makeSimDataset(object=simRates,tpts=tpts,nRep=3,NoNascent=FALSE,seed=1)
#' }
setMethod('makeSimDataset', 'INSPEcT_model', function(object
													, tpts
													, nRep
													, NoNascent = FALSE
													, seed=NULL
													, b = 0.3
													, tL = 1/6
													, noise_sd = 4.0)
{
	if(tpts[[1]]!=object@params$tpts[[1]]){stop("makeSimDataset: new and old tpts starts must coincide.")}
	if(tpts[[length(tpts)]]!=object@params$tpts[[length(object@params$tpts)]]){stop("makeSimDataset: new and old tpts ends must coincide.")}
	if(!is.null(b)|!is.null(tL))
	{
		if(!is.numeric(b)|!is.numeric(tL)){stop("makeSimDataset: b and tL must be numeric.")}
		if(b>1|b<0){stop("makeSimDataset: b must be greater than 0 and lower than 1.")}		
	}

	ratesSpecs <- object@ratesSpecs
	nGenes <- length(ratesSpecs)
	
	## I compute the corrupted data
	if(!is.null(b)&!is.null(tL))
	{
		## I solve the ode system starting from zero in order to quantify Mature and Premature RNA for the Nascent population
		cleanRates <- lapply(1:nGenes, function(i) {
			t(sapply(tpts,function(t)
			{
				outTmp <- tryCatch(
					.makeModel(tpts=c(t-tL,t)
						 , ratesSpecs[[i]][[1]]
						 # , log_shift
						 # , time_transf
						 # , deSolve::ode
						 # , .rxnrate
						 , nascent = TRUE)
					, error=function(e)
					.makeEmptyModel(tpts)
				)
				outTmp[2,]
			}))
		})

		totalSim <- t(sapply(cleanRates, function(x) unlist(x[,"total"])))
		preMRNASim <- t(sapply(cleanRates, function(x) unlist(x[,"preMRNA"])))

		L_exons <- totalSim
		L_introns <- preMRNASim

		## I solve the system for the total
		ratesSpecs <- object@ratesSpecs
		nGenes <- length(ratesSpecs)
		# log_shift <- find_tt_par(tpts)
		cleanRates <- lapply(1:nGenes, function(i) {
			tryCatch(
				.makeModel(tpts
						 , ratesSpecs[[i]][[1]]
						 # , log_shift
						 # , time_transf
						 # , deSolve::ode
						 # , .rxnrate
						 , nascent = FALSE)
				, error=function(e)
					.makeEmptyModel(tpts)
				)
			})
	
		totalSim <- t(sapply(cleanRates, function(x) x$total))
		preMRNASim <- t(sapply(cleanRates, function(x) x$preMRNA))
		
		T_exons <- totalSim
		T_introns <- preMRNASim

		U_exons <- T_exons - L_exons
		U_introns <- T_introns - L_introns

		if(is.null(rownames(L_exons)) & 
		   is.null(rownames(L_introns)) & 
		   is.null(rownames(T_exons)) & 
		   is.null(rownames(T_introns)) & 
		   is.null(rownames(U_exons)) & 
		   is.null(rownames(U_introns))){
		rownames(L_exons) <- rownames(L_introns) <- rownames(T_exons) <- 
		rownames(T_introns) <- rownames(U_exons) <- rownames(U_introns)	<- as.character(1:nrow(L_exons))
		}

		## I remove from the analysis those genes with negative expressions
		genesTmp <- apply(T_exons,1,function(r)all(r>0)) &
		apply(T_introns,1,function(r)all(r>0)) &
		apply(L_exons,1,function(r)all(r>0)) &
		apply(L_introns,1,function(r)all(r>0)) &
		apply(U_exons,1,function(r)all(r>0)) &
		apply(U_introns,1,function(r)all(r>0))

		print(paste0("makeSimDataset: ",table(genesTmp)["TRUE"]," genes suitable for the corrupted analysis."))

		T_exons <- T_exons[genesTmp,]
		T_introns <- T_introns[genesTmp,]

		L_exons <- L_exons[genesTmp,]
		L_introns <- L_introns[genesTmp,]

		U_exons <- U_exons[genesTmp,]
		U_introns <- U_introns[genesTmp,]

	### NEW CORRUPTION START ###
		set.seed(seed)

		totalFitVarianceLaw <- object@params$sim$noiseFunctions$total
		preFitVarianceLaw <- object@params$sim$noiseFunctions$pre

		## Evaluate sigma X from the sigma of b
		# X_noise_sd <- mean((L_exons+L_introns)/(U_exons+U_introns)*((1/(1-b))-(b/(1-b)^2))*noise_sd)
		## Standard deviation on X from the sd of b
		X <- (b/(1 - b))*apply(((L_exons+L_introns)/(U_exons+U_introns)),1,mean)
		# b <- X/(X+(L_exons+L_introns)/(U_exons+U_introns))
		X_noisy <- X*rnorm(length(X), 1, sd = noise_sd)#matrix(rnorm(length(X), 1, sd = noise_sd),nrow=nrow(X),ncol=ncol(X))
		
		X_noisy[X_noisy<0] <- 0
		X_noisy[X_noisy>1] <- 1

		b_noisy <- X_noisy/(X_noisy+apply((L_exons+L_introns)/(U_exons+U_introns),1,mean))

		message(paste0("b_noisy momenta: ",signif(mean(b_noisy),3)," - ",signif(sd(b_noisy),3)))

		## calculate the noise distribution
		# X_noisy <- rnorm(1000*nrow(L_exons), X, sd = X_noise_sd)

		# X_noisy <- X_noisy[X_noisy>=0&X_noisy<=1]
		# if(length(X_noisy)<nrow(L_exons))
		# {
		# 	message("makeSimDataset: warning, very large noise coefficient standard deviation! ")
		# 	X_noisy <- c(X_noisy,rep(X,length.out=nrow(L_exons) - length(X_noisy)))
		# }
		# X_noisy <- sample(X_noisy,nrow(L_exons))

		# X_noisy x<- (L_exons+L_introns)/(U_exons+U_introns)*((1/(1-b))-(b/(1-b)^2))*noise_sd
		# b_noisy <- 	X_noisy/(X_noisy+mean((L_exons+L_introns)/(U_exons+U_introns)))


		### We set a standard deviation on X, instead of b, because the boundaries for this parameter
		### are more clear and identical for all the genes (0 <= X <= 1) while for b the boundaries
		### are gene dependent (0 <= b <= 1 - (L)/(L + U))
		#
		# ## Standard deviation on b
		# b_noisy <- rnorm(1000*nrow(L_exons), b, sd = noise_sd)
		# b_noisy <- b_noisy[b_noisy>=0&b_noisy<=1]

		# if(length(b_noisy)<nrow(L_exons))
		# {
		# 	message("makeSimDataset: warning, very large noise coefficient standard deviation! ")
		# 	b_noisy <- c(b_noisy,rep(b,length.out=nrow(L_exons) - length(b_noisy)))
		# }

		# b_noisy <- sample(b_noisy,nrow(L_exons))

		# X_noisy <- (b_noisy/(1 - b_noisy))*apply(((L_exons+L_introns)/(U_exons+U_introns)),1,mean)

	### NEW CORRUPTION END ###
	
		E_exons <- L_exons + X_noisy*U_exons
		E_introns <- L_introns + X_noisy*U_introns

		E_exons_var <- t(apply(E_exons,1,function(r){totalFitVarianceLaw(r)}))
		E_introns_var <- t(apply(E_introns,1,function(r){preFitVarianceLaw(r)}))

		T_exons_var <- t(apply(T_exons,1,function(r){totalFitVarianceLaw(r)}))
		T_introns_var <- t(apply(T_introns,1,function(r){preFitVarianceLaw(r)}))
		
		## Evaluation of the new simulated data
		E_data <- list("exonsExpressions"=E_exons, "intronsExpressions"=E_introns, "exonsVariance"=E_exons_var, "intronsVariance"=E_introns_var)
		T_data <- list("exonsExpressions"=T_exons, "intronsExpressions"=T_introns, "exonsVariance"=T_exons_var, "intronsVariance"=T_introns_var)	
	
		corruptedInspObj<-newINSPEcT(tpts=tpts
									,labeling_time=tL
									,nascentExpressions=E_data
									,matureExpressions=T_data)

		totalSim <- ratesFirstGuess(corruptedInspObj,"total")
		preMRNASim <- ratesFirstGuess(corruptedInspObj,"preMRNA")
		alphaSim <- ratesFirstGuess(corruptedInspObj,"synthesis")

	### NEW CORRUPTION END ###
	}else{
		## create the clean concentrations and rates for each gene
		cleanRates <- lapply(1:nGenes, function(i) {
			tryCatch(
				.makeModel(tpts
						 , ratesSpecs[[i]][[1]]
						 # , log_shift
						 # , time_transf
						 # , deSolve::ode
						 # , .rxnrate
						 , nascent = FALSE)
				, error=function(e)
					.makeEmptyModel(tpts)
				)
			})
	
		## store total, preMRNA and alpha
		totalSim <- t(sapply(cleanRates, function(x) x$total))
		preMRNASim <- t(sapply(cleanRates, function(x) x$preMRNA))
		alphaSim <- t(sapply(cleanRates, function(x) x$alpha))

		if(is.null(rownames(totalSim))){rownames(totalSim) <- 1:nrow(totalSim)}
	}

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
	rownames(totalSimReplicates) <- 1:nrow(totalSimReplicates)
	preMRNASimReplicates <- do.call('cbind', 
		lapply(1:nRep, function(i) addNoise(preMRNASim,preMRNASim_noisevar)))
	rownames(preMRNASimReplicates) <- 1:nrow(preMRNASimReplicates)
	alphaSimReplicates <- do.call('cbind', 
		lapply(1:nRep, function(i) addNoise(alphaSim,alphaSim_noisevar)))
	rownames(alphaSimReplicates) <- 1:nrow(alphaSimReplicates)
	experimentalDesign <- rep(tpts, nRep)
	colnames(totalSimReplicates)<-colnames(preMRNASimReplicates)<-colnames(alphaSimReplicates) <- experimentalDesign

	# Required to keep track of the genes names after the dataset reduction caused by the contamination
	rownames(totalSimReplicates)<-rownames(preMRNASimReplicates)<-rownames(alphaSimReplicates)<-rownames(totalSim)

	nascentExpressions <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances = alphaSimReplicates
																			  , intronsAbundances = NULL)
														   , experimentalDesign = experimentalDesign
														   , varSamplingCondition = as.character(tpts[[1]])
														   , simulatedData = TRUE)

	matureExpressions <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances = totalSimReplicates
																			 , intronsAbundances = preMRNASimReplicates)
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

