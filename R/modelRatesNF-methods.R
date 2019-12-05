#' @rdname modelRatesNF
#'
#' @description
#' This method compute confidence intervals for the rates of synthesis, degradation and processing estimated by 
#' \code{\link{newINSPEcT}} that will be used to estimate the variability of each rate in \code{\link{ratePvals}}
#' method.
#' @param object An object of class INSPEcT
#' @param BPPARAM Parallelization parameters for bplapply. By default bpparam()
#' @return An object of class INSPEcT with modeled rates
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#' 	nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' 	## models removal
#' 	nascentInspObjThreeGenes <- removeModel(nascentInspObj10[1:3])
#' 	nascentInspObjThreeGenes <- modelRatesNF(nascentInspObjThreeGenes, 
#' 	  BPPARAM=SerialParam())
#' 	## view modeled synthesis rates
#' 	viewModelRates(nascentInspObjThreeGenes, 'synthesis')
#' 	## view gene classes
#' 	geneClass(nascentInspObjThreeGenes)
#' }
setMethod('modelRatesNF', 'INSPEcT', function(object, BPPARAM=SerialParam())
{
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	if( length(object@model@ratesSpecs) > 0 )
		stop('Remove the model before running the model again. (See "?removeModel")')

	# llConfidenceThreshold <- object@model@params$logLikelihoodConfidenceThreshold
	# if(is.null(llConfidenceThreshold)) llConfidenceThreshold <- 0.95
	# llConfidenceThreshold <- qchisq(llConfidenceThreshold,1)
	llConfidenceThreshold <- qchisq(0.95,1)

	NoNascent <- object@NoNascent

	if(NoNascent){message("No nascent RNA data mode.")}
	if(!NoNascent){message("Nascent RNA data mode.")}

	geneNames <- featureNames(object)
	tpts <- tpts(object)
	premature         <- ratesFirstGuess(object,'preMRNA')
	prematureVariance <- ratesFirstGuessVar(object,'preMRNA')
	total             <- ratesFirstGuess(object,'total')
	totalVariance     <- ratesFirstGuessVar(object,'total')

	if( NoNascent ) {

		confidenceIntervals <- classify_rates_from_priors_v5(
		      tpts,
		      premature,
		      total,
		      prematureVariance,
		      totalVariance,
		      BPPARAM=BPPARAM
		      )		

	} else {

		synthesis         <- ratesFirstGuess(object,'synthesis')
		synthesisVariance <- ratesFirstGuessVar(object,'synthesis')
		processing        <- ratesFirstGuess(object,'processing')
		degradation       <- ratesFirstGuess(object,'degradation')

		confidenceIntervals <- classify_rates_from_priors_4sU_v5(
		      tpts,
		      premature,
		      total,
		      synthesis,
		      processing,
		      degradation,
		      prematureVariance,
		      totalVariance,
		      synthesisVariance,
		      BPPARAM=BPPARAM
		      )

	}

	# make an "ExpressionSet" object containing optimized rates
	nTpts <- length(tpts)
	tptsLabels = if( is.numeric(tpts) ) signif(tpts,9) else tpts
	exprData <- cbind(total
					, premature
					, t(sapply(confidenceIntervals, function(g) g[[1]][,'opt']))
					, t(sapply(confidenceIntervals, function(g) g[[2]][,'opt']))
					, t(sapply(confidenceIntervals, function(g) g[[3]][,'opt'])))

	pData <- data.frame(feature=c(rep('total',nTpts)
								, rep('preMRNA',nTpts)
								, rep('synthesis',nTpts)
								, rep('processing',nTpts) 
								, rep('degradation',nTpts)
								)
					  , time=rep(tpts, 5))

	colnames(exprData) <- paste(pData$feature, tptsLabels, sep='_')
	rownames(pData) <- colnames(exprData)
	phenoData <- new('AnnotatedDataFrame', data=pData)

	modelRates <- ExpressionSet(assayData=exprData
								   , phenoData=phenoData)
	featureNames(modelRates) <- geneNames
	
	object@modelRates <- modelRates
	object@NF <- TRUE
	object <- setConfidenceIntervals(object=object,confidenceIntervals=confidenceIntervals)
	object <- calculateRatePvals(object)

	return(object)


})

lin2logPar <- function(par)
	c(par[1], log2(par[-1]/par[1]))

log2linPar <- function(par)
	c(par[1], par[1]*2^par[-1])

score_and_par <- function(conf_int) {

	k_score_fun <- function(k, rate_conf_int)
	{
		mean(apply(rate_conf_int, 1, function(x) {
			if( k < x[2] ) (k - x[2])^2/(x[2]-x[1])^2 
				else (k - x[2])^2/(x[2]-x[3])^2
		}), na.rm=T)
	}

	k_scores_out <- lapply(conf_int, function(gene) lapply(gene, function(rate_conf_int) {
		k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
		if(!is.finite(k_start)) return(list(par=NaN, value=NaN))
		if(!is.finite(any(rate_conf_int[,1]))) return(list(par=NaN, value=NaN))
		if(!is.finite(any(rate_conf_int[,3]))) return(list(par=NaN, value=NaN))
		optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int)	
		}))

	k_par <- t(sapply(k_scores_out, function(x) sapply(x,'[[','par')))
	k_value <- t(sapply(k_scores_out, function(x) sapply(x,'[[','value')))

	return(list(par=k_par, score=k_value))

}

mirror_low <- function(rate_conf_int) 2*rate_conf_int[,2] - rate_conf_int[,1]
mirror_hig <- function(rate_conf_int) 2*rate_conf_int[,2] - rate_conf_int[,3]

rate_conf_int_impute <- function(rate_conf_int) {
	na_low <- is.na(rate_conf_int[,1])
	if( any(na_low) ) rate_conf_int[na_low,1] <- mirror_hig(rate_conf_int)[na_low]
	na_hig <- is.na(rate_conf_int[,3])
	if( any(na_hig) ) rate_conf_int[na_hig,3] <- mirror_low(rate_conf_int)[na_hig]
	return(rate_conf_int)
}

####

errorOptim <- function(N, e) return(list(
	par = rep(NaN,N),
	value = NaN,
	counts = c("function"=NaN,"gradient"=NaN),
	convergence = NaN,
	message = e
	))

# Unit: microseconds
#                                                                               expr
#                       expressionProfilesFromPiecedFunctions(parameters2, 1:10, 10)
#    expressionProfilesFromPiecedFunctionsOneNA(parameters2[-29],      1:10, 29, 10)
#  expressionProfilesFromPiecedFunctionsManyNA(parameters, 1:10,      na_ids, 10, 3)
#     min      lq     mean  median      uq     max neval
#  71.511 72.1650 73.31720 72.4740 72.9205 120.801   100
#  78.201 79.7175 80.61715 80.0850 80.7440  98.437   100
#  80.183 81.2655 82.64758 81.5655 82.1880 129.250   100
 
###

rates_from_napars = function(parameters, na_ids, N, M) {
	parametersNew = parameters[1:na_ids[1]]
	for( i in seq_along(na_ids)[-1] ) {
		parametersNew = c(parametersNew, parameters[(na_ids[i-1]-i+2):(na_ids[i]-i+1)])
	}
	parametersNew = c(parametersNew, parameters[(na_ids[M]-M+1):(3*N-M)])
	k1 = parametersNew[1:N]
	k2 = parametersNew[(N+1):(2*N)]
	k3 = parametersNew[(2*N+1):(3*N)]
	return(list(k1,k2,k3))		
}

expressionProfilesFromPiecedFunctionsManyNA <- function(parameters, tpts, na_ids, N, M)
{

	inferPandMfromPiecedFunctions <- function(tpts, alpha, beta, gamma)
	{
		Ptmp <- seq_along(tpts)
		Ttmp <- seq_along(tpts)
	
		Ptmp[1] <- alpha[1]/gamma[1]
		Ttmp[1] <- alpha[1]/gamma[1] + alpha[1]/beta[1]

		for(t in seq_along(tpts)[-1])
		{
			mAlpha <- (alpha[t-1] - alpha[t] ) / (tpts[t-1] - tpts[t])
			qAlpha <- alpha[t-1] - mAlpha * tpts[t-1]

			Ptmp[t] <- Ptmp[t-1]*exp(-gamma[t]*(tpts[t] - tpts[t-1])) + (
			(mAlpha*tpts[t]*gamma[t] + qAlpha*gamma[t] - mAlpha ) / (gamma[t]^2) - 
			(mAlpha*tpts[t-1]*gamma[t] + qAlpha*gamma[t] - mAlpha ) * exp(-gamma[t]*(tpts[t] - tpts[t-1])) / (gamma[t]^2)
			)

			mPreMRNA <- (Ptmp[t-1] - Ptmp[t] ) / (tpts[t-1] - tpts[t] )
			qPreMRNA <- Ptmp[t-1] - mPreMRNA * tpts[t-1]

			Ttmp[t] <- Ttmp[t-1] * exp(-beta[t]*(tpts[t]-tpts[t-1])) + 
			((mAlpha*tpts[t]*beta[t] + qAlpha*beta[t] - mAlpha ) / (beta[t]^2 ) - (mAlpha*tpts[t-1]*beta[t] + qAlpha*beta[t] - mAlpha ) * exp(-beta[t]*(tpts[t]-tpts[t-1])) / (beta[t]^2 )) +
			beta[t]*((mPreMRNA*tpts[t]*beta[t] + qPreMRNA*beta[t] - mPreMRNA ) / (beta[t]^2 ) - (mPreMRNA*tpts[t-1]*beta[t] + qPreMRNA*beta[t] - mPreMRNA ) * exp(-beta[t]*(tpts[t]-tpts[t-1])) / (beta[t]^2 ))

		}

		return(c(unlist(Ttmp)-unlist(Ptmp),unlist(Ptmp)))
	}

	parametersNew = parameters[1:na_ids[1]]
	for( i in seq_along(na_ids)[-1] ) {
		parametersNew = c(parametersNew, parameters[(na_ids[i-1]-i+2):(na_ids[i]-i+1)])
	}
	parametersNew = c(parametersNew, parameters[(na_ids[M]-M+1):(3*N-M)])
	k1 = parametersNew[1:N]
	k2 = parametersNew[(N+1):(2*N)]
	k3 = parametersNew[(2*N+1):(3*N)]
	return(inferPandMfromPiecedFunctions(tpts = tpts, alpha = k1, beta = k3, gamma = k2))

}

expressionProfilesFromPiecedFunctions <- function(parameters, tpts, N)
{

	inferPandMfromPiecedFunctions <- function(tpts, alpha, beta, gamma)
	{
		Ptmp <- seq_along(tpts)
		Ttmp <- seq_along(tpts)
	
		Ptmp[1] <- alpha[1]/gamma[1]
		Ttmp[1] <- alpha[1]/gamma[1] + alpha[1]/beta[1]

		for(t in seq_along(tpts)[-1])
		{
			mAlpha <- (alpha[t-1] - alpha[t] ) / (tpts[t-1] - tpts[t])
			qAlpha <- alpha[t-1] - mAlpha * tpts[t-1]
			
			Ptmp[t] <- Ptmp[t-1]*exp(-gamma[t]*(tpts[t] - tpts[t-1])) + (
			(mAlpha*tpts[t]*gamma[t] + qAlpha*gamma[t] - mAlpha ) / (gamma[t]^2) - 
			(mAlpha*tpts[t-1]*gamma[t] + qAlpha*gamma[t] - mAlpha ) * exp(-gamma[t]*(tpts[t] - tpts[t-1])) / (gamma[t]^2)
			)

			mPreMRNA <- (Ptmp[t-1] - Ptmp[t] ) / (tpts[t-1] - tpts[t] )
			qPreMRNA <- Ptmp[t-1] - mPreMRNA * tpts[t-1]

			Ttmp[t] <- Ttmp[t-1] * exp(-beta[t]*(tpts[t]-tpts[t-1])) + 
			((mAlpha*tpts[t]*beta[t] + qAlpha*beta[t] - mAlpha ) / (beta[t]^2 ) - (mAlpha*tpts[t-1]*beta[t] + qAlpha*beta[t] - mAlpha ) * exp(-beta[t]*(tpts[t]-tpts[t-1])) / (beta[t]^2 )) +
			beta[t]*((mPreMRNA*tpts[t]*beta[t] + qPreMRNA*beta[t] - mPreMRNA ) / (beta[t]^2 ) - (mPreMRNA*tpts[t-1]*beta[t] + qPreMRNA*beta[t] - mPreMRNA ) * exp(-beta[t]*(tpts[t]-tpts[t-1])) / (beta[t]^2 ))

		}

		return(c(unlist(Ttmp)-unlist(Ptmp),unlist(Ptmp)))
	}

	k1 = parameters[1:N]
	k2 = parameters[(N+1):(2*N)]
	k3 = parameters[(2*N+1):(3*N)]
	return(inferPandMfromPiecedFunctions(tpts = tpts, alpha = k1, beta = k3, gamma = k2))

}

piecedFunctionLogLikelihoodManyNA4sU <- function(parameters, tpts, experimentalP, experimentalM, experimentalK1, varianceP, varianceM, varianceK1, na_ids, N, M)
{
	modeledProfiles <- expressionProfilesFromPiecedFunctionsManyNA(parameters = parameters, tpts = tpts, na_ids=na_ids, N=N, M=M)

	logLikelihoodFunction(experiment = c(experimentalM, experimentalP, experimentalK1),
						  model = c(modeledProfiles, parameters[1:N]),
						  variance = c(varianceM, varianceP, varianceK1))
}

piecedFunctionLogLikelihood4sU <- function(parameters, tpts, experimentalP, experimentalM, experimentalK1, varianceP, varianceM, varianceK1, N)
{
	modeledProfiles <- expressionProfilesFromPiecedFunctions(parameters=parameters, tpts = tpts, N = N)

	logLikelihoodFunction(experiment = c(experimentalM, experimentalP, experimentalK1),
						  model = c(modeledProfiles, parameters[1:N]),
						  variance = c(varianceM, varianceP, varianceK1))		
}

perturbedLogLikelihood4sU_NF <- function(parameter,name,parameters,tpts,experimentalP,experimentalM,experimentalK1,varianceP,varianceM,varianceK1,N)
{

	perturbedParameters <- parameters
	perturbedParameters[name] <- parameter

	piecedFunctionLogLikelihood4sU(parameters = perturbedParameters
							  , tpts = tpts
							  , experimentalP = experimentalP
							  , experimentalM = experimentalM
							  , experimentalK1 = experimentalK1
							  , varianceP = varianceP
							  , varianceM = varianceM
							  , varianceK1 = varianceK1
							  , N = N
							  )	
}

piecedFunctionLogLikelihoodManyNA <- function(parameters, tpts, experimentalP, experimentalM, varianceP, varianceM, na_ids, N, M)
{
	parameters <- c(log2linPar(parameters[1:N]), log2linPar(parameters[(N+1):(2*N)]), log2linPar(parameters[(2*N+1):(3*N)]))
	modeledProfiles <- expressionProfilesFromPiecedFunctionsManyNA(parameters = parameters, tpts = tpts, na_ids=na_ids, N=N, M=M)

	logLikelihoodFunction(experiment = c(experimentalM, experimentalP),
						  model = modeledProfiles,
						  variance = c(varianceM, varianceP))		
}

# piecedFunctionLogLikelihoodOneNA <- function(parameters, tpts, experimentalP, experimentalM, varianceP, varianceM, na_id, N)
# {
# 	modeledProfiles <- expressionProfilesFromPiecedFunctionsOneNA(parameters = parameters, tpts = tpts, na_id=na_id, N=N)
# 
# 	logLikelihoodFunction(experiment = c(experimentalM, experimentalP),
# 						  model = modeledProfiles,
# 						  variance = c(varianceM, varianceP))		
# }

piecedFunctionLogLikelihood <- function(parameters, tpts, experimentalP, experimentalM, varianceP, varianceM, N)
{
	parameters <- c(log2linPar(parameters[1:N]), log2linPar(parameters[(N+1):(2*N)]), log2linPar(parameters[(2*N+1):(3*N)]))
	modeledProfiles <- expressionProfilesFromPiecedFunctions(parameters=parameters, tpts = tpts, N = N)

	logLikelihoodFunction(experiment = c(experimentalM, experimentalP),
						  model = modeledProfiles,
						  variance = c(varianceM, varianceP))		
}

perturbedLogLikelihood_NF <- function(parameter,name,parameters,tpts,experimentalP,experimentalM,varianceP,varianceM,N)
{

	perturbedParameters <- parameters
	perturbedParameters[name] <- parameter

	piecedFunctionLogLikelihood(parameters = perturbedParameters
							  , tpts = tpts
							  , experimentalP = experimentalP
							  , experimentalM = experimentalM
							  , varianceP = varianceP
							  , varianceM = varianceM
							  , N = N
							  )	
}

logLikelihoodCIerror_NF <- function(parameter,name,parameters,tpts,experimentalP,experimentalM,varianceP,varianceM,N,confidenceThreshold)
{

	maximumLogLikelihoodTmp <- piecedFunctionLogLikelihood(parameters = parameters
														 , tpts = tpts
														 , experimentalP = experimentalP
														 , experimentalM = experimentalM
														 , varianceP = varianceP
														 , varianceM = varianceM
														 , N = N
														 )
	
	perturbedLogLikelihoodTmp <- perturbedLogLikelihood_NF(parameter = parameter
													  , name = name
													  , parameters = parameters
													  , tpts = tpts
													  , experimentalP = experimentalP
													  , experimentalM = experimentalM
													  , varianceP = varianceP
													  , varianceM = varianceM
													  , N = N
													  )

	return(abs(confidenceThreshold - 2*(maximumLogLikelihoodTmp - perturbedLogLikelihoodTmp)))
}

estimate_priors <- function(
	tpts,
	premature,
	mature,
	prematureVariance,
	matureVariance,
	eiGenes=NULL,
	genesFilter=TRUE,
	Dmin=1e-06,
	Dmax=10,
	BPPARAM=bpparam()
	) 
{
	k3Prior <- firstStep_NoNascent(tpts = tpts
							  ,mature = mature
							  ,premature = premature
							  ,matureVariance = matureVariance
							  ,Dmin = Dmin
							  ,Dmax = Dmax)
	rownames(k3Prior) <- eiGenes

	k1_prior_T0 <- premature[,1]^(3/4)
	k2_prior_T0 <- premature[,1]^(-1/4)
	k3_prior_T0 <- k1_prior_T0/mature[,1]

	norm_to_k3timescale <- mean(k3Prior[,1])/mean(k3_prior_T0)

	k1_prior_T0 <- k1_prior_T0 * norm_to_k3timescale
	k2_prior_T0 <- k2_prior_T0 * norm_to_k3timescale
	k3_prior_T0 <- k3_prior_T0 * norm_to_k3timescale

	tcder = function(x, y) {
		N = length(x)
		der = rep(0, N)
		der[1] = 0 # (y[2]-y[1])/(x[2]-x[1])
		der[N] = (y[N]-y[N-1])/(x[N]-x[N-1])
		der[2:(N-1)] = sapply(2:(N-1), function(j) {
			(y[j+1]-y[j-1])/(x[j+1]-x[j-1])
			})
		return(der)
	}
	prematureDer <- as.matrix(t(sapply(1:nrow(premature),function(i)
		tcder(tpts, premature[i,]))))

	matureDer <- as.matrix(t(sapply(1:nrow(mature),function(i)
		tcder(tpts, mature[i,]))))

	total <- premature + mature

	## fix k1

	ratesFromAlpha = function(alphaTC) {

		gammaTC <- ( alphaTC - prematureDer ) / premature

		#Evaluate beta as constant between intervals

		betaT0 <- alphaTC[,1] / mature[,1]	
		betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, total, premature, 
					maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
					)
		betaTC <- cbind(betaT0, 
			sapply(betaOut, function(x) sapply(x, '[[', 'root'))
		)
		betaEstimPrec <- cbind(0,
			sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
		)

		return(list(alphaTC, gammaTC, betaTC))

	}

	alphaTC <- matrix(rep(k1_prior_T0,length(tpts)),ncol=length(tpts))

	ratesOut = ratesFromAlpha(alphaTC)
	alphaTC = ratesOut[[1]]
	gammaTC = ratesOut[[2]]
	betaTC = ratesOut[[3]]

	k2varRates <- list(alphaTC, gammaTC, betaTC)

	## fix k2 

	ratesFromGamma = function(gammaTC) {

		alphaTC <- prematureDer + gammaTC * premature

		#Evaluate beta as constant between intervals

		betaT0 <- alphaTC[,1] / mature[,1]	
		betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, total, premature, 
					maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
					)
		betaTC <- cbind(betaT0, 
			sapply(betaOut, function(x) sapply(x, '[[', 'root'))
		)
		betaEstimPrec <- cbind(0,
			sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
		)

		#Evaluate gamma as constant between intervals

		gammaT0 <- alphaTC[,1] / premature[,1]
		gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, premature, 
			maxGamma=quantile(gammaT0,na.rm=TRUE,probs=.99)*10, BPPARAM=BPPARAM
			)
		gammaTC <- cbind(gammaT0, 
			sapply(gammaOut, function(x) sapply(x, '[[', 'root'))
			)
		gammaEstimPrec <- cbind(0, 
			sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec'))
			)

		return(list(alphaTC, gammaTC, betaTC))

	}

	gammaTC <- matrix(rep(k2_prior_T0,length(tpts)),ncol=length(tpts))

	ratesOut = ratesFromGamma(gammaTC)
	alphaTC = ratesOut[[1]]
	gammaTC = ratesOut[[2]]
	betaTC = ratesOut[[3]]

	k1varRates <- list(alphaTC, gammaTC, betaTC)

	## select between the two

	k1var_tot_fc <- sapply(1:nrow(mature), function(i) {
		rates <- rbind(k1varRates[[1]][i,], k1varRates[[2]][i,], k1varRates[[3]][i,])
		rates[rates<0] <- NaN
		fc <- abs(log(rates[,-1]) - log(rates[,1]))
		mean(fc, na.rm=T)
		})
	k1var_tot_fc[is.na(k1var_tot_fc)] <- max(k1var_tot_fc[!is.na(k1var_tot_fc)])

	k2var_tot_fc <- sapply(1:nrow(mature), function(i) {
		rates <- rbind(k2varRates[[1]][i,], k2varRates[[2]][i,], k2varRates[[3]][i,])
		rates[rates<0] <- NaN
		fc <- abs(log(rates[,-1]) - log(rates[,1]))
		mean(fc, na.rm=T)
		})
	k2var_tot_fc[is.na(k2var_tot_fc)] <- max(k2var_tot_fc[!is.na(k2var_tot_fc)])

	N <- length(tpts)
	select <- k1var_tot_fc < k2var_tot_fc
	alphaTC <- t(sapply(seq_along(select), function(i) {
		if( is.na(select[i]) ) return(rep(NaN, N))
		if( select[i] ) k1varRates[[1]][i,] else k2varRates[[1]][i,]
		}))
	gammaTC <- t(sapply(seq_along(select), function(i) {
		if( is.na(select[i]) ) return(rep(NaN, N))
		if( select[i] ) k1varRates[[2]][i,] else k2varRates[[2]][i,]
		}))
	betaTC <- t(sapply(seq_along(select), function(i) {
		if( is.na(select[i]) ) return(rep(NaN, N))
		if( select[i] ) k1varRates[[3]][i,] else k2varRates[[3]][i,]
		}))

	return(list(alphaTC, gammaTC, betaTC))

}

logLikelihoodCIerror4sU_NF <- function(parameter,name,parameters,tpts,experimentalP,experimentalM,experimentalK1,varianceP,varianceM,varianceK1,N,confidenceThreshold)
{

	maximumLogLikelihoodTmp <- piecedFunctionLogLikelihood4sU(parameters = parameters
														 , tpts = tpts
														 , experimentalP = experimentalP
														 , experimentalM = experimentalM
														 , experimentalK1 = experimentalK1
														 , varianceP = varianceP
														 , varianceM = varianceM
														 , varianceK1 = varianceK1
														 , N = N
														 )
	
	perturbedLogLikelihoodTmp <- perturbedLogLikelihood4sU_NF(parameter = parameter
													  , name = name
													  , parameters = parameters
													  , tpts = tpts
													  , experimentalP = experimentalP
													  , experimentalM = experimentalM
													  , experimentalK1 = experimentalK1
													  , varianceP = varianceP
													  , varianceM = varianceM
													  , varianceK1 = varianceK1
													  , N = N
													  )

	return(abs(confidenceThreshold - 2*(maximumLogLikelihoodTmp - perturbedLogLikelihoodTmp)))
}

classify_rates_from_priors_4sU_v5 <- function(
	tpts,
	premature,
	total,
	synthesis,
	gammaTC,
	betaTC,
	prematureVariance,
	totalVariance,
	synthesisVariance,
	eiGenes=NULL,
	genesFilter=TRUE,
	llConfidenceThreshold=qchisq(.95,1),
	Dmin=1e-06,
	Dmax=10,
	BPPARAM=bpparam()
	) 
{
	mature <- total - premature
	matureVariance <- totalVariance + prematureVariance

	message("Rates optimization through the minimum likelihood.")

	N <- length(tpts)
	gene_number <- nrow(premature)

	gammaTC[gammaTC<0] <- NA
	betaTC[betaTC<0]   <- NA

	N_na_genes = apply(is.na(cbind(synthesis,gammaTC,betaTC)),1,function(x) length(which(x)))

	# initalize the list of results
	newParams <- as.list(rep(NA,gene_number))
	
	# model genes without NA values
	no_na_genes = which(N_na_genes==0)
	newParams[no_na_genes] <- bplapply(no_na_genes,function(r)
	{
		parameters <- unname(c(synthesis[r,],gammaTC[r,],betaTC[r,]))
		tryCatch({
			optim(par = parameters
				, fn = piecedFunctionLogLikelihood4sU
				, tpts = tpts
				, experimentalP = premature[r,]
				, experimentalM = mature[r,]
				, experimentalK1 = synthesis[r,]
				, varianceP = prematureVariance[r,]
				, varianceM = matureVariance[r,]
				, varianceK1 = synthesisVariance[r,]
				, N = length(tpts)
				, control = list(maxit = 2000, fnscale=-1)
				, method = "BFGS")
		},error=function(e)return(errorOptim(N*3,e)))
	}, BPPARAM=BPPARAM)

	alphaTC <- t(sapply(newParams,function(g)g[[1]][1:N]))
	gammaTC <- t(sapply(newParams,function(g)g[[1]][(N+1):(2*N)]))
	betaTC  <- t(sapply(newParams,function(g)g[[1]][(2*N+1):(3*N)]))

	alphaTC[alphaTC<0] <- NA
	gammaTC[gammaTC<0] <- NA
	betaTC[betaTC<0] <- NA

	N_na_genes = apply(is.na(cbind(alphaTC,gammaTC,betaTC)),1,
		function(x) length(which(x)))

	n_iter = 10
	solved_genes = numeric(n_iter)

	if( any(N_na_genes>0 & N_na_genes<3*N) )
		for( iter in 1:n_iter ) {

			many_na_genes <- which(N_na_genes>0 & N_na_genes<3*N)
			newParams[many_na_genes] <- bplapply(many_na_genes,function(r)
			{
				parameters <- na.omit(c(alphaTC[r,],gammaTC[r,],betaTC[r,]))
				na_ids <- as.vector(attr(parameters, 'na.action'))
				M <- N_na_genes[r]
				tryCatch({
					optOut <- optim(par = as.vector(parameters)
						, fn = piecedFunctionLogLikelihoodManyNA4sU
						, tpts = tpts
						, experimentalP = premature[r,]
						, experimentalM = mature[r,]
						, experimentalK1 = synthesis[r,]
						, varianceP = prematureVariance[r,]
						, varianceM = matureVariance[r,]
						, varianceK1 = synthesisVariance[r,]
						, na_ids = na_ids
						, N = N
						, M = M
						, control = list(maxit = 2000, fnscale=-1)
						, method = "BFGS")
					## expand parameters
					optOut$par = unlist(rates_from_napars(optOut$par, na_ids, N, M))
					return(optOut)
				},error=function(e)return(errorOptim(N*3,e)))
			}, BPPARAM=BPPARAM)

			alphaTC <- t(sapply(newParams,function(g)g[[1]][1:N]))
			gammaTC <- t(sapply(newParams,function(g)g[[1]][(N+1):(2*N)]))
			betaTC  <- t(sapply(newParams,function(g)g[[1]][(2*N+1):(3*N)]))

			alphaTC[alphaTC<0] <- NA
			gammaTC[gammaTC<0] <- NA
			betaTC[betaTC<0] <- NA

			N_na_genes = apply(is.na(cbind(alphaTC,gammaTC,betaTC)),1,
				function(x) length(which(x)))

			## in case all genes are resolved or failed break the for cycle
			tab_N = table(N_na_genes)
			solved_genes[iter] = sum(tab_N[names(tab_N) %in% as.character(c(0,N*3))])
			message(paste0('Iteration ', iter, ' / ', n_iter,', solved ', round(solved_genes[iter]/gene_number*100), '% genes.'))
			if( solved_genes[iter] == gene_number ) break
			if( iter>1 ) if( solved_genes[iter] == solved_genes[iter-1] ) break

		}

	### Here I introduce the code to estimate a confidence interval for the rates first guess 

	resolved_id <- which(N_na_genes==0)
	# resolvedGenesNames <- geneNames[resolved_id]

	message("Confidece intervals estimation.")

	non_resolved_gene <- list(
			k1 = cbind(left=rep(NaN, N), opt=rep(NaN, N), right=rep(NaN, N), constant=rep(NaN, N)),
			k2 = cbind(left=rep(NaN, N), opt=rep(NaN, N), right=rep(NaN, N), constant=rep(NaN, N)),
			k3 = cbind(left=rep(NaN, N), opt=rep(NaN, N), right=rep(NaN, N), constant=rep(NaN, N))
			)
	confidence_intervals <- as.list(rep(NA,gene_number))
	for( g in which(N_na_genes>0) ) confidence_intervals[[g]] <- non_resolved_gene

	confidence_intervals[resolved_id] <- bplapply(resolved_id, function(g)
	{
		parameters <- c(alphaTC[g,],gammaTC[g,],betaTC[g,])
		allratesconfidences = sapply(seq_along(parameters), function(parname) 
		{
			par <- parameters[parname]
			suppressWarnings(capture.output(mOut <- list(
				left_1 =tryCatch(multiroot(f = logLikelihoodCIerror4sU_NF, start = 1e-2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalK1 = synthesis[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceK1 = synthesisVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				left_2 =tryCatch(multiroot(f = logLikelihoodCIerror4sU_NF, start = 1/2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalK1 = synthesis[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceK1 = synthesisVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				center =tryCatch(multiroot(f = logLikelihoodCIerror4sU_NF, start = par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalK1 = synthesis[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceK1 = synthesisVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				right_1 =tryCatch(multiroot(f = logLikelihoodCIerror4sU_NF, start = 2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalK1 = synthesis[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceK1 = synthesisVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				right_2 =tryCatch(multiroot(f = logLikelihoodCIerror4sU_NF, start = 1e2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalK1 = synthesis[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceK1 = synthesisVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN))
			)))
			precis = sapply(mOut, '[[', 'f.root')

			if( length(which(precis<1e-2))>0 )  {
				conf_int = sapply(mOut[which(precis<1e-2)], '[[', 'root')
				low_int = min(conf_int)
				high_int = max(conf_int)

				left = ifelse( low_int < par, low_int, NA)
				right = ifelse( high_int > par, high_int, NA)

				return(c(left = left, right = right))

			} else {

				return(c(left = NA, right = NA))

			}

		})

		# reintroduce NA to obtain vectors of the correct length

		k1left = allratesconfidences[1,1:N]
		k1right = allratesconfidences[2,1:N]
		k2left = allratesconfidences[1,(N+1):(2*N)]
		k2right = allratesconfidences[2,(N+1):(2*N)]
		k3left = allratesconfidences[1,(2*N+1):(3*N)]
		k3right = allratesconfidences[2,(2*N+1):(3*N)]

		rate_conf_int_k1 <- rate_conf_int_impute(cbind(left=k1left, opt=alphaTC[g,], right=k1right))
		k_start <- mean(rate_conf_int_k1[,2],na.rm=TRUE)
		if(is.finite(k_start)) {
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int_k1)$par
		} else {
			k_scores_out <- NaN
		}
		rate_conf_int_k1 <- cbind(rate_conf_int_k1, 'constant'=k_scores_out)

		rate_conf_int_k2 <- rate_conf_int_impute(cbind(left=k2left, opt=gammaTC[g,], right=k2right))
		k_start <- mean(rate_conf_int_k2[,2],na.rm=TRUE)
		if(is.finite(k_start)) {
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int_k2)$par
		} else {
			k_scores_out <- NaN
		}
		rate_conf_int_k2 <- cbind(rate_conf_int_k2, 'constant'=k_scores_out)

		rate_conf_int_k3 <- rate_conf_int_impute(cbind(left=k3left, opt=betaTC[g,], right=k3right))
		k_start <- mean(rate_conf_int_k3[,2],na.rm=TRUE)
		if(is.finite(k_start)) {
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int_k3)$par
		} else {
			k_scores_out <- NaN
		}
		rate_conf_int_k3 <- cbind(rate_conf_int_k3, 'constant'=k_scores_out)
				
		return(list(
			k1 = rate_conf_int_k1,
			k2 = rate_conf_int_k2,
			k3 = rate_conf_int_k3
			))

	}, BPPARAM=BPPARAM)

	return(confidence_intervals)

}

classify_rates_from_priors_v5 <- function(
	tpts,
	premature,
	total,
	prematureVariance,
	totalVariance,
	eiGenes=NULL,
	genesFilter=TRUE,
	llConfidenceThreshold=qchisq(.95,1),
	Dmin=1e-06,
	Dmax=10,
	BPPARAM=bpparam()
	) 
{
	mature <- total - premature
	matureVariance <- totalVariance + prematureVariance

	prOut <- estimate_priors(
		tpts,
		premature,
		mature,
		prematureVariance,
		matureVariance,
		eiGenes,
		genesFilter,
		Dmin,
		Dmax
		)

	alphaTC <- prOut[[1]]
	gammaTC <- prOut[[2]]
	betaTC  <- prOut[[3]]

	message("Rates optimization through the minimum likelihood.")

	N <- length(tpts)
	gene_number <- nrow(alphaTC)

	alphaTC[alphaTC<0] <- NA
	gammaTC[gammaTC<0] <- NA
	betaTC[betaTC<0]   <- NA

	N_na_genes = apply(is.na(cbind(alphaTC,gammaTC,betaTC)),1,function(x) length(which(x)))

	# initalize the list of results
	newParams <- as.list(rep(NA,gene_number))
	
	# model genes without NA values
	no_na_genes = which(N_na_genes==0)
	newParams[no_na_genes] <- bplapply(no_na_genes,function(r)
	{
		N <- length(tpts)
		parameters <- unname(c(lin2logPar(alphaTC[r,]),lin2logPar(gammaTC[r,]),lin2logPar(betaTC[r,])))
		tryCatch({
			optOut <- optim(par = parameters
				, fn = piecedFunctionLogLikelihood
				, tpts = tpts
				, experimentalP = premature[r,]
				, experimentalM = mature[r,]
				, varianceP = prematureVariance[r,]
				, varianceM = matureVariance[r,]
				, N = length(tpts)
				, control = list(maxit = 2000, fnscale=-1)
				, method = "BFGS")
			optOut$par <- c(log2linPar(optOut$par[1:N]), log2linPar(optOut$par[(N+1):(2*N)]), log2linPar(optOut$par[(2*N+1):(3*N)]))
			return(optOut)
		},error=function(e)return(errorOptim(N*3,e)))
	}, BPPARAM=BPPARAM)

	n_iter = length(tpts)
	solved_genes = numeric(n_iter)

	for( iter in 1:n_iter ) {

		many_na_genes <- which(N_na_genes>0 & N_na_genes<3*N)
		newParams[many_na_genes] <- bplapply(many_na_genes,function(r)
		{
			parameters <- c(lin2logPar(alphaTC[r,]),lin2logPar(gammaTC[r,]),lin2logPar(betaTC[r,]))
			parameters <- na.omit(parameters)
			na_ids <- as.vector(attr(parameters, 'na.action'))
			M <- N_na_genes[r]
			tryCatch({
				optOut <- optim(par = as.vector(parameters)
					, fn = piecedFunctionLogLikelihoodManyNA
					, tpts = tpts
					, experimentalP = premature[r,]
					, experimentalM = mature[r,]
					, varianceP = prematureVariance[r,]
					, varianceM = matureVariance[r,]
					, na_ids = na_ids
					, N = N
					, M = M
					, control = list(maxit = 2000, fnscale=-1)
					, method = "BFGS")
				## expand parameters
				optOut$par = unlist(lapply(rates_from_napars(optOut$par, na_ids, N, M), log2linPar))
				return(optOut)
			},error=function(e)return(errorOptim(N*3,e)))
		}, BPPARAM=BPPARAM)

		alphaTC <- t(sapply(newParams,function(g)g[[1]][1:N]))
		gammaTC <- t(sapply(newParams,function(g)g[[1]][(N+1):(2*N)]))
		betaTC  <- t(sapply(newParams,function(g)g[[1]][(2*N+1):(3*N)]))

		alphaTC[alphaTC<0] <- NA
		gammaTC[gammaTC<0] <- NA
		betaTC[betaTC<0] <- NA

		N_na_genes = apply(is.na(cbind(alphaTC,gammaTC,betaTC)),1,
			function(x) length(which(x)))

		## in case all genes are resolved or failed break the for cycle
		tab_N = table(N_na_genes)
		solved_genes[iter] = sum(tab_N[names(tab_N) %in% as.character(c(0,N*3))])
		message(paste0('Iteration ', iter, ' / ', n_iter,', solved ', round(solved_genes[iter]/gene_number*100), '% genes.'))
		if( solved_genes[iter] == gene_number ) {
			# message('Done.')
			break
		}
		if( iter>1 ) if( solved_genes[iter] == solved_genes[iter-1] ) {
			# message('Done.')
			break
		}

	}

	### Here I introduce the code to estimate a confidence interval for the rates first guess 

	resolved_id <- which(N_na_genes==0)
	# resolvedGenesNames <- geneNames[resolved_id]

	k1TC <- alphaTC
	k2TC <- gammaTC
	k3TC <- betaTC

	message("Confidece intervals estimation.")

	non_resolved_gene <- list(
			k1 = cbind(left=rep(NaN, N), opt=rep(NaN, N), right=rep(NaN, N), constant=rep(NaN, N)),
			k2 = cbind(left=rep(NaN, N), opt=rep(NaN, N), right=rep(NaN, N), constant=rep(NaN, N)),
			k3 = cbind(left=rep(NaN, N), opt=rep(NaN, N), right=rep(NaN, N), constant=rep(NaN, N))
			)
	confidence_intervals <- as.list(rep(NA,gene_number))
	for( g in which(N_na_genes>0) ) confidence_intervals[[g]] <- non_resolved_gene

	confidence_intervals[resolved_id] <- bplapply(resolved_id, function(g)
	{
#		parameters <- c(k1TC[g,],k2TC[g,],k3TC[g,])
		parameters <- unname(c(lin2logPar(k1TC[g,]),lin2logPar(k2TC[g,]),lin2logPar(k3TC[g,])))
		allratesconfidences = sapply(seq_along(parameters), function(parname) 
		{
			par <- parameters[parname]
			suppressWarnings(capture.output(mOut <- list(
				left_1 =tryCatch(multiroot(f = logLikelihoodCIerror_NF, start = 1e-2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				left_2 =tryCatch(multiroot(f = logLikelihoodCIerror_NF, start = 1/2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				center =tryCatch(multiroot(f = logLikelihoodCIerror_NF, start = par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				right_1 =tryCatch(multiroot(f = logLikelihoodCIerror_NF, start = 2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN)),
				right_2 =tryCatch(multiroot(f = logLikelihoodCIerror_NF, start = 1e2*par, name = parname, parameters = parameters, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], N=N, confidenceThreshold = llConfidenceThreshold), error=function(e) list(root=NaN, 'f.root'=NaN))
			)))
			precis = sapply(mOut, '[[', 'f.root')

			if( length(which(precis<1e-2))>0 )  {
				conf_int = sapply(mOut[which(precis<1e-2)], '[[', 'root')
				low_int = min(conf_int)
				high_int = max(conf_int)

				left = ifelse( low_int < par, low_int, NA)
				right = ifelse( high_int > par, high_int, NA)

				return(c(left = left, right = right))

			} else {

				return(c(left = NA, right = NA))

			}

		})

		# reintroduce NA to obtain vectors of the correct length

		k1left = log2linPar(allratesconfidences[1,1:N])
		k1right = log2linPar(allratesconfidences[2,1:N])
		k2left = log2linPar(allratesconfidences[1,(N+1):(2*N)])
		k2right = log2linPar(allratesconfidences[2,(N+1):(2*N)])
		k3left = log2linPar(allratesconfidences[1,(2*N+1):(3*N)])
		k3right = log2linPar(allratesconfidences[2,(2*N+1):(3*N)])

		# rate_conf_int_k1 <- rate_conf_int_impute(cbind(left=k1left, opt=alphaTC[g,], right=k1right))
		# rate_conf_int_k1 <- cbind(rate_conf_int_k1, 'constant'= score_and_par(rate_conf_int_k1)$par)

		# rate_conf_int_k2 <- rate_conf_int_impute(cbind(left=k2left, opt=gammaTC[g,], right=k2right))
		# rate_conf_int_k2 <- cbind(rate_conf_int_k2, 'constant'= score_and_par(rate_conf_int_k2)$par)

		# rate_conf_int_k3 <- rate_conf_int_impute(cbind(left=k3left, opt=betaTC[g,], right=k3right))
		# rate_conf_int_k3 <- cbind(rate_conf_int_k3, 'constant'= score_and_par(rate_conf_int_k3)$par)
		
		# return(list(
		# 	k1 = rate_conf_int_k1,
		# 	k2 = rate_conf_int_k2,
		# 	k3 = rate_conf_int_k3
		# 	))

		rate_conf_int_k1 <- rate_conf_int_impute(cbind(left=k1left, opt=alphaTC[g,], right=k1right))
		k_start <- mean(rate_conf_int_k1[,2],na.rm=TRUE)
		if(is.finite(k_start)) {
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int_k1)$par
		} else {
			k_scores_out <- NaN
		}
		rate_conf_int_k1 <- cbind(rate_conf_int_k1, 'constant'=k_scores_out)

		rate_conf_int_k2 <- rate_conf_int_impute(cbind(left=k2left, opt=gammaTC[g,], right=k2right))
		k_start <- mean(rate_conf_int_k2[,2],na.rm=TRUE)
		if(is.finite(k_start)) {
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int_k2)$par
		} else {
			k_scores_out <- NaN
		}
		rate_conf_int_k2 <- cbind(rate_conf_int_k2, 'constant'=k_scores_out)

		rate_conf_int_k3 <- rate_conf_int_impute(cbind(left=k3left, opt=betaTC[g,], right=k3right))
		k_start <- mean(rate_conf_int_k3[,2],na.rm=TRUE)
		if(is.finite(k_start)) {
			k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int_k3)$par
		} else {
			k_scores_out <- NaN
		}
		rate_conf_int_k3 <- cbind(rate_conf_int_k3, 'constant'=k_scores_out)
				
		return(list(
			k1 = rate_conf_int_k1,
			k2 = rate_conf_int_k2,
			k3 = rate_conf_int_k3
			))

	}, BPPARAM=BPPARAM)

	return(confidence_intervals)

}
