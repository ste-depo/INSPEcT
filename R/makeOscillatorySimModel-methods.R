#' @rdname makeOscillatorySimModel
#'
#' @description
#' This method allow the creation of synthesis, degradation and processing rates that generate
#' an oscillatory expression with a period of 24 hours. Two modes are available: one where 
#' oscillations arise just by oscillations in the synthesis of the genes (oscillatoryk3=FALSE, default) and
#' another one where both synthesis and degradation rates oscillates (oscillatoryk3=TRUE). In this latter case,
#' the oscillations of the two rates can be coupled by a cetrain delay (parametrer k3delay).
#' After the creation of the synthetic rates, a dataset with noise and contamination added can be made by \code{\link{makeSimDataset}}.
#' @param object An object of class INSPEcT
#' @param nGenes A numeric with the number of synthtic genes to be created
#' @param oscillatoryk3 A logical that enables also degradation rate to oscillate
#' @param k3delay A numeric that set the delay between synthesis and degradation oscillations. When NULL, no coupling between the two oscillations is set.
#' @param na.rm A logical that set whether missing values in the real dataset should be removed
#' @param seed A numeric to obtain reproducible results
#' @return An object of class INSPEcT_model with synthetic rates
#' @seealso \code{\link{makeSimModel}}
#' @examples
#' nascentInspObj <- readRDS(system.file(package='INSPEcT', 'nascentInspObj.rds'))
#' simRates<-makeOscillatorySimModel(nascentInspObj, 1000, seed=1)
#' table(geneClass(simRates))
setMethod('makeOscillatorySimModel', 'INSPEcT', function(object
									, nGenes
									, oscillatoryk3=FALSE
									, k3delay=NULL
									, na.rm=TRUE
									, seed=NULL)
{

	tpts <- object@tpts
	if( !is.numeric(tpts) ) stop('makeSimModel: simulated data can be created only from time-course.')

	#I remove genes without intronic signal
	genesTmp <- which(apply(ratesFirstGuess(object, 'preMRNA'),1,function(r)all(is.finite(r))&all(r>0)))		

	concentrations <- list(
		total=ratesFirstGuess(object, 'total')[genesTmp,]
		, total_var=ratesFirstGuessVar(object, 'total')[genesTmp,]
		, preMRNA=ratesFirstGuess(object, 'preMRNA')[genesTmp,]
		, preMRNA_var=ratesFirstGuessVar(object, 'preMRNA')[genesTmp,]
		)
	rates <- list(
		alpha=ratesFirstGuess(object, 'synthesis')[genesTmp,]
		, alpha_var=ratesFirstGuessVar(object, 'synthesis')[genesTmp,]
		, beta=ratesFirstGuess(object, 'degradation')[genesTmp,]
		, gamma=ratesFirstGuess(object, 'processing')[genesTmp,]
		)

	if( !is.null(seed) ) set.seed(seed)

	# read input
	alpha   <- rates$alpha
	beta    <- rates$beta
	gamma   <- rates$gamma
	total   <- concentrations$total
	preMRNA <- concentrations$preMRNA

	alpha_var <- rates$alpha_var

	total_var   <- concentrations$total_var
	preMRNA_var <- concentrations$preMRNA_var

	alphaFitVariance <- lm(formula = log(c(sqrt(alpha_var))) ~ log(c(alpha)))$coefficients
	alphaFitVarianceLaw <- function(alpha)(exp(alphaFitVariance[[1]])*alpha^(alphaFitVariance[[2]]))^2
	
	totalFitVariance <- lm(formula = log(c(sqrt(total_var))) ~ log(c(total)))$coefficients
	totalFitVarianceLaw <- function(total)(exp(totalFitVariance[[1]])*total^(totalFitVariance[[2]]))^2

	preFitVariance <- lm(formula = log(c(sqrt(preMRNA_var))) ~ log(c(preMRNA)))$coefficients
	preFitVarianceLaw <- function(pre)(exp(preFitVariance[[1]])*pre^(preFitVariance[[2]]))^2

	# sample initial timepoint
	message('sampling means from rates distribution...')
	#
	if( na.rm == TRUE ) {
		beta[beta <= 0] <- NA
		gamma[gamma <= 0] <- NA
	}
	alphaVals <- sample(alpha, nGenes, replace=TRUE)
	betaVals  <- 2^sampleNormQuantile(
		values_subject=log2(alphaVals)
		, dist_subject=log2(alpha)
		, dist_object=log2(beta)
		, na.rm=na.rm)
	gammaVals <- 2^sampleNormQuantile(
		values_subject=log2(betaVals)
		, dist_subject=log2(beta)
		, dist_object=log2(gamma)
		, na.rm=na.rm)

	if( !oscillatoryk3 ) {

		alphaParams <- generateOscillatoryParams(tpts, alphaVals, probs=c(oscillatory=.5))
		betaParams  <- generateOscillatoryParams(tpts, betaVals , probs=c(oscillatory=.0))
		gammaParams <- generateOscillatoryParams(tpts, gammaVals, probs=c(oscillatory=.0))

	} else if( oscillatoryk3 & !is.null(k3delay) ) {

		alphaParams <- generateOscillatoryParams(tpts, alphaVals, probs=c(oscillatory=.5))

		osc_ix <- which(sapply(alphaParams, '[[', 'type') == 'oscillatory')
		amplitude <- sapply(alphaParams[osc_ix], '[[', 'params')[3,]
		alpha_shift <- sapply(alphaParams[osc_ix], '[[', 'params')[4,]
		
		betaParams  <- generateParamsShift(tpts, betaVals, osc_ix, amplitude, alpha_shift+k3delay)
		gammaParams <- generateOscillatoryParams(tpts, gammaVals, probs=c(oscillatory=.0))

	} else if( oscillatoryk3 & is.null(k3delay) ) {

		alphaParams <- generateOscillatoryParams(tpts, alphaVals, probs=c(oscillatory=.5))
		betaParams  <- generateOscillatoryParams(tpts, betaVals, probs=c(oscillatory=.5))
		gammaParams <- generateOscillatoryParams(tpts, gammaVals, probs=c(oscillatory=.0))
		

	} 

	paramSpecs <- lapply(1:nGenes, 
		function(i) 
			list(alpha=alphaParams[[i]]
			   , beta=betaParams[[i]]
			   , gamma=gammaParams[[i]]))

	out <- lapply(1:nGenes, function(i){
			tryCatch(
				.makeModel(tpts, paramSpecs[[i]]),error=function(e){cbind(time=rep(NaN,length(tpts))
																		  ,preMRNA=rep(NaN,length(tpts))
																		  ,total=rep(NaN,length(tpts))
																		  ,alpha=rep(NaN,length(tpts))
																		  ,beta=rep(NaN,length(tpts))
																		  ,gamma=rep(NaN,length(tpts)))})})

	okGenes <- which(sapply(out,function(i)all(is.finite(unlist(i)))))
	out <- out[okGenes]
	paramSpecs <- paramSpecs[okGenes]

	cleanDataSet <- list(
		tpts = tpts
		, concentrations = list(
			total=t(sapply(out, function(x) x$total))
			, total_var=rep(1,nGenes)
			, preMRNA=t(sapply(out, function(x) x$preMRNA))
			, preMRNA_var=rep(1,nGenes)
			)
		, rates = list(
			alpha=t(sapply(out, function(x) x$alpha))
			, alpha_var=rep(1,nGenes)
			, beta=t(sapply(out, function(x) x$beta))
			, gamma=t(sapply(out, function(x) x$gamma))
			)
		)

	alphaSim_noisevar <- t(apply(cleanDataSet$rates$alpha,1,function(x)alphaFitVarianceLaw(x)))
	totalSim_noisevar <- t(apply(cleanDataSet$concentrations$total,1,function(x)totalFitVarianceLaw(x)))
	preSim_noisevar <- t(apply(cleanDataSet$concentrations$preMRNA,1,function(x)preFitVarianceLaw(x)))

	# select genes whose noise evaluation succeded
	okGenes <- which(
		apply(alphaSim_noisevar,1,function(r)all(is.finite(r))) &
		apply(totalSim_noisevar,1,function(r)all(is.finite(r))) &
		apply(preSim_noisevar,1,function(r)all(is.finite(r))) 
		)

	paramSpecs <- paramSpecs[okGenes]
	alphaSim_noisevar <- alphaSim_noisevar[okGenes,]
	totalSim_noisevar <- totalSim_noisevar[okGenes,]
	preSim_noisevar   <- preSim_noisevar[okGenes,]
	# add params specification
	simulatedFC <- list(
		alpha=apply(cleanDataSet$rates$alpha[okGenes, ], 1, 
			function(x) diff(log2(range(x))))
		, beta=apply(cleanDataSet$rates$beta[okGenes, ], 1, 
			function(x) diff(log2(range(x))))
		, gamma=apply(cleanDataSet$rates$gamma[okGenes, ], 1, 
			function(x) diff(log2(range(x))))
		)
	noiseVar <- list(
		alpha=alphaSim_noisevar
		, total=totalSim_noisevar
		, pre=preSim_noisevar
		)

	# arrange simdataSpecs form .makeSimData
	paramSpecs <- lapply(paramSpecs, function(x) list(x))
	#
	newObject <- new('INSPEcT_model')
	newObject@ratesSpecs <- paramSpecs
	newObject@params$sim$flag <- TRUE
	newObject@params$sim$foldchange <- simulatedFC
	newObject@params$sim$noiseVar <- noiseVar
	newObject@params$sim$noiseFunctions <- list(alpha = alphaFitVarianceLaw, preMRNA = preFitVarianceLaw, total = totalFitVarianceLaw)
	newObject@params$tpts <- tpts

	if(length(out$paramSpecs)>nGenes){return(newObject[1:nGenes])} #Return only the required number of genes or less
	return(newObject)

})

generateOscillatoryParams <- function(tpts
						 , sampled_val
						 , probs=c(oscillatory=.5))
# given a vector of absolute values and a vector of log2foldchanges
# create parametric functions (either constant, sigmoidal or impulse, 
# according to probs) and evaluate them at tpts.
{

	nGenes <- length(sampled_val)
	# 
	n_constant <- round(nGenes * (1 - probs['oscillatory']))

	# initialize
	params <- as.list(rep(NA,nGenes))

	# constant: choose the one with the lower absoulute fold change to be 
	# constant
	constant_idx <- 1:nGenes %in% sample(1:nGenes, n_constant)
	if( any(constant_idx) )
	{
		params[constant_idx] <- lapply(sampled_val[constant_idx], 
			function(val) 
				list(type='constant', fun=constantModelP , params=val, df=1)
				)
	}

	generateOscillatoryParamsFun <- function(tpts, sampled_val)
	# Given an absolute value and a value of log2fold change sample a set 
	# of parameters for the impulse function.
	{

		n <- length(sampled_val)

		sampled_freqs <- pi/12
		initial_values <- sampled_val
		amplitude <- runif(n, .05, .4)
		x_shift <- runif(n, min(tpts), max(tpts))

		oscillatorypars <- cbind(
			sampled_freqs
			, initial_values
			, amplitude
			, x_shift
			)

		return(oscillatorypars)

	}

	# sigmoid
	oscillatory_idx <- !constant_idx
	if( any(oscillatory_idx) ) {
		params[oscillatory_idx] <- lapply(
			which(oscillatory_idx)
			, function(i) list(
				type='oscillatory'
				, fun=oscillatoryModelP
				, params=generateOscillatoryParamsFun(
					tpts
					, sampled_val[i] 
					)
				, df=4
				)			
			)
	}

	return(params)

}

generateParamsShift <- function(tpts
						 , sampled_val
						 , oscillatory_idx
						 , amplitude
						 , x_shift)
# given a vector of absolute values and a vector of log2foldchanges
# create parametric functions (either constant, sigmoidal or impulse, 
# according to probs) and evaluate them at tpts.
{

	nGenes <- length(sampled_val)

	# initialize
	params <- as.list(rep(NA,nGenes))

	# constant: choose the one with the lower absoulute fold change to be 
	# constant
	constant_idx <- ! (1:nGenes %in% oscillatory_idx)
	if( any(constant_idx) )
	{
		params[constant_idx] <- lapply(sampled_val[constant_idx], 
			function(val) 
				list(type='constant', fun=constantModelP , params=val, df=1)
				)
	}

	generateOscillatoryParamsShiftFun <- function(tpts, sampled_val, amplitude, x_shift)
	# Given an absolute value and a value of log2fold change sample a set 
	# of parameters for the impulse function.
	{

		n <- length(sampled_val)

		sampled_freqs <- pi/12
		initial_values <- sampled_val
		# amplitude <- runif(n, .05, .4)

		oscillatorypars <- cbind(
			sampled_freqs
			, initial_values
			, amplitude
			, x_shift
			)

		return(oscillatorypars)

	}
	# sigmoid
	if( any(oscillatory_idx) ) {
		params[oscillatory_idx] <- lapply(
			seq_along(oscillatory_idx)
			, function(i) list(
				type='oscillatory'
				, fun=oscillatoryModelP
				, params=generateOscillatoryParamsShiftFun(
					tpts
					, sampled_val[oscillatory_idx[i]]
					, amplitude[i]
					, x_shift[i]
					)
				, df=4
				)			
			)
	}

	return(params)

}