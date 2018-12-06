find_tt_par <- function(tpts)
{
  cvLogTpts <- function(a , tpts)
  {
  	newtime <- log2(tpts + a )
    stats::sd(diff(newtime)) / mean(diff(newtime))
  }
  if(length(tpts)>2){return(optimize(f=cvLogTpts, interval=c(0,5), tpts=tpts )$minimum)}
  else{return(1)}
}

time_transf <- function(t, log_shift, c = NaN) 
{
	t[ t <= (-log_shift) ] <- NaN
	newtime <- log2(t+log_shift)
	return(newtime)
} 

time_transf_inv <- function(t, log_shift) 2^t - log_shift

time_transf_NoNascent <- function(t, log_shift, c) 
{
	newtime <- log2(t+log_shift) + c
	return(newtime)
}

chisqFunction <- function(experiment, model, variance=NULL)
{
	if( is.null(variance)) variance <- stats::var(experiment)
	sum((experiment - model )^2/variance )
}

logLikelihoodFunction <- function(experiment, model, variance=NULL)
{
    if( is.null(variance)) variance <- stats::var(experiment)
    sum(log(2*pnorm(-abs(experiment-model),mean=0,sd=sqrt(variance))))
}

.emptyGene <- function(error='')
{
	emptyRate <- function() return(list(fun=NA, type=NA, df=0, params=NaN))
	return(
		list(alpha=emptyRate(), beta=emptyRate(), gamma=emptyRate()
		, test=NaN, logLik=NaN, AIC=NaN, AICc=NaN, counts=NaN, convergence=1, message=error)
		)		
}

.makeEmptyModel <- function(tpts) {
	model <- matrix(NA, nrow=length(tpts), 5)
	colnames(model) <- c('alpha','beta','gamma','preMRNA','total')
	as.data.frame(model)
}

.makeModel <- function(tpts, hyp, log_shift, time_transf, ode, .rxnrate, c= NaN)
{
	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(time_transf(x, log_shift, c), hyp$alpha$par)
	params$beta  <- function(x) 
		hyp$beta$fun$value(time_transf(x, log_shift, c), hyp$beta$par)
	params$gamma <- function(x) 
		hyp$gamma$fun$value(time_transf(x, log_shift, c), hyp$gamma$par)
	cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1]), 
		params$alpha(tpts[1]) / params$beta(tpts[1]) + 
			params$alpha(tpts[1]) / params$gamma(tpts[1]))
	names(cinit) <- c('p', 't')
	model <- as.data.frame(
		ode(y=cinit, times=tpts, func=.rxnrate, parms=params))
	model$alpha <- params$alpha(tpts)
	model$beta  <- params$beta(tpts)
	model$gamma <- params$gamma(tpts)
	colnames(model)[2:3] <- c('preMRNA','total')
	return(model)
}

.rxnrate <- function(t,c,parms){
 
	# rate constant passed through a list called parms
	alpha <- parms$alpha
	beta  <- parms$beta
	gamma <- parms$gamma

	# derivatives dc/dt are computed below
	r=rep(0,length(c))
	r[1] <- alpha(t) - gamma(t) * c["p"]
	r[2] <- alpha(t) - beta(t) * (c["t"] - c["p"] )

	# c is the concentration of species
	
	# the computed derivatives are returned as a list
	# order of derivatives needs to be the same as the order of species in c
	return(list(r))
 
}

.makeSimpleModel <- function(tpts, hyp, log_shift, time_transf, ode, .rxnrateSimple)
{
	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(time_transf(x, log_shift), hyp$alpha$params)
	params$beta  <- function(x)
		hyp$beta$fun$value(time_transf(x, log_shift), hyp$beta$params)
	#
	cinit <- c(t = params$alpha(tpts[1]) / params$beta(tpts[1]))
	names(cinit) <- 't'
	model <- as.data.frame(
		ode(y=cinit, times=tpts, func=.rxnrateSimple, parms=params))
	model$alpha <- params$alpha(tpts)
	model$beta  <- params$beta(tpts)
	model$gamma   <- rep(NA, length(tpts))
	model$preMRNA <- rep(NA, length(tpts))
	colnames(model)[2] <- 'total'
	return(model)
}

.makeEmptySimpleModel <- function(tpts) {
	model <- matrix(NA, nrow=length(tpts), 3)
	colnames(model) <- c('alpha','beta','total')
	as.data.frame(model)
}

.rxnrateSimple <- function(t,c,parms){
 
	# rate constant passed through a list called parms
	alpha <- parms$alpha
	beta  <- parms$beta

	# derivatives dc/dt are computed below
	r=rep(0,length(c))
	r[1] <- alpha(t) - beta(t) * c["t"]

	# c is the concentration of species
	
	# the computed derivatives are returned as a list
	# order of derivatives needs to be the same as the order of species in c
	return(list(r))
 
}

.inspect.engine <- function(tpts, log_shift, concentrations, rates
	, nInit=10, nIter=300, na.rm=TRUE, BPPARAM=bpparam() #nCores=2L
	, verbose=TRUE, estimateRatesWith=c('der', 'int'), nAttempts=1
	, sigmoidDegradation=FALSE, sigmoidSynthesis=FALSE, sigmoidTotal=FALSE
	, sigmoidProcessing=FALSE, sigmoidPre=FALSE
	, testOnSmooth=TRUE, seed=NULL)
{

	chisq.test.inspect <- function(D, df)
		pchisq(D, df, lower.tail=TRUE)

	fisher.test.inspect <- function(p1, p2)
		pchisq(-2*(log(p1) + log(p2)), 4, lower.tail=FALSE)

	optimParamsSimple <- function(interpRates, tpts_exp, alpha_exp, alpha_var
		, total_exp, total_var, maxit=500
		, log_shift, time_transf, .rxnrateSimple, ode, .makeSimpleModel, logLikelihoodFunction
		, .emptyGene)
	{

		if( is.null(interpRates) ) return(.emptyGene())

		simpleModelChisq <- function(par, tpts, fun, df, alpha_exp, alpha_var
			, total_exp, total_var)
		{
			splitpar <- split(par 
				, c(rep('alpha',df[1]), rep('beta',df[2]))
				)
			#
			params <- list()
			params$alpha <- function(x) 
				fun$alpha$value(time_transf(x, log_shift), splitpar$alpha)
			params$beta <- function(x)
				fun$beta$value(time_transf(x, log_shift), splitpar$beta)
			#
			cinit <- c(params$alpha(tpts[1]) / params$beta(tpts[1]))
			names(cinit) <- 't'
			model <- ode(y=cinit, times=tpts, func=.rxnrateSimple, parms=params)
			#
			alpha_model <- params$alpha(tpts)
			total_model <- model[,'t']
			#
			D <- chisqFunction(alpha_exp, alpha_model, alpha_var) +
				chisqFunction(total_exp, total_model, total_var)
			return(D)
		}

		optOut <- tryCatch(
			optim(
				par=c(interpRates$alpha$par, interpRates$beta$par)
				, fn=simpleModelChisq, tpts=tpts_exp 
				, fun=list(alpha=interpRates$alpha$fun
					, beta=interpRates$beta$fun) 
				, df=c(interpRates$alpha$df, interpRates$beta$df)
				, alpha_exp=alpha_exp, total_exp=total_exp
				, alpha_var=alpha_var, total_var=total_var
				, control=list(maxit=maxit) #, trace=1)
				, method='Nelder-Mead'
				)
			, error=function(e)
				optim(
					par=c(interpRates$alpha$par, interpRates$beta$par)
					, fn=simpleModelChisq, tpts=tpts_exp 
					, fun=list(alpha=interpRates$alpha$fun
						, beta=interpRates$beta$fun) 
					, df=c(interpRates$alpha$df, interpRates$beta$df)
					, alpha_exp=alpha_exp, total_exp=total_exp
					, alpha_var=alpha_var, total_var=total_var
					, control=list(maxit=maxit) #, trace=1)
					, method='BFGS'
					)
				)
		#
		splitpar <- split(optOut$par, c(rep('alpha',interpRates$alpha$df)
			, rep('beta',interpRates$beta$df)))
		interpRates$alpha$params <- splitpar$alpha
		interpRates$beta$params  <- splitpar$beta
		#
		model <- .makeSimpleModel(tpts=tpts_exp
			, hyp=list(alpha=interpRates$alpha, beta=interpRates$beta)
			, log_shift, time_transf, ode, .rxnrateSimple)
		logLik <- logLikelihoodFunction(alpha_exp, model$alpha, alpha_var) + 
			logLikelihoodFunction(total_exp, model$total, total_var)
		k <- interpRates$alpha$df + interpRates$beta$df
		n <- length(alpha_exp) + length(total_exp)
		chisqTest <- log(pchisq(optOut$value, n-k, lower.tail=TRUE))
		AIC <- 2*k - 2*logLik
		AICc <- AIC + 2*k*(k+1)/(n-k-1)
		return(list(
			alpha=interpRates$alpha
			, beta=interpRates$beta
			, gamma=.emptyGene()$gamma
			, test=chisqTest
			, logLik=logLik
			, AIC=AIC
			, AICc=AICc
			, counts=optOut$counts
			, convergence=optOut$convergence
			, message=optOut$message
			))

	}

	optimParams <- function(interpRates, tpts_exp, alpha_exp, alpha_var, total_exp
		, total_var, preMRNA_exp, preMRNA_var
		# , test=c('merged', 'combined'), pval=c('lin', 'log')
		, maxit=500, log_shift, time_transf, .rxnrate, ode, .makeModel, logLikelihoodFunction)
	{
		modelChisq <- function(par, tpts, fun, df, alpha_exp, alpha_var #, pval
			, total_exp, total_var, preMRNA_exp, preMRNA_var)
		{
			splitpar <- split(par 
				, c(rep('alpha',df[1]), rep('beta',df[2]), rep('gamma',df[3])) 
				)
			#
			params <- list()
			params$alpha <- function(x) 
				fun$alpha$value(time_transf(x, log_shift), splitpar$alpha)
			params$beta  <- function(x)
				fun$beta$value(time_transf(x, log_shift), splitpar$beta)
			params$gamma <- function(x)
				fun$gamma$value(time_transf(x, log_shift), splitpar$gamma)
			#
			cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1])
				, params$alpha(tpts[1]) / params$beta(tpts[1]) + 
					params$alpha(tpts[1]) / params$gamma(tpts[1]))
			names(cinit) <- c('p', 't')
			model <- ode(y=cinit, times=tpts, func=.rxnrate, parms=params)
			#
			alpha_model <- params$alpha(tpts)
			total_model <- model[,'t']
			preMRNA_model <- model[,'p']
			#
			D <- chisqFunction(alpha_exp, alpha_model, alpha_var) +
				chisqFunction(total_exp, total_model, total_var) +
				chisqFunction(preMRNA_exp, preMRNA_model, preMRNA_var)
			df <- length(alpha_exp) + length(total_exp) + length(preMRNA_exp) - sum(df)
			testValue <- chisq.test.inspect(D, df)
			return(testValue)
		}
		optOut <- tryCatch(
			#
			# optimize the model to the experiment using 
			# Nelder-Mead method
			#
			optim(
				par=unlist(sapply(interpRates, '[[', 'params'))
				, fn=modelChisq , tpts=tpts_exp #, pval=pval 
				, fun=sapply(interpRates, '[[', 'fun')
				, df=unlist(sapply(interpRates, '[[', 'df'))
				, alpha_exp=alpha_exp, total_exp=total_exp
				, preMRNA_exp=preMRNA_exp
				, alpha_var=alpha_var, total_var=total_var
				, preMRNA_var=preMRNA_var
				, control = list(maxit = maxit)
				, method = 'Nelder-Mead'
				)
			#
			# in case Nelder-Mead method fails, try with BFGS
			#
			, error=function(e) {
				optim(
					par=unlist(sapply(interpRates, '[[', 'params'))
					, fn=modelChisq, tpts=tpts_exp #, pval=pval
					, fun=sapply(interpRates, '[[', 'fun')
					, df=unlist(sapply(interpRates, '[[', 'df'))
					, alpha_exp=alpha_exp, total_exp=total_exp
					, preMRNA_exp=preMRNA_exp
					, alpha_var=alpha_var, total_var=total_var
					, preMRNA_var=preMRNA_var
					, control = list(maxit = maxit)
					, method = 'BFGS'
					)
			})
		if( !is.null(optOut$error) ) return(optOut$error)
		#
		# assign the parameters belonging to the optimized model
		# to the parameter functions
		#
		splitpar <- split(optOut$par
			, c(rep('alpha',interpRates$alpha$df)
				, rep('beta',interpRates$beta$df)
				, rep('gamma',interpRates$gamma$df))
			)
		interpRates$alpha$params <- splitpar$alpha
		interpRates$beta$params  <- splitpar$beta
		interpRates$gamma$params <- splitpar$gamma
		# even if the minimization used the linear pvalue
		# give back the log one
		# if( pval=='lin' ) 
		optOut$value <- log(optOut$value)
		#
		# return parameter functions and other output
		# from the optimization procedure
		#
		model <- .makeModel(tpts=tpts_exp
			, hyp=list(alpha=interpRates$alpha, beta=interpRates$beta, 
				gamma=interpRates$gamma)
			, log_shift, time_transf, ode, .rxnrate)
		logLik <- logLikelihoodFunction(alpha_exp, model$alpha, alpha_var) + 
			logLikelihoodFunction(total_exp, model$total, total_var) +
			logLikelihoodFunction(preMRNA_exp, model$preMRNA, preMRNA_var)

		k <- interpRates$alpha$df + interpRates$beta$df + interpRates$gamma$df
		n <- length(alpha_exp) + length(total_exp) + length(preMRNA_exp)
		chisqTest <- optOut$value
		AIC <- 2*k - 2*logLik
		AICc <- AIC + 2*k*(k+1)/(n-k-1)
		return(list(
			alpha=interpRates$alpha
			, beta=interpRates$beta
			, gamma=interpRates$gamma
			, test=chisqTest
			, logLik=logLik
			, AIC=AIC
			, AICc=AICc
			, counts=optOut$counts
			, convergence=optOut$convergence
			, message=optOut$message
			))
	}
	
	modelOneGene <- function(i, seed=NULL,
			.chooseModel, time_transf, .DimpulseModel, .DsigmoidModel, constantModelP,
			.emptyGene, sigmoidModel, impulseModel, sigmoidModelP, impulseModelP,
			.polynomialModelP, .makeModel, .makeSimpleModel, logLikelihoodFunction, .rxnrate,
			.rxnrateSimple, optimParams, optimParamsSimple, verbose, nAttempts,
			concentrations, rates, tpts, log_shift, na.rm, sigmoidTotal,
			sigmoidSynthesis, nInit, nIter, testOnSmooth, estimateRatesWith) 
	{
		## set the mode of the gene, "only exons gene" or 
		## "introns exons gene"
		if( !is.null(seed) ) set.seed(seed)
		if( all(is.na(concentrations$preMRNA[i,])) |  
				all(is.na(rates$gamma[i,])) ) 
		{
			intExMode <- FALSE
		} else {
			intExMode <- TRUE
		}
		## start the analysis
		paramAttempts <- sapply(1:nAttempts, function(k)
			# tryCatch(
			{
				modelTotalRNAfun <- tryCatch(.chooseModel(tpts=tpts
						, log_shift=log_shift
						, experiment=concentrations$total[i,]
						, variance=concentrations$total_var[i]
						, na.rm=na.rm, sigmoid=sigmoidTotal
						, impulse=TRUE, polynomial=FALSE
						, nInit=nInit, nIter=nIter
						, time_transf=time_transf
						, sigmoidModel=sigmoidModel
						, impulseModel=impulseModel
						, sigmoidModelP=sigmoidModelP
						, impulseModelP=impulseModelP
						, .polynomialModelP=.polynomialModelP
						), error=function(e) return(.emptyGene(e)))
				modelSynthesisRatefun <- tryCatch(.chooseModel(tpts=tpts
						, log_shift=log_shift
						, experiment=rates$alpha[i,] 
						, variance=rates$alpha_var[i]
						, na.rm=na.rm, sigmoid=sigmoidSynthesis
						, impulse=TRUE, polynomial=FALSE
						, nInit=nInit, nIter=nIter
						, time_transf=time_transf
						, sigmoidModel=sigmoidModel
						, impulseModel=impulseModel
						, sigmoidModelP=sigmoidModelP
						, impulseModelP=impulseModelP
						, .polynomialModelP=.polynomialModelP
						), error=function(e) return(.emptyGene(e)))
				modelTotalRNA <- 
					if( testOnSmooth ) {
						modelTotalRNAfun$fun$value(
							time_transf(tpts, log_shift)
							, modelTotalRNAfun$params)
					} else { concentrations$total[i,] }
				modelSynthesisRate <- 
					if( testOnSmooth ) {
						modelSynthesisRatefun$fun$value(
							time_transf(tpts, log_shift)
							, modelSynthesisRatefun$params)	
					} else { rates$alpha[i,] }
				if( intExMode ) {
					modelPreMRNAfun <- tryCatch(.chooseModel(tpts=tpts
							, log_shift=log_shift
							, experiment=concentrations$preMRNA[i,]
							, variance=concentrations$preMRNA_var[i]
							, na.rm=na.rm, sigmoid=sigmoidPre
							, impulse=TRUE, polynomial=FALSE
							, nInit=nInit, nIter=nIter
							, time_transf=time_transf
							, sigmoidModel=sigmoidModel
							, impulseModel=impulseModel
							, sigmoidModelP=sigmoidModelP
							, impulseModelP=impulseModelP
							, .polynomialModelP=.polynomialModelP
							), error=function(e) return(.emptyGene(e)))
					modelPreMRNA <- 
						if( testOnSmooth ) 
							modelPreMRNAfun$fun$value(
								time_transf(tpts, log_shift)
								, modelPreMRNAfun$params) 
						else concentrations$preMRNA[i,]
				}
				#
				if( estimateRatesWith == 'der' ) {
					if( testOnSmooth ) {
						# total RNA derivative
						if( modelTotalRNAfun$type == 'impulse' )
							modelTotalRNAderivative <- .DimpulseModel(
								time_transf(tpts, log_shift)
								, modelTotalRNAfun$params)
						if( modelTotalRNAfun$type == 'sigmoid' )
							modelTotalRNAderivative <- .DsigmoidModel(
								time_transf(tpts, log_shift)
								, modelTotalRNAfun$params)
						if( modelTotalRNAfun$type == 'constant' )
							modelTotalRNAderivative <- rep(0, length(tpts))
						if( intExMode ) {
							# pre mRNA derivative
							if( modelPreMRNAfun$type == 'impulse' )
								modelPreMRNAderivative <- .DimpulseModel(
									time_transf(tpts, log_shift)
									, modelPreMRNAfun$params)
							if( modelPreMRNAfun$type == 'sigmoid' )
								modelPreMRNAderivative <- .DsigmoidModel(
									time_transf(tpts, log_shift)
									, modelPreMRNAfun$params)
							if( modelPreMRNAfun$type == 'constant' )
								modelPreMRNAderivative <- rep(0, length(tpts))
						}
					} else {
						# total RNA derivative
						spfun <- splinefun(tpts, concentrations$total[i,] )
						modelTotalRNAderivative <- spfun(tpts, deriv=1)
						if( intExMode ) {
							# pre mRNA derivative
							spfun <- splinefun(tpts, concentrations$preMRNA[i,] )
							modelPreMRNAderivative <- spfun(tpts, deriv=1)
						}
					}
					# degradation rate
					modelTotalRNAderivative[1] <- 0
					firstGuessDegrRate <- (modelSynthesisRate - 
						modelTotalRNAderivative ) / (modelTotalRNA - modelPreMRNA )
					if( intExMode ) {
						# processing rate
						modelPreMRNAderivative[1] <- 0
						firstGuessProcessingRate <- (modelSynthesisRate - 
							modelPreMRNAderivative ) / modelPreMRNA					
					}
				} else {
					# degradation rate
					firstGuessDegrRate <- rates$beta[i,]
					if( intExMode ) {
						# processing rate
						firstGuessProcessingRate <- rates$gamma[i,]
					}
				}
				idx <- firstGuessDegrRate<0 | 
					!is.finite(firstGuessDegrRate)
				firstGuessDegrRate[idx] <- NA
				if( intExMode ) {
					idx <- firstGuessProcessingRate<0 | 
						!is.finite(firstGuessProcessingRate)
					firstGuessProcessingRate[idx] <- NA
				}
				#
				constantRates <- list(
					alpha=list(fun=constantModelP
						, type='constant', df=1
						, params=mean(modelSynthesisRate, na.rm=TRUE))
					, beta=list(fun=constantModelP
						, type='constant', df=1
						, params=mean(firstGuessDegrRate, na.rm=TRUE))
					, gamma=if( intExMode ) {
						list(fun=constantModelP
							, type='constant', df=1
							, params=mean(firstGuessProcessingRate, na.rm=TRUE))
						} else { NULL }
					)
				varyingRates <- list(
						alpha=modelSynthesisRatefun
						, beta=tryCatch(
							.chooseModel(tpts=tpts
								, log_shift=log_shift
								, experiment=firstGuessDegrRate
								, variance=1
								, na.rm=na.rm, sigmoid=sigmoidDegradation
								, impulse=TRUE, polynomial=FALSE
								, nInit=nInit, nIter=nIter
								, time_transf=time_transf
								, sigmoidModel=sigmoidModel
								, impulseModel=impulseModel
								, sigmoidModelP=sigmoidModelP
								, impulseModelP=impulseModelP
								, .polynomialModelP=.polynomialModelP
								)
							, error=function(e) return(.emptyGene(e)$beta))
						, gamma=if( intExMode ) { tryCatch(
							.chooseModel(tpts=tpts
								, log_shift=log_shift
								, experiment=firstGuessProcessingRate
								, variance=1
								, na.rm=na.rm, sigmoid=sigmoidProcessing
								, impulse=TRUE, polynomial=FALSE
								, nInit=nInit, nIter=nIter
								, time_transf=time_transf
								, sigmoidModel=sigmoidModel
								, impulseModel=impulseModel
								, sigmoidModelP=sigmoidModelP
								, impulseModelP=impulseModelP
								, .polynomialModelP=.polynomialModelP
								), error=function(e) return(.emptyGene(e)$gamma))
							} else { NULL }
						)
				ratesToTest <- list('0'=constantRates
					, a=list(alpha=varyingRates$alpha, beta=constantRates$beta
						, gamma=constantRates$gamma)
					, b=list(alpha=constantRates$alpha, beta=varyingRates$beta
						, gamma=constantRates$gamma)
					, c=if( intExMode ) 
						list(alpha=constantRates$alpha, beta=constantRates$beta
						, gamma=varyingRates$gamma) else NULL
					, ab=list(alpha=varyingRates$alpha, beta=varyingRates$beta
						, gamma=constantRates$gamma)
					, bc=if( intExMode ) 
						list(alpha=constantRates$alpha, beta=varyingRates$beta
						, gamma=varyingRates$gamma) else NULL
					, ac=if( intExMode ) 
						list(alpha=varyingRates$alpha, beta=constantRates$beta
						, gamma=varyingRates$gamma) else NULL
					, abc=if( intExMode ) varyingRates else NULL
					)
				if( intExMode ) {
					results <- lapply(ratesToTest, function(interpRates) {
						tryCatch(
							optimParams(interpRates
								, tpts_exp=tpts
								, alpha_exp=modelSynthesisRate
								, alpha_var=rates$alpha_var[i]
								, total_exp=modelTotalRNA
								, total_var=concentrations$total_var[i]
								, preMRNA_exp=modelPreMRNA
								, preMRNA_var=concentrations$preMRNA_var[i]
								, maxit=nIter
								, log_shift=log_shift
								, ode=deSolve::ode
								, time_transf=time_transf
								, .rxnrate=.rxnrate
								, .makeModel=.makeModel
								, logLikelihoodFunction=logLikelihoodFunction
								)
							, error=function(e) .emptyGene(e)
							)
						})
				} else {
					results <- lapply(ratesToTest, function(interpRates) {
						tryCatch(
							optimParamsSimple(interpRates
								, tpts_exp=tpts
								, alpha_exp=modelSynthesisRate
								, alpha_var=rates$alpha_var[i]
								, total_exp=modelTotalRNA
								, total_var=concentrations$total_var[i]
								, maxit=nIter
								, log_shift=log_shift
								, ode=deSolve::ode
								, time_transf=time_transf
								, .rxnrateSimple=.rxnrateSimple
								, .makeSimpleModel=.makeSimpleModel
								, logLikelihoodFunction=logLikelihoodFunction
								, .emptyGene=.emptyGene
								)
							, error=function(e) .emptyGene(e)
							)
						})
				}
				## sometimes tryCatch is not able to return
				## the empty gene when an error occurr
				check <- sapply(results, length) == 10
				if( !all(check) ) {
					for( i in which(!check) ) {
						results[[i]] <- .emptyGene('Unknown error')
					}
				}
				return(results)
			}
			# , error=function(e) {
			# 	return(list('0'=.emptyGene(error=e), 'a'=.emptyGene(error=e)
			# 		, 'b'=.emptyGene(error=e), 'c'=.emptyGene(error=e)
			# 		, 'ab'=.emptyGene(error=e), 'bc'=.emptyGene(error=e)
			# 		, 'ac'=.emptyGene(error=e), 'abc'=.emptyGene(error=e)
			# 		))
			# 	}
			# ) ## end of tryCatch({
		) ## end of paramAttempts <- sapply(1:nAttempts, function(k)
		if( verbose ) {
			if( is.null(rownames(concentrations$total)[i]) )
				message('Gene "no_name" completed.' )
			else message(paste('Gene "', 
				rownames(concentrations$total)[i],'" completed.', sep=''))
			}
		## choose the best model for each test, out of the many attempts
		chisqPvals <- sapply(1:nAttempts, function(i) 
			sapply(paramAttempts[,i], '[[', 'test'))

		ix <- apply(chisqPvals, 1, function(x) {
			x <- na.omit(x)
			if(length(x)>0) return(which.min(x)) else return(1)
			})
		## in case the model is simple there is no minumum in all tests
		## involving 'c', therefore force to assign to them the first attempt
		if( !intExMode ) ix[grep('c', names(ix))] <- 1
		## extract the selected model, re-assign names and return
		selectedParams <- paramAttempts[cbind(1:length(ix),unlist(ix))]
		names(selectedParams) <- rownames(paramAttempts)
		return(selectedParams)
	} 

	######################
	## MAIN FUNCTION ###
	##################

	nGenes <- nrow(rates$alpha)
	tpts <- tpts
	paramSpecs <- bplapply(1:nGenes, modelOneGene, seed=seed, 
		.chooseModel=.chooseModel,
		time_transf=time_transf,
		.DimpulseModel=.DimpulseModel,
		.DsigmoidModel=.DsigmoidModel,
		constantModelP=constantModelP,
		.emptyGene=.emptyGene,
		optimParams=optimParams,
		optimParamsSimple=optimParamsSimple,
		verbose=verbose,
		nAttempts=nAttempts,
		concentrations=concentrations,
		rates=rates,
		tpts=tpts,
		log_shift=log_shift,
		na.rm=na.rm,
		sigmoidTotal=sigmoidTotal,
		nInit=nInit,
		nIter=nIter,
		sigmoidSynthesis=sigmoidSynthesis,
		testOnSmooth=testOnSmooth,
		estimateRatesWith=estimateRatesWith,
		sigmoidModel=sigmoidModel,
		impulseModel=impulseModel,
		sigmoidModelP=sigmoidModelP,
		impulseModelP=impulseModelP,
		.polynomialModelP=.polynomialModelP,
		.rxnrate=.rxnrate,
		.rxnrateSimple=.rxnrateSimple,
		.makeModel=.makeModel,
		.makeSimpleModel=.makeSimpleModel,
		logLikelihoodFunction=logLikelihoodFunction,
		BPPARAM=BPPARAM
		)
	return(paramSpecs)
}

############### pointer function

newPointer <- function(inputValue){  
	object=new.env(parent=globalenv())  
	object$value=inputValue  
	class(object)='pointer'
	return(object)  
}

############### constant

constantModel <- function(x , par ) rep(par , length(x) )
constantModelP <- newPointer(constantModel)

############### sigmoid

# cppFunction('
# NumericVector sigmoidModelC(NumericVector x, NumericVector par) {
# 	int n = x.size();
# 	NumericVector ans(n);
# 	for(int i = 0; i < n; i++) {
# 		ans[i] = par[0]+(par[1]-par[0])*(1/(1+exp(-par[3]*(x[i]-par[2]))));
# 	}
# 	return ans;
# }
# ')
# sigmoidModel <- sigmoidModelC

sigmoidModel <- function(x, par) 
{
	# h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
	par[1]+(par[2]-par[1])*(1/(1+exp(-par[4]*(x-par[3]))))
}
## compiled version
#sigmoidModel <- cmpfun(sigmoidModel)

.DsigmoidModel <- function(x, par) 
{
     h0= par[1]; h1=par[2]; t1=par[3]; b=par[4]
     S= function(b,t) 1/(1+exp(-b*(x-t)))
     dSdx= function(b,t) b/(1/exp(-b*(x-t)) + 2 + exp(-b*(x-t)) )
     s= function(x,t,h,b) h+(h1-h)*S(b,t)
     dsdx= function(x,t,h,b) (h1-h)*dSdx(b,t)
     1/h1*dsdx(x,t1,h0,b)
}

# 'pointer' for the sigmoidModel function
sigmoidModelP <- newPointer(sigmoidModel)

############### impulse

# cppFunction('
# NumericVector impulseModelC(NumericVector x, NumericVector par) {
# 	int n = x.size();
# 	NumericVector ans(n);
# 	for(int i = 0; i < n; i++) {
# 		ans[i] = 1/par[1]*(par[0]+(par[1]-par[0])*(1/(1+exp(-par[5]*(x[i]-par[3])))))*(par[2]+(par[1]-par[2])*(1/(1+exp(par[5]*(x[i]-par[4])))));
# 	}
# 	return ans;
# }
# ')
# impulseModel <- impulseModelC

impulseModel <- function(x, par) 
{
	# h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
	1/par[2]*(par[1]+(par[2]-par[1])*(1/(1+exp(-par[6]*(x-par[4])))))*
		(par[3]+(par[2]-par[3])*(1/(1+exp(par[6]*(x-par[5])))))
}
## compiled version
#impulseModel <- cmpfun(impulseModel)

.DimpulseModel <- function(x, par) 
{
     h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
     S= function(b,t) 1/(1+exp(-b*(x-t)))
     dSdx= function(b,t) b/(1/exp(-b*(x-t)) + 2 + exp(-b*(x-t)) )
     s= function(x,t,h,b) h+(h1-h)*S(b,t)
     dsdx= function(x,t,h,b) (h1-h)*dSdx(b,t)
     1/h1*(dsdx(x,t1,h0,b)*s(x,t2,h2,-b) + s(x,t1,h0,b)*dsdx(x,t2,h2,-b) )
}

# 'pointer' for the impulseModel function
impulseModelP <- newPointer(impulseModel)

.D2impulseModel <- function(t, par) {
  h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
  -(2*b^2*(h1-h0)*(h1-h2)*exp(b*(t-t2)-b*(t-t1)))/(h1*(exp(-b*(t-t1))+1)^2*(exp(b*(t-t2))+1)^2)+((h1-h2)*((2*b^2*exp(2*b*(t-t2)))/(exp(b*(t-t2))+1)^3-(b^2*exp(b*(t-t2)))/(exp(b*(t-t2))+1)^2)*((h1-h0)/(exp(-b*(t-t1))+1)+h0))/h1+((h1-h0)*((2*b^2*exp(-2*b*(t-t1)))/(exp(-b*(t-t1))+1)^3-(b^2*exp(-b*(t-t1)))/(exp(-b*(t-t1))+1)^2)*((h1-h2)/(exp(b*(t-t2))+1)+h2))/h1
}

############### polynomial

.polynomialModel <- function(x, par)
	sapply(x, function(x_i)
		sum(sapply(1:length(par), function(i) x_i^(i-1) * par[i])))

.polynomialModelP <- newPointer(.polynomialModel)

.chooseModel <- function(tpts, log_shift, experiment, variance=NULL, na.rm=TRUE
	, sigmoid=TRUE, impulse=TRUE, polynomial=TRUE, nInit=10, nIter=500
	, time_transf, impulseModel, sigmoidModel, sigmoidModelP, impulseModelP
	, .polynomialModelP)
#### choose a functional form between impulse and sigmoid according 
#### to the one that has the gratest pvalue in the chi squared test
{
	chisq.test.default <- function(experiment, model, variance=NULL, df)
	{
		if( is.null(variance) ) variance <- stats::var(experiment )
		D = chisqFunction(experiment, model, variance)
		modelDF <- max(0, length(experiment)-df)
		pchisq(D,  modelDF, lower.tail=TRUE)
	}

	optimFailOut <- function(e)
		list(par=NA, value=NA, counts=NA, convergence=1, message=e)
	#
	# impulse model functions
	#
	im.parguess <- function(tpts , values ) {
	    # values = expressions.avgd(eD)
	    # tp = tpts(eD)
	    ntp   <- length(tpts)
	    peaks <- which(diff(sign(diff(values)))!=0)+1
	    if( length(peaks) == 1 ) peak <- peaks
	    if( length(peaks)  > 1 ) peak <- sample(peaks, 1)
	    if( length(peaks) == 0 ) peak <- round(length(tpts)/2)
	    #
		initial_values <- runif( 1, min=min(values[1:3])
			, max=max(values[1:3]))
		intermediate_values <- values[peak]
		if( intermediate_values==0 ) intermediate_values <- mean(values[seq(peak-1,peak+1)])
		end_values <- runif( 1, min=min(values[(ntp-2):ntp])
			, max=max(values[(ntp-2):ntp]))
		time_of_first_response  <- tpts[peak-1]
		time_of_second_response <- tpts[peak+1]
		slope_of_response <- diff(range(tpts)) / 
			(time_of_second_response-time_of_first_response)
		#
	    return(c(h0=initial_values, h1=intermediate_values
	    	, h2=end_values, t1=time_of_first_response
	    	, t2=time_of_second_response, b=slope_of_response))
	}
	#
	im.chisq <- function(par, tpts, experiment, variance=NULL, impulseModel) 
	{
		 model <- impulseModel(tpts, par)
		 chisqFunction(experiment, model, variance)
	}
	#
	im.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
		, maxit=500) 
		sapply(1:ninit, function(x) 
	 		tryCatch(optim(
	 			par=im.parguess(tpts, experiment)
	 			, fn=im.chisq, tpts=tpts
	 			, experiment=experiment
	 			, variance=variance
	 			, impulseModel=impulseModel
	 			, control=list(maxit=maxit)
	 			), error=function(e) optimFailOut(e)))
	#
	# sigmoid model functions
	#
	sm.parguess <- function(tpts , values ) {
	    # values = expressions.avgd(eD)
	    # tp = tpts(eD)

		time_span <- diff(range(tpts))
		# sample the time uniformely
		time_of_response <- runif( 1, min=min(tpts), max=max(tpts))
		# slope of response must be high if the time of response is close to one
		# of the two boundaries
		distance_from_boundary <- min(time_of_response - min(tpts)
				, max(tpts) - time_of_response)
		slope_of_response <- time_span / distance_from_boundary
	    ntp   <- length(tpts)
		initial_values <- runif( 1, min=min(values[1:3])
			, max=max(values[1:3]))
		end_values <- runif( 1, min=min(values[(ntp-2):ntp])
			, max=max(values[(ntp-2):ntp]))
		#
	    return(c(h0=initial_values, h1=end_values, t1=time_of_response
	    	, b=slope_of_response))
	}
	#
	sm.chisq <- function(par, tpts, experiment, variance=NULL, sigmoidModel) 
	{
		 model <- sigmoidModel(tpts, par)
		 chisqFunction(experiment, model, variance)
	}
	#
	sm.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
		, maxit=500) 
		sapply(1:ninit, function(x) 
				tryCatch(optim(
					par=sm.parguess(tpts, experiment)
					, fn=sm.chisq, tpts=tpts
					, experiment=experiment
					, variance=variance
					, sigmoidModel=sigmoidModel
					, control=list(maxit=maxit)
					), error=function(e) optimFailOut(e))) 

	pn.optim.aic <- function(tpts , experiment, variance=NULL) {
		if( length(experiment) < 3 ) return(NA)
		polyOrderChisq <- function(i) {
			model <- lm(experiment~poly(tpts, i, raw=TRUE ))
			return(list(par=model$coeff, value=AIC(model)))}
		return(sapply(1:min(7,length(tpts)-2), polyOrderChisq))
	}

	# remove missing values
	if( na.rm) {
		idx <- is.finite(experiment)
		tpts <- tpts[idx]
		experiment <- experiment[idx]
	}

	## 
	if( length(experiment)==0 ) {
		stop('.chooseModel: no time points have a finite value.
			Impossible to evaluate any kind of model.')
		return(list(type='constant', fun=constantModelP
			, params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
	}
	## 
	if( length(experiment)<=2 ) {
		warning('.chooseModel: less than three time points have a finite value. 
			Impossible evaluate a variable model.
			Returning a constant model.')
		return(list(type='constant', fun=constantModelP
			, params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
	}

	## re-evaluate flags of function to evaluate according to the lenght 
	## of the experiment
	sigmoid <- sigmoid
	impulse <- impulse & length(experiment)>2
	polynomial <- polynomial & length(experiment)>2

	tptslog <- time_transf(tpts, log_shift)

	# sigmoid
	if( sigmoid ) {
		outSM  <- sm.optim.chisq(tpts=tptslog, experiment=experiment
			, variance=variance, ninit=nInit, maxit=nIter)
		bestSM <- which.min(unlist(outSM[2,]))
		pvalSM <- chisq.test.default(experiment=experiment
			, model=sigmoidModel(tptslog, outSM[,bestSM]$par)
			, variance=variance, df=length(outSM[,bestSM]$par))
		dfSM <- length(outSM[,bestSM]$par)
	} else dfSM <- NA
	# impulse
	if( impulse ) {
		outIM  <- im.optim.chisq(tpts=tptslog, experiment=experiment, 
			variance=variance, ninit=nInit, maxit=nIter)
		bestIM <- which.min(unlist(outIM[2,]))
		pvalIM <- chisq.test.default(experiment=experiment
			, model=impulseModel(tptslog, outIM[,bestIM]$par) 
			, variance=variance, df=length(outIM[,bestIM]$par))
		dfIM <- length(outIM[,bestIM]$par)
	} else dfIM <- NA

	# polynomial
	if( polynomial ) {
		outPN  <- pn.optim.aic(tptslog, experiment, variance )
		bestPN <- which.min(unlist(outPN[2,]))
		pvalPN <- chisq.test.default(experiment=experiment
			, model=.polynomialModel(tptslog, outPN[,bestPN]$par) 
			, variance=variance, df=length(outPN[,bestPN]$par))
		dfPN <- length(outPN[,bestPN]$par)
	} else dfPN <- NA

	pvals <- c(
		sigmoid=if( sigmoid ) pvalSM else NA
		, impulse=if( impulse ) pvalIM else NA
		, polynomial=if( polynomial ) pvalPN else NA
		)
	funcs <- c(sigmoidModelP, impulseModelP, .polynomialModelP)
	dfs <- c(dfSM, dfIM, dfPN)
	type <- names(pvals)[which.min(pvals)]
	df   <- dfs[which.min(pvals)]

	if( type=='sigmoid'    ) params <- outSM[,bestSM]$par
	if( type=='impulse'    ) params <- outIM[,bestIM]$par
	if( type=='polynomial' ) params <- outPN[,bestPN]$par

	pval <- pvals[which.min(pvals)]
	fun  <- funcs[[which.min(pvals)]]

	return(list(type=type, fun=fun , params=params, pval=pval, df=df))

}


####################################################################################################k1KKK_NoNascent <- function(x, par){par[1]*par[3]}

fitSmooth <- function(tpts
		            , tt_c
		            , experiment
		            , variance
		            , nInit=20
		            , nIter=500
		            , mature = FALSE
		            , seed = NULL)
{
	im_parguess <- function(tpts , values) {
    	
    	ntp   <- length(tpts)
    	peaks <- which(diff(sign(diff(values)))!=0)+1
    	if( length(peaks) == 1 ) peak <- peaks
    	if( length(peaks)  > 1 ) peak <- sample(peaks, 1)
    	if( length(peaks) == 0 ) peak <- round(length(tpts)/2)
    	
    	initial_values <- runif( 1, min=min(values[1:3]), max=max(values[1:3]))
    	
    	intermediate_values <- values[peak]
    	if( intermediate_values==0 ) intermediate_values <- mean(values[seq(peak-1,peak+1)])
    	end_values <- runif( 1, min=min(values[(ntp-2):ntp]), max=max(values[(ntp-2):ntp]))

    	time_of_first_response  <- tpts[peak-1]
    	time_of_second_response <- tpts[peak+1]
    
        slope_of_response <- 1

        par <- c(h0=initial_values
        		,h1=intermediate_values
        		,h2=end_values
        		,t1=time_of_first_response
        		,t2=time_of_second_response
        		,b=slope_of_response)

	    return(unlist(unname(par)))
	}

	im_chisq_mature <- function(par, tpts, experiment, variance=NULL, tt_c)
	{
		model <- impulseModel(tpts,par)
		if( abs(par[6]) > Inf ) return(NaN)
		if( any(model < 0) ) return(NaN)
		chisqFunction(experiment, model, variance)
	}

	im_chisq <- function(par, tpts, experiment, variance=NULL, tt_c) 
	{
		model <- impulseModel(tpts,par)
		if( any(model < 0) ) return(NaN)
		chisqFunction(experiment, model, variance)
	}
  
	if(is.numeric(seed)) set.seed(seed)
  	outIM  <- sapply(1:nInit, function(x) 
    	tryCatch(optim(
      		par=im_parguess(tpts, experiment)
		  , fn=if(mature) im_chisq_mature else im_chisq
		  , tpts=tpts
		  , experiment=experiment
		  , variance=variance
		  , tt_c = tt_c
		  , control=list(maxit=nIter)
      	), error=function(e) list(par=NA
      							, value=NA
      							, counts=NA
      							, convergence=1, message=e)))

	bestIM <- which.min(unlist(outIM[2,]))
	unlist(outIM[,bestIM])
}

prematureKKK_Int_NoNascent <- function(x, parameters)
{
  matureParameters <- parameters[1]
  k2Parameters <- parameters[2]
  k3Parameters <- parameters[3]
  
  return((k3Parameters*matureParameters)/k2Parameters)
}

k1KKK_Int_NoNascent <- function(x, par)
{
  par[1]*par[3]
}

systemSolution <- function(k1F,k2F,k3F,times)
{

  system <- function(t,c,parms)
  {
    alpha <- parms$alpha
    beta  <- parms$beta
    gamma <- parms$gamma

    r=rep(0,length(c))
    r[1] <- alpha(t) - beta(t) * c["p"]
    r[2] <- beta(t) * c["p"] - gamma(t) * c["m"]
    return(list(r))
  }

  cinit <- c(k1F(0)/k2F(0),k1F(0)/k3F(0))
  names(cinit) <- c("p","m")

  params <- list(alpha = k1F, beta = k2F, gamma = k3F)

  modData <- ode(y=cinit, times=times, func=system, parms=params)
  modData <- c(modData[,"m"],modData[,"p"])

  names(modData) <- c(rep("mature",length.out = length(times)),rep("premature",length.out = length(times)))

  return(modData) 
}

errorKKK_Int_NoNascent <- function(parameters, tpts, premature, mature, prematureVariance, matureVariance)
{

  if(parameters[1]<0)return(NaN)
  if(parameters[2]<0)return(NaN)
  if(parameters[3]<0)return(NaN)

  matureParameters <- parameters[1]

  prematureEstimated <- prematureKKK_Int_NoNascent(x = tpts, parameters = parameters)
  matureEstimated <- matureParameters

  prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
  matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

  return(sum(c(prematureChiSquare,matureChiSquare)))
}

errorVKK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) return(parameters[8])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorVKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorVVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[13:18])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorKVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[2:7])}
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorVVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

	k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
	k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
	k3F <- function(x) return(parameters[13])

	if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

	modData <- systemSolution(k1F,k2F,k3F,times)
	
	chi2 <- chisqFunction(data,modData,datavar)

	return(chi2)
}

errorKVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{
	k1F <- function(x)return(parameters[1])
	k2F <- function(x){impulseModel(log2(x+a)+c,parameters[2:7])}
	k3F <- function(x)return(parameters[8])
	
	if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)
	
	modData <- systemSolution(k1F,k2F,k3F,times)
	
	chi2 <- chisqFunction(data,modData,datavar)
	
	return(chi2)
}

errorKKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) return(parameters[2])
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[3:8])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

loglikKKK_Int_NoNascent <- function(parameters
	                    	   ,tpts
	                    	   ,premature
	                    	   ,mature
	                    	   ,prematureVariance
	                    	   ,matureVariance)
{

	matureParameters <- parameters[1]

	prematureEstimated <- prematureKKK_Int_NoNascent(x = tpts, parameters = parameters)
	matureEstimated <- matureParameters

	logLikelihoodFunction(premature, prematureEstimated, prematureVariance) + 
	logLikelihoodFunction(mature, matureEstimated, matureVariance)

}

loglikVKK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) return(parameters[8])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)

}

loglikKVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[2:7])}
  k3F <- function(x) return(parameters[8])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  logLikelihoodFunction(data, modData, datavar)
}

loglikKKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) return(parameters[2])
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[3:8])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
  k3F <- function(x) return(parameters[13])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikKVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[2:7])}
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
  k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[13:18])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}



###########################################################################

k1VKK_Der_NoNascent <- function(x, par, c)
{
	t_fact <- 2^(x-c)*log(2)
	.D2impulseModel(x, par[1:6])*t_fact^2/par[7] + .DimpulseModel(x, par[1:6])*(1+par[8]/par[7])*t_fact + par[8]*impulseModel(x, par[1:6])
}

prematureVKK_Der_NoNascent <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7]
	k3Parameters <- parameters[8]

	t_fact <- 2^(x-c)*log(2)
	(.DimpulseModel(x, matureParameters)*t_fact + k3Parameters * impulseModel(x, matureParameters))/k2Parameters

}

errorVKK_Der_NoNascent <- function(parameters
							  ,tpts
							  ,premature
							  ,mature
							  ,prematureVariance
							  ,matureVariance
							  ,c)
{

	matureParameters <- parameters[1:6]

  	if( abs(matureParameters[6]) > Inf ) return(NaN)
	D2 <- .D2impulseModel(tpts,matureParameters)
	k1 <- k1VKK_Der_NoNascent(tpts,parameters, c)

	prematureEstimated <- prematureVKK_Der_NoNascent(x = tpts, parameters = parameters, c = c)
	matureEstimated <- impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

k1VVK_Der_NoNascent <- function(x, par, c)
{
	t_fact <- 2^(x-c)*log(2)
	.D2impulseModel(x, par[1:6])/impulseModel(x, par[7:12])*t_fact^2 +
	.DimpulseModel(x, par[1:6])*t_fact*(1 - log(2)*.DimpulseModel(x, par[7:12])/impulseModel(x, par[7:12])^2 + par[13]/impulseModel(x, par[7:12])) + 
	log(2)*impulseModel(x, par[1:6])*(par[13]/log(2) - (par[13]*.DimpulseModel(x, par[7:12]))/impulseModel(x, par[7:12])^2 )
}

prematureVVK_Der_NoNascent <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7:12]
	k3Parameters <- parameters[13]
	
	t_fact <- 2^(x-c)*log(2)
	return((.DimpulseModel(x, matureParameters)*t_fact
	      + k3Parameters * impulseModel(x, matureParameters))/impulseModel(x, k2Parameters))
}

errorVVK_Der_NoNascent <- function(parameters
                    	 ,tpts
                    	 ,premature
                    	 ,mature
                    	 ,prematureVariance
                    	 ,matureVariance
                    	 ,c)
{

	matureParameters <- parameters[1:6]

	if( abs(matureParameters[6]) > Inf ) return(NaN)
	if( abs(parameters[12]) > Inf ) return(NaN)

	D2 <- .D2impulseModel(tpts,matureParameters)
	k1 <- k1VVK_Der_NoNascent(tpts,parameters, c)

	prematureEstimated <- prematureVVK_Der_NoNascent(x = tpts, parameters = parameters, c = c)
	matureEstimated <- impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

k1VKV_Der_NoNascent <- function(x, par, c)
{
	t_fact <- 2^(x-c)*log(2)
	.D2impulseModel(x, par[1:6])/par[7]*t_fact^2 +
	.DimpulseModel(x, par[1:6])*t_fact*(1 + impulseModel(x, par[8:13])/par[7]) + 
	log(2)*impulseModel(x, par[1:6])*( .DimpulseModel(x, par[8:13])/par[7] + impulseModel(x, par[8:13])/log(2) )
}

prematureVKV_Der_NoNascent <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7]
	k3Parameters <- parameters[8:13]
	
	t_fact <- 2^(x-c)*log(2)
	(.DimpulseModel(x, matureParameters)*t_fact + impulseModel(x, k3Parameters) * impulseModel(x, matureParameters))/k2Parameters
}

errorVKV_Der_NoNascent <- function(parameters
									 ,tpts
									 ,premature
									 ,mature
									 ,prematureVariance
									 ,matureVariance
									 ,c)
{
	matureParameters <- parameters[1:6]

	if( abs(matureParameters[6]) > Inf ) return(NaN)
	if( abs(parameters[13]) > Inf ) return(NaN)

	D2 <- .D2impulseModel(tpts,matureParameters)
	k1 <- k1VKV_Der_NoNascent(tpts,parameters, c)

	prematureEstimated <- prematureVKV_Der_NoNascent(x = tpts, parameters = parameters, c = c)
	matureEstimated <- impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

k1VVV_Der_NoNascent <- function(x, par, c)
{
  t_fact <- 2^(x-c)*log(2)
  .D2impulseModel(x, par[1:6])/impulseModel(x, par[7:12])*t_fact^2 +
  .DimpulseModel(x, par[1:6])*t_fact*(1 - log(2)*.DimpulseModel(x, par[7:12])/impulseModel(x, par[7:12])^2 + impulseModel(x, par[13:18])/impulseModel(x, par[7:12])) + 
  log(2)*impulseModel(x, par[1:6])*(.DimpulseModel(x, par[13:18])/impulseModel(x, par[7:12]) + impulseModel(x, par[13:18])/log(2) - (impulseModel(x, par[13:18])*.DimpulseModel(x, par[7:12]))/impulseModel(x, par[7:12])^2 )
}

prematureVVV_Der_NoNascent <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7:12]
	k3Parameters <- parameters[13:18]

	t_fact <- 2^(x-c)*log(2)
	return((.DimpulseModel(x, matureParameters)*t_fact
		+ impulseModel(x, k3Parameters)*impulseModel(x, matureParameters))/impulseModel(x, k2Parameters))
}

errorVVV_Der_NoNascent <- function(parameters
									 ,tpts
									 ,premature
									 ,mature
									 ,prematureVariance
									 ,matureVariance
									 ,c)
{

	matureParameters <- parameters[1:6]

	if( abs(matureParameters[6]) > Inf ) return(NaN)
	if( abs(parameters[12]) > Inf ) return(NaN)
	if( abs(parameters[18]) > Inf ) return(NaN)

	D2 <- .D2impulseModel(tpts,matureParameters)
	k1 <- k1VVV_Der_NoNascent(tpts,parameters, c)

	prematureEstimated <- prematureVVV_Der_NoNascent(x = tpts, parameters = parameters, c = c)
	matureEstimated <- impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

.inspect.engine_Derivative_NoNascent <- function(tptsOriginal
										   , tptsLinear
										   , a
										   , c
										   , concentrations
										   , rates
										   , BPPARAM=bpparam()
										   , na.rm=TRUE
										   , verbose=TRUE
										   , testOnSmooth=TRUE
										   , seed=NULL
										   , nInit = nInit
										   , nIter = nIter)
{

	total <- concentrations$total
	totalVariance <- concentrations$total_var

	premature <- concentrations$preMRNA
	prematureVariance <- concentrations$preMRNA_var

	mature <- concentrations$mature
	matureVariance <- concentrations$mature_var

	eiGenes <- rownames(mature)

	k2median <- median(rates$gamma,na.rm = TRUE)
	k3median <- median(rates$beta,na.rm = TRUE)

	# degreesOfFreedom <- length(tptsLinear)-6

	matureFitImpulse <- bplapply(eiGenes,function(row)
	{
  		fitSmooth(tpts = tptsLinear
        	    , tt_c = c
        	    , experiment = mature[row,]
        	    , variance = matureVariance[row,]
        	    , mature = TRUE
        	    , nInit = nInit
        	    , nIter = nIter
        	    , seed = seed)
	},BPPARAM=BPPARAM)
	names(matureFitImpulse) <- eiGenes

	if(all(sapply(matureFitImpulse,is.null)))
		stop("No genes have a mature profile possible to fit 
			with impulsive functions. Try with the options:
			'modelingParams()$estimateRatesWith <- 'int' ' and
			'modelingParams()$testOnSmooth <- FALSE'.")

	if(any(sapply(matureFitImpulse,is.null))) {
		message("Some genes have a mature profile impossible to be fitted with impulsive 
			functions therefore they will be excluded from the modelling.")

		eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null)]

		matureFitImpulse <- matureFitImpulse[eiGenes]
		total <- total[eiGenes,]
		totalVariance <- totalVariance[eiGenes,]
		premature <- premature[eiGenes,]
		prematureVariance <- prematureVariance[eiGenes,]
		mature <- mature[eiGenes,]
		matureVariance <- matureVariance[eiGenes,]
	}

	# if( degreesOfFreedom > 0 ) {
	# 	accept <- pchisq(sapply(matureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2 
	# 	if( table(accept)['FALSE']/length(accept) > .5 )
	# 		message("More than 50% did not return a good fit of their mature profiles 
	# 			with the impulsive smooth function:
	# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
	# }

	if( testOnSmooth ) {

		prematureFitImpulse <- bplapply(eiGenes,function(row)
		{
			fitSmooth(tpts = tptsLinear
					, tt_c = c
	        		, experiment = premature[row,]
	        		, variance = prematureVariance[row,]
	        		, mature = FALSE
	        		, nInit = nInit
	        		, nIter = nIter
	        		, seed = seed)     
		},BPPARAM=BPPARAM)
		names(prematureFitImpulse) <- eiGenes

		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
			stop("No genes have an expression profile possible to fit 
				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")

		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
			message("Some genes have an expression profile impossible to be fitted 
				with impulsive functions therefore they will be excluded from the modelling.")

			eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null) & 
				!sapply(prematureFitImpulse,is.null)]

			prematureFitImpulse <- prematureFitImpulse[eiGenes]
			matureFitImpulse <- matureFitImpulse[eiGenes]
			total <- total[eiGenes,]
			totalVariance <- totalVariance[eiGenes,]
			premature <- premature[eiGenes,]
			prematureVariance <- prematureVariance[eiGenes,]
			mature <- mature[eiGenes,]
			matureVariance <- matureVariance[eiGenes,]
		}

		# if( degreesOfFreedom > 0 ) {
		# 	accept <- pchisq(sapply(prematureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2 
		# 	if( table(accept)['FALSE']/length(accept) > .5 )
		# 		message("More than 50% did not return a good fit of their premature 
		# 			profiles with the impulsive smooth function:
		# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
		# }

	}

	# Equal to integrative approach

	KKK <- bplapply(eiGenes,function(row){
	
	matureParameters <- mean(impulseModel(tptsLinear, matureFitImpulse[[row]][1:6]))
	k2Parameters <- k2median
	k3Parameters <- k3median

	unlist(
		tryCatch(
			optim(c(matureParameters, k2Parameters, k3Parameters)
				,errorKKK_Int_NoNascent
				,tpts = tptsLinear
				,premature = premature[row,]
				,mature = mature[row,]
				,prematureVariance = prematureVariance[row,]
				,matureVariance = matureVariance[row,]
				,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							  , par2 = NaN
							  , par3 = NaN
							  , value = NaN
							  , counts.function = NaN
						  	  , counts.gradient = NaN
						  	  , convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KKK) <- eiGenes
	message("Model 0 finished.")

	medianAmplVKK <- sapply(eiGenes, function(row)
	{

		parameters <- unname(c(matureFitImpulse[[row]][1:6]
							 , k2median
							 , k3median))

		k1 <- k1VKK_Der_NoNascent(tptsLinear,parameters,c)
		p <- prematureVKK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)

		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)

		suppressWarnings(optimize( function(x)
		{
    		parameters <- unname(c(matureFitImpulse[[row]][1:6]
								 , k2median*x
								 , k3median*x))

			k1 <- k1VKK_Der_NoNascent(tptsLinear,parameters,c)
    		p <- prematureVKK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)

		if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(c(k1
												 ,k2median*x*length(tptsOriginal)
      									 		 ,k3median*x*length(tptsOriginal)))
  		},c(1, 1e5) ))$minimum
	})

	VKK <- bplapply(eiGenes,function(row){

		if(medianAmplVKK[[row]] > 10^2)
  		{
			return(c(par1 = NaN, par2 = NaN, par3 = NaN
				   , par4 = NaN, par5 = NaN, par6 = NaN
				   , par7 = NaN, par8 = NaN, value = NaN
				   , counts.function = NaN
				   , counts.gradient = NaN, convergence = NaN))
		}

		matureParameters <- unname(matureFitImpulse[[row]][1:6])

		k2Parameters <- k2median * unname(medianAmplVKK[row])
		k3Parameters <- k3median * unname(medianAmplVKK[row])
	  
		unlist(
			tryCatch(
	      			optim(unname(c(matureParameters, k2Parameters, k3Parameters))
	                  	 ,errorVKK_Der_NoNascent
        				 ,tpts = tptsLinear
			             ,premature = premature[row,]
			             ,mature = mature[row,]
			             ,prematureVariance = prematureVariance[row,]
			             ,matureVariance = matureVariance[row,]
	                  	 ,c = c
	                  	 ,control = list(maxit = nIter * 100)),
	    		error=function(e) c(par1 = NaN
	    						  , par2 = NaN
	    						  , par3 = NaN
	    						  , par4 = NaN
	    						  , par5 = NaN
	    						  , par6 = NaN
	    						  , par7 = NaN
	    						  , par8 = NaN
	    						  , value = NaN
	    						  , counts.function = NaN
	    						  , counts.gradient = NaN
	    						  , convergence = NaN)
	    		)
	    )
	}, BPPARAM=BPPARAM)
	names(VKK) <- eiGenes
	message("Model A finished.")

	medianAmplVVK <- sapply(eiGenes, function(row)
	{
		parameters <- c(unname(matureFitImpulse[[row]][1:6])
						, rep(k2median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1
						, k3median)
		k1 <- k1VVK_Der_NoNascent(tptsLinear,parameters, c)
		p <- prematureVVK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)
		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)
		suppressWarnings(optimize( function(x) {
			parameters <- unname(c(matureFitImpulse[[row]][1:6]
							  , c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
							  , k3median * x))
     
			k1 <- k1VVK_Der_NoNascent(tptsLinear,parameters, c)
			p <- prematureVVK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)      
     
			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(c(k1,impulseModel(tptsLinear,c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)),k3median*x*length(tptsOriginal)))
   		}, c(1, 1e5) ))$minimum
	})

	VVK <- bplapply(eiGenes, function(row)
	{

		if(medianAmplVVK[[row]] > 10^2)
		{
			return(c(par1 = NaN, par2 = NaN, par3 = NaN
				   , par4 = NaN, par5 = NaN, par6 = NaN
				   , par7 = NaN, par8 = NaN, par9 = NaN
				   , par10 = NaN, par11 = NaN, par12 = NaN
				   , par13 = NaN, value = NaN, counts.function = NaN
				   , counts.gradient = NaN, convergence = NaN))
		}
  
		matureParameters <- unname(matureFitImpulse[[row]][1:6])
		k2Parameters <- c(rep(k2median,3) * unname(medianAmplVVK[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k3Parameters <- k3median * unname(medianAmplVVK[row])
		
		unlist(
			tryCatch(
				optim(unname(c(matureParameters, k2Parameters, k3Parameters))
			        	,errorVVK_Der_NoNascent
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,c = c
						,control = list(maxit = nIter * 100)),
				error=function(e) c(par1 = NaN
								   ,par2 = NaN
								   ,par3 = NaN
								   ,par4 = NaN
								   ,par5 = NaN
								   ,par6 = NaN
								   ,par7 = NaN
								   ,par8 = NaN
								   ,par9 = NaN
								   ,par10 = NaN
								   ,par11 = NaN
								   ,par12 = NaN
								   ,par13 = NaN
								   ,value = NaN
								   ,counts.function = NaN
								   ,counts.gradient = NaN
								   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(VVK) <- eiGenes
	message("Model AC finished.")

	medianAmplVKV <- sapply(eiGenes, function(row)
	{

		parameters <- c(unname(matureFitImpulse[[row]][1:6])
							 , k2median
			  				 , rep(k3median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
  		k1 <- k1VKV_Der_NoNascent(tptsLinear,parameters, c)
		p <- prematureVKV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)

		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)

		suppressWarnings(optimize( function(x) {

			parameters <- c(unname(matureFitImpulse[[row]][1:6])
								 , k2median * x
								 , rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
      
			k1 <- k1VKV_Der_NoNascent(tptsLinear,parameters, c)
			p <- prematureVKV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)
      
			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(k1,k2median*x*length(tptsOriginal),impulseModel(tptsLinear,c(rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)))

		}, c(1, 1e5) ))$minimum
	})

	VKV <- bplapply(eiGenes, function(row)
	{

		if(medianAmplVKV[[row]] > 10^2)
		{
			return(c(par1 = NaN, par2 = NaN, par3 = NaN
				   , par4 = NaN, par5 = NaN, par6 = NaN
				   , par7 = NaN, par8 = NaN, par9 = NaN
				   , par10 = NaN, par11 = NaN, par12 = NaN
				   , par13 = NaN, value = NaN, counts.function = NaN
				   , counts.gradient = NaN, convergence = NaN))
		}

		matureParameters <- unname(matureFitImpulse[[row]][1:6])
		k2Parameters <- k2median * unname(medianAmplVKV[row])
		k3Parameters <- c(rep(k3median,3) * unname(medianAmplVKV[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
		unlist(
			tryCatch(
				optim(unname(c(matureParameters, k2Parameters, k3Parameters))
          				,errorVKV_Der_NoNascent
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,c = c
						,control = list(maxit = nIter * 100)),
				error=function(e) c(par1 = NaN
								   ,par2 = NaN
								   ,par3 = NaN
								   ,par4 = NaN
								   ,par5 = NaN
								   ,par6 = NaN
								   ,par7 = NaN
								   ,par8 = NaN
								   ,par9 = NaN
								   ,par10 = NaN
								   ,par11 = NaN
								   ,par12 = NaN
								   ,par13 = NaN
								   ,value = NaN
								   ,counts.function = NaN
								   ,counts.gradient = NaN
								   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(VKV) <- eiGenes
	message("Model AB finished.")

	medianAmplVVV <- sapply(eiGenes, function(row)
	{

		parameters <- c(unname(matureFitImpulse[[row]][1:6])
					  , rep(k2median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1
					  , rep(k3median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k1 <- k1VVV_Der_NoNascent(tptsLinear,parameters, c)
		p <- prematureVVV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)
	
		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)
	
		suppressWarnings(optimize( function(x) {
	
			parameters <- c(unname(matureFitImpulse[[row]][1:6])
						  , rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1
						  , rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
	
			k1 <- k1VVV_Der_NoNascent(tptsLinear,parameters, c)
			p <- prematureVVV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)      
	      
			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(k1
	      										, impulseModel(tptsLinear,c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1))
	      										, impulseModel(tptsLinear,c(rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)))
	
			}, c(1, 1e5) ))$minimum
	})

	VVV <- bplapply(eiGenes, function(row){

		if(medianAmplVVV[[row]] > 10^2)
		{
			return(c(par1 = NaN, par2 = NaN, par3 = NaN, par4 = NaN, par5 = NaN
				   , par6 = NaN, par7 = NaN, par8 = NaN, par9 = NaN, par10 = NaN
				   , par11 = NaN, par12 = NaN, par13 = NaN, par14 = NaN
				   , par15 = NaN, par16 = NaN, par17 = NaN, par18 = NaN
				   , value = NaN, counts.function = NaN
				   , counts.gradient = NaN, convergence = NaN))
		}

		matureParameters <- unname(matureFitImpulse[[row]][1:6])
		k2Parameters <- c(rep(k2median,3) * unname(medianAmplVVV[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k3Parameters <- c(rep(k3median,3) * unname(medianAmplVVV[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)

			unlist(
			tryCatch(
				optim(unname(c(matureParameters, k2Parameters, k3Parameters))
          				,errorVVV_Der_NoNascent
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,c = c
						,control = list(maxit = nIter * 100)),
				error=function(e) c(par1 = NaN
								   ,par2 = NaN
								   ,par3 = NaN
								   ,par4 = NaN
								   ,par5 = NaN
								   ,par6 = NaN
								   ,par7 = NaN
								   ,par8 = NaN
								   ,par9 = NaN
								   ,par10 = NaN
								   ,par11 = NaN
								   ,par12 = NaN
								   ,par13 = NaN
								   ,par14 = NaN
								   ,par15 = NaN
								   ,par16 = NaN
								   ,par17 = NaN
								   ,par18 = NaN
								   ,value = NaN
								   ,counts.function = NaN
								   ,counts.gradient = NaN
								   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(VVV) <- eiGenes
	message("Model ABC finished.")

	KKV <- bplapply(eiGenes, function(row){
 
		k1Parameters <- mean(k1VKV_Der_NoNascent(tptsLinear, VKV[[row]], c),na.rm=T)
		k2Parameters <- VKV[[row]][7]
		k3Parameters <- VKV[[row]][8:13]
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
		        	,errorKKV_Int_NoNascent
		        	,times = tptsOriginal
		        	,data = c( mature[row,], premature[row,] )
		        	,datavar = c( matureVariance[row,] , prematureVariance[row,] )
		        	,a = a
		        	,c = c
		        	,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KKV) <- eiGenes
	message("Model B finished.")

	KVK <- bplapply(eiGenes, function(row){
 
		k1Parameters <- mean(k1VVK_Der_NoNascent(tptsLinear, VVK[[row]], c),na.rm=T)
		k2Parameters <- VVK[[row]][7:12]
		k3Parameters <- VVK[[row]][13]
		
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVK_Int_NoNascent
					 ,times = tptsOriginal
					 ,data = c( mature[row,], premature[row,] )
					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
					 ,a = a
					 ,c = c
					 ,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KVK) <- eiGenes
	message("Model C finished.")

	KVV <- bplapply(eiGenes, function(row){
	
		k1Parameters <- mean(k1VVV_Der_NoNascent(tptsLinear, VVV[[row]], c),na.rm=T)
		k2Parameters <- VVV[[row]][7:12]
		k3Parameters <- VVV[[row]][13:18]
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVV_Int_NoNascent
					 ,times = tptsOriginal
					 ,data = c( mature[row,], premature[row,] )
					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
					 ,a = a
					 ,c = c
					 ,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,par9 = NaN
							   ,par10 = NaN
							   ,par11 = NaN
							   ,par12 = NaN
							   ,par13 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KVV) <- eiGenes
	message("Model BC finished.")

	if(testOnSmooth)
	{
		matureSmooth <- t(sapply(matureFitImpulse,function(i)
		{
 			impulseModel(tptsLinear,i)
		}))
		prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
		{
		 	impulseModel(tptsLinear,i)
		}))
		colnames(matureSmooth) <- tptsOriginal
		colnames(prematureSmooth) <- tptsOriginal
	}else{
		matureSmooth <- mature
		prematureSmooth <- premature
	}

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
											  ,tptsLinear
											  ,prematureSmooth[g,]
											  ,matureSmooth[g,]
											  ,prematureVariance[g,]
											  ,matureVariance[g,]),error = function(e)NaN)
	
		VKKTemp <- tryCatch(errorVKK_Der_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
											  ,tptsLinear
											  ,prematureSmooth[g,]
											  ,matureSmooth[g,]
											  ,prematureVariance[g,]
											  ,matureVariance[g,]
											  ,c),error = function(e)NaN)
	
		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)
	
		VVKTemp <- tryCatch(errorVVK_Der_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
													 ,tptsLinear
													 ,prematureSmooth[g,]
													 ,matureSmooth[g,]
													 ,prematureVariance[g,]
													 ,matureVariance[g,]
													 ,c),error = function(e)NaN)
	
		VKVTemp <- tryCatch(errorVKV_Der_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
													 ,tptsLinear
													 ,prematureSmooth[g,]
													 ,matureSmooth[g,]
													 ,prematureVariance[g,]
													 ,matureVariance[g,]
													 ,c),error = function(e)NaN)

		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)
	
		VVVTemp <- tryCatch(errorVVV_Der_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
													 ,tptsLinear
													 ,prematureSmooth[g,]
													 ,matureSmooth[g,]
													 ,prematureVariance[g,]
													 ,matureVariance[g,]
													 ,c),error = function(e)NaN)
	  
	  c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
	}, BPPARAM=BPPARAM))

	dof <- c(KKK = 3
			,VKK = 8
			,KVK = 8
			,KKV = 8
			,VVK = 13
			,VKV = 13
			,KVV = 13
			,VVV = 18)

	# P values
	pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], min(c(0,2*length(tptsOriginal)-dof['KKK'])))
						,VKK=pchisq(chi2data[,'VKK'], min(c(0,2*length(tptsOriginal)-dof['VKK'])))
						,KVK=pchisq(chi2data[,'KVK'], min(c(0,2*length(tptsOriginal)-dof['KVK'])))
						,KKV=pchisq(chi2data[,'KKV'], min(c(0,2*length(tptsOriginal)-dof['KKV'])))
						,VVK=pchisq(chi2data[,'VVK'], min(c(0,2*length(tptsOriginal)-dof['VVK'])))
						,VKV=pchisq(chi2data[,'VKV'], min(c(0,2*length(tptsOriginal)-dof['VKV'])))
						,KVV=pchisq(chi2data[,'KVV'], min(c(0,2*length(tptsOriginal)-dof['KVV'])))
						,VVV=pchisq(chi2data[,'VVV'], min(c(0,2*length(tptsOriginal)-dof['VVV']))))

	

	logLikelihood <- t(mcsapply(eiGenes,function(g)
	{
		prematureKKKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureKKK_Int_NoNascent(x = tptsLinear[t]
																						   , parameters = KKK[[g]][grep("par",names(KKK[[g]]))])))
		prematureVKKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVKK_Der_NoNascent(x = tptsLinear[t]
		                                                                         , parameters = VKK[[g]][grep("par",names(VKK[[g]]))]
		                                                                         , c = c)))
		prematureVVKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVVK_Der_NoNascent(x = tptsLinear[t]
		                                                            , parameters = VVK[[g]][grep("par",names(VVK[[g]]))]
		                                                            , c = c)))
		prematureVKVTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVKV_Der_NoNascent(x = tptsLinear[t]
		                                                            , parameters = VKV[[g]][grep("par",names(VKV[[g]]))]
		                                                            , c = c)))
		prematureVVVTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVVV_Der_NoNascent(x = tptsLinear[t]
		                                                            , parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
		                                                            , c = c)))		

		matureKKKTemp <- c(sapply(seq_along(tptsLinear),function(t)KKK[[g]][grep("par",names(KKK[[g]]))][1]))

		matureVKKTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
        		                                                               , par = VKK[[g]][grep("par",names(VKK[[g]]))][1:6])))
		matureVVKTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
                                                                      		   , par = VVK[[g]][grep("par",names(VVK[[g]]))][1:6])))
		matureVKVTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
                                                                      		   , par = VKV[[g]][grep("par",names(VKV[[g]]))][1:6])))
		matureVVVTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
                                                                      		   , par = VVV[[g]][grep("par",names(VVV[[g]]))][1:6])))

		modelKKK <- c(matureKKKTemp,prematureKKKTemp)
		modelVKK <- c(matureVKKTemp,prematureVKKTemp)
		modelVVK <- c(matureVVKTemp,prematureVVKTemp)
		modelVKV <- c(matureVKVTemp,prematureVKVTemp)
		modelVVV <- c(matureVVVTemp,prematureVVVTemp)

		KKKTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelKKK
		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
		VKKTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVKK
		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
		VVKTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVVK
		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
		VKVTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVKV
		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
		VVVTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVVV
		                               , variance = c(matureVariance[g,],prematureVariance[g,]))

		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)

		KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)

		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)

		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

	},BPPARAM=BPPARAM))

	AIC <- t(apply(logLikelihood,1,function(row)
	{
		2*(dof - row) 
	}))

	AICc <- t(apply(logLikelihood,1,function(row)
	{
	 	2*(dof - row) + (2*dof*(dof+1))/(2*length(tptsOriginal)-dof-1)
	}))

	rownames(pvaluesdata) <- rownames(logLikelihood) <- rownames(AIC) <- rownames(AICc) <- eiGenes

	ratesSpecs <- lapply(eiGenes,function(gene)
 		{
 			list(
 				"0" = list(alpha = list(fun = constantModelP
									 ,type = "constant"
									 ,df = 1
									 ,params = c(alpha = unname(KKK[[gene]]["par1"])))
						,beta = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(KKK[[gene]]["par3"])))
	 				   	,gamma = list(fun = constantModelP
									 ,type = "constant"
									 ,df = 1
									 ,params = c(gamma = unname(KKK[[gene]]["par2"])))
						,test = log(pvaluesdata[gene,"KKK"])
						,logLik = logLikelihood[gene,"KKK"]
						,AIC = AIC[gene,"KKK"]
						,AICc = AICc[gene,"KKK"]
						,counts = c("function"=unname(KKK[[gene]]["counts.function"]), gradient=unname(KKK[[gene]]["counts.gradient"]))
						,convergence = unname(KKK[[gene]]["convergence"])
						,message = NULL)

 				,a = list(alpha = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VKK[[gene]][1:6])))
 					  		,beta = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(VKK[[gene]][8])))
	 				  		,gamma = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(gamma = unname(VKK[[gene]][7])))
						,test = log(pvaluesdata[gene,"VKK"])
						,logLik = logLikelihood[gene,"VKK"]
						,AIC = AIC[gene,"VKK"]
						,AICc = AICc[gene,"VKK"]
						,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
						,convergence = unname(VKK[[gene]]["convergence"])
						,message = NULL)
 				,b = list(alpha = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KKV[[gene]][1])))
 					  		,beta = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(KKV[[gene]][3:8])))
	 				  		,gamma = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(gamma = unname(KKV[[gene]][2])))
						,test = log(pvaluesdata[gene,"KKV"])
						,logLik = logLikelihood[gene,"KKV"]
						,AIC = AIC[gene,"KKV"]
						,AICc = AICc[gene,"KKV"]
						,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
						,convergence = unname(KKV[[gene]]["convergence"])
						,message = NULL)
 				,c = list(alpha = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KVK[[gene]][1])))
 					  		,beta = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(KVK[[gene]][8])))
	 				  		,gamma = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(gamma = unname(KVK[[gene]][2:7])))
						,test = log(pvaluesdata[gene,"KVK"])
						,logLik = logLikelihood[gene,"KVK"]
						,AIC = AIC[gene,"KVK"]
						,AICc = AICc[gene,"KVK"]
						,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
						,convergence = unname(KVK[[gene]]["convergence"])
						,message = NULL)
 				,ab = list(alpha = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VKV[[gene]][1:6])))
 					  		,beta = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(VKV[[gene]][8:13])))
	 				  		,gamma = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(gamma = unname(VKV[[gene]][7])))
						,test = log(pvaluesdata[gene,"VKV"])
						,logLik = logLikelihood[gene,"VKV"]
						,AIC = AIC[gene,"VKV"]
						,AICc = AICc[gene,"VKV"]
						,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
						,convergence = unname(VKV[[gene]]["convergence"])
						,message = NULL)
 				,ac = list(alpha = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VVK[[gene]][1:6])))
 					  		,beta = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(VVK[[gene]][13])))
	 				  		,gamma = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(gamma = unname(VVK[[gene]][7:12])))
						,test = log(pvaluesdata[gene,"VVK"])
						,logLik = logLikelihood[gene,"VVK"]
						,AIC = AIC[gene,"VVK"]
						,AICc = AICc[gene,"VVK"]
						,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
						,convergence = unname(VVK[[gene]]["convergence"])
						,message = NULL)
 				,bc = list(alpha = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KVV[[gene]][1])))
 					  		,beta = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(KVV[[gene]][8:13])))
	 				  		,gamma = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(gamma = unname(KVV[[gene]][2:7])))
						,test = log(pvaluesdata[gene,"KVV"])
						,logLik = logLikelihood[gene,"KVV"]
						,AIC = AIC[gene,"KVV"]
						,AICc = AICc[gene,"KVV"]
						,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
						,convergence = unname(KVV[[gene]]["convergence"])
						,message = NULL)
 				,abc = list(alpha = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VVV[[gene]][1:6])))
 					  		,beta = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(VVV[[gene]][13:18])))
	 				  		,gamma = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(gamma = unname(VVV[[gene]][7:12])))
						,test = log(pvaluesdata[gene,"VVV"])
						,logLik = logLikelihood[gene,"VVV"]
						,AIC = AIC[gene,"VVV"]
						,AICc = AICc[gene,"VVV"]
						,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
						,convergence = unname(VVV[[gene]]["convergence"])
						,message = NULL)
			)
 		})

		return(ratesSpecs)
}

.makeModel_Derivative <- function(tpts, hyp, log_shift, time_transf, ode, .rxnrate, c= NaN, geneBestModel = NULL)
{

	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(time_transf(x, log_shift, c), hyp$alpha$par)
	params$beta  <- function(x) 
		hyp$beta$fun$value(time_transf(x, log_shift, c), hyp$beta$par)
	params$gamma <- function(x) 
		hyp$gamma$fun$value(time_transf(x, log_shift, c), hyp$gamma$par)

	matureTemp <- params$alpha(tpts)

	if(geneBestModel == "0")
	{
		prematureTemp <- sapply(tpts,function(t)prematureKKK_Int_NoNascent(t,c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params)))
		k1Temp <- sapply(tpts,function(t)k1KKK_Int_NoNascent(t,c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params)))
		
		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))
	
	}else if(geneBestModel == "a")
	{
		prematureTemp <- prematureVKK_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VKK_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))

	}else if(geneBestModel == "ac")
	{
		prematureTemp <- prematureVVK_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VVK_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$gamma$params)
		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))

	}else if(geneBestModel == "ab")
	{
		prematureTemp <- prematureVKV_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VKV_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
		k3Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$beta$params)

	}else if(geneBestModel == "abc")
	{
		prematureTemp <- prematureVVV_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VVV_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$gamma$params)
		k3Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$beta$params)

	}

	totalTemp <- matureTemp + prematureTemp

	data.frame(time = tpts, preMRNA = prematureTemp, total = totalTemp, alpha = k1Temp, beta = k3Temp, gamma = k2Temp)

}


.inspect.engine_Integrative_NoNascent <- function(tptsOriginal
										    , tptsLinear
										    , a
										    , c
										    , concentrations
										    , rates
										    , BPPARAM=bpparam()
										    , na.rm=TRUE
										    , verbose=TRUE
										    , testOnSmooth=TRUE
										    , seed = NULL
										    , nInit = nInit
										    , nIter = nIter)
{

	total <- concentrations$total
	totalVariance <- concentrations$total_var

	premature <- concentrations$preMRNA
	prematureVariance <- concentrations$preMRNA_var

	mature <- concentrations$mature
	matureVariance <- concentrations$mature_var

	eiGenes <- rownames(mature)

	k2median <- median(rates$gamma,na.rm = TRUE)
	k3median <- median(rates$beta,na.rm = TRUE)

	if( testOnSmooth ) {

		###### fit smooth functions on prematue and mature data
		###### eventually, exclude the genes that cannot be
		###### fit in either premature or mature

		matureFitImpulse <- bplapply(eiGenes,function(row)
		{
	  		fitSmooth(tpts = tptsLinear
	        	    , tt_c = c
	        	    , experiment = mature[row,]
	        	    , variance = matureVariance[row,]
	        	    , mature = TRUE
	        	    , nInit = nInit
	        	    , nIter = nIter
	        	    , seed = seed)
		},BPPARAM=BPPARAM)
		names(matureFitImpulse) <- eiGenes

		prematureFitImpulse <- bplapply(eiGenes,function(row)
		{
			fitSmooth(tpts = tptsLinear
					, tt_c = c
	        		, experiment = premature[row,]
	        		, variance = prematureVariance[row,]
	        		, mature = FALSE
	        		, nInit = nInit
	        		, nIter = nIter
	        		, seed = seed)     
		},BPPARAM=BPPARAM)
		names(prematureFitImpulse) <- eiGenes

		## in case the number of degrees of freedom allows the estimation
		## quantify how many genes have an acceptable chi2 test value (<0.2)
		## in both mature and premature
 		
 	# 	degreesOfFreedom <- length(tptsLinear)-6
		# if( degreesOfFreedom > 0 ) {
		# 	accept <- pchisq(sapply(matureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2 &
		# 		pchisq(sapply(prematureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2
		# 	if( table(accept)['FALSE']/length(accept) > .5 )
		# 		message("More than 50% did not return a good fit of their mature and
		# 			premature profiles with the impulsive smooth function:
		# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
		# }

		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
			stop("No genes have an expression profile possible to fit 
				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")

		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
			message("Some genes have an expression profile impossible to be fitted 
				with impulsive functions therefore they will be excluded from the modelling.")

			eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null) & 
				!sapply(prematureFitImpulse,is.null)]

			prematureFitImpulse <- prematureFitImpulse[eiGenes]
			matureFitImpulse <- matureFitImpulse[eiGenes]
			
			total <- total[eiGenes,]
			totalVariance <- totalVariance[eiGenes,]
			premature <- premature[eiGenes,]
			prematureVariance <- prematureVariance[eiGenes,]
			mature <- mature[eiGenes,]
			matureVariance <- matureVariance[eiGenes,]
		}


	}

	KKK <- bplapply(eiGenes,function(row){
	
			matureParameters <- mean(mature[row,])
			k2Parameters <- k2median
			k3Parameters <- k3median

			unlist(
				tryCatch(
					optim(c(matureParameters, k2Parameters, k3Parameters)
						,errorKKK_Int_NoNascent
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,control = list(maxit = nIter)),
					error=function(e) c(par1 = NaN
									  , par2 = NaN
									  , par3 = NaN
									  , value = NaN
									  , counts.function = NaN
								  	  , counts.gradient = NaN
								  	  , convergence = e)
				)
			)
		}, BPPARAM=BPPARAM)
		names(KKK) <- eiGenes
		message("Model 0 finished.")
	
	VKK <- bplapply(eiGenes,function(row){

		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k2Parameters <- KKK[[row]][2]
		k3Parameters <- KKK[[row]][3]
  
		unlist(
			tryCatch(
      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
                  	 ,errorVKK_Int_NoNascent
                  	 ,times = tptsOriginal
                  	 ,data = c( mature[row,], premature[row,] )
                  	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
                  	 ,a = a
                  	 ,c = c
                  	 ,control = list(maxit = nIter)),
    		error=function(e) c(par1 = NaN
    						  , par2 = NaN
    						  , par3 = NaN
    						  , par4 = NaN
    						  , par5 = NaN
    						  , par6 = NaN
    						  , par7 = NaN
    						  , par8 = NaN
    						  , value = NaN
    						  , counts.function = NaN
    						  , counts.gradient = NaN
    						  , convergence = e)
    		)
    	)
	}, BPPARAM=BPPARAM)
	names(VKK) <- eiGenes
	message("Model A finished.")

	KKV <- bplapply(eiGenes, function(row){

		k1Parameters <- KKK[[row]][1]
		k2Parameters <- KKK[[row]][2]
		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
		        	,errorKKV_Int_NoNascent
		        	,times = tptsOriginal
		        	,data = c( mature[row,], premature[row,] )
		        	,datavar = c( matureVariance[row,] , prematureVariance[row,] )
		        	,a = a
		        	,c = c
		        	,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KKV) <- eiGenes
	message("Model B finished.")

	KVK <- bplapply(eiGenes, function(row){

		k1Parameters <- KKK[[row]][1]
		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k3Parameters <- KKK[[row]][3]
		
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVK_Int_NoNascent
					 ,times = tptsOriginal
					 ,data = c( mature[row,], premature[row,] )
					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
					 ,a = a
					 ,c = c
					 ,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KVK) <- eiGenes
	message("Model C finished.")

	VKV <- bplapply(eiGenes, function(row){

		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k2Parameters <- KKK[[row]][2]
		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorVKV_Int_NoNascent
					 ,times = tptsOriginal
					 ,data = c( mature[row,], premature[row,] )
					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
					 ,a = a
					 ,c = c
					 ,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
		                       ,par9 = NaN
							   ,par10 = NaN
							   ,par11 = NaN
							   ,par12 = NaN
							   ,par13 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(VKV) <- eiGenes
	message("Model AB finished.")

	VVK <- bplapply(eiGenes, function(row){

		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k3Parameters <- KKK[[row]][3]

		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
    				 ,errorVVK_Int_NoNascent
    				 ,times = tptsOriginal
    				 ,data = c( mature[row,], premature[row,] )
    				 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
    				 ,a = a
    				 ,c = c
    				 ,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							  , par2 = NaN
							  , par3 = NaN
							  , par4 = NaN
							  , par5 = NaN
							  , par6 = NaN
							  , par7 = NaN
							  , par8 = NaN
                			  , par9 = NaN
                			  , par10 = NaN
                			  , par11 = NaN
                			  , par12 = NaN
                			  , par13 = NaN
                			  , value = NaN
                			  , counts.function = NaN
                			  , counts.gradient = NaN
                			  , convergence = NaN)
			)
			)
	}, BPPARAM=BPPARAM)
	names(VVK) <- eiGenes
	message("Model AC finished.")

	KVV <- bplapply(eiGenes, function(row){

		k1Parameters <- KKK[[row]][1]
		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVV_Int_NoNascent
					 ,times = tptsOriginal
					 ,data = c( mature[row,], premature[row,] )
					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
					 ,a = a
					 ,c = c
					 ,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,par9 = NaN
							   ,par10 = NaN
							   ,par11 = NaN
							   ,par12 = NaN
							   ,par13 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(KVV) <- eiGenes
	message("Model BC finished.")

	VVV <- bplapply(eiGenes, function(row){

		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					,errorVVV_Int_NoNascent
					,times = tptsOriginal
					,data = c( mature[row,], premature[row,] )
					,datavar = c( matureVariance[row,] , prematureVariance[row,] )
					,a = a
					,c = c
					,control = list(maxit = nIter)),
			error=function(e) c(par1 = NaN
							   ,par2 = NaN
							   ,par3 = NaN
							   ,par4 = NaN
							   ,par5 = NaN
							   ,par6 = NaN
							   ,par7 = NaN
							   ,par8 = NaN
							   ,par9 = NaN
							   ,par10 = NaN
							   ,par11 = NaN
							   ,par12 = NaN
							   ,par13 = NaN
							   ,par14 = NaN
							   ,par15 = NaN
							   ,par16 = NaN
							   ,par17 = NaN
							   ,par18 = NaN
							   ,value = NaN
							   ,counts.function = NaN
							   ,counts.gradient = NaN
							   ,convergence = NaN)
			)
		)
	}, BPPARAM=BPPARAM)
	names(VVV) <- eiGenes
	message("Model ABC finished.")

	if(testOnSmooth)
	{
			matureSmooth <- t(sapply(matureFitImpulse,function(i)
			{
				impulseModel(tptsLinear,i)
			}))
		
			prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
			{
		  		impulseModel(tptsLinear,i)
			}))

			colnames(matureSmooth) <- tptsOriginal
			colnames(prematureSmooth) <- tptsOriginal

	}else{
		matureSmooth <- mature
		prematureSmooth <- premature
	}

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
											  ,tptsOriginal
											  ,prematureSmooth[g,]
											  ,matureSmooth[g,]
											  ,prematureVariance[g,]
											  ,matureVariance[g,]),error = function(e)NaN)

		VKKTemp <- tryCatch(errorVKK_Int_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		VVKTemp <- tryCatch(errorVVK_Int_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		VKVTemp <- tryCatch(errorVKV_Int_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		VVVTemp <- tryCatch(errorVVV_Int_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)
	  
		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
	}, BPPARAM = BPPARAM))

	dof <- c(KKK = 3
			,VKK = 8
			,KVK = 8
			,KKV = 8
			,VVK = 13
			,VKV = 13
			,KVV = 13
			,VVV = 18)

	# P values
	pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], max(c(0,2*length(tptsOriginal)-dof['KKK'])))
						,VKK=pchisq(chi2data[,'VKK'], max(c(0,2*length(tptsOriginal)-dof['VKK'])))
						,KVK=pchisq(chi2data[,'KVK'], max(c(0,2*length(tptsOriginal)-dof['KVK'])))
						,KKV=pchisq(chi2data[,'KKV'], max(c(0,2*length(tptsOriginal)-dof['KKV'])))
						,VVK=pchisq(chi2data[,'VVK'], max(c(0,2*length(tptsOriginal)-dof['VVK'])))
						,VKV=pchisq(chi2data[,'VKV'], max(c(0,2*length(tptsOriginal)-dof['VKV'])))
						,KVV=pchisq(chi2data[,'KVV'], max(c(0,2*length(tptsOriginal)-dof['KVV'])))
						,VVV=pchisq(chi2data[,'VVV'], max(c(0,2*length(tptsOriginal)-dof['VVV']))))

	logLikelihood <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- tryCatch(loglikKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
										  ,tptsOriginal,prematureSmooth[g,]
										  ,matureSmooth[g,]
										  ,prematureVariance[g,]
										  ,matureVariance[g,]),error = function(e)NaN)

		VKKTemp <- tryCatch(loglikVKK_Int_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
										  ,tptsOriginal
										  ,c(matureSmooth[g,],prematureSmooth[g,])
										  ,c(matureVariance[g,],prematureVariance[g,])
										  ,a
										  ,c),error = function(e)NaN)

	  	KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
        				                  ,tptsOriginal
        				                  ,c(matureSmooth[g,],prematureSmooth[g,])
        				                  ,c(matureVariance[g,],prematureVariance[g,])
        				                  ,a
        				                  ,c),error = function(e)NaN)

		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
											   ,tptsOriginal
											   ,c(matureSmooth[g,],prematureSmooth[g,])
											   ,c(matureVariance[g,],prematureVariance[g,])
											   ,a
											   ,c),error = function(e)NaN)
		VVKTemp <- tryCatch(loglikVVK_Int_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
											   ,tptsOriginal
											   ,c(matureSmooth[g,],prematureSmooth[g,])
											   ,c(matureVariance[g,],prematureVariance[g,])
											   ,a
											   ,c),error = function(e)NaN)
		VKVTemp <- tryCatch(loglikVKV_Int_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
											   ,tptsOriginal
											   ,c(matureSmooth[g,],prematureSmooth[g,])
											   ,c(matureVariance[g,],prematureVariance[g,])
											   ,a
											   ,c),error = function(e)NaN)
		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
											   ,tptsOriginal
											   ,c(matureSmooth[g,],prematureSmooth[g,])
											   ,c(matureVariance[g,],prematureVariance[g,])
											   ,a
											   ,c),error = function(e)NaN)
		VVVTemp <- tryCatch(loglikVVV_Int_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
											   ,tptsOriginal
											   ,c(matureSmooth[g,],prematureSmooth[g,])
											   ,c(matureVariance[g,],prematureVariance[g,])
											   ,a
											   ,c),error = function(e)NaN)

		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

	}, BPPARAM=BPPARAM))

		AIC <- t(apply(logLikelihood,1,function(row)
		{
			2*(dof - row) 
		}))

		AICc <- t(apply(logLikelihood,1,function(row)
		{
		 	2*(dof - row) + (2*dof*(dof+1))/(2*length(tptsOriginal)-dof-1)
		}))

	rownames(pvaluesdata) <- rownames(logLikelihood) <- rownames(AIC) <- rownames(AICc) <- eiGenes

		ratesSpecs <- lapply(eiGenes,function(gene)
		{
			list(
				"0" = list(alpha = list(fun = constantModelP
								 ,type = "constant"
								 ,df = 1
								 ,params = c(alpha = k1KKK_Int_NoNascent(x = 0,par = KKK[[gene]])))
					    ,beta = list(fun = constantModelP
								,type = "constant"
								,df = 1
								,params = c(beta = unname(KKK[[gene]]["par3"])))
 				   	,gamma = list(fun = constantModelP
								 ,type = "constant"
								 ,df = 1
								 ,params = c(gamma = unname(KKK[[gene]]["par2"])))
					,test = log(pvaluesdata[gene,"KKK"])
					,logLik = logLikelihood[gene,"KKK"]
					,AIC = AIC[gene,"KKK"]
					,AICc = AICc[gene,"KKK"]
					,counts = c("function"=unname(KKK[[gene]]["counts.function"]), gradient=unname(KKK[[gene]]["counts.gradient"]))
					,convergence = unname(KKK[[gene]]["convergence"])
					,message = NULL)

				,a = list(alpha = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(alpha = unname(VKK[[gene]][1:6])))
					  		,beta = list(fun = constantModelP
								,type = "constant"
								,df = 1
								,params = c(beta = unname(VKK[[gene]][8])))
 				  		,gamma = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(gamma = unname(VKK[[gene]][7])))
					,test = log(pvaluesdata[gene,"VKK"])
					,logLik = logLikelihood[gene,"VKK"]
					,AIC = AIC[gene,"VKK"]
					,AICc = AICc[gene,"VKK"]
					,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
					,convergence = unname(VKK[[gene]]["convergence"])
					,message = NULL)
				,b = list(alpha = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(alpha = unname(KKV[[gene]][1])))
					  		,beta = list(fun = impulseModelP
								,type = "impulse"
								,df = 6
								,params = c(beta = unname(KKV[[gene]][3:8])))
 				  		,gamma = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(gamma = unname(KKV[[gene]][2])))
					,test = log(pvaluesdata[gene,"KKV"])
					,logLik = logLikelihood[gene,"KKV"]
					,AIC = AIC[gene,"KKV"]
					,AICc = AICc[gene,"KKV"]
					,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
					,convergence = unname(KKV[[gene]]["convergence"])
					,message = NULL)
				,c = list(alpha = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(alpha = unname(KVK[[gene]][1])))
					  		,beta = list(fun = constantModelP
								,type = "constant"
								,df = 1
								,params = c(beta = unname(KVK[[gene]][8])))
 				  		,gamma = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(gamma = unname(KVK[[gene]][2:7])))
					,test = log(pvaluesdata[gene,"KVK"])
					,logLik = logLikelihood[gene,"KVK"]
					,AIC = AIC[gene,"KVK"]
					,AICc = AICc[gene,"KVK"]
					,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
					,convergence = unname(KVK[[gene]]["convergence"])
					,message = NULL)
				,ab = list(alpha = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(alpha = unname(VKV[[gene]][1:6])))
					  		,beta = list(fun = impulseModelP
								,type = "impulse"
								,df = 6
								,params = c(beta = unname(VKV[[gene]][8:13])))
 				  		,gamma = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(gamma = unname(VKV[[gene]][7])))
					,test = log(pvaluesdata[gene,"VKV"])
					,logLik = logLikelihood[gene,"VKV"]
					,AIC = AIC[gene,"VKV"]
					,AICc = AICc[gene,"VKV"]
					,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
					,convergence = unname(VKV[[gene]]["convergence"])
					,message = NULL)
				,ac = list(alpha = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(alpha = unname(VVK[[gene]][1:6])))
					  		,beta = list(fun = constantModelP
								,type = "constant"
								,df = 1
								,params = c(beta = unname(VVK[[gene]][13])))
 				  		,gamma = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(gamma = unname(VVK[[gene]][7:12])))
					,test = log(pvaluesdata[gene,"VVK"])
					,logLik = logLikelihood[gene,"VVK"]
					,AIC = AIC[gene,"VVK"]
					,AICc = AICc[gene,"VVK"]
					,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
					,convergence = unname(VVK[[gene]]["convergence"])
					,message = NULL)
				,bc = list(alpha = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(alpha = unname(KVV[[gene]][1])))
					  		,beta = list(fun = impulseModelP
								,type = "impulse"
								,df = 6
								,params = c(beta = unname(KVV[[gene]][8:13])))
 				  		,gamma = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(gamma = unname(KVV[[gene]][2:7])))
					,test = log(pvaluesdata[gene,"KVV"])
					,logLik = logLikelihood[gene,"KVV"]
					,AIC = AIC[gene,"KVV"]
					,AICc = AICc[gene,"KVV"]
					,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
					,convergence = unname(KVV[[gene]]["convergence"])
					,message = NULL)
				,abc = list(alpha = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(alpha = unname(VVV[[gene]][1:6])))
					  		,beta = list(fun = impulseModelP
								,type = "impulse"
								,df = 6
								,params = c(beta = unname(VVV[[gene]][13:18])))
 				  		,gamma = list(fun = impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(gamma = unname(VVV[[gene]][7:12])))
					,test = log(pvaluesdata[gene,"VVV"])
					,logLik = logLikelihood[gene,"VVV"]
					,AIC = AIC[gene,"VVV"]
					,AICc = AICc[gene,"VVV"]
					,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
					,convergence = unname(VVV[[gene]]["convergence"])
					,message = NULL)
		)
		})

	return(ratesSpecs)
}

inferKBetaFromIntegralWithPre <- function(tpts, alpha, total, preMRNA, maxBeta=75, BPPARAM=bpparam()) 
####### accurate function for estimating the degradation rates
####### using the solution of the differential equation system under 
####### the condtion that degradation rate is constant between two 
####### consecutive time points - more stable that using derivatives
####### estimates
{
	solveBeta <- function(beta, t0, t1, alpha_t0, alpha_t1, X_t0, X_t1, P_t0, P_t1 ) 
	{
		mAlpha <- (alpha_t0 - alpha_t1 ) / (t0 - t1 )
		qAlpha <- alpha_t0 - mAlpha * t0
		#
		mPreMRNA <- (P_t0 - P_t1 ) / (t0 - t1 )
		qPreMRNA <- P_t0 - mPreMRNA * t0
		#
		X_t1 - X_t0 * exp(-beta*(t1-t0)) - 
		((mAlpha*t1*beta + qAlpha*beta - mAlpha ) / (beta^2 ) - (mAlpha*t0*beta + qAlpha*beta - mAlpha ) * exp(-beta*(t1-t0)) / (beta^2 )) -
		beta*((mPreMRNA*t1*beta + qPreMRNA*beta - mPreMRNA ) / (beta^2 ) - (mPreMRNA*t0*beta + qPreMRNA*beta - mPreMRNA ) * exp(-beta*(t1-t0)) / (beta^2 ))
	}
	bplapply(2:length(tpts), function(j)
	lapply(1:nrow(total), function(i) {
	tryCatch(
		uniroot(solveBeta
		, c(1e-5, maxBeta)
		, t0 = tpts[j-1]
		, t1 = tpts[j]
		, alpha_t0 = alpha[i,j-1]
		, alpha_t1 = alpha[i,j]
		, X_t0 = total[i,j-1]
		, X_t1 = total[i,j]
		, P_t0 = preMRNA[i,j-1]
		, P_t1 = preMRNA[i,j]
		)
		, error=function(e) return(list(root=NA, estim.prec=NA, error=e))
	)})
	, BPPARAM=BPPARAM)
}

counts2expressions <- function(counts, widths, libsize) counts*(10^9/(widths[rownames(counts)]%o%libsize))
countVar2expressions <- function(vars, widths, libsize) vars*(10^9/(widths%o%libsize))^2
inferKGammaFromIntegral <- function(tpts, alpha, preMRNA, maxGamma=150, BPPARAM=bpparam())
####### accurate function for estimating the degradation rates
####### using the solution of the differential equation system under 
####### the condtion that processing rate is constant between two 
####### consecutive time points - more stable that using derivatives
####### estimates
{
	solveFun <- function(beta, t0, t1, alpha_t0, alpha_t1, X_t0, X_t1 ) 
	{
		m <- (alpha_t0 - alpha_t1 ) / (t0 - t1 )
		q <- alpha_t0 - m * t0
		X_t1 - X_t0*exp(-beta*(t1 - t0)) - (
			(m*t1*beta + q*beta - m ) / (beta^2) - 
			(m*t0*beta + q*beta - m ) * exp(-beta*(t1 - t0)) / (beta^2)
			)
	}
	bplapply(2:length(tpts), function(j)
		lapply(1:nrow(alpha), function(i) {
			tryCatch(
				uniroot(solveFun
					, c(1e-5, maxGamma)
					, t0 = tpts[j-1]
					, t1 = tpts[j]
					, alpha_t0 = alpha[i,j-1]
					, alpha_t1 = alpha[i,j]
					, X_t0 = preMRNA[i,j-1]
					, X_t1 = preMRNA[i,j]
					)
			, error=function(e) return(list(root=NA, estim.prec=NA, error=e))
		)})
		, BPPARAM=BPPARAM)
}

impute_na_tc <- function(tpts, tcdata) {

	# impute NA values in a time course data as a linear model between non-NA values
	tc_impute_NA_linearmodel <- function(tpts, tcdata) {
		for( j in seq_along(tcdata) ) {
			if( is.na(tcdata[j]) ) {
				lower_boundary_j <- j-1
				higher_boundary_j <- j+1
				while( is.na(tcdata[higher_boundary_j]) & higher_boundary_j <= length(tcdata) ) {
					higher_boundary_j <- higher_boundary_j + 1
				}
				if( lower_boundary_j > 0 & higher_boundary_j <= length(tcdata) ) 
					if ( is.finite(tcdata[lower_boundary_j]) ) {
						x <- tpts[c(lower_boundary_j,higher_boundary_j)]
						y <- tcdata[c(lower_boundary_j,higher_boundary_j)]
						tcdata[(lower_boundary_j+1):(higher_boundary_j-1)] <- predict(lm(y ~ x), 
							data.frame(x=tpts[(lower_boundary_j+1):(higher_boundary_j-1)]))
					}
			}
		}
		return(tcdata)
	}

	# impute NA values in a time course from boundary values
	tc_impute_NA_boundaries <- function(tcdata) {

		forward_direction <- function(tcdata) {
			if( is.na(tcdata[1]) & !all(is.na(tcdata)) ) {
				lower_boundary_j <- higher_boundary_j <- 1
				while( is.na(tcdata[higher_boundary_j] & higher_boundary_j < length(tcdata) ) ) {
					higher_boundary_j <- higher_boundary_j + 1
				}
				tcdata[lower_boundary_j:(higher_boundary_j-1)] <- tcdata[higher_boundary_j]
			}
			return(tcdata)
		}

		tcdata <- forward_direction(tcdata)
		tcdata <- rev(forward_direction(rev(tcdata)))
		return(tcdata)
	}

	tcdata <- tc_impute_NA_linearmodel(tpts, tcdata)
	tcdata <- tc_impute_NA_boundaries(tcdata)
	return(tcdata)

}


########## Compare steady state no nascent

standardCurveFitFunction <- function(p,m,err)
{
	n_outliers <- function(alpha, x, y, err) {

		#Conversion
		pi_angle <- alpha * pi/180
		#Angular coefficient
		coef_ang <- tan(pi_angle)

		delta_intercept <- err/cos(pi_angle)
		intercept <- median(y,na.rm=TRUE) - coef_ang*median(x,na.rm=TRUE)

		outliers <- y > coef_ang * x + intercept + delta_intercept |
			y < coef_ang * x + intercept - delta_intercept

		length(which(outliers))
	}

	all_alphas <- seq(-89,90)
	all_alphas_outliers <- sapply(all_alphas, function(i) n_outliers(alpha = i, x=log2(p), y=log2(m), err = err))
	return(seq(-89,90)[which.min(all_alphas_outliers)])
}

classificationFunction <- function(p,m,alpha,err)
{
	standardCurveFit <- alpha
	classificationTmp <- sapply(rownames(p),function(g)
	{
		x <- log2(p[g,])
		y <- log2(m[g,])

		pi_angle <- standardCurveFit * pi/180
		coef_ang <- tan(pi_angle)
		delta_intercept <- err/cos(pi_angle)

		intercept <- median(y,na.rm=TRUE) - coef_ang*median(x,na.rm=TRUE)

		outliers <- y > coef_ang * x + intercept + delta_intercept |
				y < coef_ang * x + intercept - delta_intercept
		return(outliers)
	})

	return(t(classificationTmp))
}

# plotRegressionCurve <- function(premature, mature, alpha, err, main, xlimU, ylimU, smooth=TRUE, outliers=FALSE)
# {
# 	if(is.matrix(premature)&is.matrix(mature))
# 	{
# 		x <- log2(apply(premature,1,median,na.rm=TRUE))
# 		y <- log2(apply(mature,1,median,na.rm=TRUE))
# 	}else{
# 		x <- log2(premature)
# 		y <- log2(mature)
# 	}

# 	pi_angle <- alpha * pi/180
# 	coef_ang <- tan(pi_angle)
# 	delta_intercept <- err/cos(pi_angle)

# 	intercept <- median(y,na.rm=TRUE) - coef_ang*median(x,na.rm=TRUE)

# 	if(smooth){
# 		smoothScatter(x,y,xlab="Log2 median premature RNA",ylab="Log2 median mature RNA",main=main)

# 		abline(intercept,coef_ang,col=2,lw=3)
# 		abline(intercept + delta_intercept,coef_ang,col=2,lw=3,lty=2)	
# 		abline(intercept - delta_intercept,coef_ang,col=2,lw=3,lty=2)


# 	}else{
# 		if(!outliers)
# 		{
# 			df <- data.frame(x = x
# 						   , y = y
# 						   , d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
# 						   , cyl = rep(15,length(x)))
	
# 			p <- ggplot(df) + xlim(xlimU) + ylim(ylimU) +
#     		geom_point(aes(x, y, col = d), size = 2, shape=df$cyl) +
#     		scale_color_identity() +
#     		theme_bw() + labs(x = "Log2 premature RNA", y = "Log2 mature RNA", title = main)+
#     		geom_abline(intercept = intercept, slope = coef_ang, color="red", size=1.5)+
#     		geom_abline(intercept = intercept+ delta_intercept,linetype="dashed", slope = coef_ang, color="red", size=1.5)+
#     		geom_abline(intercept = intercept- delta_intercept,linetype="dashed", slope = coef_ang, color="red", size=1.5)
# 			message(p)
# 		}
# 		if(outliers)
# 		{
# 			x <- x[is.finite(x)&is.finite(y)]
# 			y <- y[names(x)]

# 			df <- data.frame(x = x, y = y)

# 			boolTmp <- apply(as.matrix(df),1,function(r)
# 			{
# 				as.numeric(r[1])*coef_ang + (intercept + delta_intercept) > as.numeric(r[2]) & as.numeric(r[1]) * coef_ang + (intercept- delta_intercept) < as.numeric(r[2])
# 			})	

# 			dfT <- data.frame(x=df$x[boolTmp],y=df$y[boolTmp], col =  densCols(df$x[boolTmp], df$y[boolTmp], colramp = colorRampPalette(rev(c("gray46","gray88")))))
# 			dfF <- data.frame(x=df$x[!boolTmp],y=df$y[!boolTmp], col =  densCols(df$x[!boolTmp], df$y[!boolTmp], colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

# 			p <- ggplot()
# 			p <- p + xlim(xlimU) + ylim(ylimU)
# 			p <- p + geom_point(aes(dfT$x, dfT$y, col = dfT$col), size = 2, shape=15)
# 			p <- p + geom_point(aes(dfF$x, dfF$y, col = dfF$col), size = 2, shape=15)
#     		p <- p + scale_color_identity()
#      		p <- p + labs(x = "Log2 premature RNA", y = "Log2 mature RNA", title = main) +
#     			geom_abline(intercept = intercept, slope = coef_ang, color="red", size=1.5)+
#     			geom_abline(intercept = intercept+ delta_intercept,linetype="dashed", slope = coef_ang, color="red", size=1.5)+
#     			geom_abline(intercept = intercept- delta_intercept,linetype="dashed", slope = coef_ang, color="red", size=1.5)
# 			p <- p + theme_light()
# 			message(p)
# 		}
# 	}
# }
