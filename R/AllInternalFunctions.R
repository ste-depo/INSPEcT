grepLogic <- function(text,obj,...){seq_along(obj)%in%grep(text,obj,...)}

# fromImpulseToSigmoid <- function(impulsesParameters,tpts,a,c,nIter)
# {
# 	impulseProfile <- impulseModel(x=log2(tpts+a)+c,par=impulsesParameters)

# 	internalError <- function(sigmoidsParameters,impulsesParameters,tpts,a,c)
# 	{
# 		impulseProfile <- impulseModel(x=log2(tpts+a)+c,par=impulsesParameters)
# 		sigmoidProfile <- sigmoidModel(x=log2(tpts+a)+c,par=sigmoidsParameters)

# 		chisqFunction(experiment=impulseProfile,model=sigmoidProfile,variance=impulseProfile)
# 	}

# 	optTmp <- optim(par=c(head(impulseProfile,1),tail(impulseProfile,1),max(log2(tpts+a)+c)/3,1)
# 		,internalError
# 		,impulsesParameters=impulsesParameters
# 		,tpts=tpts
# 		,a=a
# 		,c=c
# 		,control = list(maxit = nIter))

# 	return(optTmp$par)
# }

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

time_transf <- function(t, log_shift, lin_shift = 0) 
{
	t[ t <= (-log_shift) ] <- NaN
	newtime <- log2(t+log_shift) + lin_shift
	return(newtime)
} 

time_transf_inv <- function(t, log_shift, lin_shift=0) 2^(t - lin_shift) - log_shift

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

.makeModel <- function(tpts, hyp, nascent = FALSE)
{
	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(x, hyp$alpha$par)
	params$beta  <- function(x) 
		hyp$beta$fun$value(x, hyp$beta$par)
	params$gamma <- function(x) 
		hyp$gamma$fun$value(x, hyp$gamma$par)

	cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1]), 
			   params$alpha(tpts[1]) / params$beta(tpts[1]))

	if(nascent){cinit=c(0,0)}
	names(cinit) <- c('p', 'm')
	model <- as.data.frame(
		ode(y=cinit, times=tpts, func=.rxnrate, parms=params))

	model$alpha <- params$alpha(tpts)
	model$beta  <- params$beta(tpts)
	model$gamma <- params$gamma(tpts)
	colnames(model)[2:3] <- c('preMRNA','mature')
	model[,3] <- apply(model[2:3],1,sum)
	colnames(model)[2:3] <- c('preMRNA','total')

	if(nrow(model)!=length(tpts)){return(matrix(rep(NaN,2*length(tpts)),nrow=tpts,ncol=2))}
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
	r[2] <- gamma(t) * c["p"] - beta(t) * c["m"]

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
	, verbose=TRUE, limitModelComplexity=FALSE, estimateRatesWith=c('der', 'int'), nAttempts=1
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
		, .emptyGene, limitModelComplexity)
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
		if( limitModelComplexity ) {
			k <- min(interpRates$alpha$df, length(tpts_exp)) + min(interpRates$beta$df, length(tpts_exp))
		} else {
			k <- interpRates$alpha$df + interpRates$beta$df
		}
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
		, total_var, preMRNA_exp, preMRNA_var, maxit=500, log_shift, time_transf, 
		.rxnrate, ode, .makeModel, logLikelihoodFunction, limitModelComplexity)
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
			return(D)
			# df <- length(alpha_exp) + length(total_exp) + length(preMRNA_exp) - sum(df)
			# testValue <- chisq.test.inspect(D, df)
			# return(testValue)
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

		if( limitModelComplexity ) {
			k <- min(interpRates$alpha$df, length(tpts_exp)) + min(interpRates$beta$df, length(tpts_exp)) + 
				min(interpRates$gamma$df, length(tpts_exp))
		} else {
			k <- interpRates$alpha$df + interpRates$beta$df + interpRates$gamma$df
		}
		n <- length(alpha_exp) + length(total_exp) + length(preMRNA_exp)
		# chisqTest <- log(optOut$value)
		chisqTest <- log(pchisq(optOut$value, n-k, lower.tail=TRUE))
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
			sigmoidSynthesis, nInit, nIter, testOnSmooth, estimateRatesWith, limitModelComplexity) 
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
								, limitModelComplexity=limitModelComplexity
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
								, limitModelComplexity=limitModelComplexity
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
		## choose the best model for each test, out of the many attempts
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
		limitModelComplexity=limitModelComplexity,
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

	.D2sigmoidModel <- function(x, par) {
		h0= par[1]; h1=par[2]; t1=par[3]; b=par[4]
		(2*b^2*(h1-h0)*exp(-2*b*(x-t1)))/(exp(-b*(x-t1))+1)^3-(b^2*(h1-h0)*exp(-b*(x-t1)))/(exp(-b*(x-t1))+1)^2
	}

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

# .chooseModel <- function(tpts, log_shift, experiment, variance=NULL, na.rm=TRUE
# 	, sigmoid=TRUE, impulse=TRUE, polynomial=TRUE, nInit=10, nIter=500
# 	, time_transf, impulseModel, sigmoidModel, sigmoidModelP, impulseModelP
# 	, .polynomialModelP)
# #### choose a functional form between impulse and sigmoid according 
# #### to the one that has the gratest pvalue in the chi squared test
# {
# 	chisq.test.default <- function(experiment, model, variance=NULL, df)
# 	{
# 		if( is.null(variance) ) variance <- stats::var(experiment )
# 		D = chisqFunction(experiment, model, variance)
# 		modelDF <- max(0, length(experiment)-df)
# 		pchisq(D,  modelDF, lower.tail=TRUE)
# 	}

# 	optimFailOut <- function(e)
# 		list(par=NA, value=NA, counts=NA, convergence=1, message=e)
# 	#
# 	# impulse model functions
# 	#
# 	im.parguess <- function(tpts , values ) {
# 	    # values = expressions.avgd(eD)
# 	    # tp = tpts(eD)
# 	    ntp   <- length(tpts)
# 	    peaks <- which(diff(sign(diff(values)))!=0)+1
# 	    if( length(peaks) == 1 ) peak <- peaks
# 	    if( length(peaks)  > 1 ) peak <- sample(peaks, 1)
# 	    if( length(peaks) == 0 ) peak <- round(length(tpts)/2)
# 	    #
# 		initial_values <- runif( 1, min=min(values[1:3])
# 			, max=max(values[1:3]))
# 		intermediate_values <- values[peak]
# 		if( intermediate_values==0 ) intermediate_values <- mean(values[seq(peak-1,peak+1)])
# 		end_values <- runif( 1, min=min(values[(ntp-2):ntp])
# 			, max=max(values[(ntp-2):ntp]))
# 		time_of_first_response  <- tpts[peak-1]
# 		time_of_second_response <- tpts[peak+1]
# 		slope_of_response <- diff(range(tpts)) / 
# 			(time_of_second_response-time_of_first_response)
# 		#
# 	    return(c(h0=initial_values, h1=intermediate_values
# 	    	, h2=end_values, t1=time_of_first_response
# 	    	, t2=time_of_second_response, b=slope_of_response))
# 	}
# 	#
# 	im.chisq <- function(par, tpts, experiment, variance=NULL, impulseModel) 
# 	{
# 		 model <- impulseModel(tpts, par)
# 		 chisqFunction(experiment, model, variance)
# 	}
# 	#
# 	im.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
# 		, maxit=500) 
# 		sapply(1:ninit, function(x) 
# 	 		tryCatch(optim(
# 	 			par=im.parguess(tpts, experiment)
# 	 			, fn=im.chisq, tpts=tpts
# 	 			, experiment=experiment
# 	 			, variance=variance
# 	 			, impulseModel=impulseModel
# 	 			, control=list(maxit=maxit)
# 	 			), error=function(e) optimFailOut(e)))
# 	#
# 	# sigmoid model functions
# 	#
# 	sm.parguess <- function(tpts , values ) {
# 	    # values = expressions.avgd(eD)
# 	    # tp = tpts(eD)

# 		time_span <- diff(range(tpts))
# 		# sample the time uniformely
# 		time_of_response <- runif( 1, min=min(tpts), max=max(tpts))
# 		# slope of response must be high if the time of response is close to one
# 		# of the two boundaries
# 		distance_from_boundary <- min(time_of_response - min(tpts)
# 				, max(tpts) - time_of_response)
# 		slope_of_response <- time_span / distance_from_boundary
# 	    ntp   <- length(tpts)
# 		initial_values <- runif( 1, min=min(values[1:3])
# 			, max=max(values[1:3]))
# 		end_values <- runif( 1, min=min(values[(ntp-2):ntp])
# 			, max=max(values[(ntp-2):ntp]))
# 		#
# 	    return(c(h0=initial_values, h1=end_values, t1=time_of_response
# 	    	, b=slope_of_response))
# 	}
# 	#
# 	sm.chisq <- function(par, tpts, experiment, variance=NULL, sigmoidModel) 
# 	{
# 		 model <- sigmoidModel(tpts, par)
# 		 chisqFunction(experiment, model, variance)
# 	}
# 	#
# 	sm.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
# 		, maxit=500) 
# 		sapply(1:ninit, function(x) 
# 				tryCatch(optim(
# 					par=sm.parguess(tpts, experiment)
# 					, fn=sm.chisq, tpts=tpts
# 					, experiment=experiment
# 					, variance=variance
# 					, sigmoidModel=sigmoidModel
# 					, control=list(maxit=maxit)
# 					), error=function(e) optimFailOut(e))) 

# 	pn.optim.aic <- function(tpts , experiment, variance=NULL) {
# 		if( length(experiment) < 3 ) return(NA)
# 		polyOrderChisq <- function(i) {
# 			model <- lm(experiment~poly(tpts, i, raw=TRUE ))
# 			return(list(par=model$coeff, value=AIC(model)))}
# 		return(sapply(1:min(7,length(tpts)-2), polyOrderChisq))
# 	}

# 	# remove missing values
# 	if( na.rm) {
# 		idx <- is.finite(experiment)
# 		tpts <- tpts[idx]
# 		experiment <- experiment[idx]
# 	}

# 	## 
# 	if( length(experiment)==0 ) {
# 		stop('.chooseModel: no time points have a finite value.
# 			Impossible to evaluate any kind of model.')
# 		return(list(type='constant', fun=constantModelP
# 			, params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
# 	}
# 	## 
# 	if( length(experiment)<=2 ) {
# 		warning('.chooseModel: less than three time points have a finite value. 
# 			Impossible evaluate a variable model.
# 			Returning a constant model.')
# 		return(list(type='constant', fun=constantModelP
# 			, params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
# 	}

# 	## re-evaluate flags of function to evaluate according to the lenght 
# 	## of the experiment
# 	sigmoid <- sigmoid
# 	impulse <- impulse & length(experiment)>2
# 	polynomial <- polynomial & length(experiment)>2

# 	tptslog <- time_transf(tpts, log_shift)

# 	# sigmoid
# 	if( sigmoid ) {
# 		outSM  <- sm.optim.chisq(tpts=tptslog, experiment=experiment
# 			, variance=variance, ninit=nInit, maxit=nIter)
# 		bestSM <- which.min(unlist(outSM[2,]))
# 		pvalSM <- chisq.test.default(experiment=experiment
# 			, model=sigmoidModel(tptslog, outSM[,bestSM]$par)
# 			, variance=variance, df=length(outSM[,bestSM]$par))
# 		dfSM <- length(outSM[,bestSM]$par)
# 	} else dfSM <- NA
# 	# impulse
# 	if( impulse ) {
# 		outIM  <- im.optim.chisq(tpts=tptslog, experiment=experiment, 
# 			variance=variance, ninit=nInit, maxit=nIter)
# 		bestIM <- which.min(unlist(outIM[2,]))
# 		pvalIM <- chisq.test.default(experiment=experiment
# 			, model=impulseModel(tptslog, outIM[,bestIM]$par) 
# 			, variance=variance, df=length(outIM[,bestIM]$par))
# 		dfIM <- length(outIM[,bestIM]$par)
# 	} else dfIM <- NA

# 	# polynomial
# 	if( polynomial ) {
# 		outPN  <- pn.optim.aic(tptslog, experiment, variance )
# 		bestPN <- which.min(unlist(outPN[2,]))
# 		pvalPN <- chisq.test.default(experiment=experiment
# 			, model=.polynomialModel(tptslog, outPN[,bestPN]$par) 
# 			, variance=variance, df=length(outPN[,bestPN]$par))
# 		dfPN <- length(outPN[,bestPN]$par)
# 	} else dfPN <- NA

# 	pvals <- c(
# 		sigmoid=if( sigmoid ) pvalSM else NA
# 		, impulse=if( impulse ) pvalIM else NA
# 		, polynomial=if( polynomial ) pvalPN else NA
# 		)
# 	funcs <- c(sigmoidModelP, impulseModelP, .polynomialModelP)
# 	dfs <- c(dfSM, dfIM, dfPN)
# 	type <- names(pvals)[which.min(pvals)]
# 	df   <- dfs[which.min(pvals)]

# 	if( type=='sigmoid'    ) params <- outSM[,bestSM]$par
# 	if( type=='impulse'    ) params <- outIM[,bestIM]$par
# 	if( type=='polynomial' ) params <- outPN[,bestPN]$par

# 	pval <- pvals[which.min(pvals)]
# 	fun  <- funcs[[which.min(pvals)]]

# 	return(list(type=type, fun=fun , params=params, pval=pval, df=df))

# }


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

# systemSolution <- function(k1F,k2F,k3F,times)
# {

#   system <- function(t,c,parms)
#   {
#     alpha <- parms$alpha
#     beta  <- parms$beta
#     gamma <- parms$gamma

#     r=rep(0,length(c))
#     r[1] <- alpha(t) - beta(t) * c["p"]
#     r[2] <- beta(t) * c["p"] - gamma(t) * c["m"]
#     return(list(r))
#   }

#   cinit <- c(k1F(0)/k2F(0),k1F(0)/k3F(0))
#   names(cinit) <- c("p","m")

#   params <- list(alpha = k1F, beta = k2F, gamma = k3F)

#   modData <- ode(y=cinit, times=times, func=system, parms=params)
#   modData <- c(modData[,"m"],modData[,"p"])

#   names(modData) <- c(rep("mature",length.out = length(times)),rep("premature",length.out = length(times)))

#   return(modData) 
# }

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
  if(length(parameters)==8)
  {
    k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
    k2F <- function(x) return(parameters[7])
    k3F <- function(x) return(parameters[8])
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) return(parameters[5])
    k3F <- function(x) return(parameters[6])
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorVKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==13)
  {
    k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
    k2F <- function(x) return(parameters[7])
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) return(parameters[5])
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[6:9])}
  }
  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorVVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{
  if(length(parameters)==18)
  {
    k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
    k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[13:18])}
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[5:8])}
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[9:12])}
  }
  
  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorKVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{
  if(length(parameters)==13)
  {
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[2:7])}
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}
  }else{
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[2:5])}
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[6:9])}
  }
  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisqFunction(data,modData,datavar)

  return(chi2)
}

errorVVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{
  if(length(parameters)==13)
  {
	 k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
	 k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
	 k3F <- function(x) return(parameters[13])
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[5:8])}
    k3F <- function(x) return(parameters[9])
  }
	if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

	modData <- systemSolution(k1F,k2F,k3F,times)
	
	chi2 <- chisqFunction(data,modData,datavar)

	return(chi2)
}

errorKVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{
	if(length(parameters)==8)
  {
    k1F <- function(x)return(parameters[1])
  	k2F <- function(x){impulseModel(log2(x+a)+c,parameters[2:7])}
  	k3F <- function(x)return(parameters[8])
  }else{
    k1F <- function(x)return(parameters[1])
    k2F <- function(x){sigmoidModel(log2(x+a)+c,parameters[2:5])}
    k3F <- function(x)return(parameters[6])
  }	
	if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)
	
	modData <- systemSolution(k1F,k2F,k3F,times)
	
	chi2 <- chisqFunction(data,modData,datavar)
	
	return(chi2)
}

errorKKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==8)
  {
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) return(parameters[2])
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[3:8])}
  }else{
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) return(parameters[2])
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[3:6])}
  }
  
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

  if(length(parameters)==8)
  {
    k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
    k2F <- function(x) return(parameters[7])
    k3F <- function(x) return(parameters[8])
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) return(parameters[5])
    k3F <- function(x) return(parameters[6])
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)

}

loglikKVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

	if(length(parameters)==8)
  {
    k1F <- function(x)return(parameters[1])
  	k2F <- function(x){impulseModel(log2(x+a)+c,parameters[2:7])}
  	k3F <- function(x)return(parameters[8])
  }else{
    k1F <- function(x)return(parameters[1])
    k2F <- function(x){sigmoidModel(log2(x+a)+c,parameters[2:5])}
    k3F <- function(x)return(parameters[6])
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  logLikelihoodFunction(data, modData, datavar)
}

loglikKKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==8)
  {
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) return(parameters[2])
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[3:8])}
  }else{
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) return(parameters[2])
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[3:6])}
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVVK_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==13)
  {
	 k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
	 k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
	 k3F <- function(x) return(parameters[13])
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[5:8])}
    k3F <- function(x) return(parameters[9])
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVKV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==13)
  {
    k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
    k2F <- function(x) return(parameters[7])
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) return(parameters[5])
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[6:9])}
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikKVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==13)
  {
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[2:7])}
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[8:13])}
  }else{
    k1F <- function(x) return(parameters[1])
    k2F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[2:5])}
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[6:9])}
  }

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVVV_Int_NoNascent <- function(parameters, times, data, datavar, a, c)
{

  if(length(parameters)==18)
  {
    k1F <- function(x) {impulseModel(log2(x+a)+c,parameters[1:6])}
    k2F <- function(x) {impulseModel(log2(x+a)+c,parameters[7:12])}
    k3F <- function(x) {impulseModel(log2(x+a)+c,parameters[13:18])}
  }else{
    k1F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[1:4])}
    k2F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[5:8])}
    k3F <- function(x) {sigmoidModel(log2(x+a)+c,parameters[9:12])}
  }

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

# .inspect.engine_Derivative_NoNascent <- function(tptsOriginal
# 										   , tptsLinear
# 										   , a
# 										   , c
# 										   , concentrations
# 										   , rates
# 										   , BPPARAM=bpparam()
# 										   , na.rm=TRUE
# 										   , verbose=TRUE
# 										   , testOnSmooth=TRUE
# 										   , seed=NULL
# 										   , nInit = nInit
# 										   , nIter = nIter
# 										   , limitModelComplexity = FALSE)
# {

# 	total <- concentrations$total
# 	totalVariance <- concentrations$total_var

# 	premature <- concentrations$preMRNA
# 	prematureVariance <- concentrations$preMRNA_var

# 	mature <- concentrations$mature
# 	matureVariance <- concentrations$mature_var

# 	eiGenes <- rownames(mature)

# 	k2median <- median(rates$gamma,na.rm = TRUE)
# 	k3median <- median(rates$beta,na.rm = TRUE)

# 	# degreesOfFreedom <- length(tptsLinear)-6

# 	matureFitImpulse <- bplapply(eiGenes,function(row)
# 	{
#   		fitSmooth(tpts = tptsLinear
#         	    , tt_c = c
#         	    , experiment = mature[row,]
#         	    , variance = matureVariance[row,]
#         	    , mature = TRUE
#         	    , nInit = nInit
#         	    , nIter = nIter
#         	    , seed = seed)
# 	},BPPARAM=BPPARAM)
# 	names(matureFitImpulse) <- eiGenes

# 	if(all(sapply(matureFitImpulse,is.null)))
# 		stop("No genes have a mature profile possible to fit 
# 			with impulsive functions. Try with the options:
# 			'modelingParams()$estimateRatesWith <- 'int' ' and
# 			'modelingParams()$testOnSmooth <- FALSE'.")

# 	if(any(sapply(matureFitImpulse,is.null))) {
# 		message("Some genes have a mature profile impossible to be fitted with impulsive 
# 			functions therefore they will be excluded from the modelling.")

# 		eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null)]

# 		matureFitImpulse <- matureFitImpulse[eiGenes]
# 		total <- total[eiGenes,]
# 		totalVariance <- totalVariance[eiGenes,]
# 		premature <- premature[eiGenes,]
# 		prematureVariance <- prematureVariance[eiGenes,]
# 		mature <- mature[eiGenes,]
# 		matureVariance <- matureVariance[eiGenes,]
# 	}

# 	# if( degreesOfFreedom > 0 ) {
# 	# 	accept <- pchisq(sapply(matureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2 
# 	# 	if( table(accept)['FALSE']/length(accept) > .5 )
# 	# 		message("More than 50% did not return a good fit of their mature profiles 
# 	# 			with the impulsive smooth function:
# 	# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
# 	# }

# 	if( testOnSmooth ) {

# 		prematureFitImpulse <- bplapply(eiGenes,function(row)
# 		{
# 			fitSmooth(tpts = tptsLinear
# 					, tt_c = c
# 	        		, experiment = premature[row,]
# 	        		, variance = prematureVariance[row,]
# 	        		, mature = FALSE
# 	        		, nInit = nInit
# 	        		, nIter = nIter
# 	        		, seed = seed)     
# 		},BPPARAM=BPPARAM)
# 		names(prematureFitImpulse) <- eiGenes

# 		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
# 			stop("No genes have an expression profile possible to fit 
# 				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")

# 		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
# 			message("Some genes have an expression profile impossible to be fitted 
# 				with impulsive functions therefore they will be excluded from the modelling.")

# 			eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null) & 
# 				!sapply(prematureFitImpulse,is.null)]

# 			prematureFitImpulse <- prematureFitImpulse[eiGenes]
# 			matureFitImpulse <- matureFitImpulse[eiGenes]
# 			total <- total[eiGenes,]
# 			totalVariance <- totalVariance[eiGenes,]
# 			premature <- premature[eiGenes,]
# 			prematureVariance <- prematureVariance[eiGenes,]
# 			mature <- mature[eiGenes,]
# 			matureVariance <- matureVariance[eiGenes,]
# 		}

# 		# if( degreesOfFreedom > 0 ) {
# 		# 	accept <- pchisq(sapply(prematureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2 
# 		# 	if( table(accept)['FALSE']/length(accept) > .5 )
# 		# 		message("More than 50% did not return a good fit of their premature 
# 		# 			profiles with the impulsive smooth function:
# 		# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
# 		# }

# 	}

# 	# Equal to integrative approach

# 	KKK <- bplapply(eiGenes,function(row){
	
# 	matureParameters <- mean(impulseModel(tptsLinear, matureFitImpulse[[row]][1:6]))
# 	k2Parameters <- k2median
# 	k3Parameters <- k3median

# 	unlist(
# 		tryCatch(
# 			optim(c(matureParameters, k2Parameters, k3Parameters)
# 				,errorKKK_Int_NoNascent
# 				,tpts = tptsLinear
# 				,premature = premature[row,]
# 				,mature = mature[row,]
# 				,prematureVariance = prematureVariance[row,]
# 				,matureVariance = matureVariance[row,]
# 				,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							  , par2 = NaN
# 							  , par3 = NaN
# 							  , value = NaN
# 							  , counts.function = NaN
# 						  	  , counts.gradient = NaN
# 						  	  , convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KKK) <- eiGenes
# 	message("Model 0 finished.")

# 	medianAmplVKK <- sapply(eiGenes, function(row)
# 	{

# 		parameters <- unname(c(matureFitImpulse[[row]][1:6]
# 							 , k2median
# 							 , k3median))

# 		k1 <- k1VKK_Der_NoNascent(tptsLinear,parameters,c)
# 		p <- prematureVKK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)

# 		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)

# 		suppressWarnings(optimize( function(x)
# 		{
#     		parameters <- unname(c(matureFitImpulse[[row]][1:6]
# 								 , k2median*x
# 								 , k3median*x))

# 			k1 <- k1VKK_Der_NoNascent(tptsLinear,parameters,c)
#     		p <- prematureVKK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)

# 		if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(c(k1
# 												 ,k2median*x*length(tptsOriginal)
#       									 		 ,k3median*x*length(tptsOriginal)))
#   		},c(1, 1e5) ))$minimum
# 	})

# 	VKK <- bplapply(eiGenes,function(row){

# 		if(medianAmplVKK[[row]] > 10^2)
#   		{
# 			return(c(par1 = NaN, par2 = NaN, par3 = NaN
# 				   , par4 = NaN, par5 = NaN, par6 = NaN
# 				   , par7 = NaN, par8 = NaN, value = NaN
# 				   , counts.function = NaN
# 				   , counts.gradient = NaN, convergence = NaN))
# 		}

# 		matureParameters <- unname(matureFitImpulse[[row]][1:6])

# 		k2Parameters <- k2median * unname(medianAmplVKK[row])
# 		k3Parameters <- k3median * unname(medianAmplVKK[row])
	  
# 		unlist(
# 			tryCatch(
# 	      			optim(unname(c(matureParameters, k2Parameters, k3Parameters))
# 	                  	 ,errorVKK_Der_NoNascent
#         				 ,tpts = tptsLinear
# 			             ,premature = premature[row,]
# 			             ,mature = mature[row,]
# 			             ,prematureVariance = prematureVariance[row,]
# 			             ,matureVariance = matureVariance[row,]
# 	                  	 ,c = c
# 	                  	 ,control = list(maxit = nIter * 100)),
# 	    		error=function(e) c(par1 = NaN
# 	    						  , par2 = NaN
# 	    						  , par3 = NaN
# 	    						  , par4 = NaN
# 	    						  , par5 = NaN
# 	    						  , par6 = NaN
# 	    						  , par7 = NaN
# 	    						  , par8 = NaN
# 	    						  , value = NaN
# 	    						  , counts.function = NaN
# 	    						  , counts.gradient = NaN
# 	    						  , convergence = NaN)
# 	    		)
# 	    )
# 	}, BPPARAM=BPPARAM)
# 	names(VKK) <- eiGenes
# 	message("Model A finished.")

# 	medianAmplVVK <- sapply(eiGenes, function(row)
# 	{
# 		parameters <- c(unname(matureFitImpulse[[row]][1:6])
# 						, rep(k2median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1
# 						, k3median)
# 		k1 <- k1VVK_Der_NoNascent(tptsLinear,parameters, c)
# 		p <- prematureVVK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)
# 		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)
# 		suppressWarnings(optimize( function(x) {
# 			parameters <- unname(c(matureFitImpulse[[row]][1:6]
# 							  , c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 							  , k3median * x))
     
# 			k1 <- k1VVK_Der_NoNascent(tptsLinear,parameters, c)
# 			p <- prematureVVK_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)      
     
# 			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(c(k1,impulseModel(tptsLinear,c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)),k3median*x*length(tptsOriginal)))
#    		}, c(1, 1e5) ))$minimum
# 	})

# 	VVK <- bplapply(eiGenes, function(row)
# 	{

# 		if(medianAmplVVK[[row]] > 10^2)
# 		{
# 			return(c(par1 = NaN, par2 = NaN, par3 = NaN
# 				   , par4 = NaN, par5 = NaN, par6 = NaN
# 				   , par7 = NaN, par8 = NaN, par9 = NaN
# 				   , par10 = NaN, par11 = NaN, par12 = NaN
# 				   , par13 = NaN, value = NaN, counts.function = NaN
# 				   , counts.gradient = NaN, convergence = NaN))
# 		}
  
# 		matureParameters <- unname(matureFitImpulse[[row]][1:6])
# 		k2Parameters <- c(rep(k2median,3) * unname(medianAmplVVK[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k3Parameters <- k3median * unname(medianAmplVVK[row])
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(matureParameters, k2Parameters, k3Parameters))
# 			        	,errorVVK_Der_NoNascent
# 						,tpts = tptsLinear
# 						,premature = premature[row,]
# 						,mature = mature[row,]
# 						,prematureVariance = prematureVariance[row,]
# 						,matureVariance = matureVariance[row,]
# 						,c = c
# 						,control = list(maxit = nIter * 100)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(VVK) <- eiGenes
# 	message("Model AC finished.")

# 	medianAmplVKV <- sapply(eiGenes, function(row)
# 	{

# 		parameters <- c(unname(matureFitImpulse[[row]][1:6])
# 							 , k2median
# 			  				 , rep(k3median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
#   		k1 <- k1VKV_Der_NoNascent(tptsLinear,parameters, c)
# 		p <- prematureVKV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)

# 		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)

# 		suppressWarnings(optimize( function(x) {

# 			parameters <- c(unname(matureFitImpulse[[row]][1:6])
# 								 , k2median * x
# 								 , rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
      
# 			k1 <- k1VKV_Der_NoNascent(tptsLinear,parameters, c)
# 			p <- prematureVKV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)
      
# 			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(k1,k2median*x*length(tptsOriginal),impulseModel(tptsLinear,c(rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)))

# 		}, c(1, 1e5) ))$minimum
# 	})

# 	VKV <- bplapply(eiGenes, function(row)
# 	{

# 		if(medianAmplVKV[[row]] > 10^2)
# 		{
# 			return(c(par1 = NaN, par2 = NaN, par3 = NaN
# 				   , par4 = NaN, par5 = NaN, par6 = NaN
# 				   , par7 = NaN, par8 = NaN, par9 = NaN
# 				   , par10 = NaN, par11 = NaN, par12 = NaN
# 				   , par13 = NaN, value = NaN, counts.function = NaN
# 				   , counts.gradient = NaN, convergence = NaN))
# 		}

# 		matureParameters <- unname(matureFitImpulse[[row]][1:6])
# 		k2Parameters <- k2median * unname(medianAmplVKV[row])
# 		k3Parameters <- c(rep(k3median,3) * unname(medianAmplVKV[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(matureParameters, k2Parameters, k3Parameters))
#           				,errorVKV_Der_NoNascent
# 						,tpts = tptsLinear
# 						,premature = premature[row,]
# 						,mature = mature[row,]
# 						,prematureVariance = prematureVariance[row,]
# 						,matureVariance = matureVariance[row,]
# 						,c = c
# 						,control = list(maxit = nIter * 100)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(VKV) <- eiGenes
# 	message("Model AB finished.")

# 	medianAmplVVV <- sapply(eiGenes, function(row)
# 	{

# 		parameters <- c(unname(matureFitImpulse[[row]][1:6])
# 					  , rep(k2median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1
# 					  , rep(k3median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k1 <- k1VVV_Der_NoNascent(tptsLinear,parameters, c)
# 		p <- prematureVVV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)
	
# 		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)
	
# 		suppressWarnings(optimize( function(x) {
	
# 			parameters <- c(unname(matureFitImpulse[[row]][1:6])
# 						  , rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1
# 						  , rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
	
# 			k1 <- k1VVV_Der_NoNascent(tptsLinear,parameters, c)
# 			p <- prematureVVV_Der_NoNascent(x = tptsLinear, parameters = parameters, c = c)      
	      
# 			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(k1
# 	      										, impulseModel(tptsLinear,c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1))
# 	      										, impulseModel(tptsLinear,c(rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)))
	
# 			}, c(1, 1e5) ))$minimum
# 	})

# 	VVV <- bplapply(eiGenes, function(row){

# 		if(medianAmplVVV[[row]] > 10^2)
# 		{
# 			return(c(par1 = NaN, par2 = NaN, par3 = NaN, par4 = NaN, par5 = NaN
# 				   , par6 = NaN, par7 = NaN, par8 = NaN, par9 = NaN, par10 = NaN
# 				   , par11 = NaN, par12 = NaN, par13 = NaN, par14 = NaN
# 				   , par15 = NaN, par16 = NaN, par17 = NaN, par18 = NaN
# 				   , value = NaN, counts.function = NaN
# 				   , counts.gradient = NaN, convergence = NaN))
# 		}

# 		matureParameters <- unname(matureFitImpulse[[row]][1:6])
# 		k2Parameters <- c(rep(k2median,3) * unname(medianAmplVVV[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k3Parameters <- c(rep(k3median,3) * unname(medianAmplVVV[row]), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)

# 			unlist(
# 			tryCatch(
# 				optim(unname(c(matureParameters, k2Parameters, k3Parameters))
#           				,errorVVV_Der_NoNascent
# 						,tpts = tptsLinear
# 						,premature = premature[row,]
# 						,mature = mature[row,]
# 						,prematureVariance = prematureVariance[row,]
# 						,matureVariance = matureVariance[row,]
# 						,c = c
# 						,control = list(maxit = nIter * 100)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,par14 = NaN
# 								   ,par15 = NaN
# 								   ,par16 = NaN
# 								   ,par17 = NaN
# 								   ,par18 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(VVV) <- eiGenes
# 	message("Model ABC finished.")

# 	KKV <- bplapply(eiGenes, function(row){
 
# 		k1Parameters <- mean(k1VKV_Der_NoNascent(tptsLinear, VKV[[row]], c),na.rm=T)
# 		k2Parameters <- VKV[[row]][7]
# 		k3Parameters <- VKV[[row]][8:13]
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 		        	,errorKKV_Int_NoNascent
# 		        	,times = tptsOriginal
# 		        	,data = c( mature[row,], premature[row,] )
# 		        	,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 		        	,a = a
# 		        	,c = c
# 		        	,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KKV) <- eiGenes
# 	message("Model B finished.")

# 	KVK <- bplapply(eiGenes, function(row){
 
# 		k1Parameters <- mean(k1VVK_Der_NoNascent(tptsLinear, VVK[[row]], c),na.rm=T)
# 		k2Parameters <- VVK[[row]][7:12]
# 		k3Parameters <- VVK[[row]][13]
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 					 ,errorKVK_Int_NoNascent
# 					 ,times = tptsOriginal
# 					 ,data = c( mature[row,], premature[row,] )
# 					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 					 ,a = a
# 					 ,c = c
# 					 ,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KVK) <- eiGenes
# 	message("Model C finished.")

# 	KVV <- bplapply(eiGenes, function(row){
	
# 		k1Parameters <- mean(k1VVV_Der_NoNascent(tptsLinear, VVV[[row]], c),na.rm=T)
# 		k2Parameters <- VVV[[row]][7:12]
# 		k3Parameters <- VVV[[row]][13:18]
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 					 ,errorKVV_Int_NoNascent
# 					 ,times = tptsOriginal
# 					 ,data = c( mature[row,], premature[row,] )
# 					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 					 ,a = a
# 					 ,c = c
# 					 ,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,par9 = NaN
# 							   ,par10 = NaN
# 							   ,par11 = NaN
# 							   ,par12 = NaN
# 							   ,par13 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KVV) <- eiGenes
# 	message("Model BC finished.")

# 	if(testOnSmooth)
# 	{
# 		matureSmooth <- t(sapply(matureFitImpulse,function(i)
# 		{
#  			impulseModel(tptsLinear,i)
# 		}))
# 		prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
# 		{
# 		 	impulseModel(tptsLinear,i)
# 		}))
# 		colnames(matureSmooth) <- tptsOriginal
# 		colnames(prematureSmooth) <- tptsOriginal
# 	}else{
# 		matureSmooth <- mature
# 		prematureSmooth <- premature
# 	}

# 	chi2data <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 											  ,tptsLinear
# 											  ,prematureSmooth[g,]
# 											  ,matureSmooth[g,]
# 											  ,prematureVariance[g,]
# 											  ,matureVariance[g,]),error = function(e)NaN)
	
# 		VKKTemp <- tryCatch(errorVKK_Der_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
# 											  ,tptsLinear
# 											  ,prematureSmooth[g,]
# 											  ,matureSmooth[g,]
# 											  ,prematureVariance[g,]
# 											  ,matureVariance[g,]
# 											  ,c),error = function(e)NaN)
	
# 		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)
	
# 		VVKTemp <- tryCatch(errorVVK_Der_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
# 													 ,tptsLinear
# 													 ,prematureSmooth[g,]
# 													 ,matureSmooth[g,]
# 													 ,prematureVariance[g,]
# 													 ,matureVariance[g,]
# 													 ,c),error = function(e)NaN)
	
# 		VKVTemp <- tryCatch(errorVKV_Der_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
# 													 ,tptsLinear
# 													 ,prematureSmooth[g,]
# 													 ,matureSmooth[g,]
# 													 ,prematureVariance[g,]
# 													 ,matureVariance[g,]
# 													 ,c),error = function(e)NaN)

# 		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)
	
# 		VVVTemp <- tryCatch(errorVVV_Der_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
# 													 ,tptsLinear
# 													 ,prematureSmooth[g,]
# 													 ,matureSmooth[g,]
# 													 ,prematureVariance[g,]
# 													 ,matureVariance[g,]
# 													 ,c),error = function(e)NaN)
	  
# 	  c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
# 	}, BPPARAM=BPPARAM))

# 	if( limitModelComplexity ) {
# 		dof <- c(KKK = 3
# 				,VKK = min(8, length(tptsOriginal)+2)
# 				,KVK = min(8, length(tptsOriginal)+2)
# 				,KKV = min(8, length(tptsOriginal)+2)
# 				,VVK = min(13, 2*length(tptsOriginal)+1)
# 				,VKV = min(13, 2*length(tptsOriginal)+1)
# 				,KVV = min(13, 2*length(tptsOriginal)+1)
# 				,VVV = min(18, 2*length(tptsOriginal))
# 				)
# 	} else {
# 		dof <- c(KKK = 3
# 				,VKK = 8
# 				,KVK = 8
# 				,KKV = 8
# 				,VVK = 13
# 				,VKV = 13
# 				,KVV = 13
# 				,VVV = 18)
# 	}

# 	# P values
# 	pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], max(c(0,2*length(tptsOriginal)-dof['KKK'])))
# 						,VKK=pchisq(chi2data[,'VKK'], max(c(0,2*length(tptsOriginal)-dof['VKK'])))
# 						,KVK=pchisq(chi2data[,'KVK'], max(c(0,2*length(tptsOriginal)-dof['KVK'])))
# 						,KKV=pchisq(chi2data[,'KKV'], max(c(0,2*length(tptsOriginal)-dof['KKV'])))
# 						,VVK=pchisq(chi2data[,'VVK'], max(c(0,2*length(tptsOriginal)-dof['VVK'])))
# 						,VKV=pchisq(chi2data[,'VKV'], max(c(0,2*length(tptsOriginal)-dof['VKV'])))
# 						,KVV=pchisq(chi2data[,'KVV'], max(c(0,2*length(tptsOriginal)-dof['KVV'])))
# 						,VVV=pchisq(chi2data[,'VVV'], max(c(0,2*length(tptsOriginal)-dof['VVV']))))

	

# 	logLikelihood <- t(mcsapply(eiGenes,function(g)
# 	{
# 		prematureKKKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureKKK_Int_NoNascent(x = tptsLinear[t]
# 																						   , parameters = KKK[[g]][grep("par",names(KKK[[g]]))])))
# 		prematureVKKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVKK_Der_NoNascent(x = tptsLinear[t]
# 		                                                                         , parameters = VKK[[g]][grep("par",names(VKK[[g]]))]
# 		                                                                         , c = c)))
# 		prematureVVKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVVK_Der_NoNascent(x = tptsLinear[t]
# 		                                                            , parameters = VVK[[g]][grep("par",names(VVK[[g]]))]
# 		                                                            , c = c)))
# 		prematureVKVTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVKV_Der_NoNascent(x = tptsLinear[t]
# 		                                                            , parameters = VKV[[g]][grep("par",names(VKV[[g]]))]
# 		                                                            , c = c)))
# 		prematureVVVTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVVV_Der_NoNascent(x = tptsLinear[t]
# 		                                                            , parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
# 		                                                            , c = c)))		

# 		matureKKKTemp <- c(sapply(seq_along(tptsLinear),function(t)KKK[[g]][grep("par",names(KKK[[g]]))][1]))

# 		matureVKKTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
#         		                                                               , par = VKK[[g]][grep("par",names(VKK[[g]]))][1:6])))
# 		matureVVKTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
#                                                                       		   , par = VVK[[g]][grep("par",names(VVK[[g]]))][1:6])))
# 		matureVKVTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
#                                                                       		   , par = VKV[[g]][grep("par",names(VKV[[g]]))][1:6])))
# 		matureVVVTemp <- c(sapply(seq_along(tptsLinear),function(t)impulseModel(x = tptsLinear[t]
#                                                                       		   , par = VVV[[g]][grep("par",names(VVV[[g]]))][1:6])))

# 		modelKKK <- c(matureKKKTemp,prematureKKKTemp)
# 		modelVKK <- c(matureVKKTemp,prematureVKKTemp)
# 		modelVVK <- c(matureVVKTemp,prematureVVKTemp)
# 		modelVKV <- c(matureVKVTemp,prematureVKVTemp)
# 		modelVVV <- c(matureVVVTemp,prematureVVVTemp)

# 		KKKTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
# 		                               , model = modelKKK
# 		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
# 		VKKTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
# 		                               , model = modelVKK
# 		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
# 		VVKTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
# 		                               , model = modelVVK
# 		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
# 		VKVTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
# 		                               , model = modelVKV
# 		                               , variance = c(matureVariance[g,],prematureVariance[g,]))
# 		VVVTemp <- logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
# 		                               , model = modelVVV
# 		                               , variance = c(matureVariance[g,],prematureVariance[g,]))

# 		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
# 												   ,tptsOriginal
# 												   ,c(matureSmooth[g,],prematureSmooth[g,])
# 												   ,c(matureVariance[g,],prematureVariance[g,])
# 												   ,a
# 												   ,c),error = function(e)NaN)

# 		KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
# 												   ,tptsOriginal
# 												   ,c(matureSmooth[g,],prematureSmooth[g,])
# 												   ,c(matureVariance[g,],prematureVariance[g,])
# 												   ,a
# 												   ,c),error = function(e)NaN)

# 		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
# 												   ,tptsOriginal
# 												   ,c(matureSmooth[g,],prematureSmooth[g,])
# 												   ,c(matureVariance[g,],prematureVariance[g,])
# 												   ,a
# 												   ,c),error = function(e)NaN)

# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

# 	},BPPARAM=BPPARAM))

# 	AIC <- t(apply(logLikelihood,1,function(row)
# 	{
# 		2*(dof - row) 
# 	}))

# 	AICc <- t(apply(logLikelihood,1,function(row)
# 	{
# 	 	2*(dof - row) + (2*dof*(dof+1))/max(0,2*length(tptsOriginal)-dof-1)
# 	}))

# 	rownames(pvaluesdata) <- rownames(logLikelihood) <- rownames(AIC) <- rownames(AICc) <- eiGenes

# 	ratesSpecs <- lapply(eiGenes,function(gene)
#  		{
#  			list(
#  				"0" = list(alpha = list(fun = constantModelP
# 									 ,type = "constant"
# 									 ,df = 1
# 									 ,params = c(alpha = unname(KKK[[gene]]["par1"])))
# 						,beta = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(beta = unname(KKK[[gene]]["par3"])))
# 	 				   	,gamma = list(fun = constantModelP
# 									 ,type = "constant"
# 									 ,df = 1
# 									 ,params = c(gamma = unname(KKK[[gene]]["par2"])))
# 						,test = log(pvaluesdata[gene,"KKK"])
# 						,logLik = logLikelihood[gene,"KKK"]
# 						,AIC = AIC[gene,"KKK"]
# 						,AICc = AICc[gene,"KKK"]
# 						,counts = c("function"=unname(KKK[[gene]]["counts.function"]), gradient=unname(KKK[[gene]]["counts.gradient"]))
# 						,convergence = unname(KKK[[gene]]["convergence"])
# 						,message = NULL)

#  				,a = list(alpha = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(alpha = unname(VKK[[gene]][1:6])))
#  					  		,beta = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(beta = unname(VKK[[gene]][8])))
# 	 				  		,gamma = list(fun = constantModelP
# 										,type = "constant"
# 										,df = 1
# 										,params = c(gamma = unname(VKK[[gene]][7])))
# 						,test = log(pvaluesdata[gene,"VKK"])
# 						,logLik = logLikelihood[gene,"VKK"]
# 						,AIC = AIC[gene,"VKK"]
# 						,AICc = AICc[gene,"VKK"]
# 						,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
# 						,convergence = unname(VKK[[gene]]["convergence"])
# 						,message = NULL)
#  				,b = list(alpha = list(fun = constantModelP
# 										,type = "constant"
# 										,df = 1
# 										,params = c(alpha = unname(KKV[[gene]][1])))
#  					  		,beta = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(beta = unname(KKV[[gene]][3:8])))
# 	 				  		,gamma = list(fun = constantModelP
# 										,type = "constant"
# 										,df = 1
# 										,params = c(gamma = unname(KKV[[gene]][2])))
# 						,test = log(pvaluesdata[gene,"KKV"])
# 						,logLik = logLikelihood[gene,"KKV"]
# 						,AIC = AIC[gene,"KKV"]
# 						,AICc = AICc[gene,"KKV"]
# 						,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
# 						,convergence = unname(KKV[[gene]]["convergence"])
# 						,message = NULL)
#  				,c = list(alpha = list(fun = constantModelP
# 										,type = "constant"
# 										,df = 1
# 										,params = c(alpha = unname(KVK[[gene]][1])))
#  					  		,beta = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(beta = unname(KVK[[gene]][8])))
# 	 				  		,gamma = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(gamma = unname(KVK[[gene]][2:7])))
# 						,test = log(pvaluesdata[gene,"KVK"])
# 						,logLik = logLikelihood[gene,"KVK"]
# 						,AIC = AIC[gene,"KVK"]
# 						,AICc = AICc[gene,"KVK"]
# 						,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
# 						,convergence = unname(KVK[[gene]]["convergence"])
# 						,message = NULL)
#  				,ab = list(alpha = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(alpha = unname(VKV[[gene]][1:6])))
#  					  		,beta = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(beta = unname(VKV[[gene]][8:13])))
# 	 				  		,gamma = list(fun = constantModelP
# 										,type = "constant"
# 										,df = 1
# 										,params = c(gamma = unname(VKV[[gene]][7])))
# 						,test = log(pvaluesdata[gene,"VKV"])
# 						,logLik = logLikelihood[gene,"VKV"]
# 						,AIC = AIC[gene,"VKV"]
# 						,AICc = AICc[gene,"VKV"]
# 						,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
# 						,convergence = unname(VKV[[gene]]["convergence"])
# 						,message = NULL)
#  				,ac = list(alpha = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(alpha = unname(VVK[[gene]][1:6])))
#  					  		,beta = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(beta = unname(VVK[[gene]][13])))
# 	 				  		,gamma = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(gamma = unname(VVK[[gene]][7:12])))
# 						,test = log(pvaluesdata[gene,"VVK"])
# 						,logLik = logLikelihood[gene,"VVK"]
# 						,AIC = AIC[gene,"VVK"]
# 						,AICc = AICc[gene,"VVK"]
# 						,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
# 						,convergence = unname(VVK[[gene]]["convergence"])
# 						,message = NULL)
#  				,bc = list(alpha = list(fun = constantModelP
# 										,type = "constant"
# 										,df = 1
# 										,params = c(alpha = unname(KVV[[gene]][1])))
#  					  		,beta = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(beta = unname(KVV[[gene]][8:13])))
# 	 				  		,gamma = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(gamma = unname(KVV[[gene]][2:7])))
# 						,test = log(pvaluesdata[gene,"KVV"])
# 						,logLik = logLikelihood[gene,"KVV"]
# 						,AIC = AIC[gene,"KVV"]
# 						,AICc = AICc[gene,"KVV"]
# 						,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
# 						,convergence = unname(KVV[[gene]]["convergence"])
# 						,message = NULL)
#  				,abc = list(alpha = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(alpha = unname(VVV[[gene]][1:6])))
#  					  		,beta = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(beta = unname(VVV[[gene]][13:18])))
# 	 				  		,gamma = list(fun = impulseModelP
# 										,type = "impulse"
# 										,df = 6
# 										,params = c(gamma = unname(VVV[[gene]][7:12])))
# 						,test = log(pvaluesdata[gene,"VVV"])
# 						,logLik = logLikelihood[gene,"VVV"]
# 						,AIC = AIC[gene,"VVV"]
# 						,AICc = AICc[gene,"VVV"]
# 						,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
# 						,convergence = unname(VVV[[gene]]["convergence"])
# 						,message = NULL)
# 			)
#  		})

# 		return(ratesSpecs)
# }

# .makeModel_Derivative <- function(tpts, hyp, log_shift, time_transf, ode, .rxnrate, c= NaN, geneBestModel = NULL)
# {

# 	params <- list()
# 	params$alpha <- function(x) 
# 		hyp$alpha$fun$value(time_transf(x, log_shift, c), hyp$alpha$par)
# 	params$beta  <- function(x) 
# 		hyp$beta$fun$value(time_transf(x, log_shift, c), hyp$beta$par)
# 	params$gamma <- function(x) 
# 		hyp$gamma$fun$value(time_transf(x, log_shift, c), hyp$gamma$par)

# 	matureTemp <- params$alpha(tpts)

# 	if(geneBestModel == "0")
# 	{
# 		prematureTemp <- sapply(tpts,function(t)prematureKKK_Int_NoNascent(t,c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params)))
# 		k1Temp <- sapply(tpts,function(t)k1KKK_Int_NoNascent(t,c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params)))
		
# 		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
# 		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))
	
# 	}else if(geneBestModel == "a")
# 	{
# 		prematureTemp <- prematureVKK_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
# 		k1Temp <- k1VKK_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

# 		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
# 		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))

# 	}else if(geneBestModel == "ac")
# 	{
# 		prematureTemp <- prematureVVK_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
# 		k1Temp <- k1VVK_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

# 		k2Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$gamma$params)
# 		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))

# 	}else if(geneBestModel == "ab")
# 	{
# 		prematureTemp <- prematureVKV_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
# 		k1Temp <- k1VKV_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

# 		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
# 		k3Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$beta$params)

# 	}else if(geneBestModel == "abc")
# 	{
# 		prematureTemp <- prematureVVV_Der_NoNascent(time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
# 		k1Temp <- k1VVV_Der_NoNascent(time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

# 		k2Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$gamma$params)
# 		k3Temp <- impulseModel(time_transf(tpts,log_shift,c), hyp$beta$params)

# 	}

# 	totalTemp <- matureTemp + prematureTemp

# 	data.frame(time = tpts, preMRNA = prematureTemp, total = totalTemp, alpha = k1Temp, beta = k3Temp, gamma = k2Temp)

# }

# .inspect.engine_Integrative_NoNascent <- function(tptsOriginal
# 										    , tptsLinear
# 										    , a
# 										    , c
# 										    , concentrations
# 										    , rates
# 										    , BPPARAM=bpparam()
# 										    , na.rm=TRUE
# 										    , verbose=TRUE
# 										    , testOnSmooth=TRUE
# 										    , seed = NULL
# 										    , nInit = nInit
# 										    , nIter = nIter
# 										    , limitModelComplexity = FALSE
# 										    , sigmoid = FALSE)
# {
#  	total <- concentrations$total
#  	totalVariance <- concentrations$total_var

#  	premature <- concentrations$preMRNA
#  	prematureVariance <- concentrations$preMRNA_var

#  	mature <- concentrations$mature
#  	matureVariance <- concentrations$mature_var

#  	eiGenes <- rownames(mature)

#  	k2median <- median(rates$gamma,na.rm = TRUE)
#  	k3median <- median(rates$beta,na.rm = TRUE)

#  	if( testOnSmooth ) {
#  		###### fit smooth functions on prematue and mature data
#  		###### eventually, exclude the genes that cannot be
#  		###### fit in either premature or mature
#  		matureFitImpulse <- bplapply(eiGenes,function(row)
#  		{
#  	  		fitSmooth(tpts = tptsLinear
#  	        	    , tt_c = c
#  	        	    , experiment = mature[row,]
#  	        	    , variance = matureVariance[row,]
#  	        	    , mature = TRUE
#  	        	    , nInit = nInit
#  	        	    , nIter = nIter
#  	        	    , seed = seed)
#  		},BPPARAM=BPPARAM)
#  		names(matureFitImpulse) <- eiGenes
#  		prematureFitImpulse <- bplapply(eiGenes,function(row)
#  		{
#  			fitSmooth(tpts = tptsLinear
#  					, tt_c = c
#  	        		, experiment = premature[row,]
#  	        		, variance = prematureVariance[row,]
#  	        		, mature = FALSE
#  	        		, nInit = nInit
#  	        		, nIter = nIter
#  	        		, seed = seed)     
#  		},BPPARAM=BPPARAM)
#  		names(prematureFitImpulse) <- eiGenes
#  		## in case the number of degrees of freedom allows the estimation
#  		## quantify how many genes have an acceptable chi2 test value (<0.2)
#  		## in both mature and premature
		
#  		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
#  			stop("No genes have an expression profile possible to fit 
#  				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")
#  		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
#  			message("Some genes have an expression profile impossible to be fitted 
#  				with impulsive functions therefore they will be excluded from the modelling.")
#  			eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null) & 
#  				!sapply(prematureFitImpulse,is.null)]
#  			prematureFitImpulse <- prematureFitImpulse[eiGenes]
#  			matureFitImpulse <- matureFitImpulse[eiGenes]
		
#  			total <- total[eiGenes,]
#  			totalVariance <- totalVariance[eiGenes,]
#  			premature <- premature[eiGenes,]
#  			prematureVariance <- prematureVariance[eiGenes,]
#  			mature <- mature[eiGenes,]
#  			matureVariance <- matureVariance[eiGenes,]
 		
#  		### Selection of the genes which are suitable for the analysis, impulse response!
#  		}
#  	}

#  	KKK <- bplapply(eiGenes,function(row){

#  			matureParameters <- mean(mature[row,])
#  			k2Parameters <- k2median
#  			k3Parameters <- k3median
#  			unlist(
#  				tryCatch(
#  					optim(c(matureParameters, k2Parameters, k3Parameters)
#  						,errorKKK_Int_NoNascent
#  						,tpts = tptsLinear
#  						,premature = premature[row,]
#  						,mature = mature[row,]
#  						,prematureVariance = prematureVariance[row,]
#  						,matureVariance = matureVariance[row,]
#  						,control = list(maxit = nIter)),
#  					error=function(e) c(par1 = NaN
#  									  , par2 = NaN
#  									  , par3 = NaN
#  									  , value = NaN
#  									  , counts.function = NaN
#  								  	  , counts.gradient = NaN
#  								  	  , convergence = e)
#  				)
#  			)
#  		}, BPPARAM=BPPARAM)
#  		names(KKK) <- eiGenes
#  		message("Model 0 finished.")
	
# 	VKK <- bplapply(eiGenes,function(row){

# 			k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 			k2Parameters <- KKK[[row]][2]
# 			k3Parameters <- KKK[[row]][3]
	  
# 	  		#this cycle was introduced to test slopes of different sign
# 	  		#impulsesParameters <- sapply(c(1,-1),function(slope)
# 	  		impulsesParameters <- sapply(1,function(slope)
# 			{
# 				unlist(
# 					tryCatch(
# 	      				optim(unname(c(c(k1Parameters,slope), k2Parameters, k3Parameters))
# 	        	          	 ,errorVKK_Int_NoNascent
# 	        	          	 ,times = tptsOriginal
# 	        	          	 ,data = c( mature[row,], premature[row,] )
# 	        	          	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 	        	          	 ,a = a
# 	        	          	 ,c = c
# 	        	          	 ,control = list(maxit = nIter)),
# 	    			error=function(e) c(par1 = NaN
# 	    							  , par2 = NaN
# 	    							  , par3 = NaN
# 	    							  , par4 = NaN
# 	    							  , par5 = NaN
# 	    							  , par6 = NaN
# 	    							  , par7 = NaN
# 	    							  , par8 = NaN
# 	    							  , value = NaN
# 	    							  , counts.function = NaN
# 	    							  , counts.gradient = NaN
# 	    							  , convergence = e)
# 	    			)
# 	    		)
# 			})
	
# 			impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]
	
# 			if(!sigmoid){return(impulsesParameters)}
	
# 			fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters=impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
	
# 			sigmoidsParameters <- unlist(
# 				tryCatch(
# 	      			optim(unname(c(fitParameters_k1, impulsesParameters[7], impulsesParameters[8]))
# 	                  	 ,errorVKK_Int_NoNascent
# 	                  	 ,times = tptsOriginal
# 	                  	 ,data = c( mature[row,], premature[row,] )
# 	                  	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 	                  	 ,a = a
# 	                  	 ,c = c
# 	                  	 ,control = list(maxit = nIter)),
# 	    		error=function(e) c(par1 = NaN
# 	    						  , par2 = NaN
# 	    						  , par3 = NaN
# 	    						  , par4 = NaN
# 	    						  , par5 = NaN
# 	    						  , par6 = NaN
# 	    						  , par7 = NaN
# 	    						  , par8 = NaN
# 	    						  , value = NaN
# 	    						  , counts.function = NaN
# 	    						  , counts.gradient = NaN
# 	    						  , convergence = e)
# 	    		)
# 	    	)

# 			return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(VKK) <- eiGenes

# 	VKK_sigmoid <- lapply(VKK,"[[","sigmoid")
# 	VKK_impulse <- lapply(VKK,"[[","impulse")

# 	message("Model A finished.")

# 	KKV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
		
# 		#impulsesParameters <- sapply(c(1,-1),function(slope)
# 		impulsesParameters <- sapply(1,function(slope)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(k1Parameters, k2Parameters, c(k3Parameters,slope)))
# 			        	,errorKKV_Int_NoNascent
# 			        	,times = tptsOriginal
# 			        	,data = c( mature[row,], premature[row,] )
# 			        	,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 			        	,a = a
# 			        	,c = c
# 			        	,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)		
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 		if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[3:8],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(impulsesParameters[1], impulsesParameters[2], fitParameters_k3))
#                   	 ,errorKKV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(KKV) <- eiGenes

# 	KKV_sigmoid <- lapply(KKV,"[[","sigmoid")
# 	KKV_impulse <- lapply(KKV,"[[","impulse")

# 	message("Model B finished.")

# 	KVK <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- KKK[[row]][3]
		
#   		#impulsesParameters <- sapply(c(1,-1),function(slope)
#   		impulsesParameters <- sapply(1,function(slope)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(k1Parameters, c(k2Parameters,slope), k3Parameters))
# 						 ,errorKVK_Int_NoNascent
# 						 ,times = tptsOriginal
# 						 ,data = c( mature[row,], premature[row,] )
# 						 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						 ,a = a
# 						 ,c = c
# 						 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[2:7],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(impulsesParameters[1], fitParameters_k2, impulsesParameters[8]))
#                   	 ,errorKVK_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(KVK) <- eiGenes

# 	KVK_sigmoid <- lapply(KVK,"[[","sigmoid")
# 	KVK_impulse <- lapply(KVK,"[[","impulse")

# 	message("Model C finished.")

# 	VKV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
		
#   		#impulsesParameters <- apply(matrix(c(1,1,-1,-1,1,-1,1,-1),nrow=2),2,function(slopes)
#   		impulsesParameters <- apply(matrix(c(1,1),nrow=2),2,function(slopes)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(c(k1Parameters,slopes[[1]]), k2Parameters, c(k3Parameters,slopes[[2]])))
# 						 ,errorVKV_Int_NoNascent
# 						 ,times = tptsOriginal
# 						 ,data = c( mature[row,], premature[row,] )
# 						 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						 ,a = a
# 						 ,c = c
# 						 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 			                       ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]
	
# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[8:13],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(fitParameters_k1, impulsesParameters[7], fitParameters_k3))
#                   	 ,errorVKV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)
# 		return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(VKV) <- eiGenes

# 	VKV_sigmoid <- lapply(VKV,"[[","sigmoid")
# 	VKV_impulse <- lapply(VKV,"[[","impulse")

# 	message("Model AB finished.")

# 	VVK <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- KKK[[row]][3]

#   		#impulsesParameters <- apply(matrix(c(1,1,-1,-1,1,-1,1,-1),nrow=2),2,function(slopes)
#   		impulsesParameters <- apply(matrix(c(1,1),nrow=2),2,function(slopes)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(c(k1Parameters,slopes[[1]]), c(k2Parameters,slopes[[2]]), k3Parameters))
#     					 ,errorVVK_Int_NoNascent
#     					 ,times = tptsOriginal
#     					 ,data = c( mature[row,], premature[row,] )
#     					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#     					 ,a = a
#     					 ,c = c
#     					 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								  , par2 = NaN
# 								  , par3 = NaN
# 								  , par4 = NaN
# 								  , par5 = NaN
# 								  , par6 = NaN
# 								  , par7 = NaN
# 								  , par8 = NaN
#         	        			  , par9 = NaN
#         	        			  , par10 = NaN
#         	        			  , par11 = NaN
#         	        			  , par12 = NaN
#         	        			  , par13 = NaN
#         	        			  , value = NaN
#         	        			  , counts.function = NaN
#         	        			  , counts.gradient = NaN
#         	        			  , convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[7:12],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(fitParameters_k1, fitParameters_k2, impulsesParameters[13]))
#                   	 ,errorVVK_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(VVK) <- eiGenes

# 	VVK_sigmoid <- lapply(VVK,"[[","sigmoid")
# 	VVK_impulse <- lapply(VVK,"[[","impulse")

# 	message("Model AC finished.")

# 	KVV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)

#   		#impulsesParameters <- apply(matrix(c(1,1,-1,-1,1,-1,1,-1),nrow=2),2,function(slopes)
#   		impulsesParameters <- apply(matrix(c(1,1),nrow=2),2,function(slopes)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(k1Parameters, c(k2Parameters,slopes[[1]]), c(k3Parameters,slopes[[2]])))/3
# 						 ,errorKVV_Int_NoNascent
# 						 ,times = tptsOriginal
# 						 ,data = c( mature[row,], premature[row,] )
# 						 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						 ,a = a
# 						 ,c = c
# 						 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[2:7],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[8:13],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(impulsesParameters[1], fitParameters_k2, fitParameters_k3))
#                   	 ,errorKVV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(KVV) <- eiGenes

# 	KVV_sigmoid <- lapply(KVV,"[[","sigmoid")
# 	KVV_impulse <- lapply(KVV,"[[","impulse")

# 	message("Model BC finished.")

# 	VVV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
		
#   		#impulsesParameters <- apply(matrix(c(1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,-1,1,-1),nrow=3),2,function(slopes)
# 		impulsesParameters <- apply(matrix(c(1,1,1),nrow=3),2,function(slopes)		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(c(k1Parameters,slopes[[1]]), c(k2Parameters,slopes[[2]]), c(k3Parameters,slopes[[3]])))
# 						,errorVVV_Int_NoNascent
# 						,times = tptsOriginal
# 						,data = c( mature[row,], premature[row,] )
# 						,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						,a = a
# 						,c = c
# 						,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,par14 = NaN
# 								   ,par15 = NaN
# 								   ,par16 = NaN
# 								   ,par17 = NaN
# 								   ,par18 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[7:12],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[13:18],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(fitParameters_k1, fitParameters_k2, fitParameters_k3))
#                   	 ,errorVVV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		return(list("sigmoid"=sigmoidsParameters,"impulse"=impulsesParameters))

# 	}, BPPARAM=BPPARAM)
# 	names(VVV) <- eiGenes

# 	VVV_sigmoid <- lapply(VVV,"[[","sigmoid")
# 	VVV_impulse <- lapply(VVV,"[[","impulse")

# 	message("Model ABC finished.")

# 	if(testOnSmooth)
# 	{
# 			matureSmooth <- t(sapply(matureFitImpulse,function(i)
# 			{
# 				impulseModel(tptsLinear,i)
# 			}))
		
# 			prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
# 			{
# 		  		impulseModel(tptsLinear,i)
# 			}))

# 			colnames(matureSmooth) <- tptsOriginal
# 			colnames(prematureSmooth) <- tptsOriginal

# 	}else{
# 		matureSmooth <- mature
# 		prematureSmooth <- premature
# 	}

# 	## Chisquare - sigmoid
# 	chi2data_sigmoid <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 											  ,tptsOriginal
# 											  ,prematureSmooth[g,]
# 											  ,matureSmooth[g,]
# 											  ,prematureVariance[g,]
# 											  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(errorVKK_Int_NoNascent(VKK_sigmoid[[g]][grep("par",names(VKK_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK_sigmoid[[g]][grep("par",names(KVK_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV_sigmoid[[g]][grep("par",names(KKV_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVKTemp <- tryCatch(errorVVK_Int_NoNascent(VVK_sigmoid[[g]][grep("par",names(VVK_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VKVTemp <- tryCatch(errorVKV_Int_NoNascent(VKV_sigmoid[[g]][grep("par",names(VKV_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV_sigmoid[[g]][grep("par",names(KVV_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVVTemp <- tryCatch(errorVVV_Int_NoNascent(VVV_sigmoid[[g]][grep("par",names(VVV_sigmoid[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)
	  
# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
# 	}, BPPARAM = BPPARAM))

# 	## Chisquare - impulse
# 	chi2data_impulse <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 											  ,tptsOriginal
# 											  ,prematureSmooth[g,]
# 											  ,matureSmooth[g,]
# 											  ,prematureVariance[g,]
# 											  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(errorVKK_Int_NoNascent(VKK_impulse[[g]][grep("par",names(VKK_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK_impulse[[g]][grep("par",names(KVK_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV_impulse[[g]][grep("par",names(KKV_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVKTemp <- tryCatch(errorVVK_Int_NoNascent(VVK_impulse[[g]][grep("par",names(VVK_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VKVTemp <- tryCatch(errorVKV_Int_NoNascent(VKV_impulse[[g]][grep("par",names(VKV_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV_impulse[[g]][grep("par",names(KVV_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVVTemp <- tryCatch(errorVVV_Int_NoNascent(VVV_impulse[[g]][grep("par",names(VVV_impulse[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)
	  
# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
# 	}, BPPARAM = BPPARAM))

# 	if( limitModelComplexity )
# 	{
# 		complexity_sigmoid <- c('KKK'=min(3,length(tptsOriginal))
# 					,'VKK'=min(6,length(tptsOriginal)+2)
# 					,'KVK'=min(6,length(tptsOriginal)+2)
# 					,'KKV'=min(6,length(tptsOriginal)+2)
# 					,'VVK'=min(9,length(tptsOriginal)+1)
# 					,'VKV'=min(9,length(tptsOriginal)+1)
# 					,'KVV'=min(9,length(tptsOriginal)+1)
# 					,'VVV'=min(12,length(tptsOriginal)))
# 		complexity_sigmoid[complexity_sigmoid<0] <- 0

# 		complexity_impulse <- c('KKK'=min(3,length(tptsOriginal))
# 						,'VKK'=min(8,length(tptsOriginal)+2)
# 						,'KVK'=min(8,length(tptsOriginal)+2)
# 						,'KKV'=min(8,length(tptsOriginal)+2)
# 						,'VVK'=min(13,length(tptsOriginal)+1)
# 						,'VKV'=min(13,length(tptsOriginal)+1)
# 						,'KVV'=min(13,length(tptsOriginal)+1)
# 						,'VVV'=min(18,length(tptsOriginal)))
# 		complexity_impulse[complexity_impulse<0] <- 0
# 	}else
# 	{
# 		complexity_sigmoid <- c('KKK'=3
# 					,'VKK'=6
# 					,'KVK'=6
# 					,'KKV'=6
# 					,'VVK'=9
# 					,'VKV'=9
# 					,'KVV'=9
# 					,'VVV'=12)
# 		complexity_sigmoid[complexity_sigmoid<0] <- 0

# 		complexity_impulse <- c('KKK'=3
# 					,'VKK'=8
# 					,'KVK'=8
# 					,'KKV'=8
# 					,'VVK'=13
# 					,'VKV'=13
# 					,'KVV'=13
# 					,'VVV'=18)
# 		complexity_impulse[complexity_impulse<0] <- 0

# 	}

# 	pvaluesdata_sigmoid <- cbind(KKK=pchisq(chi2data_sigmoid[,'KKK'], 2*length(tptsOriginal)-complexity_sigmoid['KKK'])
# 						,VKK=pchisq(chi2data_sigmoid[,'VKK'], 2*length(tptsOriginal)-complexity_sigmoid['VKK'])
# 						,KVK=pchisq(chi2data_sigmoid[,'KVK'], 2*length(tptsOriginal)-complexity_sigmoid['KVK'])
# 						,KKV=pchisq(chi2data_sigmoid[,'KKV'], 2*length(tptsOriginal)-complexity_sigmoid['KKV'])
# 						,VVK=pchisq(chi2data_sigmoid[,'VVK'], 2*length(tptsOriginal)-complexity_sigmoid['VVK'])
# 						,VKV=pchisq(chi2data_sigmoid[,'VKV'], 2*length(tptsOriginal)-complexity_sigmoid['VKV'])
# 						,KVV=pchisq(chi2data_sigmoid[,'KVV'], 2*length(tptsOriginal)-complexity_sigmoid['KVV'])
# 						,VVV=pchisq(chi2data_sigmoid[,'VVV'], 2*length(tptsOriginal)-complexity_sigmoid['VVV']))

# 	pvaluesdata_impulse <- cbind(KKK=pchisq(chi2data_impulse[,'KKK'], 2*length(tptsOriginal)-complexity_impulse['KKK'])
# 						,VKK=pchisq(chi2data_impulse[,'VKK'], 2*length(tptsOriginal)-complexity_impulse['VKK'])
# 						,KVK=pchisq(chi2data_impulse[,'KVK'], 2*length(tptsOriginal)-complexity_impulse['KVK'])
# 						,KKV=pchisq(chi2data_impulse[,'KKV'], 2*length(tptsOriginal)-complexity_impulse['KKV'])
# 						,VVK=pchisq(chi2data_impulse[,'VVK'], 2*length(tptsOriginal)-complexity_impulse['VVK'])
# 						,VKV=pchisq(chi2data_impulse[,'VKV'], 2*length(tptsOriginal)-complexity_impulse['VKV'])
# 						,KVV=pchisq(chi2data_impulse[,'KVV'], 2*length(tptsOriginal)-complexity_impulse['KVV'])
# 						,VVV=pchisq(chi2data_impulse[,'VVV'], 2*length(tptsOriginal)-complexity_impulse['VVV']))

# 	logLikelihood_sigmoid <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(loglikKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 										  ,tptsOriginal,prematureSmooth[g,]
# 										  ,matureSmooth[g,]
# 										  ,prematureVariance[g,]
# 										  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(loglikVKK_Int_NoNascent(VKK_sigmoid[[g]][grep("par",names(VKK_sigmoid[[g]]))]
# 										  ,tptsOriginal
# 										  ,c(matureSmooth[g,],prematureSmooth[g,])
# 										  ,c(matureVariance[g,],prematureVariance[g,])
# 										  ,a
# 										  ,c),error = function(e)NaN)

# 	  	KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK_sigmoid[[g]][grep("par",names(KVK_sigmoid[[g]]))]
#         				                  ,tptsOriginal
#         				                  ,c(matureSmooth[g,],prematureSmooth[g,])
#         				                  ,c(matureVariance[g,],prematureVariance[g,])
#         				                  ,a
#         				                  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV_sigmoid[[g]][grep("par",names(KKV_sigmoid[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVKTemp <- tryCatch(loglikVVK_Int_NoNascent(VVK_sigmoid[[g]][grep("par",names(VVK_sigmoid[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VKVTemp <- tryCatch(loglikVKV_Int_NoNascent(VKV_sigmoid[[g]][grep("par",names(VKV_sigmoid[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV_sigmoid[[g]][grep("par",names(KVV_sigmoid[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVVTemp <- tryCatch(loglikVVV_Int_NoNascent(VVV_sigmoid[[g]][grep("par",names(VVV_sigmoid[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)

# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

# 	}, BPPARAM=BPPARAM))

# 	logLikelihood_impulse <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(loglikKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 										  ,tptsOriginal,prematureSmooth[g,]
# 										  ,matureSmooth[g,]
# 										  ,prematureVariance[g,]
# 										  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(loglikVKK_Int_NoNascent(VKK_impulse[[g]][grep("par",names(VKK_impulse[[g]]))]
# 										  ,tptsOriginal
# 										  ,c(matureSmooth[g,],prematureSmooth[g,])
# 										  ,c(matureVariance[g,],prematureVariance[g,])
# 										  ,a
# 										  ,c),error = function(e)NaN)

# 	  	KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK_impulse[[g]][grep("par",names(KVK_impulse[[g]]))]
#         				                  ,tptsOriginal
#         				                  ,c(matureSmooth[g,],prematureSmooth[g,])
#         				                  ,c(matureVariance[g,],prematureVariance[g,])
#         				                  ,a
#         				                  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV_impulse[[g]][grep("par",names(KKV_impulse[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVKTemp <- tryCatch(loglikVVK_Int_NoNascent(VVK_impulse[[g]][grep("par",names(VVK_impulse[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VKVTemp <- tryCatch(loglikVKV_Int_NoNascent(VKV_impulse[[g]][grep("par",names(VKV_impulse[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV_impulse[[g]][grep("par",names(KVV_impulse[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVVTemp <- tryCatch(loglikVVV_Int_NoNascent(VVV_impulse[[g]][grep("par",names(VVV_impulse[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)

# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

# 	}, BPPARAM=BPPARAM))

# 	AIC_sigmoid <- t(sapply(1:nrow(logLikelihood_sigmoid),function(r)
# 	{
# 		2*(complexity_sigmoid - logLikelihood_sigmoid[r,])
# 	}))

# 	AIC_impulse <- t(sapply(1:nrow(logLikelihood_impulse),function(r)
# 	{
# 		2*(complexity_impulse - logLikelihood_impulse[r,])
# 	}))

# 	rownames(chi2data_sigmoid) <- rownames(chi2data_impulse) <- 
# 	rownames(pvaluesdata_sigmoid) <- rownames(pvaluesdata_impulse) <- 
# 	rownames(logLikelihood_sigmoid) <- rownames(logLikelihood_impulse) <- 
# 	rownames(AIC_sigmoid) <- rownames(AIC_impulse) <- 
# 	eiGenes

# 	## Gene labeling according to the best model, sigmoidal or impulsive

# 	chi2data <- matrix(NaN,nrow=nrow(chi2data_sigmoid),ncol=ncol(chi2data_sigmoid))
# 	pvaluesdata <- matrix(NaN,nrow=nrow(pvaluesdata_sigmoid),ncol=ncol(pvaluesdata_sigmoid))
# 	logLikelihood <- matrix(NaN,nrow=nrow(logLikelihood_sigmoid),ncol=ncol(logLikelihood_sigmoid))
# 	AIC <- matrix(NaN,nrow=nrow(AIC_sigmoid),ncol=ncol(AIC_sigmoid))

# 	rownames(chi2data) <- rownames(chi2data_sigmoid)
# 	rownames(pvaluesdata) <- rownames(pvaluesdata_sigmoid)
# 	rownames(logLikelihood) <- rownames(logLikelihood_sigmoid)
# 	rownames(AIC) <- rownames(AIC_sigmoid)

# 	colnames(chi2data) <- colnames(chi2data_sigmoid)
# 	colnames(pvaluesdata) <- colnames(pvaluesdata_sigmoid)
# 	colnames(logLikelihood) <- colnames(logLikelihood_sigmoid)
# 	colnames(AIC) <- colnames(AIC_sigmoid)

# 	for(g in eiGenes)
# 	{
# 		if(min(AIC_sigmoid[g,],na.rm=TRUE)<=min(AIC_impulse[g,],na.rm=TRUE))
# 		{
# 			VKK[[g]] <- VKK_sigmoid[[g]]
# 			KVK[[g]] <- KVK_sigmoid[[g]]
# 			KKV[[g]] <- KKV_sigmoid[[g]]
# 			VKV[[g]] <- VKV_sigmoid[[g]]
# 			VVK[[g]] <- VVK_sigmoid[[g]]
# 			KVV[[g]] <- KVV_sigmoid[[g]]
# 			VVV[[g]] <- VVV_sigmoid[[g]]

# 			chi2data[g,] <- chi2data_sigmoid[g,]
# 			pvaluesdata[g,] <- pvaluesdata_sigmoid[g,]
# 			logLikelihood[g,] <- logLikelihood_sigmoid[g,]
# 			AIC[g,] <- AIC_sigmoid[g,]			
# 		}else{
# 			VKK[[g]] <- VKK_impulse[[g]]
# 			KVK[[g]] <- KVK_impulse[[g]]
# 			KKV[[g]] <- KKV_impulse[[g]]
# 			VKV[[g]] <- VKV_impulse[[g]]
# 			VVK[[g]] <- VVK_impulse[[g]]
# 			KVV[[g]] <- KVV_impulse[[g]]
# 			VVV[[g]] <- VVV_impulse[[g]]

# 			chi2data[g,] <- chi2data_impulse[g,]
# 			pvaluesdata[g,] <- pvaluesdata_impulse[g,]
# 			logLikelihood[g,] <- logLikelihood_impulse[g,]
# 			AIC[g,] <- AIC_impulse[g,]			

# 		}
# 	}
# 		ratesSpecs <- lapply(eiGenes,function(gene)
# 		{

# 			paramsTmp <- c("KKK"=length(grep("par",names(KKK[[gene]])))
# 						  ,"VKK"=length(grep("par",names(VKK[[gene]])))
# 						  ,"KVK"=length(grep("par",names(KVK[[gene]])))
# 						  ,"KKV"=length(grep("par",names(KKV[[gene]])))
# 						  ,"VVK"=length(grep("par",names(VVK[[gene]])))
# 						  ,"VKV"=length(grep("par",names(VKV[[gene]])))
# 						  ,"KVV"=length(grep("par",names(KVV[[gene]])))
# 						  ,"VVV"=length(grep("par",names(VVV[[gene]]))))

# 				zero_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = k1KKK_Int_NoNascent(x = 0,par = KKK[[gene]])))
# 					    ,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KKK[[gene]]["par3"])))
#  				   		,gamma = list(fun = constantModelP
# 									 ,type = "constant"
# 									 ,df = 1
# 									 ,params = c(gamma = unname(KKK[[gene]]["par2"])))
# 					,test = log(pvaluesdata[gene,"KKK"])
# 					,logLik = logLikelihood[gene,"KKK"]
# 					,AIC = AIC[gene,"KKK"]
# 					,counts = c("function"=unname(KKK[[gene]]["counts.function"]), gradient=unname(KKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKK[[gene]]["convergence"])
# 					,message = NULL)

# 			if(paramsTmp[["VKK"]]==8)
# 			{
# 				a_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VKK[[gene]][1:6])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VKK[[gene]][8])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKK[[gene]][7])))
# 					,test = log(pvaluesdata[gene,"VKK"])
# 					,logLik = logLikelihood[gene,"VKK"]
# 					,AIC = AIC[gene,"VKK"]
# 					,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKK[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				a_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VKK[[gene]][1:4])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VKK[[gene]][6])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKK[[gene]][5])))
# 					,test = log(pvaluesdata[gene,"VKK"])
# 					,logLik = logLikelihood[gene,"VKK"]
# 					,AIC = AIC[gene,"VKK"]
# 					,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKK[[gene]]["convergence"])
# 					,message = NULL)				
# 			}

# 			if(paramsTmp[["KKV"]]==8)
# 			{
# 				b_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KKV[[gene]][1])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(KKV[[gene]][3:8])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(KKV[[gene]][2])))
# 					,test = log(pvaluesdata[gene,"KKV"])
# 					,logLik = logLikelihood[gene,"KKV"]
# 					,AIC = AIC[gene,"KKV"]
# 					,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				b_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KKV[[gene]][1])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(KKV[[gene]][3:6])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(KKV[[gene]][2])))
# 					,test = log(pvaluesdata[gene,"KKV"])
# 					,logLik = logLikelihood[gene,"KKV"]
# 					,AIC = AIC[gene,"KKV"]
# 					,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			if(paramsTmp[["KVK"]]==8)
# 			{
# 				c_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVK[[gene]][1])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KVK[[gene]][8])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(KVK[[gene]][2:7])))
# 					,test = log(pvaluesdata[gene,"KVK"])
# 					,logLik = logLikelihood[gene,"KVK"]
# 					,AIC = AIC[gene,"KVK"]
# 					,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVK[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				c_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVK[[gene]][1])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KVK[[gene]][6])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(KVK[[gene]][2:5])))
# 					,test = log(pvaluesdata[gene,"KVK"])
# 					,logLik = logLikelihood[gene,"KVK"]
# 					,AIC = AIC[gene,"KVK"]
# 					,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVK[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			if(paramsTmp[["VKV"]]==13)
# 			{
# 				ab_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VKV[[gene]][1:6])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(VKV[[gene]][8:13])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKV[[gene]][7])))
# 					,test = log(pvaluesdata[gene,"VKV"])
# 					,logLik = logLikelihood[gene,"VKV"]
# 					,AIC = AIC[gene,"VKV"]
# 					,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				ab_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VKV[[gene]][1:4])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(VKV[[gene]][6:9])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKV[[gene]][5])))
# 					,test = log(pvaluesdata[gene,"VKV"])
# 					,logLik = logLikelihood[gene,"VKV"]
# 					,AIC = AIC[gene,"VKV"]
# 					,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			if(paramsTmp[["VVK"]]==13)
# 			{
# 				ac_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VVK[[gene]][1:6])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VVK[[gene]][13])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(VVK[[gene]][7:12])))
# 					,test = log(pvaluesdata[gene,"VVK"])
# 					,logLik = logLikelihood[gene,"VVK"]
# 					,AIC = AIC[gene,"VVK"]
# 					,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVK[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				ac_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VVK[[gene]][1:4])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VVK[[gene]][9])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(VVK[[gene]][5:8])))
# 					,test = log(pvaluesdata[gene,"VVK"])
# 					,logLik = logLikelihood[gene,"VVK"]
# 					,AIC = AIC[gene,"VVK"]
# 					,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVK[[gene]]["convergence"])
# 					,message = NULL)				
# 			}

# 			if(paramsTmp[["KVV"]]==13)
# 			{
# 				bc_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVV[[gene]][1])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(KVV[[gene]][8:13])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(KVV[[gene]][2:7])))
# 					,test = log(pvaluesdata[gene,"KVV"])
# 					,logLik = logLikelihood[gene,"KVV"]
# 					,AIC = AIC[gene,"KVV"]
# 					,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				bc_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVV[[gene]][1])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(KVV[[gene]][6:9])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(KVV[[gene]][2:5])))
# 					,test = log(pvaluesdata[gene,"KVV"])
# 					,logLik = logLikelihood[gene,"KVV"]
# 					,AIC = AIC[gene,"KVV"]
# 					,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVV[[gene]]["convergence"])
# 					,message = NULL)				
# 			}

# 			if(paramsTmp[["VVV"]]==18)
# 			{
# 				abc_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VVV[[gene]][1:6])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(VVV[[gene]][13:18])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(VVV[[gene]][7:12])))
# 					,test = log(pvaluesdata[gene,"VVV"])
# 					,logLik = logLikelihood[gene,"VVV"]
# 					,AIC = AIC[gene,"VVV"]
# 					,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				abc_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VVV[[gene]][1:4])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(VVV[[gene]][9:12])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(VVV[[gene]][5:8])))
# 					,test = log(pvaluesdata[gene,"VVV"])
# 					,logLik = logLikelihood[gene,"VVV"]
# 					,AIC = AIC[gene,"VVV"]
# 					,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVV[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			list("0"=zero_out,"a"=a_out,"b"=b_out,"c"=c_out,"ab"=ab_out,"ac"=ac_out,"bc"=bc_out,"abc"=abc_out)
# 		}
# 		)
# 	return(ratesSpecs)
# }

# .inspect.engine_Integrative_NoNascent_MixedImpulseAndSigmoid <- function(tptsOriginal
# 										    , tptsLinear
# 										    , a
# 										    , c
# 										    , concentrations
# 										    , rates
# 										    , BPPARAM=bpparam()
# 										    , na.rm=TRUE
# 										    , verbose=TRUE
# 										    , testOnSmooth=TRUE
# 										    , seed = NULL
# 										    , nInit = nInit
# 										    , nIter = nIter
# 										    , limitModelComplexity = FALSE
# 										    , sigmoid = FALSE)
# {
#  	total <- concentrations$total
#  	totalVariance <- concentrations$total_var

#  	premature <- concentrations$preMRNA
#  	prematureVariance <- concentrations$preMRNA_var

#  	mature <- concentrations$mature
#  	matureVariance <- concentrations$mature_var

#  	eiGenes <- rownames(mature)

#  	k2median <- median(rates$gamma,na.rm = TRUE)
#  	k3median <- median(rates$beta,na.rm = TRUE)

#  	if( testOnSmooth ) {
#  		###### fit smooth functions on prematue and mature data
#  		###### eventually, exclude the genes that cannot be
#  		###### fit in either premature or mature
#  		matureFitImpulse <- bplapply(eiGenes,function(row)
#  		{
#  	  		fitSmooth(tpts = tptsLinear
#  	        	    , tt_c = c
#  	        	    , experiment = mature[row,]
#  	        	    , variance = matureVariance[row,]
#  	        	    , mature = TRUE
#  	        	    , nInit = nInit
#  	        	    , nIter = nIter
#  	        	    , seed = seed)
#  		},BPPARAM=BPPARAM)
#  		names(matureFitImpulse) <- eiGenes
#  		prematureFitImpulse <- bplapply(eiGenes,function(row)
#  		{
#  			fitSmooth(tpts = tptsLinear
#  					, tt_c = c
#  	        		, experiment = premature[row,]
#  	        		, variance = prematureVariance[row,]
#  	        		, mature = FALSE
#  	        		, nInit = nInit
#  	        		, nIter = nIter
#  	        		, seed = seed)     
#  		},BPPARAM=BPPARAM)
#  		names(prematureFitImpulse) <- eiGenes
#  		## in case the number of degrees of freedom allows the estimation
#  		## quantify how many genes have an acceptable chi2 test value (<0.2)
#  		## in both mature and premature
		
#  		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
#  			stop("No genes have an expression profile possible to fit 
#  				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")
#  		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
#  			message("Some genes have an expression profile impossible to be fitted 
#  				with impulsive functions therefore they will be excluded from the modelling.")
#  			eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null) & 
#  				!sapply(prematureFitImpulse,is.null)]
#  			prematureFitImpulse <- prematureFitImpulse[eiGenes]
#  			matureFitImpulse <- matureFitImpulse[eiGenes]
		
#  			total <- total[eiGenes,]
#  			totalVariance <- totalVariance[eiGenes,]
#  			premature <- premature[eiGenes,]
#  			prematureVariance <- prematureVariance[eiGenes,]
#  			mature <- mature[eiGenes,]
#  			matureVariance <- matureVariance[eiGenes,]
#  		}
#  	}

#  	KKK <- bplapply(eiGenes,function(row){

#  			matureParameters <- mean(mature[row,])
#  			k2Parameters <- k2median
#  			k3Parameters <- k3median
#  			unlist(
#  				tryCatch(
#  					optim(c(matureParameters, k2Parameters, k3Parameters)
#  						,errorKKK_Int_NoNascent
#  						,tpts = tptsLinear
#  						,premature = premature[row,]
#  						,mature = mature[row,]
#  						,prematureVariance = prematureVariance[row,]
#  						,matureVariance = matureVariance[row,]
#  						,control = list(maxit = nIter)),
#  					error=function(e) c(par1 = NaN
#  									  , par2 = NaN
#  									  , par3 = NaN
#  									  , value = NaN
#  									  , counts.function = NaN
#  								  	  , counts.gradient = NaN
#  								  	  , convergence = e)
#  				)
#  			)
#  		}, BPPARAM=BPPARAM)
#  		names(KKK) <- eiGenes
#  		message("Model 0 finished.")
	
# 	VKK <- bplapply(eiGenes,function(row){

# 			k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 			k2Parameters <- KKK[[row]][2]
# 			k3Parameters <- KKK[[row]][3]
	  
# 	  		#this cycle was introduced to test slopes of different sign
# 	  		#impulsesParameters <- sapply(c(1,-1),function(slope)
# 	  		impulsesParameters <- sapply(1,function(slope)
# 			{
# 				unlist(
# 					tryCatch(
# 	      				optim(unname(c(c(k1Parameters,slope), k2Parameters, k3Parameters))
# 	        	          	 ,errorVKK_Int_NoNascent
# 	        	          	 ,times = tptsOriginal
# 	        	          	 ,data = c( mature[row,], premature[row,] )
# 	        	          	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 	        	          	 ,a = a
# 	        	          	 ,c = c
# 	        	          	 ,control = list(maxit = nIter)),
# 	    			error=function(e) c(par1 = NaN
# 	    							  , par2 = NaN
# 	    							  , par3 = NaN
# 	    							  , par4 = NaN
# 	    							  , par5 = NaN
# 	    							  , par6 = NaN
# 	    							  , par7 = NaN
# 	    							  , par8 = NaN
# 	    							  , value = NaN
# 	    							  , counts.function = NaN
# 	    							  , counts.gradient = NaN
# 	    							  , convergence = e)
# 	    			)
# 	    		)
# 			})
	
# 			impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]
	
# 		if(!sigmoid){return(impulsesParameters)}
	
# 			fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters=impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
	
# 			sigmoidsParameters <- unlist(
# 				tryCatch(
# 	      			optim(unname(c(fitParameters_k1, impulsesParameters[7], impulsesParameters[8]))
# 	                  	 ,errorVKK_Int_NoNascent
# 	                  	 ,times = tptsOriginal
# 	                  	 ,data = c( mature[row,], premature[row,] )
# 	                  	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 	                  	 ,a = a
# 	                  	 ,c = c
# 	                  	 ,control = list(maxit = nIter)),
# 	    		error=function(e) c(par1 = NaN
# 	    						  , par2 = NaN
# 	    						  , par3 = NaN
# 	    						  , par4 = NaN
# 	    						  , par5 = NaN
# 	    						  , par6 = NaN
# 	    						  , par7 = NaN
# 	    						  , par8 = NaN
# 	    						  , value = NaN
# 	    						  , counts.function = NaN
# 	    						  , counts.gradient = NaN
# 	    						  , convergence = e)
# 	    		)
# 	    	)
	
# 			models <- list(sigmoidsParameters, impulsesParameters)
# 			## Selection of the best model
	
# 			bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 							pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))
	
# 			model <- models[[bestModel]]
# 			model <- model[-which(names(model)=="value")]
# 			model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])
	
# 			return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(VKK) <- eiGenes
# 	message("Model A finished.")

# 	KKV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
		
# 		#impulsesParameters <- sapply(c(1,-1),function(slope)
# 		impulsesParameters <- sapply(1,function(slope)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(k1Parameters, k2Parameters, c(k3Parameters,slope)))
# 			        	,errorKKV_Int_NoNascent
# 			        	,times = tptsOriginal
# 			        	,data = c( mature[row,], premature[row,] )
# 			        	,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 			        	,a = a
# 			        	,c = c
# 			        	,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)		
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 		if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[3:8],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(impulsesParameters[1], impulsesParameters[2], fitParameters_k3))
#                   	 ,errorKKV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		models <- list(sigmoidsParameters, impulsesParameters)
# 		## Selection of the best model
# 		bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 						pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))

# 		model <- models[[bestModel]]
# 		model <- model[-which(names(model)=="value")]
# 		model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])

# 		return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(KKV) <- eiGenes
# 	message("Model B finished.")

# 	KVK <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- KKK[[row]][3]
		
#   		#impulsesParameters <- sapply(c(1,-1),function(slope)
#   		impulsesParameters <- sapply(1,function(slope)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(k1Parameters, c(k2Parameters,slope), k3Parameters))
# 						 ,errorKVK_Int_NoNascent
# 						 ,times = tptsOriginal
# 						 ,data = c( mature[row,], premature[row,] )
# 						 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						 ,a = a
# 						 ,c = c
# 						 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[2:7],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(impulsesParameters[1], fitParameters_k2, impulsesParameters[8]))
#                   	 ,errorKVK_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		models <- list(sigmoidsParameters, impulsesParameters)
# 		## Selection of the best model
# 		bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 						pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))

# 		model <- models[[bestModel]]
# 		model <- model[-which(names(model)=="value")]
# 		model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])

# 		return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(KVK) <- eiGenes
# 	message("Model C finished.")

# 	VKV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
		
#   		#impulsesParameters <- apply(matrix(c(1,1,-1,-1,1,-1,1,-1),nrow=2),2,function(slopes)
#   		impulsesParameters <- apply(matrix(c(1,1),nrow=2),2,function(slopes)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(c(k1Parameters,slopes[[1]]), k2Parameters, c(k3Parameters,slopes[[2]])))
# 						 ,errorVKV_Int_NoNascent
# 						 ,times = tptsOriginal
# 						 ,data = c( mature[row,], premature[row,] )
# 						 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						 ,a = a
# 						 ,c = c
# 						 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 			                       ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]
	
# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[8:13],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(fitParameters_k1, impulsesParameters[7], fitParameters_k3))
#                   	 ,errorVKV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)
# 		models <- list(sigmoidsParameters, impulsesParameters)
# 		## Selection of the best model

# 		bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 						pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))

# 		model <- models[[bestModel]]
# 		model <- model[-which(names(model)=="value")]
# 		model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])

# 		return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(VKV) <- eiGenes
# 	message("Model AB finished.")

# 	VVK <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- KKK[[row]][3]

#   		#impulsesParameters <- apply(matrix(c(1,1,-1,-1,1,-1,1,-1),nrow=2),2,function(slopes)
#   		impulsesParameters <- apply(matrix(c(1,1),nrow=2),2,function(slopes)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(c(k1Parameters,slopes[[1]]), c(k2Parameters,slopes[[2]]), k3Parameters))
#     					 ,errorVVK_Int_NoNascent
#     					 ,times = tptsOriginal
#     					 ,data = c( mature[row,], premature[row,] )
#     					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#     					 ,a = a
#     					 ,c = c
#     					 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								  , par2 = NaN
# 								  , par3 = NaN
# 								  , par4 = NaN
# 								  , par5 = NaN
# 								  , par6 = NaN
# 								  , par7 = NaN
# 								  , par8 = NaN
#         	        			  , par9 = NaN
#         	        			  , par10 = NaN
#         	        			  , par11 = NaN
#         	        			  , par12 = NaN
#         	        			  , par13 = NaN
#         	        			  , value = NaN
#         	        			  , counts.function = NaN
#         	        			  , counts.gradient = NaN
#         	        			  , convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[7:12],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(fitParameters_k1, fitParameters_k2, impulsesParameters[13]))
#                   	 ,errorVVK_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		models <- list(sigmoidsParameters, impulsesParameters)
# 		## Selection of the best model
# 		bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 								 pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))

# 		model <- models[[bestModel]]
# 		model <- model[-which(names(model)=="value")]
# 		model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])

# 		return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(VVK) <- eiGenes
# 	message("Model AC finished.")

# 	KVV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)

#   		#impulsesParameters <- apply(matrix(c(1,1,-1,-1,1,-1,1,-1),nrow=2),2,function(slopes)
#   		impulsesParameters <- apply(matrix(c(1,1),nrow=2),2,function(slopes)
# 		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(k1Parameters, c(k2Parameters,slopes[[1]]), c(k3Parameters,slopes[[2]])))/3
# 						 ,errorKVV_Int_NoNascent
# 						 ,times = tptsOriginal
# 						 ,data = c( mature[row,], premature[row,] )
# 						 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						 ,a = a
# 						 ,c = c
# 						 ,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[2:7],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[8:13],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(impulsesParameters[1], fitParameters_k2, fitParameters_k3))
#                   	 ,errorKVV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		models <- list(sigmoidsParameters, impulsesParameters)
# 		## Selection of the best model
# 		bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 								 pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))

# 		model <- models[[bestModel]]
# 		model <- model[-which(names(model)=="value")]
# 		model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])

# 		return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(KVV) <- eiGenes
# 	message("Model BC finished.")

# 	VVV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2)
		
#   		#impulsesParameters <- apply(matrix(c(1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,-1,1,-1),nrow=3),2,function(slopes)
# 		impulsesParameters <- apply(matrix(c(1,1,1),nrow=3),2,function(slopes)		{
# 			unlist(
# 				tryCatch(
# 					optim(unname(c(c(k1Parameters,slopes[[1]]), c(k2Parameters,slopes[[2]]), c(k3Parameters,slopes[[3]])))
# 						,errorVVV_Int_NoNascent
# 						,times = tptsOriginal
# 						,data = c( mature[row,], premature[row,] )
# 						,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 						,a = a
# 						,c = c
# 						,control = list(maxit = nIter)),
# 				error=function(e) c(par1 = NaN
# 								   ,par2 = NaN
# 								   ,par3 = NaN
# 								   ,par4 = NaN
# 								   ,par5 = NaN
# 								   ,par6 = NaN
# 								   ,par7 = NaN
# 								   ,par8 = NaN
# 								   ,par9 = NaN
# 								   ,par10 = NaN
# 								   ,par11 = NaN
# 								   ,par12 = NaN
# 								   ,par13 = NaN
# 								   ,par14 = NaN
# 								   ,par15 = NaN
# 								   ,par16 = NaN
# 								   ,par17 = NaN
# 								   ,par18 = NaN
# 								   ,value = NaN
# 								   ,counts.function = NaN
# 								   ,counts.gradient = NaN
# 								   ,convergence = NaN)
# 				)
# 			)
# 		})

# 		impulsesParameters <- impulsesParameters[,which.min(impulsesParameters["value",])]

# 	if(!sigmoid){return(impulsesParameters)}

# 		fitParameters_k1 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[1:6],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k2 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[7:12],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})
# 		fitParameters_k3 <- tryCatch(fromImpulseToSigmoid(impulsesParameters[13:18],tpts=tptsOriginal,a=a,c=c,nIter=nIter),error=function(e){rep(NaN,4)})

# 		sigmoidsParameters <- unlist(
# 			tryCatch(
#       			optim(unname(c(fitParameters_k1, fitParameters_k2, fitParameters_k3))
#                   	 ,errorVVV_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)

# 		models <- list(sigmoidsParameters, impulsesParameters)
# 		## Selection of the best model
# 		bestModel <- which.min(c(pchisq(sigmoidsParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(sigmoidsParameters))))),
# 								 pchisq(impulsesParameters["value"],max(0,2*length(tptsOriginal)-length(grep("par",names(impulsesParameters)))))))

# 		model <- models[[bestModel]]
# 		model <- model[-which(names(model)=="value")]
# 		model <- c(model,"sigmoid"=sigmoidsParameters["value"],"impulse"=impulsesParameters["value"])

# 		return(model)

# 	}, BPPARAM=BPPARAM)
# 	names(VVV) <- eiGenes
# 	message("Model ABC finished.")

# 	if(testOnSmooth)
# 	{
# 			matureSmooth <- t(sapply(matureFitImpulse,function(i)
# 			{
# 				impulseModel(tptsLinear,i)
# 			}))
		
# 			prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
# 			{
# 		  		impulseModel(tptsLinear,i)
# 			}))

# 			colnames(matureSmooth) <- tptsOriginal
# 			colnames(prematureSmooth) <- tptsOriginal

# 	}else{
# 		matureSmooth <- mature
# 		prematureSmooth <- premature
# 	}

# 	## Chisquare
# 	chi2data <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 											  ,tptsOriginal
# 											  ,prematureSmooth[g,]
# 											  ,matureSmooth[g,]
# 											  ,prematureVariance[g,]
# 											  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(errorVKK_Int_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVKTemp <- tryCatch(errorVVK_Int_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VKVTemp <- tryCatch(errorVKV_Int_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVVTemp <- tryCatch(errorVVV_Int_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)
	  
# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
# 	}, BPPARAM = BPPARAM))

# 	if( limitModelComplexity )
# 	{
# 		dof <- cbind('KKK'=2*length(tptsOriginal)-sapply(KKK,function(i){min(length(grep("par",names(i))), length(tptsOriginal))})
# 					,'VKK'=2*length(tptsOriginal)-sapply(VKK,function(i){min(length(grep("par",names(i))), length(tptsOriginal)+2)})
# 					,'KVK'=2*length(tptsOriginal)-sapply(KVK,function(i){min(length(grep("par",names(i))), length(tptsOriginal)+2)})
# 					,'KKV'=2*length(tptsOriginal)-sapply(KKV,function(i){min(length(grep("par",names(i))), length(tptsOriginal)+2)})
# 					,'VVK'=2*length(tptsOriginal)-sapply(VVK,function(i){min(length(grep("par",names(i))), length(tptsOriginal)+1)})
# 					,'VKV'=2*length(tptsOriginal)-sapply(VKV,function(i){min(length(grep("par",names(i))), length(tptsOriginal)+1)})
# 					,'KVV'=2*length(tptsOriginal)-sapply(KVV,function(i){min(length(grep("par",names(i))), length(tptsOriginal)+1)})
# 					,'VVV'=2*length(tptsOriginal)-sapply(VVV,function(i){min(length(grep("par",names(i))), length(tptsOriginal))}))
# 		dof[dof<0] <- 0
# 	}else
# 	{
# 		dof <- cbind('KKK'=2*length(tptsOriginal)-sapply(KKK,function(i)length(grep("par",names(i))))
# 					,'VKK'=2*length(tptsOriginal)-sapply(VKK,function(i)length(grep("par",names(i))))
# 					,'KVK'=2*length(tptsOriginal)-sapply(KVK,function(i)length(grep("par",names(i))))
# 					,'KKV'=2*length(tptsOriginal)-sapply(KKV,function(i)length(grep("par",names(i))))
# 					,'VVK'=2*length(tptsOriginal)-sapply(VVK,function(i)length(grep("par",names(i))))
# 					,'VKV'=2*length(tptsOriginal)-sapply(VKV,function(i)length(grep("par",names(i))))
# 					,'KVV'=2*length(tptsOriginal)-sapply(KVV,function(i)length(grep("par",names(i))))
# 					,'VVV'=2*length(tptsOriginal)-sapply(VVV,function(i)length(grep("par",names(i)))))
# 		dof[dof<0] <- 0
# 	}

# 	pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], dof[,'KKK'])
# 						,VKK=pchisq(chi2data[,'VKK'], dof[,'VKK'])
# 						,KVK=pchisq(chi2data[,'KVK'], dof[,'KVK'])
# 						,KKV=pchisq(chi2data[,'KKV'], dof[,'KKV'])
# 						,VVK=pchisq(chi2data[,'VVK'], dof[,'VVK'])
# 						,VKV=pchisq(chi2data[,'VKV'], dof[,'VKV'])
# 						,KVV=pchisq(chi2data[,'KVV'], dof[,'KVV'])
# 						,VVV=pchisq(chi2data[,'VVV'], dof[,'VVV']))

# 	logLikelihood <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(loglikKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 										  ,tptsOriginal,prematureSmooth[g,]
# 										  ,matureSmooth[g,]
# 										  ,prematureVariance[g,]
# 										  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(loglikVKK_Int_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
# 										  ,tptsOriginal
# 										  ,c(matureSmooth[g,],prematureSmooth[g,])
# 										  ,c(matureVariance[g,],prematureVariance[g,])
# 										  ,a
# 										  ,c),error = function(e)NaN)

# 	  	KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
#         				                  ,tptsOriginal
#         				                  ,c(matureSmooth[g,],prematureSmooth[g,])
#         				                  ,c(matureVariance[g,],prematureVariance[g,])
#         				                  ,a
#         				                  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVKTemp <- tryCatch(loglikVVK_Int_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VKVTemp <- tryCatch(loglikVKV_Int_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVVTemp <- tryCatch(loglikVVV_Int_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)

# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

# 	}, BPPARAM=BPPARAM))

# 		AIC <- t(sapply(1:nrow(logLikelihood),function(r)
# 		{
# 			2*(dof[r,] - logLikelihood[r,])
# 		}))

# 		AICc <- t(sapply(1:nrow(logLikelihood),function(r)
# 		{
# 			2*(dof[r,] - logLikelihood[r,]) + (2*dof[r,]*(dof[r,]+1))/max(0,2*length(tptsOriginal)-dof[r,]-1)
# 		}))
		
# 		rownames(pvaluesdata) <- rownames(logLikelihood) <- rownames(AIC) <- rownames(AICc) <- eiGenes

# 		ratesSpecs <- lapply(eiGenes,function(gene)
# 		{

# 			paramsTmp <- c("KKK"=length(grep("par",names(KKK[[gene]])))
# 						  ,"VKK"=length(grep("par",names(VKK[[gene]])))
# 						  ,"KVK"=length(grep("par",names(KVK[[gene]])))
# 						  ,"KKV"=length(grep("par",names(KKV[[gene]])))
# 						  ,"VVK"=length(grep("par",names(VVK[[gene]])))
# 						  ,"VKV"=length(grep("par",names(VKV[[gene]])))
# 						  ,"KVV"=length(grep("par",names(KVV[[gene]])))
# 						  ,"VVV"=length(grep("par",names(VVV[[gene]]))))

# 				zero_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = k1KKK_Int_NoNascent(x = 0,par = KKK[[gene]])))
# 					    ,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KKK[[gene]]["par3"])))
#  				   		,gamma = list(fun = constantModelP
# 									 ,type = "constant"
# 									 ,df = 1
# 									 ,params = c(gamma = unname(KKK[[gene]]["par2"])))
# 					,test = log(pvaluesdata[gene,"KKK"])
# 					,logLik = logLikelihood[gene,"KKK"]
# 					,AIC = AIC[gene,"KKK"]
# 					,AICc = AICc[gene,"KKK"]
# 					,counts = c("function"=unname(KKK[[gene]]["counts.function"]), gradient=unname(KKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKK[[gene]]["convergence"])
# 					,message = NULL)

# 			if(paramsTmp[["VKK"]]==8)
# 			{
# 				a_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VKK[[gene]][1:6])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VKK[[gene]][8])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKK[[gene]][7])))
# 					,test = log(pvaluesdata[gene,"VKK"])
# 					,logLik = logLikelihood[gene,"VKK"]
# 					,AIC = AIC[gene,"VKK"]
# 					,AICc = AICc[gene,"VKK"]
# 					,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKK[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				a_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VKK[[gene]][1:4])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VKK[[gene]][6])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKK[[gene]][5])))
# 					,test = log(pvaluesdata[gene,"VKK"])
# 					,logLik = logLikelihood[gene,"VKK"]
# 					,AIC = AIC[gene,"VKK"]
# 					,AICc = AICc[gene,"VKK"]
# 					,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKK[[gene]]["convergence"])
# 					,message = NULL)				
# 			}

# 			if(paramsTmp[["KKV"]]==8)
# 			{
# 				b_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KKV[[gene]][1])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(KKV[[gene]][3:8])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(KKV[[gene]][2])))
# 					,test = log(pvaluesdata[gene,"KKV"])
# 					,logLik = logLikelihood[gene,"KKV"]
# 					,AIC = AIC[gene,"KKV"]
# 					,AICc = AICc[gene,"KKV"]
# 					,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				b_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KKV[[gene]][1])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(KKV[[gene]][3:6])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(KKV[[gene]][2])))
# 					,test = log(pvaluesdata[gene,"KKV"])
# 					,logLik = logLikelihood[gene,"KKV"]
# 					,AIC = AIC[gene,"KKV"]
# 					,AICc = AICc[gene,"KKV"]
# 					,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			if(paramsTmp[["KVK"]]==8)
# 			{
# 				c_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVK[[gene]][1])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KVK[[gene]][8])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(KVK[[gene]][2:7])))
# 					,test = log(pvaluesdata[gene,"KVK"])
# 					,logLik = logLikelihood[gene,"KVK"]
# 					,AIC = AIC[gene,"KVK"]
# 					,AICc = AICc[gene,"KVK"]
# 					,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVK[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				c_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVK[[gene]][1])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KVK[[gene]][6])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(KVK[[gene]][2:5])))
# 					,test = log(pvaluesdata[gene,"KVK"])
# 					,logLik = logLikelihood[gene,"KVK"]
# 					,AIC = AIC[gene,"KVK"]
# 					,AICc = AICc[gene,"KVK"]
# 					,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVK[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			if(paramsTmp[["VKV"]]==13)
# 			{
# 				ab_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VKV[[gene]][1:6])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(VKV[[gene]][8:13])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKV[[gene]][7])))
# 					,test = log(pvaluesdata[gene,"VKV"])
# 					,logLik = logLikelihood[gene,"VKV"]
# 					,AIC = AIC[gene,"VKV"]
# 					,AICc = AICc[gene,"VKV"]
# 					,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				ab_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VKV[[gene]][1:4])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(VKV[[gene]][6:9])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKV[[gene]][5])))
# 					,test = log(pvaluesdata[gene,"VKV"])
# 					,logLik = logLikelihood[gene,"VKV"]
# 					,AIC = AIC[gene,"VKV"]
# 					,AICc = AICc[gene,"VKV"]
# 					,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKV[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			if(paramsTmp[["VVK"]]==13)
# 			{
# 				ac_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VVK[[gene]][1:6])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VVK[[gene]][13])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(VVK[[gene]][7:12])))
# 					,test = log(pvaluesdata[gene,"VVK"])
# 					,logLik = logLikelihood[gene,"VVK"]
# 					,AIC = AIC[gene,"VVK"]
# 					,AICc = AICc[gene,"VVK"]
# 					,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVK[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				ac_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VVK[[gene]][1:4])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VVK[[gene]][9])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(VVK[[gene]][5:8])))
# 					,test = log(pvaluesdata[gene,"VVK"])
# 					,logLik = logLikelihood[gene,"VVK"]
# 					,AIC = AIC[gene,"VVK"]
# 					,AICc = AICc[gene,"VVK"]
# 					,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVK[[gene]]["convergence"])
# 					,message = NULL)				
# 			}

# 			if(paramsTmp[["KVV"]]==13)
# 			{
# 				bc_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVV[[gene]][1])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(KVV[[gene]][8:13])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(KVV[[gene]][2:7])))
# 					,test = log(pvaluesdata[gene,"KVV"])
# 					,logLik = logLikelihood[gene,"KVV"]
# 					,AIC = AIC[gene,"KVV"]
# 					,AICc = AICc[gene,"KVV"]
# 					,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				bc_out = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVV[[gene]][1])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(KVV[[gene]][6:9])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(KVV[[gene]][2:5])))
# 					,test = log(pvaluesdata[gene,"KVV"])
# 					,logLik = logLikelihood[gene,"KVV"]
# 					,AIC = AIC[gene,"KVV"]
# 					,AICc = AICc[gene,"KVV"]
# 					,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVV[[gene]]["convergence"])
# 					,message = NULL)				
# 			}

# 			if(paramsTmp[["VVV"]]==18)
# 			{
# 				abc_out = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VVV[[gene]][1:6])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(VVV[[gene]][13:18])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(VVV[[gene]][7:12])))
# 					,test = log(pvaluesdata[gene,"VVV"])
# 					,logLik = logLikelihood[gene,"VVV"]
# 					,AIC = AIC[gene,"VVV"]
# 					,AICc = AICc[gene,"VVV"]
# 					,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVV[[gene]]["convergence"])
# 					,message = NULL)
# 			}else{
# 				abc_out = list(alpha = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(alpha = unname(VVV[[gene]][1:4])))
# 					  		,beta = list(fun = sigmoidModelP
# 								,type = "sigmoid"
# 								,df = 4
# 								,params = c(beta = unname(VVV[[gene]][9:12])))
#  				  		,gamma = list(fun = sigmoidModelP
# 									,type = "sigmoid"
# 									,df = 4
# 									,params = c(gamma = unname(VVV[[gene]][5:8])))
# 					,test = log(pvaluesdata[gene,"VVV"])
# 					,logLik = logLikelihood[gene,"VVV"]
# 					,AIC = AIC[gene,"VVV"]
# 					,AICc = AICc[gene,"VVV"]
# 					,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVV[[gene]]["convergence"])
# 					,message = NULL)
# 			}

# 			list("0"=zero_out,"a"=a_out,"b"=b_out,"c"=c_out,"ab"=ab_out,"ac"=ac_out,"bc"=bc_out,"abc"=abc_out)
# 		}
# 		)
# 	return(ratesSpecs)
# }




# .inspect.engine_Integrative_NoNascent_NoSigmoid <- function(tptsOriginal
# 										    , tptsLinear
# 										    , a
# 										    , c
# 										    , concentrations
# 										    , rates
# 										    , BPPARAM=bpparam()
# 										    , na.rm=TRUE
# 										    , verbose=TRUE
# 										    , testOnSmooth=TRUE
# 										    , seed = NULL
# 										    , nInit = nInit
# 										    , nIter = nIter
# 										    , limitModelComplexity = FALSE
# 										    , sigmoid = FALSE)
# {

# 	total <- concentrations$total
# 	totalVariance <- concentrations$total_var

# 	premature <- concentrations$preMRNA
# 	prematureVariance <- concentrations$preMRNA_var

# 	mature <- concentrations$mature
# 	matureVariance <- concentrations$mature_var

# 	eiGenes <- rownames(mature)

# 	k2median <- median(rates$gamma,na.rm = TRUE)
# 	k3median <- median(rates$beta,na.rm = TRUE)

# 	if( testOnSmooth ) {

# 		###### fit smooth functions on prematue and mature data
# 		###### eventually, exclude the genes that cannot be
# 		###### fit in either premature or mature

# 		matureFitImpulse <- bplapply(eiGenes,function(row)
# 		{
# 	  		fitSmooth(tpts = tptsLinear
# 	        	    , tt_c = c
# 	        	    , experiment = mature[row,]
# 	        	    , variance = matureVariance[row,]
# 	        	    , mature = TRUE
# 	        	    , nInit = nInit
# 	        	    , nIter = nIter
# 	        	    , seed = seed)
# 		},BPPARAM=BPPARAM)
# 		names(matureFitImpulse) <- eiGenes

# 		prematureFitImpulse <- bplapply(eiGenes,function(row)
# 		{
# 			fitSmooth(tpts = tptsLinear
# 					, tt_c = c
# 	        		, experiment = premature[row,]
# 	        		, variance = prematureVariance[row,]
# 	        		, mature = FALSE
# 	        		, nInit = nInit
# 	        		, nIter = nIter
# 	        		, seed = seed)     
# 		},BPPARAM=BPPARAM)
# 		names(prematureFitImpulse) <- eiGenes

# 		## in case the number of degrees of freedom allows the estimation
# 		## quantify how many genes have an acceptable chi2 test value (<0.2)
# 		## in both mature and premature
 		
#  	# 	degreesOfFreedom <- length(tptsLinear)-6
# 		# if( degreesOfFreedom > 0 ) {
# 		# 	accept <- pchisq(sapply(matureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2 &
# 		# 		pchisq(sapply(prematureFitImpulse, '[[', 'value'), degreesOfFreedom) < 0.2
# 		# 	if( table(accept)['FALSE']/length(accept) > .5 )
# 		# 		message("More than 50% did not return a good fit of their mature and
# 		# 			premature profiles with the impulsive smooth function:
# 		# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
# 		# }

# 		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
# 			stop("No genes have an expression profile possible to fit 
# 				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")

# 		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
# 			message("Some genes have an expression profile impossible to be fitted 
# 				with impulsive functions therefore they will be excluded from the modelling.")

# 			eiGenes <- eiGenes[!sapply(matureFitImpulse,is.null) & 
# 				!sapply(prematureFitImpulse,is.null)]

# 			prematureFitImpulse <- prematureFitImpulse[eiGenes]
# 			matureFitImpulse <- matureFitImpulse[eiGenes]
			
# 			total <- total[eiGenes,]
# 			totalVariance <- totalVariance[eiGenes,]
# 			premature <- premature[eiGenes,]
# 			prematureVariance <- prematureVariance[eiGenes,]
# 			mature <- mature[eiGenes,]
# 			matureVariance <- matureVariance[eiGenes,]
# 		}


# 	}

# 	KKK <- bplapply(eiGenes,function(row){
	
# 			matureParameters <- mean(mature[row,])
# 			k2Parameters <- k2median
# 			k3Parameters <- k3median

# 			unlist(
# 				tryCatch(
# 					optim(c(matureParameters, k2Parameters, k3Parameters)
# 						,errorKKK_Int_NoNascent
# 						,tpts = tptsLinear
# 						,premature = premature[row,]
# 						,mature = mature[row,]
# 						,prematureVariance = prematureVariance[row,]
# 						,matureVariance = matureVariance[row,]
# 						,control = list(maxit = nIter)),
# 					error=function(e) c(par1 = NaN
# 									  , par2 = NaN
# 									  , par3 = NaN
# 									  , value = NaN
# 									  , counts.function = NaN
# 								  	  , counts.gradient = NaN
# 								  	  , convergence = e)
# 				)
# 			)
# 		}, BPPARAM=BPPARAM)
# 		names(KKK) <- eiGenes
# 		message("Model 0 finished.")
	
# 	VKK <- bplapply(eiGenes,function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- KKK[[row]][3]
  
# 		unlist(
# 			tryCatch(
#       			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
#                   	 ,errorVKK_Int_NoNascent
#                   	 ,times = tptsOriginal
#                   	 ,data = c( mature[row,], premature[row,] )
#                   	 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#                   	 ,a = a
#                   	 ,c = c
#                   	 ,control = list(maxit = nIter)),
#     		error=function(e) c(par1 = NaN
#     						  , par2 = NaN
#     						  , par3 = NaN
#     						  , par4 = NaN
#     						  , par5 = NaN
#     						  , par6 = NaN
#     						  , par7 = NaN
#     						  , par8 = NaN
#     						  , value = NaN
#     						  , counts.function = NaN
#     						  , counts.gradient = NaN
#     						  , convergence = e)
#     		)
#     	)
# 	}, BPPARAM=BPPARAM)
# 	names(VKK) <- eiGenes
# 	message("Model A finished.")

# 	KKV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 		        	,errorKKV_Int_NoNascent
# 		        	,times = tptsOriginal
# 		        	,data = c( mature[row,], premature[row,] )
# 		        	,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 		        	,a = a
# 		        	,c = c
# 		        	,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KKV) <- eiGenes
# 	message("Model B finished.")

# 	KVK <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k3Parameters <- KKK[[row]][3]
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 					 ,errorKVK_Int_NoNascent
# 					 ,times = tptsOriginal
# 					 ,data = c( mature[row,], premature[row,] )
# 					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 					 ,a = a
# 					 ,c = c
# 					 ,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KVK) <- eiGenes
# 	message("Model C finished.")

# 	VKV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k2Parameters <- KKK[[row]][2]
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 					 ,errorVKV_Int_NoNascent
# 					 ,times = tptsOriginal
# 					 ,data = c( mature[row,], premature[row,] )
# 					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 					 ,a = a
# 					 ,c = c
# 					 ,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 		                       ,par9 = NaN
# 							   ,par10 = NaN
# 							   ,par11 = NaN
# 							   ,par12 = NaN
# 							   ,par13 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(VKV) <- eiGenes
# 	message("Model AB finished.")

# 	VVK <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k3Parameters <- KKK[[row]][3]

# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
#     				 ,errorVVK_Int_NoNascent
#     				 ,times = tptsOriginal
#     				 ,data = c( mature[row,], premature[row,] )
#     				 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
#     				 ,a = a
#     				 ,c = c
#     				 ,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							  , par2 = NaN
# 							  , par3 = NaN
# 							  , par4 = NaN
# 							  , par5 = NaN
# 							  , par6 = NaN
# 							  , par7 = NaN
# 							  , par8 = NaN
#                 			  , par9 = NaN
#                 			  , par10 = NaN
#                 			  , par11 = NaN
#                 			  , par12 = NaN
#                 			  , par13 = NaN
#                 			  , value = NaN
#                 			  , counts.function = NaN
#                 			  , counts.gradient = NaN
#                 			  , convergence = NaN)
# 			)
# 			)
# 	}, BPPARAM=BPPARAM)
# 	names(VVK) <- eiGenes
# 	message("Model AC finished.")

# 	KVV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- KKK[[row]][1]
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 					 ,errorKVV_Int_NoNascent
# 					 ,times = tptsOriginal
# 					 ,data = c( mature[row,], premature[row,] )
# 					 ,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 					 ,a = a
# 					 ,c = c
# 					 ,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,par9 = NaN
# 							   ,par10 = NaN
# 							   ,par11 = NaN
# 							   ,par12 = NaN
# 							   ,par13 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(KVV) <- eiGenes
# 	message("Model BC finished.")

# 	VVV <- bplapply(eiGenes, function(row){

# 		k1Parameters <- c(rep(k1KKK_Int_NoNascent(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
# 		k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		
# 		unlist(
# 			tryCatch(
# 				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
# 					,errorVVV_Int_NoNascent
# 					,times = tptsOriginal
# 					,data = c( mature[row,], premature[row,] )
# 					,datavar = c( matureVariance[row,] , prematureVariance[row,] )
# 					,a = a
# 					,c = c
# 					,control = list(maxit = nIter)),
# 			error=function(e) c(par1 = NaN
# 							   ,par2 = NaN
# 							   ,par3 = NaN
# 							   ,par4 = NaN
# 							   ,par5 = NaN
# 							   ,par6 = NaN
# 							   ,par7 = NaN
# 							   ,par8 = NaN
# 							   ,par9 = NaN
# 							   ,par10 = NaN
# 							   ,par11 = NaN
# 							   ,par12 = NaN
# 							   ,par13 = NaN
# 							   ,par14 = NaN
# 							   ,par15 = NaN
# 							   ,par16 = NaN
# 							   ,par17 = NaN
# 							   ,par18 = NaN
# 							   ,value = NaN
# 							   ,counts.function = NaN
# 							   ,counts.gradient = NaN
# 							   ,convergence = NaN)
# 			)
# 		)
# 	}, BPPARAM=BPPARAM)
# 	names(VVV) <- eiGenes
# 	message("Model ABC finished.")

# 	if(testOnSmooth)
# 	{
# 			matureSmooth <- t(sapply(matureFitImpulse,function(i)
# 			{
# 				impulseModel(tptsLinear,i)
# 			}))
		
# 			prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
# 			{
# 		  		impulseModel(tptsLinear,i)
# 			}))

# 			colnames(matureSmooth) <- tptsOriginal
# 			colnames(prematureSmooth) <- tptsOriginal

# 	}else{
# 		matureSmooth <- mature
# 		prematureSmooth <- premature
# 	}

# 	chi2data <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(errorKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 											  ,tptsOriginal
# 											  ,prematureSmooth[g,]
# 											  ,matureSmooth[g,]
# 											  ,prematureVariance[g,]
# 											  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(errorVKK_Int_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVKTemp <- tryCatch(errorKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(errorKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVKTemp <- tryCatch(errorVVK_Int_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VKVTemp <- tryCatch(errorVKV_Int_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		KVVTemp <- tryCatch(errorKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)

# 		VVVTemp <- tryCatch(errorVVV_Int_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
# 											  ,tptsOriginal
# 											  ,c(matureSmooth[g,],prematureSmooth[g,])
# 											  ,c(matureVariance[g,],prematureVariance[g,])
# 											  ,a
# 											  ,c),error = function(e)NaN)
	  
# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
# 	}, BPPARAM = BPPARAM))

# 	if( limitModelComplexity ) {
# 		dof <- c(KKK = 3
# 				,VKK = min(8, length(tptsOriginal)+2)
# 				,KVK = min(8, length(tptsOriginal)+2)
# 				,KKV = min(8, length(tptsOriginal)+2)
# 				,VVK = min(13, 2*length(tptsOriginal)+1)
# 				,VKV = min(13, 2*length(tptsOriginal)+1)
# 				,KVV = min(13, 2*length(tptsOriginal)+1)
# 				,VVV = min(18, 2*length(tptsOriginal))
# 				)
# 	} else {
# 	if( limitModelComplexity ) {
# 		dof <- c(KKK = 3
# 				,VKK = min(8, length(tptsOriginal)+2)
# 				,KVK = min(8, length(tptsOriginal)+2)
# 				,KKV = min(8, length(tptsOriginal)+2)
# 				,VVK = min(13, 2*length(tptsOriginal)+1)
# 				,VKV = min(13, 2*length(tptsOriginal)+1)
# 				,KVV = min(13, 2*length(tptsOriginal)+1)
# 				,VVV = min(18, 2*length(tptsOriginal))
# 				)
# 	} else {
# 		dof <- c(KKK = 3
# 				,VKK = 8
# 				,KVK = 8
# 				,KKV = 8
# 				,VVK = 13
# 				,VKV = 13
# 				,KVV = 13
# 				,VVV = 18)
# 	}
# 	}

# 	# P values
# 	pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], max(c(0,2*length(tptsOriginal)-dof['KKK'])))
# 						,VKK=pchisq(chi2data[,'VKK'], max(c(0,2*length(tptsOriginal)-dof['VKK'])))
# 						,KVK=pchisq(chi2data[,'KVK'], max(c(0,2*length(tptsOriginal)-dof['KVK'])))
# 						,KKV=pchisq(chi2data[,'KKV'], max(c(0,2*length(tptsOriginal)-dof['KKV'])))
# 						,VVK=pchisq(chi2data[,'VVK'], max(c(0,2*length(tptsOriginal)-dof['VVK'])))
# 						,VKV=pchisq(chi2data[,'VKV'], max(c(0,2*length(tptsOriginal)-dof['VKV'])))
# 						,KVV=pchisq(chi2data[,'KVV'], max(c(0,2*length(tptsOriginal)-dof['KVV'])))
# 						,VVV=pchisq(chi2data[,'VVV'], max(c(0,2*length(tptsOriginal)-dof['VVV']))))

# 	logLikelihood <- t(mcsapply(eiGenes,function(g)
# 	{
# 		KKKTemp <- tryCatch(loglikKKK_Int_NoNascent(KKK[[g]][grep("par",names(KKK[[g]]))]
# 										  ,tptsOriginal,prematureSmooth[g,]
# 										  ,matureSmooth[g,]
# 										  ,prematureVariance[g,]
# 										  ,matureVariance[g,]),error = function(e)NaN)

# 		VKKTemp <- tryCatch(loglikVKK_Int_NoNascent(VKK[[g]][grep("par",names(VKK[[g]]))]
# 										  ,tptsOriginal
# 										  ,c(matureSmooth[g,],prematureSmooth[g,])
# 										  ,c(matureVariance[g,],prematureVariance[g,])
# 										  ,a
# 										  ,c),error = function(e)NaN)

# 	  	KVKTemp <- tryCatch(loglikKVK_Int_NoNascent(KVK[[g]][grep("par",names(KVK[[g]]))]
#         				                  ,tptsOriginal
#         				                  ,c(matureSmooth[g,],prematureSmooth[g,])
#         				                  ,c(matureVariance[g,],prematureVariance[g,])
#         				                  ,a
#         				                  ,c),error = function(e)NaN)

# 		KKVTemp <- tryCatch(loglikKKV_Int_NoNascent(KKV[[g]][grep("par",names(KKV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVKTemp <- tryCatch(loglikVVK_Int_NoNascent(VVK[[g]][grep("par",names(VVK[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VKVTemp <- tryCatch(loglikVKV_Int_NoNascent(VKV[[g]][grep("par",names(VKV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		KVVTemp <- tryCatch(loglikKVV_Int_NoNascent(KVV[[g]][grep("par",names(KVV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)
# 		VVVTemp <- tryCatch(loglikVVV_Int_NoNascent(VVV[[g]][grep("par",names(VVV[[g]]))]
# 											   ,tptsOriginal
# 											   ,c(matureSmooth[g,],prematureSmooth[g,])
# 											   ,c(matureVariance[g,],prematureVariance[g,])
# 											   ,a
# 											   ,c),error = function(e)NaN)

# 		c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)

# 	}, BPPARAM=BPPARAM))

# 		AIC <- t(apply(logLikelihood,1,function(row)
# 		{
# 			2*(dof - row) 
# 		}))

# 		AICc <- t(apply(logLikelihood,1,function(row)
# 		{
# 		 	2*(dof - row) + (2*dof*(dof+1))/max(0,2*length(tptsOriginal)-dof-1)
# 		}))

# 	rownames(pvaluesdata) <- rownames(logLikelihood) <- rownames(AIC) <- rownames(AICc) <- eiGenes

# 		ratesSpecs <- lapply(eiGenes,function(gene)
# 		{
# 			list(
# 				"0" = list(alpha = list(fun = constantModelP
# 								 ,type = "constant"
# 								 ,df = 1
# 								 ,params = c(alpha = k1KKK_Int_NoNascent(x = 0,par = KKK[[gene]])))
# 					    ,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KKK[[gene]]["par3"])))
#  				   	,gamma = list(fun = constantModelP
# 								 ,type = "constant"
# 								 ,df = 1
# 								 ,params = c(gamma = unname(KKK[[gene]]["par2"])))
# 					,test = log(pvaluesdata[gene,"KKK"])
# 					,logLik = logLikelihood[gene,"KKK"]
# 					,AIC = AIC[gene,"KKK"]
# 					,AICc = AICc[gene,"KKK"]
# 					,counts = c("function"=unname(KKK[[gene]]["counts.function"]), gradient=unname(KKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKK[[gene]]["convergence"])
# 					,message = NULL)

# 				,a = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VKK[[gene]][1:6])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VKK[[gene]][8])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKK[[gene]][7])))
# 					,test = log(pvaluesdata[gene,"VKK"])
# 					,logLik = logLikelihood[gene,"VKK"]
# 					,AIC = AIC[gene,"VKK"]
# 					,AICc = AICc[gene,"VKK"]
# 					,counts = c("function"=unname(VKK[[gene]]["counts.function"]), gradient=unname(VKK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKK[[gene]]["convergence"])
# 					,message = NULL)
# 				,b = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KKV[[gene]][1])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(KKV[[gene]][3:8])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(KKV[[gene]][2])))
# 					,test = log(pvaluesdata[gene,"KKV"])
# 					,logLik = logLikelihood[gene,"KKV"]
# 					,AIC = AIC[gene,"KKV"]
# 					,AICc = AICc[gene,"KKV"]
# 					,counts = c("function"=unname(KKV[[gene]]["counts.function"]), gradient=unname(KKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KKV[[gene]]["convergence"])
# 					,message = NULL)
# 				,c = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVK[[gene]][1])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(KVK[[gene]][8])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(KVK[[gene]][2:7])))
# 					,test = log(pvaluesdata[gene,"KVK"])
# 					,logLik = logLikelihood[gene,"KVK"]
# 					,AIC = AIC[gene,"KVK"]
# 					,AICc = AICc[gene,"KVK"]
# 					,counts = c("function"=unname(KVK[[gene]]["counts.function"]), gradient=unname(KVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVK[[gene]]["convergence"])
# 					,message = NULL)
# 				,ab = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VKV[[gene]][1:6])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(VKV[[gene]][8:13])))
#  				  		,gamma = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(gamma = unname(VKV[[gene]][7])))
# 					,test = log(pvaluesdata[gene,"VKV"])
# 					,logLik = logLikelihood[gene,"VKV"]
# 					,AIC = AIC[gene,"VKV"]
# 					,AICc = AICc[gene,"VKV"]
# 					,counts = c("function"=unname(VKV[[gene]]["counts.function"]), gradient=unname(VKV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VKV[[gene]]["convergence"])
# 					,message = NULL)
# 				,ac = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VVK[[gene]][1:6])))
# 					  		,beta = list(fun = constantModelP
# 								,type = "constant"
# 								,df = 1
# 								,params = c(beta = unname(VVK[[gene]][13])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(VVK[[gene]][7:12])))
# 					,test = log(pvaluesdata[gene,"VVK"])
# 					,logLik = logLikelihood[gene,"VVK"]
# 					,AIC = AIC[gene,"VVK"]
# 					,AICc = AICc[gene,"VVK"]
# 					,counts = c("function"=unname(VVK[[gene]]["counts.function"]), gradient=unname(VVK[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVK[[gene]]["convergence"])
# 					,message = NULL)
# 				,bc = list(alpha = list(fun = constantModelP
# 									,type = "constant"
# 									,df = 1
# 									,params = c(alpha = unname(KVV[[gene]][1])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(KVV[[gene]][8:13])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(KVV[[gene]][2:7])))
# 					,test = log(pvaluesdata[gene,"KVV"])
# 					,logLik = logLikelihood[gene,"KVV"]
# 					,AIC = AIC[gene,"KVV"]
# 					,AICc = AICc[gene,"KVV"]
# 					,counts = c("function"=unname(KVV[[gene]]["counts.function"]), gradient=unname(KVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(KVV[[gene]]["convergence"])
# 					,message = NULL)
# 				,abc = list(alpha = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(alpha = unname(VVV[[gene]][1:6])))
# 					  		,beta = list(fun = impulseModelP
# 								,type = "impulse"
# 								,df = 6
# 								,params = c(beta = unname(VVV[[gene]][13:18])))
#  				  		,gamma = list(fun = impulseModelP
# 									,type = "impulse"
# 									,df = 6
# 									,params = c(gamma = unname(VVV[[gene]][7:12])))
# 					,test = log(pvaluesdata[gene,"VVV"])
# 					,logLik = logLikelihood[gene,"VVV"]
# 					,AIC = AIC[gene,"VVV"]
# 					,AICc = AICc[gene,"VVV"]
# 					,counts = c("function"=unname(VVV[[gene]]["counts.function"]), gradient=unname(VVV[[gene]]["counts.gradient"]))
# 					,convergence = unname(VVV[[gene]]["convergence"])
# 					,message = NULL)
# 		)
# 		})

# 	return(ratesSpecs)
# }



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





































### Choose among constant, sigmoid and impulsive function ###
	.chooseModel <- function(tpts
						   , experiment
						   , variance=NULL
						   , na.rm=TRUE
						   , sigmoid=TRUE
						   , impulse=TRUE
						   , polynomial=TRUE
						   , nInit=10
						   , nIter=500
						   , impulseModel
						   , sigmoidModel
						   , sigmoidModelP
						   , impulseModelP
						   , .polynomialModelP
						   , seed = 1
						   , computeDerivatives = TRUE)
	{
		if(is.null(seed))seed <- 1
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
			 
			 D <- abs(.DimpulseModel(tpts, par))
			 D2 <- abs(.D2impulseModel(tpts, par))
			 
			 D0 <- .DimpulseModel(0, par)

			 if(!is.finite(D0))return(NaN)

			 if(computeDerivatives)
			 {
				 if(all(is.finite(D)) & all(is.finite(D2))) return(chisqFunction(experiment, model, variance)+abs(D0)+max(D)+max(D2))
				 else(return(NaN))	 	
			 }else{
			 	return(chisqFunction(experiment, model, variance)+abs(D0))
			 }
		}
		#
		im.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
			, maxit=500){
			set.seed(seed) 
			sapply(1:ninit, function(x) 
		 		tryCatch(optim(
		 			par=im.parguess(tpts, experiment)
		 			, fn=im.chisq, tpts=tpts
		 			, experiment=experiment
		 			, variance=variance
		 			, impulseModel=impulseModel
		 			, control=list(maxit=maxit)
		 			), error=function(e) optimFailOut(e)))}
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
			 D <- abs(.DsigmoidModel(tpts, par))
			 D2 <- abs(.D2sigmoidModel(tpts, par))

	 		 D0 <- .DsigmoidModel(0, par)

			 if(!is.finite(D0))return(NaN)

			 if(computeDerivatives)
			 {
				 if(all(is.finite(D)) & all(is.finite(D2))) return(chisqFunction(experiment, model, variance)+abs(D0)+max(D)+max(D2))
				 else(return(NaN))
			 }else{
			 	return(chisqFunction(experiment, model, variance)+abs(D0))
			 }
		}
		#
		sm.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
			, maxit=500){ 
			set.seed(seed)
			sapply(1:ninit, function(x) 
					tryCatch(optim(
						par=sm.parguess(tpts, experiment)
						, fn=sm.chisq, tpts=tpts
						, experiment=experiment
						, variance=variance
						, sigmoidModel=sigmoidModel
						, control=list(maxit=maxit)
						), error=function(e) optimFailOut(e)))
		}
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
		## re-evaluate flags of function to evaluate according to the length 
		## of the experiment
		sigmoid <- sigmoid
		impulse <- impulse & length(experiment)>2
		polynomial <- polynomial & length(experiment)>2
		# sigmoid
		if( sigmoid ) {
			outSM  <- sm.optim.chisq(tpts=tpts, experiment=experiment
				, variance=variance, ninit=nInit, maxit=nIter)
			bestSM <- which.min(unlist(outSM[2,]))
			pvalSM <- tryCatch(chisq.test.default(experiment=experiment
				, model=sigmoidModel(tpts, outSM[,bestSM]$par)
				, variance=variance, df=length(outSM[,bestSM]$par)),error=function(e)NaN)
			dfSM <- length(outSM[,bestSM]$par)
		} else dfSM <- NA
		# impulse
		if( impulse ) {
			outIM  <- im.optim.chisq(tpts=tpts, experiment=experiment, 
				variance=variance, ninit=nInit, maxit=nIter)
			bestIM <- which.min(unlist(outIM[2,]))
			pvalIM <- tryCatch(chisq.test.default(experiment=experiment
				, model=impulseModel(tpts, outIM[,bestIM]$par) 
				, variance=variance, df=length(outIM[,bestIM]$par)),error=function(e)NaN)
			dfIM <- length(outIM[,bestIM]$par)
		} else dfIM <- NA
		# polynomial
		if( polynomial ) {
			outPN  <- pn.optim.aic(tpts, experiment, variance )
			bestPN <- which.min(unlist(outPN[2,]))
			pvalPN <- tryCatch(chisq.test.default(experiment=experiment
				, model=.polynomialModel(tpts, outPN[,bestPN]$par) 
				, variance=variance, df=length(outPN[,bestPN]$par)),error=function(e)NaN)
			dfPN <- length(outPN[,bestPN]$par)
		} else dfPN <- NA
		pvals <- c(
			sigmoid=if( sigmoid ) pvalSM else NA
			, impulse=if( impulse ) pvalIM else NA
			, polynomial=if( polynomial ) pvalPN else NA
			)

		# if(pvals["impulse"]>0.05&is.finite(pvals["sigmoid"])){pvals["impulse"] <- NA} # I prefer the sigmoid function if the impulse is not good enough

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

### Numerical solution of the complete ODE system both for unlabled and labled RNA ###
	systemSolution <- function(k1F,k2F,k3F,times,nascent=FALSE)
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

	  if(nascent){cinit <- c(0,0)}else{cinit <- c(k1F(0)/k2F(0),k1F(0)/k3F(0))}
	  names(cinit) <- c("p","m")

	  params <- list(alpha = k1F, beta = k2F, gamma = k3F)

	  modData <- ode(y=cinit, times=times, func=system, parms=params)
	  modData <- c(modData[,"m"],modData[,"p"])
	  
	  if(length(modData)!=2*length(times)){modData <- rep(NaN,2*length(times))}

	  names(modData) <- c(rep("mature",length.out = length(times)),rep("premature",length.out = length(times)))

	  return(modData) 
	}

	systemSolution_Simple <- function(k1F,k3F,times,nascent=FALSE)
	{

	  system <- function(t,c,parms)
	  {
	    alpha <- parms$alpha
	    gamma <- parms$gamma

	    r=rep(0,length(c))
	    r[1] <- alpha(t) - gamma(t) * c["t"]
	    return(list(r))
	  }
		
	  if(nascent){cinit <- 0}else{cinit <- k1F(0)/k3F(0)}
	  names(cinit) <- c("t")

	  params <- list(alpha = k1F, gamma = k3F)

	  modData <- ode(y=cinit, times=times, func=system, parms=params)
	  modData <- c(modData[,"t"])

	  ### POSSIBLE REQUEST FOR A CHECK LIKE:
	  #
	  # if(length(modData)!=2*length(times)){modData <- rep(NaN,2*length(times))}

	  names(modData) <- rep("total",length.out = length(times))

	  return(modData) 
	}

### Genesis of an object reporting expression data, kinetic rates and the temporal profile ###
	.makeModel_Derivative <- function(tpts, hyp, geneBestModel)
	{
		if(!any(geneBestModel %in% c("b","c","bc")))
		{
			params <- list()
			params$mature <- function(x) 
				hyp$mature$fun$value(x, hyp$mature$par)
			params$beta  <- function(x) 
				hyp$beta$fun$value(x, hyp$beta$par)
			params$gamma <- function(x) 
				hyp$gamma$fun$value(x, hyp$gamma$par)

			matureTemp <- params$mature(tpts)
			k2Temp <- params$gamma(tpts)
			k3Temp <- params$beta(tpts)

			if(geneBestModel == "0")
			{
				prematureTemp <- sapply(tpts,function(t)prematureKKK_Der(t,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params)))
				k1Temp <- sapply(tpts,function(t)k1KKK_Der(t,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params)))
			}else if(geneBestModel == "a")
			{
				prematureTemp <- prematureVKK_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
				k1Temp <- k1VKK_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
			}else if(geneBestModel == "ac")
			{
				prematureTemp <- prematureVVK_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
				k1Temp <- k1VVK_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
			}else if(geneBestModel == "ab")
			{
				prematureTemp <- prematureVKV_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
				k1Temp <- k1VKV_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
			}else if(geneBestModel == "abc")
			{
				prematureTemp <- prematureVVV_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
				k1Temp <- k1VVV_Der(tpts,c(hyp$mature$params,hyp$gamma$params,hyp$beta$params))
			}
			totalTemp <- matureTemp + prematureTemp
		}

		if(geneBestModel == "b")
		{
			params <- list()
			params$total <- function(x) 
				hyp$total$fun$value(x, hyp$total$par)
			params$alpha  <- function(x) 
				hyp$alpha$fun$value(x, hyp$alpha$par)
			params$gamma <- function(x) 
				hyp$gamma$fun$value(x, hyp$gamma$par)

			totalTemp <- params$total(tpts)
			prematureTemp <- sapply(tpts, function(t) prematureKKV_Der(t, c(hyp$total$params,hyp$alpha$params,hyp$gamma$params)))

			k1Temp <- params$alpha(tpts)
			k2Temp <- params$gamma(tpts)
			k3Temp <- sapply(tpts, function(t) k3KKV_Der(t, c(hyp$total$params,hyp$alpha$params,hyp$gamma$params)))
		}else if(geneBestModel == "c"){
			params <- list()
			params$total <- function(x) 
				hyp$total$fun$value(x, hyp$total$par)
			params$alpha  <- function(x) 
				hyp$alpha$fun$value(x, hyp$alpha$par)
			params$beta <- function(x) 
				hyp$beta$fun$value(x, hyp$beta$par)

			totalTemp <- params$total(tpts)
			prematureTemp <- sapply(tpts, function(t) prematureKVK_Der(t, c(hyp$total$params,hyp$alpha$params,hyp$beta$params)))

			k1Temp <- params$alpha(tpts)
			k3Temp <- params$beta(tpts)
			k2Temp <- sapply(tpts, function(t) k2KVK_Der(t, c(hyp$total$params,hyp$alpha$params,hyp$beta$params)))
		}else if(geneBestModel == "bc"){
			params <- list()
			params$total <- function(x) 
				hyp$total$fun$value(x, hyp$total$par)
			params$alpha  <- function(x) 
				hyp$alpha$fun$value(x, hyp$alpha$par)
			params$beta <- function(x) 
				hyp$beta$fun$value(x, hyp$beta$par)

			totalTemp <- params$total(tpts)
			prematureTemp <- sapply(tpts, function(t) prematureKVV_Der(t, c(hyp$total$params,hyp$alpha$params,hyp$beta$params)))

			k1Temp <- params$alpha(tpts)
			k3Temp <- params$beta(tpts)
			k2Temp <- sapply(tpts, function(t) k2KVV_Der(t, c(hyp$total$params,hyp$alpha$params,hyp$beta$params)))
		}

		data.frame(time = tpts
				 , preMRNA = prematureTemp
				 , total = totalTemp
				 , alpha = k1Temp
				 , beta = k3Temp
				 , gamma = k2Temp)
	}

	.makeModel_Derivative_Simple <- function(tpts, hyp, geneBestModel)
	{
		if(geneBestModel != "b")
		{
			params <- list()
			params$total <- function(x) 
				hyp$total$fun$value(x, hyp$total$par)
			params$beta  <- function(x) 
				hyp$beta$fun$value(x, hyp$beta$par)

			totalTemp <- params$total(tpts)
			k3Temp <- params$beta(tpts)

			if(geneBestModel == "0")
			{
				k1Temp <- sapply(tpts,function(t)k1KKK_Der_Simple(t,c(hyp$total$params,hyp$beta$params)))
			}else if(geneBestModel == "a")
			{
				k1Temp <- k1VKK_Der_Simple(tpts,c(hyp$total$params,hyp$beta$params))
			}else if(geneBestModel == "ab")
			{
				k1Temp <- k1VKV_Der_Simple(tpts,c(hyp$total$params,hyp$beta$params))
			}		
		}else{
			params <- list()
			params$total <- function(x) 
				hyp$total$fun$value(x, hyp$total$par)
			params$alpha  <- function(x) 
				hyp$alpha$fun$value(x, hyp$alpha$par)

			totalTemp <- params$total(tpts)

			k1Temp <- params$alpha(tpts)
			k3Temp <- sapply(tpts, function(t) k3KKV_Der_Simple(t, c(hyp$total$params,hyp$alpha$params)))
		}

		data.frame(time = tpts
		 , preMRNA = unname(rep(NaN,length(tpts)))
		 , total = totalTemp
		 , alpha = k1Temp
		 , beta = k3Temp
		 , gamma = unname(rep(NaN,length(tpts))))
	}

#############################################################################
### NEW INSPEcT FUNCTIONS FOR THE COMPUTATION ON THE CONFIDENCE INTERVALS ###
#############################################################################

logLikelihoodCIerror <- function(parameter,name,parameters,class,tpts,experimentalP,experimentalM,experimentalA=NULL,varianceP,varianceM,varianceA=NULL,confidenceThreshold,derivative=TRUE)
{
	if(derivative)
	{
		maximumLogLikelihoodTmp <- logLikelihood_derivativeModels(tpts=tpts
																 ,class=class
																 ,parameters=parameters
																 ,premature=experimentalP
																 ,mature=experimentalM
																 ,alpha=experimentalA
																 ,prematureVariance=varianceP
																 ,matureVariance=varianceM
																 ,alphaVariance=varianceA)

		perturbedParameters <- parameters
		perturbedParameters[name] <- parameter
		
		perturbedLogLikelihoodTmp <- logLikelihood_derivativeModels(tpts=tpts
																   ,class=class
																   ,parameters=perturbedParameters
																   ,premature=experimentalP
																   ,mature=experimentalM
																   ,alpha=experimentalA
																   ,prematureVariance=varianceP
																   ,matureVariance=varianceM
																   ,alphaVariance=varianceA)
	}else{
		maximumLogLikelihoodTmp <- logLikelihood_integrativeModels(tpts=tpts
																 ,class=class
																 ,parameters=parameters
																 ,premature=experimentalP
																 ,mature=experimentalM
																 ,alpha=experimentalA
																 ,prematureVariance=varianceP
																 ,matureVariance=varianceM
																 ,alphaVariance=varianceA)

		perturbedParameters <- parameters
		perturbedParameters[name] <- parameter
		
		perturbedLogLikelihoodTmp <- logLikelihood_integrativeModels(tpts=tpts
																   ,class=class
																   ,parameters=perturbedParameters
																   ,premature=experimentalP
																   ,mature=experimentalM
																   ,alpha=experimentalA
																   ,prematureVariance=varianceP
																   ,matureVariance=varianceM
																   ,alphaVariance=varianceA)
	}

	return(abs(confidenceThreshold - 2*(maximumLogLikelihoodTmp - perturbedLogLikelihoodTmp)))
}

logLikelihood_derivativeModels <- function(tpts,class,parameters,premature,mature,alpha,prematureVariance,matureVariance,alphaVariance)
{
	if(class=="0")
	{
		prematureKKKTemp <- c(sapply(seq_along(tpts),function(t)prematureKKK_Der(x = tpts[t], parameters = parameters)))
		matureKKKTemp <- rep(parameters[[1]],length(tpts))
		
		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelKKK <- c(matureKKKTemp,prematureKKKTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelKKK, variance = c(matureVariance,prematureVariance)))
		}

		alphaKKKTemp <- c(sapply(seq_along(tpts),function(t)k1KKK_Der(x = tpts[t], parameters = parameters)))
		modelKKK <- c(matureKKKTemp,prematureKKKTemp,alphaKKKTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha), model = modelKKK, variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else if(class=="a")
	{
		prematureVKKTemp <- c(sapply(seq_along(tpts),function(t)prematureVKK_Der(x = tpts[t], parameters = parameters)))
		matureVKKTemp <- c(sapply(seq_along(tpts),function(t)matureVKK_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelVKK <- c(matureVKKTemp,prematureVKKTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelVKK, variance = c(matureVariance,prematureVariance)))
		}

		alphaVKKTemp <- c(sapply(seq_along(tpts),function(t)k1VKK_Der(x = tpts[t], parameters = parameters)))
		modelVKK <- c(matureVKKTemp,prematureVKKTemp,alphaVKKTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelVKK , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else if(class=="c")
	{
		prematureKVKTemp <- c(sapply(seq_along(tpts),function(t)prematureKVK_Der(x = tpts[t], parameters = parameters)))
		matureKVKTemp <- c(sapply(seq_along(tpts),function(t)matureKVK_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelKVK <- c(matureKVKTemp,prematureKVKTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelKVK, variance = c(matureVariance,prematureVariance)))
		}

		alphaKVKTemp <- c(sapply(seq_along(tpts),function(t)k1KVK_Der(x = tpts[t], parameters = parameters)))
		modelKVK <- c(matureKVKTemp,prematureKVKTemp,alphaKVKTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelKVK , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else if(class=="b")
	{
		prematureKKVTemp <- c(sapply(seq_along(tpts),function(t)prematureKKV_Der(x = tpts[t], parameters = parameters)))
		matureKKVTemp <- c(sapply(seq_along(tpts),function(t)matureKKV_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelKKV <- c(matureKKVTemp,prematureKKVTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelKKV, variance = c(matureVariance,prematureVariance)))
		}

		alphaKKVTemp <- c(sapply(seq_along(tpts),function(t)k1KKV_Der(x = tpts[t], parameters = parameters)))
		modelKKV <- c(matureKKVTemp,prematureKKVTemp,alphaKKVTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelKKV , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else if(class=="ac")
	{
		prematureVVKTemp <- c(sapply(seq_along(tpts),function(t)prematureVVK_Der(x = tpts[t], parameters = parameters)))
		matureVVKTemp <- c(sapply(seq_along(tpts),function(t)matureVVK_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelVVK <- c(matureVVKTemp,prematureVVKTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelVVK, variance = c(matureVariance,prematureVariance)))
		}

		alphaVVKTemp <- c(sapply(seq_along(tpts),function(t)k1VVK_Der(x = tpts[t], parameters = parameters)))
		modelVVK <- c(matureVVKTemp,prematureVVKTemp,alphaVVKTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelVVK , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else if(class=="ab")
	{
		prematureVKVTemp <- c(sapply(seq_along(tpts),function(t)prematureVKV_Der(x = tpts[t], parameters = parameters)))
		matureVKVTemp <- c(sapply(seq_along(tpts),function(t)matureVKV_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelVKV <- c(matureVKVTemp,prematureVKVTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelVKV, variance = c(matureVariance,prematureVariance)))
		}

		alphaVKVTemp <- c(sapply(seq_along(tpts),function(t)k1VKV_Der(x = tpts[t], parameters = parameters)))
		modelVKV <- c(matureVKVTemp,prematureVKVTemp,alphaVKVTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelVKV , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else if(class=="bc")
	{
		prematureKVVTemp <- c(sapply(seq_along(tpts),function(t)prematureKVV_Der(x = tpts[t], parameters = parameters)))
		matureKVVTemp <- c(sapply(seq_along(tpts),function(t)matureKVV_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelKVV <- c(matureKVVTemp,prematureKVVTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelKVV, variance = c(matureVariance,prematureVariance)))
		}

		alphaKVVTemp <- c(sapply(seq_along(tpts),function(t)k1KVV_Der(x = tpts[t], parameters = parameters)))
		modelKVV <- c(matureKVVTemp,prematureKVVTemp,alphaKVVTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelKVV , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}else{
		prematureVVVTemp <- c(sapply(seq_along(tpts),function(t)prematureVVV_Der(x = tpts[t], parameters = parameters)))
		matureVVVTemp <- c(sapply(seq_along(tpts),function(t)matureVVV_Der(x = tpts[t], parameters = parameters)))

		if(is.null(alpha)&is.null(alphaVariance))
		{
			modelVVV <- c(matureVVVTemp,prematureVVVTemp)
			return(logLikelihoodFunction(experiment = c(mature,premature), model = modelVVV, variance = c(matureVariance,prematureVariance)))
		}

		alphaVVVTemp <- c(sapply(seq_along(tpts),function(t)k1VVV_Der(x = tpts[t], parameters = parameters)))
		modelVVV <- c(matureVVVTemp,prematureVVVTemp,alphaVVVTemp)
		return(logLikelihoodFunction(experiment = c(mature,premature,alpha) , model = modelVVV , variance = c(matureVariance,prematureVariance,alphaVariance)))
	}
}

logLikelihood_integrativeModels <- function(tpts,class,parameters,premature,mature,alpha,prematureVariance,matureVariance,alphaVariance)
{
	if(class=="KKK"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3]

		k1F <- function(x) {k1Parameters}
		k2F <- function(x) {k2Parameters}
		k3F <- function(x) {k3Parameters}
	}else if(class=="VKK"){
		k1Parameters <- parameters[1:(length(parameters)-2)]
		k2Parameters <- parameters[length(parameters)-1]
		k3Parameters <- parameters[length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		k2F <- function(x) {k2Parameters}
		k3F <- function(x) {k3Parameters}

	}else if(class=="KKV"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:length(parameters)]

		k1F <- function(x) {k1Parameters}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		k3F <- function(x) {k3Parameters}
	}else if(class=="KVK"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2:(length(parameters)-1)]
		k3Parameters <- parameters[length(parameters)]

		k1F <- function(x) {k1Parameters}
		k2F <- function(x) {k2Parameters}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}else if(class=="VKV"){
		k1Parameters <- parameters[1:((length(parameters)-1)/2)]
		k2Parameters <- parameters[1+(length(parameters)-1)/2]
		k3Parameters <- parameters[(2+(length(parameters)-1)/2):length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		k3F <- function(x) {k3Parameters}
	}else if(class=="VVK"){
		k1Parameters <- parameters[1:((length(parameters)-1)/2)]
		k2Parameters <- parameters[(1+(length(parameters)-1)/2):(length(parameters)-1)]
		k3Parameters <- parameters[length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		k2F <- function(x) {k2Parameters}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}else if(class=="KVV"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2:(1+(length(parameters)-1)/2)]
		k3Parameters <- parameters[(1+(1+(length(parameters)-1)/2)):length(parameters)]

		k1F <- function(x) {k1Parameters}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}else{
		k1Parameters <- parameters[1:(length(parameters)/3)]
		k2Parameters <- parameters[(1+length(parameters)/3):(2*(length(parameters)/3))]
		k3Parameters <- parameters[(1+2*(length(parameters)/3)):length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}

	modData <- systemSolution(k1F,k2F,k3F,tpts)

	if(!is.null(alpha))modData <- c(modData,sapply(tpts,k1F))

	return(logLikelihoodFunction(experiment = c(mature,premature,alpha), model = modData, variance = c(matureVariance,prematureVariance,alphaVariance)))
}

rates_derivativeModels <- function(tpts,class,parameters)
{
	if(class=="0")
	{
		alphaKKKTemp <- c(sapply(seq_along(tpts),function(t)k1KKK_Der(x = tpts[t], parameters = parameters)))
		betaKKKTemp <- c(sapply(seq_along(tpts),function(t)k3KKK_Der(x = tpts[t], parameters = parameters)))
		gammaKKKTemp <- c(sapply(seq_along(tpts),function(t)k2KKK_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaKKKTemp), beta=unname(betaKKKTemp), gamma=unname(gammaKKKTemp)))

	}else if(class=="a")
	{
		alphaVKKTemp <- c(sapply(seq_along(tpts),function(t)k1VKK_Der(x = tpts[t], parameters = parameters)))
		betaVKKTemp <- c(sapply(seq_along(tpts),function(t)k3VKK_Der(x = tpts[t], parameters = parameters)))
		gammaVKKTemp <- c(sapply(seq_along(tpts),function(t)k2VKK_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaVKKTemp), beta=unname(betaVKKTemp), gamma=unname(gammaVKKTemp)))
	
	}else if(class=="c")
	{
		alphaKVKTemp <- c(sapply(seq_along(tpts),function(t)k1KVK_Der(x = tpts[t], parameters = parameters)))
		betaKVKTemp <- c(sapply(seq_along(tpts),function(t)k3KVK_Der(x = tpts[t], parameters = parameters)))
		gammaKVKTemp <- c(sapply(seq_along(tpts),function(t)k2KVK_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaKVKTemp), beta=unname(betaKVKTemp), gamma=unname(gammaKVKTemp)))
	
	}else if(class=="b")
	{
		alphaKKVTemp <- c(sapply(seq_along(tpts),function(t)k1KKV_Der(x = tpts[t], parameters = parameters)))
		betaKKVTemp <- c(sapply(seq_along(tpts),function(t)k3KKV_Der(x = tpts[t], parameters = parameters)))
		gammaKKVTemp <- c(sapply(seq_along(tpts),function(t)k2KKV_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaKKVTemp), beta=unname(betaKKVTemp), gamma=unname(gammaKKVTemp)))
	
	}else if(class=="ac")
	{
		alphaVVKTemp <- c(sapply(seq_along(tpts),function(t)k1VVK_Der(x = tpts[t], parameters = parameters)))
		betaVVKTemp <- c(sapply(seq_along(tpts),function(t)k3VVK_Der(x = tpts[t], parameters = parameters)))
		gammaVVKTemp <- c(sapply(seq_along(tpts),function(t)k2VVK_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaVVKTemp), beta=unname(betaVVKTemp), gamma=unname(gammaVVKTemp)))
	
	}else if(class=="ab")
	{
		alphaVKVTemp <- c(sapply(seq_along(tpts),function(t)k1VKV_Der(x = tpts[t], parameters = parameters)))
		betaVKVTemp <- c(sapply(seq_along(tpts),function(t)k3VKV_Der(x = tpts[t], parameters = parameters)))
		gammaVKVTemp <- c(sapply(seq_along(tpts),function(t)k2VKV_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaVKVTemp), beta=unname(betaVKVTemp), gamma=unname(gammaVKVTemp)))
	
	}else if(class=="bc")
	{
		alphaKVVTemp <- c(sapply(seq_along(tpts),function(t)k1KVV_Der(x = tpts[t], parameters = parameters)))
		betaKVVTemp <- c(sapply(seq_along(tpts),function(t)k3KVV_Der(x = tpts[t], parameters = parameters)))
		gammaKVVTemp <- c(sapply(seq_along(tpts),function(t)k2KVV_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaKVVTemp), beta=unname(betaKVVTemp), gamma=unname(gammaKVVTemp)))
	
	}else
	{
		alphaVVVTemp <- c(sapply(seq_along(tpts),function(t)k1VVV_Der(x = tpts[t], parameters = parameters)))
		betaVVVTemp <- c(sapply(seq_along(tpts),function(t)k3VVV_Der(x = tpts[t], parameters = parameters)))
		gammaVVVTemp <- c(sapply(seq_along(tpts),function(t)k2VVV_Der(x = tpts[t], parameters = parameters)))

		return(c(alpha=unname(alphaVVVTemp), beta=unname(betaVVVTemp), gamma=unname(gammaVVVTemp)))
	
	}
}

rates_integrativeModels <- function(tpts, class, parameters)
{
	rate <- function(tpts,parameters){if(length(parameters)==6){return(impulseModel(tpts,parameters))}
									  else if(length(parameters)==4){return(sigmoidModel(tpts,parameters))}
									  else{return(rep(parameters,length(tpts)))}}
	if(class=="KKK")
	{
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3]
	}else if(class=="VKK")
	{
		k1Parameters <- parameters[1:(length(parameters)-2)]
		k2Parameters <- parameters[length(parameters)-1]
		k3Parameters <- parameters[length(parameters)]
	}else if(class=="KKV")
	{
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:length(parameters)]
	}else if(class=="KVK")
	{
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2:(length(parameters)-1)]
		k3Parameters <- parameters[length(parameters)]
	}else if(class=="VKV")
	{
		k1Parameters <- parameters[1:((length(parameters)-1)/2)]
		k2Parameters <- parameters[1+(length(parameters)-1)/2]
		k3Parameters <- parameters[(2+(length(parameters)-1)/2):length(parameters)]
	}else if(class=="VVK")
	{
		k1Parameters <- parameters[1:((length(parameters)-1)/2)]
		k2Parameters <- parameters[(1+(length(parameters)-1)/2):(length(parameters)-1)]
		k3Parameters <- parameters[length(parameters)]
	}else if(class=="KVV")
	{
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2:(1+(length(parameters)-1)/2)]
		k3Parameters <- parameters[(1+(1+(length(parameters)-1)/2)):length(parameters)]
	}else if(class=="VVV")
	{
		k1Parameters <- parameters[1:(length(parameters)/3)]
		k2Parameters <- parameters[(1+length(parameters)/3):(2*(length(parameters)/3))]
		k3Parameters <- parameters[(1+2*(length(parameters)/3)):length(parameters)]
	}

	return(c("alpha"=rate(tpts,k1Parameters),"gamma"=rate(tpts,k2Parameters),"beta"=rate(tpts,k3Parameters)))
}

expressionData_integrativeModels <- function(tpts,class,parameters)
{
	if(class=="KKK"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3]

		k1F <- function(x) {k1Parameters}
		k2F <- function(x) {k2Parameters}
		k3F <- function(x) {k3Parameters}
	}else if(class=="VKK"){
		k1Parameters <- parameters[1:(length(parameters)-2)]
		k2Parameters <- parameters[length(parameters)-1]
		k3Parameters <- parameters[length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		k2F <- function(x) {k2Parameters}
		k3F <- function(x) {k3Parameters}

	}else if(class=="KKV"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:length(parameters)]

		k1F <- function(x) {k1Parameters}
		k2F <- function(x) {k2Parameters}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}else if(class=="KVK"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2:(length(parameters)-1)]
		k3Parameters <- parameters[length(parameters)]

		k1F <- function(x) {k1Parameters}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		k3F <- function(x) {k3Parameters}
	}else if(class=="VKV"){
		k1Parameters <- parameters[1:((length(parameters)-1)/2)]
		k2Parameters <- parameters[1+(length(parameters)-1)/2]
		k3Parameters <- parameters[(2+(length(parameters)-1)/2):length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		k2F <- function(x) {k2Parameters}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}else if(class=="VVK"){
		k1Parameters <- parameters[1:((length(parameters)-1)/2)]
		k2Parameters <- parameters[(1+(length(parameters)-1)/2):(length(parameters)-1)]
		k3Parameters <- parameters[length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		k3F <- function(x) {k3Parameters}
	}else if(class=="KVV"){
		k1Parameters <- parameters[1]
		k2Parameters <- parameters[2:(1+(length(parameters)-1)/2)]
		k3Parameters <- parameters[(1+(1+(length(parameters)-1)/2)):length(parameters)]

		k1F <- function(x) {k1Parameters}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}else{
		k1Parameters <- parameters[1:(length(parameters)/3)]
		k2Parameters <- parameters[(1+length(parameters)/3):(2*(length(parameters)/3))]
		k3Parameters <- parameters[(1+2*(length(parameters)/3)):length(parameters)]

		if(length(k1Parameters)==6){k1F <- function(x){impulseModel(x, k1Parameters)}}else{k1F <- function(x){sigmoidModel(x, k1Parameters)}}
		if(length(k2Parameters)==6){k2F <- function(x){impulseModel(x, k2Parameters)}}else{k2F <- function(x){sigmoidModel(x, k2Parameters)}}
		if(length(k3Parameters)==6){k3F <- function(x){impulseModel(x, k3Parameters)}}else{k3F <- function(x){sigmoidModel(x, k3Parameters)}}
	}

	modData <- systemSolution(k1F,k2F,k3F,tpts)

	return(modData)
}

emptyList <- list(root = NaN, f.root = NaN, iter = NaN, estim.precis = NaN)

k_score_fun <- function(k, rate_conf_int)
{
	sum(apply(rate_conf_int, 1, function(x) {
		if( k < x[2] ) (k - x[2])^2/(x[2]-x[1])^2 
			else (k - x[2])^2/(x[2]-x[3])^2
	}), na.rm=T)
}

#####################################################################
################ NEW INSPEcT INTEGRATIVE APPROACHES #################
#####################################################################

.inspect.engine_Integrative_Nascent <- function(tpts
											 , concentrations
											 , rates
											 , BPPARAM
											 , na.rm
											 , verbose
											 , testOnSmooth
											 , seed
											 , nInit
											 , nIter
											 , limitModelComplexity
											 , computeDerivatives = TRUE
											 , useSigmoidFun = TRUE
											 , initialPenalityRelevance = 1
											 , derivativePenalityRelevance = 10^-50
											 , llConfidenceThreshold)
{
	total <- concentrations$total
 	totalVariance <- concentrations$total_var

 	premature <- concentrations$preMRNA
 	prematureVariance <- concentrations$preMRNA_var

 	mature <- concentrations$mature
 	matureVariance <- concentrations$mature_var

 	alpha <- rates$alpha
 	alphaVariance <- rates$alpha_var

 	beta <- rates$beta
 	gamma <- rates$gamma

	prematureSmooth <- premature
	matureSmooth <- mature

	eiGenes <- rownames(mature)

 	KKK <- bplapply(eiGenes,function(row){

 			k1Parameters <- mean(alpha[row,])
 			k2Parameters <- mean(gamma[row,])
 			k3Parameters <- mean(beta[row,])

 			unlist(
 				tryCatch(
 					optim(c(k1Parameters, k2Parameters, k3Parameters)
 						,errorKKK_Int
 						,tpts = tpts
 						,premature = premature[row,]
 						,mature = mature[row,]
 						,alpha = alpha[row,]
 						,prematureVariance = prematureVariance[row,]
 						,matureVariance = matureVariance[row,]
 						,alphaVariance = alphaVariance[row,]
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

 	VVV <- bplapply(eiGenes, function(row){
 		if(useSigmoidFun)
 		{
			k1Parameters <- c(rep(KKK[[row]][1],2), max(tpts)/3,1)
			k2Parameters <- c(rep(KKK[[row]][2],2), max(tpts)/3,1)
			k3Parameters <- c(rep(KKK[[row]][3],2), max(tpts)/3,1)

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorVVV_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = alpha[row,]
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = alphaVariance[row,]	
						 ,KKK = NULL
						 ,initialChisquare = NULL
						 ,initialDistances = NULL
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
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
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[1:4],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)})
			k2Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[5:8],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)})
			k3Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[9:12],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)})

 		}else{
 			k1Parameters <- c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)
 		}

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorVVV_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = alpha[row,]
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = alphaVariance[row,]	
					 ,KKK = NULL
					 ,initialChisquare = NULL
					 ,initialDistances = NULL
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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
							  , par14 = NaN
							  , par15 = NaN
							  , par16 = NaN
							  , par17 = NaN
							  , par18 = NaN
							  , value = NaN
							  , counts.function = NaN
							  , counts.gradient = NaN
							  , convergence = e)
			)
		)

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,3*length(tpts) - 18)) < pchisq(sigmoidsParameters["value"],max(0,3*length(tpts) - 12))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)

	names(VVV) <- eiGenes

	### Standard outputs

	# Log likelihood
	logLikelihood <- t(sapply(eiGenes,function(g)
	{
		modelVVV <- expressionData_integrativeModels(tpts, class = "VVV", parameters = VVV[[g]][grep("par",names(VVV[[g]]))])
		ratesVVV <- rates_integrativeModels(tpts, class = "VVV", parameters = VVV[[g]][grep("par",names(VVV[[g]]))])
		
		matureModel <- modelVVV[grep("^mature",names(modelVVV))]
		prematureModel <- modelVVV[grep("^premature",names(modelVVV))]
		alphaModel <- ratesVVV[grep("alpha",names(ratesVVV))]

		modelVVV <- c(matureModel,prematureModel,alphaModel)

		VVVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,],alpha[g,])
		                               , model = modelVVV
		                               , variance = c(matureVariance[g,],prematureVariance[g,],alphaVariance[g,])),error=function(e)NaN)

		c("KKK" = NaN,"VKK" = NaN,"KVK" = NaN,"KKV" = NaN,"VVK" = NaN,"VKV" = NaN,"KVV" = NaN,"VVV" = VVVTemp)
	}))

	rownames(logLikelihood) <- eiGenes

	dof <- cbind(KKK = NaN
				,VKK = NaN
				,KVK = NaN
				,KKV = NaN
				,VVK = NaN
				,VKV = NaN
				,KVV = NaN
				,VVV = sapply(VVV,function(m)length(grep("par",names(m)))))

 	AIC <- 2*(dof - logLikelihood)
	AICc <- 2*(dof - logLikelihood) + (2*dof*(dof+1))/max(0,2*length(tpts)-dof-1)

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- NaN
		VKKTemp <- NaN
		KVKTemp <- NaN
		KKVTemp <- NaN	
		VVKTemp <- NaN
		VKVTemp <- NaN
		KVVTemp <- NaN
	
		VVVTemp <- tryCatch(errorVVV_Int(parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = alpha[g,]
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = alphaVariance[g,]
										,clean = TRUE),error = function(e)NaN)

	  
	  c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
	}, BPPARAM=BPPARAM))

	rownames(chi2data) <- eiGenes

 	# P values
 	pvaluesdata <- cbind(KKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKK'], max(c(0,2*length(tpts)-dof[g,'KKK']))))
 						,VKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKK'], max(c(0,2*length(tpts)-dof[g,'VKK']))))
 						,KVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVK'], max(c(0,2*length(tpts)-dof[g,'KVK']))))
 						,KKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKV'], max(c(0,2*length(tpts)-dof[g,'KKV']))))
 						,VVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVK'], max(c(0,2*length(tpts)-dof[g,'VVK']))))
 						,VKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKV'], max(c(0,2*length(tpts)-dof[g,'VKV']))))
 						,KVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVV'], max(c(0,2*length(tpts)-dof[g,'KVV']))))
 						,VVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVV'], max(c(0,2*length(tpts)-dof[g,'VVV'])))))

	ratesSpecs <- lapply(eiGenes,function(g)
	{
		list("0" = list(alpha = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = NaN)
					   ,beta = list(fun = constantModelP
								   ,type = "constant"
								   ,df = 1
								   ,params = NaN)
					   ,gamma = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = NaN)
					   ,test = NaN
					   ,logLik = NaN
					   ,AIC = NaN
					   ,AICc = NaN
					   ,counts = NaN
					   ,convergence = NaN
					   ,message = NaN)
			,"a" = if(length(grep("par",names(VVV[[g]])))==18)
				   {
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,beta = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,beta = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }
			,"b" = if(length(grep("par",names(VVV[[g]])))==18)
				   {
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }else{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = sigmoidModelP
										,type = "sigmoid"
										,df = 4
										,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }
			,"c" = if(length(grep("par",names(VVV[[g]])))==18)
				   {
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }else{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = sigmoidModelP
										,type = "sigmoid"
										,df = 4
										,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }
			,"ab" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}
			,"ac" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}
	 		,"bc" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}else{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}
			,"abc" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(alpha = unname(VVV[[g]][1:6])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(VVV[[g]][13:18])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(VVV[[g]][7:12])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(alpha = unname(VVV[[g]][1:4])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(VVV[[g]][9:12])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(VVV[[g]][5:8])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}
			)
 	})

	names(ratesSpecs) <- eiGenes

 	print("Confidence intervals.")
	confidenceIntervals <- bplapply(eiGenes,function(g)
	{
		classTmp <- "VVV"

		parameters <- VVV[[g]][grep("par",names(VVV[[g]]))]
		optTmp <- rates_integrativeModels(tpts=tpts, class=classTmp, parameters=parameters)

		foe <- capture.output({ # Just to capture the output of multiroot function
			suppressWarnings({
				intervals <- sapply(names(parameters),function(parname)
				{
					par <- parameters[parname]

						mOut = list(
							left_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e-2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
							left_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1/2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
							center = tryCatch(multiroot(f = logLikelihoodCIerror, start = par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
							right_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1.5*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
							right_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = FALSE),error=function(e)return(emptyList))
						)
					precis = sapply(mOut, '[[', 'f.root')

					if( length(which(precis<1e-2))>0 )  {
						conf_int = sapply(mOut[which(precis<1e-2)], '[[', 'root')
						low_int = min(conf_int)
						high_int = max(conf_int)

						left = ifelse( low_int < par, low_int, NA)
						right = ifelse( high_int > par, high_int, NA)

						left = unname(left)
						right = unname(right)

					} else {
						left = NA
						right = NA
					}
					return(c(left,right))
				})
				intervals[1,!is.finite(intervals[2,])] <- NaN
				intervals[2,!is.finite(intervals[1,])] <- NaN
			})
		})

		perturbedRates <- matrix(rep(NaN,3*length(tpts)),ncol=1)
		for(parname in names(parameters))
		{
			for(extremePar in intervals[,parname])
			{
				perturbedParameters <- parameters
				perturbedParameters[parname] <- extremePar

				perturbedRates <- cbind(perturbedRates,rates_derivativeModels(tpts=tpts, class=classTmp, parameters=perturbedParameters))
			}
		};perturbedRates <- perturbedRates[,-1]
		perturbedRates[perturbedRates<0] <- 0

		k1left <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,min,na.rm=TRUE)
		k1TC <- optTmp[grep("alpha",names(optTmp))]
		k1right <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,max,na.rm=TRUE)

		k2left <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,min,na.rm=TRUE)
		k2TC <- optTmp[grep("gamma",names(optTmp))]
		k2right <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,max,na.rm=TRUE)

		k3left <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,min,na.rm=TRUE)
		k3TC <- optTmp[grep("beta",names(optTmp))]
		k3right <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,max,na.rm=TRUE)

		return(list(
			k1 = cbind(left=k1left, opt=k1TC, right=k1right),
			k2 = cbind(left=k2left, opt=k2TC, right=k2right),
			k3 = cbind(left=k3left, opt=k3TC, right=k3right)
			))
	},BPPARAM=BPPARAM)

	names(confidenceIntervals) <- eiGenes

	for(g in seq_along(confidenceIntervals))
	{
		for(r in 1:3)
		{
			confidenceIntervals[[g]][[r]] <- t(apply(confidenceIntervals[[g]][[r]],1,function(row)
			{
				if((!is.finite(row[1])|row[1]==row[2])&(is.finite(row[3])&row[3]!=row[2])) row[1] <- row[2] - (row[3]-row[2])
				if((!is.finite(row[3])|row[3]==row[2])&(is.finite(row[1])&row[1]!=row[2])) row[3] <- row[2] + (row[2]-row[1])
				row
			}))
		}
	}

	k1_low <- median(sapply(confidenceIntervals,function(g){abs(g[[1]][,2] - g[[1]][,1])/g[[1]][,1]}),na.rm=TRUE)
	k1_high <- median(sapply(confidenceIntervals,function(g){abs(g[[1]][,3] - g[[1]][,1])/g[[1]][,1]}),na.rm=TRUE)

	k2_low <- median(sapply(confidenceIntervals,function(g){abs(g[[2]][,2] - g[[2]][,1])/g[[2]][,1]}),na.rm=TRUE)
	k2_high <- median(sapply(confidenceIntervals,function(g){abs(g[[2]][,3] - g[[2]][,1])/g[[2]][,1]}),na.rm=TRUE)

	k3_low <- median(sapply(confidenceIntervals,function(g){abs(g[[3]][,2] - g[[3]][,1])/g[[3]][,1]}),na.rm=TRUE)
	k3_high <- median(sapply(confidenceIntervals,function(g){abs(g[[3]][,3] - g[[3]][,1])/g[[3]][,1]}),na.rm=TRUE)

	median_low <- c(k1=k1_low,k2=k2_low,k3=k3_low)
	median_high <- c(k1=k1_high,k2=k2_high,k3=k3_high)

	for(g in seq_along(confidenceIntervals))
	{
		for(r in 1:3)
		{
			confidenceIntervals[[g]][[r]] <- t(apply(confidenceIntervals[[g]][[r]],1,function(row)
			{
				if(row[1]==row[2] & row[1]==row[3]) row[1] <- row[2]*(1-median_low[[r]]); row[3] <- row[2]*(1+median_high[[r]])
				row
			}))
		}
	}

	# Removal of not modeled genes
	eiGenes <- eiGenes[sapply(confidenceIntervals,function(g)all(is.finite(g[[1]]))&all(is.finite(g[[2]]))&all(is.finite(g[[3]])))]
	names(confidenceIntervals) <- eiGenes
	confidenceIntervals <- confidenceIntervals[eiGenes]

	out <- list(ratesSpecs=ratesSpecs[eiGenes],
				confidenceIntervals=confidenceIntervals)

	return(out)
}

.inspect.engine_Integrative_NoNascent <- function(tpts
												, concentrations
												, rates
												, BPPARAM
												, na.rm
												, verbose
												, testOnSmooth
												, seed
												, nInit
												, nIter
												, limitModelComplexity
												, computeDerivatives = TRUE
												, useSigmoidFun = TRUE
												, initialPenalityRelevance = 1
												, derivativePenalityRelevance = 10^-50
												, llConfidenceThreshold)
{
	total <- concentrations$total
 	totalVariance <- concentrations$total_var

 	premature <- concentrations$preMRNA
 	prematureVariance <- concentrations$preMRNA_var

 	mature <- concentrations$mature
 	matureVariance <- concentrations$mature_var

 	alpha <- rates$alpha
 	alphaVariance <- rates$alpha_var

 	beta <- rates$beta
 	gamma <- rates$gamma

	prematureSmooth <- premature
	matureSmooth <- mature

	eiGenes <- rownames(mature)

 	KKK <- bplapply(eiGenes,function(row){

 			k1Parameters <- mean(alpha[row,])
 			k2Parameters <- mean(gamma[row,])
 			k3Parameters <- mean(beta[row,])

 			unlist(
 				tryCatch(
 					optim(c(k1Parameters, k2Parameters, k3Parameters)
 						,errorKKK_Int
 						,tpts = tpts
 						,premature = premature[row,]
 						,mature = mature[row,]
 						,alpha = NULL
 						,prematureVariance = prematureVariance[row,]
 						,matureVariance = matureVariance[row,]
 						,alphaVariance = NULL
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

	VKK <- bplapply(eiGenes, function(row){
		
		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
 			k1Parameters <- c(rep(KKK[[row]][1],2), max(tpts)/3,1)
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- KKK[[row]][3]

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelVKK <- expressionData_integrativeModels(tpts, class = "VKK", parameters = parameters)
			ratesVKK <- rates_integrativeModels(0, class = "VKK", parameters = parameters)

			prematureEstimated <- modelVKK[grep("^premature",names(modelVKK))]
			matureEstimated <- modelVKK[grep("^mature",names(modelVKK))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesVKK[grep("alpha",names(ratesVKK))]
			betaEstimated <- ratesVKK[grep("beta",names(ratesVKK))]
			gammaEstimated <- ratesVKK[grep("gamma",names(ratesVKK))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorVKK_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
				error=function(e) c(par1 = NaN
								  , par2 = NaN
								  , par3 = NaN
								  , par4 = NaN
								  , par5 = NaN
								  , par6 = NaN
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[1:4],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)})
			k2Parameters <- tryCatch(sigmoidsParameters[2],error=function(e)KKK[[row]][2])
			k3Parameters <- tryCatch(sigmoidsParameters[6],error=function(e)KKK[[row]][3])
 		}else{
			k1Parameters <- c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)
			k2Parameters <- KKK[[row]][2]
 			k3Parameters <- KKK[[row]][3]
 		}

		modelVKK <- expressionData_integrativeModels(tpts, class = "VKK", parameters = parameters)
		ratesVKK <- rates_integrativeModels(0, class = "VKK", parameters = parameters)

		prematureEstimated <- modelVKK[grep("^mature",names(modelVKK))]
		matureEstimated <- modelVKK[grep("^premature",names(modelVKK))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesVKK[grep("alpha",names(ratesVKK))]
		betaEstimated <- ratesVKK[grep("beta",names(ratesVKK))]
		gammaEstimated <- ratesVKK[grep("gamma",names(ratesVKK))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorVKK_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 8)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 6))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
 	names(VKK) <- eiGenes
 	message("Model A finished.")

	KKV <- bplapply(eiGenes, function(row){

		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
			k1Parameters <- KKK[[row]][1]
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- c(rep(KKK[[row]][3],2), max(tpts)/3,1)

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelKKV <- expressionData_integrativeModels(tpts, class = "KKV", parameters = parameters)
			ratesKKV <- rates_integrativeModels(0, class = "KKV", parameters = parameters)

			prematureEstimated <- modelKKV[grep("^premature",names(modelKKV))]
			matureEstimated <- modelKKV[grep("^mature",names(modelKKV))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesKKV[grep("alpha",names(ratesKKV))]
			betaEstimated <- ratesKKV[grep("beta",names(ratesKKV))]
			gammaEstimated <- ratesKKV[grep("gamma",names(ratesKKV))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorKKV_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
				error=function(e) c(par1 = NaN
								  , par2 = NaN
								  , par3 = NaN
								  , par4 = NaN
								  , par5 = NaN
								  , par6 = NaN
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(sigmoidsParameters[1],error=function(e)KKK[[row]][1])
			k2Parameters <- tryCatch(sigmoidsParameters[2],error=function(e)KKK[[row]][2])
			k3Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[3:6],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)})
 		}else{
 			k1Parameters <- KKK[[row]][1]
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)
 		}

		modelKKV <- expressionData_integrativeModels(tpts, class = "KKV", parameters = parameters)
		ratesKKV <- rates_integrativeModels(0, class = "KKV", parameters = parameters)

		prematureEstimated <- modelKKV[grep("^mature",names(modelKKV))]
		matureEstimated <- modelKKV[grep("^premature",names(modelKKV))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesKKV[grep("alpha",names(ratesKKV))]
		betaEstimated <- ratesKKV[grep("beta",names(ratesKKV))]
		gammaEstimated <- ratesKKV[grep("gamma",names(ratesKKV))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKKV_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 8)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 6))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
 	names(KKV) <- eiGenes
 	message("Model B finished.")

	KVK <- bplapply(eiGenes, function(row){

		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
			k1Parameters <- KKK[[row]][1]
			k2Parameters <- c(rep(KKK[[row]][2],2), max(tpts)/3,1)
			k3Parameters <- KKK[[row]][3]

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelKVK <- expressionData_integrativeModels(tpts, class = "KVK", parameters = parameters)
			ratesKVK <- rates_integrativeModels(0, class = "KVK", parameters = parameters)

			prematureEstimated <- modelKVK[grep("^premature",names(modelKVK))]
			matureEstimated <- modelKVK[grep("^mature",names(modelKVK))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesKVK[grep("alpha",names(ratesKVK))]
			betaEstimated <- ratesKVK[grep("beta",names(ratesKVK))]
			gammaEstimated <- ratesKVK[grep("gamma",names(ratesKVK))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorKVK_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
				error=function(e) c(par1 = NaN
								  , par2 = NaN
								  , par3 = NaN
								  , par4 = NaN
								  , par5 = NaN
								  , par6 = NaN
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(sigmoidsParameters[1],error=function(e){KKK[[row]][1]})
			k2Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[2:5],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)})
			k3Parameters <- tryCatch(sigmoidsParameters[6],error=function(e){KKK[[row]][3]})

 		}else{
 			k1Parameters <- KKK[[row]][1]
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)
			k3Parameters <- KKK[[row]][3]
 		}

		modelKVK <- expressionData_integrativeModels(tpts, class = "KVK", parameters = parameters)
		ratesKVK <- rates_integrativeModels(0, class = "KVK", parameters = parameters)

		prematureEstimated <- modelKVK[grep("^mature",names(modelKVK))]
		matureEstimated <- modelKVK[grep("^premature",names(modelKVK))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesKVK[grep("alpha",names(ratesKVK))]
		betaEstimated <- ratesKVK[grep("beta",names(ratesKVK))]
		gammaEstimated <- ratesKVK[grep("gamma",names(ratesKVK))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVK_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 8)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 6))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
 	names(KVK) <- eiGenes
 	message("Model C finished.")

	VKV <- bplapply(eiGenes, function(row){

		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
			k1Parameters <- c(rep(KKK[[row]][1],2), max(tpts)/3,1)
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- c(rep(KKK[[row]][3],2), max(tpts)/3,1)

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelVKV <- expressionData_integrativeModels(tpts, class = "VKV", parameters = parameters)
			ratesVKV <- rates_integrativeModels(0, class = "VKV", parameters = parameters)

			prematureEstimated <- modelVKV[grep("^premature",names(modelVKV))]
			matureEstimated <- modelVKV[grep("^mature",names(modelVKV))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesVKV[grep("alpha",names(ratesVKV))]
			betaEstimated <- ratesVKV[grep("beta",names(ratesVKV))]
			gammaEstimated <- ratesVKV[grep("gamma",names(ratesVKV))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorVKV_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
				error=function(e) c(par1 = NaN
								  , par2 = NaN
								  , par3 = NaN
								  , par4 = NaN
								  , par5 = NaN
								  , par6 = NaN
								  , par7 = NaN
								  , par8 = NaN
								  , par9 = NaN
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[1:4],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)})
			k2Parameters <- tryCatch(sigmoidsParameters[5],error=function(e)KKK[[row]][2])
			k3Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[6:9],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)})

 		}else{
 			k1Parameters <- c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)
 		}

		modelVKV <- expressionData_integrativeModels(tpts, class = "VKV", parameters = parameters)
		ratesVKV <- rates_integrativeModels(0, class = "VKV", parameters = parameters)

		prematureEstimated <- modelVKV[grep("^mature",names(modelVKV))]
		matureEstimated <- modelVKV[grep("^premature",names(modelVKV))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesVKV[grep("alpha",names(ratesVKV))]
		betaEstimated <- ratesVKV[grep("beta",names(ratesVKV))]
		gammaEstimated <- ratesVKV[grep("gamma",names(ratesVKV))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorVKV_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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
							  , convergence = e)
			)
		)

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 13)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 9))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
	names(VKV) <- eiGenes
	message("Model AB finished.")

	VVK <- bplapply(eiGenes, function(row){

		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
			k1Parameters <- c(rep(KKK[[row]][1],2), max(tpts)/3,1)
			k2Parameters <- c(rep(KKK[[row]][2],2), max(tpts)/3,1)
			k3Parameters <- KKK[[row]][3]

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelVVK <- expressionData_integrativeModels(tpts, class = "VVK", parameters = parameters)
			ratesVVK <- rates_integrativeModels(0, class = "VVK", parameters = parameters)

			prematureEstimated <- modelVVK[grep("^premature",names(modelVVK))]
			matureEstimated <- modelVVK[grep("^mature",names(modelVVK))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesVVK[grep("alpha",names(ratesVVK))]
			betaEstimated <- ratesVVK[grep("beta",names(ratesVVK))]
			gammaEstimated <- ratesVVK[grep("gamma",names(ratesVVK))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorVVK_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
				error=function(e) c(par1 = NaN
								  , par2 = NaN
								  , par3 = NaN
								  , par4 = NaN
								  , par5 = NaN
								  , par6 = NaN
								  , par7 = NaN
								  , par8 = NaN
								  , par9 = NaN
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[1:4],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)})
			k2Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[5:8],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)})
			k3Parameters <- tryCatch(sigmoidsParameters[9],error=function(e)KKK[[row]][3])

 		}else{
 			k1Parameters <- c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)
			k3Parameters <- KKK[[row]][3]
 		}

		modelVVK <- expressionData_integrativeModels(tpts, class = "VVK", parameters = parameters)
		ratesVVK <- rates_integrativeModels(0, class = "VVK", parameters = parameters)

		prematureEstimated <- modelVVK[grep("^mature",names(modelVVK))]
		matureEstimated <- modelVVK[grep("^premature",names(modelVVK))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesVVK[grep("alpha",names(ratesVVK))]
		betaEstimated <- ratesVVK[grep("beta",names(ratesVVK))]
		gammaEstimated <- ratesVVK[grep("gamma",names(ratesVVK))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorVVK_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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
							  , convergence = e)
			)
		)

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 13)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 9))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
	names(VVK) <- eiGenes
	message("Model AC finished.")

	KVV <- bplapply(eiGenes, function(row){

		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
			k1Parameters <- KKK[[row]][1]
			k2Parameters <- c(rep(KKK[[row]][2],2), max(tpts)/3,1)
			k3Parameters <- c(rep(KKK[[row]][3],2), max(tpts)/3,1)

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelKVV <- expressionData_integrativeModels(tpts, class = "KVV", parameters = parameters)
			ratesKVV <- rates_integrativeModels(0, class = "KVV", parameters = parameters)

			prematureEstimated <- modelKVV[grep("^premature",names(modelKVV))]
			matureEstimated <- modelKVV[grep("^mature",names(modelKVV))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesKVV[grep("alpha",names(ratesKVV))]
			betaEstimated <- ratesKVV[grep("beta",names(ratesKVV))]
			gammaEstimated <- ratesKVV[grep("gamma",names(ratesKVV))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorKVV_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
				error=function(e) c(par1 = NaN
								  , par2 = NaN
								  , par3 = NaN
								  , par4 = NaN
								  , par5 = NaN
								  , par6 = NaN
								  , par7 = NaN
								  , par8 = NaN
								  , par9 = NaN
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(sigmoidsParameters[1],error=function(e)KKK[[row]][1])
			k2Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[2:5],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)})
			k3Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[6:9],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)})

 		}else{
 			k1Parameters <- KKK[[row]][1]
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)
 		}

		modelKVV <- expressionData_integrativeModels(tpts, class = "KVV", parameters = parameters)
		ratesKVV <- rates_integrativeModels(0, class = "KVV", parameters = parameters)

		prematureEstimated <- modelKVV[grep("^mature",names(modelKVV))]
		matureEstimated <- modelKVV[grep("^premature",names(modelKVV))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesKVV[grep("alpha",names(ratesKVV))]
		betaEstimated <- ratesKVV[grep("beta",names(ratesKVV))]
		gammaEstimated <- ratesKVV[grep("gamma",names(ratesKVV))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVV_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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
							  , convergence = e)
			)
		)

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 13)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 9))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
	names(KVV) <- eiGenes
	message("Model BC finished.")

 	VVV <- bplapply(eiGenes, function(row){

		ratesKKK <- rates_integrativeModels(0, class = "KKK", parameters = KKK[[row]][1:3])

		alphaKKK <- ratesKKK[grep("alpha",names(ratesKKK))]
		betaKKK <- ratesKKK[grep("beta",names(ratesKKK))]
		gammaKKK <- ratesKKK[grep("gamma",names(ratesKKK))]

 		if(useSigmoidFun)
 		{
			k1Parameters <- c(rep(KKK[[row]][1],2), max(tpts)/3,1)
			k2Parameters <- c(rep(KKK[[row]][2],2), max(tpts)/3,1)
			k3Parameters <- c(rep(KKK[[row]][3],2), max(tpts)/3,1)

			parameters <- c(k1Parameters, k2Parameters, k3Parameters)

			modelVVV <- expressionData_integrativeModels(tpts, class = "VVV", parameters = parameters)
			ratesVVV <- rates_integrativeModels(0, class = "VVV", parameters = parameters)

			prematureEstimated <- modelVVV[grep("^premature",names(modelVVV))]
			matureEstimated <- modelVVV[grep("^mature",names(modelVVV))]

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			alphaEstimated <- ratesVVV[grep("alpha",names(ratesVVV))]
			betaEstimated <- ratesVVV[grep("beta",names(ratesVVV))]
			gammaEstimated <- ratesVVV[grep("gamma",names(ratesVVV))]

			sigmoidsParameters <- unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorVVV_Int
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
							 					 , (betaKKK-betaEstimated)^2
							 					 , (gammaKKK-gammaEstimated)^2))
						 ,initialPenalityRelevance = 1
						 ,derivativePenalityRelevance = 10^-50
						 ,clean = FALSE
						 ,control = list(maxit = nIter*10)),
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
								  , value = NaN
								  , counts.function = NaN
								  , counts.gradient = NaN
								  , convergence = e)
				)
			)

			k1Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[1:4],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)})
			k2Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[5:8],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)})
			k3Parameters <- tryCatch(fromSigmoidToImpulse(sigmoidsParameters[9:12],tpts=tpts,nIter=nIter),error=function(e){c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)})

 		}else{
 			k1Parameters <- c(rep(KKK[[row]][1],3), max(tpts)/3, max(tpts)/3*2,1)
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tpts)/3, max(tpts)/3*2,1)
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tpts)/3, max(tpts)/3*2,1)
 		}

		modelVVV <- expressionData_integrativeModels(tpts, class = "VVV", parameters = parameters)
		ratesVVV <- rates_integrativeModels(0, class = "VVV", parameters = parameters)

		prematureEstimated <- modelVVV[grep("^mature",names(modelVVV))]
		matureEstimated <- modelVVV[grep("^premature",names(modelVVV))]

		prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
		matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

		alphaEstimated <- ratesVVV[grep("alpha",names(ratesVVV))]
		betaEstimated <- ratesVVV[grep("beta",names(ratesVVV))]
		gammaEstimated <- ratesVVV[grep("gamma",names(ratesVVV))]

		impulsesParameters <- unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorVVV_Int
					 ,tpts = tpts
					 ,premature = premature[row,]
					 ,mature = mature[row,]
					 ,alpha = NULL
					 ,prematureVariance = prematureVariance[row,]
					 ,matureVariance = matureVariance[row,]
					 ,alphaVariance = NULL
					 ,KKK = KKK[[row]]
					 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
					 ,initialDistances = sum(c((alphaKKK-alphaEstimated)^2
						 					 , (betaKKK-betaEstimated)^2
						 					 , (gammaKKK-gammaEstimated)^2))
					 ,initialPenalityRelevance = 1
					 ,derivativePenalityRelevance = 10^-50
					 ,clean = FALSE
					 ,control = list(maxit = nIter*10)),
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
							  , par14 = NaN
							  , par15 = NaN
							  , par16 = NaN
							  , par17 = NaN
							  , par18 = NaN
							  , value = NaN
							  , counts.function = NaN
							  , counts.gradient = NaN
							  , convergence = e)
			)
		)

		if(!useSigmoidFun)return(impulsesParameters)

		if(!is.finite(sigmoidsParameters["value"])){sigmoidsParameters["value"] <- Inf}
		if(!is.finite(impulsesParameters["value"])){impulsesParameters["value"] <- Inf}

		if(pchisq(impulsesParameters["value"],max(0,2*length(tpts) - 18)) < pchisq(sigmoidsParameters["value"],max(0,2*length(tpts) - 12))){return(impulsesParameters)}else{return(sigmoidsParameters)}

	}, BPPARAM=BPPARAM)
	names(VVV) <- eiGenes
	message("Model ABC finished.")

	# Log likelihood
	logLikelihood <- t(sapply(eiGenes,function(g)
	{
		modelKKK <- expressionData_integrativeModels(tpts, class = "KKK", parameters = KKK[[g]][grep("par",names(KKK[[g]]))])
		
		matureModel <- modelKKK[grep("^mature",names(modelKKK))]
		prematureModel <- modelKKK[grep("^premature",names(modelKKK))]

		modelKKK <- c(matureModel,prematureModel)

		KKKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelKKK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelVKK <- expressionData_integrativeModels(tpts, class = "VKK", parameters = VKK[[g]][grep("par",names(VKK[[g]]))])
		
		matureModel <- modelVKK[grep("^mature",names(modelVKK))]
		prematureModel <- modelVKK[grep("^premature",names(modelVKK))]

		modelVKK <- c(matureModel,prematureModel)

		VKKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVKK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelKVK <- expressionData_integrativeModels(tpts, class = "KVK", parameters = KVK[[g]][grep("par",names(KVK[[g]]))])
		
		matureModel <- modelKVK[grep("^mature",names(modelKVK))]
		prematureModel <- modelKVK[grep("^premature",names(modelKVK))]

		modelKVK <- c(matureModel,prematureModel)

		KVKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelKVK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelKKV <- expressionData_integrativeModels(tpts, class = "KKV", parameters = KKV[[g]][grep("par",names(KKV[[g]]))])
		
		matureModel <- modelKKV[grep("^mature",names(modelKKV))]
		prematureModel <- modelKKV[grep("^premature",names(modelKKV))]

		modelKKV <- c(matureModel,prematureModel)

		KKVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelKKV
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelVVK <- expressionData_integrativeModels(tpts, class = "VVK", parameters = VVK[[g]][grep("par",names(VVK[[g]]))])
		
		matureModel <- modelVVK[grep("^mature",names(modelVVK))]
		prematureModel <- modelVVK[grep("^premature",names(modelVVK))]

		modelVVK <- c(matureModel,prematureModel)

		VVKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVVK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelVKV <- expressionData_integrativeModels(tpts, class = "VKV", parameters = VKV[[g]][grep("par",names(VKV[[g]]))])
		
		matureModel <- modelVKV[grep("^mature",names(modelVKV))]
		prematureModel <- modelVKV[grep("^premature",names(modelVKV))]

		modelVKV <- c(matureModel,prematureModel)

		VKVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVKV
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelKVV <- expressionData_integrativeModels(tpts, class = "KVV", parameters = KVV[[g]][grep("par",names(KVV[[g]]))])
		
		matureModel <- modelKVV[grep("^mature",names(modelKVV))]
		prematureModel <- modelKVV[grep("^premature",names(modelKVV))]

		modelKVV <- c(matureModel,prematureModel)

		KVVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelKVV
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		modelVVV <- expressionData_integrativeModels(tpts, class = "VVV", parameters = VVV[[g]][grep("par",names(VVV[[g]]))])
		
		matureModel <- modelVVV[grep("^mature",names(modelVVV))]
		prematureModel <- modelVVV[grep("^premature",names(modelVVV))]

		modelVVV <- c(matureModel,prematureModel)

		VVVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVVV
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		c("KKK" = KKKTemp,"VKK" = VKKTemp,"KVK" = KVKTemp,"KKV" = KKVTemp,"VVK" = VVKTemp,"VKV" = VKVTemp,"KVV" = KVVTemp,"VVV" = VVVTemp)
	}))

	rownames(logLikelihood) <- eiGenes

	dof <- cbind(KKK = sapply(KKK,function(m)length(grep("par",names(m))))
				,VKK = sapply(VKK,function(m)length(grep("par",names(m))))
				,KVK = sapply(KVK,function(m)length(grep("par",names(m))))
				,KKV = sapply(KKV,function(m)length(grep("par",names(m))))
				,VVK = sapply(VVK,function(m)length(grep("par",names(m))))
				,VKV = sapply(VKV,function(m)length(grep("par",names(m))))
				,KVV = sapply(KVV,function(m)length(grep("par",names(m))))
				,VVV = sapply(VVV,function(m)length(grep("par",names(m)))))

 	AIC <- 2*(dof - logLikelihood)
	AICc <- 2*(dof - logLikelihood) + (2*dof*(dof+1))/max(0,2*length(tpts)-dof-1)

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- tryCatch(errorKKK_Int(parameters = KKK[[g]][grep("par",names(KKK[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL),error = function(e)NaN)

		VKKTemp <- tryCatch(errorVKK_Int(parameters = VKK[[g]][grep("par",names(VKK[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

		KVKTemp <- tryCatch(errorKVK_Int(parameters = KVK[[g]][grep("par",names(KVK[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

		KKVTemp <- tryCatch(errorKKV_Int(parameters = KKV[[g]][grep("par",names(KKV[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

		VVKTemp <- tryCatch(errorVVK_Int(parameters = VVK[[g]][grep("par",names(VVK[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

		VKVTemp <- tryCatch(errorVKV_Int(parameters = VKV[[g]][grep("par",names(VKV[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

		KVVTemp <- tryCatch(errorKVV_Int(parameters = KVV[[g]][grep("par",names(KVV[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

		VVVTemp <- tryCatch(errorVVV_Int(parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
										,tpts = tpts
										,premature = prematureSmooth[g,]
										,mature = matureSmooth[g,]
										,alpha = NULL
										,prematureVariance = prematureVariance[g,]
										,matureVariance = matureVariance[g,]
										,alphaVariance = NULL
										,clean = TRUE),error = function(e)NaN)

	  
	  c(KKK = KKKTemp, VKK = VKKTemp, KVK = KVKTemp, KKV = KKVTemp, VVK = VVKTemp, VKV = VKVTemp, KVV = KVVTemp, VVV = VVVTemp)
	
	}, BPPARAM=BPPARAM))

	rownames(chi2data) <- eiGenes

 	# P values
 	pvaluesdata <- cbind(KKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKK'], max(c(0,2*length(tpts)-dof[g,'KKK']))))
 						,VKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKK'], max(c(0,2*length(tpts)-dof[g,'VKK']))))
 						,KVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVK'], max(c(0,2*length(tpts)-dof[g,'KVK']))))
 						,KKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKV'], max(c(0,2*length(tpts)-dof[g,'KKV']))))
 						,VVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVK'], max(c(0,2*length(tpts)-dof[g,'VVK']))))
 						,VKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKV'], max(c(0,2*length(tpts)-dof[g,'VKV']))))
 						,KVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVV'], max(c(0,2*length(tpts)-dof[g,'KVV']))))
 						,VVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVV'], max(c(0,2*length(tpts)-dof[g,'VVV'])))))

	ratesSpecs <- lapply(eiGenes,function(g)
	{
		list("0" = list(alpha = list(fun = constantModelP
								   ,type = "constant"
								   ,df = 1
								   ,params = c(alpha = unname(KKK[[g]][1])))
					  ,beta = list(fun = constantModelP
								  ,type = "constant"
								  ,df = 1
								  ,params = c(beta = unname(KKK[[g]][3])))
					  ,gamma = list(fun = constantModelP
								   ,type = "constant"
								   ,df = 1
								   ,params = c(gamma = unname(KKK[[g]][2])))
					  ,test = log(pvaluesdata[g,"KKK"])
					  ,logLik = logLikelihood[g,"KKK"]
					  ,AIC = AIC[g,"KKK"]
					  ,AICc = AICc[g,"KKK"]
					  ,counts = c("function"=unname(KKK[[g]]["counts.function"]), gradient=unname(KKK[[g]]["counts.gradient"]))
					  ,convergence = unname(KKK[[g]]["convergence"])
					  ,message = NULL)
			,"a" = if(length(grep("par",names(VKK[[g]])))==8)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(alpha = unname(VKK[[g]][1:6])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(VKK[[g]][8])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKK[[g]][7])))
							,test = log(pvaluesdata[g,"VKK"])
							,logLik = logLikelihood[g,"VKK"]
							,AIC = AIC[g,"VKK"]
							,AICc = AICc[g,"VKK"]
							,counts = c("function"=unname(VKK[[g]]["counts.function"]), gradient=unname(VKK[[g]]["counts.gradient"]))
							,convergence = unname(VKK[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(alpha = unname(VKK[[g]][1:4])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(VKK[[g]][6])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKK[[g]][5])))
							,test = log(pvaluesdata[g,"VKK"])
							,logLik = logLikelihood[g,"VKK"]
							,AIC = AIC[g,"VKK"]
							,AICc = AICc[g,"VKK"]
							,counts = c("function"=unname(VKK[[g]]["counts.function"]), gradient=unname(VKK[[g]]["counts.gradient"]))
							,convergence = unname(VKK[[g]]["convergence"])
							,message = NULL)
					}
			,"b" = if(length(grep("par",names(KKV[[g]])))==8)
					{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KKV[[g]][1])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(KKV[[g]][3:8])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(KKV[[g]][2])))
							,test = log(pvaluesdata[g,"KKV"])
							,logLik = logLikelihood[g,"KKV"]
							,AIC = AIC[g,"KKV"]
							,AICc = AICc[g,"KKV"]
							,counts = c("function"=unname(KKV[[g]]["counts.function"]), gradient=unname(KKV[[g]]["counts.gradient"]))
							,convergence = unname(KKV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KKV[[g]][1])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(KKV[[g]][3:6])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(KKV[[g]][2])))
							,test = log(pvaluesdata[g,"KKV"])
							,logLik = logLikelihood[g,"KKV"]
							,AIC = AIC[g,"KKV"]
							,AICc = AICc[g,"KKV"]
							,counts = c("function"=unname(KKV[[g]]["counts.function"]), gradient=unname(KKV[[g]]["counts.gradient"]))
							,convergence = unname(KKV[[g]]["convergence"])
							,message = NULL)
					}
			,"c" = if(length(grep("par",names(KVK[[g]])))==8)
					{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVK[[g]][1])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(KVK[[g]][8])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(KVK[[g]][2:7])))
							,test = log(pvaluesdata[g,"KVK"])
							,logLik = logLikelihood[g,"KVK"]
							,AIC = AIC[g,"KVK"]
							,AICc = AICc[g,"KVK"]
							,counts = c("function"=unname(KVK[[g]]["counts.function"]), gradient=unname(KVK[[g]]["counts.gradient"]))
							,convergence = unname(KVK[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVK[[g]][1])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(KVK[[g]][6])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(KVK[[g]][2:5])))
							,test = log(pvaluesdata[g,"KVK"])
							,logLik = logLikelihood[g,"KVK"]
							,AIC = AIC[g,"KVK"]
							,AICc = AICc[g,"KVK"]
							,counts = c("function"=unname(KVK[[g]]["counts.function"]), gradient=unname(KVK[[g]]["counts.gradient"]))
							,convergence = unname(KVK[[g]]["convergence"])
							,message = NULL)
					}
			,"ab" = if(length(grep("par",names(VKV[[g]])))==13)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(alpha = unname(VKV[[g]][1:6])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(VKV[[g]][8:13])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKV[[g]][7])))
							,test = log(pvaluesdata[g,"VKV"])
							,logLik = logLikelihood[g,"VKV"]
							,AIC = AIC[g,"VKV"]
							,AICc = AICc[g,"VKV"]
							,counts = c("function"=unname(VKV[[g]]["counts.function"]), gradient=unname(VKV[[g]]["counts.gradient"]))
							,convergence = unname(VKV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(alpha = unname(VKV[[g]][1:4])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(VKV[[g]][6:9])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKV[[g]][5])))
							,test = log(pvaluesdata[g,"VKV"])
							,logLik = logLikelihood[g,"VKV"]
							,AIC = AIC[g,"VKV"]
							,AICc = AICc[g,"VKV"]
							,counts = c("function"=unname(VKV[[g]]["counts.function"]), gradient=unname(VKV[[g]]["counts.gradient"]))
							,convergence = unname(VKV[[g]]["convergence"])
							,message = NULL)
					}
			,"ac" = if(length(grep("par",names(VVK[[g]])))==13)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(alpha = unname(VVK[[g]][1:6])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(VVK[[g]][13])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(VVK[[g]][7:12])))
							,test = log(pvaluesdata[g,"VVK"])
							,logLik = logLikelihood[g,"VVK"]
							,AIC = AIC[g,"VVK"]
							,AICc = AICc[g,"VVK"]
							,counts = c("function"=unname(VVK[[g]]["counts.function"]), gradient=unname(VVK[[g]]["counts.gradient"]))
							,convergence = unname(VVK[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(alpha = unname(VVK[[g]][1:4])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(VVK[[g]][9])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(VVK[[g]][5:8])))
							,test = log(pvaluesdata[g,"VVK"])
							,logLik = logLikelihood[g,"VVK"]
							,AIC = AIC[g,"VVK"]
							,AICc = AICc[g,"VVK"]
							,counts = c("function"=unname(VVK[[g]]["counts.function"]), gradient=unname(VVK[[g]]["counts.gradient"]))
							,convergence = unname(VVK[[g]]["convergence"])
							,message = NULL)
					}
			,"bc" = if(length(grep("par",names(KVV[[g]])))==13)
					{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVV[[g]][1])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(KVV[[g]][8:12])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(KVV[[g]][2:7])))
							,test = log(pvaluesdata[g,"KVV"])
							,logLik = logLikelihood[g,"KVV"]
							,AIC = AIC[g,"KVV"]
							,AICc = AICc[g,"KVV"]
							,counts = c("function"=unname(KVV[[g]]["counts.function"]), gradient=unname(KVV[[g]]["counts.gradient"]))
							,convergence = unname(KVV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVV[[g]][1])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(KVV[[g]][6:9])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(KVV[[g]][2:5])))
							,test = log(pvaluesdata[g,"KVV"])
							,logLik = logLikelihood[g,"KVV"]
							,AIC = AIC[g,"KVV"]
							,AICc = AICc[g,"KVV"]
							,counts = c("function"=unname(KVV[[g]]["counts.function"]), gradient=unname(KVV[[g]]["counts.gradient"]))
							,convergence = unname(KVV[[g]]["convergence"])
							,message = NULL)
					}
			,"abc" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(alpha = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(alpha = unname(VVV[[g]][1:6])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(VVV[[g]][13:18])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(VVV[[g]][7:12])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(alpha = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(alpha = unname(VVV[[g]][1:4])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(VVV[[g]][9:12])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(VVV[[g]][5:8])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}
			)
 	})
	names(ratesSpecs) <- eiGenes
	return(ratesSpecs)
}

#####################################
### Errors integrative functions ####
#####################################

genericError_Int <- function(k1F,k2F,k3F
							,Dk1F,Dk2F,Dk3F
							,tpts
							,premature, mature, alpha
							,prematureVariance, matureVariance, alphaVariance
							,KKK = NULL
							,initialChisquare = NULL
							,initialDistances = NULL
							,initialPenalityRelevance = 1
							,derivativePenalityRelevance = 10^-50
							,clean)
{
	D0_k1 <- Dk1F(0)
	D0_k2 <- Dk2F(0)
	D0_k3 <- Dk3F(0)

	modData <- systemSolution(k1F,k2F,k3F,tpts)

	prematureEstimated <- modData[grep("^premature",names(modData))]
	matureEstimated <- modData[grep("^mature",names(modData))]

	alphaEstimated <- k1F(tpts)
	gammaEstimated <- k2F(tpts)
	betaEstimated <- k3F(tpts)

	modData[modData<0] <- NaN
	alphaEstimated[alphaEstimated<0] <- NaN
	gammaEstimated[gammaEstimated<0] <- NaN
	betaEstimated[betaEstimated<0] <- NaN

	if(any(!is.finite(modData)) | 
	   any(!is.finite(alphaEstimated)) | 
	   any(!is.finite(gammaEstimated)) | 
	   any(!is.finite(betaEstimated)) | 
	   !is.finite(D0_k1) | 
	   !is.finite(D0_k2) | 
	   !is.finite(D0_k3)
	   ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
	{
		alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
		initialPenality <- 0
	}else{
		if(clean)
		{
			initialPenality <- 0
		}else{
			initialPenality <- initialPenalityRelevance*initialChisquare*((k1KKK_Int(0,KKK)-alphaEstimated[1])^2
																		+ (k2KKK_Int(0,KKK)-gammaEstimated[1])^2
																		+ (k3KKK_Int(0,KKK)-betaEstimated[1])^2)
		}
		alphaChiSquare <- 0
	}

	chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
	penalty <- abs(D0_k1)+abs(D0_k2)+abs(D0_k3)

	if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}
	
	if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
}

errorKKK_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance)
{

	if(parameters[1]<0)return(NaN)
	if(parameters[2]<0)return(NaN)
	if(parameters[3]<0)return(NaN)

	k1Parameters <- parameters[1]
	k2Parameters <- parameters[2]
	k3Parameters <- parameters[3]

	prematureEstimated <- rep(k1Parameters/k2Parameters,length(tpts))
	matureEstimated <- rep(k1Parameters/k3Parameters,length(tpts))
	alphaEstimated <- rep(k1Parameters,length(tpts))

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)
	if(is.null(alpha)&is.null(alphaVariance)){alphaChiSquare <- 0}else{alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)}

	return(sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare)))
}

errorVKK_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==8)
	{
		k1F <- function(x) {impulseModel(x,parameters[1:6])}
		Dk1F <- function(x){.DimpulseModel(x,parameters[1:6])}
		k2F <- function(x) {parameters[7]}
		Dk2F <- function(x){0}
		k3F <- function(x) {parameters[8]}
		Dk3F <- function(x){0}
	}else{
		k1F <- function(x) {sigmoidModel(x,parameters[1:4])}
		Dk1F <- function(x){.DsigmoidModel(x,parameters[1:4])}
		k2F <- function(x) {parameters[5]}
		Dk2F <- function(x){0}
		k3F <- function(x) {parameters[6]}
		Dk3F <- function(x){0}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

errorKVK_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==8)
	{
		k1F <- function(x) {parameters[1]}
		Dk1F <- function(x){0}
		k2F <- function(x) {impulseModel(x,parameters[2:7])}
		Dk2F <- function(x){.DimpulseModel(x,parameters[2:7])}
		k3F <- function(x) {parameters[8]}
		Dk3F <- function(x){0}
	}else{
		k1F <- function(x) {parameters[1]}
		Dk1F <- function(x){0}
		k2F <- function(x) {sigmoidModel(x,parameters[2:5])}
		Dk2F <- function(x){.DsigmoidModel(x,parameters[2:5])}
		k3F <- function(x) {parameters[6]}
		Dk3F <- function(x){0}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

errorKKV_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==8)
	{
		k1F <- function(x) {parameters[1]}
		Dk1F <- function(x){0}
		k2F <- function(x) {parameters[2]}
		Dk2F <- function(x){0}
		k3F <- function(x) {impulseModel(x,parameters[3:8])}
		Dk3F <- function(x){.DimpulseModel(x,parameters[3:8])}
	}else{
		k1F <- function(x) {parameters[1]}
		Dk1F <- function(x){0}
		k2F <- function(x) {parameters[2]}
		Dk2F <- function(x){0}
		k3F <- function(x) {sigmoidModel(x,parameters[3:6])}
		Dk3F <- function(x){.DsigmoidModel(x,parameters[3:6])}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

errorVVK_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==13)
	{
		k1F <- function(x) {impulseModel(x,parameters[1:6])}
		Dk1F <- function(x){.DimpulseModel(x,parameters[1:6])}
		k2F <- function(x) {impulseModel(x,parameters[7:12])}
		Dk2F <- function(x){.DimpulseModel(x,parameters[7:12])}
		k3F <- function(x) {parameters[13]}
		Dk3F <- function(x){0}
	}else{
		k1F <- function(x) {sigmoidModel(x,parameters[1:4])}
		Dk1F <- function(x){.DsigmoidModel(x,parameters[1:4])}
		k2F <- function(x) {sigmoidModel(x,parameters[5:8])}
		Dk2F <- function(x){.DsigmoidModel(x,parameters[5:8])}
		k3F <- function(x) {parameters[9]}
		Dk3F <- function(x){0}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

errorVKV_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==13)
	{
		k1F <- function(x) {impulseModel(x,parameters[1:6])}
		Dk1F <- function(x){.DimpulseModel(x,parameters[1:6])}
		k2F <- function(x) {parameters[7]}
		Dk2F <- function(x){0}
		k3F <- function(x) {impulseModel(x,parameters[8:13])}
		Dk3F <- function(x){.DimpulseModel(x,parameters[8:13])}		
	}else{
		k1F <- function(x) {sigmoidModel(x,parameters[1:4])}
		Dk1F <- function(x){.DsigmoidModel(x,parameters[1:4])}
		k2F <- function(x) {parameters[5]}
		Dk2F <- function(x){0}
		k3F <- function(x) {sigmoidModel(x,parameters[6:9])}
		Dk3F <- function(x){.DsigmoidModel(x,parameters[6:9])}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

errorKVV_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==13)
	{
		k1F <- function(x) {parameters[1]}
		Dk1F <- function(x){0}
		k2F <- function(x) {impulseModel(x,parameters[2:7])}
		Dk2F <- function(x){.DimpulseModel(x,parameters[2:7])}
		k3F <- function(x) {impulseModel(x,parameters[8:13])}
		Dk3F <- function(x){.DimpulseModel(x,parameters[8:13])}		
	}else{
		k1F <- function(x) {parameters[1]}
		Dk1F <- function(x){0}
		k2F <- function(x) {sigmoidModel(x,parameters[2:5])}
		Dk2F <- function(x){.DsigmoidModel(x,parameters[2:5])}
		k3F <- function(x) {sigmoidModel(x,parameters[6:9])}
		Dk3F <- function(x){.DsigmoidModel(x,parameters[6:9])}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

errorVVV_Int <- function(parameters, tpts
					   , premature, mature, alpha
					   , prematureVariance, matureVariance, alphaVariance
					   , KKK = NULL
					   , initialChisquare = NULL
					   , initialDistances = NULL
					   , initialPenalityRelevance = 1
					   , derivativePenalityRelevance = 10^-50
					   , clean)
{
	if(length(parameters)==18)
	{
		k1F <- function(x) {impulseModel(x,parameters[1:6])}
		Dk1F <- function(x){.DimpulseModel(x,parameters[1:6])}
		k2F <- function(x) {impulseModel(x,parameters[7:12])}
		Dk2F <- function(x){.DimpulseModel(x,parameters[7:12])}
		k3F <- function(x) {impulseModel(x,parameters[13:18])}
		Dk3F <- function(x){.DimpulseModel(x,parameters[13:18])}		
	}else{
		k1F <- function(x) {sigmoidModel(x,parameters[1:4])}
		Dk1F <- function(x){.DsigmoidModel(x,parameters[1:4])}
		k2F <- function(x) {sigmoidModel(x,parameters[5:8])}
		Dk2F <- function(x){.DsigmoidModel(x,parameters[5:8])}
		k3F <- function(x) {sigmoidModel(x,parameters[9:12])}
		Dk3F <- function(x){.DsigmoidModel(x,parameters[9:12])}
	}
	return(genericError_Int(k1F=k1F,k2F=k2F,k3F=k3F,Dk1F=Dk1F,Dk2F=Dk2F,Dk3F=Dk3F
						   ,tpts=tpts
						   ,premature=premature, mature=mature, alpha=alpha
						   ,prematureVariance=prematureVariance, matureVariance=matureVariance, alphaVariance=alphaVariance
						   ,KKK=KKK
						   ,initialChisquare=initialChisquare
						   ,initialDistances=initialDistances
						   ,initialPenalityRelevance=initialPenalityRelevance
						   ,derivativePenalityRelevance=derivativePenalityRelevance
						   ,clean=clean))
}

fromSigmoidToImpulse <- function(sigmoidsParameters,tpts,nIter)
{
	sigmoidProfile <- sigmoidModel(x=tpts,par=sigmoidsParameters)

	internalError <- function(impulsesParameters,sigmoidsParameters,tpts)
	{
		impulseProfile <- impulseModel(x=tpts,par=impulsesParameters)
		sigmoidProfile <- sigmoidModel(x=tpts,par=sigmoidsParameters)

		if(any(impulseProfile<0))return(NaN)

		chisqFunction(experiment=impulseProfile,model=sigmoidProfile,variance=impulseProfile)
	}

	optTmp <- optim(par=c(head(sigmoidProfile,1),median(sigmoidProfile),tail(sigmoidProfile,1),max(tpts)/3,2*max(tpts)/3,1)
		,internalError
		,sigmoidsParameters=sigmoidsParameters
		,tpts=tpts
		,control = list(maxit = nIter))

	return(optTmp$par)
}

####################################################################
################ NEW INSPEcT DERIVATIVE APPROACHES #################
####################################################################

.inspect.engine_Derivative_Nascent <- function(tpts
											 , concentrations
											 , rates
											 , BPPARAM
											 , na.rm
											 , verbose
											 , testOnSmooth
											 , seed
											 , nInit
											 , nIter
											 , limitModelComplexity
											 , computeDerivatives = TRUE
											 , useSigmoidFun = TRUE
											 , initialPenalityRelevance = 1
											 , derivativePenalityRelevance = 10^-50
											 , llConfidenceThreshold)
{
	total <- concentrations$total
 	totalVariance <- concentrations$total_var

 	premature <- concentrations$preMRNA
 	prematureVariance <- concentrations$preMRNA_var

 	mature <- concentrations$mature
 	matureVariance <- concentrations$mature_var

 	alpha <- rates$alpha
 	alphaVariance <- rates$alpha_var

 	beta <- rates$beta
 	gamma <- rates$gamma

	prematureSmooth <- premature
	matureSmooth <- mature

	eiGenes <- rownames(premature)

 	print("Mature RNA fit.")
	modelMatureRNAfun <- bplapply(eiGenes,function(i)
	{
		tryCatch(.chooseModel(tpts=tpts
				, experiment=mature[i,]
				, variance=matureVariance[i,]
				, na.rm=na.rm
				, sigmoid=useSigmoidFun
				, impulse=TRUE
				, polynomial=FALSE
				, nInit=nInit
				, nIter=nIter
				, sigmoidModel=sigmoidModel
				, impulseModel=impulseModel
				, sigmoidModelP=sigmoidModelP
				, impulseModelP=impulseModelP
				, .polynomialModelP=.polynomialModelP
				, seed = seed
				, computeDerivatives = computeDerivatives
				), error=function(e) return(.emptyGene(e)))
	},BPPARAM=BPPARAM)
	names(modelMatureRNAfun) <- eiGenes

	accelerationCoefficient <- sapply(eiGenes, function(row)
	{
		matureParameters <- unname(modelMatureRNAfun[[row]]$params)

		if(is.null(matureParameters)) return(NaN)

		if(length(matureParameters)==6)
		{
			matureEstimated <- impulseModel(x = tpts, par = matureParameters)

			k2Parameters <- c(rep(mean(gamma[row,]),3), max(tpts)/3, max(tpts)/3*2, 1)
			k3Parameters <- c(rep(mean(beta[row,]),3), max(tpts)/3, max(tpts)/3*2, 1)

			parameters <- c(matureParameters,k2Parameters,k3Parameters)

			D0_M <- .DimpulseModel(0,parameters[1:6])
			D0_k2 <- .DimpulseModel(0,parameters[7:12])
			D0_k3 <- .DimpulseModel(0,parameters[13:18])

		} else {
			matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)

			k2Parameters <- c(rep(mean(gamma[row,]),2), max(tpts)/3, 1)
			k3Parameters <- c(rep(mean(beta[row,]),2), max(tpts)/3, 1)

			parameters <- c(matureParameters,k2Parameters,k3Parameters)

			D0_M <- .DsigmoidModel(0,parameters[1:4])
			D0_k2 <- .DsigmoidModel(0,parameters[5:8])
			D0_k3 <- .DsigmoidModel(0,parameters[9:12])
		}

		D0_P <- .DprematureVVV_Der(0, parameters)

		prematureEstimated <- prematureVVV_Der(x = tpts, parameters = parameters)
		alphaEstimated <- k1VVV_Der(x = tpts, parameters = parameters)

		alphaEstimated[alphaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(all(is.finite(alphaEstimated)) & 
		   all(is.finite(prematureEstimated)) & 
		   all(is.finite(matureEstimated)) & 
		   all(is.finite(c(D0_M, D0_P, D0_k2, D0_k3)))
		) return(1)

		suppressWarnings(optimize( function(x)
		{
			if(length(matureParameters)==6)
			{
				k2Parameters[1:3] <- k2Parameters[1:3]*x
				k3Parameters[1:3] <- k3Parameters[1:3]*x

				parameters <- c(matureParameters,k2Parameters,k3Parameters)

				D0_M <- .DimpulseModel(0,parameters[1:6])
				D0_k2 <- .DimpulseModel(0,parameters[7:12])
				D0_k3 <- .DimpulseModel(0,parameters[13:18])

			}else{
				k2Parameters[1:2] <- k2Parameters[1:2]*x
				k3Parameters[1:2] <- k3Parameters[1:2]*x

				parameters <- c(matureParameters,k2Parameters,k3Parameters)

				D0_M <- .DsigmoidModel(0,parameters[1:4])
				D0_k2 <- .DsigmoidModel(0,parameters[5:8])
				D0_k3 <- .DsigmoidModel(0,parameters[9:12])
			}

			D0_P <- .DprematureVVV_Der(0, parameters)

			prematureEstimated <- prematureVVV_Der(x = tpts, parameters = parameters)
			alphaEstimated <- k1VVV_Der(x = tpts, parameters = parameters)

			alphaEstimated[alphaEstimated<0] <- NaN
			prematureEstimated[prematureEstimated<0] <- NaN
			matureEstimated[matureEstimated<0] <- NaN

			if(all(is.finite(alphaEstimated)) & 
			   all(is.finite(prematureEstimated)) & 
			   all(is.finite(matureEstimated)) & 
			   all(is.finite(c(D0_M, D0_P, D0_k2, D0_k3)))
			) return(x) else NaN
 		},c(1, 1e5) ))$minimum
	})

	eiGenes <- rownames(mature)

 	print("VVV modeling.")
	VVV <- bplapply(eiGenes,function(row){

		matureParameters <- tryCatch(unname(modelMatureRNAfun[[row]]$params),error=function(e) rep(NaN, length(tpts)))
		if(length(matureParameters)==6)
		{
			k2Parameters <- c(rep(mean(gamma[row,]),3)*accelerationCoefficient[row], max(tpts)/3, max(tpts)/3*2, 1)
			k3Parameters <- c(rep(mean(beta[row,]),3)*accelerationCoefficient[row], max(tpts)/3, max(tpts)/3*2, 1)
		}else{
			k2Parameters <- c(rep(mean(gamma[row,]),2)*accelerationCoefficient[row], max(tpts)/3, 1)
			k3Parameters <- c(rep(mean(beta[row,]),2)*accelerationCoefficient[row], max(tpts)/3, 1)
		}

		unlist(
			tryCatch(
	      			optim(unname(c(matureParameters, k2Parameters, k3Parameters))
	                  	 ,errorVVV_Der
        				 ,tpts = tpts
			             ,premature = premature[row,]
			             ,mature = mature[row,]
						 ,alpha = alpha[row,]
			             ,prematureVariance = prematureVariance[row,]
			             ,matureVariance = matureVariance[row,]
			             ,alphaVariance = alphaVariance[row,]
			             ,KKK = NULL
						 ,initialChisquare = NULL
						 ,initialDistances = NULL
						 ,initialPenalityRelevance = initialPenalityRelevance
						 ,derivativePenalityRelevance = derivativePenalityRelevance
						 ,clean = FALSE
	                  	 ,control = list(maxit = nIter * 1000)),
				error=function(e) list("par"=rep(NaN,length(c(matureParameters, k2Parameters, k3Parameters)))
									 , "value" = NaN
									 , "counts" = c("function" = NaN, "gradient" = NaN)
									 , "convergence" = NaN)
	    		)
	    )
	}, BPPARAM=BPPARAM)
	names(VVV) <- eiGenes

	### Standard outputs

	# Log likelihood
	logLikelihood <- t(sapply(eiGenes,function(g)
	{
		prematureVVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureVVV_Der(x = tpts[t], parameters = VVV[[g]][grep("par",names(VVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		matureVVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureVVV_Der(x = tpts[t], parameters = VVV[[g]][grep("par",names(VVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		alphaVVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)k1VVV_Der(x = tpts[t], parameters = VVV[[g]][grep("par",names(VVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		
		modelVVV <- c(matureVVVTemp,prematureVVVTemp,alphaVVVTemp)

		VVVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,],alpha[g,])
		                               , model = modelVVV
		                               , variance = c(matureVariance[g,],prematureVariance[g,],alphaVariance[g,])),error=function(e)NaN)

		c("KKK" = NaN,"VKK" = NaN,"KVK" = NaN,"KKV" = NaN,"VVK" = NaN,"VKV" = NaN,"KVV" = NaN,"VVV" = VVVTemp)
	}))

	rownames(logLikelihood) <- eiGenes

	### Common code for confidence bars computation
	# dof
 		dof <- cbind(KKK = NaN
					,VKK = NaN
					,KVK = NaN
					,KKV = NaN
					,VVK = NaN
					,VKV = NaN
					,KVV = NaN
					,VVV = sapply(VVV,function(m)length(grep("par",names(m)))))

 	AIC <- 2*(dof - logLikelihood)
	AICc <- 2*(dof - logLikelihood) + (2*dof*(dof+1))/max(0,2*length(tpts)-dof-1)

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- NaN
		VKKTemp <- NaN
		KVKTemp <- NaN
		KKVTemp <- NaN	
		VVKTemp <- NaN
		VKVTemp <- NaN
		KVVTemp <- NaN
	
		VVVTemp <- tryCatch(errorVVV_Der(parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = alpha[g,]
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = alphaVariance[g,]
											  , clean = TRUE),error = function(e)NaN)

	  
	  c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
	}, BPPARAM=BPPARAM))

	rownames(chi2data) <- eiGenes

 	# P values
 	pvaluesdata <- cbind(KKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKK'], max(c(0,2*length(tpts)-dof[g,'KKK']))))
 						,VKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKK'], max(c(0,2*length(tpts)-dof[g,'VKK']))))
 						,KVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVK'], max(c(0,2*length(tpts)-dof[g,'KVK']))))
 						,KKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKV'], max(c(0,2*length(tpts)-dof[g,'KKV']))))
 						,VVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVK'], max(c(0,2*length(tpts)-dof[g,'VVK']))))
 						,VKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKV'], max(c(0,2*length(tpts)-dof[g,'VKV']))))
 						,KVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVV'], max(c(0,2*length(tpts)-dof[g,'KVV']))))
 						,VVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVV'], max(c(0,2*length(tpts)-dof[g,'VVV'])))))

	ratesSpecs <- lapply(eiGenes,function(g)
	{
		list("0" = list(mature = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = NaN)
					   ,beta = list(fun = constantModelP
								   ,type = "constant"
								   ,df = 1
								   ,params = NaN)
					   ,gamma = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = NaN)
					   ,test = NaN
					   ,logLik = NaN
					   ,AIC = NaN
					   ,AICc = NaN
					   ,counts = NaN
					   ,convergence = NaN
					   ,message = NaN)
			,"a" = if(length(grep("par",names(VVV[[g]])))==18)
				   {
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,beta = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,beta = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }
			,"b" = if(length(grep("par",names(VVV[[g]])))==18)
				   {
						list(total = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = NaN)
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }else{
						list(total = list(fun = sigmoidModelP
										,type = "sigmoid"
										,df = 4
										,params = NaN)
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }
			,"c" = if(length(grep("par",names(VVV[[g]])))==18)
				   {
						list(total = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = NaN)
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }else{
						list(total = list(fun = sigmoidModelP
										,type = "sigmoid"
										,df = 4
										,params = NaN)
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
				   }
			,"ab" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}
			,"ac" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}
	 		,"bc" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(total = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}else{
						list(total = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = NaN)
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = NaN)
							,test = NaN
							,logLik = NaN
							,AIC = NaN
							,AICc = NaN
							,counts = NaN
							,convergence = NaN
							,message = NaN)
					}
			,"abc" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(mature = unname(VVV[[g]][1:6])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(VVV[[g]][13:18])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(VVV[[g]][7:12])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(mature = unname(VVV[[g]][1:4])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(VVV[[g]][9:12])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(VVV[[g]][5:8])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}
			)
 	})

	names(ratesSpecs) <- eiGenes

 	print("Confidence intervals.")
	confidenceIntervals <- bplapply(eiGenes,function(g)
	{
		classTmp <- "VVV"

		parameters <- VVV[[g]][grep("par",names(VVV[[g]]))]
		optTmp <- rates_derivativeModels(tpts=tpts, class=classTmp, parameters=parameters)

		foe <- capture.output({ # Just to capture the output of multiroot function
			suppressWarnings({
				intervals <- sapply(names(parameters),function(parname)
				{
					par <- parameters[parname]

						mOut = list(
							left_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e-2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
							left_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1/2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
							center = tryCatch(multiroot(f = logLikelihoodCIerror, start = par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
							right_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1.5*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
							right_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = alpha[g,], varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = alphaVariance[g,], confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList))
						)
					precis = sapply(mOut, '[[', 'f.root')

					if( length(which(precis<1e-2))>0 )  {
						conf_int = sapply(mOut[which(precis<1e-2)], '[[', 'root')
						low_int = min(conf_int)
						high_int = max(conf_int)

						left = ifelse( low_int < par, low_int, NA)
						right = ifelse( high_int > par, high_int, NA)

						left = unname(left)
						right = unname(right)

					} else {
						left = NA
						right = NA
					}
					return(c(left,right))
				})
				intervals[1,!is.finite(intervals[2,])] <- NaN
				intervals[2,!is.finite(intervals[1,])] <- NaN
			})
		})

		perturbedRates <- matrix(rep(NaN,3*length(tpts)),ncol=1)
		for(parname in names(parameters))
		{
			for(extremePar in intervals[,parname])
			{
				perturbedParameters <- parameters
				perturbedParameters[parname] <- extremePar

				perturbedRates <- cbind(perturbedRates,rates_derivativeModels(tpts=tpts, class=classTmp, parameters=perturbedParameters))
			}
		};perturbedRates <- perturbedRates[,-1]
		perturbedRates[perturbedRates<0] <- 0

		k1left <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,min,na.rm=TRUE)
		k1TC <- optTmp[grep("alpha",names(optTmp))]
		k1right <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,max,na.rm=TRUE)

		k2left <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,min,na.rm=TRUE)
		k2TC <- optTmp[grep("gamma",names(optTmp))]
		k2right <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,max,na.rm=TRUE)

		k3left <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,min,na.rm=TRUE)
		k3TC <- optTmp[grep("beta",names(optTmp))]
		k3right <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,max,na.rm=TRUE)

		return(list(
			k1 = cbind(left=k1left, opt=k1TC, right=k1right),
			k2 = cbind(left=k2left, opt=k2TC, right=k2right),
			k3 = cbind(left=k3left, opt=k3TC, right=k3right)
			))
	},BPPARAM=BPPARAM)

	names(confidenceIntervals) <- eiGenes

	for(g in seq_along(confidenceIntervals))
	{
		for(r in 1:3)
		{
			confidenceIntervals[[g]][[r]] <- t(apply(confidenceIntervals[[g]][[r]],1,function(row)
			{
				if((!is.finite(row[1])|row[1]==row[2])&(is.finite(row[3])&row[3]!=row[2])) row[1] <- row[2] - (row[3]-row[2])
				if((!is.finite(row[3])|row[3]==row[2])&(is.finite(row[1])&row[1]!=row[2])) row[3] <- row[2] + (row[2]-row[1])
				row
			}))
		}
	}

	k1_low <- median(sapply(confidenceIntervals,function(g){abs(g[[1]][,2] - g[[1]][,1])/g[[1]][,1]}),na.rm=TRUE)
	k1_high <- median(sapply(confidenceIntervals,function(g){abs(g[[1]][,3] - g[[1]][,1])/g[[1]][,1]}),na.rm=TRUE)

	k2_low <- median(sapply(confidenceIntervals,function(g){abs(g[[2]][,2] - g[[2]][,1])/g[[2]][,1]}),na.rm=TRUE)
	k2_high <- median(sapply(confidenceIntervals,function(g){abs(g[[2]][,3] - g[[2]][,1])/g[[2]][,1]}),na.rm=TRUE)

	k3_low <- median(sapply(confidenceIntervals,function(g){abs(g[[3]][,2] - g[[3]][,1])/g[[3]][,1]}),na.rm=TRUE)
	k3_high <- median(sapply(confidenceIntervals,function(g){abs(g[[3]][,3] - g[[3]][,1])/g[[3]][,1]}),na.rm=TRUE)

	median_low <- c(k1=k1_low,k2=k2_low,k3=k3_low)
	median_high <- c(k1=k1_high,k2=k2_high,k3=k3_high)

	for(g in seq_along(confidenceIntervals))
	{
		for(r in 1:3)
		{
			confidenceIntervals[[g]][[r]] <- t(apply(confidenceIntervals[[g]][[r]],1,function(row)
			{
				if(is.finite(row[2]))
				{
					if(row[1]==row[2] & row[1]==row[3]) row[1] <- row[2]*(1-median_low[[r]]); row[3] <- row[2]*(1+median_high[[r]])
				}
				row
			}))
		}
	}

	# Removal of not modeled genes
	eiGenes <- eiGenes[sapply(confidenceIntervals,function(g)all(is.finite(g[[1]]))&all(is.finite(g[[2]]))&all(is.finite(g[[3]])))]
	confidenceIntervals <- confidenceIntervals[eiGenes]

	out <- list(ratesSpecs=ratesSpecs[eiGenes], confidenceIntervals=confidenceIntervals)

	return(out)
}

.inspect.engine_Derivative_NoNascent <- function(tpts
											 , concentrations
											 , rates
											 , BPPARAM
											 , na.rm
											 , verbose
											 , testOnSmooth
											 , seed
											 , nInit
											 , nIter
											 , limitModelComplexity
											 , computeDerivatives = TRUE
											 , useSigmoidFun = TRUE
											 , initialPenalityRelevance = 1
											 , derivativePenalityRelevance = 10^-50
											 , llConfidenceThreshold)
{
	#$# print("newSave")
	total <- concentrations$total
	totalVariance <- concentrations$total_var

	premature <- concentrations$preMRNA
	prematureVariance <- concentrations$preMRNA_var

	mature <- concentrations$mature
	matureVariance <- concentrations$mature_var

	alpha <- rates$alpha
	beta <- rates$beta
	gamma <- rates$gamma

	prematureSmooth <- premature
	matureSmooth <- mature

	eiGenes <- rownames(total)

 	message(".inspect.engine_Derivative_NoNascent: experimental data functional form acquisition:")
 	message("          mature RNA")
	modelMatureRNAfun <- bplapply(eiGenes,function(i)
	{
		tryCatch(.chooseModel(tpts=tpts
				, experiment=mature[i,]
				, variance=matureVariance[i,]
				, na.rm=na.rm
				, sigmoid=useSigmoidFun
				, impulse=TRUE 
				, polynomial=FALSE
				, nInit=nInit
				, nIter=nIter
				, sigmoidModel=sigmoidModel
				, impulseModel=impulseModel
				, sigmoidModelP=sigmoidModelP
				, impulseModelP=impulseModelP
				, .polynomialModelP=.polynomialModelP
				, seed = seed
				, computeDerivatives = computeDerivatives
				), error=function(e) return(.emptyGene(e)))
	},BPPARAM=BPPARAM)
	names(modelMatureRNAfun) <- eiGenes
	#$# saveRDS(modelMatureRNAfun,"modelMatureRNAfun.rds")

 	message("          total RNA")
	modelTotalRNAfun <- bplapply(eiGenes,function(i)
	{
		tryCatch(.chooseModel(tpts=tpts
				, experiment=total[i,]
				, variance=totalVariance[i,]
				, na.rm=na.rm
				, sigmoid=useSigmoidFun
				, impulse=TRUE
				, polynomial=FALSE
				, nInit=nInit
				, nIter=nIter
				, sigmoidModel=sigmoidModel
				, impulseModel=impulseModel
				, sigmoidModelP=sigmoidModelP
				, impulseModelP=impulseModelP
				, .polynomialModelP=.polynomialModelP
				, seed = seed
				, computeDerivatives = computeDerivatives
				), error=function(e) return(.emptyGene(e)))
	},BPPARAM=BPPARAM)
	names(modelTotalRNAfun) <- eiGenes
	#$# saveRDS(modelTotalRNAfun,"modelTotalRNAfun.rds")

	impulseGenes <- as.numeric(sapply(modelMatureRNAfun,"[[","type")=="impulse") + 
					as.numeric(sapply(modelTotalRNAfun,"[[","type")=="impulse")

	impulseGenes <- eiGenes[impulseGenes!=2&impulseGenes!=0]

	message("          mixed genes coercion: impulse")
	for(i in impulseGenes)
	{
		modelMatureRNAfun[[i]] <- tryCatch(.chooseModel(tpts=tpts
										 , experiment=mature[i,]
										 , variance=matureVariance[i,]
										 , na.rm=na.rm
										 , sigmoid=FALSE
										 , impulse=TRUE
										 , polynomial=FALSE
										 , nInit=nInit
										 , nIter=nIter
										 , sigmoidModel=sigmoidModel
										 , impulseModel=impulseModel
										 , sigmoidModelP=sigmoidModelP
										 , impulseModelP=impulseModelP
										 , .polynomialModelP=.polynomialModelP
										 , seed = seed
										 , computeDerivatives = computeDerivatives
										 ), error=function(e) return(.emptyGene(e)))

		modelTotalRNAfun[[i]] <- tryCatch(.chooseModel(tpts=tpts
										, experiment=total[i,]
										, variance=totalVariance[i,]
										, na.rm=na.rm
										, sigmoid=FALSE
										, impulse=TRUE
										, polynomial=FALSE
										, nInit=nInit
										, nIter=nIter
										, sigmoidModel=sigmoidModel
										, impulseModel=impulseModel
										, sigmoidModelP=sigmoidModelP
										, impulseModelP=impulseModelP
										, .polynomialModelP=.polynomialModelP
										, seed = seed
										, computeDerivatives = computeDerivatives
										), error=function(e) return(.emptyGene(e)))

	}

	message("          mixed genes coercion: sigmoid")
	sigmoidGenes <- eiGenes[!sapply(sapply(modelMatureRNAfun,"[[","message"),is.null)|!sapply(sapply(modelTotalRNAfun,"[[","message"),is.null)]
	for(i in sigmoidGenes)
	{
		modelMatureRNAfun[[i]] <- tryCatch(.chooseModel(tpts=tpts
										 , experiment=mature[i,]
										 , variance=matureVariance[i,]
										 , na.rm=na.rm
										 , sigmoid=TRUE
										 , impulse=FALSE
										 , polynomial=FALSE
										 , nInit=nInit
										 , nIter=nIter
										 , sigmoidModel=sigmoidModel
										 , impulseModel=impulseModel
										 , sigmoidModelP=sigmoidModelP
										 , impulseModelP=impulseModelP
										 , .polynomialModelP=.polynomialModelP
										 , seed = seed
										 , computeDerivatives = computeDerivatives
										 ), error=function(e) return(.emptyGene(e)))

		modelTotalRNAfun[[i]] <- tryCatch(.chooseModel(tpts=tpts
										, experiment=total[i,]
										, variance=totalVariance[i,]
										, na.rm=na.rm
										, sigmoid=TRUE
										, impulse=FALSE
										, polynomial=FALSE
										, nInit=nInit
										, nIter=nIter
										, sigmoidModel=sigmoidModel
										, impulseModel=impulseModel
										, sigmoidModelP=sigmoidModelP
										, impulseModelP=impulseModelP
										, .polynomialModelP=.polynomialModelP
										, seed = seed
										, computeDerivatives = computeDerivatives
										), error=function(e) return(.emptyGene(e)))
	}
	#$# saveRDS(modelMatureRNAfun,"modelMatureRNAfun.rds")
	#$# saveRDS(modelTotalRNAfun,"modelTotalRNAfun.rds")

	### KKK
		KKK <- bplapply(eiGenes,function(row)
		{
			matureParameters <- tryCatch(mean(modelMatureRNAfun[[row]]$fun$value(tpts, modelMatureRNAfun[[row]]$params)),error=function(e)list(NaN))
			k2Parameters <- mean(gamma,na.rm=TRUE)
			k3Parameters <- mean(beta,na.rm=TRUE)

			unlist(
				tryCatch(
					optim(c(matureParameters, k2Parameters, k3Parameters)
						,errorKKK_Der
						,tpts = tpts
						,premature = premature[row,]
						,mature = mature[row,]
						,alpha = NULL
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,alphaVariance = NULL
						,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=c(NaN,NaN,NaN)
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
				)
			)
		}, BPPARAM=BPPARAM)
		names(KKK) <- eiGenes
		message("Model 0 finished.")
	#$# saveRDS(KKK,"KKK.rds")

	### Acceleration coefficients
	accelerationCoefficient <- sapply(eiGenes, function(row)
	{
		matureParameters <- unname(modelMatureRNAfun[[row]]$params)

		if(is.null(matureParameters)) return(NaN)

		k2Parameters <- KKK[[row]]["par2"] 
		k3Parameters <- KKK[[row]]["par3"] 

		parameters <- c(matureParameters,k2Parameters,k3Parameters)

		prematureEstimated <- prematureVKK_Der(x = tpts, parameters = parameters)
		alphaEstimated <- k1VKK_Der(x = tpts, parameters = parameters)

		prematureEstimated[prematureEstimated<0] <- NaN
		alphaEstimated[alphaEstimated<0] <- NaN

		if(all(is.finite(alphaEstimated)) & all(is.finite(prematureEstimated))) return(1)

		suppressWarnings(optimize( function(x)
		{
			k2Parameters <- k2Parameters*x
			k2Tmp <- k2Parameters*length(tpts)

			k3Parameters <- k3Parameters*x
			k3Tmp <- k3Parameters*length(tpts)

			parameters <- c(matureParameters,k2Parameters,k3Parameters)	
			
			prematureEstimated <- prematureVKK_Der(x = tpts, parameters = parameters)
			alphaEstimated <- k1VKK_Der(x = tpts, parameters = parameters)

			prematureEstimated[prematureEstimated<0] <- NaN
			alphaEstimated[alphaEstimated<0] <- NaN

			if(any(!is.finite(prematureEstimated))|any(!is.finite(alphaEstimated))) NaN else x
  		},c(1, 1e5) ))$minimum
	})
	#$# saveRDS(accelerationCoefficient,"accelerationCoefficient.rds")

	accelerationCoefficient_constantSynthesis_variableDegradation <- sapply(eiGenes, function(row)
	{
		totalParameters <- modelTotalRNAfun[[row]]$params

		if(is.null(totalParameters)) return(NaN)

		k1Parameters <- k1KKK_Der(0,KKK[[row]][1:3]) 
		k3Parameters <- KKK[[row]]["par3"] 

		parameters <- c(totalParameters,k1Parameters,k3Parameters)

		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))
			alphaEstimated <- rep(parameters[7], length(tpts))
			betaEstimated <- rep(parameters[8], length(tpts))
		} else {
			totalParameters <- parameters[1:4]
			totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))
			alphaEstimated <- rep(parameters[5], length(tpts))
			betaEstimated <- rep(parameters[6], length(tpts))
		}

		prematureEstimated <- sapply(tpts,function(t)prematureKVK_Der(x = t, parameters = parameters))
		matureEstimated <- totalEstimated - prematureEstimated

		gammaEstimated <- sapply(tpts,function(t)k2KVK_Der(t, parameters))

		totalEstimated[totalEstimated<0] <- NaN
		alphaEstimated[alphaEstimated<0] <- NaN
		betaEstimated[betaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(all(is.finite(c(totalEstimated,alphaEstimated,betaEstimated,prematureEstimated,matureEstimated,gammaEstimated)))) return(1)

		suppressWarnings(optimize(function(x)
		{
			parameters <- c(totalParameters,k1Parameters,k3Parameters*x)
			
			if(length(parameters)==8)
			{
				totalParameters <- parameters[1:6]
				totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))
				alphaEstimated <- rep(parameters[7], length(tpts))
				betaEstimated <- rep(parameters[8], length(tpts))
			} else {
				totalParameters <- parameters[1:4]
				totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))
				alphaEstimated <- rep(parameters[5], length(tpts))
				betaEstimated <- rep(parameters[6], length(tpts))
			}
	
			prematureEstimated <- sapply(tpts,function(t)prematureKVK_Der(x = t, parameters = parameters))
			matureEstimated <- totalEstimated - prematureEstimated
	
			gammaEstimated <- sapply(tpts,function(t)k2KVK_Der(t, parameters))

			totalEstimated[totalEstimated<0] <- NaN
			alphaEstimated[alphaEstimated<0] <- NaN
			betaEstimated[betaEstimated<0] <- NaN
			prematureEstimated[prematureEstimated<0] <- NaN
			matureEstimated[matureEstimated<0] <- NaN


			if(all(is.finite(c(totalEstimated
							  ,alphaEstimated
							  ,betaEstimated
							  ,prematureEstimated
							  ,matureEstimated)))) x else NaN
 		},c(1,1e5)))$minimum
	})
	#$# saveRDS(accelerationCoefficient_constantSynthesis_variableDegradation,"accelerationCoefficient_constantSynthesis_variableDegradation.rds")

	accelerationCoefficient_constantSynthesis_variableProcessing <- sapply(eiGenes, function(row)
	{
		totalParameters <- modelTotalRNAfun[[row]]$params

		if(is.null(totalParameters)) return(NaN)

		k1Parameters <- k1KKK_Der(0,KKK[[row]][1:3]) 
		k2Parameters <- KKK[[row]]["par2"] 

		parameters <- c(totalParameters,k1Parameters,k2Parameters)

		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))
			alphaEstimated <- rep(parameters[7], length(tpts))
			gammaEstimated <- rep(parameters[8], length(tpts))
		} else {
			totalParameters <- parameters[1:4]
			totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))
			alphaEstimated <- rep(parameters[5], length(tpts))
			gammaEstimated <- rep(parameters[6], length(tpts))
		}

		prematureEstimated <- sapply(tpts,function(t)prematureKKV_Der(x = t, parameters = parameters))
		matureEstimated <- totalEstimated - prematureEstimated

		betaEstimated <- sapply(tpts,function(t)k3KKV_Der(t, parameters))

		totalEstimated[totalEstimated<0] <- NaN
		alphaEstimated[alphaEstimated<0] <- NaN
		betaEstimated[betaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(all(is.finite(c(totalEstimated,alphaEstimated,betaEstimated,prematureEstimated,matureEstimated,gammaEstimated)))) return(1)

		suppressWarnings(optimize( function(x)
		{
			parameters <- c(totalParameters,k1Parameters,k2Parameters*x)
			
			if(length(parameters)==8)
			{
				totalParameters <- parameters[1:6]
				totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))
				alphaEstimated <- rep(parameters[7], length(tpts))
				gammaEstimated <- rep(parameters[8], length(tpts))
			} else {
				totalParameters <- parameters[1:4]
				totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))
				alphaEstimated <- rep(parameters[5], length(tpts))
				gammaEstimated <- rep(parameters[6], length(tpts))
			}
	
			prematureEstimated <- sapply(tpts,function(t)prematureKKV_Der(x = t, parameters = parameters))
			matureEstimated <- totalEstimated - prematureEstimated
	
			betaEstimated <- sapply(tpts,function(t)k3KKV_Der(t, parameters))

			totalEstimated[totalEstimated<0] <- NaN
			alphaEstimated[alphaEstimated<0] <- NaN
			betaEstimated[betaEstimated<0] <- NaN
			prematureEstimated[prematureEstimated<0] <- NaN
			matureEstimated[matureEstimated<0] <- NaN


			if(all(is.finite(c(totalEstimated
							  ,alphaEstimated
							  ,betaEstimated
							  ,prematureEstimated
							  ,matureEstimated)))) x else NaN
 		},c(1,1e5)))$minimum
	})
	#$# saveRDS(accelerationCoefficient_constantSynthesis_variableProcessing,"accelerationCoefficient_constantSynthesis_variableProcessing.rds")

	# VKK
		VKK <- bplapply(eiGenes,function(row){

			matureParameters <- tryCatch(unname(modelMatureRNAfun[[row]]$params),error=function(e) rep(NaN, length(tpts)))
			k2Parameters <- accelerationCoefficient[row]*KKK[[row]]["par2"] 
			k3Parameters <- accelerationCoefficient[row]*KKK[[row]]["par3"] 

			parameters <- unname(c(matureParameters,k2Parameters,k3Parameters))

			prematureEstimated <- prematureVKK_Der(x = tpts, parameters = parameters)
			matureEstimated <- matureVKK_Der(x = tpts, parameters = parameters)
			alphaEstimated <- k1VKK_Der(x = tpts, parameters = parameters)

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])
		  
			unlist(
				tryCatch(
					optim(parameters
						 ,errorVKK_Der
						 ,tpts = tpts
						 ,premature = premature[row,]
						 ,mature = mature[row,]
						 ,alpha = NULL
						 ,prematureVariance = prematureVariance[row,]
						 ,matureVariance = matureVariance[row,]
						 ,alphaVariance = NULL
						 ,KKK = KKK[[row]]
						 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
						 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1VKK_Der(0,parameters))^2
							 					  ,(k2KKK_Der(0,KKK[[row]])-k2VKK_Der(0,parameters))^2
							 					  ,(k3KKK_Der(0,KKK[[row]])-k3VKK_Der(0,parameters))^2))
						 ,initialPenalityRelevance = initialPenalityRelevance
						 ,derivativePenalityRelevance = derivativePenalityRelevance
						 ,clean = FALSE
						 ,control = list(maxit = nIter * 1000)),
				error=function(e) list("par"=rep(NaN,length(c(matureParameters, k2Parameters, k3Parameters)))
									  ,"value" = NaN
									  ,"counts" = c("function" = NaN, "gradient" = NaN)
									  ,"convergence" = NaN)
				)
		    )
		}, BPPARAM=BPPARAM)
		names(VKK) <- eiGenes
		message("Model A finished.")
	#$# saveRDS(VKK,"VKK.rds")

	# KKV
		KKV <- bplapply(eiGenes,function(row){

			totalParameters <- modelTotalRNAfun[[row]]$params
			k1Parameters <- k1KKK_Der(0,KKK[[row]][1:3])
			k2Parameters <- accelerationCoefficient_constantSynthesis_variableProcessing[row]*KKK[[row]]["par2"] 

			parameters <- unname(c(totalParameters, k1Parameters, k2Parameters))

			if(length(parameters)==8)
			{
				totalParameters <- parameters[1:6]
				totalEstimated <- sapply(tpts, function(t)impulseModel(t, totalParameters))
			}else{
				totalParameters <- parameters[1:4]
				totalEstimated <- sapply(tpts, function(t)sigmoidModel(t, totalParameters))
			}

			prematureEstimated <- sapply(tpts, function(x)prematureKKV_Der(x, parameters))
			matureEstimated <- totalEstimated - prematureEstimated

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])
		  
			unlist(
				tryCatch(
		     			optim(parameters
							 ,errorKKV_Der
							 ,tpts = tpts
				             ,premature = premature[row,]
				             ,mature = mature[row,]
							 ,alpha = NULL
				             ,prematureVariance = prematureVariance[row,]
				             ,matureVariance = matureVariance[row,]
				             ,alphaVariance = NULL
							 ,KKK = KKK[[row]]
							 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
							 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1KKV_Der(0,parameters))^2
							 						 , (k2KKK_Der(0,KKK[[row]])-k2KKV_Der(0,parameters))^2
							 						 , (k3KKK_Der(0,KKK[[row]])-k3KKV_Der(0,parameters))^2))
							 ,initialPenalityRelevance = initialPenalityRelevance
							 ,derivativePenalityRelevance = derivativePenalityRelevance
							 ,clean = FALSE
		                  	 ,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=rep(NaN,length(c(totalParameters, k1Parameters, k2Parameters)))
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
		   		)
		    )
		}, BPPARAM=BPPARAM)
		names(KKV) <- eiGenes
		message("Model B finished.")
	#$# saveRDS(KKV,"KKV.rds")

	# KVK
		KVK <- bplapply(eiGenes,function(row){

			totalParameters <- modelTotalRNAfun[[row]]$params
			k1Parameters <- k1KKK_Der(0,KKK[[row]][1:3]) 
			k3Parameters <- accelerationCoefficient_constantSynthesis_variableDegradation[row]*KKK[[row]]["par3"]

			parameters <- unname(c(totalParameters, k1Parameters, k3Parameters))

			if(length(parameters)==8)
			{
				totalParameters <- parameters[1:6]
				totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))
			} else {
				totalParameters <- parameters[1:4]
				totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))
			}

			prematureEstimated <- sapply(tpts,function(t)prematureKVK_Der(x = t, parameters = parameters))
			matureEstimated <- totalEstimated - prematureEstimated

			prematureChiSquare <- sum((premature[row,] - prematureEstimated)^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])
		  
			unlist(
				tryCatch(
		     			optim(parameters
							 ,errorKVK_Der
							 ,tpts = tpts
				             ,premature = premature[row,]
				             ,mature = mature[row,]
							 ,alpha = NULL
				             ,prematureVariance = prematureVariance[row,]
				             ,matureVariance = matureVariance[row,]
				             ,alphaVariance = NULL
							 ,KKK = KKK[[row]]
							 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
							 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1KVK_Der(0,parameters))^2
							 						 , (k2KKK_Der(0,KKK[[row]])-k2KVK_Der(0,parameters))^2
							 						 , (k3KKK_Der(0,KKK[[row]])-k3KVK_Der(0,parameters))^2))
							 ,initialPenalityRelevance = initialPenalityRelevance
							 ,derivativePenalityRelevance = derivativePenalityRelevance
							 ,clean = FALSE
		                  	 ,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=rep(NaN,length(c(totalParameters, k1Parameters, k3Parameters)))
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
		   		)
		    )
		}, BPPARAM=BPPARAM)
		names(KVK) <- eiGenes
		message("Model C finished.")
	#$# saveRDS(KVK,"KVK.rds")

	# VKV
		VKV <- bplapply(eiGenes,function(row){

			matureParameters <- tryCatch(unname(modelMatureRNAfun[[row]]$params),error=function(e) rep(NaN, length(tpts)))
			k2Parameters <- accelerationCoefficient[row]*KKK[[row]]["par2"] 
			
			if(length(matureParameters)==6)
			{
				k3Parameters <- c(rep(KKK[[row]]["par3"],3)*accelerationCoefficient[row], max(tpts)/3, max(tpts)/3*2, 1)
			}else{
				k3Parameters <- c(rep(KKK[[row]]["par3"],2)*accelerationCoefficient[row], max(tpts)/3, 1)
			}

			parameters <- unname(c(matureParameters, k2Parameters, k3Parameters))

			if(length(parameters)==13)
			{
				matureParameters <- parameters[1:6]
				matureEstimated <- impulseModel(x = tpts, par = matureParameters)
			} else {
				matureParameters <- parameters[1:4]
				matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)
			}

			prematureEstimated <- prematureVKV_Der(x = tpts, parameters = parameters)

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			unlist(
				tryCatch(
		      			optim(parameters
		                  	 ,errorVKV_Der
							 ,tpts = tpts
				             ,premature = premature[row,]
				             ,mature = mature[row,]
							 ,alpha = NULL
				             ,prematureVariance = prematureVariance[row,]
				             ,matureVariance = matureVariance[row,]
				             ,alphaVariance = NULL
							 ,KKK = KKK[[row]]
							 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
							 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1VKV_Der(0,parameters))^2
							 						 , (k2KKK_Der(0,KKK[[row]])-k2VKV_Der(0,parameters))^2
							 						 , (k3KKK_Der(0,KKK[[row]])-k3VKV_Der(0,parameters))^2))
							 ,initialPenalityRelevance = initialPenalityRelevance
							 ,derivativePenalityRelevance = derivativePenalityRelevance
							 ,clean = FALSE
		                  	 ,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=rep(NaN,length(c(matureParameters, k2Parameters, k3Parameters)))
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
		    		)
		    )
		}, BPPARAM=BPPARAM)
		names(VKV) <- eiGenes
		message("Model AB finished.")
	#$# saveRDS(VKV,"VKV.rds")

	# VVK
		VVK <- bplapply(eiGenes,function(row){

			matureParameters <- tryCatch(unname(modelMatureRNAfun[[row]]$params),error=function(e) rep(NaN, length(tpts)))
			
			if(length(matureParameters)==6)
			{
				k2Parameters <- c(rep(KKK[[row]]["par2"],3)*accelerationCoefficient[row], max(tpts)/3, max(tpts)/3*2, 1)
			}else{
				k2Parameters <- c(rep(KKK[[row]]["par2"],2)*accelerationCoefficient[row], max(tpts)/3, 1)
			}

			k3Parameters <- accelerationCoefficient[row]*KKK[[row]]["par3"] 

			parameters <- unname(c(matureParameters, k2Parameters, k3Parameters))

			if(length(parameters)==13)
			{
				matureParameters <- parameters[1:6]
				matureEstimated <- impulseModel(x = tpts, par = matureParameters)
			} else {
				matureParameters <- parameters[1:4]
				matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)
			}

			prematureEstimated <- prematureVVK_Der(x = tpts, parameters = parameters)

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			unlist(
				tryCatch(
		      			optim(parameters
		                  	 ,errorVVK_Der
							 ,tpts = tpts
				             ,premature = premature[row,]
				             ,mature = mature[row,]
							 ,alpha = NULL
				             ,prematureVariance = prematureVariance[row,]
				             ,matureVariance = matureVariance[row,]
				             ,alphaVariance = NULL
							 ,KKK = KKK[[row]]
							 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
							 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1VVK_Der(0,parameters))^2
							 						 , (k2KKK_Der(0,KKK[[row]])-k2VVK_Der(0,parameters))^2
							 						 , (k3KKK_Der(0,KKK[[row]])-k3VVK_Der(0,parameters))^2))
							 ,initialPenalityRelevance = initialPenalityRelevance
							 ,derivativePenalityRelevance = derivativePenalityRelevance
							 ,clean = FALSE
		                  	 ,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=rep(NaN,length(c(matureParameters, k2Parameters, k3Parameters)))
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
		    		)
		    )
		}, BPPARAM=BPPARAM)
		names(VVK) <- eiGenes
		message("Model AC finished.")
	#$# saveRDS(VVK,"VVK.rds")

	# KVV
		KVV <- bplapply(eiGenes,function(row){

			totalParameters <- modelTotalRNAfun[[row]]$params
			k1Parameters <- k1KKK_Der(0,KKK[[row]][1:3]) 
			if(length(totalParameters)==6)
			{
				k3Parameters <- c(rep(KKK[[row]]["par3"],3)*accelerationCoefficient_constantSynthesis_variableDegradation[row], max(tpts)/3, max(tpts)/3*2, 1)
			}else{
				k3Parameters <- c(rep(KKK[[row]]["par3"],2)*accelerationCoefficient_constantSynthesis_variableDegradation[row], max(tpts)/3, 1)
			}

			parameters <- unname(c(totalParameters, k1Parameters, k3Parameters))

			prematureEstimated <- sapply(tpts,function(t)prematureKVV_Der(x = t, parameters = parameters))
			matureEstimated <- matureKVV_Der(tpts,parameters)

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])
  
			unlist(
				tryCatch(
		     			optim(parameters
		                  	 ,errorKVV_Der
							 ,tpts = tpts
				             ,premature = premature[row,]
				             ,mature = mature[row,]
							 ,alpha = NULL
				             ,prematureVariance = prematureVariance[row,]
				             ,matureVariance = matureVariance[row,]
				             ,alphaVariance = NULL
							 ,KKK = KKK[[row]]
							 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
							 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1KVV_Der(0,parameters))^2
							 						 , (k2KKK_Der(0,KKK[[row]])-k2KVV_Der(0,parameters))^2
							 						 , (k3KKK_Der(0,KKK[[row]])-k3KVV_Der(0,parameters))^2))
							 ,initialPenalityRelevance = initialPenalityRelevance
							 ,derivativePenalityRelevance = derivativePenalityRelevance
							 ,clean = FALSE
		                  	 ,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=rep(NaN,length(c(totalParameters, k1Parameters, k3Parameters)))
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
		   		)
		    )
		}, BPPARAM=BPPARAM)
		names(KVV) <- eiGenes
		message("Model BC finished.")
	#$# saveRDS(KVV,"KVV.rds")

	# VVV
		VVV <- bplapply(eiGenes,function(row){

			matureParameters <- tryCatch(unname(modelMatureRNAfun[[row]]$params),error=function(e) rep(NaN, length(tpts)))
			if(length(matureParameters)==6)
			{
				k2Parameters <- c(rep(KKK[[row]]["par2"],3)*accelerationCoefficient[row], max(tpts)/3, max(tpts)/3*2, 1)
				k3Parameters <- c(rep(KKK[[row]]["par3"],3)*accelerationCoefficient[row], max(tpts)/3, max(tpts)/3*2, 1)
			}else{
				k2Parameters <- c(rep(KKK[[row]]["par2"],2)*accelerationCoefficient[row], max(tpts)/3, 1)
				k3Parameters <- c(rep(KKK[[row]]["par3"],2)*accelerationCoefficient[row], max(tpts)/3, 1)
			}

			parameters <- unname(c(matureParameters, k2Parameters, k3Parameters))

			if(length(parameters)==18)
			{
				matureParameters <- parameters[1:6]
				matureEstimated <- impulseModel(x = tpts, par = matureParameters)
			} else {
				matureParameters <- parameters[1:4]
				matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)
			}

			prematureEstimated <- prematureVVV_Der(x = tpts, parameters = parameters)

			prematureChiSquare <- sum((premature[row,] - prematureEstimated )^2/prematureVariance[row,])
			matureChiSquare <- sum((mature[row,] - matureEstimated)^2/matureVariance[row,])

			unlist(
				tryCatch(
		      			optim(parameters
		                  	 ,errorVVV_Der
							 ,tpts = tpts
				             ,premature = premature[row,]
				             ,mature = mature[row,]
							 ,alpha = NULL
				             ,prematureVariance = prematureVariance[row,]
				             ,matureVariance = matureVariance[row,]
				             ,alphaVariance = NULL
							 ,KKK = KKK[[row]]
							 ,initialChisquare = sum(c(prematureChiSquare,matureChiSquare))
							 ,initialDistances = sum(c((k1KKK_Der(0,KKK[[row]])-k1VVV_Der(0,parameters))^2
							 						 , (k2KKK_Der(0,KKK[[row]])-k2VVV_Der(0,parameters))^2
							 						 , (k3KKK_Der(0,KKK[[row]])-k3VVV_Der(0,parameters))^2))
							 ,initialPenalityRelevance = initialPenalityRelevance
							 ,derivativePenalityRelevance = derivativePenalityRelevance
							 ,clean = FALSE
		                  	 ,control = list(maxit = nIter * 1000)),
					error=function(e) list("par"=rep(NaN,length(c(matureParameters, k2Parameters, k3Parameters)))
										 , "value" = NaN
										 , "counts" = c("function" = NaN, "gradient" = NaN)
										 , "convergence" = NaN)
		    		)
		    )
		}, BPPARAM=BPPARAM)
		names(VVV) <- eiGenes
		message("Model ABC finished.")
	#$# saveRDS(VVV,"VVV.rds")

	# Log likelihood
	logLikelihood <- t(sapply(eiGenes,function(g)
	{
		prematureKKKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureKKK_Der(x = tpts[t], parameters = KKK[[g]][grep("par",names(KKK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
			prematureKVKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureKVK_Der(x = tpts[t], parameters = KVK[[g]][grep("par",names(KVK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
			prematureKKVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureKKV_Der(x = tpts[t], parameters = KKV[[g]][grep("par",names(KKV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		prematureVKKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureVKK_Der(x = tpts[t], parameters = VKK[[g]][grep("par",names(VKK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		prematureVVKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureVVK_Der(x = tpts[t], parameters = VVK[[g]][grep("par",names(VVK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		prematureVKVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureVKV_Der(x = tpts[t], parameters = VKV[[g]][grep("par",names(VKV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
			prematureKVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureKVV_Der(x = tpts[t], parameters = KVV[[g]][grep("par",names(KVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		prematureVVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)prematureVVV_Der(x = tpts[t], parameters = VVV[[g]][grep("par",names(VVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))

		matureKKKTemp <- tryCatch(rep(KKK[[g]][[1]],length(tpts)),error=function(e)rep(NaN,length(tpts)))
		matureVKKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureVKK_Der(x = tpts[t], parameters = VKK[[g]][grep("par",names(VKK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
			matureKVKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureKVK_Der(x = tpts[t], parameters = KVK[[g]][grep("par",names(KVK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
			matureKKVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureKKV_Der(x = tpts[t], parameters = KKV[[g]][grep("par",names(KKV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		matureVVKTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureVVK_Der(x = tpts[t], parameters = VVK[[g]][grep("par",names(VVK[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		matureVKVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureVKV_Der(x = tpts[t], parameters = VKV[[g]][grep("par",names(VKV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
			matureKVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureKVV_Der(x = tpts[t], parameters = KVV[[g]][grep("par",names(KVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))
		matureVVVTemp <- tryCatch(c(sapply(seq_along(tpts),function(t)matureVVV_Der(x = tpts[t], parameters = VVV[[g]][grep("par",names(VVV[[g]]))]))),error=function(e)rep(NaN,length(tpts)))

		modelKKK <- c(matureKKKTemp,prematureKKKTemp)
			modelKVK <- c(matureKVKTemp,prematureKVKTemp)
			modelKKV <- c(matureKKVTemp,prematureKKVTemp)
		modelVKK <- c(matureVKKTemp,prematureVKKTemp)
		modelVVK <- c(matureVVKTemp,prematureVVKTemp)
		modelVKV <- c(matureVKVTemp,prematureVKVTemp)
			modelKVV <- c(matureKVVTemp,prematureKVVTemp)
		modelVVV <- c(matureVVVTemp,prematureVVVTemp)

		KKKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelKKK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
			KVKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
			                               , model = modelKVK
			                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
			KKVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
			                               , model = modelKKV
			                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
		VKKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVKK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
		VVKTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVVK
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
		VKVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVKV
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
			KVVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
			                               , model = modelKVV
			                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)
		VVVTemp <- tryCatch(logLikelihoodFunction(experiment = c(matureSmooth[g,],prematureSmooth[g,])
		                               , model = modelVVV
		                               , variance = c(matureVariance[g,],prematureVariance[g,])),error=function(e)NaN)

		c("KKK" = KKKTemp,"VKK" = VKKTemp,"KVK" = KVKTemp,"KKV" = KKVTemp,"VVK" = VVKTemp,"VKV" = VKVTemp,"KVV" = KVVTemp,"VVV" = VVVTemp)
	}))

	rownames(logLikelihood) <- eiGenes
	#$# saveRDS(logLikelihood,"logLikelihood.rds")

	### Common code for confidence bars computation
	# dof
 		dof <- cbind(KKK = sapply(KKK,function(m)length(grep("par",names(m))))
					,VKK = sapply(VKK,function(m)length(grep("par",names(m))))
					,KVK = sapply(KVK,function(m)length(grep("par",names(m))))
					,KKV = sapply(KKV,function(m)length(grep("par",names(m))))
					,VVK = sapply(VVK,function(m)length(grep("par",names(m))))
					,VKV = sapply(VKV,function(m)length(grep("par",names(m))))
					,KVV = sapply(KVV,function(m)length(grep("par",names(m))))
					,VVV = sapply(VVV,function(m)length(grep("par",names(m)))))

 	AIC <- 2*(dof - logLikelihood)
	AICc <- 2*(dof - logLikelihood) + (2*dof*(dof+1))/max(0,2*length(tpts)-dof-1)

	#$# saveRDS(AIC,"AIC.rds")
	#$# saveRDS(AICc,"AICc.rds")

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- tryCatch(errorKKK_Der(parameters = KKK[[g]][grep("par",names(KKK[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL),error = function(e)NaN)

		VKKTemp <- tryCatch(errorVKK_Der(parameters = VKK[[g]][grep("par",names(VKK[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)

		KVKTemp <- tryCatch(errorKVK_Der(parameters = KVK[[g]][grep("par",names(KVK[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)

		KKVTemp <- tryCatch(errorKKV_Der(parameters = KKV[[g]][grep("par",names(KKV[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)
	
		VVKTemp <- tryCatch(errorVVK_Der(parameters = VVK[[g]][grep("par",names(VVK[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)

		VKVTemp <- tryCatch(errorVKV_Der(parameters = VKV[[g]][grep("par",names(VKV[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)

		KVVTemp <- tryCatch(errorKVV_Der(parameters = KVV[[g]][grep("par",names(KVV[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)
	
		VVVTemp <- tryCatch(errorVVV_Der(parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
											  , tpts = tpts
											  , premature = prematureSmooth[g,]
											  , mature = matureSmooth[g,]
											  , alpha = NULL
											  , prematureVariance = prematureVariance[g,]
											  , matureVariance = matureVariance[g,]
											  , alphaVariance = NULL
											  , clean = TRUE),error = function(e)NaN)

	  c(KKK = KKKTemp,VKK = VKKTemp,KVK = KVKTemp,KKV = KKVTemp,VVK = VVKTemp,VKV = VKVTemp,KVV = KVVTemp,VVV = VVVTemp)
	
	}, BPPARAM=BPPARAM))

	rownames(chi2data) <- eiGenes
	#$# saveRDS(chi2data,"chi2data.rds")

 	# P values
 	pvaluesdata <- cbind(KKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKK'], max(c(0,2*length(tpts)-dof[g,'KKK']))))
 						,VKK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKK'], max(c(0,2*length(tpts)-dof[g,'VKK']))))
 						,KVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVK'], max(c(0,2*length(tpts)-dof[g,'KVK']))))
 						,KKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KKV'], max(c(0,2*length(tpts)-dof[g,'KKV']))))
 						,VVK=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVK'], max(c(0,2*length(tpts)-dof[g,'VVK']))))
 						,VKV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VKV'], max(c(0,2*length(tpts)-dof[g,'VKV']))))
 						,KVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'KVV'], max(c(0,2*length(tpts)-dof[g,'KVV']))))
 						,VVV=sapply(eiGenes,function(g)pchisq(chi2data[g,'VVV'], max(c(0,2*length(tpts)-dof[g,'VVV'])))))	
	#$# saveRDS(pvaluesdata,"pvaluesdata.rds")

	ratesSpecs <- lapply(eiGenes,function(g)
	{
		list("0" = list(mature = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(mature = unname(KKK[[g]]["par1"])))
					   ,beta = list(fun = constantModelP
								   ,type = "constant"
								   ,df = 1
								   ,params = c(beta = unname(KKK[[g]]["par3"])))
					   ,gamma = list(fun = constantModelP
									,type = "constant"
									,df = 1
									,params = c(gamma = unname(KKK[[g]]["par2"])))
					   ,test = log(pvaluesdata[g,"KKK"])
					   ,logLik = logLikelihood[g,"KKK"]
					   ,AIC = AIC[g,"KKK"]
					   ,AICc = AICc[g,"KKK"]
					   ,counts = c("function"=unname(KKK[[g]]["counts.function"]), gradient=unname(KKK[[g]]["counts.gradient"]))
					   ,convergence = unname(KKK[[g]]["convergence"])
					   ,message = NULL)
			,"a" = if(length(grep("par",names(VKK[[g]])))==8)
				   {
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(mature = unname(VKK[[g]][1:6])))
							,beta = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(beta = unname(VKK[[g]][8])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKK[[g]][7])))
							,test = log(pvaluesdata[g,"VKK"])
							,logLik = logLikelihood[g,"VKK"]
							,AIC = AIC[g,"VKK"]
							,AICc = AICc[g,"VKK"]
							,counts = c("function"=unname(VKK[[g]]["counts.function"]), gradient=unname(VKK[[g]]["counts.gradient"]))
							,convergence = unname(VKK[[g]]["convergence"])
							,message = NULL)
				   }else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(mature = unname(VKK[[g]][1:4])))
							,beta = list(fun = constantModelP
										,type = "constant"
										,df = 1
										,params = c(beta = unname(VKK[[g]][6])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKK[[g]][5])))
							,test = log(pvaluesdata[g,"VKK"])
							,logLik = logLikelihood[g,"VKK"]
							,AIC = AIC[g,"VKK"]
							,AICc = AICc[g,"VKK"]
							,counts = c("function"=unname(VKK[[g]]["counts.function"]), gradient=unname(VKK[[g]]["counts.gradient"]))
							,convergence = unname(VKK[[g]]["convergence"])
							,message = NULL)
				   }
			,"b" = if(length(grep("par",names(KKV[[g]])))==8)
				   {
						list(total = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(total = unname(KKV[[g]][1:6])))
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KKV[[g]][7])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(KKV[[g]][8])))
							,test = log(pvaluesdata[g,"KKV"])
							,logLik = logLikelihood[g,"KKV"]
							,AIC = AIC[g,"KKV"]
							,AICc = AICc[g,"KKV"]
							,counts = c("function"=unname(KKV[[g]]["counts.function"]), gradient=unname(KKV[[g]]["counts.gradient"]))
							,convergence = unname(KKV[[g]]["convergence"])
							,message = NULL)
				   }else{
						list(total = list(fun = sigmoidModelP
										,type = "sigmoid"
										,df = 4
										,params = c(total = unname(KKV[[g]][1:4])))
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KKV[[g]][5])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(KKV[[g]][6])))
							,test = log(pvaluesdata[g,"KKV"])
							,logLik = logLikelihood[g,"KKV"]
							,AIC = AIC[g,"KKV"]
							,AICc = AICc[g,"KKV"]
							,counts = c("function"=unname(KKV[[g]]["counts.function"]), gradient=unname(KKV[[g]]["counts.gradient"]))
							,convergence = unname(KKV[[g]]["convergence"])
							,message = NULL)
				   }
			,"c" = if(length(grep("par",names(KVK[[g]])))==8)
				   {
						list(total = list(fun = impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(total = unname(KVK[[g]][1:6])))
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVK[[g]][7])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(KVK[[g]][8])))
							,test = log(pvaluesdata[g,"KVK"])
							,logLik = logLikelihood[g,"KVK"]
							,AIC = AIC[g,"KVK"]
							,AICc = AICc[g,"KVK"]
							,counts = c("function"=unname(KVK[[g]]["counts.function"]), gradient=unname(KVK[[g]]["counts.gradient"]))
							,convergence = unname(KVK[[g]]["convergence"])
							,message = NULL)
				   }else{
						list(total = list(fun = sigmoidModelP
										,type = "sigmoid"
										,df = 4
										,params = c(total = unname(KVK[[g]][1:4])))
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVK[[g]][5])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(KVK[[g]][6])))
							,test = log(pvaluesdata[g,"KVK"])
							,logLik = logLikelihood[g,"KVK"]
							,AIC = AIC[g,"KVK"]
							,AICc = AICc[g,"KVK"]
							,counts = c("function"=unname(KVK[[g]]["counts.function"]), gradient=unname(KVK[[g]]["counts.gradient"]))
							,convergence = unname(KVK[[g]]["convergence"])
							,message = NULL)
				   }
			,"ab" = if(length(grep("par",names(VKV[[g]])))==13)
					{
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(mature = unname(VKV[[g]][1:6])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(VKV[[g]][8:13])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKV[[g]][7])))
							,test = log(pvaluesdata[g,"VKV"])
							,logLik = logLikelihood[g,"VKV"]
							,AIC = AIC[g,"VKV"]
							,AICc = AICc[g,"VKV"]
							,counts = c("function"=unname(VKV[[g]]["counts.function"]), gradient=unname(VKV[[g]]["counts.gradient"]))
							,convergence = unname(VKV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(mature = unname(VKV[[g]][1:4])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(VKV[[g]][6:9])))
							,gamma = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(gamma = unname(VKV[[g]][5])))
							,test = log(pvaluesdata[g,"VKV"])
							,logLik = logLikelihood[g,"VKV"]
							,AIC = AIC[g,"VKV"]
							,AICc = AICc[g,"VKV"]
							,counts = c("function"=unname(VKV[[g]]["counts.function"]), gradient=unname(VKV[[g]]["counts.gradient"]))
							,convergence = unname(VKV[[g]]["convergence"])
							,message = NULL)
					}
			,"ac" = if(length(grep("par",names(VVK[[g]])))==13)
					{
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(mature = unname(VVK[[g]][1:6])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(VVK[[g]][13])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(VVK[[g]][7:12])))
							,test = log(pvaluesdata[g,"VVK"])
							,logLik = logLikelihood[g,"VVK"]
							,AIC = AIC[g,"VVK"]
							,AICc = AICc[g,"VVK"]
							,counts = c("function"=unname(VVK[[g]]["counts.function"]), gradient=unname(VVK[[g]]["counts.gradient"]))
							,convergence = unname(VVK[[g]]["convergence"])
							,message = NULL)
					}else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(mature = unname(VVK[[g]][1:4])))
							,beta = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(beta = unname(VVK[[g]][9])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(VVK[[g]][5:8])))
							,test = log(pvaluesdata[g,"VVK"])
							,logLik = logLikelihood[g,"VVK"]
							,AIC = AIC[g,"VVK"]
							,AICc = AICc[g,"VVK"]
							,counts = c("function"=unname(VVK[[g]]["counts.function"]), gradient=unname(VVK[[g]]["counts.gradient"]))
							,convergence = unname(VVK[[g]]["convergence"])
							,message = NULL)
					}
	 		,"bc" = if(length(grep("par",names(KVV[[g]])))==13)
					{
						list(total = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(total = unname(KVV[[g]][1:6])))
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVV[[g]][7])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(KVV[[g]][8:13])))
							,test = log(pvaluesdata[g,"KVV"])
							,logLik = logLikelihood[g,"KVV"]
							,AIC = AIC[g,"KVV"]
							,AICc = AICc[g,"KVV"]
							,counts = c("function"=unname(KVV[[g]]["counts.function"]), gradient=unname(KVV[[g]]["counts.gradient"]))
							,convergence = unname(KVV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(total = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(total = unname(KVV[[g]][1:4])))
							,alpha = list(fun = constantModelP
										 ,type = "constant"
										 ,df = 1
										 ,params = c(alpha = unname(KVV[[g]][5])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(KVV[[g]][6:9])))
							,test = log(pvaluesdata[g,"KVV"])
							,logLik = logLikelihood[g,"KVV"]
							,AIC = AIC[g,"KVV"]
							,AICc = AICc[g,"KVV"]
							,counts = c("function"=unname(KVV[[g]]["counts.function"]), gradient=unname(KVV[[g]]["counts.gradient"]))
							,convergence = unname(KVV[[g]]["convergence"])
							,message = NULL)
					}
			,"abc" = if(length(grep("par",names(VVV[[g]])))==18)
					{
						list(mature = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(mature = unname(VVV[[g]][1:6])))
							,beta = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(beta = unname(VVV[[g]][13:18])))
							,gamma = list(fun = impulseModelP
										 ,type = "impulse"
										 ,df = 6
										 ,params = c(gamma = unname(VVV[[g]][7:12])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}else{
						list(mature = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(mature = unname(VVV[[g]][1:4])))
							,beta = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(beta = unname(VVV[[g]][9:12])))
							,gamma = list(fun = sigmoidModelP
										 ,type = "sigmoid"
										 ,df = 4
										 ,params = c(gamma = unname(VVV[[g]][5:8])))
							,test = log(pvaluesdata[g,"VVV"])
							,logLik = logLikelihood[g,"VVV"]
							,AIC = AIC[g,"VVV"]
							,AICc = AICc[g,"VVV"]
							,counts = c("function"=unname(VVV[[g]]["counts.function"]), gradient=unname(VVV[[g]]["counts.gradient"]))
							,convergence = unname(VVV[[g]]["convergence"])
							,message = NULL)
					}
			)
 	})
	names(ratesSpecs) <- eiGenes
	#$# saveRDS(ratesSpecs,"ratesSpecs.rds")

	return(ratesSpecs)
}

####################################
### Errors derivative functions ####
####################################

	errorKKK_Der <- function(parameters, tpts
						   , premature, mature, alpha
						   , prematureVariance, matureVariance, alphaVariance)
	{
	
		if(parameters[1]<0)return(NaN)
		if(parameters[2]<0)return(NaN)
		if(parameters[3]<0)return(NaN)
	
		matureParameters <- parameters[1]
	
		prematureEstimated <- prematureKKK_Der(x = tpts, parameters = parameters)
		matureEstimated <- rep(matureParameters,length(tpts))
		alphaEstimated <- k1KKK_Der(x = tpts, parameters = parameters)
	
		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)
		if(is.null(alpha)&is.null(alphaVariance)){alphaChiSquare <- 0}else{alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)}
	
		return(sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare)))
	}

	errorVKK_Der <- function(parameters, tpts
									   , premature, mature, alpha
									   , prematureVariance, matureVariance, alphaVariance
									   , KKK = NULL
									   , initialChisquare = NULL
									   , initialDistances = NULL
									   , initialPenalityRelevance = 1
									   , derivativePenalityRelevance = 1
									   , clean
									   )
	{
		if(length(parameters)==8)
		{
			matureParameters <- parameters[1:6]
			matureEstimated <- impulseModel(x = tpts, par = matureParameters)

			D0_M <- .DimpulseModel(0,parameters[1:6])
			D0_k2 <- 0
			D0_k3 <- 0

		} else {
			matureParameters <- parameters[1:4]
			matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)

			D0_M <- .DsigmoidModel(0,parameters[1:4])
			D0_k2 <- 0
			D0_k3 <- 0
		}

		D0_P <- .DprematureVKK_Der(0, parameters)

		prematureEstimated <- prematureVKK_Der(x = tpts, parameters = parameters)
		alphaEstimated <- k1VKK_Der(x = tpts, parameters = parameters)

		alphaEstimated[alphaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(any(!is.finite(alphaEstimated)) | 
		   any(!is.finite(prematureEstimated)) | 
		   any(!is.finite(matureEstimated)) | 
		   !is.finite(D0_M) | 
		   !is.finite(D0_k2) | 
		   !is.finite(D0_k3) | 
		   !is.finite(D0_P)
		   ) return(NaN)

		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

		if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
		{
			alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
			initialPenality <- 0
		}else{
			if(clean){initialPenality <- 0}else{
			initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1VKK_Der(0,parameters))^2
																						   + (k2KKK_Der(0,KKK)-k2VKK_Der(0,parameters))^2
																						   + (k3KKK_Der(0,KKK)-k3VKK_Der(0,parameters))^2)
			}
			alphaChiSquare <- 0
		}

		chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
		penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)

		if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}

		if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
	}

	errorKVK_Der <- function(parameters, tpts
						   , premature, mature, alpha
						   , prematureVariance, matureVariance, alphaVariance
						   , KKK = NULL
						   , initialChisquare = NULL
						   , initialDistances = NULL
						   , initialPenalityRelevance = 1
						   , derivativePenalityRelevance = 1
						   , clean
						   )
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k3Parameters <- parameters[8]

			totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))

			D0_T <- .DimpulseModel(0,parameters[1:6])
		} else {
			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k3Parameters <- parameters[6]

			totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))

			D0_T <- .DsigmoidModel(0,parameters[1:4])
		}

		D0_k1 <- 0
		D0_k3 <- 0

		alphaEstimated <- rep(k1Parameters, length(tpts))
		betaEstimated <- rep(k3Parameters, length(tpts))

		gammaEstimated <- sapply(tpts,function(t)k2KVK_Der(t, parameters))
		prematureEstimated <- sapply(tpts,function(t)prematureKVK_Der(x = t, parameters = parameters))
		matureEstimated <- totalEstimated - prematureEstimated

		D0_P <- alphaEstimated[[1]] - gammaEstimated[[1]]*prematureEstimated[[1]]

		totalEstimated[totalEstimated<0] <- NaN
		alphaEstimated[alphaEstimated<0] <- NaN
		betaEstimated[betaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN
		gammaEstimated[gammaEstimated<0] <- NaN

		if(any(!is.finite(alphaEstimated)) | 
		   any(!is.finite(betaEstimated)) | 
		   any(!is.finite(gammaEstimated)) | 
		   any(!is.finite(prematureEstimated)) | 
		   any(!is.finite(matureEstimated)) | 
		   any(!is.finite(totalEstimated)) | 
		   !is.finite(D0_T) | 
		   !is.finite(D0_k1) | 
		   !is.finite(D0_k3) | 
		   !is.finite(D0_P)) return(NaN)

		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

		if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
		{
			alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
			initialPenality <- 0
		}else{
			if(clean){initialPenality <- 0}else{
		initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1KVK_Der(0,parameters))^2
																					   + (k2KKK_Der(0,KKK)-k2KVK_Der(0,parameters))^2
																					   + (k3KKK_Der(0,KKK)-k3KVK_Der(0,parameters))^2)
			}
			alphaChiSquare <- 0
		}

		chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
		penalty <- (abs(D0_T)+abs(D0_P)+abs(D0_k1)+abs(D0_k3))

		if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}

		if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}

	}

	errorKKV_Der <- function(parameters, tpts
						   , premature, mature, alpha
						   , prematureVariance, matureVariance, alphaVariance
						   , KKK = NULL
						   , initialChisquare = NULL
						   , initialDistances = NULL
						   , initialPenalityRelevance = 1
						   , derivativePenalityRelevance = 1
						   , clean
						   )
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k2Parameters <- parameters[8]

			totalEstimated <- sapply(tpts, function(t)impulseModel(t, totalParameters))

			D0_T <- .DimpulseModel(0,parameters[1:6])

		}else{

			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k2Parameters <- parameters[6]

			totalEstimated <- sapply(tpts, function(t)sigmoidModel(t, totalParameters))
			
			D0_T <- .DsigmoidModel(0,parameters[1:4])
		}

		D0_k1 <- 0
		D0_k2 <- 0

		alphaEstimated <- rep(k1Parameters, length(tpts))
		gammaEstimated <- rep(k2Parameters, length(tpts))

		betaEstimated <- sapply(tpts, function(x)k3KKV_Der(x, parameters))
		prematureEstimated <- sapply(tpts, function(x)prematureKKV_Der(x, parameters))
		matureEstimated <- totalEstimated - prematureEstimated

		D0_P <- 0

		totalEstimated[totalEstimated<0] <- NaN
		alphaEstimated[alphaEstimated<0] <- NaN
		gammaEstimated[gammaEstimated<0] <- NaN
		betaEstimated[betaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(any(!is.finite(alphaEstimated)) | 
		   any(!is.finite(betaEstimated)) | 
		   any(!is.finite(gammaEstimated)) | 
		   any(!is.finite(prematureEstimated)) | 
		   any(!is.finite(matureEstimated)) | 
		   !is.finite(D0_T) | 
		   !is.finite(D0_k1) | 
		   !is.finite(D0_k2) | 
		   !is.finite(D0_P)
		   ) return(NaN)

		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

		if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
		{
			alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
			initialPenality <- 0
		}else{
			if(clean){initialPenality <- 0}else{
		initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1KKV_Der(0,parameters))^2
																					 + (k2KKK_Der(0,KKK)-k2KKV_Der(0,parameters))^2
																					 + (k3KKK_Der(0,KKK)-k3KKV_Der(0,parameters))^2)
			}
			alphaChiSquare <- 0
		}

		chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
		penalty <- abs(D0_T)+abs(D0_P)+abs(D0_k1)+abs(D0_k2)

		if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}

		if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
	}

	errorVVK_Der <- function(parameters, tpts
						   , premature, mature, alpha
						   , prematureVariance, matureVariance, alphaVariance
						   , KKK = NULL
						   , initialChisquare = NULL
						   , initialDistances = NULL
						   , initialPenalityRelevance = 1
						   , derivativePenalityRelevance = 1
						   , clean
						   )
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			matureEstimated <- impulseModel(x = tpts, par = matureParameters)

			D0_M <- .DimpulseModel(0,parameters[1:6])
			D0_k2 <- .DimpulseModel(0,parameters[7:12])
			D0_k3 <- 0

		} else {
			matureParameters <- parameters[1:4]
			matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)

			D0_M <- .DsigmoidModel(0,parameters[1:4])
			D0_k2 <- .DsigmoidModel(0,parameters[5:8])
			D0_k3 <- 0
		}

		D0_P <- .DprematureVVK_Der(0, parameters)

		prematureEstimated <- prematureVVK_Der(x = tpts, parameters = parameters)
		alphaEstimated <- k1VVK_Der(x = tpts, parameters = parameters)

		alphaEstimated[alphaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(any(!is.finite(alphaEstimated)) | 
		   any(!is.finite(prematureEstimated)) | 
		   any(!is.finite(matureEstimated)) | 
		   !is.finite(D0_M) | 
		   !is.finite(D0_k2) | 
		   !is.finite(D0_k3) | 
		   !is.finite(D0_P)
		   ) return(NaN)

		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

		if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
		{
			alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
			initialPenality <- 0
		}else{
			if(clean){initialPenality <- 0}else{
		initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1VVK_Der(0,parameters))^2
																					 + (k2KKK_Der(0,KKK)-k2VVK_Der(0,parameters))^2
																					 + (k3KKK_Der(0,KKK)-k3VVK_Der(0,parameters))^2)
			}
			alphaChiSquare <- 0
		}

		chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
		penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)

		if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}

		if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
	}

	errorVKV_Der <- function(parameters, tpts
						   , premature, mature, alpha
						   , prematureVariance, matureVariance, alphaVariance
						   , KKK = NULL
						   , initialChisquare = NULL
						   , initialDistances = NULL
						   , initialPenalityRelevance = 1
						   , derivativePenalityRelevance = 1
						   , clean
						   )
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			matureEstimated <- impulseModel(x = tpts, par = matureParameters)

			D0_M <- .DimpulseModel(0,parameters[1:6])
			D0_k2 <- 0
			D0_k3 <- .DimpulseModel(0,parameters[8:13])

		} else {
			matureParameters <- parameters[1:4]
			matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)

			D0_M <- .DsigmoidModel(0,parameters[1:4])
			D0_k2 <- 0
			D0_k3 <- .DsigmoidModel(0,parameters[6:9])
		}

		D0_P <- .DprematureVKV_Der(0, parameters)

		prematureEstimated <- prematureVKV_Der(x = tpts, parameters = parameters)
		alphaEstimated <- k1VKV_Der(x = tpts, parameters = parameters)

		alphaEstimated[alphaEstimated<0] <- NaN
		prematureEstimated[prematureEstimated<0] <- NaN
		matureEstimated[matureEstimated<0] <- NaN

		if(any(!is.finite(alphaEstimated)) | 
		   any(!is.finite(prematureEstimated)) | 
		   any(!is.finite(matureEstimated)) | 
		   !is.finite(D0_M) | 
		   !is.finite(D0_k2) | 
		   !is.finite(D0_k3) | 
		   !is.finite(D0_P)
		   ) return(NaN)

		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

		if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
		{
			alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
			initialPenality <- 0
		}else{
			if(clean){initialPenality <- 0}else{
		initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1VKV_Der(0,parameters))^2
																					 + (k2KKK_Der(0,KKK)-k2VKV_Der(0,parameters))^2
																					 + (k3KKK_Der(0,KKK)-k3VKV_Der(0,parameters))^2)
			}
			alphaChiSquare <- 0
		}

		chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
		penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)

		if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}

		if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
	}

	errorKVV_Der <- function(parameters, tpts
						   , premature, mature, alpha
						   , prematureVariance, matureVariance, alphaVariance
						   , KKK = NULL
						   , initialChisquare = NULL
						   , initialDistances = NULL
						   , initialPenalityRelevance = 1
						   , derivativePenalityRelevance = 1
						   , clean
						   )
	{
		if(length(parameters)==13)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k3Parameters <- parameters[8:13]

			totalEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = totalParameters))
			alphaEstimated <- rep(k1Parameters, length(tpts))
			betaEstimated <- sapply(tpts,function(t)impulseModel(x = t, par = k3Parameters))
		
			D0_T <- .DimpulseModel(0, totalParameters)
			D0_k1 <- 0
			D0_k3 <- .DimpulseModel(0, k3Parameters)
		} else {
			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k3Parameters <- parameters[6:9]

			totalEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = totalParameters))
			alphaEstimated <- rep(k1Parameters, length(tpts))
			betaEstimated <- sapply(tpts,function(t)sigmoidModel(x = t, par = k3Parameters))
		
			D0_T <- .DsigmoidModel(0, totalParameters)
			D0_k1 <- 0
			D0_k3 <- .DsigmoidModel(0, k3Parameters)
		}

		gammaEstimated <- sapply(tpts,function(t)k2KVV_Der(t, parameters))
		prematureEstimated <- sapply(tpts,function(t)prematureKVV_Der(x = t, parameters = parameters))
		matureEstimated <- totalEstimated - prematureEstimated

		D0_P <- alphaEstimated[[1]] - gammaEstimated[[1]]*prematureEstimated[[1]]

		betaEstimated[betaEstimated<0] <- NaN
		gammaEstimated[gammaEstimated<0] <- NaN

		if(is.null(alpha)&is.null(alphaVariance))
		{
			alphaEstimated[alphaEstimated<0] <- NaN						
		}

		if(any(!is.finite(alphaEstimated)) | 
		   any(!is.finite(betaEstimated)) | 
		   any(!is.finite(gammaEstimated)) | 
		   any(!is.finite(prematureEstimated)) | 
		   any(!is.finite(matureEstimated)) | 
		   any(!is.finite(totalEstimated)) | 
		   !is.finite(D0_T) | 
		   !is.finite(D0_k1) | 
		   !is.finite(D0_k3) | 
		   !is.finite(D0_P)) return(NaN)

		prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
		matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

		if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
		{
			alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
			initialPenality <- 0
		}else{
			if(clean){initialPenality <- 0}else{
		initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1KVV_Der(0,parameters))^2
																					 + (k2KKK_Der(0,KKK)-k2KVV_Der(0,parameters))^2
																					 + (k3KKK_Der(0,KKK)-k3KVV_Der(0,parameters))^2)
			}
			alphaChiSquare <- 0
		}

		chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
		penalty <- abs(D0_T)+abs(D0_P)+abs(D0_k1)+abs(D0_k3)

		if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}

		if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
	}

		errorVVV_Der <- function(parameters, tpts
							   , premature, mature, alpha
							   , prematureVariance, matureVariance, alphaVariance
							   , KKK = NULL
							   , initialChisquare = NULL
							   , initialDistances = NULL
							   , initialPenalityRelevance = 1
							   , derivativePenalityRelevance = 10^-50
							   , clean)
		{
			if(length(parameters)==18)
			{
				matureParameters <- parameters[1:6]
				matureEstimated <- impulseModel(x = tpts, par = matureParameters)

				D0_M <- .DimpulseModel(0,parameters[1:6])
				D0_k2 <- .DimpulseModel(0,parameters[7:12])
				D0_k3 <- .DimpulseModel(0,parameters[13:18])

			} else {
				matureParameters <- parameters[1:4]
				matureEstimated <- sigmoidModel(x = tpts, par = matureParameters)

				D0_M <- .DsigmoidModel(0,parameters[1:4])
				D0_k2 <- .DsigmoidModel(0,parameters[5:8])
				D0_k3 <- .DsigmoidModel(0,parameters[9:12])
			}

			D0_P <- .DprematureVVV_Der(0, parameters)

			prematureEstimated <- prematureVVV_Der(x = tpts, parameters = parameters)
			alphaEstimated <- k1VVV_Der(x = tpts, parameters = parameters)

			alphaEstimated[alphaEstimated<0] <- NaN
			prematureEstimated[prematureEstimated<0] <- NaN
			matureEstimated[matureEstimated<0] <- NaN

			if(any(!is.finite(alphaEstimated)) | 
			   any(!is.finite(prematureEstimated)) | 
			   any(!is.finite(matureEstimated)) | 
			   !is.finite(D0_M) | 
			   !is.finite(D0_k2) | 
			   !is.finite(D0_k3) | 
			   !is.finite(D0_P)
			   ) return(NaN)

			prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
			matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

			if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
			{
				alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
				initialPenality <- 0
			}else{
				if(clean){initialPenality <- 0}else{
				initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1VVV_Der(0,parameters))^2
																							   + (k2KKK_Der(0,KKK)-k2VVV_Der(0,parameters))^2
																							   + (k3KKK_Der(0,KKK)-k3VVV_Der(0,parameters))^2)
				}
				alphaChiSquare <- 0
			}

			chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
			penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)

			if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}
			
			if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
		}

############################################
### Kinetic rates integrative functions ####
############################################

	k1KKK_Int <- function(x, parameters)
	{
	  parameters[1]
	}

	k2KKK_Int <- function(x, parameters)
	{
	  parameters[2]
	}

	k3KKK_Int <- function(x, parameters)
	{
	  parameters[3]
	}

###########################################
### Kinetic rates derivative functions ####
###########################################

	k1KKK_Der <- function(x, parameters)
	{
	  parameters[1]*parameters[3]
	}

	k2KKK_Der <- function(x, parameters)
	{
	  parameters[2]
	}

	k3KKK_Der <- function(x, parameters)
	{
	  parameters[3]
	}

	k1VKK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7]
			k3Parameters <- parameters[8]

			return(.D2impulseModel(x, matureParameters)/k2Parameters + .DimpulseModel(x, matureParameters)*(1+k3Parameters/k2Parameters) + k3Parameters*impulseModel(x, matureParameters))
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5]
			k3Parameters <- parameters[6]

			return(.D2sigmoidModel(x, matureParameters)/k2Parameters + .DsigmoidModel(x, matureParameters)*(1+k3Parameters/k2Parameters) + k3Parameters*sigmoidModel(x, matureParameters))
		}
	}

	k2VKK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			k2Parameters <- parameters[7]

		}else{
			k2Parameters <- parameters[5]
		}
	}

	k3VKK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			k3Parameters <- parameters[8]

		}else{
			k3Parameters <- parameters[6]
		}
	}	

	k1KVK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			alphaParameters <- parameters[7]
			return(alphaParameters)
		}else{
			alphaParameters <- parameters[5]
			return(alphaParameters)
		}
	}

	k2KVK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k3Parameters <- parameters[8]

			return((k1Parameters - (.DimpulseModel(x, totalParameters) + .D2impulseModel(x, totalParameters)/k3Parameters))/prematureKVK_Der(x = x, parameters = parameters))
		}else{
			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k3Parameters <- parameters[6]

			return((k1Parameters - (.DsigmoidModel(x, totalParameters) + .D2sigmoidModel(x, totalParameters)/k3Parameters))/prematureKVK_Der(x = x, parameters = parameters))
		}
	}

	k3KVK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			k3Parameters <- parameters[8]

		}else{
			k3Parameters <- parameters[6]
		}
	}

	k1KKV_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			alphaParameters <- parameters[7]
			return(alphaParameters)
		}else{
			alphaParameters <- parameters[5]
			return(alphaParameters)
		}
	}

	k2KKV_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			k2Parameters <- parameters[8]
		} else {
			k2Parameters <- parameters[6]
		}
	}

	k3KKV_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k2Parameters <- parameters[8]

			return((k1Parameters - .DimpulseModel(x, totalParameters))/(impulseModel(x, totalParameters) - k1Parameters/k2Parameters))
		} else {
			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k2Parameters <- parameters[6]

			return((k1Parameters - .DsigmoidModel(x, totalParameters))/(sigmoidModel(x, totalParameters) - k1Parameters/k2Parameters))
		}
	}

	k1VVK_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7:12]
			k3Parameters <- parameters[13]

			return(.D2impulseModel(x, matureParameters)/impulseModel(x, k2Parameters) + 
				   .DimpulseModel(x, matureParameters)*(1 - .DimpulseModel(x, k2Parameters)/impulseModel(x, k2Parameters)^2 + k3Parameters/impulseModel(x, k2Parameters)) + 
				   impulseModel(x, matureParameters)*(k3Parameters - (k3Parameters*.DimpulseModel(x, k2Parameters))/impulseModel(x, k2Parameters)^2))
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5:8]
			k3Parameters <- parameters[9]

			return(.D2sigmoidModel(x, matureParameters)/sigmoidModel(x, k2Parameters) + 
				   .DsigmoidModel(x, matureParameters)*(1 - .DsigmoidModel(x, k2Parameters)/sigmoidModel(x, k2Parameters)^2 + k3Parameters/sigmoidModel(x, k2Parameters)) + 
				   sigmoidModel(x, matureParameters)*(k3Parameters - (k3Parameters*.DsigmoidModel(x, k2Parameters))/sigmoidModel(x, k2Parameters)^2))
		}
	}

	k2VVK_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			k2Parameters <- parameters[7:12]

			return(impulseModel(x, k2Parameters))
		}else{
			k2Parameters <- parameters[5:8]

			return(sigmoidModel(x, k2Parameters))

		}
	}

	k3VVK_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			k3Parameters <- parameters[13]
		}else{
			k3Parameters <- parameters[9]
		}
	}

	k1VKV_Der <- function(x, parameters)
	{

		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7]
			k3Parameters <- parameters[8:13]

			return(.D2impulseModel(x, matureParameters)/k2Parameters + 
				   .DimpulseModel(x, matureParameters)*(1 + impulseModel(x, k3Parameters)/k2Parameters) + 
				   impulseModel(x, matureParameters)*(.DimpulseModel(x, k3Parameters)/k2Parameters + impulseModel(x, k3Parameters)))
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5]
			k3Parameters <- parameters[6:9]

			return(.D2sigmoidModel(x, matureParameters)/k2Parameters + 
				   .DsigmoidModel(x, matureParameters)*(1 + sigmoidModel(x, k3Parameters)/k2Parameters) + 
				   sigmoidModel(x, matureParameters)*(.DsigmoidModel(x, k3Parameters)/k2Parameters + sigmoidModel(x, k3Parameters)))
		}
	}

	k2VKV_Der <- function(x, parameters)
	{

		if(length(parameters)==13)
		{
			k2Parameters <- parameters[7]
		}else{
			k2Parameters <- parameters[5]
		}
	}

	k3VKV_Der <- function(x, parameters)
	{

		if(length(parameters)==13)
		{
			k3Parameters <- parameters[8:13]

			return(impulseModel(x, k3Parameters))
		}else{
			k3Parameters <- parameters[6:9]

			return(sigmoidModel(x, k3Parameters))
		}
	}

	k1KVV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			alphaParameters <- parameters[7]
			return(alphaParameters)
		}else{
			alphaParameters <- parameters[5]
			return(alphaParameters)
		}
	}

	k2KVV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k3Parameters <- parameters[8:13]

			return((k1Parameters - (.DimpulseModel(x, totalParameters) + .D2impulseModel(x, totalParameters)/impulseModel(x, k3Parameters) - ((.DimpulseModel(x, totalParameters) - k1Parameters)*.DimpulseModel(x, k3Parameters))/impulseModel(x, k3Parameters)^2))/(impulseModel(x, totalParameters) + (.DimpulseModel(x, totalParameters) - k1Parameters)/impulseModel(x, k3Parameters)))
		}else{
			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k3Parameters <- parameters[6:9]

			return((k1Parameters - (.DsigmoidModel(x, totalParameters) + .D2sigmoidModel(x, totalParameters)/sigmoidModel(x, k3Parameters) - ((.DsigmoidModel(x, totalParameters) - k1Parameters)*.DsigmoidModel(x, k3Parameters))/sigmoidModel(x, k3Parameters)^2))/(sigmoidModel(x, totalParameters) + (.DsigmoidModel(x, totalParameters) - k1Parameters)/sigmoidModel(x, k3Parameters)))
		}
	}

	k3KVV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			k3Parameters <- parameters[8:13]

			return(impulseModel(x, k3Parameters))
		}else{
			k3Parameters <- parameters[6:9]

			return(sigmoidModel(x, k3Parameters))
		}
	}

	k1VVV_Der <- function(x, parameters)
	{

		if(length(parameters)==18)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7:12]
			k3Parameters <- parameters[13:18]

			return(.D2impulseModel(x, matureParameters)/impulseModel(x, k2Parameters) +
				   .DimpulseModel(x, matureParameters)*(1 - .DimpulseModel(x, k2Parameters)/impulseModel(x, k2Parameters)^2 + impulseModel(x, k3Parameters)/impulseModel(x, k2Parameters)) + 
				   impulseModel(x, matureParameters)*(.DimpulseModel(x, k3Parameters)/impulseModel(x, k2Parameters) + impulseModel(x, k3Parameters) - (impulseModel(x, k3Parameters)*.DimpulseModel(x, k2Parameters))/impulseModel(x, k2Parameters)^2 ))
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5:8]
			k3Parameters <- parameters[9:12]

			return(.D2sigmoidModel(x, matureParameters)/sigmoidModel(x, k2Parameters) +
				   .DsigmoidModel(x, matureParameters)*(1 - .DsigmoidModel(x, k2Parameters)/sigmoidModel(x, k2Parameters)^2 + sigmoidModel(x, k3Parameters)/sigmoidModel(x, k2Parameters)) + 
				   sigmoidModel(x, matureParameters)*(.DsigmoidModel(x, k3Parameters)/sigmoidModel(x, k2Parameters) + sigmoidModel(x, k3Parameters) - (sigmoidModel(x, k3Parameters)*.DsigmoidModel(x, k2Parameters))/sigmoidModel(x, k2Parameters)^2 ))
		}
	}

	k2VVV_Der <- function(x, parameters)
	{

		if(length(parameters)==18)
		{
			k2Parameters <- parameters[7:12]
			return(impulseModel(x, k2Parameters))
		}else{
			k2Parameters <- parameters[5:8]
			return(sigmoidModel(x, k2Parameters))
		}
	}

	k3VVV_Der <- function(x, parameters)
	{

		if(length(parameters)==18)
		{
			k3Parameters <- parameters[13:18]
			return(impulseModel(x, k3Parameters))
		}else{
			k3Parameters <- parameters[9:12]
			return(sigmoidModel(x, k3Parameters))
		}
	}

#######################################
### Premature derivative functions ####
#######################################

	prematureKKK_Der <- function(x, parameters)
	{
	  matureParameters <- parameters[1]
	  k2Parameters <- parameters[2]
	  k3Parameters <- parameters[3]
	  
	  return((k3Parameters*matureParameters)/k2Parameters)
	}

	prematureVKK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7]
			k3Parameters <- parameters[8]

			return((.DimpulseModel(x, matureParameters) + k3Parameters * impulseModel(x, matureParameters))/k2Parameters)
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5]
			k3Parameters <- parameters[6]

			return((.DsigmoidModel(x, matureParameters) + k3Parameters * sigmoidModel(x, matureParameters))/k2Parameters)
		}	
	}

	.DprematureVKK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7]
			k3Parameters <- parameters[8]

			M <- impulseModel(x, matureParameters)
			DM <- .DimpulseModel(x, matureParameters)
			D2M <- .D2impulseModel(x, matureParameters)
			
			k3 <- k3Parameters
			Dk3 <- 0
			
			k2 <- k2Parameters
			Dk2 <- 0

		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5]
			k3Parameters <- parameters[6]

			M <- sigmoidModel(x, matureParameters)
			DM <- .DsigmoidModel(x, matureParameters)
			D2M <- .D2sigmoidModel(x, matureParameters)
			
			k3 <- k3Parameters
			Dk3 <- 0
			
			k2 <- k2Parameters
			Dk2 <- 0
		}
		return((k2*(M*Dk3+k3*DM+D2M)-Dk2*(k3*M+DM))/k2^2)		
	}

	prematureKVK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k3Parameters <- parameters[8]

			return(impulseModel(x, totalParameters) + (.DimpulseModel(x, totalParameters) - k1Parameters)/k3Parameters)
		}else{

			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k3Parameters <- parameters[6]

			return(sigmoidModel(x, totalParameters) + (.DsigmoidModel(x, totalParameters) - k1Parameters)/k3Parameters)
		}
	}

	prematureKKV_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k2Parameters <- parameters[8]

			return(k1Parameters/k2Parameters)
		}else{

			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k2Parameters <- parameters[6]

			return(k1Parameters/k2Parameters)
		}
	}

	prematureVVK_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7:12]
			k3Parameters <- parameters[13]

			return((.DimpulseModel(x, matureParameters) + k3Parameters * impulseModel(x, matureParameters))/impulseModel(x, k2Parameters))
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5:8]
			k3Parameters <- parameters[9]

			return((.DsigmoidModel(x, matureParameters) + k3Parameters * sigmoidModel(x, matureParameters))/sigmoidModel(x, k2Parameters))
		}	
	}

	.DprematureVVK_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7:12]
			k3Parameters <- parameters[13]

			M <- impulseModel(x, matureParameters)
			DM <- .DimpulseModel(x, matureParameters)
			D2M <- .D2impulseModel(x, matureParameters)
			
			k2 <- impulseModel(x, k2Parameters)
			Dk2 <- .DimpulseModel(x, k2Parameters)

			k3 <- k3Parameters
			Dk3 <- 0

		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5:8]
			k3Parameters <- parameters[9]

			M <- sigmoidModel(x, matureParameters)
			DM <- .DsigmoidModel(x, matureParameters)
			D2M <- .D2sigmoidModel(x, matureParameters)
			
			k2 <- sigmoidModel(x, k2Parameters)
			Dk2 <- .DsigmoidModel(x, k2Parameters)
			
			k3 <- k3Parameters
			Dk3 <- 0
		}
		return((k2*(M*Dk3+k3*DM+D2M)-Dk2*(k3*M+DM))/k2^2)
	}

	prematureVKV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7]
			k3Parameters <- parameters[8:13]

			return((.DimpulseModel(x, matureParameters) + impulseModel(x, k3Parameters) * impulseModel(x, matureParameters))/k2Parameters)
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5]
			k3Parameters <- parameters[6:9]

			return((.DsigmoidModel(x, matureParameters) + sigmoidModel(x, k3Parameters) * sigmoidModel(x, matureParameters))/k2Parameters)
		}	
	}

	.DprematureVKV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7]
			k3Parameters <- parameters[8:13]

			M <- impulseModel(x, matureParameters)
			DM <- .DimpulseModel(x, matureParameters)
			D2M <- .D2impulseModel(x, matureParameters)
			
			k2 <- k2Parameters
			Dk2 <- 0
			
			k3 <- impulseModel(x, k3Parameters)
			Dk3 <- .DimpulseModel(x, k3Parameters)

		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5]
			k3Parameters <- parameters[6:9]

			M <- sigmoidModel(x, matureParameters)
			DM <- .DsigmoidModel(x, matureParameters)
			D2M <- .D2sigmoidModel(x, matureParameters)
						
			k2 <- k2Parameters
			Dk2 <- 0

			k3 <- sigmoidModel(x, k3Parameters)
			Dk3 <- .DsigmoidModel(x, k3Parameters)

		}
		return((k2*(M*Dk3+k3*DM+D2M)-Dk2*(k3*M+DM))/k2^2)
	}

	prematureKVV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			totalParameters <- parameters[1:6]
			k1Parameters <- parameters[7]
			k3Parameters <- parameters[8:13]

			return(impulseModel(x, totalParameters) + (.DimpulseModel(x, totalParameters) - k1Parameters)/impulseModel(x, k3Parameters))
		}else{

			totalParameters <- parameters[1:4]
			k1Parameters <- parameters[5]
			k3Parameters <- parameters[6:9]

			return(sigmoidModel(x, totalParameters) + (.DsigmoidModel(x, totalParameters) - k1Parameters)/sigmoidModel(x, k3Parameters))
		}
	}

	prematureVVV_Der <- function(x, parameters)
	{
		if(length(parameters)==18)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7:12]
			k3Parameters <- parameters[13:18]

			return((.DimpulseModel(x, matureParameters) + impulseModel(x, k3Parameters) * impulseModel(x, matureParameters))/impulseModel(x, k2Parameters))
		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5:8]
			k3Parameters <- parameters[9:12]

			return((.DsigmoidModel(x, matureParameters) + sigmoidModel(x, k3Parameters) * sigmoidModel(x, matureParameters))/sigmoidModel(x, k2Parameters))
		}	
	}

	.DprematureVVV_Der <- function(x, parameters)
	{
		if(length(parameters)==18)
		{
			matureParameters <- parameters[1:6]
			k2Parameters <- parameters[7:12]
			k3Parameters <- parameters[13:18]

			M <- impulseModel(x, matureParameters)
			DM <- .DimpulseModel(x, matureParameters)
			D2M <- .D2impulseModel(x, matureParameters)
			
			k3 <- impulseModel(x, k3Parameters)
			Dk3 <- .DimpulseModel(x, k3Parameters)
			
			k2 <- impulseModel(x, k2Parameters)
			Dk2 <- .DimpulseModel(x, k2Parameters)


		}else{
			matureParameters <- parameters[1:4]
			k2Parameters <- parameters[5:8]
			k3Parameters <- parameters[9:12]

			M <- sigmoidModel(x, matureParameters)
			DM <- .DsigmoidModel(x, matureParameters)
			D2M <- .D2sigmoidModel(x, matureParameters)
			
			k3 <- sigmoidModel(x, k3Parameters)
			Dk3 <- .DsigmoidModel(x, k3Parameters)
			
			k2 <- sigmoidModel(x, k2Parameters)
			Dk2 <- .DsigmoidModel(x, k2Parameters)
		}
		return((k2*(M*Dk3+k3*DM+D2M)-Dk2*(k3*M+DM))/k2^2)		
	}

####################################
### Mature derivative functions ####
####################################

	matureKKK_Der <- function(x, parameters)
	{
	  matureParameters <- parameters[1]
	  return(matureParameters)
	}

	matureVKK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			matureParameters <- parameters[1:6]
			return(impulseModel(x, matureParameters))
		}else{
			matureParameters <- parameters[1:4]
			return(sigmoidModel(x, matureParameters))
		}	
	}

	matureKVK_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			return(impulseModel(x, totalParameters) - prematureKVK_Der(x, parameters))
		}else{

			totalParameters <- parameters[1:4]
			return(sigmoidModel(x, totalParameters) - prematureKVK_Der(x, parameters))
		}
	}

	matureKKV_Der <- function(x, parameters)
	{
		if(length(parameters)==8)
		{
			totalParameters <- parameters[1:6]
			return(impulseModel(x, parameters) - prematureKKV_Der(x, parameters))
		}else{

			totalParameters <- parameters[1:4]
			return(sigmoidModel(x, parameters) - prematureKKV_Der(x, parameters))
		}
	}

	matureVVK_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			return(impulseModel(x, matureParameters))
		}else{
			matureParameters <- parameters[1:4]
			return(sigmoidModel(x, matureParameters))
		}	
	}

	matureVKV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			matureParameters <- parameters[1:6]
			return(impulseModel(x, matureParameters))
		}else{
			matureParameters <- parameters[1:4]
			return(sigmoidModel(x, matureParameters))
		}	
	}

	matureKVV_Der <- function(x, parameters)
	{
		if(length(parameters)==13)
		{
			totalParameters <- parameters[1:6]
			return(impulseModel(x, totalParameters) - prematureKVV_Der(x, parameters))
		}else{
			totalParameters <- parameters[1:4]
			return(sigmoidModel(x, totalParameters) - prematureKVV_Der(x, parameters))
		}
	}

	matureVVV_Der <- function(x, parameters)
	{
		if(length(parameters)==18)
		{
			matureParameters <- parameters[1:6]
			return(impulseModel(x, matureParameters))
		}else{
			matureParameters <- parameters[1:4]
			return(sigmoidModel(x, matureParameters))
		}	
	}