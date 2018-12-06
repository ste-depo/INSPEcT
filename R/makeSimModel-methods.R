#' @rdname makeSimModel
#'
#' @description
#' This method allow the creation of synthesis, degradation and processing rates for a certain number of genes.
#' The rates are created according to the distributions of the real data-set which is given as an input of the
#' method. Different proportions of constant varying rates can be set and a new vector of time points can be
#' provided. This method has to be used before the \code{\link{makeSimDataset}} method.
#' @param object An object of class INSPEcT
#' @param nGenes A numeric with the number of synthtic genes to be created
#' @param newTpts A numeric verctor with time points of the synthtic dataset, if NULL the time points of the real dataset will be used
#' @param probs A numeric vector wich describes the probability of a rate to be constant, shaped like a sigmoid or like an impulse model
#' @param na.rm A logical that set whether missing values in the real dataset should be removed
#' @param seed A numeric to obtain reproducible results
#' @details The method \code{\link{makeSimModel}} generates an object of class INSPEcT_model that stores the parametric functions to genrate clean rates of a time-course. To any of the rates also a noise variance is associate but not used yet. In a typical workflow the output of \code{\link{makeSimModel}} is the input of the method \code{\link{makeSimDataset}}, that build the noisy rates and concentrations, given a specified number of replicates.
#' @return An object of class INSPEcT_model with synthetic rates
#' @seealso \code{\link{makeSimDataset}}
#' @examples
#' nascentInspObj <- readRDS(system.file(package='INSPEcT', 'nascentInspObj.rds'))
#' simRates<-makeSimModel(nascentInspObj, 1000, seed=1)

#' tpts <- simRates@params$tpts
#' 
#' simData2rep_Nascent <- makeSimDataset(object=simRates,tpts=tpts,nRep=3,NoNascent=FALSE,seed=1)
#' simData2rep_Nascent <- modelRates(simData2rep_Nascent[1:10], seed=1)
setMethod('makeSimModel', 'INSPEcT', function(object
											, nGenes
											, probs=c(constant=.5
													, sigmoid=.3
													, impulse=.2)
											, na.rm=TRUE
											, seed=NULL)
{

	tpts <- object@tpts

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
	#
	suppressWarnings(
		out <- .makeSimData(nGenes
						  , tpts
						  , concentrations
						  , rates
						  , probs = probs
						  , na.rm = na.rm
						  , seed = seed)
	)

	# arrange simdataSpecs form .makeSimData
	simdataSpecs <- out$simdataSpecs
	simdataSpecs <- lapply(simdataSpecs, function(x) list(x))
	#
	newObject <- new('INSPEcT_model')
	newObject@ratesSpecs <- simdataSpecs
	newObject@params$sim$flag <- TRUE
	newObject@params$sim$foldchange <- out$simulatedFC
	newObject@params$sim$noiseVar <- out$noiseVar
	newObject@params$sim$noiseFunctions <- out$noiseFunctions
	newObject@params$tpts <- tpts

	if(length(out$simdataSpecs)>nGenes){return(newObject[1:nGenes])} #Return only the required number of genes or less
	return(newObject)

	})

# .rowVars <- function(data, na.rm=FALSE)
# {
# 	n <- ncol(data)
# 	(rowMeans(data^2, na.rm=na.rm) - rowMeans(data, na.rm=na.rm)^2) * n / (n-1)
# }

# .rowVars2 <- function(data)
# {
# 	n <- ncol(data)
# 	rowMeans((data-rowMeans(data))^2) * n / (n-1)
# }

.which.quantile <- function(values, distribution=values, na.rm=FALSE, 
	quantiles=100)
# given a number of quantiles, returns the quantile each element of 'values'
# belongs to. By default the quantiles are evaluested on 'values' itself, 
# otherwise can be calculated on a specified 'distribution'
{
	if( is.null(distribution)) distribution <- values
	# calculate the boundaries of the quantiles
	qtlBoundaries <- quantile(distribution, probs=seq(0, 1, by=1/quantiles), 
		na.rm=na.rm)
	qtlBoundaries <- cbind(
		qtlBoundaries[-(length(qtlBoundaries))], 
		qtlBoundaries[-1]
		)
	# assign -Inf to the lower boundary
	# and +Inf to the upper one
	qtlBoundaries[1, 1] <- -Inf
	qtlBoundaries[nrow(qtlBoundaries), 2] <- Inf
	ix <- sapply(values, function(x) 
		min(which(x < qtlBoundaries[,2])))

	return(ix)

}

.makeSimData <- function(nGenes
					   , tpts
					   , concentrations
					   , rates
					   , probs=c(constant=.5,sigmoid=.3,impulse=.2)
					   , na.rm=FALSE
					   , seed=NULL)
{

	######################
	# internal functions ###
	##########################

	sampleNormQuantile <- function(values_subject
								 , dist_subject
								 , dist_object, na.rm=FALSE
								 , quantiles=100)
	# sample values from the distribution OBJECT given that some values of the 
	# distribution SUBJECT are known.
	{

		quantileMeanVar <- function(dist_subject, dist_object=NULL, na.rm=FALSE, quantiles)
		# for each quantile of the distribution SUBJECT gives back
		# the mean and the standard deviation of distribution OBJECT
		{
			if( is.null(dist_object)) 
				dist_object <- dist_subject
			idx <- .which.quantile(values=dist_subject, na.rm=na.rm, 
				quantiles=quantiles)
			distMean <- tapply(dist_object, idx, mean)
			distVar <- tapply(dist_object, idx, stats::var)
			return(cbind(mean=distMean, var=distVar))
		}
		## linearize the time-course matrices into vectors of values
		dist_subject <- c(dist_subject)
		dist_object  <- c(dist_object)
		if( na.rm ) {
			tokeep <- is.finite(dist_subject) & is.finite(dist_object)
			dist_subject <- dist_subject[tokeep]
			dist_object  <- dist_object[tokeep]
		}
		## number of quantile can't be too large in order that each quantile
		## can host al least 4 elements
		quantiles <- min(quantiles, floor(length(dist_subject)/4))
		## sample the values
		idx <- .which.quantile(
			values         = values_subject 
			, distribution = dist_subject 
			, quantiles    = quantiles
			, na.rm        = na.rm
			)
		qmv <- quantileMeanVar(
			dist_subject  = dist_subject 
			, dist_object = dist_object
			, quantiles   = quantiles
			, na.rm       = na.rm
			)
		values_object <- rep(NA, length(values_subject))

		for(i in 1:quantiles)
		{
			nobjects <- length(which(idx==i))
			if(nobjects!=0){
				values_object[idx==i] <- rnorm(
					nobjects
					, mean=qmv[as.character(i),'mean'] 
					, sd=sqrt(qmv[as.character(i),'var']) 
					)
			}
		}

		# M.F. cambio perché si piantava con il gamma, utilizzo gli indici di riga per ovviare a problemi con quantili non definiti
		# for(i in 1:quantiles)
		# {
		# 	nobjects <- length(which(idx==i))
		# 	values_object[idx==i] <- rnorm(
		# 		nobjects
		# 		, mean=qmv[i,'mean'] 
		# 		, sd=sqrt(qmv[i,'var']) 
		# 		)
		# }
		return(values_object)
	}

	sampleNorm2DQuantile <- function(values_subject1
								   , values_subject2
								   , dist_subject1
								   , dist_subject2
								   , dist_object
								   , na.rm=FALSE
								   , quantiles=10)
	# sample values from the distribution OBJECT given that some values odf the 
	# distribution SUBJECT are known.
	{
		dist_subject1 <- c(dist_subject1)
		dist_subject2 <- c(dist_subject2)
		dist_object   <- c(dist_object)
		if( na.rm ) {
			tokeep <- is.finite(dist_subject1) & is.finite(dist_subject2) & 
				is.finite(dist_object)
			dist_subject1 <- dist_subject1[tokeep]
			dist_subject2 <- dist_subject2[tokeep]
			dist_object <- dist_object[tokeep]
		}
		## number of quantile can't be too large in order that each quantile
		## can host al least 4 elements
		quantiles <- min(quantiles, floor(sqrt(length(dist_subject1)/4)))
		##
		idx1 <- .which.quantile(values_subject1, dist_subject1, 
			na.rm=na.rm, quantiles=quantiles)
		idx2 <- .which.quantile(values_subject2, dist_subject2, 
			na.rm=na.rm, quantiles=quantiles)

		quantile2DMeanVar <- function(dist_subject1
									, dist_subject2
									, dist_object
									, na.rm=FALSE
									, quantiles=100)
		# for each quantile of the distribution SUBJECT1 and SUBJECT2 gives
		# back the mean and the standard deviation of distribution OBJECT. 
		# Returns the two square matrices of mean and variance corresponding 
		# to each pair of quantiles of SUBJECT1 and SUBJECT2.
		{
			idx1 <- .which.quantile(dist_subject1, na.rm=na.rm, quantiles=quantiles)
			idx2 <- .which.quantile(dist_subject2, na.rm=na.rm, quantiles=quantiles)
			meansTab <- matrix(NA, nrow=quantiles, ncol=quantiles)
			varsTab  <- matrix(NA, nrow=quantiles, ncol=quantiles)
			for(i1 in unique(idx1))
			{
				for(i2 in unique(idx2))
				{
					# belonging to either quantiles
					ix <- idx1 == i1 & idx2 == i2
					meansTab[i1,i2] <- mean(dist_object[ix])
					varsTab[i1,i2]  <- stats::var(dist_object[ix])
				}
			}
			# fill the missing values
			na.fill <- function(mat)
			# Fill the NA values of a matrix with the mean of the surroundings.
			# Iterates until all the missing values are filled.
			{
				if( all(is.na(mat))) return(mat)
				nRow <- nrow(mat)
				nCol <- ncol(mat)
				while(length(which(is.na(mat))) > 0){
					for(i in 1:nrow(mat)){
						for(j in 1:ncol(mat)){
							if( is.na(mat[i,j]))
							{
								idx_top <- max(1,i-1)
								idx_bottom <- min(nRow,i+1)
								idx_left <- max(1,j-1)
								idx_right <- min(nCol,j+1)
								surroundingRows <- idx_top:idx_bottom
								surroundingCols <- idx_left:idx_right
								mat[i,j] <- mean(
									mat[surroundingRows,surroundingCols], 
									na.rm=TRUE
									)
							} } } }
				return(mat)
			}
			meansTab <- na.fill(meansTab)
			varsTab  <- na.fill(varsTab)
			return(list(mean=meansTab,var=varsTab))
		}

		q2dmv <- quantile2DMeanVar(dist_subject1, dist_subject2, dist_object, 
			na.rm=na.rm, quantiles=quantiles)
		sampledValues <- sapply(1:length(idx1), function(i) {
			qtMean <- q2dmv$mean[idx1[i], idx2[i]]
			qtVar <- q2dmv$var[idx1[i], idx2[i]]
			########## Why not sqrt(qtVar) ?????????????
			# changed to sd=sqrt(qtVar), previously was:
			# sd=qtVar
			return(rnorm(1, mean=qtMean, sd=sqrt(qtVar)))
			})
		return(sampledValues)		
	}

	generateParams <- function(tpts
							 , sampled_val
							 , log2foldchange
							 , probs=c(constant=.5,sigmoid=.3,impulse=.2))
	# given a vector of absolute values and a vector of log2foldchanges
	# create parametric functions (either constant, sigmoidal or impulse, 
	# according to probs) and evaluate them at tpts.
	{

		##############################
		#### define local functions ####
		##################################

		generateImpulseParams <- function(tpts, sampled_val, log2foldchange)
		# Given an absolute value and a value of log2fold change sample a set 
		# of parameters for the impulse function.
		{

			n <- length(sampled_val)

			# sample the delta of the two responses between a range that 
			# is considered valid to reproduce the expected fold change
			# (intervals that are too small or too large compared to the 
			# length of the dynamics can lead to a reduced fold change)
			time_span <- diff(range(tpts))
			delta_max <- time_span / 1.5
			delta_min <- time_span / 6.5

			# sample the delta of the response (difference between first and 
			# second response) uniformly over the confidence interval
			sampled_deltas <- runif( n, min=delta_min, max=delta_max)

			# the time of first response is sampled in order to include the 
			# whole response within the time course
			time_of_first_response <- sapply(
				max(tpts) - sampled_deltas
				, function(max_first_response) 
					runif( 1,min=min(tpts),max=max_first_response)
				)
			# second response is then trivial
			time_of_second_response <- time_of_first_response + sampled_deltas
			
			# the slope of the response is inversely proportional to the delta
			# sampled (the shorter is the response the fastest it has to be, 
			# in order to satisfy the fold change)
			slope_of_response <- time_span / sampled_deltas

			initial_values      <- sampled_val * 2^(-log2foldchange)
			intermediate_values <- sampled_val
			end_values          <- initial_values

			impulsepars <- cbind(
				initial_values
				, intermediate_values
				, end_values
				, time_of_first_response
				, time_of_second_response
				, slope_of_response
				)

			return(impulsepars)

		}

		generateSigmoidParams <- function(tpts, sampled_val, log2foldchange)
		# Given an absolute value and a value of log2fold change sample a set 
		# of parameters for the sigmoid function.
		{

			n <- length(sampled_val)
			time_span <- diff(range(tpts))

			# sample the time uniformely
			time_of_response <- runif( n, min=min(tpts), max=max(tpts))

			# slope of response must be high if the time of response is close 
			# to one of the two boundaries
			distance_from_boundary <- apply(
				cbind(
					time_of_response - min(tpts)
					, max(tpts) - time_of_response
				),1,min)
			slope_of_response <- time_span / distance_from_boundary

			initial_values <- sampled_val * 2^(-log2foldchange)
			end_values     <- sampled_val

			sigmoidpars <- cbind(
				initial_values
				, end_values
				, time_of_response
				, slope_of_response
				)

			return(sigmoidpars)

		}

		#####################################################
		# body of the 'generateParams' function starts here ###
		#########################################################

		nGenes <- length(sampled_val)
		# 
		n_constant <- round(nGenes * probs['constant'])
		n_sigmoid  <- round(nGenes * probs['sigmoid'])
		n_impulse  <- nGenes - (n_constant + n_sigmoid)

		# initialize
		params <- as.list(rep(NA,nGenes))

		# constant: choose the one with the lower absoulute fold change to be 
		# constant
		constant_idx <- 1:nGenes %in% order(abs(log2foldchange))[1:n_constant]
		if( any(constant_idx) )
		{
			params[constant_idx] <- lapply(sampled_val[constant_idx], 
				function(val) 
					list(type='constant', fun=constantModelP , params=val, df=1)
					)
		}
		# impulse varying, fist guess
		impulse_idx <- 1:nGenes %in% sample(which(!constant_idx), n_impulse)
		if( any(impulse_idx) ) {
			impulseParamGuess <- lapply(which(impulse_idx), 
				function(i)
					generateImpulseParams(tpts, sampled_val[i], log2foldchange[i])
					)
			valuesImpulse <- do.call('rbind', lapply(impulseParamGuess, 
				function(par) impulseModel(tpts,par)
				))
			expectedFC  <- abs(log2foldchange[impulse_idx])
			simulatedFC <- apply(valuesImpulse, 1, 
				function(x) diff(log2(range(x)))
				)
			# due to the nature of the impulse function, by average the real fold
			# change (the one generated by the sampled data) is lower than the
			# one expected. For this reason, we calculate the factor of scale
			# between the real and expected fold changes and we generate new 
			# data

		 	factor_of_correction <- lm(simulatedFC ~ expectedFC)$coefficients[2]
			params[impulse_idx] <- lapply(
				which(impulse_idx)
				, function(i) list(
					type='impulse'
					, fun=impulseModelP
					, params=generateImpulseParams(
						tpts 
						, sampled_val[i]
						, log2foldchange[i] * factor_of_correction
						)
					, df=6
					)
				)
		}
		
		# sigmoid
		sigmoid_idx <- !constant_idx & !impulse_idx
		if( any(sigmoid_idx) ) {
			params[sigmoid_idx] <- lapply(
				which(sigmoid_idx)
				, function(i) list(
					type='sigmoid'
					, fun=sigmoidModelP
					, params=generateSigmoidParams(
						tpts
						, sampled_val[i] 
						, log2foldchange[i] 
						)
					, df=4
					)			
				)
		}

		# # report true foldchanges
		# simulatedFC <- apply(values, 1, function(x) diff(log2(range(x))))

		return(params)

	}

	#########################################
	# body of the main function starts here ###
	#############################################

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

	# make twice the number of genes and then select only the valid ones
	nGenes.bkp <- nGenes
	nGenes <- nGenes * 2
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
	# fold change
	message('sampling fold changes from rates distribution...')
	# get log2 fc distribution
	alphaLog2FC <- log2(alpha[,-1]) - log2(alpha[,1])
	betaLog2FC  <- log2(beta[,-1])  - log2(beta[,1])
	gammaLog2FC <- log2(gamma[,-1]) - log2(gamma[,1])
	# sample log2 fc
	alphaSimLog2FC <- sampleNormQuantile(
		values_subject = log2(alphaVals) 
		, dist_subject = log2(alpha[,-1]) 
		, dist_object  = alphaLog2FC
		, na.rm=na.rm
		)
	betaSimLog2FC  <- sampleNorm2DQuantile(
		values_subject1   = log2(betaVals)
		, values_subject2 = alphaSimLog2FC
		, dist_subject1   = log2(beta[,-1])
		, dist_subject2   = alphaLog2FC
		, dist_object     = betaLog2FC
		, na.rm=na.rm, quantiles=50
		)
	gammaSimLog2FC  <- sampleNorm2DQuantile(
		values_subject1   = log2(gammaVals)
		, values_subject2 = alphaSimLog2FC
		, dist_subject1   = log2(gamma[,-1])
		, dist_subject2   = alphaLog2FC
		, dist_object     = gammaLog2FC
		, na.rm=na.rm, quantiles=50
		)

	# transform time in log scale
	log_shift <- find_tt_par(tpts)
	x <- time_transf(tpts, log_shift)
	# generate alpha, beta, gamma
	message('generating rates time course...')
	alphaParams <- generateParams(x, alphaVals, alphaSimLog2FC, probs)
	betaParams  <- generateParams(x, betaVals , betaSimLog2FC , probs)
	gammaParams <- generateParams(x, gammaVals, gammaSimLog2FC, probs)
	# generate total and preMRNA from alpha,beta,gamma
	paramSpecs <- lapply(1:nGenes, 
		function(i) 
			list(alpha=alphaParams[[i]]
			   , beta=betaParams[[i]]
			   , gamma=gammaParams[[i]]))

	out <- lapply(1:nGenes, function(i){
			tryCatch(
				.makeModel(tpts, paramSpecs[[i]], log_shift, 
					time_transf, ode, .rxnrate),error=function(e){cbind(time=rep(NaN,length(tpts))
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
	nGenes <- nGenes.bkp

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

	return(list(
		simdataSpecs=paramSpecs
		, simulatedFC=simulatedFC
		, noiseVar=noiseVar
		, noiseFunctions = list(alpha = alphaFitVarianceLaw, preMRNA = preFitVarianceLaw, total = totalFitVarianceLaw)
		))
}
