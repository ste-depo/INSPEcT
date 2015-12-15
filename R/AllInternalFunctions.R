.find_tt_par <- function(tpts)
{
	cvLogTpts <- function(a , tpts) {
		newtime <- log2(tpts + a )
		stats::sd(diff(newtime)) / mean(diff(newtime))}
	optimize(f=cvLogTpts, interval=c(0,5), tpts=tpts )$minimum
}

.time_transf <- function(t, log_shift) 
{
	newtime <- log2(t+log_shift)
	return(newtime)
} 

.getRatesAndConcentrationsFromRpkms <- function(totRpkms, labeledRpkms, tpts
	, tL=NULL, degDuringPulse=FALSE, simulatedData=FALSE, parallelize=TRUE
	, totalMedianNorm=TRUE, labeledMedianNorm=FALSE, totalSF=NULL, labeledSF=NULL
	, steadyStateMode=0) 
{

	## steadyStateMode
	## 0 -> derivative at time zero is zero
	## 1 -> derivative at time zero is non zero only for labeledSF
	## 2 -> derivative at time zero is non zero only for labeledSF and rates

	## the control during the workflow will be based on the 
	## negation of simulatedData (that is realData)
	realData <- !simulatedData

	## set parameters for the parallel mode
	if( parallelize ) {
		maxcores <- parallel::detectCores()
		if( Sys.info()[['sysname']] == 'Windows' ) {
			applyfun <- bplapply
			options(MulticoreParam=SnowParam(workers=maxcores))
			# print(bpparam())
		} else {
			applyfun <- parallel::mclapply
			options(mc.cores=maxcores)
			# print(bpparam())
		}
	} else applyfun <- lapply

	## eventually set the 'simple' mode
	if( is.null(totRpkms$introns) | is.null(labeledRpkms$introns) ) 
		onlyExonsMode <- TRUE
	else
		onlyExonsMode <- FALSE

	## retrieve gene names from rownames of exon total rpkms
	## (all the matrices have the same dimensions and rownames)
	geneNames <- rownames(totRpkms$exons)

	###########
	# internal functions
	########
	speedyVar <- function(x) sum((x - mean.default(x))^2)/length(x[-1])

	##### total fraction
	## rename total exons
	Texo  <- totRpkms$exons
	## assign total exons variance
	if( !is.null(totRpkms$exons_var) ) {
		Texo_var <- totRpkms$exons_var
		if( !identical(dim(Texo_var), dim(Texo)) )
			stop('Variance must be a matrix with same dimensions of rpkms.')
		else
			totalVariance <- TRUE
	} else 
		totalVariance <- FALSE
	if( !onlyExonsMode ) {
		## rename total introns
		Tint  <- totRpkms$introns
		## assign total introns variance
		if( totalVariance )
			Tint_var <- totRpkms$introns_var
	}
	##### labeled fraction
	## rename labeled exons
	Lexo <- labeledRpkms$exons
	## assign total exons variance
	if( !is.null(labeledRpkms$exons_var) ) {
		Lexo_var <- labeledRpkms$exons_var
		if( !identical(dim(Lexo_var), dim(Lexo)) )
			stop('Variance must be a matrix with same dimensions of rpkms.')
		else
			labeledVarince <- TRUE
	} else 
		labeledVarince <- FALSE
	if( !onlyExonsMode ) {
		## rename total introns
		Lint <- labeledRpkms$introns
		## assign total introns variance
		if( labeledVarince )
			Lint_var <- labeledRpkms$introns_var
	}

	if( realData )
	##############
	# NORMALIZE THE DATA  - only needed with real data, synthetic data are 
	# 				generated already scaled
	#########################
	{

		#########################
		#### scale total rpkms ###
		#########################
		# eventually inter normalize total with a user defined scaling factor
		# and adjust variance accordingly
		if( !is.null(totalSF) & totalMedianNorm )
			stop('If totalSF is provided totalMedianNorm must be set to FALSE')

		if( is.null(totalSF) & !totalMedianNorm ) {
			totalSF <- rep(1, length(tpts))
		}

		if( totalMedianNorm ) {
			ssq.residuals <- function(a, x, y) 
				stats::median(x - a*y, na.rm=TRUE)^2
			totalSF <- sapply(1:ncol(Texo), 
				function(i) 
					optimize(ssq.residuals, c(-10,10), x=Texo[,1]
						, y=Texo[,i] )$minimum)				
		}

		if( !is.null(totalSF) ) {
			Texo <- t(t(Texo)*totalSF)
			if( totalVariance )
				Texo_var <- apply(t(t(Texo_var)*totalSF^2), 1, mean, na.rm=TRUE)
			if( !onlyExonsMode ) {
				Tint <- t(t(Tint)*totalSF)
				if( totalVariance )
					Tint_var <- apply(t(t(Tint_var)*totalSF^2), 1, mean, na.rm=TRUE)
			}
		}

		if( !totalVariance ) {
			## evaluate the variance now
			Texo_var <- apply( Texo , 1 , speedyVar )
			if( !onlyExonsMode )
				Tint_var <- apply( Tint , 1 , speedyVar )

		}

		## in case varince is still null (meaning that the user
		## didn't provided one) calculate it now

		# ## scale the 4sU according to the labelling time
		# Lexo <- Lexo/tL
		# if( !onlyExonsMode )
		# 	Tint <- Tint/tL
		# if( labeledVarince ) {
		# 	Lexo_var <- Lexo_var/tL^2
		# 	if( !onlyExonsMode )
		# 		Lint_var <- Lint_var/tL^2
		# }

		#########################
		#### scale 4sU rpkms ###
		#########################

		## only exons mode
		if( onlyExonsMode ) 
		{

			if( !is.null(labeledSF) & labeledMedianNorm )
				stop('If labeledSF is provided labeledMedianNorm must be set to FALSE')

			# normalize 4su according to median
			if( labeledMedianNorm ) {
				ssq.residuals <- function(a, x, y) stats::median(x - a*y, na.rm=TRUE)^2
				# calculate a scaling factor for each time point
				labeledSF <- sapply(1:ncol(Lexo), function(i)
						optimize(ssq.residuals, c(-10,10), x=Lexo[,1], y=Lexo[,i])$minimum)
			}

			# normalize and estimate variance
			if( !is.null(labeledSF) ) {
				Lexo <- t(t(Lexo)*labeledSF)
				if( labeledVarince )
					Lexo_var <- apply(t(t(Lexo_var)*labeledSF^2), 1, mean, na.rm=TRUE)
				else
					Lexo_var <- apply(Lexo, 1, speedyVar)					
			} else {
			# only estimate variance
				if( labeledVarince )
					Lexo_var <- apply(Lexo_var, 1, mean, na.rm=TRUE)
				else
					Lexo_var <- apply(Lexo, 1, speedyVar)
			}

		## introns and exons mode
		} else {

			# derivatives from time course
			TintDer <- as.matrix(sapply(1:nrow(Tint), 
				function(i) {
					if( all(is.finite(Tint[i,] ) ) ) {
						spfun <- splinefun(tpts, Tint[i,])
						return(spfun(tpts, deriv=1) )
					} else return(rep(NA, length(tpts)) )
				}) )
			if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
			if( steadyStateMode == 0 ) TintDer[, 1] <- 0 

			TexoDer <- as.matrix(sapply(1:nrow(Texo), 
				function(i) {
					if( all(is.finite(Texo[i,] ) ) ) {
						spfun <- splinefun(tpts, Texo[i,])
						return(spfun(tpts, deriv=1) )
					} else return(rep(NA, length(tpts)) )
				}) )
			if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
			if( steadyStateMode == 0 ) TexoDer[, 1] <- 0

			if( labeledMedianNorm )
				warning('When introns of labelled and total fraction are provided normalization of 
					labelled fraction is computed internally. "labeledMedianNorm" argument is ignored.')

			if( !is.null(labeledSF) )
				warning('When introns of labelled and total fraction are provided normalization of 
					labelled fraction is computed internally. "labeledSF" argument is ignored.')

			message('Calculating scaling factor between total and 4su libraries...')

			#########################
			##### local functions ######
			########################
			fLint <- function(Lexo,gamma,tL) Lexo/(tL*gamma)*(1-exp(-gamma*tL))
				# given the number of labeled molecules, gamma and tL gives back the 
				# number of processed molecules
			fGamma <- function(Lint, Lexo, tL, maxGamma=1e3) {
				# given the number of labeled molecules and the number of processed
				# molecules gives back the processing rate. It's an inverse function,
				# therefore, the precision (step of the interval evaluated) and the
				# max value where gamma is possibly evaluated should be provided
				if( is.na(Lint) | is.na(Lexo) ) return(NA)
				if( Lint >= Lexo ) return(0)
				if( Lint == 0 ) return(Inf)
				errorfun <- function(gamma, Lexo, Lint, tL ) 
					(Lint-fLint(Lexo, gamma, tL))^2
				optimize(errorfun, c(0,maxGamma), Lexo=Lexo, Lint=Lint
					, tL=tL)$minimum
			}
			# calculate the factor which bring the median of the residuals between
			# the modeled preMRNA levels and the measured to zero
			sq.median.resids <- function(sf, P, dP, alpha, gamma) 
				sapply(sf, function(i) {
					t1 <- dP + gamma*P
					t2 <- i*alpha
					idx <- is.finite(t1) & is.finite(t2) & t1 > 0 & t2 > 0
					resids <- t1[idx] - t2[idx]
					stats::median(resids , na.rm=TRUE)^2
				})

			##################
			#### scale data ###
			##################
			# preMRNA derivative as splines 
			# (force the first time point to have derivative zero )
			# estimate of alpha and gamma from 4sU data

			gammaTC <- sapply(
				1:length(tpts), function(j) 
					unlist(applyfun(1:nrow(Lint), 
						function(i, Lint, Lexo, tL) fGamma(Lint[i,j] , Lexo[i,j] , tL)
						, Lint=Lint, Lexo=Lexo, tL=tL ))
					) 
			# scale factor 
			labeledSF_prior <- sapply(1:ncol(Tint), function(j)
				optimize(sq.median.resids, c(0.01,100), P=Tint[,j], dP=TintDer[,j], 
					alpha=Lexo[,j]/tL, gamma=gammaTC[,j] )$minimum)

			if( degDuringPulse ) {

				if( labeledMedianNorm )
					warning('When introns of labelled and total fraction are provided normalization of 
						labelled fraction is computed internally. "labeledMedianNorm" argument is ignored.')

				if( !is.null(labeledSF) )
					warning('When introns of labelled and total fraction are provided normalization of 
						labelled fraction is computed internally. "labeledSF" argument is ignored.')

				# message('Calculating scaling factor between total and 4su libraries...')

				## select the 500 most synthesized genes at time zero
				## to calculate the scaling factor
				if( nrow(Lexo) > 500 ) {
					geneSubset <- order(Texo[,1],decreasing=TRUE)[1:500]
				} else {
					geneSubset <- 1:nrow(Lexo)
				}

				## re-calculate scaling factor

				## degradation during pulse
				# system with 4su scaling factor as a 4th variable to be identified
				initOneGene <- function(i, j) {
					c(
						Lexo[i,j]/tL
						,(Lexo[i,j]/tL-TexoDer[i,j])/(Texo[i,j]-Tint[i,j])
						,(Lexo[i,j]/tL-TintDer[i,j])/Tint[i,j]
						)	
				}
				sys4suScale <- function(x) {
					y <- numeric(4)
					y[1] <- TintDer[i,j] - x[1] + x[3]*Tint[i,j]
					y[2] <- TexoDer[i,j] - x[1] + x[2]*(Texo[i,j]-Tint[i,j])
					y[3] <- x[4]*Lint[i,j] - x[1]/x[3]*(1-exp(-x[3]*tL))
					y[4] <- x[4]*Lexo[i,j] - (x[1]*exp(-x[2]*tL)*(x[3]^2-x[2]^2*exp((x[2]-x[3])*tL)+
						x[2]^2*exp(x[2]*tL)-x[3]^2*exp(x[2]*tL)))/(x[2]*(x[2]-x[3])*x[3])
					y/c(x[1],x[1],x[4]*Lint[i,j],x[4]*Lexo[i,j])
				}
				## get only 4sU scales
				nGenes <- nrow(Lexo)
				nTpts <- ncol(Lexo)
				labeledSfep <- matrix(NA, nrow=nGenes, ncol=nTpts)
				labeledSf <- matrix(NA, nrow=nGenes, ncol=nTpts)
				capture.output(suppressWarnings({
					for( i in geneSubset ) {
						for( j in 1:nTpts ) {
							init <- initOneGene(i,j)*labeledSF_prior[j]
							# if( all(is.finite(init)) ) {
							mrOut <- tryCatch(
								multiroot(
									sys4suScale
									, c(init,labeledSF_prior[j])
									# , control = list(maxit=1e3)
									# , positive=TRUE
									)
								, error=function(e)
									list(
										root=rep(NA, 4)
										, estim.precis=NA
										)
									)
							labeledSf[i,j] <- mrOut$root[4]
							labeledSfep[i,j] <- mrOut$estim.precis
							# } else {
							# 	labeledSf[i,j] <- NA
							# 	labeledSfep[i,j] <- NA
							# }
						}
					}
				}))
				
				## chose the best resolved genes to estimate 
				## from them the scale factor
				epTsh <- apply(labeledSfep,2,quantile,probs=.75,na.rm=TRUE)
				ix <- t(apply(labeledSfep,1, function(x) x>epTsh))
				labeledSf[ix] <- NA
				labeledSF <- apply(labeledSf, 2, stats::median, na.rm=TRUE)

				# Lexo <- t(t(Lexo)*labeledSF)
				# if( labeledVarince )
				# 	Lexo_var <- apply(t(t(Lexo_var)*labeledSF^2), 1, mean, na.rm=TRUE)
				# else
				# 	Lexo_var <- apply(Lexo, 1, speedyVar)

			} else {

				labeledSF <- labeledSF_prior

				# if( labeledMedianNorm )
				# 	warning('When introns of labelled and total fraction are provided normalization of 
				# 		labelled fraction is computed internally. "labeledMedianNorm" argument is ignored.')

				# if( !is.null(labeledSF) )
				# 	warning('When introns of labelled and total fraction are provided normalization of 
				# 		labelled fraction is computed internally. "labeledSF" argument is ignored.')

				# message('Calculating scaling factor between total and 4su libraries...')

				# #########################
				# ##### local functions ######
				# ########################
				# fLint <- function(Lexo,gamma,tL) Lexo/(tL*gamma)*(1-exp(-gamma*tL))
				# 	# given the number of labeled molecules, gamma and tL gives back the 
				# 	# number of processed molecules
				# fGamma <- function(Lint, Lexo, tL, maxGamma=1e3) {
				# 	# given the number of labeled molecules and the number of processed
				# 	# molecules gives back the processing rate. It's an inverse function,
				# 	# therefore, the precision (step of the interval evaluated) and the
				# 	# max value where gamma is possibly evaluated should be provided
				# 	if( is.na(Lint) | is.na(Lexo) ) return(NA)
				# 	if( Lint >= Lexo ) return(0)
				# 	if( Lint == 0 ) return(Inf)
				# 	errorfun <- function(gamma, Lexo, Lint, tL ) 
				# 		(Lint-fLint(Lexo, gamma, tL))^2
				# 	optimize(errorfun, c(0,maxGamma), Lexo=Lexo, Lint=Lint
				# 		, tL=tL)$minimum
				# }
				# # calculate the factor which bring the median of the residuals between
				# # the modeled preMRNA levels and the measured to zero
				# sq.median.resids <- function(sf, P, dP, alpha, gamma) 
				# 	sapply(sf, function(i) {
				# 		t1 <- dP + gamma*P
				# 		t2 <- i*alpha
				# 		idx <- is.finite(t1) & is.finite(t2) & t1 > 0 & t2 > 0
				# 		resids <- t1[idx] - t2[idx]
				# 		stats::median(resids , na.rm=TRUE)^2
				# 	})

				# ##################
				# #### scale data ###
				# ##################
				# # preMRNA derivative as splines 
				# # (force the first time point to have derivative zero )
				# # estimate of alpha and gamma from 4sU data

				# gammaTC <- sapply(
				# 	1:length(tpts), function(j) 
				# 		unlist(applyfun(1:nrow(Lint), 
				# 			function(i, Lint, Lexo, tL) fGamma(Lint[i,j] , Lexo[i,j] , tL)
				# 			, Lint=Lint, Lexo=Lexo, tL=tL ))
				# 		) 
				# # scale factor 
				# labeledSF <- sapply(1:ncol(Tint), function(j)
				# 	optimize(sq.median.resids, c(0.01,100), P=Tint[,j], dP=TintDer[,j], 
				# 		alpha=Lexo[,j]/tL, gamma=gammaTC[,j] )$minimum)
				# scaling (only on Lexo, gammaTC is independent to scaling factor )

				# Lexo <- t(t(Lexo)*labeledSF)
				# if( labeledVarince )
				# 	Lexo_var <- apply(t(t(Lexo_var)*labeledSF^2), 1, mean, na.rm=TRUE)
				# else
				# 	Lexo_var <- apply(Lexo, 1, speedyVar)

			}

		}

	## simulated data
	} else {

		# in case of synthetic data the time course is already scaled
		# therefore just rename the variables and compute varince in 
		# case is not provided
		alphaTC <- Lexo
		if( labeledVarince )
			alphaTC_var <- rowMeans(Lexo_var)
		else alphaTC_var <- apply(alphaTC, 1, speedyVar)
		if( totalVariance ) {
			Texo_var <- rowMeans(Texo_var)
			Tint_var <- rowMeans(Tint_var)
		} else {
			Texo_var <- apply( Texo , 1 , speedyVar )
			Tint_var <- apply( Tint , 1 , speedyVar )
		}
		# scaling factor for synthetic dataset is meaningless
		labeledSF <- rep(1, length(tpts))

		# derivatives from time course
		TintDer <- as.matrix(sapply(1:nrow(Tint), 
			function(i) {
				if( all(is.finite(Tint[i,] ) ) ) {
					spfun <- splinefun(tpts, Tint[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}) )
		if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
		if( steadyStateMode == 0 ) TintDer[, 1] <- 0 

		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}) )
		if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		if( steadyStateMode == 0 ) TexoDer[, 1] <- 0


	}

	###################################
	## estimate degradation rates #########
	#####################################

	# only exons mode
	if( onlyExonsMode ) {

		# in case introns are not provided calculate degradation rates with
		# a reduced function that doesn't takes preMRNAs into account and
		# return results without estimating processing rates

		# derivatives from time course
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}) )
		if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		if( steadyStateMode == 0 ) TexoDer[, 1] <- 0


		## degradation during pulse
		## estimate both scaling and rates at the same time
		if( degDuringPulse ) {

			message('Estimating all rates...')

			## degradation during pulse
			initOneGene <- function(i, j) {
				c(
					Lexo[i,j]*labeledSF[j]/tL
					,(Lexo[i,j]*labeledSF[j]/tL-TexoDer[i,j])/Texo[i,j]
					)	
			}
			sys4suSmall <- function(x) {
				y <- numeric(2)
				y[1] <- TexoDer[i,j] - x[1] + x[2]*Texo[i,j]
				y[2] <- Lexo[i,j] - x[1]/x[2]*(1-exp(-x[2]*tL))
				y
			}
			## get only the rates
			nGenes <- nrow(Lexo)
			nTpts <- ncol(Lexo)
			ratesEstimPrec <- matrix(NA, nrow=nGenes, ncol=nTpts)
			alphaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
			betaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
			capture.output(suppressWarnings({
				for( i in 1:nGenes ) {
					for( j in 1:nTpts ) {
						init <- initOneGene(i,j)
						# if( all(is.finite(init)) ) {
						mrOut <- tryCatch(
							multiroot(
								sys4suSmall
								, initOneGene(i,j)
								# , labeledSF=labeledSF
								# , control = list(maxit=1e3)
								# , positive=TRUE
								)
							, error=function(e)
								list(
									root=rep(NA, 2)
									, estim.precis=NA
									)
								)
						ratesEstimPrec[i,j] <- mrOut$estim.precis
						alphaTC[i,j] <- mrOut$root[1]
						betaTC[i,j] <- mrOut$root[2]
						# } else {
						# 	labeledSfep[i,j] <- NA
						# 	alphaTC[i,j] <- NA
						# 	betaTC[i,j] <- NA
						# 	gammaTC[i,j] <- NA
						# }
					}
				}
			}))

			# ## keep only the best resolved rates (put to NA the worst 10%
			# ## of rates in time points later than zero)
			# epTsh <- apply(labeledSfep,2,quantile,probs=.9,na.rm=TRUE)
			# ix <- t(apply(labeledSfep,1, function(x) x>epTsh))
			# ix[,1] <- FALSE
			# alphaTC[ix] <- NA
			# betaTC[ix] <- NA
			# gammaTC[ix] <- NA

			## put negative values to NA and rise a 
			## warning if they are more than 20% of a specific rate
			## synthesis
			ix <- alphaTC<0
			negativePerc <- length(which(ix))/length(ix)*100
			if( negativePerc>20 ) {
				warning(paste(round(negativePerc), 
					'% of the synthesis rates are negative. Putting them to NA.'))
			}
			alphaTC[ix] <- NA
			## degradation
			ix <- betaTC<0
			negativePerc <- length(which(ix))/length(ix)*100
			if( negativePerc>20 ) {
				warning(paste(round(negativePerc), 
					'% of the degradation rates are negative. Putting them to NA.'))
			}
			betaTC[ix] <- NA

			# ## some rates are unbelivably large, respect to the assumption
			# ## that no degradation occurs during the pulse
			# ## based on an aprroximated assumtion of what is the expected
			# ## ratio between the two we discard the ones that have a ratio 
			# ## superior than the threshold
			# approxExpectedRatio <- function(beta, tL) {
			# 	beta*tL/(1-exp(-beta*tL))
			# }
			# ## ratio between median accounts for the conversion between the scaling factors
			# ## of the noDegr mode and the Degr mode
			# noDegrAlphaTC <- Lexo/tL
			# noDegrAlphaTC <- noDegrAlphaTC*
			# 	stats::median(alphaTC, na.rm=TRUE)/stats::median(noDegrAlphaTC, na.rm=TRUE)
			# measuredRatio <- alphaTC/noDegrAlphaTC
			# qt <- seq(.9,1,by=.001)
			# res <- numeric(length(qt))
			# for(i in 1:length(qt)) {
			# 	q <- qt[i]
			# 	tOut <- table(
			# 		measuredRatio>
			# 			approxExpectedRatio(
			# 				quantile(betaTC, na.rm=TRUE, probs=q)
			# 				, tL)
			# 			)
			# 	res[i] <- tOut['FALSE']/sum(tOut)
			# }
			# if( all(res>qt) ) {
			# 	# plot(qt, res, type='l')
			# 	# abline(0,1,col='red')
			# 	warning('Too many rates are over the expected using degDuringPulse mode.')
			# }
			# tsh <- qt[min(which(res>qt))]
			# tshRatio <- approxExpectedRatio(
			# 	quantile(betaTC, na.rm=TRUE, probs=tsh)
			# 	, tL)

			# alphaTC[measuredRatio>tshRatio] <- NA
			# betaTC[measuredRatio>tshRatio] <- NA

			## recalculate the variance
			if( labeledVarince )
				alphaTC_var <- Lexo_var/tL^2
			else
				alphaTC_var <- apply(alphaTC, 1, speedyVar)

			## set preMRNA and gamma to NA
			Tint <- matrix(NA, nrow(Texo), ncol(Texo))
			Tint_var <- rep(NA, length(Texo_var))
			gammaTC <- matrix(NA, nrow(betaTC), ncol(betaTC))


		# assume that no degradation occur during pulse
		} else {

			inferKBetaFromIntegral <- function(tpts, alpha, totalRNA, maxBeta=150) 
			{
				solveFun <- function(beta, t0, t1, alpha_t0, alpha_t1, X_t0, X_t1 ) 
				{
					intsol <- function(t , m , q , beta )
						(m*t*beta + q*beta - m ) * exp(beta*t ) / (beta^2 )
					defintsol <- function(t0, t1, m, q, beta ) 
						intsol(t1, m, q, beta ) - intsol(t0, m, q, beta )
					#
					mAlpha <- (alpha_t0 - alpha_t1 ) / (t0 - t1 )
					qAlpha <- alpha_t0 - mAlpha * t0
					#
					X_t1*exp(beta*t1 ) - X_t0*exp(beta*t0 ) - 
						defintsol(t0, t1, mAlpha, qAlpha, beta )
				}
				maxBetas <- seq(5, maxBeta, by=5)
				# applyfun <- if( parallel ) bplapply else lapply
				lapply(2:length(tpts), function(j) {
					unlist(applyfun(1:nrow(alpha),
						function(i) {
							k <- 1
							solutionFound <- FALSE
							while((!solutionFound) & k <= length(maxBetas)) {
								try({
									solution <- uniroot(solveFun
										, c(1e-5, maxBetas[k])
										, t0 = tpts[j-1]
										, t1 = tpts[j]
										, alpha_t0 = alpha[i,j-1]
										, alpha_t1 = alpha[i,j]
										, X_t0 = totalRNA[i,j-1]
										, X_t1 = totalRNA[i,j]
										)#$root
									solutionFound = TRUE
									}, silent=TRUE)
								k <- k + 1
							}
							if( solutionFound ) 
								return(solution) 
							else 
								return(list(root=NA, estim.prec=NA))
						}))
					})
			}
			## calculate alpha and recalculate the variance
			alphaTC <- Lexo/tL
			if( labeledVarince )
				alphaTC_var <- Lexo_var/tL^2
			else
				alphaTC_var <- apply(alphaTC, 1, speedyVar)
			## calculate beta
			message('Estimating degradation rates...')
			if( steadyStateMode == 1 ) TexoDer[,1] <- 0
			betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / Texo[,1]
			betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
			if( length(tpts)>1 ) {
				betaOut <- inferKBetaFromIntegral(tpts, alphaTC, Texo)
				betaTC <- cbind(betaT0, 
					sapply(betaOut, function(x) x[names(x)=='root'])
					)
				ratesEstimPrec <- cbind(0,
					sapply(betaOut, function(x) x[names(x)=='estim.prec'])
					)
			} else {
				betaTC <- as.matrix(betaT0)
				ratesEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
			}
			## set preMRNA and gamma to NA
			Tint <- matrix(NA, nrow(Texo), ncol(Texo))
			Tint_var <- rep(NA, length(Texo_var))
			gammaTC <- matrix(NA, nrow(betaTC), ncol(betaTC))
		}

	# introns and exons mode
	} else {

		# in case introns are provided calculate degradation rates with
		# a function that takes preMRNAs into account and
		# claculate processing rates

		## degradation during pulse
		## estimate both scaling and rates at the same time
		if( degDuringPulse ) {

			# Lexo <- alphaTC
			# Lexo_var <- alphaTC_var

			message('Estimating all rates...')

			## once set the scale factor calculate the rates
			initOneGene <- function(i, j) {
				c(
					Lexo[i,j]*labeledSF[j]/tL
					,(Lexo[i,j]*labeledSF[j]/tL-TexoDer[i,j])/(Texo[i,j]-Tint[i,j])
					,(Lexo[i,j]*labeledSF[j]/tL-TintDer[i,j])/Tint[i,j]
					)	
			}
			sys4su <- function(x, labeledSF) {
				y <- numeric(3)
				y[1] <- TintDer[i,j] - x[1] + x[3]*Tint[i,j]
				y[2] <- TexoDer[i,j] - x[1] + x[2]*(Texo[i,j]-Tint[i,j])
				y[3] <- labeledSF[j]*Lexo[i,j] - (x[1]*exp(-x[2]*tL)*(x[3]^2-x[2]^2*exp((x[2]-x[3])*tL)+
					x[2]^2*exp(x[2]*tL)-x[3]^2*exp(x[2]*tL)))/(x[2]*(x[2]-x[3])*x[3])
				y
			}
			## get only the rates
			nGenes <- nrow(Lexo)
			nTpts <- ncol(Lexo)
			ratesEstimPrec <- matrix(NA, nrow=nGenes, ncol=nTpts)
			alphaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
			betaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
			gammaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
			capture.output(suppressWarnings({
				for( i in 1:nGenes ) {
					for( j in 1:nTpts ) {
						# if( all(is.finite(init)) ) {
						mrOut <- tryCatch(
							multiroot(
								sys4su
								, initOneGene(i,j)
								, labeledSF=labeledSF
								# , control = list(maxit=1e3)
								# , positive=TRUE
								)
							, error=function(e)
								list(
									root=rep(NA, 4)
									, estim.precis=NA
									)
								)
						ratesEstimPrec[i,j] <- mrOut$estim.precis
						alphaTC[i,j] <- mrOut$root[1]
						betaTC[i,j] <- mrOut$root[2]
						gammaTC[i,j] <- mrOut$root[3]
						# } else {
						# 	labeledSfep[i,j] <- NA
						# 	alphaTC[i,j] <- NA
						# 	betaTC[i,j] <- NA
						# 	gammaTC[i,j] <- NA
						# }
					}
				}
			}))

			# ## keep only the best resolved rates (put to NA the worst 10%
			# ## of rates in time points later than zero)
			# epTsh <- apply(labeledSfep,2,quantile,probs=.9,na.rm=TRUE)
			# ix <- t(apply(labeledSfep,1, function(x) x>epTsh))
			# ix[,1] <- FALSE
			# alphaTC[ix] <- NA
			# betaTC[ix] <- NA
			# gammaTC[ix] <- NA

			## put negative values to NA and rise a 
			## warning if they are more than 20% of a specific rate
			ix <- alphaTC<0 | betaTC<0 | gammaTC<0
			negativePerc <- length(which(ix))/length(ix)*100
			if( negativePerc>20 ) {
				warning(paste(round(negativePerc), 
					'% of the genes contains negative rates. Putting them to NA.'))
			}
			alphaTC[ix] <- NA
			betaTC[ix] <- NA
			gammaTC[ix] <- NA
			ratesEstimPrec[ix] <- NA

			# ## some rates are unbelivably large, respect to the assumption
			# ## that no degradation occurs during the pulse
			# ## based on an aprroximated assumtion of what is the expected
			# ## ratio between the two we discard the ones that have a ratio 
			# ## superior than the threshold
			# approxExpectedRatio <- function(beta, tL) {
			# 	beta*tL/(1-exp(-beta*tL))
			# }
			# ## ratio between median accounts for the conversion between the scaling factors
			# ## of the noDegr mode and the Degr mode
			# noDegrAlphaTC <- Lexo/tL
			# noDegrAlphaTC <- noDegrAlphaTC*
			# 	stats::median(alphaTC, na.rm=TRUE)/stats::median(noDegrAlphaTC, na.rm=TRUE)
			# measuredRatio <- alphaTC/noDegrAlphaTC
			# qt <- seq(.9,1,by=.001)
			# res <- numeric(length(qt))
			# for(i in 1:length(qt)) {
			# 	q <- qt[i]
			# 	tOut <- table(
			# 		measuredRatio>
			# 			approxExpectedRatio(
			# 				quantile(betaTC, na.rm=TRUE, probs=q)
			# 				, tL)
			# 			)
			# 	res[i] <- tOut['FALSE']/sum(tOut)
			# }
			# if( all(res>qt) ) {
			# 	# plot(qt, res, type='l')
			# 	# abline(0,1,col='red')
			# 	warning('Too many rates are over the expected using degDuringPulse mode.')
			# }
			# tsh <- qt[min(which(res>qt))]
			# tshRatio <- approxExpectedRatio(
			# 	quantile(betaTC, na.rm=TRUE, probs=tsh)
			# 	, tL)

			# alphaTC[measuredRatio>tshRatio] <- NA
			# betaTC[measuredRatio>tshRatio] <- NA
			# gammaTC[measuredRatio>tshRatio] <- NA

			## recalculate the variance
			if( labeledVarince )
				alphaTC_var <- apply(t(t(Lexo_var)*(labeledSF/tL)^2), 1, mean, na.rm=TRUE)
			else
				alphaTC_var <- apply(alphaTC, 1, speedyVar)

		# assume that no degradation occur during pulse
		} else {

			inferKBetaFromIntegral <- function(tpts, alpha, totalRNA, maxBeta=150) 
			{
				solveFun <- function(beta, t0, t1, alpha_t0, alpha_t1, X_t0, X_t1 ) 
				{
					intsol <- function(t , m , q , beta )
						(m*t*beta + q*beta - m ) * exp(beta*t ) / (beta^2 )
					defintsol <- function(t0, t1, m, q, beta ) 
						intsol(t1, m, q, beta ) - intsol(t0, m, q, beta )
					#
					mAlpha <- (alpha_t0 - alpha_t1 ) / (t0 - t1 )
					qAlpha <- alpha_t0 - mAlpha * t0
					#
					X_t1*exp(beta*t1 ) - X_t0*exp(beta*t0 ) - 
						defintsol(t0, t1, mAlpha, qAlpha, beta )
				}
				maxBetas <- seq(5, maxBeta, by=5)
				# applyfun <- if( parallel ) bplapply else lapply
				lapply(2:length(tpts), function(j) {
					unlist(applyfun(1:nrow(alpha),
						function(i) {
							k <- 1
							solutionFound <- FALSE
							while((!solutionFound) & k <= length(maxBetas)) {
								try({
									solution <- uniroot(solveFun
										, c(1e-5, maxBetas[k])
										, t0 = tpts[j-1]
										, t1 = tpts[j]
										, alpha_t0 = alpha[i,j-1]
										, alpha_t1 = alpha[i,j]
										, X_t0 = totalRNA[i,j-1]
										, X_t1 = totalRNA[i,j]
										)#$root
									solutionFound = TRUE
									}, silent=TRUE)
								k <- k + 1
							}
							if( solutionFound ) 
								return(solution) 
							else 
								return(list(root=NA, estim.prec=NA))
						}))
					})
			}

			# calculate alpha and recalculate the variance
			alphaTC <- Lexo/tL
			if( labeledVarince )
				alphaTC_var <- Lexo_var/tL^2
			else
				alphaTC_var <- apply(alphaTC, 1, speedyVar)
			alphaTC <- t(t(alphaTC)*labeledSF)
			if( labeledVarince )
				alphaTC_var <- apply(t(t(alphaTC_var)*labeledSF^2), 1, mean, na.rm=TRUE)
			else
				alphaTC_var <- apply(alphaTC, 1, speedyVar)

			# calculate beta
			message('Estimating degradation rates...')
			if( steadyStateMode == 1 ) TexoDer[,1] <- 0
			betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / (Texo[,1] - Tint[,1] )
			betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
			if( length(tpts)>1 ) {
				betaOut <- inferKBetaFromIntegral(tpts, alphaTC, Texo)
				betaTC <- cbind(betaT0, 
					sapply(betaOut, function(x) x[names(x)=='root'])
					)
				betaEstimPrec <- cbind(0,
					sapply(betaOut, function(x) x[names(x)=='estim.prec'])
					)
			} else {
				betaTC <- as.matrix(betaT0)
				betaEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
			}

			###################################
			## estimate processing rates #########
			#####################################

			inferKGammaFromIntegral <- function(tpts, alpha, preMRNA, maxGamma=150)
			####### accurate function for estimating the degradation rates
			####### using the solution of the differential equation system under 
			####### the condtion that processing rate is constant between two 
			####### consecutive time points - more stable that using derivatives
			####### estimates
			{
				solveFun <- function(beta, t0, t1, alpha_t0, alpha_t1, X_t0, X_t1 ) 
				{
					intsol <- function(t , m , q , beta )
						(m*t*beta + q*beta - m ) * exp(beta*t ) / (beta^2 )
					defintsol <- function(t0, t1, m, q, beta ) 
						intsol(t1, m, q, beta ) - intsol(t0, m, q, beta )
					#
					mAlpha <- (alpha_t0 - alpha_t1 ) / (t0 - t1 )
					qAlpha <- alpha_t0 - mAlpha * t0
					#
					X_t1*exp(beta*t1 ) - X_t0*exp(beta*t0 ) - 
						defintsol(t0, t1, mAlpha, qAlpha, beta )
				}
				maxGammas <- seq(5, maxGamma, by=5)
				# applyfun <- if( parallel ) bplapply else lapply
				lapply(2:length(tpts), function(j)
					{
						unlist(applyfun(1:nrow(alpha), function(i) {
							k <- 1
							solutionFound <- FALSE
							while((!solutionFound) & k <= length(maxGammas)) {
								try({
									solution <- uniroot(solveFun
										, c(1e-5, maxGammas[k])
										, t0 = tpts[j-1]
										, t1 = tpts[j]
										, alpha_t0 = alpha[i,j-1]
										, alpha_t1 = alpha[i,j]
										, X_t0 = preMRNA[i,j-1]
										, X_t1 = preMRNA[i,j]
										)
									solutionFound = TRUE
									}, silent=TRUE)
								k <- k + 1
							}
							if( solutionFound )
								return(solution)
							else
								return(list(root=NA, estim.prec=NA))
						}))
					})
			}
			# calculate gamma (from  total RNA introns and alphas )
			message('Estimating processing rates...')
			if( steadyStateMode == 1 ) TintDer[,1] <- 0
			gammaT0 <- ( alphaTC[,1] - TintDer[,1] ) / Tint[,1]
			gammaT0[gammaT0 < 0 | !is.finite(gammaT0)] <- NA
			if( length(tpts)>1 ) {
				gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, Tint)
				gammaTC <- cbind(gammaT0, 
					sapply(gammaOut, function(x) x[names(x)=='root'])
					)
				gammaEstimPrec <- cbind(0,
					sapply(gammaOut, function(x) x[names(x)=='estim.prec'])
					)
			} else {
				gammaTC <- as.matrix(gammaT0)
				gammaEstimPrec <- matrix(0, nrow=nrow(gammaTC), ncol=ncol(gammaTC))
			}

			ratesEstimPrec <- betaEstimPrec + gammaEstimPrec

		}

	# end of intron-exons mode
	}

	## remove dimnames from matrices and names from vectors 
	## that will be output
	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(Texo, 'dimnames') <- NULL
	attr(Tint, 'dimnames') <- NULL
	attr(alphaTC_var, 'names') <- NULL
	attr(Texo_var, 'names') <- NULL
	attr(Tint_var, 'names') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL

	## return scaled results and rates (adding NA values 
	## for preMRNA and gamma)
	return(list(
		concentrations=list(
			total=Texo
			, total_var=Texo_var
			, preMRNA=Tint
			, preMRNA_var=Tint_var)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, geneNames=geneNames
		, tpts=tpts
		, labeledSF=labeledSF
		, totalSF=totalSF
		, tL=tL))

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

.makeModel <- function(tpts, hyp, log_shift, .time_transf, ode, .rxnrate)
{
	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(.time_transf(x, log_shift), hyp$alpha$par)
	params$beta  <- function(x) 
		hyp$beta$fun$value(.time_transf(x, log_shift), hyp$beta$par)
	params$gamma <- function(x) 
		hyp$gamma$fun$value(.time_transf(x, log_shift), hyp$gamma$par)
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

.makeSimpleModel <- function(tpts, hyp, log_shift, .time_transf, ode, .rxnrateSimple)
{
	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(.time_transf(x, log_shift), hyp$alpha$params)
	params$beta  <- function(x)
		hyp$beta$fun$value(.time_transf(x, log_shift), hyp$beta$params)
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

.logLikelihood <- function(experiment, model, variance=NULL)
{
    if( is.null(variance)) variance <- stats::var(experiment)
    sum(log(2*pnorm(-abs(experiment-model),mean=0,sd=sqrt(variance))))
}

.inspect.engine <- function(tpts, log_shift, concentrations, rates
	, nInit=10, nIter=300, na.rm=TRUE, nCores=2L
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
		, log_shift, .time_transf, .rxnrateSimple, ode, .makeSimpleModel, .logLikelihood
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
				fun$alpha$value(.time_transf(x, log_shift), splitpar$alpha)
			params$beta <- function(x)
				fun$beta$value(.time_transf(x, log_shift), splitpar$beta)
			#
			cinit <- c(params$alpha(tpts[1]) / params$beta(tpts[1]))
			names(cinit) <- 't'
			model <- ode(y=cinit, times=tpts, func=.rxnrateSimple, parms=params)
			#
			alpha_model <- params$alpha(tpts)
			total_model <- model[,'t']
			#
			D <- chisq(alpha_exp, alpha_model, alpha_var) +
				chisq(total_exp, total_model, total_var)
			# df <- length(alpha_exp) + length(total_exp) - sum(df[1:2])
			# return(log(chisq.test.inspect(D, df)))
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
			, log_shift, .time_transf, ode, .rxnrateSimple)
		logLik <- .logLikelihood(alpha_exp, model$alpha, alpha_var) + 
			.logLikelihood(total_exp, model$total, total_var)
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
		, maxit=500, log_shift, .time_transf, .rxnrate, ode, .makeModel, .logLikelihood)
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
				fun$alpha$value(.time_transf(x, log_shift), splitpar$alpha)
			params$beta  <- function(x)
				fun$beta$value(.time_transf(x, log_shift), splitpar$beta)
			params$gamma <- function(x)
				fun$gamma$value(.time_transf(x, log_shift), splitpar$gamma)
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
			D <- chisq(alpha_exp, alpha_model, alpha_var) +
				chisq(total_exp, total_model, total_var) +
				chisq(preMRNA_exp, preMRNA_model, preMRNA_var)
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
			, log_shift, .time_transf, ode, .rxnrate)
		logLik <- .logLikelihood(alpha_exp, model$alpha, alpha_var) + 
			.logLikelihood(total_exp, model$total, total_var) +
			.logLikelihood(preMRNA_exp, model$preMRNA, preMRNA_var)

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
			.chooseModel, .time_transf, .DimpulseModel, .DsigmoidModel, .constantModelP,
			.emptyGene, .sigmoidModel, .impulseModel, .sigmoidModelP, .impulseModelP,
			.polynomialModelP, .makeModel, .makeSimpleModel, .logLikelihood, .rxnrate,
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
						, .time_transf=.time_transf
						, .sigmoidModel=.sigmoidModel
						, .impulseModel=.impulseModel
						, .sigmoidModelP=.sigmoidModelP
						, .impulseModelP=.impulseModelP
						, .polynomialModelP=.polynomialModelP
						), error=function(e) return(.emptyGene(e)))
				modelSynthesisRatefun <- tryCatch(.chooseModel(tpts=tpts
						, log_shift=log_shift
						, experiment=rates$alpha[i,] 
						, variance=rates$alpha_var[i]
						, na.rm=na.rm, sigmoid=sigmoidSynthesis
						, impulse=TRUE, polynomial=FALSE
						, nInit=nInit, nIter=nIter
						, .time_transf=.time_transf
						, .sigmoidModel=.sigmoidModel
						, .impulseModel=.impulseModel
						, .sigmoidModelP=.sigmoidModelP
						, .impulseModelP=.impulseModelP
						, .polynomialModelP=.polynomialModelP
						), error=function(e) return(.emptyGene(e)))
				modelTotalRNA <- 
					if( testOnSmooth ) {
						modelTotalRNAfun$fun$value(
							.time_transf(tpts, log_shift)
							, modelTotalRNAfun$params)
					} else { concentrations$total[i,] }
				modelSynthesisRate <- 
					if( testOnSmooth ) {
						modelSynthesisRatefun$fun$value(
							.time_transf(tpts, log_shift)
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
							, .time_transf=.time_transf
							, .sigmoidModel=.sigmoidModel
							, .impulseModel=.impulseModel
							, .sigmoidModelP=.sigmoidModelP
							, .impulseModelP=.impulseModelP
							, .polynomialModelP=.polynomialModelP
							), error=function(e) return(.emptyGene(e)))
					modelPreMRNA <- 
						if( testOnSmooth ) 
							modelPreMRNAfun$fun$value(
								.time_transf(tpts, log_shift)
								, modelPreMRNAfun$params) 
						else concentrations$preMRNA[i,]
				}
				#
				if( estimateRatesWith == 'der' ) {
					if( testOnSmooth ) {
						# total RNA derivative
						if( modelTotalRNAfun$type == 'impulse' )
							modelTotalRNAderivative <- .DimpulseModel(
								.time_transf(tpts, log_shift)
								, modelTotalRNAfun$params)
						if( modelTotalRNAfun$type == 'sigmoid' )
							modelTotalRNAderivative <- .DsigmoidModel(
								.time_transf(tpts, log_shift)
								, modelTotalRNAfun$params)
						if( modelTotalRNAfun$type == 'constant' )
							modelTotalRNAderivative <- rep(0, length(tpts))
						if( intExMode ) {
							# pre mRNA derivative
							if( modelPreMRNAfun$type == 'impulse' )
								modelPreMRNAderivative <- .DimpulseModel(
									.time_transf(tpts, log_shift)
									, modelPreMRNAfun$params)
							if( modelPreMRNAfun$type == 'sigmoid' )
								modelPreMRNAderivative <- .DsigmoidModel(
									.time_transf(tpts, log_shift)
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
					alpha=list(fun=.constantModelP
						, type='constant', df=1
						, params=mean(modelSynthesisRate, na.rm=TRUE))
					, beta=list(fun=.constantModelP
						, type='constant', df=1
						, params=mean(firstGuessDegrRate, na.rm=TRUE))
					, gamma=if( intExMode ) {
						list(fun=.constantModelP
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
								, .time_transf=.time_transf
								, .sigmoidModel=.sigmoidModel
								, .impulseModel=.impulseModel
								, .sigmoidModelP=.sigmoidModelP
								, .impulseModelP=.impulseModelP
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
								, .time_transf=.time_transf
								, .sigmoidModel=.sigmoidModel
								, .impulseModel=.impulseModel
								, .sigmoidModelP=.sigmoidModelP
								, .impulseModelP=.impulseModelP
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
								, .time_transf=.time_transf
								, .rxnrate=.rxnrate
								, .makeModel=.makeModel
								, .logLikelihood=.logLikelihood
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
								, .time_transf=.time_transf
								, .rxnrateSimple=.rxnrateSimple
								, .makeSimpleModel=.makeSimpleModel
								, .logLikelihood=.logLikelihood
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

	if( nCores > 1 ) {
		nCores <- min(nCores, parallel::detectCores())
		if( Sys.info()[['sysname']] == 'Windows' ) {
			applyfun <- bplapply
			options(MulticoreParam=SnowParam(workers=nCores))
			# print(bpparam())
		} else {
			applyfun <- parallel::mclapply
			options(mc.cores=nCores)
			# print(bpparam())
		}
	} else applyfun <- lapply
	nGenes <- nrow(rates$alpha)
	tpts <- tpts
	paramSpecs <- applyfun(1:nGenes, modelOneGene, seed=seed, 
		.chooseModel=.chooseModel,
		.time_transf=.time_transf,
		.DimpulseModel=.DimpulseModel,
		.DsigmoidModel=.DsigmoidModel,
		.constantModelP=.constantModelP,
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
		.sigmoidModel=.sigmoidModel,
		.impulseModel=.impulseModel,
		.sigmoidModelP=.sigmoidModelP,
		.impulseModelP=.impulseModelP,
		.polynomialModelP=.polynomialModelP,
		.rxnrate=.rxnrate,
		.rxnrateSimple=.rxnrateSimple,
		.makeModel=.makeModel,
		.makeSimpleModel=.makeSimpleModel,
		.logLikelihood=.logLikelihood
		)
	return(paramSpecs)
}

.brown_method <- function(y, na.rm=FALSE) {
	## initialize answer as a vector of NA
	ans <- rep(NA, nrow(y))
	names(ans) <- rownames(y)
	## count per each row (gene) the number of tests that
	## could actually being performed
	k <- apply(y,1,function(x) length(na.omit(x)))
	## if it's zero for every gene return the NA vector
	if( all(k==0) ) return(ans)
	## otherwise calculate the brown's test for the availbale genes
	## dividing them in cathegories, according to the availbale tests
	y_logical <- matrix(FALSE, nrow(y), ncol(y))
	y_logical[!is.na(y)] <- TRUE
	groups <- unique(as.data.frame(y_logical))
	if( nrow(groups)==1 )
		belongs_to <- rep(1, nrow(y_logical))
	else
		belongs_to <- apply(apply(y_logical,1,function(x) 
			apply(groups, 1, function(y) all(x==y))),2,which)
	fisherSum <- function(x) 
		apply(x,1,function(z) -2*sum(log(z), na.rm=na.rm))
	tmp <- c(-2.59,-2.382,-2.17,-1.946,-1.709,-1.458,-1.194,
		-.916,-.625,-.320,0,.334,.681,1.044,1.421,1.812,
		2.219,2.641,3.079,3.531,4)
	for( i in 1:nrow(groups) ) {
		acceptedTests <- as.logical(groups[i,])
		if( any(acceptedTests) ) {
			ix <- belongs_to==i
			y_group <- y[,acceptedTests, drop=FALSE]
			fX2 <- fisherSum(y_group[ix, , drop=FALSE])
			k_group <- length(which(acceptedTests))
			xout_matrix <- stats::cor(y_group, use='complete.obs')
			xout_vector <- xout_matrix[which(as.vector(lower.tri(xout_matrix)))]
			s2X2 <- 4 * k_group + 2 * sum(approx(seq(-1,1,.1),tmp,xout=xout_vector)$y)
			f <- 2 * (2 * k_group)^2 /s2X2
			correction <- s2X2/(2*2*k_group)
			ans[ix] <- pchisq(fX2/correction, f, lower.tail=FALSE)
		}
	}
	return(ans)
}

.brown_method_mask <- function(y, mask) {
	# if there is only one gene, no assumption about
	# the correlation of the tests (unless it is set and fixed)
	# can be done, therefore fisher is done
	if( nrow(y)==1 ) {
		y_filtered <- y[mask]
		fX2 <- -2*sum(log(y_filtered))
		f <- 2*length(y_filtered)
		out <- pchisq(fX2, f, lower.tail=FALSE)
		names(out) <- rownames(y)
	} else {
		y[c(!mask)] <- NA
		out <- .brown_method(y, na.rm=TRUE)
	}
	return(out)
}

.logLikRatioTest <- function(null, alt)
{
	D <- - 2*null$logLik + 2*alt$logLik
	df_null <- null$alpha$df + null$beta$df +
		ifelse(!is.null(null$gamma$df), null$gamma$df , 0)
	df_alt <- alt$alpha$df + alt$beta$df +
		ifelse(!is.null(alt$gamma$df), alt$gamma$df , 0)
	#chisq.test.inspect(D,  df_alt - df_null)
	pchisq(D, df_alt-df_null, lower.tail=FALSE)
}

############### pointer function

.newPointer <- function(inputValue){  
	object=new.env(parent=globalenv())  
	object$value=inputValue  
	class(object)='pointer'
	return(object)  
}

############### constant

.constantModel <- function(x , par ) rep(par , length(x) )
.constantModelP <- .newPointer(.constantModel)

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
# .sigmoidModel <- sigmoidModelC

.sigmoidModel <- function(x, par) 
{
	# h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
	par[1]+(par[2]-par[1])*(1/(1+exp(-par[4]*(x-par[3]))))
}
## compiled version
.sigmoidModel <- cmpfun(.sigmoidModel)

.DsigmoidModel <- function(x, par) 
{
     h0= par[1]; h1=par[2]; t1=par[3]; b=par[4]
     S= function(b,t) 1/(1+exp(-b*(x-t)))
     dSdx= function(b,t) b/(1/exp(-b*(x-t)) + 2 + exp(-b*(x-t)) )
     s= function(x,t,h,b) h+(h1-h)*S(b,t)
     dsdx= function(x,t,h,b) (h1-h)*dSdx(b,t)
     1/h1*dsdx(x,t1,h0,b)
}

# 'pointer' for the .sigmoidModel function
.sigmoidModelP <- .newPointer(.sigmoidModel)

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
# .impulseModel <- impulseModelC

.impulseModel <- function(x, par) 
{
	# h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
	1/par[2]*(par[1]+(par[2]-par[1])*(1/(1+exp(-par[6]*(x-par[4])))))*
		(par[3]+(par[2]-par[3])*(1/(1+exp(par[6]*(x-par[5])))))
}
## compiled version
.impulseModel <- cmpfun(.impulseModel)

.DimpulseModel <- function(x, par) 
{
     h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
     S= function(b,t) 1/(1+exp(-b*(x-t)))
     dSdx= function(b,t) b/(1/exp(-b*(x-t)) + 2 + exp(-b*(x-t)) )
     s= function(x,t,h,b) h+(h1-h)*S(b,t)
     dsdx= function(x,t,h,b) (h1-h)*dSdx(b,t)
     1/h1*(dsdx(x,t1,h0,b)*s(x,t2,h2,-b) + s(x,t1,h0,b)*dsdx(x,t2,h2,-b) )
}

# 'pointer' for the .impulseModel function
.impulseModelP <- .newPointer(.impulseModel)

############### polynomial

.polynomialModel <- function(x, par)
	sapply(x, function(x_i)
		sum(sapply(1:length(par), function(i) x_i^(i-1) * par[i])))

.polynomialModelP <- .newPointer(.polynomialModel)

chisq <- function(experiment, model, variance=NULL)
{
	if( is.null(variance) )  variance <- stats::var(experiment )
	else if( variance==0 ) variance <- stats::var(experiment )
	sum((experiment - model )^2/variance )
}

.chooseModel <- function(tpts, log_shift, experiment, variance=NULL, na.rm=TRUE
	, sigmoid=TRUE, impulse=TRUE, polynomial=TRUE, nInit=10, nIter=500
	, .time_transf, .impulseModel, .sigmoidModel, .sigmoidModelP, .impulseModelP
	, .polynomialModelP)
#### choose a functional form between impulse and sigmoid according 
#### to the one that has the gratest pvalue in the chi squared test
{



	chisq.test.default <- function(experiment, model, variance=NULL, df)
	{
		if( is.null(variance) ) variance <- stats::var(experiment )
		D = chisq(experiment, model, variance)
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
	im.chisq <- function(par, tpts, experiment, variance=NULL, .impulseModel) 
	{
		 model <- .impulseModel(tpts, par)
		 chisq(experiment, model, variance)
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
	 			, .impulseModel=.impulseModel
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
	sm.chisq <- function(par, tpts, experiment, variance=NULL, .sigmoidModel) 
	{
		 model <- .sigmoidModel(tpts, par)
		 chisq(experiment, model, variance)
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
					, .sigmoidModel=.sigmoidModel
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
		return(list(type='constant', fun=.constantModelP
			, params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
	}
	## 
	if( length(experiment)<=2 ) {
		warning('.chooseModel: less than three time points have a finite value. 
			Impossible evaluate a variable model.
			Returning a constant model.')
		return(list(type='constant', fun=.constantModelP
			, params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
	}

	## re-evaluate flags of function to evaluate according to the lenght 
	## of the experiment
	sigmoid <- sigmoid
	impulse <- impulse & length(experiment)>2
	polynomial <- polynomial & length(experiment)>2

	tptslog <- .time_transf(tpts, log_shift)

	# sigmoid
	if( sigmoid ) {
		outSM  <- sm.optim.chisq(tpts=tptslog, experiment=experiment
			, variance=variance, ninit=nInit, maxit=nIter)
		bestSM <- which.min(unlist(outSM[2,]))
		pvalSM <- chisq.test.default(experiment=experiment
			, model=.sigmoidModel(tptslog, outSM[,bestSM]$par)
			, variance=variance, df=length(outSM[,bestSM]$par))
		dfSM <- length(outSM[,bestSM]$par)
	} else dfSM <- NA
	# impulse
	if( impulse ) {
		outIM  <- im.optim.chisq(tpts=tptslog, experiment=experiment, 
			variance=variance, ninit=nInit, maxit=nIter)
		bestIM <- which.min(unlist(outIM[2,]))
		pvalIM <- chisq.test.default(experiment=experiment
			, model=.impulseModel(tptslog, outIM[,bestIM]$par) 
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
	funcs <- c(.sigmoidModelP, .impulseModelP, .polynomialModelP)
	dfs <- c(dfSM, dfIM, dfPN)
	type <- names(pvals)[which.max(pvals)]
	df   <- dfs[which.max(pvals)]

	if( type=='sigmoid'    ) params <- outSM[,bestSM]$par
	if( type=='impulse'    ) params <- outIM[,bestIM]$par
	if( type=='polynomial' ) params <- outPN[,bestPN]$par

	pval <- pvals[which.max(pvals)]
	fun  <- funcs[[which.max(pvals)]]

	return(list(type=type, fun=fun , params=params, pval=pval, df=df))

}

########################
#### SYNTHETIC DATA ######
########################

.rowVars <- function(data, na.rm=FALSE)
{
	n <- ncol(data)
	(rowMeans(data^2, na.rm=na.rm) - rowMeans(data, na.rm=na.rm)^2) * n / (n-1)
}

.rowVars2 <- function(data)
{
	n <- ncol(data)
	rowMeans((data-rowMeans(data))^2) * n / (n-1)
}

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

.makeSimData <- function(nGenes, tpts, concentrations, rates, newTpts=NULL, 
	probs=c(constant=.5,sigmoid=.3,impulse=.2), na.rm=FALSE, seed=NULL)
{

	######################
	# internal functions ###
	##########################

	sampleNormQuantile <- function(values_subject, dist_subject, dist_object, 
		na.rm=FALSE, quantiles=100)
	# sample values from the distribution OBJECT given that some values of the 
	# distribution SUBJECT are known.
	{

		quantileMeanVar <- function(dist_subject, dist_object=NULL, na.rm=FALSE, 
			quantiles)
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
			values_object[idx==i] <- rnorm(
				nobjects
				, mean=qmv[i,'mean'] 
				, sd=sqrt(qmv[i,'var']) 
				)
		}
		return(values_object)
	}

	sampleNorm2DQuantile <- function(values_subject1, values_subject2, 
		dist_subject1, dist_subject2, dist_object, na.rm=FALSE, quantiles=10)
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

		quantile2DMeanVar <- function(dist_subject1, dist_subject2, dist_object, 
			na.rm=FALSE, quantiles=100)
		# for each quantile of the distribution SUBJECT1 and SUBJECT2 gives
		# back the mean and the standard deviation of distribution OBJECT. 
		# Returns the two square matrices of mean and variance corresponding 
		# to each pair of quantiles of SUBJECT1 and SUBJECT2.
		{
			idx1 <- .which.quantile(dist_subject1, na.rm=na.rm, 
				quantiles=quantiles)
			idx2 <- .which.quantile(dist_subject2, na.rm=na.rm, 
				quantiles=quantiles)
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

	generateParams <- function(tpts, sampled_val, log2foldchange, 
		probs=c(constant=.5,sigmoid=.3,impulse=.2))
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
		params[constant_idx] <- lapply(sampled_val[constant_idx], 
			function(val) 
				list(type='constant', fun=.constantModelP , params=val, df=1)
				)
		
		# impulse varying, fist guess
		impulse_idx <- 1:nGenes %in% sample(which(!constant_idx), n_impulse)
		impulseParamGuess <- lapply(which(impulse_idx), 
			function(i)
				generateImpulseParams(tpts, sampled_val[i], log2foldchange[i])
				)
		valuesImpulse <- do.call('rbind', lapply(impulseParamGuess, 
			function(par) .impulseModel(tpts,par)
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
				, fun=.impulseModelP
				, params=generateImpulseParams(
					tpts 
					, sampled_val[i]
					, log2foldchange[i] * factor_of_correction
					)
				, df=6
				)
			)
		
		# sigmoid
		sigmoid_idx <- !constant_idx & !impulse_idx
		params[sigmoid_idx] <- lapply(
			which(sigmoid_idx)
			, function(i) list(
				type='sigmoid'
				, fun=.sigmoidModelP
				, params=generateSigmoidParams(
					tpts
					, sampled_val[i] 
					, log2foldchange[i] 
					)
				, df=4
				)			
			)

		# # report true foldchanges
		# simulatedFC <- apply(values, 1, function(x) diff(log2(range(x))))

		return(params)

	}

	noiseEval <- function(sim, real, plotout = FALSE)
	# plotout implemented in a previous version
	{
		quantiles <- 10
		sim.bkp <- sim
		# calculate mean and variance for genes of the real data
		rdMeans <- rowMeans(real, na.rm=TRUE)
		rdVars  <- .rowVars(real, na.rm=TRUE)
		## remove NA
		ix <- !is.na(rdMeans) & !is.na(rdVars)
		rdMeans <- rdMeans[ix]
		rdVars <- rdVars[ix]
		##
		rdQt <- .which.quantile(rdMeans, quantiles = quantiles)
		qtCenters <- quantile(rdMeans, 
			probs=seq(1/quantiles/2, 1, by=1/quantiles))
		# calculate mean and variance for genes of the sim data
		sdMeans <- rowMeans(sim)
		sdVars  <- .rowVars(sim)
		sdQt <- .which.quantile(sdMeans, rdMeans, quantiles = quantiles)
		# assign a positive variace to the constant sim genes
		# that is lower than the other observed variances
		nullVar <- sdVars < 1e-9
		sdVars[nullVar] <- seq(0,min(sdVars[!nullVar]),
			length.out=length(which(nullVar)))
		# identify genes whose variance before the addition of noise
		# is above the 9th decile of the real varince of genes within 
		# the same quantile
		thresholds <- tapply(rdVars, rdQt, quantile, probs=.9)
		sdOverVar  <- sdVars > thresholds[sdQt]
		# exclude those genes
		sim <- sim[!sdOverVar,]
		sdMeans <- sdMeans[!sdOverVar]
		sdVars  <- sdVars[!sdOverVar]
		sdQt <- sdQt[!sdOverVar] 
		# .which.quantile(sdMeans, rdMeans, quantiles = 10)
		#
		newVar <- sapply(
			1:length(sdVars) 
			, function(i) 
				quantile(rdVars[rdQt==sdQt[i]], 
					probs=ecdf(sdVars[sdQt==sdQt[i]])(sdVars[i]))
				)
		# noiseVar <- newVar
		# ??????? why the noise variance is not the difference between
		# the new variance evaluated and the old variance?
		noiseVar <- newVar-sdVars
		noiseVar[noiseVar <= 0] <- NA
		#
		outNoise <- rep(NA,nrow(sim.bkp))
		outNoise[!sdOverVar] <- noiseVar
		return(outNoise)
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
	if( !is.null(newTpts) ) tpts <- newTpts
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

	# # time of max response
	# # alpha
	# alphaLog2FC <- alphaLog2FC[apply(alphaLog2FC, 1, function(x) any(!is.na(x))),]
	# alphaabsfc <- abs(alphaLog2FC)
	# alphatmaxresponse <- apply(alphaabsfc, 1, which.max)
	# alphatmaxresponsePDF <- table(alphatmaxresponse)/length(alphatmaxresponse)
	# # beta
	# betaSimLog2FC <- betaSimLog2FC[apply(betaSimLog2FC, 1, function(x) any(!is.na(x))),]
	# betaabsfc <- abs(betaSimLog2FC)
	# betatmaxresponse <- apply(betaabsfc, 1, which.max)
	# betatmaxresponsePDF <- table(betatmaxresponse)/length(betatmaxresponse)
	# # gamma
	# gammaSimLog2FC <- gammaSimLog2FC[apply(gammaSimLog2FC, 1, function(x) any(!is.na(x))),]
	# gammaabsfc <- abs(gammaSimLog2FC)
	# gammatmaxresponse <- apply(gammaabsfc, 1, which.max)
	# gammatmaxresponsePDF <- table(gammatmaxresponse)/length(gammatmaxresponse)
	# possibly reduce variance prior to the addition of noise

	# if fc_scale > 1
	fc_scale <- 0
	if( fc_scale != 0)
	{
		alphaSimLog2FC <- sign(alphaSimLog2FC) * 
			(abs(alphaSimLog2FC) - fc_scale)
		betaSimLog2FC  <- sign(betaSimLog2FC)  * 
			(abs(betaSimLog2FC)  - fc_scale)
		gammaSimLog2FC <- sign(gammaSimLog2FC) * 
			(abs(gammaSimLog2FC) - fc_scale)
	}
	# transform time in log scale
	log_shift <- .find_tt_par(tpts)
	x <- .time_transf(tpts, log_shift)
	# generate alpha, beta, gamma
	message('generating rates time course...')
	alphaParams <- generateParams(x, alphaVals, alphaSimLog2FC, probs)
	betaParams  <- generateParams(x, betaVals , betaSimLog2FC , probs)
	gammaParams <- generateParams(x, gammaVals, gammaSimLog2FC, probs)
	# evaluate noise
	message('evaluating noise for simulated alpha, total and pre...')
	# generate total and preMRNA from alpha,beta,gamma
	paramSpecs <- lapply(1:nGenes, 
		function(i) 
			list(alpha=alphaParams[[i]], beta=betaParams[[i]], 
				gamma=gammaParams[[i]]))
	out <- lapply(1:nGenes, function(i) 
		.makeModel(tpts, paramSpecs[[i]], log_shift, 
			.time_transf, ode, .rxnrate))
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
	alphaSim_noisevar <- noiseEval(cleanDataSet$rates$alpha, alpha)
	totalSim_noisevar <- noiseEval(cleanDataSet$concentrations$total, total)
	preSim_noisevar <- noiseEval(cleanDataSet$concentrations$preMRNA, preMRNA)
	# select genes whose noise evaluation succeded
	okGenes <- which(
		!is.na(alphaSim_noisevar) & 
		!is.na(totalSim_noisevar) & 
		!is.na(preSim_noisevar) 
		)
	nGenes <- nGenes.bkp
	# select randomly among the okGenes the genes of the sample
	if( length(okGenes) >= nGenes ) 
		okGenes <- sample(okGenes, nGenes)
	else warning(paste('makeSimData: Only',length(which(okGenes)), 
		'genes were be generated, instead of',nGenes) )
	# keep okGenes only
	paramSpecs <- paramSpecs[okGenes]
	alphaSim_noisevar <- alphaSim_noisevar[okGenes]
	totalSim_noisevar <- totalSim_noisevar[okGenes]
	preSim_noisevar   <- preSim_noisevar[okGenes]
	# cleanDataSet <- subsetSnOut(cleanDataSet, okGenes)
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
		))
}

.bestModel <- function(object, bTsh=NULL, cTsh=NULL) {
		## in case bTsh or bTsh are provided set them as
		# permanent for the object
		if( is.null(bTsh) )
			bTsh <- object@params$thresholds$brown
		if( is.null(cTsh) )
			cTsh <- object@params$thresholds$chisquare
		## calculate ratePvals
		ratePvals <- ratePvals(object, cTsh)
		## give a discrete classification per each rate per each gene
		# according to the brown's threshold for the pvalues
		acceptedVarModels <- sapply(1:3, function(i) ratePvals[,i]<bTsh[i])
		if( !is.matrix(acceptedVarModels) )
			acceptedVarModels <- t(as.matrix(acceptedVarModels))
		# nonResolvedGenes <- apply(acceptedVarModels, 1, 
		# 	function(x) all(is.na(x)))
		acceptedVarModels[is.na(acceptedVarModels)] <- FALSE
		rownames(acceptedVarModels) <- rownames(ratePvals)
		colnames(acceptedVarModels) <- colnames(ratePvals)
		geneClass <- apply(acceptedVarModels, 1, 
			function(accepted) paste(c('a','b','c')[accepted],collapse=''))
		geneClass[geneClass==''] <- '0'
		# geneClass[nonResolvedGenes] <- NA
		## retrive all the models
		ratesSpecs <- object@ratesSpecs
		## select the best model (according to geneClass) per gene
		nGenes <- length(ratesSpecs)
		object@ratesSpecs <- lapply(1:nGenes, 
			function(i) ratesSpecs[[i]][geneClass[i]])
		return(object)
	}











