### minimization of the b,c,bc models but 
### 100 times less than other models

.time_transf <- function(t, log_shift, c = NaN) 
{
	newtime <- log2(t+log_shift)
	return(newtime)
} 

.time_transf_No4sU <- function(t, log_shift, c) 
{
	newtime <- log2(t+log_shift) + c
	return(newtime)
} 

.D2impulseModel <- function(t, par) {
  h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
  -(2*b^2*(h1-h0)*(h1-h2)*exp(b*(t-t2)-b*(t-t1)))/(h1*(exp(-b*(t-t1))+1)^2*(exp(b*(t-t2))+1)^2)+((h1-h2)*((2*b^2*exp(2*b*(t-t2)))/(exp(b*(t-t2))+1)^3-(b^2*exp(b*(t-t2)))/(exp(b*(t-t2))+1)^2)*((h1-h0)/(exp(-b*(t-t1))+1)+h0))/h1+((h1-h0)*((2*b^2*exp(-2*b*(t-t1)))/(exp(-b*(t-t1))+1)^3-(b^2*exp(-b*(t-t1)))/(exp(-b*(t-t1))+1)^2)*((h1-h2)/(exp(b*(t-t2))+1)+h2))/h1
}

chisq <- function(experiment, model, variance=NULL)
{
	sum((experiment - model )^2/variance )
}

.getRatesAndConcentrationsFromRpkms <- function(totRpkms, labeledRpkms, tpts
, tL=NULL, degDuringPulse=FALSE, simulatedData=FALSE, BPPARAM=bpparam() #parallelize=TRUE
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
	
	## eventually set the 'simple' mode
	if( is.null(totRpkms$introns) | is.null(labeledRpkms$introns) ){onlyExonsMode <- TRUE}else{onlyExonsMode <- FALSE}
	
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
	
	##############
	# NORMALIZE THE DATA  - only needed with real data, synthetic data are 
	# generated already scaled
	#########################
	if( realData ){
	
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
					Texo_var <- t(t(Texo_var)*totalSF^2)
				if( !onlyExonsMode ) {
					Tint <- t(t(Tint)*totalSF)
					if( totalVariance )
						Tini_var <- t(t(Tint_var)*totalSF^2)
				}
		}
		
		if( !totalVariance ) {
				## evaluate the variance now
				Texo_var <- apply( Texo , 1 , speedyVar )
				if( !onlyExonsMode )
					Tint_var <- apply( Tint , 1 , speedyVar )
		}
		
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
			t(t(Lexo_var)*labeledSF^2)
			#Lexo_var <- apply(t(t(Lexo_var)*labeledSF^2), 1, mean, na.rm=TRUE)
			else
			Lexo_var <- apply(Lexo, 1, speedyVar)
			} else {
			# only estimate variance
			if( labeledVarince )
			Lexo_var <- Lexo_var
			#Lexo_var <- apply(Lexo_var, 1, mean, na.rm=TRUE)
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
				}
			))
	
			if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
			if( steadyStateMode == 0 ) TintDer[, 1] <- 0 
		
			TexoDer <- as.matrix(sapply(1:nrow(Texo), 
				function(i) {
					if( all(is.finite(Texo[i,] ) ) ) {
						spfun <- splinefun(tpts, Texo[i,])
						return(spfun(tpts, deriv=1) )
					} else return(rep(NA, length(tpts)) )
				}
			))
	
			if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
			if( steadyStateMode == 0 ) TexoDer[, 1] <- 0
		
			if( labeledMedianNorm )
				warning('When introns of labelled and total fraction are provided normalization of 
				labelled fraction is computed internally. "labeledMedianNorm" argument is ignored.')
		
			if( !is.null(labeledSF) ){
				warning('When introns of labelled and total fraction are provided normalization of 
				labelled fraction could be computed internally.')
				labeledSF_prior <- labeledSF
			}else{
		
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
					errorfun <- function(gamma, Lexo, Lint, tL ) (Lint-fLint(Lexo, gamma, tL))^2
					
					optimize(errorfun, c(0,maxGamma), Lexo=Lexo, Lint=Lint, tL=tL)$minimum
				}
				# calculate the factor which bring the median of the residuals between
				# the modeled preMRNA levels and the measured to zero
				sq.median.resids <- function(sf, P, dP, alpha, gamma) sapply(sf, function(i) {
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
		
				gammaTC <- do.call('cbind',bplapply(1:length(tpts), function(j) 
					sapply(1:nrow(Lint), function(i, Lint, Lexo, tL) 
						fGamma(Lint[i,j] , Lexo[i,j] , tL)
					, Lint=Lint, Lexo=Lexo, tL=tL)
				,BPPARAM=BPPARAM))
		
				# scale factor 
				labeledSF_prior <- sapply(1:ncol(Tint), function(j)
					optimize(sq.median.resids, c(0.01,100), P=Tint[,j], dP=TintDer[,j], 
					alpha=Lexo[,j]/tL, gamma=gammaTC[,j] )$minimum)
			}
		
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
					c(Lexo[i,j]/tL
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
							mrOut <- tryCatch(
								multiroot(sys4suScale, c(init,labeledSF_prior[j]))
							
							, error=function(e)
							list(root=rep(NA, 4), estim.precis=NA)
							)
							labeledSf[i,j] <- mrOut$root[4]
							labeledSfep[i,j] <- mrOut$estim.precis
						}
					}
				}))
		
				## chose the best resolved genes to estimate 
				## from them the scale factor
				epTsh <- apply(labeledSfep,2,quantile,probs=.75,na.rm=TRUE)
				ix <- t(apply(labeledSfep,1, function(x) x>epTsh))
				labeledSf[ix] <- NA
				labeledSF <- apply(labeledSf, 2, stats::median, na.rm=TRUE)
		
			} else {
		
				labeledSF <- labeledSF_prior
		
			}
		
		}
		
		## simulated data
	} else {
	
		# in case of synthetic data the time course is already scaled
		# therefore just rename the variables and compute varince in 
		# case is not provided
		alphaTC <- Lexo
		if( labeledVarince )
			alphaTC_var <- Lexo_var ################################## M.F. varianza alpha
		else {
			# calculate variance from the time course and replicate it
			# for each time point
			alphaTC_var <- apply(alphaTC, 1, speedyVar)
			alphaTC_var <- sapply(1:ncol(alphaTC), function(i) alphaTC_var)
		}
		if( totalVariance ) {
			Texo_var <- Texo_var
			Tint_var <- Tint_var
		} else {
			Texo_var <- apply( Texo , 1 , speedyVar )
			Texo_var <- sapply(1:ncol(Texo), function(i) Texo_var)
			Tint_var <- apply( Tint , 1 , speedyVar )
			Tint_var <- sapply(1:ncol(Tint), function(i) Tint_var)
		}
		# scaling factor for synthetic dataset is meaningless
		totalSF <- labeledSF <- rep(1, length(tpts))
	
		# derivatives from time course
		TintDer <- as.matrix(sapply(1:nrow(Tint), 
			function(i) {
				if( all(is.finite(Tint[i,] ) ) ) {
					spfun <- splinefun(tpts, Tint[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
		if( steadyStateMode == 0 ) TintDer[, 1] <- 0 
	
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
				spfun <- splinefun(tpts, Texo[i,])
				return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
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
						# labeledSfep[i,j] <- NA
						# alphaTC[i,j] <- NA
						# betaTC[i,j] <- NA
						# gammaTC[i,j] <- NA
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
			
			## recalculate the variance
			if( labeledVarince )
				alphaTC_var <- Lexo_var/tL^2
			else {
				alphaTC_var <- apply(alphaTC, 1, speedyVar)
				alphaTC_var <- sapply(1:ncol(alphaTC), function(i) alphaTC_var)
			}
			
			## set preMRNA and gamma to NA
			Tint <- matrix(NA, nrow(Texo), ncol(Texo))
			Tint_var <- matrix(NA, nrow(Texo), ncol(Texo))
			gammaTC <- matrix(NA, nrow(betaTC), ncol(betaTC))
		
		# assume that no degradation occur during pulse
		} else {
		
			## calculate alpha and recalculate the variance
			alphaTC <- Lexo/tL
			if( labeledVarince )
				alphaTC_var <- Lexo_var/tL^2
			else {
				alphaTC_var <- apply(alphaTC, 1, speedyVar)
				alphaTC_var <- sapply(1:ncol(alphaTC), function(i) alphaTC_var)
			}
			
			## calculate beta
			message('Estimating degradation rates...')
			inferKBetaFromIntegral <- function(tpts, alpha, totalRNA, maxBeta=75) 
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
							, c(1e-5, maxBeta)
							, t0 = tpts[j-1]
							, t1 = tpts[j]
							, alpha_t0 = alpha[i,j-1]
							, alpha_t1 = alpha[i,j]
							, X_t0 = totalRNA[i,j-1]
							, X_t1 = totalRNA[i,j]
							)
						, error=function(e) return(list(root=NA, estim.prec=NA, error=e))
					)})
				, BPPARAM=BPPARAM)
			}
			if( steadyStateMode == 1 ) TexoDer[,1] <- 0
			betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / Texo[,1]
			betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
			if( length(tpts)>1 ) {
				betaOut <- inferKBetaFromIntegral(tpts, alphaTC, Texo, 
					maxBeta=quantile(betaT0, na.rm=TRUE, probs=.99)*10
					)
				betaTC <- cbind(betaT0, 
					sapply(betaOut, function(x) sapply(x, '[[', 'root'))
					)
				ratesEstimPrec <- cbind(0,
					sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
					)
			} else {
				betaTC <- as.matrix(betaT0)
				ratesEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
			}

			# ## impute NA values
			betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
				function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=bpparam()))

			## set preMRNA and gamma to NA
			Tint <- matrix(NA, nrow(Texo), ncol(Texo))
			Tint_var <- matrix(NA, nrow(Texo), ncol(Texo))
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
				# labeledSfep[i,j] <- NA
				# alphaTC[i,j] <- NA
				# betaTC[i,j] <- NA
				# gammaTC[i,j] <- NA
				# }
			}
		}
		}))
		
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
		
		## recalculate the variance
		if( labeledVarince )
			alphaTC_var <- Lexo_var*(labeledSF/tL)^2
		else {
			alphaTC_var <- apply(alphaTC, 1, speedyVar)
			alphaTC_var <- sapply(1:ncol(alphaTC), function(i) alphaTC_var)
		}
		
		# assume that no degradation occur during pulse
		} else {
		
		# calculate alpha and recalculate the variance
		alphaTC <- Lexo/tL
		if( labeledVarince )
			alphaTC_var <- Lexo_var/tL^2
		else {
			alphaTC_var <- apply(alphaTC, 1, speedyVar)
			alphaTC_var <- sapply(1:ncol(alphaTC), function(i) alphaTC_var)
		}

		alphaTC <- t(t(alphaTC)*labeledSF)
		
		if( labeledVarince )
			alphaTC_var <- t(t(alphaTC_var)*labeledSF^2)
		else {
			alphaTC_var <- apply(alphaTC, 1, speedyVar)
			alphaTC_var <- sapply(1:ncol(alphaTC), function(i) alphaTC_var)
		}
		
		# calculate beta
		
		message('Estimating degradation rates...')
		if( steadyStateMode == 1 ) TexoDer[,1] <- 0
		betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / (Texo[,1] - Tint[,1] )
		betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
		if( length(tpts)>1 ) {
			betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, Texo, Tint, 
				maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
				)
			betaTC <- cbind(betaT0, 
				sapply(betaOut, function(x) sapply(x, '[[', 'root'))
				)
			betaEstimPrec <- cbind(0,
				sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
				)
		} else {
			betaTC <- as.matrix(betaT0)
			betaEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
		}

		###################################
		## estimate processing rates #########
		#####################################
		
		# calculate gamma (from  total RNA introns and alphas )
		message('Estimating processing rates...')
		if( steadyStateMode == 1 ) TintDer[,1] <- 0
		gammaT0 <- ( alphaTC[,1] - TintDer[,1] ) / Tint[,1]
		gammaT0[gammaT0 < 0 | !is.finite(gammaT0)] <- NA
		if( length(tpts)>1 ) {
			gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, Tint, 
				maxGamma=quantile(gammaT0,na.rm=TRUE,probs=.99)*10, BPPARAM=BPPARAM
				)
			gammaTC <- cbind(gammaT0, 
				sapply(gammaOut, function(x) sapply(x, '[[', 'root'))
				)
			gammaEstimPrec <- cbind(0, 
				sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec'))
				)
		} else {
			gammaTC <- as.matrix(gammaT0)
			gammaEstimPrec <- matrix(0, nrow=nrow(gammaTC), ncol=ncol(gammaTC))
		}

		# ## impute NA values
		betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
			function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=bpparam()))
		gammaTC <- do.call('rbind',bplapply(1:nrow(gammaTC), 
			function(i) impute_na_tc(tpts, gammaTC[i,]), BPPARAM=bpparam()))

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

	if(!is.matrix(Texo_var)){Texo_var <- matrix(rep(Texo_var,length(tpts)),nrow=length(Texo_var),ncol=length(tpts))}
	if(!is.matrix(Tint_var)){Tint_var <- matrix(rep(Tint_var,length(tpts)),nrow=length(Tint_var),ncol=length(tpts))}

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

.makeModel <- function(tpts, hyp, log_shift, .time_transf, ode, .rxnrate, c= NaN)
{
	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(.time_transf(x, log_shift, c), hyp$alpha$par)
	params$beta  <- function(x) 
		hyp$beta$fun$value(.time_transf(x, log_shift, c), hyp$beta$par)
	params$gamma <- function(x) 
		hyp$gamma$fun$value(.time_transf(x, log_shift, c), hyp$gamma$par)
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
.find_tt_par <- function(tpts)
{
  cvLogTpts <- function(a , tpts)
  {
  	newtime <- log2(tpts + a )
    stats::sd(diff(newtime)) / mean(diff(newtime))
  }
  if(length(tpts)>2){return(optimize(f=cvLogTpts, interval=c(0,5), tpts=tpts )$minimum)}
  else{return(1)}
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

	nGenes <- nrow(rates$alpha)
	tpts <- tpts
	paramSpecs <- bplapply(1:nGenes, modelOneGene, seed=seed, 
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
		.logLikelihood=.logLikelihood,
		BPPARAM=BPPARAM
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
	## otherwise calculate the brown's test for the available genes
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
			xout_matrix <- suppressWarnings(stats::cor(y_group, use='complete.obs'))
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
# NumericVector .impulseModelC(NumericVector x, NumericVector par) {
# 	int n = x.size();
# 	NumericVector ans(n);
# 	for(int i = 0; i < n; i++) {
# 		ans[i] = 1/par[1]*(par[0]+(par[1]-par[0])*(1/(1+exp(-par[5]*(x[i]-par[3])))))*(par[2]+(par[1]-par[2])*(1/(1+exp(par[5]*(x[i]-par[4])))));
# 	}
# 	return ans;
# }
# ')
# .impulseModel <- .impulseModelC

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
	type <- names(pvals)[which.min(pvals)]
	df   <- dfs[which.min(pvals)]

	if( type=='sigmoid'    ) params <- outSM[,bestSM]$par
	if( type=='impulse'    ) params <- outIM[,bestIM]$par
	if( type=='polynomial' ) params <- outPN[,bestPN]$par

	pval <- pvals[which.min(pvals)]
	fun  <- funcs[[which.min(pvals)]]

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
					list(type='constant', fun=.constantModelP , params=val, df=1)
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
		}
		
		# sigmoid
		sigmoid_idx <- !constant_idx & !impulse_idx
		if( any(sigmoid_idx) ) {
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
	log_shift <- .find_tt_par(tpts)
	x <- .time_transf(tpts, log_shift)
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
					.time_transf, ode, .rxnrate),error=function(e){cbind(time=rep(NaN,length(tpts))
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

.bestModel <- function(object, bTsh=NULL, cTsh=NULL) {

	preferPValue <- object@params$preferPValue
	
	## in case bTsh or bTsh are provided set them as
	# permanent for the object
	if( is.null(bTsh) )
		bTsh <- object@params$thresholds$brown
	if( is.null(cTsh) )
		cTsh <- object@params$thresholds$chisquare
	## calculate ratePvals
	ratePvals <- ratePvals(object, cTsh)
	ratePvals <- replace(ratePvals,is.na(ratePvals),1)

	if(preferPValue)
	{
		pValues <- chisqtest(object)
		rownames(pValues) <- rownames(ratePvals)

		geneClass <- sapply(rownames(ratePvals),function(i)
		{

			acceptableModelsTemp <- which(pValues[i,] <= cTsh)
			if(length(acceptableModelsTemp)==0){return("0")}
			if(length(acceptableModelsTemp)==1){return(colnames(pValues)[acceptableModelsTemp])}    

			nameTemp <- rep("K",3)
			if(ratePvals[i,"synthesis"] <= bTsh["synthesis"]){nameTemp[[1]] <- "V"}
			if(ratePvals[i,"processing"] <= bTsh["processing"]){nameTemp[[2]] <- "V"}
			if(ratePvals[i,"degradation"] <= bTsh["degradation"]){nameTemp[[3]] <- "V"}

			nameTemp <- paste0(nameTemp[[1]],nameTemp[[2]],nameTemp[[3]])
			nameTemp <- c("KKK" = "0"
						 ,"VKK" = "a"
						 ,"KVK" = "c"
						 ,"KKV" = "b"
						 ,"VVK" = "ac"
						 ,"VKV" = "ab"
						 ,"KVV" = "bc"
						 ,"VVV" = "abc")[nameTemp]

			if(pValues[i,nameTemp]<=cTsh|!is.finite(pValues[i,nameTemp])){return(nameTemp)}else{
				return(names(which.min(pValues[i,])))}
    	})

		ratesSpecs <- object@ratesSpecs
		## select the best model (according to geneClass) per gene
		nGenes <- length(ratesSpecs)
		object@ratesSpecs <- lapply(1:nGenes, 
			function(i) ratesSpecs[[i]][geneClass[i]])
		return(object)

	}else{

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
}


####################################################################################################k1KKK_No4sU <- function(x, par){par[1]*par[3]}

mcsapply <- function( X, FUN, ... ) do.call('cbind', bplapply( X, FUN, ... ))

lineCoefficients_No4sU <- function(xInitial
								  ,xFinal
								  ,yInitial
								  ,yFinal)
{
  return(c(m = (yFinal-yInitial)/(xFinal-xInitial)
          ,q = (yInitial*xFinal-yFinal*xInitial)/(xFinal-xInitial)))
}

firstStep_No4sU <- function(tpts
						   ,mature
						   ,premature
						   ,matureVariance
						   ,Dmin
						   ,Dmax) 
{
	fits <- lapply(1:nrow(mature), function(row)
  	{
	    optimize(firstStepError_No4sU
	            ,c(Dmin, Dmax)
	            ,k2K3Ratio = mature[row,1]/premature[row,1]
	            ,tpts = tpts
	            ,premature = premature[row,]
	            ,mature = mature[row,]
	            ,matureVariance = matureVariance[row,])
   })
   
	out <- t(sapply(fits,unlist))
	out[,2] <- pchisq(out[,2], length(tpts)-1)
	colnames(out) <- c('k3','p')
	return(out)
}

firstStepError_No4sU <- function(tpts
								,mature
								,premature
								,matureVariance
								,k3
								,k2K3Ratio)
{
	matureEstimated <- numeric(length = length(mature))
	matureEstimated[1] <- mature[1]

	for(t in 2:length(mature))
	{
		tInitial <- tpts[t-1]
		tFinal <- tpts[t]

		prematureInitial <- premature[t-1]
		prematureFinal <- premature[t]

		matureInitial <- matureEstimated[t-1]

		coefficients <- lineCoefficients_No4sU(xInitial = tInitial
		                                	  ,xFinal = tFinal
		                                	  ,yInitial = prematureInitial
		                                	  ,yFinal = prematureFinal)
		m <- coefficients[1]
		q <- coefficients[2]

    	matureEstimated[t] <- (exp(-k3*(tFinal-tInitial))*matureInitial
                        	+ exp(-k3*tFinal)*k2K3Ratio*(q*(exp(k3*tFinal)-exp(k3*tInitial))
                        	+ m*(exp(k3*tFinal)*(k3*tFinal-1)/k3
                        	- exp(k3*tInitial)*(k3*tInitial-1)/k3)))
  	}

   return(sum((mature[-1] - matureEstimated[-1])^2/matureVariance[-1]))
}

secondStepError_No4sU <- function(tpts
								 ,mature
								 ,premature
								 ,matureVariance
								 ,k2k3)
{
	k2 <- k2k3[1]
	k3 <- k2k3[2]
	
	matureEstimated <- numeric(length = length(mature))
	matureEstimated[1] <- mature[1]
  
	for(t in 2:length(mature))
	{
		tInitial <- tpts[t-1]
		tFinal <- tpts[t]

		prematureInitial <- premature[t-1]
		prematureFinal <- premature[t]

		matureInitial <- matureEstimated[t-1]

		coefficients <- lineCoefficients_No4sU(xInitial = tInitial
											  ,xFinal = tFinal
											  ,yInitial = prematureInitial
											  ,yFinal = prematureFinal)
		m <- coefficients[1]
		q <- coefficients[2]

		matureEstimated[t] <- (exp(-k3*(tFinal-tInitial))*matureInitial
							 + exp(-k3*tFinal)*k2/k3*(q*(exp(k3*tFinal)-exp(k3*tInitial))
							 + m*(exp(k3*tFinal)*(k3*tFinal-1)/k3
							 - exp(k3*tInitial)*(k3*tInitial-1)/k3)))
	}

	return(mean((mature[-1] - matureEstimated[-1])^2/matureVariance[-1]))
}

.getRatesAndConcentrationsFromRpkms_No4sU <- function(totRpkms
													, tpts
													, BPPARAM=bpparam()
													, modellingParameters=list(Dmin = 1e-6
																			 , Dmax = 10)
													, genesFilter
													)
{

	Dmin <- modellingParameters$Dmin
	Dmax <- modellingParameters$Dmax

	eiGenes <- rownames(totRpkms$exons)

	# Mature, premature and total rpkms
	mature <- totRpkms$exons - totRpkms$introns
	premature <- totRpkms$introns
	total <- totRpkms$exons

	# Mature, premature and total variance, rpkms
	matureVariance <- totRpkms$introns_var+totRpkms$exons_var
	prematureVariance <- totRpkms$introns_var
	totalVariance <- totRpkms$exons_var

	boolExpressionsAndVariances <- apply(mature,1,function(row)any(row < 0)) |
	                               apply(premature,1,function(row)any(row < 0)) |
	                               apply(matureVariance,1,function(row)any(row < 0)) |
	                               apply(prematureVariance,1,function(row)any(row < 0))

	boolExpressionsAndVariances[!is.finite(boolExpressionsAndVariances)] <- TRUE

	if(any(boolExpressionsAndVariances))
		{print("newINSPEcT: for some genes mRNA or preMRNA expressions or variances are negative. Those genes are excluded from the analysis.")}

	eiGenes <- eiGenes[!boolExpressionsAndVariances]

	mature <- mature[eiGenes,]
	premature <- premature[eiGenes,]
	matureVariance <- matureVariance[eiGenes,]
	prematureVariance <- prematureVariance[eiGenes,]
	total <- total[eiGenes,]
	totalVariance <- totalVariance[eiGenes,]

	# Constant post transcriptional rates and fixed post transcriptional ratio 
	k3Prior <- firstStep_No4sU(tpts = tpts
							  ,mature = mature
							  ,premature = premature
							  ,matureVariance = matureVariance
							  ,Dmin = Dmin
							  ,Dmax = Dmax)
	rownames(k3Prior) <- eiGenes

	# Constant post transcriptional rates and variable post transcriptiona ratio
	fits <- t(mcsapply(1:nrow(mature), function(row)
	{
		unlist(
			tryCatch(
	    		optim(par = c(mature[row,1]/premature[row,1]*k3Prior[row,'k3'], k3Prior[row,'k3'])
	    				   ,fn = secondStepError_No4sU
	        			   ,tpts = tpts
	        			   ,premature = premature[row,]
	        			   ,mature = mature[row,]
	        			   ,matureVariance = matureVariance[row,])
			,error=function(e)list(par = c(NaN,NaN), value = NaN, counts = NaN, convergence = NaN, message = "Optimization error."))[1:4])
	},BPPARAM = BPPARAM))
	
	fits[,3] <- pchisq(fits[,3], length(tpts)-3)
	colnames(fits) <- c('k2','k3','p','counts','gradient','convergence')
	rownames(fits) <- eiGenes

	# Correction of negative priors with the median
	if(genesFilter){
		fits[fits[,'k2']<0,'k2'] <- NaN
		fits[fits[,'k3']<0,'k3'] <- NaN

		notFiniteRates <- !is.finite(fits[,'k2']) | !is.finite(fits[,'k3'])

		fits[notFiniteRates,'k2'] <- median(fits[is.finite(fits[,'k2']),'k2'])
		fits[notFiniteRates,'k3'] <- median(fits[is.finite(fits[,'k3']),'k3'])

		fits[notFiniteRates,'p'] <- NaN
	}

	# Data formatting
	constantModels <- list(models = fits
						 , premature = premature
						 , mature = mature
						 , prematureVariance = prematureVariance
						 , matureVariance = matureVariance)

	ratesConstantPriors <- constantModels$models  

	# alphaTC <- t(sapply(seq_along(ratesConstantPriors[,'k3']),function(g)
	# {
	# 	sapply(tpts,function(t){k1KKK_No4sU(t,par = c(mean(mature[g,],na.rm = T),ratesConstantPriors[g,'k2'],ratesConstantPriors[g,'k3']))})
	# }))

	# betaTC <- matrix(rep(ratesConstantPriors[,'k3'],length(tpts)),ncol=length(tpts))
	gammaTC <- matrix(rep(ratesConstantPriors[,'k2'],length(tpts)),ncol=length(tpts))

	prematureDer <- as.matrix(t(sapply(1:nrow(premature),function(i){
		if(all(is.finite(premature[i,]))){
			spfun <- splinefun(tpts, premature[i,])
			return(spfun(tpts, deriv=1))
		} else return(rep(NA,length(tpts)))
	})))
	prematureDer[,1] <- 0 # force the steady state at time 0

	alphaTC <- prematureDer + gammaTC * premature
	alphaTC[alphaTC<0] <- NaN

	#Evaluate beta as constant between intervals

	betaT0 <- alphaTC[,1] / mature[,1]
	betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
	
	betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, total, premature, 
				maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
				)
	betaTC <- cbind(betaT0, 
		sapply(betaOut, function(x) sapply(x, '[[', 'root'))
	)
	betaEstimPrec <- cbind(0,
		sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
	)

<<<<<<< HEAD
	#Evaluate gamma as constant between intervals

	gammaT0 <- alphaTC[,1] / premature[,1]
	gammaT0[gammaT0 < 0 | !is.finite(gammaT0)] <- NA

	gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, premature, 
		maxGamma=quantile(gammaT0,na.rm=TRUE,probs=.99)*10, BPPARAM=BPPARAM
		)
	gammaTC <- cbind(gammaT0, 
		sapply(gammaOut, function(x) sapply(x, '[[', 'root'))
		)
	gammaEstimPrec <- cbind(0, 
		sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec'))
		)

	print('imputing NA values')
	# ## impute NA values
	alphaTC <- do.call('rbind',bplapply(1:nrow(alphaTC), 
		function(i) impute_na_tc(tpts, alphaTC[i,]), BPPARAM=bpparam()))
	betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
		function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=bpparam()))
	gammaTC <- do.call('rbind',bplapply(1:nrow(gammaTC), 
		function(i) impute_na_tc(tpts, gammaTC[i,]), BPPARAM=bpparam()))

	## caluculate error through integration of alphaTC, betaTC, gammaTC?

=======
>>>>>>> f7d7f2c07ed09f1c3246f3f5954906b07f324835
	pModel <- fits[,"p"]
	pModel[apply(alphaTC,1,function(row)any(!is.finite(row)))] <- NaN

	alphaTC_var <- rep(1, length(eiGenes))
	ratesEstimPrec <- betaEstimPrec + gammaEstimPrec


	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(total, 'dimnames') <- NULL
	attr(premature, 'dimnames') <- NULL
	attr(alphaTC_var, 'names') <- NULL
	attr(totalVariance, 'names') <- NULL
	attr(prematureVariance, 'names') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL

	return(list(concentrations=list(total=total
						  		  , total_var=totalVariance
						  		  , preMRNA=premature
						  		  , preMRNA_var=prematureVariance)
			  , rates=list(alpha=alphaTC
						 , alpha_var=alphaTC_var
						 , beta=betaTC
						 , gamma=gammaTC)
			  , ratesFirstGuessP = pModel
	  		  , ratesEstimPrec=ratesEstimPrec
			  , geneNames=eiGenes
			  , tpts=tpts
			  , labeledSF=NULL
			  , totalSF=NULL
			  , tL = NULL
	))

}

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
		model <- .impulseModel(tpts,par)
		if( abs(par[6]) > Inf ) return(NaN)
		if( any(model < 0) ) return(NaN)
		chisq(experiment, model, variance)
	}

	im_chisq <- function(par, tpts, experiment, variance=NULL, tt_c) 
	{
		model <- .impulseModel(tpts,par)
		if( any(model < 0) ) return(NaN)
		chisq(experiment, model, variance)
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

prematureKKK_Int_No4sU <- function(x, parameters)
{
  matureParameters <- parameters[1]
  k2Parameters <- parameters[2]
  k3Parameters <- parameters[3]
  
  return((k3Parameters*matureParameters)/k2Parameters)
}

k1KKK_Int_No4sU <- function(x, par)
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

errorKKK_Int_No4sU <- function(parameters, tpts, premature, mature, prematureVariance, matureVariance)
{

  if(parameters[1]<0)return(NaN)
  if(parameters[2]<0)return(NaN)
  if(parameters[3]<0)return(NaN)

  matureParameters <- parameters[1]

  prematureEstimated <- prematureKKK_Int_No4sU(x = tpts, parameters = parameters)
  matureEstimated <- matureParameters

  prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
  matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

  return(sum(c(prematureChiSquare,matureChiSquare)))
}

errorVKK_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) return(parameters[8])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  chi2 <- chisq(data,modData,datavar)

  return(chi2)
}

errorVKV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisq(data,modData,datavar)

  return(chi2)
}

errorVVV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[7:12])}
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[13:18])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisq(data,modData,datavar)

  return(chi2)
}

errorKVV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[2:7])}
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  chi2 <- chisq(data,modData,datavar)

  return(chi2)
}

errorVVK_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

	k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
	k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[7:12])}
	k3F <- function(x) return(parameters[13])

	if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

	modData <- systemSolution(k1F,k2F,k3F,times)
	
	chi2 <- chisq(data,modData,datavar)

	return(chi2)
}

errorKVK_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{
	k1F <- function(x)return(parameters[1])
	k2F <- function(x){.impulseModel(log2(x+a)+c,parameters[2:7])}
	k3F <- function(x)return(parameters[8])
	
	if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)
	
	modData <- systemSolution(k1F,k2F,k3F,times)
	
	chi2 <- chisq(data,modData,datavar)
	
	return(chi2)
}

errorKKV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) return(parameters[2])
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[3:8])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  chi2 <- chisq(data,modData,datavar)

  return(chi2)
}

logLikelihoodFunction <- function(experiment, model, variance=NULL)
{
    if( is.null(variance)) variance <- stats::var(experiment)
    sum(log(2*pnorm(-abs(experiment-model),mean=0,sd=sqrt(variance))))
}

loglikKKK_Int_No4sU <- function(parameters
	                    	   ,tpts
	                    	   ,premature
	                    	   ,mature
	                    	   ,prematureVariance
	                    	   ,matureVariance)
{

	matureParameters <- parameters[1]

	prematureEstimated <- prematureKKK_Int_No4sU(x = tpts, parameters = parameters)
	matureEstimated <- matureParameters

	logLikelihoodFunction(premature, prematureEstimated, prematureVariance) + 
	logLikelihoodFunction(mature, matureEstimated, matureVariance)

}

loglikVKK_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) return(parameters[8])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)

}

loglikKVK_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[2:7])}
  k3F <- function(x) return(parameters[8])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)
  # modelLm <- sapply(times,function(t)k1F(t)/k3F(t)*(1 - exp(-k3F(t)*1/6))+k1F(t)/(k3F(t) - k2F(t))*(exp(-k3F(t)*1/6) - exp(-k2F(t)*1/6)))

  # modData <- c(modData, modelLm)

  logLikelihoodFunction(data, modData, datavar)
}

loglikKKV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) return(parameters[2])
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[3:8])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVVK_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[7:12])}
  k3F <- function(x) return(parameters[13])

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVKV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) return(parameters[7])
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikKVV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) return(parameters[1])
  k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[2:7])}
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[8:13])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}

loglikVVV_Int_No4sU <- function(parameters, times, data, datavar, a, c)
{

  k1F <- function(x) {.impulseModel(log2(x+a)+c,parameters[1:6])}
  k2F <- function(x) {.impulseModel(log2(x+a)+c,parameters[7:12])}
  k3F <- function(x) {.impulseModel(log2(x+a)+c,parameters[13:18])}

  if( any(c(k1F(times),k2F(times),k3F(times))<0) ) return(NaN)

  modData <- systemSolution(k1F,k2F,k3F,times)

  logLikelihoodFunction(data, modData, datavar)
}



###########################################################################

k1VKK_Der_No4sU <- function(x, par, c)
{
	t_fact <- 2^(x-c)*log(2)
	.D2impulseModel(x, par[1:6])*t_fact^2/par[7] + .DimpulseModel(x, par[1:6])*(1+par[8]/par[7])*t_fact + par[8]*.impulseModel(x, par[1:6])
}

prematureVKK_Der_No4sU <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7]
	k3Parameters <- parameters[8]

	t_fact <- 2^(x-c)*log(2)
	(.DimpulseModel(x, matureParameters)*t_fact + k3Parameters * .impulseModel(x, matureParameters))/k2Parameters

}

errorVKK_Der_No4sU <- function(parameters
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
	k1 <- k1VKK_Der_No4sU(tpts,parameters, c)

	prematureEstimated <- prematureVKK_Der_No4sU(x = tpts, parameters = parameters, c = c)
	matureEstimated <- .impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

k1VVK_Der_No4sU <- function(x, par, c)
{
	t_fact <- 2^(x-c)*log(2)
	.D2impulseModel(x, par[1:6])/.impulseModel(x, par[7:12])*t_fact^2 +
	.DimpulseModel(x, par[1:6])*t_fact*(1 - log(2)*.DimpulseModel(x, par[7:12])/.impulseModel(x, par[7:12])^2 + par[13]/.impulseModel(x, par[7:12])) + 
	log(2)*.impulseModel(x, par[1:6])*(par[13]/log(2) - (par[13]*.DimpulseModel(x, par[7:12]))/.impulseModel(x, par[7:12])^2 )
}

prematureVVK_Der_No4sU <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7:12]
	k3Parameters <- parameters[13]
	
	t_fact <- 2^(x-c)*log(2)
	return((.DimpulseModel(x, matureParameters)*t_fact
	      + k3Parameters * .impulseModel(x, matureParameters))/.impulseModel(x, k2Parameters))
}

errorVVK_Der_No4sU <- function(parameters
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
	k1 <- k1VVK_Der_No4sU(tpts,parameters, c)

	prematureEstimated <- prematureVVK_Der_No4sU(x = tpts, parameters = parameters, c = c)
	matureEstimated <- .impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

k1VKV_Der_No4sU <- function(x, par, c)
{
	t_fact <- 2^(x-c)*log(2)
	.D2impulseModel(x, par[1:6])/par[7]*t_fact^2 +
	.DimpulseModel(x, par[1:6])*t_fact*(1 + .impulseModel(x, par[8:13])/par[7]) + 
	log(2)*.impulseModel(x, par[1:6])*( .DimpulseModel(x, par[8:13])/par[7] + .impulseModel(x, par[8:13])/log(2) )
}

prematureVKV_Der_No4sU <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7]
	k3Parameters <- parameters[8:13]
	
	t_fact <- 2^(x-c)*log(2)
	(.DimpulseModel(x, matureParameters)*t_fact + .impulseModel(x, k3Parameters) * .impulseModel(x, matureParameters))/k2Parameters
}

errorVKV_Der_No4sU <- function(parameters
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
	k1 <- k1VKV_Der_No4sU(tpts,parameters, c)

	prematureEstimated <- prematureVKV_Der_No4sU(x = tpts, parameters = parameters, c = c)
	matureEstimated <- .impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

k1VVV_Der_No4sU <- function(x, par, c)
{
  t_fact <- 2^(x-c)*log(2)
  .D2impulseModel(x, par[1:6])/.impulseModel(x, par[7:12])*t_fact^2 +
  .DimpulseModel(x, par[1:6])*t_fact*(1 - log(2)*.DimpulseModel(x, par[7:12])/.impulseModel(x, par[7:12])^2 + .impulseModel(x, par[13:18])/.impulseModel(x, par[7:12])) + 
  log(2)*.impulseModel(x, par[1:6])*(.DimpulseModel(x, par[13:18])/.impulseModel(x, par[7:12]) + .impulseModel(x, par[13:18])/log(2) - (.impulseModel(x, par[13:18])*.DimpulseModel(x, par[7:12]))/.impulseModel(x, par[7:12])^2 )
}

prematureVVV_Der_No4sU <- function(x, parameters, c)
{
	matureParameters <- parameters[1:6]
	k2Parameters <- parameters[7:12]
	k3Parameters <- parameters[13:18]

	t_fact <- 2^(x-c)*log(2)
	return((.DimpulseModel(x, matureParameters)*t_fact
		+ .impulseModel(x, k3Parameters)*.impulseModel(x, matureParameters))/.impulseModel(x, k2Parameters))
}

errorVVV_Der_No4sU <- function(parameters
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
	k1 <- k1VVV_Der_No4sU(tpts,parameters, c)

	prematureEstimated <- prematureVVV_Der_No4sU(x = tpts, parameters = parameters, c = c)
	matureEstimated <- .impulseModel(x = tpts, par = matureParameters)

	if( any(is.na(D2)) | any(k1<0) | any(prematureEstimated<0) | any(matureEstimated<0) ) return(NaN)

	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)

	return(sum(c(prematureChiSquare,matureChiSquare)))
}

.inspect.engine_Derivative_No4sU <- function(tptsOriginal
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
		print("Some genes have a mature profile impossible to be fitted with impulsive 
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
	# 		print("More than 50% did not return a good fit of their mature profiles 
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
			print("Some genes have an expression profile impossible to be fitted 
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
		# 		print("More than 50% did not return a good fit of their premature 
		# 			profiles with the impulsive smooth function:
		# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
		# }

	}

	# Equal to integrative approach

	KKK <- bplapply(eiGenes,function(row){
	
	matureParameters <- mean(.impulseModel(tptsLinear, matureFitImpulse[[row]][1:6]))
	k2Parameters <- k2median
	k3Parameters <- k3median

	unlist(
		tryCatch(
			optim(c(matureParameters, k2Parameters, k3Parameters)
				,errorKKK_Int_No4sU
				,tpts = tptsLinear
				,premature = premature[row,]
				,mature = mature[row,]
				,prematureVariance = prematureVariance[row,]
				,matureVariance = matureVariance[row,]
				,control = list(maxit = 100*(nIter))),
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
	print("Model 0 finished.")

	medianAmplVKK <- sapply(eiGenes, function(row)
	{

		parameters <- unname(c(matureFitImpulse[[row]][1:6]
							 , k2median
							 , k3median))

		k1 <- k1VKK_Der_No4sU(tptsLinear,parameters,c)
		p <- prematureVKK_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)

		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)

		suppressWarnings(optimize( function(x)
		{
    		parameters <- unname(c(matureFitImpulse[[row]][1:6]
								 , k2median*x
								 , k3median*x))

			k1 <- k1VKK_Der_No4sU(tptsLinear,parameters,c)
    		p <- prematureVKK_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)

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
	                  	 ,errorVKK_Der_No4sU
        				 ,tpts = tptsLinear
			             ,premature = premature[row,]
			             ,mature = mature[row,]
			             ,prematureVariance = prematureVariance[row,]
			             ,matureVariance = matureVariance[row,]
	                  	 ,c = c
<<<<<<< HEAD
	                  	 ,control = list(maxit = nIter * 100)),
=======
	                  	 ,control = list(maxit = 100*(nIter))),
>>>>>>> f7d7f2c07ed09f1c3246f3f5954906b07f324835
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
	print("Model A finished.")

	medianAmplVVK <- sapply(eiGenes, function(row)
	{
		parameters <- c(unname(matureFitImpulse[[row]][1:6])
						, rep(k2median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1
						, k3median)
		k1 <- k1VVK_Der_No4sU(tptsLinear,parameters, c)
		p <- prematureVVK_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)
		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)
		suppressWarnings(optimize( function(x) {
			parameters <- unname(c(matureFitImpulse[[row]][1:6]
							  , c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
							  , k3median * x))
     
			k1 <- k1VVK_Der_No4sU(tptsLinear,parameters, c)
			p <- prematureVVK_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)      
     
			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(c(k1,.impulseModel(tptsLinear,c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)),k3median*x*length(tptsOriginal)))
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
			        	,errorVVK_Der_No4sU
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,c = c
<<<<<<< HEAD
						,control = list(maxit = nIter * 100)),
=======
						,control = list(maxit = 100*(nIter))),
>>>>>>> f7d7f2c07ed09f1c3246f3f5954906b07f324835
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
	print("Model AC finished.")

	medianAmplVKV <- sapply(eiGenes, function(row)
	{

		parameters <- c(unname(matureFitImpulse[[row]][1:6])
							 , k2median
			  				 , rep(k3median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
  		k1 <- k1VKV_Der_No4sU(tptsLinear,parameters, c)
		p <- prematureVKV_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)

		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)

		suppressWarnings(optimize( function(x) {

			parameters <- c(unname(matureFitImpulse[[row]][1:6])
								 , k2median * x
								 , rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
      
			k1 <- k1VKV_Der_No4sU(tptsLinear,parameters, c)
			p <- prematureVKV_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)
      
			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(k1,k2median*x*length(tptsOriginal),.impulseModel(tptsLinear,c(rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)))

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
          				,errorVKV_Der_No4sU
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,c = c
<<<<<<< HEAD
						,control = list(maxit = nIter * 100)),
=======
						,control = list(maxit = 100*(nIter))),
>>>>>>> f7d7f2c07ed09f1c3246f3f5954906b07f324835
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
	print("Model AB finished.")

	medianAmplVVV <- sapply(eiGenes, function(row)
	{

		parameters <- c(unname(matureFitImpulse[[row]][1:6])
					  , rep(k2median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1
					  , rep(k3median,3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
		k1 <- k1VVV_Der_No4sU(tptsLinear,parameters, c)
		p <- prematureVVV_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)
	
		if(all(k1>0,na.rm=TRUE) & all(p>0,na.rm=TRUE) & all(is.finite(k1)) & all(is.finite(p))) return(1)
	
		suppressWarnings(optimize( function(x) {
	
			parameters <- c(unname(matureFitImpulse[[row]][1:6])
						  , rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1
						  , rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
	
			k1 <- k1VVV_Der_No4sU(tptsLinear,parameters, c)
			p <- prematureVVV_Der_No4sU(x = tptsLinear, parameters = parameters, c = c)      
	      
			if(!all(k1>0,na.rm=TRUE) | !all(p>0,na.rm=TRUE) | !all(is.finite(k1)) | !all(is.finite(p))) NaN else sum(k1
	      										, .impulseModel(tptsLinear,c(rep(k2median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1))
	      										, .impulseModel(tptsLinear,c(rep(k3median,3) * x, max(tptsLinear)/3, max(tptsLinear)/3*2, 1)))
	
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
          				,errorVVV_Der_No4sU
						,tpts = tptsLinear
						,premature = premature[row,]
						,mature = mature[row,]
						,prematureVariance = prematureVariance[row,]
						,matureVariance = matureVariance[row,]
						,c = c
<<<<<<< HEAD
						,control = list(maxit = nIter * 100)),
=======
						,control = list(maxit = 100*(nIter))),
>>>>>>> f7d7f2c07ed09f1c3246f3f5954906b07f324835
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
	print("Model ABC finished.")

	KKV <- bplapply(eiGenes, function(row){
 
		k1Parameters <- mean(k1VKV_Der_No4sU(tptsLinear, VKV[[row]], c),na.rm=T)
		k2Parameters <- VKV[[row]][7]
		k3Parameters <- VKV[[row]][8:13]
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
		        	,errorKKV_Int_No4sU
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
	print("Model B finished.")

	KVK <- bplapply(eiGenes, function(row){
 
		k1Parameters <- mean(k1VVK_Der_No4sU(tptsLinear, VVK[[row]], c),na.rm=T)
		k2Parameters <- VVK[[row]][7:12]
		k3Parameters <- VVK[[row]][13]
		
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVK_Int_No4sU
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
	print("Model C finished.")

	KVV <- bplapply(eiGenes, function(row){
	
		k1Parameters <- mean(k1VVV_Der_No4sU(tptsLinear, VVV[[row]], c),na.rm=T)
		k2Parameters <- VVV[[row]][7:12]
		k3Parameters <- VVV[[row]][13:18]
		unlist(
			tryCatch(
				optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
					 ,errorKVV_Int_No4sU
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
	print("Model BC finished.")

	if(testOnSmooth)
	{
		matureSmooth <- t(sapply(matureFitImpulse,function(i)
		{
 			.impulseModel(tptsLinear,i)
		}))
		prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
		{
		 	.impulseModel(tptsLinear,i)
		}))
		colnames(matureSmooth) <- tptsOriginal
		colnames(prematureSmooth) <- tptsOriginal
	}else{
		matureSmooth <- mature
		prematureSmooth <- premature
	}

	chi2data <- t(mcsapply(eiGenes,function(g)
	{
		KKKTemp <- tryCatch(errorKKK_Int_No4sU(KKK[[g]][grep("par",names(KKK[[g]]))]
											  ,tptsLinear
											  ,prematureSmooth[g,]
											  ,matureSmooth[g,]
											  ,prematureVariance[g,]
											  ,matureVariance[g,]),error = function(e)NaN)
	
		VKKTemp <- tryCatch(errorVKK_Der_No4sU(VKK[[g]][grep("par",names(VKK[[g]]))]
											  ,tptsLinear
											  ,prematureSmooth[g,]
											  ,matureSmooth[g,]
											  ,prematureVariance[g,]
											  ,matureVariance[g,]
											  ,c),error = function(e)NaN)
	
		KVKTemp <- tryCatch(errorKVK_Int_No4sU(KVK[[g]][grep("par",names(KVK[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		KKVTemp <- tryCatch(errorKKV_Int_No4sU(KKV[[g]][grep("par",names(KKV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)
	
		VVKTemp <- tryCatch(errorVVK_Der_No4sU(VVK[[g]][grep("par",names(VVK[[g]]))]
													 ,tptsLinear
													 ,prematureSmooth[g,]
													 ,matureSmooth[g,]
													 ,prematureVariance[g,]
													 ,matureVariance[g,]
													 ,c),error = function(e)NaN)
	
		VKVTemp <- tryCatch(errorVKV_Der_No4sU(VKV[[g]][grep("par",names(VKV[[g]]))]
													 ,tptsLinear
													 ,prematureSmooth[g,]
													 ,matureSmooth[g,]
													 ,prematureVariance[g,]
													 ,matureVariance[g,]
													 ,c),error = function(e)NaN)

		KVVTemp <- tryCatch(errorKVV_Int_No4sU(KVV[[g]][grep("par",names(KVV[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)
	
		VVVTemp <- tryCatch(errorVVV_Der_No4sU(VVV[[g]][grep("par",names(VVV[[g]]))]
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
	pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], 2*length(tptsOriginal)-dof['KKK'])
						,VKK=pchisq(chi2data[,'VKK'], 2*length(tptsOriginal)-dof['VKK'])
						,KVK=pchisq(chi2data[,'KVK'], 2*length(tptsOriginal)-dof['KVK'])
						,KKV=pchisq(chi2data[,'KKV'], 2*length(tptsOriginal)-dof['KKV'])
						,VVK=pchisq(chi2data[,'VVK'], 2*length(tptsOriginal)-dof['VVK'])
						,VKV=pchisq(chi2data[,'VKV'], 2*length(tptsOriginal)-dof['VKV'])
						,KVV=pchisq(chi2data[,'KVV'], 2*length(tptsOriginal)-dof['KVV'])
						,VVV=pchisq(chi2data[,'VVV'], 2*length(tptsOriginal)-dof['VVV']))

	logLikelihood <- t(mcsapply(eiGenes,function(g)
	{
		prematureKKKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureKKK_Int_No4sU(x = tptsLinear[t]
																						   , parameters = KKK[[g]][grep("par",names(KKK[[g]]))])))
		prematureVKKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVKK_Der_No4sU(x = tptsLinear[t]
		                                                                         , parameters = VKK[[g]][grep("par",names(VKK[[g]]))]
		                                                                         , c = c)))
		prematureVVKTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVVK_Der_No4sU(x = tptsLinear[t]
		                                                            , parameters = VVK[[g]][grep("par",names(VVK[[g]]))]
		                                                            , c = c)))
		prematureVKVTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVKV_Der_No4sU(x = tptsLinear[t]
		                                                            , parameters = VKV[[g]][grep("par",names(VKV[[g]]))]
		                                                            , c = c)))
		prematureVVVTemp <- c(sapply(seq_along(tptsLinear),function(t)prematureVVV_Der_No4sU(x = tptsLinear[t]
		                                                            , parameters = VVV[[g]][grep("par",names(VVV[[g]]))]
		                                                            , c = c)))		

		matureKKKTemp <- c(sapply(seq_along(tptsLinear),function(t)KKK[[g]][grep("par",names(KKK[[g]]))][1]))

		matureVKKTemp <- c(sapply(seq_along(tptsLinear),function(t).impulseModel(x = tptsLinear[t]
        		                                                               , par = VKK[[g]][grep("par",names(VKK[[g]]))][1:6])))
		matureVVKTemp <- c(sapply(seq_along(tptsLinear),function(t).impulseModel(x = tptsLinear[t]
                                                                      		   , par = VVK[[g]][grep("par",names(VVK[[g]]))][1:6])))
		matureVKVTemp <- c(sapply(seq_along(tptsLinear),function(t).impulseModel(x = tptsLinear[t]
                                                                      		   , par = VKV[[g]][grep("par",names(VKV[[g]]))][1:6])))
		matureVVVTemp <- c(sapply(seq_along(tptsLinear),function(t).impulseModel(x = tptsLinear[t]
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

		KKVTemp <- tryCatch(loglikKKV_Int_No4sU(KKV[[g]][grep("par",names(KKV[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)

		KVKTemp <- tryCatch(loglikKVK_Int_No4sU(KVK[[g]][grep("par",names(KVK[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)

		KVVTemp <- tryCatch(loglikKVV_Int_No4sU(KVV[[g]][grep("par",names(KVV[[g]]))]
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
 				"0" = list(alpha = list(fun = .constantModelP
									 ,type = "constant"
									 ,df = 1
									 ,params = c(alpha = unname(KKK[[gene]]["par1"])))
						,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(KKK[[gene]]["par3"])))
	 				   	,gamma = list(fun = .constantModelP
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

 				,a = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VKK[[gene]][1:6])))
 					  		,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(VKK[[gene]][8])))
	 				  		,gamma = list(fun = .constantModelP
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
 				,b = list(alpha = list(fun = .constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KKV[[gene]][1])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(KKV[[gene]][3:8])))
	 				  		,gamma = list(fun = .constantModelP
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
 				,c = list(alpha = list(fun = .constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KVK[[gene]][1])))
 					  		,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(KVK[[gene]][8])))
	 				  		,gamma = list(fun = .impulseModelP
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
 				,ab = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VKV[[gene]][1:6])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(VKV[[gene]][8:13])))
	 				  		,gamma = list(fun = .constantModelP
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
 				,ac = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VVK[[gene]][1:6])))
 					  		,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(VVK[[gene]][13])))
	 				  		,gamma = list(fun = .impulseModelP
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
 				,bc = list(alpha = list(fun = .constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KVV[[gene]][1])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(KVV[[gene]][8:13])))
	 				  		,gamma = list(fun = .impulseModelP
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
 				,abc = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VVV[[gene]][1:6])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(VVV[[gene]][13:18])))
	 				  		,gamma = list(fun = .impulseModelP
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

.makeModel_Derivative <- function(tpts, hyp, log_shift, .time_transf, ode, .rxnrate, c= NaN, geneBestModel = NULL)
{

	params <- list()
	params$alpha <- function(x) 
		hyp$alpha$fun$value(.time_transf(x, log_shift, c), hyp$alpha$par)
	params$beta  <- function(x) 
		hyp$beta$fun$value(.time_transf(x, log_shift, c), hyp$beta$par)
	params$gamma <- function(x) 
		hyp$gamma$fun$value(.time_transf(x, log_shift, c), hyp$gamma$par)

	matureTemp <- params$alpha(tpts)

	if(geneBestModel == "0")
	{
		prematureTemp <- sapply(tpts,function(t)prematureKKK_Int_No4sU(t,c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params)))
		k1Temp <- sapply(tpts,function(t)k1KKK_Int_No4sU(t,c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params)))
		
		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))
	
	}else if(geneBestModel == "a")
	{
		prematureTemp <- prematureVKK_Der_No4sU(.time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VKK_Der_No4sU(.time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))

	}else if(geneBestModel == "ac")
	{
		prematureTemp <- prematureVVK_Der_No4sU(.time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VVK_Der_No4sU(.time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- .impulseModel(.time_transf(tpts,log_shift,c), hyp$gamma$params)
		k3Temp <- rep(hyp$beta$params, length.out = length(tpts))

	}else if(geneBestModel == "ab")
	{
		prematureTemp <- prematureVKV_Der_No4sU(.time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VKV_Der_No4sU(.time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- rep(hyp$gamma$params, length.out = length(tpts))
		k3Temp <- .impulseModel(.time_transf(tpts,log_shift,c), hyp$beta$params)

	}else if(geneBestModel == "abc")
	{
		prematureTemp <- prematureVVV_Der_No4sU(.time_transf(tpts, log_shift, c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)
		k1Temp <- k1VVV_Der_No4sU(.time_transf(tpts,log_shift,c),c(hyp$alpha$params,hyp$gamma$params,hyp$beta$params),c)

		k2Temp <- .impulseModel(.time_transf(tpts,log_shift,c), hyp$gamma$params)
		k3Temp <- .impulseModel(.time_transf(tpts,log_shift,c), hyp$beta$params)

	}

	totalTemp <- matureTemp + prematureTemp

	data.frame(time = tpts, preMRNA = prematureTemp, total = totalTemp, alpha = k1Temp, beta = k3Temp, gamma = k2Temp)

}


.inspect.engine_Integrative_No4sU <- function(tptsOriginal
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
		# 		print("More than 50% did not return a good fit of their mature and
		# 			premature profiles with the impulsive smooth function:
		# 			consider using the option 'modelingParams()$testOnSmooth <- FALSE'.")
		# }

		if(all(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null)))
			stop("No genes have an expression profile possible to fit 
				with impulsive functions. Try with the option 'modelingParams()$testOnSmooth <- FALSE'.")

		if(any(sapply(matureFitImpulse,is.null) | sapply(prematureFitImpulse,is.null))) {
			print("Some genes have an expression profile impossible to be fitted 
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
						,errorKKK_Int_No4sU
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
		print("Model 0 finished.")
	
		VKK <- bplapply(eiGenes,function(row){

			k1Parameters <- c(rep(k1KKK_Int_No4sU(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- KKK[[row]][3]
	  
			unlist(
				tryCatch(
	      			optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
	                  	 ,errorVKK_Int_No4sU
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
		print("Model A finished.")

		KKV <- bplapply(eiGenes, function(row){
  
			k1Parameters <- KKK[[row]][1]
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			
			unlist(
				tryCatch(
					optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
			        	,errorKKV_Int_No4sU
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
		print("Model B finished.")

		KVK <- bplapply(eiGenes, function(row){
  
			k1Parameters <- KKK[[row]][1]
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k3Parameters <- KKK[[row]][3]
			
			unlist(
				tryCatch(
					optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorKVK_Int_No4sU
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
		print("Model C finished.")

		VKV <- bplapply(eiGenes, function(row){

			k1Parameters <- c(rep(k1KKK_Int_No4sU(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k2Parameters <- KKK[[row]][2]
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			
			unlist(
				tryCatch(
					optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorVKV_Int_No4sU
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
		print("Model AB finished.")

		VVK <- bplapply(eiGenes, function(row){

			k1Parameters <- c(rep(k1KKK_Int_No4sU(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k3Parameters <- KKK[[row]][3]

			unlist(
				tryCatch(
					optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
        				 ,errorVVK_Int_No4sU
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
		print("Model AC finished.")

		KVV <- bplapply(eiGenes, function(row){
  
			k1Parameters <- KKK[[row]][1]
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			unlist(
				tryCatch(
					optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						 ,errorKVV_Int_No4sU
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
		print("Model BC finished.")

		VVV <- bplapply(eiGenes, function(row){

			k1Parameters <- c(rep(k1KKK_Int_No4sU(0,KKK[[row]]),3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k2Parameters <- c(rep(KKK[[row]][2],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			k3Parameters <- c(rep(KKK[[row]][3],3), max(tptsLinear)/3, max(tptsLinear)/3*2, 1)
			
			unlist(
				tryCatch(
					optim(unname(c(k1Parameters, k2Parameters, k3Parameters))
						,errorVVV_Int_No4sU
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
		print("Model ABC finished.")

		if(testOnSmooth)
		{
				matureSmooth <- t(sapply(matureFitImpulse,function(i)
				{
  					.impulseModel(tptsLinear,i)
				}))
			
				prematureSmooth <- t(sapply(prematureFitImpulse,function(i)
				{
			  		.impulseModel(tptsLinear,i)
				}))

				colnames(matureSmooth) <- tptsOriginal
				colnames(prematureSmooth) <- tptsOriginal

		}else{
			matureSmooth <- mature
			prematureSmooth <- premature
		}

		chi2data <- t(mcsapply(eiGenes,function(g)
		{
			KKKTemp <- tryCatch(errorKKK_Int_No4sU(KKK[[g]][grep("par",names(KKK[[g]]))]
												  ,tptsOriginal
												  ,prematureSmooth[g,]
												  ,matureSmooth[g,]
												  ,prematureVariance[g,]
												  ,matureVariance[g,]),error = function(e)NaN)

			VKKTemp <- tryCatch(errorVKK_Int_No4sU(VKK[[g]][grep("par",names(VKK[[g]]))]
												  ,tptsOriginal
												  ,c(matureSmooth[g,],prematureSmooth[g,])
												  ,c(matureVariance[g,],prematureVariance[g,])
												  ,a
												  ,c),error = function(e)NaN)

			KVKTemp <- tryCatch(errorKVK_Int_No4sU(KVK[[g]][grep("par",names(KVK[[g]]))]
												  ,tptsOriginal
												  ,c(matureSmooth[g,],prematureSmooth[g,])
												  ,c(matureVariance[g,],prematureVariance[g,])
												  ,a
												  ,c),error = function(e)NaN)

			KKVTemp <- tryCatch(errorKKV_Int_No4sU(KKV[[g]][grep("par",names(KKV[[g]]))]
												  ,tptsOriginal
												  ,c(matureSmooth[g,],prematureSmooth[g,])
												  ,c(matureVariance[g,],prematureVariance[g,])
												  ,a
												  ,c),error = function(e)NaN)

			VVKTemp <- tryCatch(errorVVK_Int_No4sU(VVK[[g]][grep("par",names(VVK[[g]]))]
												  ,tptsOriginal
												  ,c(matureSmooth[g,],prematureSmooth[g,])
												  ,c(matureVariance[g,],prematureVariance[g,])
												  ,a
												  ,c),error = function(e)NaN)

			VKVTemp <- tryCatch(errorVKV_Int_No4sU(VKV[[g]][grep("par",names(VKV[[g]]))]
												  ,tptsOriginal
												  ,c(matureSmooth[g,],prematureSmooth[g,])
												  ,c(matureVariance[g,],prematureVariance[g,])
												  ,a
												  ,c),error = function(e)NaN)

			KVVTemp <- tryCatch(errorKVV_Int_No4sU(KVV[[g]][grep("par",names(KVV[[g]]))]
												  ,tptsOriginal
												  ,c(matureSmooth[g,],prematureSmooth[g,])
												  ,c(matureVariance[g,],prematureVariance[g,])
												  ,a
												  ,c),error = function(e)NaN)

			VVVTemp <- tryCatch(errorVVV_Int_No4sU(VVV[[g]][grep("par",names(VVV[[g]]))]
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
		pvaluesdata <- cbind(KKK=pchisq(chi2data[,'KKK'], 2*length(tptsOriginal)-dof['KKK'])
							,VKK=pchisq(chi2data[,'VKK'], 2*length(tptsOriginal)-dof['VKK'])
							,KVK=pchisq(chi2data[,'KVK'], 2*length(tptsOriginal)-dof['KVK'])
							,KKV=pchisq(chi2data[,'KKV'], 2*length(tptsOriginal)-dof['KKV'])
							,VVK=pchisq(chi2data[,'VVK'], 2*length(tptsOriginal)-dof['VVK'])
							,VKV=pchisq(chi2data[,'VKV'], 2*length(tptsOriginal)-dof['VKV'])
							,KVV=pchisq(chi2data[,'KVV'], 2*length(tptsOriginal)-dof['KVV'])
							,VVV=pchisq(chi2data[,'VVV'], 2*length(tptsOriginal)-dof['VVV']))

		logLikelihood <- t(mcsapply(eiGenes,function(g)
		{
			KKKTemp <- tryCatch(loglikKKK_Int_No4sU(KKK[[g]][grep("par",names(KKK[[g]]))]
											  ,tptsOriginal,prematureSmooth[g,]
											  ,matureSmooth[g,]
											  ,prematureVariance[g,]
											  ,matureVariance[g,]),error = function(e)NaN)

			VKKTemp <- tryCatch(loglikVKK_Int_No4sU(VKK[[g]][grep("par",names(VKK[[g]]))]
											  ,tptsOriginal
											  ,c(matureSmooth[g,],prematureSmooth[g,])
											  ,c(matureVariance[g,],prematureVariance[g,])
											  ,a
											  ,c),error = function(e)NaN)

		  	KVKTemp <- tryCatch(loglikKVK_Int_No4sU(KVK[[g]][grep("par",names(KVK[[g]]))]
            				                  ,tptsOriginal
            				                  ,c(matureSmooth[g,],prematureSmooth[g,])
            				                  ,c(matureVariance[g,],prematureVariance[g,])
            				                  ,a
            				                  ,c),error = function(e)NaN)

			KKVTemp <- tryCatch(loglikKKV_Int_No4sU(KKV[[g]][grep("par",names(KKV[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)
			VVKTemp <- tryCatch(loglikVVK_Int_No4sU(VVK[[g]][grep("par",names(VVK[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)
			VKVTemp <- tryCatch(loglikVKV_Int_No4sU(VKV[[g]][grep("par",names(VKV[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)
			KVVTemp <- tryCatch(loglikKVV_Int_No4sU(KVV[[g]][grep("par",names(KVV[[g]]))]
												   ,tptsOriginal
												   ,c(matureSmooth[g,],prematureSmooth[g,])
												   ,c(matureVariance[g,],prematureVariance[g,])
												   ,a
												   ,c),error = function(e)NaN)
			VVVTemp <- tryCatch(loglikVVV_Int_No4sU(VVV[[g]][grep("par",names(VVV[[g]]))]
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
 				"0" = list(alpha = list(fun = .constantModelP
									 ,type = "constant"
									 ,df = 1
									 ,params = c(alpha = k1KKK_Int_No4sU(x = 0,par = KKK[[gene]])))
 					    ,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(KKK[[gene]]["par3"])))
	 				   	,gamma = list(fun = .constantModelP
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

 				,a = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VKK[[gene]][1:6])))
 					  		,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(VKK[[gene]][8])))
	 				  		,gamma = list(fun = .constantModelP
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
 				,b = list(alpha = list(fun = .constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KKV[[gene]][1])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(KKV[[gene]][3:8])))
	 				  		,gamma = list(fun = .constantModelP
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
 				,c = list(alpha = list(fun = .constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KVK[[gene]][1])))
 					  		,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(KVK[[gene]][8])))
	 				  		,gamma = list(fun = .impulseModelP
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
 				,ab = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VKV[[gene]][1:6])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(VKV[[gene]][8:13])))
	 				  		,gamma = list(fun = .constantModelP
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
 				,ac = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VVK[[gene]][1:6])))
 					  		,beta = list(fun = .constantModelP
									,type = "constant"
									,df = 1
									,params = c(beta = unname(VVK[[gene]][13])))
	 				  		,gamma = list(fun = .impulseModelP
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
 				,bc = list(alpha = list(fun = .constantModelP
										,type = "constant"
										,df = 1
										,params = c(alpha = unname(KVV[[gene]][1])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(KVV[[gene]][8:13])))
	 				  		,gamma = list(fun = .impulseModelP
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
 				,abc = list(alpha = list(fun = .impulseModelP
										,type = "impulse"
										,df = 6
										,params = c(alpha = unname(VVV[[gene]][1:6])))
 					  		,beta = list(fun = .impulseModelP
									,type = "impulse"
									,df = 6
									,params = c(beta = unname(VVV[[gene]][13:18])))
	 				  		,gamma = list(fun = .impulseModelP
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


# 11/5/18

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
