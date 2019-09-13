#' @rdname computeConfidenceIntervals
#'
#' @description
#' This function is used to compute the confidence intervals for a given set of modeled genes in the NoNascent scenario.
#' @param object An object of class INSPEcT_model
#' @return An object of class INSPEcT.

setMethod(f='computeConfidenceIntervals', 'INSPEcT', definition=function(object, singleGeneClass=NULL, BPPARAM=bpparam())
{
	if(!object@NoNascent){

		message("+ Nascent objects have confidence intervals by default!")
		return(object)
		
	} else {

		if(nGenes(object)==1&!is.null(singleGeneClass))
		{
			gc <- singleGeneClass
		}else{
			gc <- geneClass(object)
			names(gc) <- featureNames(object)
		}
		
		llConfidenceThreshold <- object@model@params$logLikelihoodConfidenceThreshold
		if(is.null(llConfidenceThreshold)) llConfidenceThreshold <- 0.95
		llConfidenceThreshold <- qchisq(llConfidenceThreshold,1)

		eiGenes <- featureNames(object)

		tpts <- tpts(object)

		total <- ratesFirstGuess(object,'total')
		premature <- ratesFirstGuess(object,'preMRNA')
		mature <- total - premature

		totalVariance <- ratesFirstGuessVar(object,'total')
		prematureVariance <- ratesFirstGuessVar(object,'preMRNA')
		matureVariance <- totalVariance + prematureVariance

		### Confidence intervals
 		print("Confidence intervals.")

		confidenceIntervals <- bplapply(eiGenes,function(g)
		{
		  classTmp <- gc[g]

			parameters <- unlist(object@model@ratesSpecs[[g]][[classTmp]])
			parameters <- unlist(parameters[grep("params",names(parameters))])
			
			parameters <- c(parameters[grep("total",names(parameters))],
							parameters[grep("mature",names(parameters))],
							parameters[grep("alpha",names(parameters))],
							parameters[grep("gamma",names(parameters))],
							parameters[grep("beta",names(parameters))])

			optTmp <- rates_derivativeModels(tpts=tpts, class=classTmp, parameters=parameters)

			foe <- capture.output({ # Just to capture the output of multiroot function
				suppressWarnings({
					intervals <- sapply(names(parameters),function(parname)
					{
						par <- parameters[parname]

							mOut = list(
								left_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e-2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = NULL, varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = NULL, confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
								left_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1/2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = NULL, varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = NULL, confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
								center = tryCatch(multiroot(f = logLikelihoodCIerror, start = par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = NULL, varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = NULL, confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
								right_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1.5*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = NULL, varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = NULL, confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList)),
								right_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = premature[g,], experimentalM = mature[g,], experimentalA = NULL, varianceP = prematureVariance[g,], varianceM = matureVariance[g,], varianceA = NULL, confidenceThreshold = llConfidenceThreshold, derivative = TRUE),error=function(e)return(emptyList))
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

		# Usefull for single gene (no nascent)
		if(!is.finite(k1_low)|k1_low==0) k1_low <- NaN
		if(!is.finite(k1_high)|k1_high==0) k1_high <- NaN
		if(!is.finite(k2_low)|k2_low==0) k2_low <- NaN
		if(!is.finite(k2_high)|k2_high==0) k2_high <- NaN
		if(!is.finite(k3_low)|k3_low==0) k3_low <- NaN
		if(!is.finite(k3_high)|k3_high==0) k3_high <- NaN

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
		# Check for single modeled genes (no nascent)
		if(!object@NoNascent)
		{
			eiGenes <- eiGenes[sapply(confidenceIntervals,function(g)all(is.finite(g[[1]]))&all(is.finite(g[[2]]))&all(is.finite(g[[3]])))]
			confidenceIntervals <- confidenceIntervals[eiGenes]
		}

		# I compute che constant rates
		fitResults_synthesis <- unlist(lapply(eiGenes,function(g)
		{
			rate_conf_int <- confidenceIntervals[[g]][["k1"]]
			k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
			if(!is.finite(k_start)) NaN #return(list(par=NaN, value=NaN))
			k_scores_out <- tryCatch(optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int),error=function(e)list('par'=NaN))
			return(k_scores_out$par)
		}))

		fitResults_processing <- unlist(lapply(eiGenes,function(g)
		{
			rate_conf_int <- confidenceIntervals[[g]][["k2"]]
			k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
			if(!is.finite(k_start)) NaN #return(list(par=NaN, value=NaN))
			k_scores_out <- tryCatch(optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int),error=function(e)list('par'=NaN))
			return(k_scores_out$par)
		}))

		fitResults_degradation <- unlist(lapply(eiGenes,function(g)
		{
			rate_conf_int <- confidenceIntervals[[g]][["k3"]]
			k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
			if(!is.finite(k_start)) NaN #return(list(par=NaN, value=NaN))
			k_scores_out <- tryCatch(optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int),error=function(e)list('par'=NaN))
			return(k_scores_out$par)
		}))

		names(fitResults_synthesis) <- 
		names(fitResults_processing) <- 
		names(fitResults_degradation) <- eiGenes

		confidenceIntervals <- lapply(eiGenes,function(g)
		{
			confidenceIntervals[[g]][['k1']] <- cbind(confidenceIntervals[[g]][['k1']],'constant'=rep(fitResults_synthesis[[g]],length(tpts)))
			confidenceIntervals[[g]][['k2']] <- cbind(confidenceIntervals[[g]][['k2']],'constant'=rep(fitResults_processing[[g]],length(tpts)))
			confidenceIntervals[[g]][['k3']] <- cbind(confidenceIntervals[[g]][['k3']],'constant'=rep(fitResults_degradation[[g]],length(tpts)))
			confidenceIntervals[[g]]
		})

		confidenceIntervals <- lapply(confidenceIntervals,function(i)lapply(i,function(j){j[!is.finite(j)]<-NaN;j}))

		object <- setConfidenceIntervals(object=object,confidenceIntervals=confidenceIntervals)

		return(object)
	}
})
