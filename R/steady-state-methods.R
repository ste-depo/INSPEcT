#' @rdname compareSteady
#'
#' @description
#' This method compares two object of class INSPEcT in order to identify differential usage
#' of synthesis, processing or degradation rates in two different steady-state conditions. 
#' The two INSPEcT objects must have been profiled with replicates in order to provide
#' a statistical significance to the differences between their rates.
#' @param inspectIds An object of calss INSPEcT with two conditions
#' @param BPPARAM Configuration for BiocParallel parallelization. By default is set to SerialParam()
#' @return An object of class INSPEcT_diffsteady which contains both the absolute 
#' quantification of the rates as well as the comparison with the statistical significance
#' associated for each gene and rate. (See \code{\link{INSPEcT_diffsteady-class}})
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   data('allcounts', package='INSPEcT')
#'   data('featureWidths', package='INSPEcT')
#'   data('libsizes', package='INSPEcT')
#'   
#'   nascentCounts<-allcounts$nascent
#'   matureCounts<-allcounts$mature
#'   conditions<-letters[1:11]
#'   expDes<-rep(conditions,3)
#'   tL<-1/6
#'   
#'   nasExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'         allcounts=nascentCounts
#'         ,libsize=nascentLS
#'         ,exonsWidths=exWdths
#'         ,intronsWidths=intWdths
#'         ,experimentalDesign=expDes)
#'   
#'   matExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'         allcounts=matureCounts
#'         ,libsize=totalLS
#'         ,exonsWidths=exWdths
#'         ,intronsWidths=intWdths
#'         ,experimentalDesign=expDes)
#'  
#'   nasFullObj <- newINSPEcT(
#'         tpts=conditions
#'         ,labeling_time=tL
#'         ,nascentExpressions=nasExp_DESeq2
#'         ,matureExpressions=matExp_DESeq2)
#'   
#'   diffrates = compareSteady(nasFullObj[,c(1,11)])
#' }
setMethod('compareSteady', signature('INSPEcT'), function(inspectIds, BPPARAM=SerialParam()) 
{

	if( length(tpts(inspectIds))!=2 ) stop('compareSteady: two conditions at the time can be compared.')

	## get parameters for the analysis from the INSPEcT object

	ddp <- inspectIds@degDuringPulse
	cTsh <- 1 # accept all models
	sf <- inspectIds@labeledSF
	tL <- inspectIds@tL

	## divide the dataset in 2: full (eiGenes) and simple (eGenes)

	allGenes <- featureNames(inspectIds)
	eGenes <- allGenes[apply(is.na(ratesFirstGuess(inspectIds, 'processing')),1,all)]
	eiGenes <- setdiff(allGenes, eGenes)

	#########
	## FULL #####
	##############

	message(paste('Comparative analysis of the',length(eiGenes),'with both intronic and exonic signals.'))

	## concentrations

	wt_preT <- ratesFirstGuess(inspectIds,'preMRNA')[eiGenes,1]
	wt_matT <- ratesFirstGuess(inspectIds,'total')[eiGenes,1] - wt_preT
	wt_preTvar <- ratesFirstGuessVar(inspectIds,'preMRNA')[eiGenes,1]
	wt_matTvar <- ratesFirstGuessVar(inspectIds,'total')[eiGenes,1] + wt_preTvar

	sh_preT <- ratesFirstGuess(inspectIds,'preMRNA')[eiGenes,2]
	sh_matT <- ratesFirstGuess(inspectIds,'total')[eiGenes,2] - sh_preT
	sh_preTvar <- ratesFirstGuessVar(inspectIds,'preMRNA')[eiGenes,2]
	sh_matTvar <- ratesFirstGuessVar(inspectIds,'total')[eiGenes,2] + sh_preTvar

	## rates

	wt_rates <- cbind(
		k1=ratesFirstGuess(inspectIds,'synthesis')[eiGenes,1],
		k2=ratesFirstGuess(inspectIds,'processing')[eiGenes,1],
		k3=ratesFirstGuess(inspectIds,'degradation')[eiGenes,1])

	sh_rates <- cbind(
		k1=ratesFirstGuess(inspectIds,'synthesis')[eiGenes,2],
		k2=ratesFirstGuess(inspectIds,'processing')[eiGenes,2],
		k3=ratesFirstGuess(inspectIds,'degradation')[eiGenes,2])

	# give the rates corresponding to the other condition in case of NA

	# wt_rates[!is.finite(wt_rates)] <- sh_rates[!is.finite(wt_rates)]
	# sh_rates[!is.finite(sh_rates)] <- wt_rates[!is.finite(sh_rates)]

	# wt_rates[wt_rates==0] <- sh_rates[wt_rates==0]
	# sh_rates[sh_rates==0] <- wt_rates[sh_rates==0]

	wt_rates[wt_rates[,1]==0|is.na(wt_rates[,1]),1] <- min(wt_rates[wt_rates[,1]!=0 & !is.na(wt_rates[,1]),1])
	wt_rates[wt_rates[,2]==0|is.na(wt_rates[,2]),2] <- min(wt_rates[wt_rates[,2]!=0 & !is.na(wt_rates[,2]),2])
	wt_rates[wt_rates[,3]==0|is.na(wt_rates[,3]),3] <- min(wt_rates[wt_rates[,3]!=0 & !is.na(wt_rates[,3]),3])
	wt_rates[wt_rates[,1]==Inf,1] <- max(wt_rates[wt_rates[,1]!=Inf,1])
	wt_rates[wt_rates[,2]==Inf,2] <- max(wt_rates[wt_rates[,2]!=Inf,2])
	wt_rates[wt_rates[,3]==Inf,3] <- max(wt_rates[wt_rates[,3]!=Inf,3])

	sh_rates[sh_rates[,1]==0|is.na(sh_rates[,1]),1] <- min(sh_rates[sh_rates[,1]!=0 & !is.na(sh_rates[,1]),1])
	sh_rates[sh_rates[,2]==0|is.na(sh_rates[,2]),2] <- min(sh_rates[sh_rates[,2]!=0 & !is.na(sh_rates[,2]),2])
	sh_rates[sh_rates[,3]==0|is.na(sh_rates[,3]),3] <- min(sh_rates[sh_rates[,3]!=0 & !is.na(sh_rates[,3]),3])
	sh_rates[sh_rates[,1]==Inf,1] <- max(sh_rates[sh_rates[,1]!=Inf,1])
	sh_rates[sh_rates[,2]==Inf,2] <- max(sh_rates[sh_rates[,2]!=Inf,2])
	sh_rates[sh_rates[,3]==Inf,3] <- max(sh_rates[sh_rates[,3]!=Inf,3])

	mean_rates <- apply(array(cbind(wt_rates[,1:3], sh_rates[,1:3]), dim=c(nrow(wt_rates), 3,2)), 1:2, mean, na.rm=TRUE)

	if( !ddp ) {

		wt_k1 <- ratesFirstGuess(inspectIds,'synthesis')[eiGenes,1]
		wt_k1var <- ratesFirstGuessVar(inspectIds,'synthesis')[eiGenes,1]

		wt_data <- cbind(wt_matT, wt_preT, wt_k1)
		wt_datavar <- cbind(wt_matTvar, wt_preTvar, wt_k1var)

		sh_k1 <- ratesFirstGuess(inspectIds,'synthesis')[eiGenes,2]
		sh_k1var <- ratesFirstGuessVar(inspectIds,'synthesis')[eiGenes,2]

		sh_data <- cbind(sh_matT, sh_preT, sh_k1)
		sh_datavar <- cbind(sh_matTvar, sh_preTvar, sh_k1var)

		# give the variance corresponding to the other condition in case of Zeros

		wt_datavar[wt_datavar==0] <- sh_datavar[wt_datavar==0]
		sh_datavar[sh_datavar==0] <- wt_datavar[sh_datavar==0]

		wt_datavar[wt_datavar==0] <- min(wt_datavar[wt_datavar!=0])
		sh_datavar[sh_datavar==0] <- min(sh_datavar[sh_datavar!=0])

		# model!

		message('Evaluating model no-reg [1/8]...')
		outKKK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim(mean_rates[i,], errorKKK_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,],
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 7))
			,BPPARAM=BPPARAM))
		outKKK[apply(outKKK[,1:3]<0,1,any),] <- NaN

		message('Evaluating model s [2/8]...')
		outVKK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], mean_rates[i,2], mean_rates[i,3] ), 
					errorVKK_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,],  
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 8))
			,BPPARAM=BPPARAM))
		outVKK[apply(outVKK[,1:4]<0,1,any),] <- NaN

		message('Evaluating model p [3/8]...')
		outKVK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( mean_rates[i,1], wt_rates[i,2], sh_rates[i,2], mean_rates[i,3] ), 
					errorKVK_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 8))
			,BPPARAM=BPPARAM))
		outKVK[apply(outKVK[,1:4]<0,1,any),] <- NaN

		message('Evaluating model d [4/8]...')
		outKKV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( mean_rates[i,1], mean_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorKKV_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 8))
			,BPPARAM=BPPARAM))
		outKKV[apply(outKKV[,1:4]<0,1,any),] <- NaN

		message('Evaluating model sp [5/8]...')
		outVVK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2], mean_rates[i,3] ), 
					errorVVK_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 9))
			,BPPARAM=BPPARAM))
		outVVK[apply(outVVK[,1:5]<0,1,any),] <- NaN

		message('Evaluating model sd [6/8]...')
		outVKV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], mean_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorVKV_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 9))
			,BPPARAM=BPPARAM))
		outVKV[apply(outVKV[,1:5]<0,1,any),] <- NaN

		message('Evaluating model pd [7/8]...')
		outKVV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( mean_rates[i,1], wt_rates[i,2], sh_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorKVV_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 9))
			,BPPARAM=BPPARAM))
		outKVV[apply(outKVV[,1:5]<0,1,any),] <- NaN

		message('Evaluating model spd [8/8]...')
		outVVV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorVVV_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,],
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
				, error=function(e) rep(NaN, 10))
			,BPPARAM=BPPARAM))
		outVVV[apply(outVVV[,1:6]<0,1,any),] <- NaN

	} else { ## ddp

		datader <- c(0,0)

		wt_preL <- sf[1] * ratesFirstGuess(inspectIds,'labeled_preMRNA')[eiGenes,1]
		wt_matL <- sf[1] * ratesFirstGuess(inspectIds,'labeled_total')[eiGenes,1] - wt_preL
		wt_preLvar <- sf[1]^2 * ratesFirstGuessVar(inspectIds,'labeled_preMRNA')[eiGenes,1]
		wt_matLvar <- sf[1]^2 * ratesFirstGuessVar(inspectIds,'labeled_total')[eiGenes,1] + wt_preLvar

		wt_data <- cbind(wt_matT, wt_preT, wt_matL)
		wt_datavar <- cbind(wt_matTvar, wt_preTvar, wt_matLvar)

		sh_preL <- sf[2] * ratesFirstGuess(inspectIds,'labeled_preMRNA')[eiGenes,2]
		sh_matL <- sf[2] * ratesFirstGuess(inspectIds,'labeled_total')[eiGenes,2] - sh_preL
		sh_preLvar <- sf[2]^2 * ratesFirstGuessVar(inspectIds,'labeled_preMRNA')[eiGenes,2]
		sh_matLvar <- sf[2]^2 * ratesFirstGuessVar(inspectIds,'labeled_total')[eiGenes,2] + sh_preLvar

		sh_data <- cbind(sh_matT, sh_preT, sh_matL)
		sh_datavar <- cbind(sh_matTvar, sh_preTvar, sh_matLvar)

		# give the variance corresponding to the other condition in case of Zeros

		wt_datavar[wt_datavar==0] <- sh_datavar[wt_datavar==0]
		sh_datavar[sh_datavar==0] <- wt_datavar[sh_datavar==0]

		wt_datavar[wt_datavar==0] <- min(wt_datavar[wt_datavar!=0])
		sh_datavar[sh_datavar==0] <- min(sh_datavar[sh_datavar!=0])

		# model!

		message('Evaluating model no-reg [1/8]...')
		outKKK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim(mean_rates[i,], errorKKK_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader,
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 7))
			,BPPARAM=BPPARAM))
		outKKK[apply(outKKK[,1:3]<0,1,any),] <- NaN

		message('Evaluating model s [2/8]...')
		outVKK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], mean_rates[i,2], mean_rates[i,3] ), 
					errorVKK_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 8))
			,BPPARAM=BPPARAM))
		outVKK[apply(outVKK[,1:4]<0,1,any),] <- NaN

		message('Evaluating model p [3/8]...')
		outKVK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( mean_rates[i,1], wt_rates[i,2], sh_rates[i,2], mean_rates[i,3] ), 
					errorKVK_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 8))
			,BPPARAM=BPPARAM))
		outKVK[apply(outKVK[,1:4]<0,1,any),] <- NaN

		message('Evaluating model d [4/8]...')
		outKKV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( mean_rates[i,1], mean_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorKKV_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 8))
			,BPPARAM=BPPARAM))
		outKKV[apply(outKKV[,1:4]<0,1,any),] <- NaN

		message('Evaluating model sp [5/8]...')
		outVVK <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2], mean_rates[i,3] ), 
					errorVVK_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 9))
			,BPPARAM=BPPARAM))
		outVVK[apply(outVVK[,1:5]<0,1,any),] <- NaN

		message('Evaluating model sd [6/8]...')
		outVKV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], mean_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorVKV_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 9))
			,BPPARAM=BPPARAM))
		outVKV[apply(outVKV[,1:5]<0,1,any),] <- NaN

		message('Evaluating model pd [7/8]...')
		outKVV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( mean_rates[i,1], wt_rates[i,2], sh_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorKVV_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 9))
			,BPPARAM=BPPARAM))
		outKVV[apply(outKVV[,1:5]<0,1,any),] <- NaN

		message('Evaluating model spd [8/8]...')
		outVVV <- do.call('rbind', bplapply(seq_along(eiGenes), function(i)
			tryCatch(
				unlist(optim( c( wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2], wt_rates[i,3], sh_rates[i,3] ), 
					errorVVV_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
					sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
				, error=function(e) rep(NaN, 10))
			,BPPARAM=BPPARAM))
		outVVV[apply(outVVV[,1:6]<0,1,any),] <- NaN

	}

	models <- list(
		KKK=outKKK,
		VKK=outVKK,
		KVK=outKVK,
		KKV=outKKV,
		VVK=outVVK,
		VKV=outVKV,
		KVV=outKVV,
		VVV=outVVV
		)

	p_vals <- cbind(
		KKK=pchisq(models[['KKK']][,'value'], 3),
		VKK=pchisq(models[['VKK']][,'value'], 2),
		KVK=pchisq(models[['KVK']][,'value'], 2),
		KKV=pchisq(models[['KKV']][,'value'], 2),
		VVK=pchisq(models[['VVK']][,'value'], 1),
		VKV=pchisq(models[['VKV']][,'value'], 1),
		KVV=pchisq(models[['KVV']][,'value'], 1),
		VVV=pchisq(models[['VVV']][,'value'], 0)
		)
	p_vals[is.na(p_vals)] <- 1

	if( !ddp ) {
		log_liks <- cbind(
			KKK=sapply(seq_along(eiGenes), function(i) logLikKKK_ddpFALSE(models[['KKK']][i,1:3], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			VKK=sapply(seq_along(eiGenes), function(i) logLikVKK_ddpFALSE(models[['VKK']][i,1:4], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			KVK=sapply(seq_along(eiGenes), function(i) logLikKVK_ddpFALSE(models[['KVK']][i,1:4], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			KKV=sapply(seq_along(eiGenes), function(i) logLikKKV_ddpFALSE(models[['KKV']][i,1:4], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			VVK=sapply(seq_along(eiGenes), function(i) logLikVVK_ddpFALSE(models[['VVK']][i,1:5], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			VKV=sapply(seq_along(eiGenes), function(i) logLikVKV_ddpFALSE(models[['VKV']][i,1:5], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			KVV=sapply(seq_along(eiGenes), function(i) logLikKVV_ddpFALSE(models[['KVV']][i,1:5], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
			VVV=sapply(seq_along(eiGenes), function(i) logLikVVV_ddpFALSE(models[['VVV']][i,1:6], 
				wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,]))
			)		
	} else { ## ddp
		log_liks <- cbind(
			KKK=sapply(seq_along(eiGenes), function(i) logLikKKK_ddpTRUE(models[['KKK']][i,1:3], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			VKK=sapply(seq_along(eiGenes), function(i) logLikVKK_ddpTRUE(models[['VKK']][i,1:4], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			KVK=sapply(seq_along(eiGenes), function(i) logLikKVK_ddpTRUE(models[['KVK']][i,1:4], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			KKV=sapply(seq_along(eiGenes), function(i) logLikKKV_ddpTRUE(models[['KKV']][i,1:4], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			VVK=sapply(seq_along(eiGenes), function(i) logLikVVK_ddpTRUE(models[['VVK']][i,1:5], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			VKV=sapply(seq_along(eiGenes), function(i) logLikVKV_ddpTRUE(models[['VKV']][i,1:5], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			KVV=sapply(seq_along(eiGenes), function(i) logLikKVV_ddpTRUE(models[['KVV']][i,1:5], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
			VVV=sapply(seq_along(eiGenes), function(i) logLikVVV_ddpTRUE(models[['VVV']][i,1:6], 
				wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL))
			)		
	}

	## filter genes with no resolved models (probably because they were solved with a negative rate in all
	## the models
	allNA <- apply(is.na(log_liks),1,all)
	if( any(allNA) ) {
		NallNA <- length(which(allNA))
		message(paste('Removing', NallNA, 'genes with unresolved models'))
	}
	eiOut <- list(
		eiGenes=eiGenes[!allNA],
		models=lapply(models, function(x) x[!allNA,]),
		p_vals=p_vals[!allNA,],
		log_liks=log_liks[!allNA,]
		)

	if( length( eGenes ) > 0 ) 
	{
		message(paste('Comparative analysis of the',length(eGenes),'without intronic signals.'))

		wt_matT <- ratesFirstGuess(inspectIds,'total')[eGenes,1]
		wt_matTvar <- ratesFirstGuessVar(inspectIds,'total')[eGenes,1]

		sh_matT <- ratesFirstGuess(inspectIds,'total')[eGenes,2]
		sh_matTvar <- ratesFirstGuessVar(inspectIds,'total')[eGenes,2]

		wt_rates <- cbind(
			k1=ratesFirstGuess(inspectIds,'synthesis')[eGenes,1],
			k3=ratesFirstGuess(inspectIds,'degradation')[eGenes,1])

		sh_rates <- cbind(
			k1=ratesFirstGuess(inspectIds,'synthesis')[eGenes,2],
			k3=ratesFirstGuess(inspectIds,'degradation')[eGenes,2])

		# give the rates corresponding to the other condition in case of NA or zeros

		# wt_rates[!is.finite(wt_rates)] <- sh_rates[!is.finite(wt_rates)]
		# sh_rates[!is.finite(sh_rates)] <- wt_rates[!is.finite(sh_rates)]

		# wt_rates[wt_rates==0] <- sh_rates[wt_rates==0]
		# sh_rates[sh_rates==0] <- wt_rates[sh_rates==0]

		# wt_rates[wt_rates[,1]==0,1] <- median(wt_rates[,1])
		# wt_rates[wt_rates[,2]==0,2] <- median(wt_rates[,2])
		# sh_rates[sh_rates[,1]==0,1] <- median(sh_rates[,1])
		# sh_rates[sh_rates[,2]==0,2] <- median(sh_rates[,2])

		wt_rates[wt_rates[,1]==0|is.na(wt_rates[,1]),1] <- min(wt_rates[wt_rates[,1]!=0 & !is.na(wt_rates[,1]),1])
		wt_rates[wt_rates[,2]==0|is.na(wt_rates[,2]),2] <- min(wt_rates[wt_rates[,2]!=0 & !is.na(wt_rates[,2]),2])
		wt_rates[wt_rates[,1]==Inf,1] <- max(wt_rates[wt_rates[,1]!=Inf,1])
		wt_rates[wt_rates[,2]==Inf,2] <- max(wt_rates[wt_rates[,2]!=Inf,2])

		sh_rates[sh_rates[,1]==0|is.na(sh_rates[,1]),1] <- min(sh_rates[sh_rates[,1]!=0 & !is.na(sh_rates[,1]),1])
		sh_rates[sh_rates[,2]==0|is.na(sh_rates[,2]),2] <- min(sh_rates[sh_rates[,2]!=0 & !is.na(sh_rates[,2]),2])
		sh_rates[sh_rates[,1]==Inf,1] <- max(sh_rates[sh_rates[,1]!=Inf,1])
		sh_rates[sh_rates[,2]==Inf,2] <- max(sh_rates[sh_rates[,2]!=Inf,2])

		mean_rates <- apply(array(cbind(wt_rates[,1:2], sh_rates[,1:2]), dim=c(nrow(wt_rates),2,2)), 1:2, mean, na.rm=TRUE)

		if( !ddp ) {

			wt_k1 <- ratesFirstGuess(inspectIds,'synthesis')[eGenes,1]
			wt_k1var <- ratesFirstGuessVar(inspectIds,'synthesis')[eGenes,1]

			wt_data <- cbind(wt_matT, wt_k1)
			wt_datavar <- cbind(wt_matTvar, wt_k1var)

			sh_k1 <- ratesFirstGuess(inspectIds,'synthesis')[eGenes,2]
			sh_k1var <- ratesFirstGuessVar(inspectIds,'synthesis')[eGenes,2]

			sh_data <- cbind(sh_matT, sh_k1)
			sh_datavar <- cbind(sh_matTvar, sh_k1var)

			rownames(wt_datavar) <- rownames(sh_datavar) <- eGenes

			# give the variance corresponding to the other condition in case of Zeros

			wt_datavar[wt_datavar==0] <- sh_datavar[wt_datavar==0]
			sh_datavar[sh_datavar==0] <- wt_datavar[sh_datavar==0]

			wt_datavar[wt_datavar==0] <- min(wt_datavar[wt_datavar!=0])
			sh_datavar[sh_datavar==0] <- min(sh_datavar[sh_datavar!=0])

			# model!

			message('Evaluating model no-reg [1/4]...')
			outK_K <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim(mean_rates[i,], errorK_K_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,],
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
					, error=function(e) c(mean_rates[i,], NA, NA, NA, e[1]))
				,BPPARAM=BPPARAM))
			outK_K[apply(outK_K[,1:2]<0,1,any),] <- NaN

			message('Evaluating model s [2/4]...')
			outV_K <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim( c( wt_rates[i,1], sh_rates[i,1], mean_rates[i,2] ), 
						errorV_K_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,],  
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
					, error=function(e) c(wt_rates[i,1], sh_rates[i,1], mean_rates[i,2], NA, NA, NA, e[1]))
				,BPPARAM=BPPARAM))
			outV_K[apply(outV_K[,1:3]<0,1,any),] <- NaN

			message('Evaluating model d [3/4]...')
			outK_V <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim( c( mean_rates[i,1], wt_rates[i,2], sh_rates[i,2] ), 
						errorK_V_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
					, error=function(e) c(mean_rates[i,1], wt_rates[i,2], sh_rates[i,2], NA, NA, NA, e[1]))
				,BPPARAM=BPPARAM))
			outK_V[apply(outK_V[,1:3]<0,1,any),] <- NaN
			
			message('Evaluating model sd [4/4]...')
			outV_V <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim( c( wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2] ), 
						errorV_V_ddpFALSE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], 
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,])[1:4])
					, error=function(e) c(wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2], NA, NA, NA, e[1]))
				,BPPARAM=BPPARAM))
			outV_V[apply(outV_V[,1:4]<0,1,any),] <- NaN

		} else { ## ddp

			wt_matL <- sf[1] * ratesFirstGuess(inspectIds,'labeled_total')[eGenes,1]
			wt_matLvar <- sf[1]^2 * ratesFirstGuessVar(inspectIds,'labeled_total')[eGenes,1]

			wt_data <- cbind(wt_matT, wt_matL)
			wt_datavar <- cbind(wt_matTvar, wt_matLvar)

			sh_matL <- sf[2] * ratesFirstGuess(inspectIds,'labeled_total')[eGenes,1]
			sh_matLvar <- sf[2]^2 * ratesFirstGuessVar(inspectIds,'labeled_total')[eGenes,1]

			sh_data <- cbind(sh_matT, sh_matL)
			sh_datavar <- cbind(sh_matTvar, sh_matLvar)

			datader <- c(0)

			# give the variance corresponding to the other condition in case of zeros

			wt_datavar[wt_datavar==0] <- sh_datavar[wt_datavar==0]
			sh_datavar[sh_datavar==0] <- wt_datavar[sh_datavar==0]

			wt_datavar[wt_datavar==0] <- min(wt_datavar[wt_datavar!=0])
			sh_datavar[sh_datavar==0] <- min(sh_datavar[sh_datavar!=0])

			# model!

			message('Evaluating model no-reg [1/4]...')
			outK_K <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim(mean_rates[i,], errorKK_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader,
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
					, error=function(e) rep(NaN, 6))
				,BPPARAM=BPPARAM))

			message('Evaluating model s [2/4]...')
			outV_K <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim( c( wt_rates[i,1], sh_rates[i,1], mean_rates[i,2] ), 
						errorVK_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
					, error=function(e) rep(NaN, 7))
				,BPPARAM=BPPARAM))

			message('Evaluating model d [3/4]...')
			outK_V <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim( c( mean_rates[i,1], wt_rates[i,2], sh_rates[i,2] ), 
						errorKV_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
					, error=function(e) rep(NaN, 7))
				,BPPARAM=BPPARAM))

			message('Evaluating model sd [4/4]...')
			outV_V <- do.call('rbind', bplapply(seq_along(eGenes), function(i)
				tryCatch(
					unlist(optim( c( wt_rates[i,1], sh_rates[i,1], wt_rates[i,2], sh_rates[i,2] ), 
						errorVV_ddpTRUE, wt_data=wt_data[i,], wt_datavar=wt_datavar[i,], wt_datader=datader, 
						sh_data=sh_data[i,], sh_datavar=sh_datavar[i,], sh_datader=datader, tL=tL)[1:4])
					, error=function(e) rep(NaN, 8))
				,BPPARAM=BPPARAM))

		}

		models <- list(
			'K_K'=outK_K,
			'V_K'=outV_K,
			'K_V'=outK_V,
			'V_V'=outV_V
			)

		p_vals <- cbind(
			K_K=pchisq(models[['K_K']][,'value'], 2),
			V_K=pchisq(models[['V_K']][,'value'], 1),
			K_V=pchisq(models[['K_V']][,'value'], 1),
			V_V=pchisq(models[['V_V']][,'value'], 0)
			)
		p_vals[is.na(p_vals)] <- 1

		if( !ddp ) {
			log_liks <- cbind(
				K_K=sapply(seq_along(eGenes), function(i) logLikK_K_ddpFALSE(models[['K_K']][i,1:3], 
					wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
				V_K=sapply(seq_along(eGenes), function(i) logLikV_K_ddpFALSE(models[['V_K']][i,1:5], 
					wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
				K_V=sapply(seq_along(eGenes), function(i) logLikK_V_ddpFALSE(models[['K_V']][i,1:5], 
					wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,])),
				V_V=sapply(seq_along(eGenes), function(i) logLikV_V_ddpFALSE(models[['V_V']][i,1:6], 
					wt_data[i,] , wt_datavar[i,], sh_data[i,], sh_datavar[i,]))
				)
		} else {
			log_liks <- cbind(
				K_K=sapply(seq_along(eGenes), function(i) logLikKK_ddpTRUE(models[['K_K']][i,1:3], 
					wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
				V_K=sapply(seq_along(eGenes), function(i) logLikVK_ddpTRUE(models[['V_K']][i,1:5], 
					wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
				K_V=sapply(seq_along(eGenes), function(i) logLikKV_ddpTRUE(models[['K_V']][i,1:5], 
					wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL)),
				V_V=sapply(seq_along(eGenes), function(i) logLikVV_ddpTRUE(models[['V_V']][i,1:6], 
					wt_data[i,] , wt_datavar[i,], datader, sh_data[i,], sh_datavar[i,], datader, tL))
				)			
		}

		## filter genes with no resolved models (probably because they were solved with a negative rate in all
		## the models
		allNA <- apply(is.na(log_liks),1,all)
		if( any(allNA) ) {
			NallNA <- length(which(allNA))
			message(paste('Removing', NallNA, 'genes with unresolved models'))
		}
		eOut <- list(
			eGenes=eGenes[!allNA],
			models=lapply(models, function(x) x[!allNA,]),
			p_vals=p_vals[!allNA,],
			log_liks=log_liks[!allNA,]
			)

	} else {
		message('No intronless genes found.')
		eOut <- NULL
	}

	model_steady_states_out= list(eiOut, eOut)

	##################
	### make the differential analysis based on the 
	### chi squared and brown test (and the respective thresholds cTsh and bTsh)
	##############################

	ie_model_steady_states_out <- model_steady_states_out[[1]]

	if( ddp ) tL <- ie_model_steady_states_out$tL
	eiGenes <- ie_model_steady_states_out$eiGenes
	models <- ie_model_steady_states_out$models
	log_liks <- ie_model_steady_states_out$log_liks

	k_pars <- c(3,4,4,4,5,5,5,6)
	AIC <- t(2*k_pars - 2*t(log_liks))
	diz <- c('KKK', 'VKK', 'KVK', 'KKV', 'VVK', 'VKV', 'KVV', 'VVV') 
	gene_class <- diz[apply(AIC, 1, which.min)]
	rates_pvals <- t(sapply(seq_along(gene_class), function(i) {
		switch(gene_class[i],
			'KKK' = c(
				k1=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'VKK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KKV'], 1, lower.tail=FALSE)
				),
			'VKK' = c(
				k1=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'VKK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE)
				),
			'KVK' = c(
				k1=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE)
				),
			'KKV' = c(
				k1=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KKK'] + 2*log_liks[i,'KKV'], 1, lower.tail=FALSE)
				),
			'VVK' = c(
				k1=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VVK'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VVK'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE)
				),
			'VKV' = c(
				k1=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VKK'] + 2*log_liks[i,'VKV'], 1, lower.tail=FALSE)
				),
			'KVV' = c(
				k1=pchisq(- 2*log_liks[i,'KVV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'KKV'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'KVK'] + 2*log_liks[i,'KVV'], 1, lower.tail=FALSE)
				),
			'VVV' = c(
				k1=pchisq(- 2*log_liks[i,'KVV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k2=pchisq(- 2*log_liks[i,'VKV'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE),
				k3=pchisq(- 2*log_liks[i,'VVK'] + 2*log_liks[i,'VVV'], 1, lower.tail=FALSE)
				)
			)
		}))
	
	# ei_p_vals <- p_vals <- ie_model_steady_states_out$p_vals

	# synthesis_tests <- cbind(
	# 	logLikRatioTest( log_liks[,1], log_liks[,2], 1 ),
	# 	logLikRatioTest( log_liks[,3], log_liks[,5], 1 ),
	# 	logLikRatioTest( log_liks[,4], log_liks[,6], 1 ),
	# 	logLikRatioTest( log_liks[,7], log_liks[,8], 1 )
	# 	)

	# synthesis_mask <- cbind(
	# 	p_vals[,1] <= cTsh | p_vals[,2] <= cTsh,
	# 	p_vals[,3] <= cTsh | p_vals[,5] <= cTsh,
	# 	p_vals[,4] <= cTsh | p_vals[,6] <= cTsh,
	# 	p_vals[,7] <= cTsh | p_vals[,8] <= cTsh
	# 	)

	# processing_tests <- cbind(
	# 	logLikRatioTest( log_liks[,1], log_liks[,3], 1 ),
	# 	logLikRatioTest( log_liks[,2], log_liks[,5], 1 ),
	# 	logLikRatioTest( log_liks[,4], log_liks[,7], 1 ),
	# 	logLikRatioTest( log_liks[,6], log_liks[,8], 1 )
	# 	)

	# processing_mask <- cbind(
	# 	p_vals[,1] <= cTsh | p_vals[,3] <= cTsh,
	# 	p_vals[,2] <= cTsh | p_vals[,5] <= cTsh,
	# 	p_vals[,4] <= cTsh | p_vals[,7] <= cTsh,
	# 	p_vals[,6] <= cTsh | p_vals[,8] <= cTsh
	# 	)

	# degradation_tests <- cbind(
	# 	logLikRatioTest( log_liks[,1], log_liks[,4], 1 ),
	# 	logLikRatioTest( log_liks[,2], log_liks[,6], 1 ),
	# 	logLikRatioTest( log_liks[,3], log_liks[,7], 1 ),
	# 	logLikRatioTest( log_liks[,5], log_liks[,8], 1 )
	# 	)

	# degradation_mask <- cbind(
	# 	p_vals[,1] <= cTsh | p_vals[,4] <= cTsh,
	# 	p_vals[,2] <= cTsh | p_vals[,6] <= cTsh,
	# 	p_vals[,3] <= cTsh | p_vals[,7] <= cTsh,
	# 	p_vals[,5] <= cTsh | p_vals[,8] <= cTsh
	# 	)

	# rates_pvals <- cbind(
	# 	k1=brown_method_mask(synthesis_tests, synthesis_mask),
	# 	k2=brown_method_mask(processing_tests, processing_mask),
	# 	k3=brown_method_mask(degradation_tests, degradation_mask)
	# 	)
	# rates_pvals[is.na(rates_pvals)] <- 1

	# gene_class <- apply( rates_pvals , 1 , function(x) {vec <- x <= bTsh; paste(ifelse(vec, 'V', 'K'), collapse='')} )

	# diz <- c('KKK'=1, 'VKK'=2, 'KVK'=3, 'KKV'=4, 'VVK'=5, 'VKV'=6, 'KVV'=7, 'VVV'=8) 
	# numeric_gene_class <- diz[gene_class]
	# p_vals_geneclass <- p_vals[cbind(seq_along(eiGenes),numeric_gene_class)]

	# rates_c <- t(sapply( seq_along(gene_class) , function(i) {

	# 	switch(
	# 		gene_class[i],
	# 		'KKK'=models[['KKK']][i, c(1,2,3)],
	# 		'VKK'=models[['VKK']][i, c(1,3,4)],
	# 		'KVK'=models[['KVK']][i, c(1,2,4)],
	# 		'KKV'=models[['KKV']][i, c(1,2,3)],
	# 		'VVK'=models[['VVK']][i, c(1,3,5)],
	# 		'VKV'=models[['VKV']][i, c(1,3,4)],
	# 		'KVV'=models[['KVV']][i, c(1,2,4)],
	# 		'VVV'=models[['VVV']][i, c(1,3,5)]
	# 		)

	# 	}))

	# if( ddp ) {
	# 	concentrations_c <- t(apply(rates_c, 1, sys4suModel_ddpTRUE, datader=c(0,0), tL=tL))[,1:2]	
	# } else {
	# 	concentrations_c <- t(apply(rates_c, 1, sys4suModel_ddpFALSE))[,1:2]
	# }
	
	# rates_t <- t(sapply( seq_along(gene_class) , function(i) {

	# 	switch(
	# 		gene_class[i],
	# 		'KKK'=models[['KKK']][i, c(1,2,3)],
	# 		'VKK'=models[['VKK']][i, c(2,3,4)],
	# 		'KVK'=models[['KVK']][i, c(1,3,4)],
	# 		'KKV'=models[['KKV']][i, c(1,2,4)],
	# 		'VVK'=models[['VVK']][i, c(2,4,5)],
	# 		'VKV'=models[['VKV']][i, c(2,3,5)],
	# 		'KVV'=models[['KVV']][i, c(1,3,5)],
	# 		'VVV'=models[['VVV']][i, c(2,4,6)]
	# 		)

	# 	}))

	# if( ddp ) {
	# 	concentrations_t <- t(apply(rates_t, 1, sys4suModel_ddpTRUE, datader=c(0,0), tL=tL))[,1:2]
	# } else {
	# 	concentrations_t <- t(apply(rates_t, 1, sys4suModel_ddpFALSE))[,1:2]
	# }

	rates_c <- models[['VVV']][, c(1,3,5)]
	rates_t <- models[['VVV']][, c(2,4,6)]

	ei_mat <- data.frame(#concentrations_c, 
		rates_c, 
		#concentrations_t, 
		rates_t, 
		rates_pvals#, FAL#SE 
		#, p_vals_geneclass, gene_class
		)
	colnames(ei_mat) <- c(#'mat_c','pre_c',
		'k1_c','k2_c','k3_c',
		# 'mat_t','pre_t',
		'k1_t','k2_t','k3_t',
		'p_k1','p_k2','p_k3'##, 'intronless'
		# ,'p','class'
		)
	rownames(ei_mat) <- eiGenes

	e_model_steady_states_out <- model_steady_states_out[[2]]

	if( !is.null(e_model_steady_states_out) )
	{

		models <- e_model_steady_states_out$models
		log_liks <- e_model_steady_states_out$log_liks

		k_pars <- c(2,3,3,4)
		AIC <- t(2*k_pars - 2*t(log_liks))
		diz <- c('K_K', 'V_K', 'K_V', 'V_V')
		gene_class <- diz[apply(AIC, 1, which.min)]
		rates_pvals <- t(sapply(seq_along(gene_class), function(i) {
			switch(gene_class[i],
				'K_K' = c(
					k1=pchisq(- 2*log_liks[i,'K_K'] + 2*log_liks[i,'V_K'], 1, lower.tail=FALSE),
					k3=pchisq(- 2*log_liks[i,'K_K'] + 2*log_liks[i,'K_V'], 1, lower.tail=FALSE)
					),
				'V_K' = c(
					k1=pchisq(- 2*log_liks[i,'K_K'] + 2*log_liks[i,'V_K'], 1, lower.tail=FALSE),
					k3=pchisq(- 2*log_liks[i,'V_K'] + 2*log_liks[i,'V_V'], 1, lower.tail=FALSE)
					),
				'K_V' = c(
					k1=pchisq(- 2*log_liks[i,'K_V'] + 2*log_liks[i,'V_V'], 1, lower.tail=FALSE),
					k3=pchisq(- 2*log_liks[i,'K_K'] + 2*log_liks[i,'K_V'], 1, lower.tail=FALSE)
					),
				'V_V' = c(
					k1=pchisq(- 2*log_liks[i,'K_V'] + 2*log_liks[i,'V_V'], 1, lower.tail=FALSE),
					k3=pchisq(- 2*log_liks[i,'V_K'] + 2*log_liks[i,'V_V'], 1, lower.tail=FALSE)
					)
				)
			}))

		# e_pvals <- p_vals <- e_model_steady_states_out$p_vals
		
		# synthesis_tests <- cbind(
		# 	logLikRatioTest( log_liks[,1], log_liks[,2], 1 ),
		# 	logLikRatioTest( log_liks[,3], log_liks[,4], 1 )
		# 	)

		# synthesis_mask <- cbind(
		# 	p_vals[,1] <= cTsh | p_vals[,2] <= cTsh,
		# 	p_vals[,3] <= cTsh | p_vals[,4] <= cTsh
		# 	)

		# degradation_tests <- cbind(
		# 	logLikRatioTest( log_liks[,1], log_liks[,3], 1 ),
		# 	logLikRatioTest( log_liks[,2], log_liks[,4], 1 )
		# 	)

		# degradation_mask <- cbind(
		# 	p_vals[,1] <= cTsh | p_vals[,3] <= cTsh,
		# 	p_vals[,2] <= cTsh | p_vals[,4] <= cTsh
		# 	)

		# rates_pvals <- cbind(
		# 	k1=brown_method_mask(synthesis_tests, synthesis_mask),
		# 	k3=brown_method_mask(degradation_tests, degradation_mask)
		# 	)
		# rates_pvals[is.na(rates_pvals)] <- 1

		# gene_class <- apply( rates_pvals, 1, function(x) {vec <- x <= bTsh[1:2]; paste(ifelse(vec, 'V', 'K'), collapse='_')} )

		# diz <- c('K_K'=1, 'V_K'=2, 'K_V'=3, 'V_V'=4) 
		# numeric_gene_class <- diz[gene_class]
		# p_vals_geneclass <- p_vals[cbind(seq_along(eGenes),numeric_gene_class)]

		# rates_c <- t(sapply( seq_along(gene_class) , function(i) {

		# 	switch(
		# 		gene_class[i],
		# 		'K_K'=models[['K_K']][i, c(1,2)],
		# 		'V_K'=models[['V_K']][i, c(1,3)],
		# 		'K_V'=models[['K_V']][i, c(1,2)],
		# 		'V_V'=models[['V_V']][i, c(1,3)]
		# 		)

		# 	}))

		# if( ddp ) {
		# 	concentrations_c <- t(apply(rates_c, 1, simplesys4suModel_ddpTRUE, datader=c(0), tL=tL))[,1]
		# } else {
		# 	concentrations_c <- t(apply(rates_c, 1, sys4suSimpleModel_ddpFALSE))[,1]
		# }

		# rates_t <- t(sapply( seq_along(gene_class) , function(i) {

		# 	switch(
		# 		gene_class[i],
		# 		'K_K'=models[['K_K']][i, c(1,2)],
		# 		'V_K'=models[['V_K']][i, c(2,3)],
		# 		'K_V'=models[['K_V']][i, c(1,3)],
		# 		'V_V'=models[['V_V']][i, c(2,4)]
		# 		)

		# 	}))

		# if( ddp ) {
		# 	concentrations_t <- t(apply(rates_t, 1, simplesys4suModel_ddpTRUE, datader=c(0), tL=tL))[,1]
		# } else {
		# 	concentrations_t <- t(apply(rates_t, 1, sys4suSimpleModel_ddpFALSE))[,1]
		# }

		rates_c <- models[['V_V']][, c(1,3)]
		rates_t <- models[['V_V']][, c(2,4)]

		e_mat <- data.frame(#concentrations_c, 
			rates_c, 
			#concentrations_t, 
			rates_t, 
			rates_pvals#, TRUE
			#, p_vals_geneclass, gene_class
			)
		colnames(e_mat) <- c(#'mat_c',
			'k1_c','k3_c',
			# 'mat_t',
			'k1_t','k3_t',
			'p_k1','p_k3'#, 'intronless'
			# ,'p','class'
			)
		rownames(e_mat) <- eGenes

		new_e_class <- as.character(e_mat$class)
		new_e_class <- paste(substr(new_e_class,1,1), '-', substr(new_e_class,2,2), sep='')
		e_mat <- data.frame(
			# mat_c = e_mat$mat_c,
			# pre_c = NA,
			k1_c = e_mat$k1_c,
			k2_c = NA,
			k3_c = e_mat$k3_c,
			# mat_t = e_mat$mat_t,
			# pre_t = NA,
			k1_t = e_mat$k1_t,
			k2_t = NA,
			k3_t = e_mat$k3_t,
			p_k1 = e_mat$p_k1,
			p_k2 = NA,
			p_k3 = e_mat$p_k3,
			# intronless = e_mat$intronless,
			# p = e_mat$p,
			# class = new_e_class,
			row.names = rownames(e_mat)
			)
		# e_pvals <- cbind(
		# 	KKK=e_pvals[,'K_K'],
		# 	VKK=e_pvals[,'V_K'],
		# 	KVK=NaN,
		# 	KKV=e_pvals[,'K_V'],
		# 	VVK=NaN,
		# 	VKV=e_pvals[,'V_V'],
		# 	KVV=NaN,
		# 	VVV=NaN
		# 	)
	} else {
		e_mat <- NULL
	}

	mat <- rbind(ei_mat, e_mat)

	synthesis_res <- data.frame(
		condition1=mat[,'k1_c'],
		# variance1=s1logvar,
		condition2=mat[,'k1_t'],
		# variance2=s2logvar,
		# samplesize1=s1n,
		# samplesize2=s2n,
		log2mean=log2(sqrt(mat[,'k1_c']*mat[,'k1_t'])),
		log2fc=log2(mat[,'k1_t']/mat[,'k1_c']),
		pval=mat[,'p_k1'],
		padj=p.adjust(mat[,'p_k1'], method='BH'),
		# intronless=mat[,'intronless'],
		row.names=c(eiGenes, eGenes)
		)

	processing_res <- data.frame(
		condition1=mat[,'k2_c'],
		# variance1=s1logvar,
		condition2=mat[,'k2_t'],
		# variance2=s2logvar,
		# samplesize1=s1n,
		# samplesize2=s2n,
		log2mean=log2(sqrt(mat[,'k2_c']*mat[,'k2_t'])),
		log2fc=log2(mat[,'k2_t']/mat[,'k2_c']),
		pval=mat[,'p_k2'],
		padj=p.adjust(mat[,'p_k2'], method='BH'),
		# intronless=mat[,'intronless'],
		row.names=c(eiGenes, eGenes)
		)

	degradation_res <- data.frame(
		condition1=mat[,'k3_c'],
		# variance1=s1logvar,
		condition2=mat[,'k3_t'],
		# variance2=s2logvar,
		# samplesize1=s1n,
		# samplesize2=s2n,
		log2mean=log2(sqrt(mat[,'k3_c']*mat[,'k3_t'])),
		log2fc=log2(mat[,'k3_t']/mat[,'k3_c']),
		pval=mat[,'p_k3'],
		padj=p.adjust(mat[,'p_k3'], method='BH'),
		# intronless=mat[,'intronless'],
		row.names=c(eiGenes, eGenes)
		)

	colnames(synthesis_res)[1:2] <- colnames(processing_res)[1:2] <- 
		colnames(degradation_res)[1:2] <- tpts(inspectIds)

	# modeling_res <- data.frame(
	# 	cbind(rbind(ei_p_vals, e_pvals)), # intronless=mat[,'intronless']),
	# 	row.names=c(eiGenes, eGenes)
	# 	)

	new_object <- new('INSPEcT_diffsteady')
	new_object@synthesis <- synthesis_res
	new_object@degradation <- degradation_res
	new_object@processing <- processing_res
	# new_object@modeling_res <- modeling_res

	return(new_object)

})

#' @rdname INSPEcT_diffsteady-class
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   data('allcounts', package='INSPEcT')
#'   data('featureWidths', package='INSPEcT')
#'   data('libsizes', package='INSPEcT')
#'   
#'   nascentCounts<-allcounts$nascent
#'   matureCounts<-allcounts$mature
#'   conditions<-letters[1:11]
#'   expDes<-rep(conditions,3)
#'   tL<-1/6
#'   
#'   nasExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'         allcounts=nascentCounts
#'         ,libsize=nascentLS
#'         ,exonsWidths=exWdths
#'         ,intronsWidths=intWdths
#'         ,experimentalDesign=expDes)
#'   
#'   matExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'         allcounts=matureCounts
#'         ,libsize=totalLS
#'         ,exonsWidths=exWdths
#'         ,intronsWidths=intWdths
#'         ,experimentalDesign=expDes)
#'  
#'   nasFullObj <- newINSPEcT(tpts=conditions,labeling_time=tL
#'         ,nascentExpressions=nasExp_DESeq2,matureExpressions=matExp_DESeq2)
#'   
#'   diffrates = compareSteady(nasFullObj[,c(1,11)])
#'   head(synthesis(diffrates))
#' }
setMethod('synthesis', 'INSPEcT_diffsteady', function(object) slot(object, 'synthesis'))
#' @rdname INSPEcT_diffsteady-class
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   head(processing(diffrates))
#' }
setMethod('processing', 'INSPEcT_diffsteady', function(object) slot(object, 'processing'))
#' @rdname INSPEcT_diffsteady-class
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   head(degradation(diffrates))
#' }
setMethod('degradation', 'INSPEcT_diffsteady', function(object) slot(object, 'degradation'))
#' @rdname INSPEcT_diffsteady-class
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   featureNames(diffrates)
#' }
setMethod('featureNames', 'INSPEcT_diffsteady', function(object) rownames(slot(object, 'synthesis')))
#' @rdname geneClass
#' @param ... specify the threshold for rate variability 'bTsh' in case of 'INSPEcT_diffsteady' objects (default = .1)
setMethod('geneClass', 'INSPEcT_diffsteady', function(object, ...) {
	arguments <- list(...)
	if( any(names(arguments) == 'bTsh') ) {
		bTsh <- arguments$bTsh
	} else {
		bTsh <- .1	
	}
	ratePvalsTmp = cbind(
		synthesis(object)$padj,
		processing(object)$padj,
		degradation(object)$padj
		)
	# gc = apply(padj_Vals, 1, function(x) {
	# 	classNum = 1+1*(x<bTsh)
	# 	classNum[is.na(classNum)] = 1
	# 	classLetter = c('K','V')[classNum]
	# 	paste(classLetter, collapse='')
	# 	})
	# return(factor(gc))
	geneClass <- apply(ratePvalsTmp,1,function(r)
	{
		r <- r < unlist(bTsh)
		if(all(is.na(r))) return(NA)
		if(!r[1]&!r[2]&!r[3]) return("no-reg") # 0
		if(r[1]&!r[2]&!r[3]) return("s") # a
		if(!r[1]&r[2]&!r[3]) return("p") # c
		if(!r[1]&!r[2]&r[3]) return("d") # b
		if(r[1]&!r[2]&r[3]) return("sd") # ab
		if(r[1]&r[2]&!r[3]) return("sp") # ac
		if(!r[1]&r[2]&r[3]) return("pd") # bc
		if(r[1]&r[2]&r[3]) return("spd") # abc
	})
	return(geneClass)
	})


#' @name plotMA
#' @title MA-plot from base means and log fold changes
NULL

#' @rdname plotMA
#' @description Visualize the comparison between the rates calculated from two different INSPEcT objects
#' profiled in steady-state conditions.
#' @param object An object of calss INSPEcT_diffsteady
#' @param ... Additional parameters, see Details section
#' @details
#' Possible arguments to "plotMA":
#' \itemize{
#' \item "rate" - A character, which represent the rate to be visualized, either "synthesis", "processing" or "degradation". By default, "synthesis" is chosen.
#' \item "padj" - A numeric, The p-adjusted threshold for significance. Genes with p-adjusted lower than the threshold will be depicted as orange triangles. By default set to -Inf, meaning that no genes will be highlighted.
#' \item "xlim" - A numeric vector of length 2, limits of x-axis, by default the range of the data.
#' \item "xlab" - A character, the label of x-axis, by default "log2 geometric mean"
#' \item "ylim" - A numeric vector of length 2, limits of y-axis, by default the range of the data.
#' \item "ylab" - A character, the label of y-axis, by default "log2 fold change"
#' \item "main" - A character, the title of the plot, by default the name of the visualized rate.
#' }
#' @seealso \url{http://en.wikipedia.org/wiki/MA_plot}
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   data('allcounts', package='INSPEcT')
#'   data('featureWidths', package='INSPEcT')
#'   data('libsizes', package='INSPEcT')
#'   
#'   nascentCounts<-allcounts$nascent
#'   matureCounts<-allcounts$mature
#'   conditions<-letters[1:11]
#'   expDes<-rep(conditions,3)
#'   tL<-1/6
#'   
#'   nasExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'         allcounts=nascentCounts
#'         ,libsize=nascentLS
#'         ,exonsWidths=exWdths
#'         ,intronsWidths=intWdths
#'         ,experimentalDesign=expDes)
#'   
#'   matExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'         allcounts=matureCounts
#'         ,libsize=totalLS
#'         ,exonsWidths=exWdths
#'         ,intronsWidths=intWdths
#'         ,experimentalDesign=expDes)
#'  
#'   nasFullObj <- newINSPEcT(tpts=conditions
#'         ,labeling_time=tL
#'         ,nascentExpressions=nasExp_DESeq2
#'         ,matureExpressions=matExp_DESeq2)
#'   
#'   diffrates = compareSteady(nasFullObj[,c(1,11)])
#'  
#'   plotMA(diffrates, padj=.01)
#' }
setMethod('plotMA', 'INSPEcT_diffsteady', function(object, ...) {
	addargs <- list(...)

	## argument "rate"
	feasiblerates <- c('synthesis','processing','degradation')
	if( any(names(addargs) == 'rate') ) {
		rate <- addargs[['rate']]
		if( ! rate %in% feasiblerates )
			stop('plotMA: Unrecognized "rate" argument.')
	} else rate <- "synthesis"
	
	## argument "padj"
	if( any(names(addargs) == 'padj') ) {
		padj <- addargs[['padj']]
		if( !is.numeric(padj) | padj<0 | padj>1 )
			stop('plotMA: "padj" must be numeric and between 0 and 1.')
	} else padj <- -Inf

	data <- slot(object, rate)
	x <- data$log2mean
	y <- data$log2fc
	signif_genes <- data$padj<padj
	ix <- is.na(x) | is.na(y)
	if( any(ix) ) {
		x <- x[!ix]
		y <- y[!ix]
		signif_genes <- signif_genes[!ix]
	}

	## argument "xlim"
	if( any(names(addargs) == 'xlim') ) {
		xlim <- addargs[['xlim']]
		if( !(is.numeric(xlim) & length(xlim)==2) )
			stop('plotMA: "xlim" must be a numeric of length 2.')
		x[x<xlim[1]] <- xlim[1]
		x[x>xlim[2]] <- xlim[2]
	} else xlim <- range(x, na.rm=TRUE)

	## argument "ylim"
	if( any(names(addargs) == 'ylim') ) {
		ylim <- addargs[['ylim']]
		if( !(is.numeric(ylim) & length(ylim)==2) )
			stop('plotMA: "ylim" must be a numeric of length 2.')
		y[y<ylim[1]] <- ylim[1]
		y[y>ylim[2]] <- ylim[2]
	} else ylim <- range(y, na.rm=TRUE)

	## argument "xlab"
	if( any(names(addargs) == 'xlab') ) {
		xlab <- addargs[['xlab']]
		if( !is.character(xlab) )
			stop('plotMA: "xlab" must be a character.')	 		
	} else xlab <- 'log2 geometric mean'

	## argument "ylab"
	if( any(names(addargs) == 'ylab') ) {
		ylab <- addargs[['ylab']]
		if( !is.character(ylab) )
			stop('plotMA: "ylab" must be a character.')	 		
	} else {
		condition1 <- colnames(data)[1]
		condition2 <- colnames(data)[2]
		ylab <- paste('log2 fold change',condition2,'vs',condition1)
	}

	## argument "main"
	if( any(names(addargs) == 'main') ) {
		main <- addargs[['main']]
		if( !is.character(main) )
			stop('plotMA: "main" must be a character.')	 		
	} else main <- rate

	smoothScatter(x, y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
	points(x[signif_genes], y[signif_genes], col='orange', pch=2)
	abline(h=0, lty=3)

	})


##############################
### LOW-LEVEL FUNCTIONS #########
############################

## complete model error function

sys4suModel_ddpFALSE <- function(par)
{
	k1 <- par[1]; k2 <- par[2]; k3 <- par[3]

	preTmod <- k1 / k2
	matTmod <- k1 / k3
	k1mod <- k1

	return(c(matTmod, preTmod, k1mod))

}

sys4suChisq_ddpFALSE <- function(par, data, datavar)
{

	model <- sys4suModel_ddpFALSE(par)
	sum( ( model - data )^2 / datavar )
		
}

sys4suLoglik_ddpFALSE <- function(par, data, datavar)
{

	model <- sys4suModel_ddpFALSE(par)
	sum(log(2*pnorm(-abs(data-model),mean=0,sd=sqrt(datavar))))

}

errorKKK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( par , wt_data, wt_datavar) + 
		sys4suChisq_ddpFALSE( par , sh_data, sh_datavar)
}

errorVKK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c(par[1], par[3], par[4]) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c(par[2], par[3], par[4]) , sh_data, sh_datavar)
}

errorKVK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c( par[1], par[3], par[4] ) , sh_data, sh_datavar )
}

errorKKV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c( par[1], par[2], par[3] ) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c( par[1], par[2], par[4] ) , sh_data, sh_datavar )
}

errorVVK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c( par[2], par[4], par[5] ) , sh_data, sh_datavar )
}

errorVKV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c( par[1], par[3], par[4] ) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c( par[2], par[3], par[5] ) , sh_data, sh_datavar )
}

errorKVV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c( par[1], par[3], par[5] ) , sh_data, sh_datavar )
}

errorVVV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suChisq_ddpFALSE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar ) + 
		sys4suChisq_ddpFALSE( c( par[2], par[4], par[6] ) , sh_data, sh_datavar )
}

logLikKKK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( par , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( par , sh_data, sh_datavar)
}

logLikVKK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c(par[1], par[3], par[4]) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c(par[2], par[3], par[4]) , sh_data, sh_datavar)
}

logLikKVK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c( par[1], par[3], par[4] ) , sh_data, sh_datavar )
}

logLikKKV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c( par[1], par[2], par[3] ) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c( par[1], par[2], par[4] ) , sh_data, sh_datavar )
}

logLikVVK_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c( par[2], par[4], par[5] ) , sh_data, sh_datavar )
}

logLikVKV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c( par[1], par[3], par[4] ) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c( par[2], par[3], par[5] ) , sh_data, sh_datavar )
}

logLikKVV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c( par[1], par[3], par[5] ) , sh_data, sh_datavar )
}

logLikVVV_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suLoglik_ddpFALSE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar ) + 
		sys4suLoglik_ddpFALSE( c( par[2], par[4], par[6] ) , sh_data, sh_datavar )
}

## simple model error function

sys4suSimpleModel_ddpFALSE <- function(par)
{
	k1 <- par[1]; k3 <- par[2]
	matTmod <- k1 / k3
	k1mod <- k1

	return(c(matTmod, k1mod))

}

sys4suSimpleChisq_ddpFALSE <- function(par, data, datavar)
{

	model <- sys4suSimpleModel_ddpFALSE(par)
	sum( ( model - data )^2 / datavar )
		
}

sys4suSimpleLoglik_ddpFALSE <- function(par, data, datavar)
{

	model <- sys4suSimpleModel_ddpFALSE(par)
	sum(log(2*pnorm(-abs(data-model),mean=0,sd=sqrt(datavar))))

}

errorK_K_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleChisq_ddpFALSE( par , wt_data, wt_datavar) + 
		sys4suSimpleChisq_ddpFALSE( par , sh_data, sh_datavar)
}

errorV_K_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleChisq_ddpFALSE( c(par[1], par[3]) , wt_data, wt_datavar ) + 
		sys4suSimpleChisq_ddpFALSE( c(par[2], par[3]) , sh_data, sh_datavar)
}

errorK_V_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleChisq_ddpFALSE( c( par[1], par[2] ) , wt_data, wt_datavar ) + 
		sys4suSimpleChisq_ddpFALSE( c( par[1], par[3] ) , sh_data, sh_datavar )
}

errorV_V_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleChisq_ddpFALSE( c( par[1], par[3] ) , wt_data, wt_datavar ) + 
		sys4suSimpleChisq_ddpFALSE( c( par[2], par[4] ) , sh_data, sh_datavar )
}

logLikK_K_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleLoglik_ddpFALSE( par , wt_data, wt_datavar) + 
		sys4suSimpleLoglik_ddpFALSE( par , sh_data, sh_datavar)
}

logLikV_K_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleLoglik_ddpFALSE( c(par[1], par[3]) , wt_data, wt_datavar ) + 
		sys4suSimpleLoglik_ddpFALSE( c(par[2], par[3]) , sh_data, sh_datavar)
}

logLikK_V_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleLoglik_ddpFALSE( c( par[1], par[2] ) , wt_data, wt_datavar ) + 
		sys4suSimpleLoglik_ddpFALSE( c( par[1], par[3] ) , sh_data, sh_datavar )
}

logLikV_V_ddpFALSE <- function( par , wt_data , wt_datavar, sh_data, sh_datavar )
{ 
	sys4suSimpleLoglik_ddpFALSE( c( par[1], par[3] ) , wt_data, wt_datavar ) + 
		sys4suSimpleLoglik_ddpFALSE( c( par[2], par[4] ) , sh_data, sh_datavar )
}

logLikRatioTest <- function(null, alt, deltadf)
{
	D <- - 2*null + 2*alt
	pchisq(D, deltadf, lower.tail=FALSE)
}


#################
### DDP ##############
###############

sysScaleChisq_ddpTRUE <- function(par, data, datavar, datader, tL)
{
	k1 <- par[1]; k2 <- par[2]; k3 <- par[3]; sf <- par[4]
	matTder <- datader[1]; preTder <- datader[1]

	preTmod <- ( k1 - preTder ) / k2
	matTmod <- ( k2 * preTmod - matTder ) / k3
	preLmod <- sf * k1 / k2 * ( 1 - exp( -k2 * tL ) )
	matLmod <- sf * k1 / k3 * ( k2 * ( 1 - exp( -k3 * tL ) ) - k3 * ( 1 - exp( -k2 * tL ) ) ) / ( k2 - k3 )

	model <- c(matTmod, preTmod, matLmod, preLmod)
	sum( ( model - data )^2 / datavar )

}

sys4suModel_ddpTRUE <- function(par, datader, tL)
{
	k1 <- par[1]; k2 <- par[2]; k3 <- par[3]
	matTder <- datader[1]; preTder <- datader[1]

	preTmod <- ( k1 - preTder ) / k2
	matTmod <- ( k2 * preTmod - matTder ) / k3
	# preLmod <- k1 / k2 * ( 1 - exp( -k2 * tL ) )
	matLmod <- k1 / k3 * ( k2 * ( 1 - exp( -k3 * tL ) ) - k3 * ( 1 - exp( -k2 * tL ) ) ) / ( k2 - k3 )

	return(c(matTmod, preTmod, matLmod))

}

sys4suChisq_ddpTRUE <- function(par, data, datavar, datader, tL)
{

	model <- sys4suModel_ddpTRUE(par, datader, tL)
	sum( ( model - data )^2 / datavar )
		
}

sys4suLoglik_ddpTRUE <- function(par, data, datavar, datader, tL)
{

	model <- sys4suModel_ddpTRUE(par, datader, tL)
	sum(log(2*pnorm(-abs(data-model),mean=0,sd=sqrt(datavar))))

}

priorRates_ddpTRUE <- function(data, datader, tL)
{
	data <- unname(data)
	matT <- data[1]; preT <- data[2]; matL <- data[3]; preL <- data[4]
	matTder <- datader[1]; preTder <- datader[1]
	k1 <- ( matL + preL ) / tL
	k2 <- ( k1 - preTder ) / preT
	k3 <- ( k2 * preT - matTder ) / matT
	c( k1 , k2 , k3 )
}

errorKKK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( par , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( par , sh_data, sh_datavar, sh_datader, tL)
}

errorVKK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c(par[1], par[3], par[4]) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c(par[2], par[3], par[4]) , sh_data, sh_datavar, sh_datader, tL)
}

errorKVK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c( par[1], par[3], par[4] ) , sh_data, sh_datavar, sh_datader, tL )
}

errorKKV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c( par[1], par[2], par[3] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c( par[1], par[2], par[4] ) , sh_data, sh_datavar, sh_datader, tL )
}

errorVVK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c( par[2], par[4], par[5] ) , sh_data, sh_datavar, sh_datader, tL )
}

errorVKV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c( par[1], par[3], par[4] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c( par[2], par[3], par[5] ) , sh_data, sh_datavar, sh_datader, tL )
}

errorKVV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c( par[1], par[3], par[5] ) , sh_data, sh_datavar, sh_datader, tL )
}

errorVVV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suChisq_ddpTRUE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suChisq_ddpTRUE( c( par[2], par[4], par[6] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikKKK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( par , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( par , sh_data, sh_datavar, sh_datader, tL)
}

logLikVKK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c(par[1], par[3], par[4]) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c(par[2], par[3], par[4]) , sh_data, sh_datavar, sh_datader, tL)
}

logLikKVK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c( par[1], par[3], par[4] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikKKV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c( par[1], par[2], par[3] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c( par[1], par[2], par[4] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikVVK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c( par[2], par[4], par[5] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikVKV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c( par[1], par[3], par[4] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c( par[2], par[3], par[5] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikKVV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c( par[1], par[2], par[4] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c( par[1], par[3], par[5] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikVVV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	sys4suLoglik_ddpTRUE( c( par[1], par[3], par[5] ) , wt_data, wt_datavar, wt_datader, tL) + 
		sys4suLoglik_ddpTRUE( c( par[2], par[4], par[6] ) , sh_data, sh_datavar, sh_datader, tL )
}

simplesys4suModel_ddpTRUE <- function(par, datader, tL)
{
	k1 <- par[1]; k3 <- par[2]
	matTder <- datader[1]

	matTmod <- ( k1 - matTder ) / k3
	matLmod <- k1 / k3 * ( 1 - exp( -k3 * tL ) )

	return(c(matTmod, matLmod))

}

simplesys4suChisq_ddpTRUE <- function(par, data, datavar, datader, tL)
{

	model <- simplesys4suModel_ddpTRUE(par, datader, tL)
	sum( ( model - data )^2 / datavar )
		
}

simplesys4suLoglik_ddpTRUE <- function(par, data, datavar, datader, tL)
{

	model <- simplesys4suModel_ddpTRUE(par, datader, tL)
	sum(log(2*pnorm(-abs(data-model),mean=0,sd=sqrt(datavar))))

}

simplepriorRates_ddpTRUE <- function(data, datader, tL)
{
	data <- unname(data)
	matT <- data[1]; matL <- data[2]
	matTder <- datader
	k1 <- matL / tL
	k3 <- ( k1 - matTder ) / matT
	c( k1 , k3 )
}

errorKK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suChisq_ddpTRUE( par , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suChisq_ddpTRUE( par , sh_data, sh_datavar, sh_datader, tL)
}

errorVK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suChisq_ddpTRUE( c( par[1], par[3] ) , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suChisq_ddpTRUE( c( par[2], par[3] ) , sh_data, sh_datavar, sh_datader, tL)
}

errorKV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suChisq_ddpTRUE( c( par[1], par[2] ) , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suChisq_ddpTRUE( c( par[1], par[3] ) , sh_data, sh_datavar, sh_datader, tL )
}

errorVV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suChisq_ddpTRUE( c( par[1], par[3] ) , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suChisq_ddpTRUE( c( par[2], par[4] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikKK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suLoglik_ddpTRUE( par , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suLoglik_ddpTRUE( par , sh_data, sh_datavar, sh_datader, tL)
}

logLikVK_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suLoglik_ddpTRUE( c( par[1], par[3] ) , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suLoglik_ddpTRUE( c( par[2], par[3] ) , sh_data, sh_datavar, sh_datader, tL)
}

logLikKV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suLoglik_ddpTRUE( c( par[1], par[2] ) , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suLoglik_ddpTRUE( c( par[1], par[3] ) , sh_data, sh_datavar, sh_datader, tL )
}

logLikVV_ddpTRUE <- function( par , wt_data , wt_datavar, wt_datader, sh_data, sh_datavar, sh_datader, tL )
{ 
	simplesys4suLoglik_ddpTRUE( c( par[1], par[3] ) , wt_data, wt_datavar, wt_datader, tL) + 
		simplesys4suLoglik_ddpTRUE( c( par[2], par[4] ) , sh_data, sh_datavar, sh_datader, tL )
}


