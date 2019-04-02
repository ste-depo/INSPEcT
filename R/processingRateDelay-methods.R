#' @rdname processingDelay
#'
#' @description
#' These functions calculates the tau and delta metrics for all genes with introns and exons in an oblect of class INSPEcT.
#' If the INSPEcT dataset was obtained with nascent RNA the metrics are caluclated using RNA dynamics and solving numerically the system of equations.
#' If the INSPEcT dataset was obtained without nascent RNA the metrics are approximated using premature and mature levels.
#' @param inspectIds An object of class INSPEcT.
#' @param tauThreshold A numeric representing the tau threshold to define a gene affected by processing. Default: 1.2
#' @param deltaThreshold A numeric representing the delta threshold to define a gene affected by processing. Default: 1.0
#' @param silent A logical indicating whether informaiton about the procedure should be printed or not.
#' @examples
#' data('allcounts', package='INSPEcT')
#' data('featureWidths', package='INSPEcT')
#' data('libsizes', package='INSPEcT')
#' 
#' nascentCounts<-allcounts$nascent
#' matureCounts<-allcounts$mature
#' conditions<-c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
#' expDes<-rep(conditions,3)
#' tL <- 1/6
#'
#' nasExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'       allcounts=matureCounts
#'       ,libsize=totalLS
#'       ,exonsWidths=exWdths
#'       ,intronsWidths=intWdths
#'       ,experimentalDesign=expDes)
#' 
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'       allcounts=matureCounts
#'       ,libsize=totalLS
#'       ,exonsWidths=exWdths
#'       ,intronsWidths=intWdths
#'       ,experimentalDesign=expDes)
#'
#' matureInspObj <- newINSPEcT(
#'       tpts=conditions
#'       ,labeling_time=tL
#'       ,nascentExpressions=nasExp_DESeq2
#'       ,matureExpressions=matExp_DESeq2)
#' 
#' procDelay<- processingDelay(inspectIds=matureInspObj
#'       ,tauThreshold=1.2
#'       ,deltaThreshold=1.0)
#'
#' head(procDelay)
#' table(procDelay)
#'
setMethod('processingDelay', signature('INSPEcT'), function(inspectIds, tauThreshold=1.2, deltaThreshold=1.0, silent=TRUE) {

	tau <- calculateTau(inspectIds, silent=silent)
	delta <- calculateDelta(inspectIds, silent=silent)
	tau > tauThreshold & delta > deltaThreshold

	})

#' @rdname processingDelay
#' @examples
#' head(calculateTau(matureInspObj))
#'
setMethod('calculateTau', signature('INSPEcT'), function(inspectIds, silent=FALSE) {

	eiGenes <- featureNames(inspectIds)[!apply(is.na(ratesFirstGuess(inspectIds, 'preMRNA')),1,all)]

	if( inspectIds@NoNascent ) {

		if( !silent ) 
			message(paste('No nascent dataset, approximated tau will be calculated on', 
				length(eiGenes),'genes with both intronic and exonic signals.'))

		tot <- ratesFirstGuess(inspectIds[eiGenes,], 'total')
		pre <- ratesFirstGuess(inspectIds[eiGenes,], 'preMRNA')
		mat <- tot - pre

		tau <- 1+pre/mat

	} else {

		if( !silent )
			message(paste('Tau will be calculated on', length(eiGenes),
				'genes with both intronic and exonic signals.'))

		k2 <- ratesFirstGuess(inspectIds[eiGenes,], 'processing')
		k3 <- ratesFirstGuess(inspectIds[eiGenes,], 'degradation')

		tau <- as.matrix(sapply(1:ncol(k2), function(j) sapply(1:nrow(k2), function(i) tau_fun(k2[i,j], k3[i,j]) )))
		rownames(tau) <- rownames(k2)

	}

	colnames(tau) <- tpts(inspectIds)
	return(tau)

	})

#' @rdname processingDelay
#' @examples
#' head(calculateDelta(matureInspObj))
#'
setMethod('calculateDelta', signature('INSPEcT'), function(inspectIds, silent=FALSE) {

	eiGenes <- featureNames(inspectIds)[!apply(is.na(ratesFirstGuess(inspectIds, 'preMRNA')),1,all)]

	if( inspectIds@NoNascent ) {

		if( !silent )
			message(paste('No nascent dataset, approximated delta will be calculated on', 
				length(eiGenes),'genes with both intronic and exonic signals.'))

		pre <- ratesFirstGuess(inspectIds[eiGenes,], 'preMRNA')
		delta <- pre/2

	} else {

		if( !silent )
			message(paste('Delta will be calculated on', length(eiGenes),
				'genes with both intronic and exonic signals.'))

		k1 <- ratesFirstGuess(inspectIds[eiGenes,], 'synthesis')
		k2 <- ratesFirstGuess(inspectIds[eiGenes,], 'processing')
		k3 <- ratesFirstGuess(inspectIds[eiGenes,], 'degradation')

		delta <- as.matrix(sapply(1:ncol(k1), function(j) sapply(1:nrow(k1), function(i) delta_fun(k1[i,j], k2[i,j], k3[i,j]) )))
		rownames(delta) <- rownames(k1)

	}

	colnames(delta) <- tpts(inspectIds)
	return(delta)

	})

