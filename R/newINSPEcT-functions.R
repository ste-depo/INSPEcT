#' Create a new INSPEcT object
#'
#' @description
#' The function newINSPEcT creates a new instance of the class INSPEcT provided the experimental time points, RPKMs of total RNA
#' and eventually the 4sU-seq experiments. For the 4sU analysis, it is also requires the labeling time and the scaling factor to
#' normalize the 4sU-seq libraries. This latter parameter can also be calculated by the function itself if both exonic and intronic
#' RPKMs are provided; otherwise it must be given as an input and it is essential to guarantee the robustness of the analysis.
#' The function new INSPEcT also evaluates the variance associated with the RPKMs. It can be done through DESeq2, if the package has
#' been previously used in \code{\link{makeRPKMsFromBams}} or \code{\link{makeRPKMsFromCounts}}, or from the experimental data in a
#' specific experimental condition.
#' @param tpts A vector of time points, one for each sample
#' @param labeling_time A number, lenght of the 4sU pulse
#' @param rpkms_4su_exons A matrix containing expression levels of 4su exons
#' @param rpkms_total_exons A matrix containing expression levels of total exons
#' @param rpkms_4su_introns A matrix containing expression levels of 4su introns
#' @param rpkms_total_introns A matrix containing expression levels of total introns
#' @param BPPARAM Configuration for BiocParallel parallelization. By default is set to bpparam()
#' @param totalMedianNorm A logical to perform median normalization over total RNA exons rpkms, it will apply also on introns
#' @param labeledMedianNorm A logical to perform median normalization over 4sU RNA exons rpkms, it will apply also on introns
#' @param totalSF A vector storing user defined normalization scale over Total RNA exons and introns rpkms
#' @param labeledSF A vector storing user defined normalization scale over 4sU RNA exons and introns rpkms
#' @param totalQuantileNorm A logical to perform to perform median normalization over total RNA exons rpkms, it will apply also on introns
#' @param labeledQuantileNorm A logical to perform to perform median normalization over 4sU RNA exons rpkms, it will apply also on introns
#' @param simulatedData A logical, set to TRUE in case the analysis is on simulated data
#' @param degDuringPulse A logical, set to TRUE in case of a long labelling time. Also degradation of newly synthesized transcripts will be taken into account
#' @param Dmin A numerical, it is the lower bound of the degradation rate domain for the prior optimization 
#' @param Dmax A numerical, it is the upper bound of the degradation rate domain for the prior optimization
#' @param genesFilter, A logical which, if TRUE, filters out genes which have no signal in at least 2/3 of the time points in each feature
#' @param dispersion_parameters_DESeq2, A list of parameters for total and, eventually, 4sU data produced by the DESeq2 analysis and exploited by newINSPEcT
#' for the estimation of the variance in the DESeq2 configuration 
#' @param varSamplingCondition, a character containing the name of the time points selected for the evaluation of the variance if DESeq2 is not used
#' @return An object of class INSPEcT with rates guessed, rates can be accessed by \code{\link{ratesFirstGuess}} method.
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' data('allcounts4su', package='INSPEcT')
#' 
#' tpts <- c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
#' temporalDesign <- rep(tpts,3)
#' tL <- 1/6
#' 
#' #4sU analysis with DESeq2
#' makeRPKMsOut_4sU <- makeRPKMsFromCounts(txdb=txdb,allcounts=allcounts4su,temporalDesign=temporalDesign,DESeq2=TRUE)
#' 
#' rpkms_4sU <- makeRPKMsOut_4sU$rpkms
#' counts_4sU <- makeRPKMsOut_4sU$counts
#' annotations_4sU <- makeRPKMsOut_4sU$annotation
#' dispersion_parameters_DESeq2_4sU <- makeRPKMsOut_4sU$dispersion_parameters_DESeq2
#' 
#' #Analysis with 4sU and DESeq2
#' 
#' mycerIds_4sU <- newINSPEcT(tpts=temporalDesign,labeling_time=tL,rpkms_4su_exons=rpkms_4sU$foursu_exons,rpkms_total_exons=rpkms_4sU$total_exons,rpkms_4su_introns=rpkms_4sU$foursu_introns,rpkms_total_introns=rpkms_4sU$total_introns,dispersion_parameters_DESeq2=dispersion_parameters_DESeq2_4sU,varSamplingCondition=NULL)

newINSPEcT <- function(tpts
					 , labeling_time = NULL
					 , rpkms_4su_exons = NULL
					 , rpkms_total_exons
					 , rpkms_4su_introns = NULL
					 , rpkms_total_introns
					 , BPPARAM = bpparam()
					 , totalMedianNorm = TRUE
					 , labeledMedianNorm = FALSE
					 , totalSF = NULL
					 , labeledSF = NULL
					 , totalQuantileNorm = FALSE
					 , labeledQuantileNorm = FALSE
					 , simulatedData = FALSE
					 , degDuringPulse = FALSE
					 , Dmin = 1e-6
					 , Dmax = 10
					 , genesFilter = TRUE
					 , dispersion_parameters_DESeq2 = NULL
					 , varSamplingCondition = NULL)
{

	if(all(table(tpts)==1))stop("newINSPEcT: at least one replicate is required.")
	speedyVar <- function(x) sum((x - mean.default(x))^2)/length(x[-1])
	countVar2rpkmVar <- function(vars, widths, libsize) vars*(10^9/(widths%o%libsize))^2

	if(is.null(dispersion_parameters_DESeq2)){DESeq2 <- FALSE}else{DESeq2 <- TRUE}
	if(is.null(labeling_time) & is.null(rpkms_4su_exons)){No4sU <- TRUE}else{No4sU <- FALSE}
	if(!DESeq2)
	{
		noRep <- TRUE
		if(is.null(varSamplingCondition)){stop('newINSPEcT: if some time points do not have at least one replica, and the analysis is done without DESeq2, it is mandatory to select a timepoint for the sampling of the variance.')}

		options(warn=-1)
		options(menu.graphics=FALSE)

		# Fit of a powerlaw needed to evaluate the variance of timepoint without replicates
		exprData <- rpkms_total_exons
		pData <- data.frame(feature=tpts)
		colnames(exprData) <- paste0(pData$feature,"_",seq_along(tpts))
		rownames(pData) <- paste0(pData$feature,"_",seq_along(tpts))
		phenoData <- new('AnnotatedDataFrame', data=pData)
		fData <- data.frame(exprData)
		featureData <- new('AnnotatedDataFrame', data=fData)
		foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
		totalFitVariance <- plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,plot.file=FALSE,fittingEval=FALSE, verbose=FALSE)
		totalFitVarianceLaw <- function(expression){(exp(totalFitVariance$INTERCEPT)*expression^(totalFitVariance$SLOPE))^2}

		if( !is.null(rpkms_total_introns) )
		{
			exprData <- rpkms_total_introns
			pData <- data.frame(feature=tpts)
			colnames(exprData) <- paste0(pData$feature,"_",seq_along(tpts))
			rownames(pData) <- paste0(pData$feature,"_",seq_along(tpts))
			phenoData <- new('AnnotatedDataFrame', data=pData)
			fData <- data.frame(exprData)
			featureData <- new('AnnotatedDataFrame', data=fData)
			foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
			preMRNAFitVariance <- plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,plot.file=FALSE,fittingEval=FALSE, verbose=FALSE)
			preMRNAFitVarianceLaw <- function(expression){(exp(preMRNAFitVariance$INTERCEPT)*expression^(preMRNAFitVariance$SLOPE))^2}
		}
		if(!No4sU)
		{
			exprData <- rpkms_4su_exons
			pData <- data.frame(feature=tpts)
			colnames(exprData) <- paste0(pData$feature,"_",seq_along(tpts))
			rownames(pData) <- paste0(pData$feature,"_",seq_along(tpts))
			phenoData <- new('AnnotatedDataFrame', data=pData)
			fData <- data.frame(exprData)
			featureData <- new('AnnotatedDataFrame', data=fData)
			foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
			totalFitVariance_4sU <- plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,plot.file=FALSE,fittingEval=FALSE, verbose=FALSE)
			totalFitVarianceLaw_4sU <- function(expression){(exp(totalFitVariance_4sU$INTERCEPT)*expression^(totalFitVariance_4sU$SLOPE))^2}
	
			if( !is.null(rpkms_4su_introns) )
			{
				exprData <- rpkms_4su_introns
				pData <- data.frame(feature=tpts)
				colnames(exprData) <- paste0(pData$feature,"_",seq_along(tpts))
				rownames(pData) <- paste0(pData$feature,"_",seq_along(tpts))
				phenoData <- new('AnnotatedDataFrame', data=pData)
				fData <- data.frame(exprData)
				featureData <- new('AnnotatedDataFrame', data=fData)
				foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
				preMRNAFitVariance_4sU <- plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,plot.file=FALSE,fittingEval=FALSE, verbose=FALSE)
				preMRNAFitVarianceLaw_4sU <- function(expression){(exp(preMRNAFitVariance_4sU$INTERCEPT)*expression^(preMRNAFitVariance_4sU$SLOPE))^2}
			}
		}

		options(warn=0)

#		dev.off()

	}else{noRep <- FALSE}

	if(No4sU)
	{
		eiGenes <- intersect(rownames(rpkms_total_exons),rownames(rpkms_total_introns))

		rpkms_total_exons <- rpkms_total_exons[eiGenes,]
		rpkms_total_introns <- rpkms_total_introns[eiGenes,]

		if( !is.numeric(Dmin) )
			stop('newINSPEcT: Dmin must be a numeric')
	
		if( !is.numeric(Dmax) )
			stop('newINSPEcT: Dmax must be a numeric')
		
		if( !is.numeric(tpts) )
			stop('newINSPEcT: tpts must be a numeric')
	
		if( class(rpkms_total_exons)=='data.frame' ) 
			rpkms_total_exons <- as.matrix(rpkms_total_exons)
		if( !is.numeric(rpkms_total_exons) )
			stop('newINSPEcT: rpkms_total_exons must be a numeric or a data.frame coercible to numeric')
	
		if( class(rpkms_total_introns)=='data.frame' ) 
			rpkms_total_introns <- as.matrix(rpkms_total_introns)
		if( !is.numeric(rpkms_total_introns) )
			stop('newINSPEcT: rpkms_total_introns must be a numeric or a data.frame coercible to numeric')
		
		if( !class(BPPARAM) %in% sapply(registered(),class) )
			stop('newINSPEcT: BPPARAM argument not registered.')
	
		# check that matrices have rownames
		if( is.null(rownames(rpkms_total_exons))
		  | is.null(rownames(rpkms_total_introns)))
		stop('newINSPEcT: rpkms matrices must have rownames')
		
		##### CALCOLO DELLA VARIANZA CAMPIONARIA
		message('Estimating variance from replicates...')

		if(DESeq2)
		{
			varCountsExons<-t(sapply(rownames(rpkms_total_exons),function(gene){dispersion_parameters_DESeq2$total$muExons[gene,]+dispersion_parameters_DESeq2$total$alphaExons[gene]*dispersion_parameters_DESeq2$total$muExons[gene,]^2}))
			rpkms_total_exons_variances<-t(apply(varCountsExons,1,function(x)tapply(x, tpts, mean)))
  			rpkms_total_exons_variances<-countVar2rpkmVar(rpkms_total_exons_variances
  												,dispersion_parameters_DESeq2$total$exonsWidths[rownames(rpkms_total_exons_variances)]
  												,tapply(dispersion_parameters_DESeq2$total$libsizes,tpts,mean))

  			rownames(rpkms_total_exons_variances) <- rownames(rpkms_total_exons)
  		}else{
			
			rpkms_total_exons_variances <- t(sapply(1:nrow(rpkms_total_exons),function(r){
				sapply(1:ncol(rpkms_total_exons),function(c)
				{
					totalFitVarianceLaw(rpkms_total_exons[r,c])
				})
			}))
			rpkms_total_exons_variances <- t(as.matrix(apply(rpkms_total_exons_variances, 1, function(x) tapply(x, tpts, mean.default))))
		}
		#
		if(DESeq2)
		{
			varCountsIntrons<-t(sapply(rownames(rpkms_total_introns),function(gene){dispersion_parameters_DESeq2$total$muIntrons[gene,]+dispersion_parameters_DESeq2$total$alphaIntrons[gene]*dispersion_parameters_DESeq2$total$muIntrons[gene,]^2}))
			rpkms_total_introns_variances<-t(apply(varCountsIntrons,1,function(x) tapply(x, tpts, mean)))
	
			rpkms_total_introns_variances<-countVar2rpkmVar(rpkms_total_introns_variances
  											,dispersion_parameters_DESeq2$total$intronsWidths[rownames(rpkms_total_introns_variances)]
  											,tapply(dispersion_parameters_DESeq2$total$libsizes,tpts,mean))

			rownames(rpkms_total_introns_variances) <- rownames(rpkms_total_introns)
  		}else{
			rpkms_total_introns_variances <- t(sapply(1:nrow(rpkms_total_introns),function(r){
				sapply(1:ncol(rpkms_total_introns),function(c)
				{
					totalFitVarianceLaw(rpkms_total_introns[r,c])
				})
			}))
			rpkms_total_introns_variances <- t(as.matrix(apply(rpkms_total_introns_variances, 1, function(x) tapply(x, tpts, mean.default))))
		}

		rownames(rpkms_total_introns_variances) <- rownames(rpkms_total_introns)
		rownames(rpkms_total_exons_variances) <- rownames(rpkms_total_exons)

		rpkms_total_exons <- as.matrix(apply(rpkms_total_exons, 1, function(x) 
			tapply(x, tpts, mean.default)))
		if( ncol(rpkms_total_exons)>1 ) rpkms_total_exons <- t(rpkms_total_exons)
		
		rpkms_total_introns <- as.matrix(apply(rpkms_total_introns, 1, function(x) 
			tapply(x, tpts, mean.default)))
		if( ncol(rpkms_total_introns)>1 ) rpkms_total_introns <- t(rpkms_total_introns)
			tpts <- sort(unique(tpts))

		if(genesFilter)
		{
	
			########
			### filter out genes which have no signal in at least 2/3 of the time points 
			### in each feature
			#######
	
			ix1 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
			ix2 <- apply(rpkms_total_exons_variances, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons_variances)
			ix3 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
			ix4 <- apply(rpkms_total_introns_variances, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns_variances)
	
			filteroutGenes <- rownames(rpkms_total_exons)[ix1>2/3 | ix2>2/3 | ix3>2/3 | ix4>2/3 ]
			if( length(filteroutGenes)>0 ) {
				message(
					'Some genes have been filtered out because they have zero valued features in more than 2/3 of the time points in their exonic or intronic features: '
					, paste(filteroutGenes[1:min(10, length(filteroutGenes))], collapse='; ')
					, sep= ''
					)
				if( length(filteroutGenes)>10 )
					message('... and other ', length(filteroutGenes)-10)
				rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons_variances <- rpkms_total_exons_variances[!rownames(rpkms_total_exons_variances) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_introns_variances <- rpkms_total_introns_variances[!rownames(rpkms_total_introns_variances) %in% filteroutGenes, ,drop=FALSE]
			}
	
		}

		# exons and introns rpkms must be provided for every genes, if not: stop
		if( length(setdiff(rownames(rpkms_total_introns), rownames(rpkms_total_exons)))>0 )
			stop('newINSPEcT: for some genes introns RPKMs are provided but not the corresponding exons RPKMs.')
		if( length(setdiff(rownames(rpkms_total_introns), rownames(rpkms_total_exons)))<0 )
			stop('newINSPEcT: for some genes exons RPKMs are provided but not the corresponding introns RPKMs.')
	
		totRpkms <- list()
	
		## subset the "introns and exons genes" from the total
		totRpkmsIntEx <- list(
				exons=rpkms_total_exons
				, exons_var=rpkms_total_exons_variances
				, introns=rpkms_total_introns
				, introns_var=rpkms_total_introns_variances)

		tpts <- sort(unique(tpts))

		## estimate the rates
		outIntEx <- .getRatesAndConcentrationsFromRpkms_No4sU(totRpkmsIntEx
															, tpts
															, BPPARAM=BPPARAM
															, modellingParameters=list(Dmin = Dmin
																					 , Dmax = Dmax
																					 , BPPARAM = BPPARAM)
															, genesFilter = genesFilter)
	
		colnames(outIntEx$concentrations$total_var) <- colnames(outIntEx$concentrations$preMRNA_var) <- NULL
	
		out <- list(concentrations=list(total=outIntEx$concentrations$total
									  , total_var=outIntEx$concentrations$total_var
									  , preMRNA=outIntEx$concentrations$preMRNA
									  , preMRNA_var=outIntEx$concentrations$preMRNA_var )
				  , rates=list(alpha=outIntEx$rates$alpha
							 , alpha_var=outIntEx$rates$alpha_var
				  			 , beta=outIntEx$rates$beta
				  			 , gamma=outIntEx$rates$gamma )
				  , ratesEstimPrec=outIntEx$ratesEstimPrec
				  , geneNames=outIntEx$geneNames
				  , tpts=outIntEx$tpts
				  , labeledSF=outIntEx$labeledSF
				  , totalSF=outIntEx$totalSF
				  , tL=outIntEx$tL)
	
		# make an "ExpressionSet" object containing all the information
		nTpts <- length(tpts)
		exprData <- cbind(out$concentrations$total
						, out$concentrations$preMRNA
						, out$rates$alpha
						, out$rates$beta
						, out$rates$gamma)
	
		pData <- data.frame(feature=c(rep('total',nTpts)
									, rep('preMRNA',nTpts)
									, rep('synthesis',nTpts)
									, rep('degradation',nTpts)
									, rep('processing',nTpts) )
						  , time=rep(tpts, 5))
	
		colnames(exprData) <- paste(pData$feature, signif(pData$time,2), sep='_')
		rownames(pData) <- colnames(exprData)
		phenoData <- new('AnnotatedDataFrame', data=pData)
		
		fData <- data.frame(total=out$concentrations$total_var
						  , preMRNA=out$concentrations$preMRNA_var
						  , synthesis=out$rates$alpha_var
						  , degradation=1
						  , processing=1
						  , synthesis_t0=NA
						  , degradation_t0=NA
						  , processing_t0=NA
						  , nRepT0=NA)

		rownames(fData) <- out$geneNames
	
		featureData <- new('AnnotatedDataFrame', data=fData)
		ratesFirstGuess <- ExpressionSet(assayData=exprData
									   , phenoData=phenoData
									   , featureData=featureData)
	
		# Controls necessary to keep the output of newINSPEcT_No4sU and newINSPEcT equal
	
		if( is.null(out$totalSF) ) out$totalSF <- numeric(0)
		if( is.null(out$labeledSF) ) out$labeledSF <- numeric(0)
		if( is.null(out$tL) ) out$tL <- numeric(0)
	
		## update the object and return it
		object <- new('INSPEcT')
		object@tpts <- tpts
		object@totalSF <- out$totalSF
		object@labeledSF <- out$labeledSF
		object@tL <- out$tL	
		object@ratesFirstGuess <- ratesFirstGuess
		object@precision <- out$ratesEstimPrec
		object@model@simple <- TRUE
		object@ratesFirstGuessP <- outIntEx$ratesFirstGuessP
	
		return(object)

	}else{

		if( class(rpkms_total_exons)=='data.frame' ) 
			rpkms_total_exons <- as.matrix(rpkms_total_exons)
		if( class(rpkms_4su_exons)=='data.frame' ) 
			rpkms_4su_exons <- as.matrix(rpkms_4su_exons)
		if( class(rpkms_4su_exons)=='data.frame' ) 
			rpkms_4su_exons <- as.matrix(rpkms_4su_exons)
	
		## check the input arguments
		if( !is.numeric(tpts) )
			stop('newINSPEcT: tpts must be a numeric')
		if( !is.numeric(labeling_time) )
			stop('newINSPEcT: labeling_time must be a numeric')
	
		if( class(rpkms_4su_exons)=='data.frame' ) 
			rpkms_4su_exons <- as.matrix(rpkms_4su_exons)
		if( !is.numeric(rpkms_4su_exons) )
			stop('newINSPEcT: rpkms_4su_exons must be a numeric or a data.frame coercible to numeric')
	
		if( class(rpkms_total_exons)=='data.frame' ) 
			rpkms_total_exons <- as.matrix(rpkms_total_exons)
		if( !is.numeric(rpkms_total_exons) )
			stop('newINSPEcT: rpkms_total_exons must be a numeric or a data.frame coercible to numeric')
	
		if( !is.null(rpkms_4su_introns) ) {
			if( class(rpkms_4su_introns)=='data.frame' ) 
				rpkms_4su_introns <- as.matrix(rpkms_4su_introns)
			if( !is.numeric(rpkms_4su_introns) )
				stop('newINSPEcT: rpkms_4su_introns must be a numeric or a data.frame coercible to numeric')		
		}
	
		if( !is.null(rpkms_total_introns) ) {
			if( class(rpkms_total_introns)=='data.frame' ) 
				rpkms_total_introns <- as.matrix(rpkms_total_introns)
			if( !is.numeric(rpkms_total_introns) )
				stop('newINSPEcT: rpkms_total_introns must be a numeric or a data.frame coercible to numeric')		
		}
	
		if( !(is.numeric(totalSF) || is.null(totalSF)) )
			stop('newINSPEcT: totalSF must be either numeric or NULL')
		if( !(is.numeric(labeledSF) || is.null(labeledSF)) )
			stop('newINSPEcT: labeledSF must be either numeric or NULL')
	
		if( !is.null(totalSF) )
			if( length(totalSF) != length(unique(tpts)) )
				stop('newINSPEcT: length of totalSF must equal to length(unique(tpts))')
		if( !is.null(labeledSF) )
			if( length(labeledSF) != length(unique(tpts)) )
				stop('newINSPEcT: length of labeledSF must equal to length(unique(tpts))')
	
		if( !class(BPPARAM) %in% sapply(registered(),class) )
			stop('newINSPEcT: BPPARAM argument not registered.')
		if( !is.logical(totalMedianNorm) )
			stop('newINSPEcT: totalMedianNorm must be a logical.')
		if( !is.logical(labeledMedianNorm) )
			stop('newINSPEcT: labeledMedianNorm must be a logical.')
		if( !is.logical(totalQuantileNorm) )
			stop('newINSPEcT: totalQuantileNorm must be a logical.')
		if( !is.logical(labeledQuantileNorm) )
			stop('newINSPEcT: labeledQuantileNorm must be a logical.')
		if( !is.logical(simulatedData) )
			stop('newINSPEcT: simulatedData must be a logical.')
		if( !is.logical(degDuringPulse) )
			stop('newINSPEcT: degDuringPulse must be a logical.')
	
		if( length(tpts) != ncol(rpkms_4su_exons) )
			stop('newINSPEcT: length of tpts is not equal to rpkms matrices coulumns number')

		originalTpts <- tpts

		if( simulatedData ) {

			intExGenes <- rownames(rpkms_4su_exons)
			onlyExGenes <- character(0)
			rpkms_4su_introns <- matrix(0, nrow(rpkms_4su_exons), ncol(rpkms_4su_exons))
			rownames(rpkms_4su_introns) <- rownames(rpkms_4su_exons)
	
		} else {

			## in case either introns of 4sU or introns of total fraction are not
			# provided (both are required) proceed in the simple mode
			if( is.null(rpkms_4su_introns) | is.null(rpkms_total_introns) )
			{
				message('Intronic RPKMs have not been provided, only degradation rates and synthesis rates will be evaluated.')
	
				## in the simple mode the scaling factor between the two library is
				# strongly reccomended
				if( is.null(labeledSF) )
					warning('Without Intronic RPKMs nor labeledSF the analysis could be meaningless.')
	
				## check that matrices have rownames
				if( is.null(rownames(rpkms_4su_exons)) 
					| is.null(rownames(rpkms_total_exons)) )
					stop('newINSPEcT: rpkms matrices must have rownames')
	
				## check that the rownames correspond
				if( !identical(rownames(rpkms_4su_exons), rownames(rpkms_total_exons)) )
					stop('newINSPEcT: rpkms matrices rownames do not correspond')		
	
				## check that matrices have same dimensions
				if( !identical(dim(rpkms_4su_exons), dim(rpkms_total_exons)) )
					stop("newINSPEcT: input matrices don't have same dimensions")
	
				########
				### filter out genes which have no signal in at least 2/3 of the time points 
				### in each feature
				#######
	
				ix1 <- apply(rpkms_4su_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_4su_exons)
				ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
				filteroutGenes <- rownames(rpkms_4su_exons)[ix1>2/3 | ix2>2/3]
	
				if( length(filteroutGenes)>0 ) {
					message(
						'Some genes have been filtered out because they have zero valued features in more than 2/3 of the time points: '
						, paste(filteroutGenes, collapse='; ')
						, sep= ''
						)
					rpkms_4su_exons <- rpkms_4su_exons[!rownames(rpkms_4su_exons) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				}
	
				## assign all the genes to the "only exons" mode.
				intExGenes <- character(0)
				onlyExGenes <- rownames(rpkms_4su_exons)
	
			} else { # analysis with introns
				## check that matrices have rownames
				if( is.null(rownames(rpkms_4su_exons)) 
					| is.null(rownames(rpkms_total_exons)) 
					| is.null(rownames(rpkms_4su_introns)) 
					| is.null(rownames(rpkms_total_introns)) )
					stop('newINSPEcT: rpkms matrices must have rownames')

				## check that corresponding matrices have same rownames (exons with exons, introns with introns)
				if( !identical(rownames(rpkms_4su_exons), rownames(rpkms_total_exons)) )
					stop('newINSPEcT: rpkms_4su_exons and rpkms_total_exons rownames does not correspond.')
				if( !identical(rownames(rpkms_4su_introns), rownames(rpkms_total_introns)) )
					stop('newINSPEcT: rpkms_4su_introns and rpkms_total_introns rownames does not correspond.')

				## check that corresponding matrices have same dimensions
				if( !identical(dim(rpkms_4su_exons), dim(rpkms_total_exons)) )
					stop("newINSPEcT: rpkms exons matrices don't have same dimensions")
				if( !identical(dim(rpkms_4su_introns), dim(rpkms_total_introns)) )
					stop("newINSPEcT: rpkms introns matrices don't have same dimensions")

				## in case, apply quantile normalization
				if( totalQuantileNorm ) {
					normMat <- normalize.quantiles(rbind(rpkms_total_exons, rpkms_total_introns))
					exonFeatures <- rownames(rpkms_total_exons)
					rpkms_total_exons <- normMat[1:nrow(rpkms_total_exons),]
					rownames(rpkms_total_exons) <- exonFeatures
					intronFeatures <- rownames(rpkms_total_introns)
					rpkms_total_introns <- normMat[(nrow(rpkms_total_exons)+1):nrow(normMat),]
					rownames(rpkms_total_introns) <- intronFeatures
				}
				if( labeledQuantileNorm ) {
					normMat <- normalize.quantiles(rbind(rpkms_4su_exons, rpkms_4su_introns))
					exonFeatures <- rownames(rpkms_4su_exons)
					rpkms_4su_exons <- normMat[1:nrow(rpkms_4su_exons),]
					rownames(rpkms_4su_exons) <- exonFeatures
					intronFeatures <- rownames(rpkms_4su_introns)
					rpkms_4su_introns <- normMat[(nrow(rpkms_4su_exons)+1):nrow(normMat),]
					rownames(rpkms_4su_introns) <- intronFeatures
				}

				########
				### filter out genes which have no signal in at least 2/3 of the time points 
				### in each feature
				#######
	
				ix1 <- apply(rpkms_4su_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_4su_exons)
				ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
				filteroutGenes <- rownames(rpkms_4su_exons)[ix1>2/3 | ix2>2/3] # | ix3>2/3 | ix4>2/3]
				if( length(filteroutGenes)>0 ) {
					message(
						'Some genes have been filtered out because they have zero valued features in more than 2/3 of the time points in their exonic features: '
						, paste(filteroutGenes[1:min(10, length(filteroutGenes))], collapse='; ')
						, sep= ''
						)
					if( length(filteroutGenes)>10 )
						message('... and other ', length(filteroutGenes)-10)
					rpkms_4su_exons <- rpkms_4su_exons[!rownames(rpkms_4su_exons) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
					rpkms_4su_introns <- rpkms_4su_introns[!rownames(rpkms_4su_introns) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
				}
	
				ix3 <- apply(rpkms_4su_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_4su_introns)
				ix4 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
				filteroutGenes <- rownames(rpkms_4su_exons)[ix3>2/3 | ix4>2/3] # | ix3>2/3 | ix4>2/3]
				if( length(filteroutGenes)>0 ) {
					message(
						'For some genes only synthesis and degradation will be evaluated because they have zero valued features in more than 2/3 of the time points in their intronic features: '
						, paste(filteroutGenes[1:min(10, length(filteroutGenes))], collapse='; ')
						, sep= ''
						)
					if( length(filteroutGenes)>10 )
						message('... and other ', length(filteroutGenes)-10)
					rpkms_4su_introns <- rpkms_4su_introns[!rownames(rpkms_4su_introns) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
				}
	
				## assign genes to "only exons" or "introns and exons" wheter they have introns or not
				if( identical(rownames(rpkms_4su_exons), rownames(rpkms_4su_introns)) ) {
					intExGenes <- rownames(rpkms_4su_exons)
					onlyExGenes <- character(0)
				} else {
					## exons rpkms must be provided for every genes, if not: stop
					if( length(setdiff(rownames(rpkms_4su_introns), rownames(rpkms_4su_exons)))>0 )
						stop('newINSPEcT: for some genes introns RPKMs are provided but not the corresponding exons RPKMs.')
					intExGenes <- intersect(rownames(rpkms_4su_exons), rownames(rpkms_4su_introns))
					onlyExGenes <- setdiff(rownames(rpkms_4su_exons), rownames(rpkms_4su_introns))
					message('Some genes have only exons RPKMs, on them only synthesis and degradation will be evaluated.')
					
				}
			}
		}

		totRpkms <- list()
		labeledRpkms <- list()

		##### CALCOLO DELLA VARIANZA CAMPIONARIA
		message('Estimating variance from replicates...')

# browser()

		#EXONS TOTAL
		if(DESeq2)
		{
			varCountsExons<-t(sapply(rownames(rpkms_total_exons),function(gene){dispersion_parameters_DESeq2$total$muExons[gene,]+dispersion_parameters_DESeq2$total$alphaExons[gene]*dispersion_parameters_DESeq2$total$muExons[gene,]^2}))
			totRpkms$exons_var<-t(apply(varCountsExons,1,function(x)tapply(x, tpts, mean)))
		  	totRpkms$exons_var<-countVar2rpkmVar(totRpkms$exons_var
		  										  ,dispersion_parameters_DESeq2$total$exonsWidths[rownames(totRpkms$exons_var)]
		  										  ,tapply(dispersion_parameters_DESeq2$total$libsizes,tpts,mean))
		
			rownames(totRpkms$exons_var) <- rownames(rpkms_total_exons)
		}else{
			totRpkms$exons_var <- t(sapply(1:nrow(rpkms_total_exons),function(r){
				sapply(1:ncol(rpkms_total_exons),function(c)
				{
					totalFitVarianceLaw(rpkms_total_exons[r,c])
				})
			}))
			totRpkms$exons_var <- t(as.matrix(apply(totRpkms$exons_var, 1, function(x) tapply(x, tpts, mean.default))))

		}
		#
		rownames(totRpkms$exons_var) <- rownames(rpkms_total_exons)
		totRpkms$exons <- as.matrix(apply(rpkms_total_exons, 1, function(x) tapply(x, tpts, mean.default)))
		if( ncol(totRpkms$exons)>1 ) totRpkms$exons <- t(totRpkms$exons)
		
		#INTRONS TOTAL
		if( !is.null(rpkms_total_introns) )
		{
			if(DESeq2)
			{
				varCountsIntrons<-t(sapply(rownames(rpkms_total_introns),function(gene){dispersion_parameters_DESeq2$total$muIntrons[gene,]+dispersion_parameters_DESeq2$total$alphaIntrons[gene]*dispersion_parameters_DESeq2$total$muIntrons[gene,]^2}))
				totRpkms$introns_var<-t(apply(varCountsIntrons,1,function(x)tapply(x, tpts, mean)))
			  	totRpkms$introns_var<-countVar2rpkmVar(totRpkms$introns_var
			  										  ,dispersion_parameters_DESeq2$total$intronsWidths[rownames(totRpkms$introns_var)]
			  										  ,tapply(dispersion_parameters_DESeq2$total$libsizes,tpts,mean))
			
				rownames(totRpkms$introns_var) <- rownames(rpkms_total_introns)
			}else{
				totRpkms$introns_var <- t(sapply(1:nrow(rpkms_total_introns),function(r){
					sapply(1:ncol(rpkms_total_introns),function(c)
					{
						totalFitVarianceLaw(rpkms_total_introns[r,c])
					})
				}))
				totRpkms$introns_var <- t(as.matrix(apply(totRpkms$introns_var, 1, function(x) tapply(x, tpts, mean.default))))

			}
			#
			rownames(totRpkms$introns_var) <- rownames(rpkms_total_introns)
		
			totRpkms$introns <- as.matrix(apply(rpkms_total_introns, 1, function(x) tapply(x, tpts, mean.default)))
			if( ncol(totRpkms$introns)>1 ) totRpkms$introns <- t(totRpkms$introns)
		}
		
		#EXONS LABELED
		if(DESeq2)
		{
		
			varCountsExons<-t(sapply(rownames(rpkms_4su_exons),function(gene){dispersion_parameters_DESeq2$foursu$muExons[gene,]+dispersion_parameters_DESeq2$foursu$alphaExons[gene]*dispersion_parameters_DESeq2$foursu$muExons[gene,]^2}))
			labeledRpkms$exons_var<-t(apply(varCountsExons,1,function(x)tapply(x, tpts, mean)))
		  	labeledRpkms$exons_var<-countVar2rpkmVar(labeledRpkms$exons_var
		  										  ,dispersion_parameters_DESeq2$foursu$exonsWidths[rownames(labeledRpkms$exons_var)]
		  										  ,tapply(dispersion_parameters_DESeq2$foursu$libsizes,tpts,mean))
		
			rownames(labeledRpkms$exons_var) <- rownames(rpkms_4su_exons)
		}else{
			labeledRpkms$exons_var <- t(sapply(1:nrow(rpkms_4su_exons),function(r){
				sapply(1:ncol(rpkms_4su_exons),function(c)
				{
					totalFitVarianceLaw(rpkms_4su_exons[r,c])
				})
			}))
			labeledRpkms$exons_var <- t(as.matrix(apply(labeledRpkms$exons_var, 1, function(x) tapply(x, tpts, mean.default))))

		}
		#
		rownames(labeledRpkms$exons_var) <- rownames(rpkms_4su_exons)
		labeledRpkms$exons <- as.matrix(apply(rpkms_4su_exons, 1, function(x) tapply(x, tpts, mean.default)))
		if( ncol(labeledRpkms$exons)>1 ) labeledRpkms$exons <- t(labeledRpkms$exons)
		
		#INTRONS LABELED
		if( !is.null(rpkms_4su_introns) )
		{
			if(DESeq2)
			{
			
				varCountsIntrons<-t(sapply(rownames(rpkms_4su_introns),function(gene){dispersion_parameters_DESeq2$foursu$muIntrons[gene,]+dispersion_parameters_DESeq2$foursu$alphaIntrons[gene]*dispersion_parameters_DESeq2$foursu$muIntrons[gene,]^2}))
				labeledRpkms$introns_var<-t(apply(varCountsIntrons,1,function(x)tapply(x, tpts, mean)))
			  	labeledRpkms$introns_var<-countVar2rpkmVar(labeledRpkms$introns_var
			  										  ,dispersion_parameters_DESeq2$foursu$intronsWidths[rownames(labeledRpkms$introns_var)]
			  										  ,tapply(dispersion_parameters_DESeq2$foursu$libsizes,tpts,mean))
			
				rownames(labeledRpkms$introns_var) <- rownames(rpkms_4su_introns)
			}else{
				labeledRpkms$introns_var <- t(sapply(1:nrow(rpkms_4su_introns),function(r){
					sapply(1:ncol(rpkms_4su_introns),function(c)
					{
						totalFitVarianceLaw(rpkms_4su_introns[r,c])
					})
				}))
				labeledRpkms$introns_var <- t(as.matrix(apply(labeledRpkms$introns_var, 1, function(x) tapply(x, tpts, mean.default))))
			}
			#
			rownames(labeledRpkms$introns_var) <- rownames(rpkms_4su_introns)
			labeledRpkms$introns <- as.matrix(apply(rpkms_4su_introns, 1, function(x) tapply(x, tpts, mean.default)))
			if( ncol(labeledRpkms$introns)>1 ) labeledRpkms$introns <- t(labeledRpkms$introns)
		}

		tpts <- sort(unique(tpts))
	
		if( length(intExGenes)>0 ) {
			message(paste('Number of genes with introns and exons: ', length(intExGenes)))
			## subset the "introns and exons genes" from the total
			totRpkmsIntEx <- list(
					exons=totRpkms$exons[intExGenes, , drop=FALSE]
					, exons_var=totRpkms$exons_var[intExGenes, , drop=FALSE]
					, introns=totRpkms$introns[intExGenes, , drop=FALSE]
					, introns_var=totRpkms$introns_var[intExGenes, , drop=FALSE]
			)
			labeledRpkmsIntEx <- list(
					exons=labeledRpkms$exons[intExGenes, , drop=FALSE]
					, exons_var=labeledRpkms$exons_var[intExGenes, , drop=FALSE]
					, introns=labeledRpkms$introns[intExGenes, , drop=FALSE]
					, introns_var=labeledRpkms$introns_var[intExGenes, , drop=FALSE]
				)

			## estimate the rates
			outIntEx <- .getRatesAndConcentrationsFromRpkms(totRpkmsIntEx
														  , labeledRpkmsIntEx
														  , tpts
														  , tL=labeling_time
														  , simulatedData=simulatedData
														  , BPPARAM=BPPARAM
														  , totalMedianNorm=totalMedianNorm
														  , labeledMedianNorm=labeledMedianNorm
														  , totalSF=totalSF
														  , labeledSF=labeledSF
														  , degDuringPulse=degDuringPulse)
			## set the labeledSF and totalSF, they can be used by "only exons genes" (if present)
			labeledSF <- outIntEx$labeledSF
			totalSF <- outIntEx$totalSF
			totalMedianNorm <- FALSE
			labeledMedianNorm <- FALSE
		}

		if( length(onlyExGenes)>0 ) {
			message(paste('Number of genes with only exons: ', length(onlyExGenes)))
			## subset the "only exons genes" from the total
			totRpkmsOnlyEx <- list(
					exons=totRpkms$exons[onlyExGenes, , drop=FALSE]
					, exons_var=totRpkms$exons_var[onlyExGenes, , drop=FALSE]
				)
			labeledRpkmsOnlyEx <- list(
					exons=labeledRpkms$exons[onlyExGenes, , drop=FALSE]
					, exons_var=labeledRpkms$exons_var[onlyExGenes, , drop=FALSE]
				)
			## estimate the rates, eventually using "labeledSF" and "totalSF" calculated from "introns and exons genes"
			outOnlyEx <- .getRatesAndConcentrationsFromRpkms(totRpkmsOnlyEx
														   , labeledRpkmsOnlyEx
														   , tpts
				, tL=labeling_time, simulatedData=simulatedData, BPPARAM=BPPARAM
				, totalMedianNorm=totalMedianNorm, labeledMedianNorm=labeledMedianNorm
				, totalSF=totalSF, labeledSF=labeledSF, degDuringPulse=degDuringPulse)
		}

		## merge the output of "introns and exons genes" and "only exons genes" (if both present)
		if( length(intExGenes)>0 & length(onlyExGenes)>0 ) {
			## internal control, these data must be identical in the two cathegories
			## if not: something unpredicted happened
			stopifnot(identical(outIntEx$tpts, outOnlyEx$tpts))
			stopifnot(identical(outIntEx$labeledSF, outOnlyEx$labeledSF))
			stopifnot(identical(outIntEx$totalSF, outOnlyEx$totalSF))
			stopifnot(identical(outIntEx$tL, outOnlyEx$tL))
			## merge the two lists

			out <- list(
				concentrations=list(
					total=rbind(outIntEx$concentrations$total, outOnlyEx$concentrations$total)
					#, total_var=c(outIntEx$concentrations$total_var, outOnlyEx$concentrations$total_var)
					, total_var=rbind(outIntEx$concentrations$total_var, outOnlyEx$concentrations$total_var)
					, preMRNA=rbind(outIntEx$concentrations$preMRNA, outOnlyEx$concentrations$preMRNA)
					#, preMRNA_var=c(outIntEx$concentrations$preMRNA_var, outOnlyEx$concentrations$preMRNA_var)
					, preMRNA_var=rbind(outIntEx$concentrations$preMRNA_var, outOnlyEx$concentrations$preMRNA_var)
					)
				, rates=list(
					alpha=rbind(outIntEx$rates$alpha, outOnlyEx$rates$alpha)
					, alpha_var=rbind(outIntEx$rates$alpha_var, outOnlyEx$rates$alpha_var)
					, beta=rbind(outIntEx$rates$beta, outOnlyEx$rates$beta)
					, gamma=rbind(outIntEx$rates$gamma, outOnlyEx$rates$gamma)
					)
				, ratesEstimPrec=rbind(outIntEx$ratesEstimPrec, outOnlyEx$ratesEstimPrec)
				, geneNames=c(outIntEx$geneNames, outOnlyEx$geneNames)
				## from here both lists have the same information, pick it up from the first one
				, tpts=outIntEx$tpts
				, labeledSF=outIntEx$labeledSF
				, totalSF=outIntEx$totalSF
				, tL=outIntEx$tL
				)
		} else if( length(intExGenes)>0 ) {
			out <- outIntEx
		} else {
			out <- outOnlyEx
		}

		## if replicates are available at time zero, compute rate variances
		t0ix <- originalTpts==originalTpts[1]
		nRepT0 <- length(which(t0ix))
		evaluateRateVarAtT0 <- nRepT0>1 # & !simulatedData
		if( evaluateRateVarAtT0 ) {
	
			message('Estimating rate variances at time zero...')
	
			a <- matrix(NA, ncol=nRepT0, nrow=length(out$geneNames))
			rownames(a) <- out$geneNames
			p <- t <- a
			a[rownames(rpkms_4su_exons),] <- rpkms_4su_exons[,t0ix]/out$tL*out$labeledSF[1]
			t[rownames(rpkms_total_exons),] <- rpkms_total_exons[,t0ix]*out$totalSF[1]
			p[rownames(rpkms_total_introns),] <- rpkms_total_introns[,t0ix]*out$totalSF[1]
			mean_a <- apply(a, 1, mean, na.rm=TRUE)
			mean_t <- apply(t, 1, mean, na.rm=TRUE)
			mean_p <- apply(p, 1, mean, na.rm=TRUE)
			mean_1_over_p <- apply(1/p, 1, mean, na.rm=TRUE)
			mean_t_minus_p <- apply(t-p, 1, mean, na.rm=TRUE)
			mean_1_over_t <- apply(1/t, 1, mean, na.rm=TRUE)
			mean_1_over_t_minus_p <- apply(1/(t-p), 1, mean, na.rm=TRUE)
	
			matCov <- function(x, y)
				sapply(1:nrow(x), function(i) 
					tryCatch(cov(x[i,], y[i,], use="complete.obs"), error=function(e) NA))	
	
			var_a <- apply(a, 1, stats::var, na.rm=TRUE)
			var_t <- apply(t, 1, stats::var, na.rm=TRUE)
			var_p <- apply(p, 1, stats::var, na.rm=TRUE)
			var_t_minus_p <- var_t + var_p - 2*matCov(t, p)
			var_1_over_p <- apply(1/p, 1, stats::var, na.rm=TRUE)
			var_c <- var_a_over_p <- matCov(a^2,1/p^2) + ( var_a + mean_a^2 )*( var_1_over_p + mean_1_over_p^2 ) -
				( matCov(a, 1/p) + mean_a*mean_1_over_p )^2
			var_1_over_t <- apply(1/t, 1, stats::var, na.rm=TRUE)
			var_1_over_t_minus_p <- apply(1/(t-p), 1, stats::var, na.rm=TRUE)
			var_a_over_t_minus_p <- matCov(a^2, 1/(t-p)^2) + ( var_a + mean_a^2 )*( var_1_over_t_minus_p + mean_1_over_t_minus_p^2 ) -
				( matCov(a, 1/(t-p)) + mean_a*mean_1_over_t_minus_p )^2
			var_a_over_t <- matCov(a^2, 1/t^2) + ( var_a + mean_a^2 )*( var_1_over_t + mean_1_over_t^2 ) -
				( matCov(a, 1/t) + mean_a*mean_1_over_t )^2
			var_b <- ifelse(is.na(mean_p), var_a_over_t, var_a_over_t_minus_p)
			# mask (non-sensed) negative variances
			var_b[var_b<0] <- NaN
			var_c[var_c<0] <- NaN
	
		}		

		## make an "ExpressionSet" object containing all the information
		nTpts <- length(tpts)
		exprData <- cbind(out$concentrations$total, out$concentrations$preMRNA
			, out$rates$alpha, out$rates$beta, out$rates$gamma)
		pData <- data.frame(
			feature=c(
				rep('total',nTpts)
				, rep('preMRNA',nTpts)
				, rep('synthesis',nTpts)
				, rep('degradation',nTpts)
				, rep('processing',nTpts)
				)
			, time=rep(tpts, 5))
		colnames(exprData) <- paste(pData$feature, 
			signif(pData$time,2), sep='_')
		rownames(pData) <- colnames(exprData)
		phenoData <- new('AnnotatedDataFrame', data=pData)
		fData <- data.frame(
			total=out$concentrations$total_var
			, preMRNA=out$concentrations$preMRNA_var
			, synthesis=out$rates$alpha_var
			, degradation=1
			, processing=1
			, synthesis_t0=if( evaluateRateVarAtT0 ) var_a else NA
			, degradation_t0=if( evaluateRateVarAtT0 ) var_b else NA
			, processing_t0=if( evaluateRateVarAtT0 ) var_c else NA
			, nRepT0=if( evaluateRateVarAtT0 ) nRepT0 else NA
			)
		rownames(fData) <- out$geneNames
		## adjust variance of "gamma" of "onlyExGenes", assigning NA
		if( length(onlyExGenes)>0 ) fData[onlyExGenes,]$processing <- NA
		featureData <- new('AnnotatedDataFrame', data=fData)
		ratesFirstGuess <- ExpressionSet(
			assayData=exprData
			, phenoData=phenoData
			, featureData=featureData
			)
	
		if( is.null(out$totalSF) ) out$totalSF <- numeric(0)
		if( is.null(out$labeledSF) ) out$labeledSF <- numeric(0)
		if( is.null(out$tL) ) out$tL <- numeric(0)
	
		## update the object and return it
		object <- new('INSPEcT')
		object@tpts <- tpts
		object@totalSF <- out$totalSF
		object@labeledSF <- out$labeledSF
		object@tL <- out$tL	
		object@ratesFirstGuess <- ratesFirstGuess
		object@precision <- out$ratesEstimPrec
		object@params$No4sU <- No4sU
		if( is.null(rpkms_4su_introns) | is.null(rpkms_total_introns) )
			object@model@simple <- TRUE
	
		return(object)
	}
}

