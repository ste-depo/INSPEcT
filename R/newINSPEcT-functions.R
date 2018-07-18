#' Create a new INSPEcT object
#'
#' @description
#' The function newINSPEcT creates a new instance of the class INSPEcT provided the experimental time points, expression data (like RPKMs) of mature
#' and eventually nascent RNA. For the nascent analysis, it is also requires a collecting time and the scaling factor to
#' normalize the nascent RNA-seq libraries. This latter parameter can also be calculated by the function itself if both exonic and intronic
#' expression data are provided; otherwise it must be given as an input and it is essential to guarantee the robustness of the analysis.
#' @param tpts A vector of time points, one for each sample
#' @param labeling_time A number, lenght of the Nascent pulse
#' @param nascentExpressions A list which contains exons and introns expression matrices and variances for the nascent RNA
#' @param matureExpressions A list which contains exons and introns expression matrices and variances for the mature RNA
#' @param BPPARAM Configuration for BiocParallel parallelization. By default is set to bpparam()
#' @param totalMedianNorm A logical to perform median normalization over total RNA exons rpkms, it will apply also on introns
#' @param labeledMedianNorm A logical to perform median normalization over Nascent RNA exons rpkms, it will apply also on introns
#' @param totalSF A vector storing user defined normalization scale over Total RNA exons and introns rpkms
#' @param labeledSF A vector storing user defined normalization scale over Nascent RNA exons and introns rpkms
#' @param totalQuantileNorm A logical to perform to perform median normalization over total RNA exons rpkms, it will apply also on introns
#' @param labeledQuantileNorm A logical to perform to perform median normalization over Nascent RNA exons rpkms, it will apply also on introns
#' @param simulatedData A logical, set to TRUE in case the analysis is on simulated data
#' @param degDuringPulse A logical, set to TRUE in case of a long labelling time. Also degradation of newly synthesized transcripts will be taken into account
#' @param Dmin A numerical, it is the lower bound of the degradation rate domain for the prior optimization 
#' @param Dmax A numerical, it is the upper bound of the degradation rate domain for the prior optimization
#' @param genesFilter, A logical which, if TRUE, filters out genes which have no signal in at least 2/3 of the time points in each feature
#' @return An object of class INSPEcT with a first estimation of the rates which can be accessed by the method \code{\link{ratesFirstGuess}}.
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' data('allcounts', package='INSPEcT')
#' 
#' matureCounts<-allcounts$mature
#' 
#' testGenes<-rownames(matureCounts$exonsCounts)
#' 
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb<-TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' exonsDB<-reduce(exonsBy(txdb ,'gene'))
#' exonsDB<-exonsDB[elementNROWS(range(exonsDB))==1]
#' intronsDB<-psetdiff(unlist(range(exonsDB)),exonsDB)
#' intronsDB<-intronsDB[elementNROWS(intronsDB)>0]
#' 
#' exWdths<-sapply(width(exonsDB),sum)
#' intWdths<-sapply(width(intronsDB),sum)
#' 
#' exWdths<-exWdths[testGenes]
#' intWdths<-intWdths[testGenes]
#'     
#' totalLS<-colSums(matureCounts$
#'   stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
#' 
#' tpts<-c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
#' expDes<-rep(tpts,3)
#'     
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(libsize=totalLS
#'                                               ,exonsWidths=exWdths
#'                                               ,intronsWidths=intWdths
#'                                               ,allcounts=matureCounts
#'                                               ,DESeq2=TRUE
#'                                               ,experimentalDesign=expDes)
#' 
#' matureInspObj<-newINSPEcT(tpts=tpts
#'                          ,labeling_time=NULL
#'                          ,nascentExpressions=NULL
#'                          ,matureExpressions=matExp_DESeq2)

newINSPEcT <- function(tpts
					 , labeling_time = NULL
					 , nascentExpressions = NULL
					 , matureExpressions
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
					 , genesFilter = TRUE)
{

	if(is.null(labeling_time) & is.null(nascentExpressions)){NoNascent <- TRUE}
	else if(!is.null(labeling_time) & !is.null(nascentExpressions)){NoNascent <- FALSE}
	else(stop('newINSPEcT: labeling_time and nascentExpressions must be both NULL or defined'))

	rpkms_total_exons <- matureExpressions$exonsExpressions
	rpkms_total_exons_variances <- matureExpressions$exonsVariance
	rpkms_total_introns <- matureExpressions$intronsExpressions
	rpkms_total_introns_variances <- matureExpressions$intronsVariance

	if(NoNascent)
	{
		#Only genes with exons and introns can be modeled without nascent RNA 
		eiGenes <- intersect(rownames(rpkms_total_exons),rownames(rpkms_total_introns))

		rpkms_total_exons <- rpkms_total_exons[eiGenes,]
		rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes,]
		rpkms_total_introns <- rpkms_total_introns[eiGenes,]
		rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes,]

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
		outIntEx <- .getRatesAndConcentrationsFromRpkms_NoNascent(totRpkmsIntEx
															, tpts
															, BPPARAM=BPPARAM
															, modellingParameters=list(Dmin = Dmin
																					 , Dmax = Dmax)
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
	
		# Controls necessary to keep the output of newINSPEcT_NoNascent and newINSPEcT equal
	
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
	
	}else{ #Nascent RNA mode

		rpkms_Nascent_exons <- nascentExpressions$exonsExpressions
		rpkms_Nascent_exons_variances <- nascentExpressions$exonsVariance
		rpkms_Nascent_introns <- nascentExpressions$intronsExpressions
		rpkms_Nascent_introns_variances <- nascentExpressions$intronsVariance

		originalTpts <- tpts

		if( simulatedData ) {

			intExGenes <- rownames(rpkms_Nascent_exons)
			onlyExGenes <- character(0)
			rpkms_Nascent_introns <- matrix(0, nrow(rpkms_Nascent_exons), ncol(rpkms_Nascent_exons))
			rownames(rpkms_Nascent_introns) <- rownames(rpkms_Nascent_exons)
	
		} else {

			## in case either introns of Nascent or introns of total fraction are not
			# provided (both are required) proceed in the simple mode
			if( is.null(rpkms_Nascent_introns) | is.null(rpkms_total_introns) )
			{
				message('Intronic RPKMs have not been provided, only degradation rates and synthesis rates will be evaluated.')
	
				## in the simple mode the scaling factor between the two library is
				# strongly reccomended
				if( is.null(labeledSF) )
					warning('Without Intronic RPKMs nor labeledSF the analysis could be meaningless.')
	
				## check that matrices have rownames
				if( is.null(rownames(rpkms_Nascent_exons)) 
					| is.null(rownames(rpkms_total_exons)) )
					stop('newINSPEcT: rpkms matrices must have rownames')
	
				## check that the rownames correspond
				if( !identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons)) )
					stop('newINSPEcT: rpkms matrices rownames do not correspond')		
	
				## check that matrices have same dimensions
				if( !identical(dim(rpkms_Nascent_exons), dim(rpkms_total_exons)) )
					stop("newINSPEcT: input matrices don't have same dimensions")
	
				########
				### filter out genes which have no signal in at least 2/3 of the time points 
				### in each feature
				#######
	
				ix1 <- apply(rpkms_Nascent_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_exons)
				ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
				filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1>2/3 | ix2>2/3]
	
				if( length(filteroutGenes)>0 ) {
					message(
						'Some genes have been filtered out because they have zero valued features in more than 2/3 of the time points: '
						, paste(filteroutGenes, collapse='; ')
						, sep= ''
						)
					rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				}
	
				## assign all the genes to the "only exons" mode.
				intExGenes <- character(0)
				onlyExGenes <- rownames(rpkms_Nascent_exons)
	
			} else {

				# analysis with introns
				## check that matrices have rownames
				if( is.null(rownames(rpkms_Nascent_exons)) 
					| is.null(rownames(rpkms_total_exons)) 
					| is.null(rownames(rpkms_Nascent_introns)) 
					| is.null(rownames(rpkms_total_introns)) )
				{stop('newINSPEcT: rpkms matrices must have rownames')}

				## check that corresponding matrices have same rownames (exons with exons, introns with introns)
				if( !identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons)) )
					stop('newINSPEcT: rpkms_Nascent_exons and rpkms_total_exons rownames does not correspond.')
				if( !identical(rownames(rpkms_Nascent_introns), rownames(rpkms_total_introns)) )
					stop('newINSPEcT: rpkms_Nascent_introns and rpkms_total_introns rownames does not correspond.')

				## check that corresponding matrices have same dimensions
				if( !identical(dim(rpkms_Nascent_exons), dim(rpkms_total_exons)) )
					stop("newINSPEcT: rpkms exons matrices don't have same dimensions")
				if( !identical(dim(rpkms_Nascent_introns), dim(rpkms_total_introns)) )
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
					normMat <- normalize.quantiles(rbind(rpkms_Nascent_exons, rpkms_Nascent_introns))
					exonFeatures <- rownames(rpkms_Nascent_exons)
					rpkms_Nascent_exons <- normMat[1:nrow(rpkms_Nascent_exons),]
					rownames(rpkms_Nascent_exons) <- exonFeatures
					intronFeatures <- rownames(rpkms_Nascent_introns)
					rpkms_Nascent_introns <- normMat[(nrow(rpkms_Nascent_exons)+1):nrow(normMat),]
					rownames(rpkms_Nascent_introns) <- intronFeatures
				}

				########
				### filter out genes which have no signal in at least 2/3 of the time points 
				### in each feature
				#######
	
				ix1 <- apply(rpkms_Nascent_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_exons)
				ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
				filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1>2/3 | ix2>2/3] # | ix3>2/3 | ix4>2/3]
				if( length(filteroutGenes)>0 ) {
					message(
						'Some genes have been filtered out because they have zero valued features in more than 2/3 of the time points in their exonic features: '
						, paste(filteroutGenes[1:min(10, length(filteroutGenes))], collapse='; ')
						, sep= ''
						)
					if( length(filteroutGenes)>10 )
						message('... and other ', length(filteroutGenes)-10)
					rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
					rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
				}
	
				ix3 <- apply(rpkms_Nascent_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_introns)
				ix4 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
				filteroutGenes <- rownames(rpkms_Nascent_exons)[ix3>2/3 | ix4>2/3] # | ix3>2/3 | ix4>2/3]
				if( length(filteroutGenes)>0 ) {
					message(
						'For some genes only synthesis and degradation will be evaluated because they have zero valued features in more than 2/3 of the time points in their intronic features: '
						, paste(filteroutGenes[1:min(10, length(filteroutGenes))], collapse='; ')
						, sep= ''
						)
					if( length(filteroutGenes)>10 )
						message('... and other ', length(filteroutGenes)-10)
					rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
					rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
				}
	
				## assign genes to "only exons" or "introns and exons" wheter they have introns or not
				if( identical(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_introns)) ) {
					intExGenes <- rownames(rpkms_Nascent_exons)
					onlyExGenes <- character(0)
				} else {
					## exons rpkms must be provided for every genes, if not: stop
					if( length(setdiff(rownames(rpkms_Nascent_introns), rownames(rpkms_Nascent_exons)))>0 )
						stop('newINSPEcT: for some genes introns RPKMs are provided but not the corresponding exons RPKMs.')
					intExGenes <- intersect(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_introns))
					onlyExGenes <- setdiff(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_introns))
					message('Some genes have only exons RPKMs, on them only synthesis and degradation will be evaluated.')
					
				}
			}
		
		}

		totRpkms <- list()
		labeledRpkms <- list()
	
		tpts <- sort(unique(tpts))

		if( length(intExGenes)>0 ) {
			message(paste('Number of genes with introns and exons: ', length(intExGenes)))
			## subset the "introns and exons genes" from the total
			totRpkmsIntEx <- list(
					exons=rpkms_total_exons[intExGenes, , drop=FALSE]
					, exons_var=rpkms_total_exons_variances[intExGenes, , drop=FALSE]
					, introns=rpkms_total_introns[intExGenes, , drop=FALSE]
					, introns_var=rpkms_total_introns_variances[intExGenes, , drop=FALSE]
			)
			labeledRpkmsIntEx <- list(
					exons=rpkms_Nascent_exons[intExGenes, , drop=FALSE]
					, exons_var=rpkms_Nascent_exons_variances[intExGenes, , drop=FALSE]
					, introns=rpkms_Nascent_introns[intExGenes, , drop=FALSE]
					, introns_var=rpkms_Nascent_introns_variances[intExGenes, , drop=FALSE]
				)
			## estimate the rates
			outIntEx <- .getRatesAndConcentrationsFromRpkms(totRpkms=totRpkmsIntEx
														  , labeledRpkms=labeledRpkmsIntEx
														  , tpts=tpts
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
					exons=rpkms_total_exons[onlyExGenes, , drop=FALSE]
					, exons_var=rpkms_total_exons_variances[onlyExGenes, , drop=FALSE]
				)
			labeledRpkmsOnlyEx <- list(
					exons=rpkms_Nascent_exons[onlyExGenes, , drop=FALSE]
					, exons_var=rpkms_Nascent_exons_variances[onlyExGenes, , drop=FALSE]
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
#		evaluateRateVarAtT0 <- nRepT0>1 # & !simulatedData
		evaluateRateVarAtT0 <- TRUE
		if( evaluateRateVarAtT0 ) {
	
			message('Estimating rate variances at time zero...')

			a <- rpkms_Nascent_exons[,1]/out$tL*out$labeledSF[1]
			t <- rpkms_total_exons[,1]*out$totalSF[1]
			p <- rpkms_total_introns[,1]*out$totalSF[1]

			p[setdiff(names(t),names(p))] <- NaN

			c <- a/p
			b <- a/(t - p)

			var_a <- rpkms_Nascent_exons_variances[,1]/out$tL*out$labeledSF[1]
			var_t <- rpkms_total_exons_variances[,1]*out$totalSF[1]
			var_p <- rpkms_total_introns_variances[,1]*out$totalSF[1]

			var_b <- var_a/(t-p)^2 + (a^2)/(t-p)^4*(var_t + var_p)
			var_c <- var_a/p^2 + (p^4*var_p)/a^2

			var_a[var_a<0] <- NaN
			var_t[var_t<0] <- NaN
			var_p[var_p<0] <- NaN
			var_b[var_b<0] <- NaN
			var_c[var_c<0] <- NaN

#			a <- matrix(NA, ncol=nRepT0, nrow=length(out$geneNames))
#			rownames(a) <- out$geneNames
#			p <- t <- a
#			a[rownames(rpkms_Nascent_exons),] <- rpkms_Nascent_exons[,t0ix]/out$tL*out$labeledSF[1]
#			t[rownames(rpkms_total_exons),] <- rpkms_total_exons[,t0ix]*out$totalSF[1]
#			p[rownames(rpkms_total_introns),] <- rpkms_total_introns[,t0ix]*out$totalSF[1]
#			mean_a <- apply(a, 1, mean, na.rm=TRUE)
#			mean_t <- apply(t, 1, mean, na.rm=TRUE)
#			mean_p <- apply(p, 1, mean, na.rm=TRUE)
#			mean_1_over_p <- apply(1/p, 1, mean, na.rm=TRUE)
#			mean_t_minus_p <- apply(t-p, 1, mean, na.rm=TRUE)
#			mean_1_over_t <- apply(1/t, 1, mean, na.rm=TRUE)
#			mean_1_over_t_minus_p <- apply(1/(t-p), 1, mean, na.rm=TRUE)
#	
#			matCov <- function(x, y)
#				sapply(1:nrow(x), function(i) 
#					tryCatch(cov(x[i,], y[i,], use="complete.obs"), error=function(e) NA))	
#	
#			var_a <- apply(a, 1, stats::var, na.rm=TRUE)
#			var_t <- apply(t, 1, stats::var, na.rm=TRUE)
#			var_p <- apply(p, 1, stats::var, na.rm=TRUE)
#			var_t_minus_p <- var_t + var_p - 2*matCov(t, p)
#			var_1_over_p <- apply(1/p, 1, stats::var, na.rm=TRUE)
#			var_c <- var_a_over_p <- matCov(a^2,1/p^2) + ( var_a + mean_a^2 )*( var_1_over_p + mean_1_over_p^2 ) -
#				( matCov(a, 1/p) + mean_a*mean_1_over_p )^2
#			var_1_over_t <- apply(1/t, 1, stats::var, na.rm=TRUE)
#			var_1_over_t_minus_p <- apply(1/(t-p), 1, stats::var, na.rm=TRUE)
#			var_a_over_t_minus_p <- matCov(a^2, 1/(t-p)^2) + ( var_a + mean_a^2 )*( var_1_over_t_minus_p + mean_1_over_t_minus_p^2 ) -
#				( matCov(a, 1/(t-p)) + mean_a*mean_1_over_t_minus_p )^2
#			var_a_over_t <- matCov(a^2, 1/t^2) + ( var_a + mean_a^2 )*( var_1_over_t + mean_1_over_t^2 ) -
#				( matCov(a, 1/t) + mean_a*mean_1_over_t )^2
#			var_b <- ifelse(is.na(mean_p), var_a_over_t, var_a_over_t_minus_p)
#			# mask (non-sensed) negative variances
#			var_b[var_b<0] <- NaN
#			var_c[var_c<0] <- NaN
	
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
			, nRepT0=NA
#			, nRepT0=if( evaluateRateVarAtT0 ) nRepT0 else NA
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
		object@params$NoNascent <- NoNascent
		if( is.null(rpkms_Nascent_introns) | is.null(rpkms_total_introns) )
			object@model@simple <- TRUE
	
		return(object)


	}
}