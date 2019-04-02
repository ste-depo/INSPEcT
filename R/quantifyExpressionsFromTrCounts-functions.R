#' Evaluates introns and exons RPKMs, per gene, from counts data.
#' @param allcounts A named list containing "exonsCounts" and "intronsCounts".
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates.
#' Time points must be ordered according to the sequence of files submitted for the analysis, these labels characterize
#' different files as replicates of a given condition.
#' @param exonsWidths A numeric containing the exons widths.
#' @param intronsWidths A numeric containing the intorns widths.
#' @param libsize A numeric containing the library size.
#' @param DESeq2 A logical, if TRUE the RPKMs variances are evaluated through the package DESeq2, if FALSE plgem is used.
#' @param varSamplingCondition A character reporting which experimental condition should be used to sample the variance if DESeq2 = FALSE. By default, the first element of "experimentalDesign" with replicates.
#' @return A list containing RPKMs and associated variances for exons and introns.
#' @examples
#' data('allcounts', package='INSPEcT')
#' data('featureWidths', package='INSPEcT')
#' data('libsizes', package='INSPEcT')
#' 
#' nascentCounts<-allcounts$nascent
#' matureCounts<-allcounts$mature
#' 
#' expDes<-rep(c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16),3)
#' 
#' nasExp_DESeq2<-quantifyExpressionsFromTrCounts(libsize=nascentLS
#'                                               ,exonsWidths=exWdths
#'                                               ,intronsWidths=intWdths
#'                                               ,allcounts=nascentCounts
#'                                               ,experimentalDesign=expDes)
#' 
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(libsize=totalLS
#'                                               ,exonsWidths=exWdths
#'                                               ,intronsWidths=intWdths
#'                                               ,allcounts=matureCounts
#'                                               ,experimentalDesign=expDes)
#' 
#' nasExp_plgem<-quantifyExpressionsFromTrCounts(libsize=nascentLS
#'                                              ,exonsWidths=exWdths
#'                                              ,intronsWidths=intWdths
#'                                              ,allcounts=nascentCounts
#'                                              ,DESeq2=FALSE
#'                                              ,experimentalDesign=expDes)
#' 
#' matExp_plgem<-quantifyExpressionsFromTrCounts(libsize=totalLS
#'                                              ,exonsWidths=exWdths
#'                                              ,intronsWidths=intWdths
#'                                              ,allcounts=matureCounts
#'                                              ,DESeq2=FALSE
#'                                              ,experimentalDesign=expDes)

quantifyExpressionsFromTrCounts <- function(allcounts
							  , experimentalDesign
							  , exonsWidths
							  , intronsWidths
							  , libsize = NULL
							  , DESeq2 = TRUE
							  , varSamplingCondition = NULL)
{

	############################################
	### CHECK ARGUMENTS ########################
	############################################
	# allcounts
	if( ! is.list(allcounts) )
		stop('quantifyExpressionsFromTrCounts: "allcounts" must be a list with elements "exonsCounts" and "intronsCounts"')
	if( ! all(c('exonsCounts','intronsCounts') %in% names(allcounts)) )
		stop('quantifyExpressionsFromTrCounts: "allcounts" must be a list with elements "exonsCounts" and "intronsCounts"')
	exonsCounts <- allcounts$exonsCounts
	intronsCounts <- allcounts$intronsCounts
	if( !( is.matrix(exonsCounts) & is.matrix(intronsCounts) ) )
		stop('quantifyExpressionsFromTrCounts: the elements "exonsCounts" and "intronsCounts" of "allcounts" must be matrices with the same numebr of columns.')
	if( ncol(exonsCounts) != ncol(intronsCounts) )
		stop('quantifyExpressionsFromTrCounts: the elements "exonsCounts" and "intronsCounts" of "allcounts" must be matrices with the same numebr of columns.')
	# experimentalDesign
	if(all(table(experimentalDesign)==1))
		stop("quantifyExpressionsFromTrCounts: at least one condition with replicates is required.")
	if(length(experimentalDesign)!=ncol(exonsCounts))
		stop('quantifyExpressionsFromTrCounts: each counts column must be accounted in the experimentalDesign')
	# exonsWidths
	if( !is.numeric(exonsWidths) )
		stop('quantifyExpressionsFromTrCounts: "exonsWidths" must be a numeric of length equal to the number of rows of the element "exonsCounts" of "allcounts".')
	if( nrow(exonsCounts)!=length(exonsWidths) )
		stop('quantifyExpressionsFromTrCounts: "exonsWidths" must be a numeric of length equal to the number of rows of the element "exonsCounts" of "allcounts".')
	# intronsWidths
	if( !is.numeric(intronsWidths) )
		stop('quantifyExpressionsFromTrCounts: "intronsWidths" must be a numeric of length equal to the number of rows of the element "intronsCounts" of "allcounts".')
	if( nrow(intronsCounts)!=length(intronsWidths) )
		stop('quantifyExpressionsFromTrCounts: "intronsWidths" must be a numeric of length equal to the number of rows of the element "intronsCounts" of "allcounts".')
	# libsize
	if( !is.null(libsize) ) {
		if( !is.numeric(libsize) )
			stop('quantifyExpressionsFromTrCounts: "libsize" must be either NULL or a numeric of length equal to the number of columns of the element "exonsCounts" or "intronsCounts" of "allcounts".')
		if( length(libsize) != ncol(exonsCounts) )
			stop('quantifyExpressionsFromTrCounts: "libsize" must be either NULL or a numeric of length equal to the number of columns of the element "exonsCounts" or "intronsCounts" of "allcounts".')
	} else {
		libsize <- colSums(exonsCounts) + colSums(intronsCounts)
	}
	# DESeq2
	if( !is.logical(DESeq2) )
		stop('quantifyExpressionsFromTrCounts: "DESeq2" must be a logical.')
	# varSamplingCondition
	if( !DESeq2 ) {
		if( is.null(varSamplingCondition) ) {
			varSamplingCondition <- names(which(table(experimentalDesign)>1)[1])
		} else {
			if( length(which(as.character(experimentalDesign) == varSamplingCondition)) < 2 )
				stop('quantifyExpressionsFromTrCounts: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
		}
	} 

	############################################
	### ESTIMATE EXPRESSION AND VARIANCE #######
	############################################


	############## with DESeq2 #################
	if(DESeq2)
	{

		message('Estimation of expressions and variances using DESeq2...')

		if( is.numeric(experimentalDesign) ) 
			experimentalDesign= signif(experimentalDesign,9)
		sampleTptsNames <- factor(experimentalDesign)
		colData <- data.frame(tpts=sampleTptsNames)
		countsTemp <- list(exonsCounts,intronsCounts)

		singleConditionAnalysis <- length(levels(sampleTptsNames)) == 1

		if( singleConditionAnalysis ) {

			ddsList <- lapply(countsTemp,function(countData){
				ori <- DESeqDataSetFromMatrix(countData = countData
											,colData = colData
											,design = ~ 1)
				suppressWarnings(suppressMessages(ori <- DESeq(ori)))
				return(ori)})

		} else {

			ddsList <- lapply(countsTemp,function(countData){
				ori <- DESeqDataSetFromMatrix(countData = countData
											,colData = colData
											,design = ~ tpts)
				suppressMessages(ori <- DESeq(ori))
				return(ori)})

		}

		muExons<-assays(ddsList[[1]])[['mu']]
		alphaExons<-dispersions(ddsList[[1]])
		names(alphaExons) <- rownames(muExons)

		muIntrons<-assays(ddsList[[2]])[['mu']]
		alphaIntrons<-dispersions(ddsList[[2]])
		names(alphaIntrons) <- rownames(muIntrons)

		varCountsExons<-t(sapply(rownames(muExons),function(gene)
			{muExons[gene,]+alphaExons[gene]*muExons[gene,]^2}))

		varCountsIntrons<-t(sapply(rownames(muIntrons),function(gene)
			{muIntrons[gene,]+alphaIntrons[gene]*muIntrons[gene,]^2}))

		#from counts to RPKMs

		expressionExonsDESeq<-counts2expressions(exonsCounts,exonsWidths,libsize)
		expressionIntronsDESeq<-counts2expressions(intronsCounts,intronsWidths,libsize)

  		varianceExpressionsExons<-countVar2expressions(varCountsExons,exonsWidths,libsize)
  		varianceExpressionsIntrons<-countVar2expressions(varCountsIntrons,intronsWidths,libsize)

  		if( singleConditionAnalysis ) {

			expressionsExons <- as.matrix(rowMeans(expressionExonsDESeq))
			expressionsIntrons <- as.matrix(rowMeans(expressionIntronsDESeq))

			varianceExpressionsExons <- as.matrix(rowMeans(varianceExpressionsExons))
			varianceExpressionsIntrons <- as.matrix(rowMeans(varianceExpressionsIntrons))

		} else {

			expressionsExons <- t(apply(expressionExonsDESeq,1,function(x)tapply(x, experimentalDesign, mean)))
			expressionsIntrons <- t(apply(expressionIntronsDESeq,1,function(x)tapply(x, experimentalDesign, mean)))

			varianceExpressionsExons <- t(apply(varianceExpressionsExons,1,function(x)tapply(x, experimentalDesign, mean)))
			varianceExpressionsIntrons <- t(apply(varianceExpressionsIntrons,1,function(x)tapply(x, experimentalDesign, mean)))

		}

		rownames(expressionsExons) <- rownames(varianceExpressionsExons) <- rownames(muExons)
		rownames(expressionsIntrons) <- rownames(varianceExpressionsIntrons) <- rownames(muIntrons)

		colnames(expressionsExons) <- colnames(varianceExpressionsExons) <- sort(unique(experimentalDesign))
		colnames(expressionsIntrons) <- colnames(varianceExpressionsIntrons) <- sort(unique(experimentalDesign))

		return(list(exonsExpressions = expressionsExons
				  , intronsExpressions = expressionsIntrons
				  , exonsVariance = varianceExpressionsExons
				  , intronsVariance = varianceExpressionsIntrons
				  ))

	} else {

		message('Estimation of expressions and variances using plgem...')

		expressionExons<-counts2expressions(exonsCounts,exonsWidths,libsize)
		expressionIntrons<-counts2expressions(intronsCounts,intronsWidths,libsize)

		############## with PLGEM ##################
		return(quantifyExpressionsFromTrAbundance(
								trAbundaces = list(
									exonsAbundances = expressionExons
									, intronsAbundances = expressionIntrons
									)
								, experimentalDesign = experimentalDesign
								, varSamplingCondition = varSamplingCondition))

	}
}
