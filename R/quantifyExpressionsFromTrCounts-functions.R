#' Evaluates introns and exons RPKMs, per gene, from counts data.
#' @param libsize A numeric reporting the number of assigned reads.
#' @param exonsWidths A numeric containing the exons widths.
#' @param intronsWidths A numeric containing the intorns widths.
#' @param allcounts A list object containing introns and exons counts.
#' @param by A character, either "gene" or "tx", indicating if rpkms and counts should be summarized at the levels of genes or transcripts. "gene" by default
#' @param DESeq2 A logical, if TRUE the RPKMs variances are evaluated through the package DESeq2, if FALSE plgem is used.
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the columns of the count matrices submitted for the analysis; these labels define conditions and replicates.
#' @param varSamplingCondition A character reporting which experimental condition should be used to sample the fatiance if DESeq2 = FALSE.
#' @return A list containing RPKMs and associated variances for exons and introns.
#' @examples
#' data('allcounts', package='INSPEcT')
#' 
#' nascentCounts<-allcounts$nascent
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
#' nascentLS<-colSums(nascentCounts$
#'   stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
#' 
#' totalLS<-colSums(matureCounts$
#'   stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
#' 
#' expDes<-rep(c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16),3)
#' 
#' nasExp_DESeq2<-quantifyExpressionsFromTrCounts(libsize=nascentLS
#'                                               ,exonsWidths=exWdths
#'                                               ,intronsWidths=intWdths
#'                                               ,allcounts=nascentCounts
#'                                               ,DESeq2=TRUE
#'                                               ,experimentalDesign=expDes)
#' 
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(libsize=totalLS
#'                                               ,exonsWidths=exWdths
#'                                               ,intronsWidths=intWdths
#'                                               ,allcounts=matureCounts
#'                                               ,DESeq2=TRUE
#'                                               ,experimentalDesign=expDes)
#' 
#' vsc<-as.character(expDes[[1]])
#' 
#' nasExp_plgem<-quantifyExpressionsFromTrCounts(libsize=nascentLS
#'                                              ,exonsWidths=exWdths
#'                                              ,intronsWidths=intWdths
#'                                              ,allcounts=nascentCounts
#'                                              ,DESeq2=FALSE
#'                                              ,experimentalDesign=expDes
#'                                              ,varSamplingCondition=vsc)
#' 
#' matExp_plgem<-quantifyExpressionsFromTrCounts(libsize=totalLS
#'                                              ,exonsWidths=exWdths
#'                                              ,intronsWidths=intWdths
#'                                              ,allcounts=matureCounts
#'                                              ,DESeq2=FALSE
#'                                              ,experimentalDesign=expDes
#'                                              ,varSamplingCondition=vsc)

quantifyExpressionsFromTrCounts <- function(libsize
							  , exonsWidths
							  , intronsWidths
							  , allcounts
							  , by = c('gene','tx')
							  , DESeq2 = TRUE
							  , experimentalDesign
							  , varSamplingCondition = NULL)
{
	exonsCounts <- allcounts$exonsCounts
	intronsCounts <- allcounts$intronsCounts

	############################################
	### CHECK ARGUMENTS ########################
	############################################
	if( !is.logical(DESeq2) & !any(as.character(experimentalDesign)==varSamplingCondition))
		stop('makeExpressions: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
	if(all(table(experimentalDesign)==1))
		stop("makeExpressions: at least one replicate is required.")
	if(length(experimentalDesign)!=ncol(exonsCounts))
		stop('makeExpressions: each counts column must be accounted in the experimentalDesign')
	if(ncol(exonsCounts)!=ncol(intronsCounts)
	 | ncol(exonsCounts)!=length(libsize)
	 | nrow(exonsCounts)!=length(exonsWidths)
	 | nrow(intronsCounts)!=length(intronsWidths))
		stop('makeExpressions: dimensionality issue in input data.')


	############################################
	### ESTIMATE EXPRESSION AND VARIANCE #######
	############################################

	############## with DESeq2 #################
	if(DESeq2)
	{
		sampleTptsNames <- factor(signif(experimentalDesign,2))
		colData <- data.frame(tpts=sampleTptsNames)
		countsTemp <- list(exonsCounts,intronsCounts)

		ddsList <- lapply(countsTemp,function(countData){ori <- DESeqDataSetFromMatrix(countData = countData
																,colData = colData
																,design = ~ tpts)
					return(DESeq(ori))})

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

		expressionsExons <- t(apply(expressionExonsDESeq,1,function(x)tapply(x, experimentalDesign, mean)))
		expressionsIntrons <- t(apply(expressionIntronsDESeq,1,function(x)tapply(x, experimentalDesign, mean)))

		varianceExpressionsExons <- t(apply(varianceExpressionsExons,1,function(x)tapply(x, experimentalDesign, mean)))
		varianceExpressionsIntrons <- t(apply(varianceExpressionsIntrons,1,function(x)tapply(x, experimentalDesign, mean)))

		rownames(expressionsExons) <- rownames(varianceExpressionsExons) <- rownames(muExons)
		rownames(expressionsIntrons) <- rownames(varianceExpressionsIntrons) <- rownames(muIntrons)

		colnames(expressionsExons) <- colnames(varianceExpressionsExons) <- sort(unique(experimentalDesign))
		colnames(expressionsIntrons) <- colnames(varianceExpressionsIntrons) <- sort(unique(experimentalDesign))

		return(list(exonsExpressions = expressionsExons
				  , intronsExpressions = expressionsIntrons
				  , exonsVariance = varianceExpressionsExons
				  , intronsVariance = varianceExpressionsIntrons))

	} else {

		expressionExons<-counts2expressions(exonsCounts,exonsWidths,libsize)
		expressionIntrons<-counts2expressions(intronsCounts,intronsWidths,libsize)

		############## with PLGEM ##################
		return(quantifyExpressionsFromTrAbundance(exonsAbundances = expressionExons
										       , intronsAbundances = expressionIntrons
											   , experimentalDesign = experimentalDesign
											   , varSamplingCondition = varSamplingCondition))

	}
}