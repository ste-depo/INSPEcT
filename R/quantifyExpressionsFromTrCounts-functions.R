#' Calculate RPKM and count values on introns and exons from bam/sam files
#' @description Given a TranscriptDb object and a matrix of counts for total and eventually RNA 
#' experiments, "quantifyExpressionsFromTrCounts" function calculates RPKM on exonic and intronic 
#' features per each gene. Reads that fall where intronic and exonic features overlaps are 
#' univoquely assigned to exons.
#' @param txdb A TranscriptDB object
#' @param allcounts A list object containing intronic and exonic counts from total and eventually Nascent experiments with the associated statistics redarding the assigned reads.
#' @param by A character, either "gene" or "tx", indicating if rpkms and counts should be summarized at the levels of genes or transcripts. "gene" by default
#' @param DESeq2 A logical, if TRUE the RPKMs from exons and introns and associated variances are evaluated through the package DESeq2
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the columns of the count matrices submitted for the analysis; these labels define conditions and replicates.
#' @return A list containing rpkms, counts and the annotation extracted from TxDB for exons and introns, if DESeq2 = TRUE the output also contains a set of data
#' needed to estimate rpkms variances.
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' data('allcountsNascent', package='INSPEcT')
#' 
#' tpts <- c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
#' experimentalDesign <- rep(tpts,3)
#' 
#' #Nascent analysis with DESeq2
#' makeRPKMsOut_Nascent <- quantifyExpressionsFromTrCounts(txdb=txdb,allcounts=allcountsNascent,experimentalDesign=experimentalDesign,DESeq2=TRUE)
#' 
#' rpkms_Nascent <- makeRPKMsOut_Nascent$rpkms
#' counts_Nascent <- makeRPKMsOut_Nascent$counts
#' annotations_Nascent <- makeRPKMsOut_Nascent$annotation
#' dispersion_parameters_DESeq2_Nascent <- makeRPKMsOut_Nascent$dispersion_parameters_DESeq2
#' 
#' #Nascent analysis without DESeq2
#' makeRPKMsOut_Nascent <- quantifyExpressionsFromTrCounts(txdb=txdb,allcounts=allcountsNascent,experimentalDesign=experimentalDesign,DESeq2=FALSE)
#' 
#' rpkms_Nascent <- makeRPKMsOut_Nascent$rpkms
#' counts_Nascent <- makeRPKMsOut_Nascent$counts
#' annotations_Nascent <- makeRPKMsOut_Nascent$annotation
#' 
#' #NoNascent analysis with DESeq2
#' allcounts <- allcountsNascent$total
#' 
#' makeRPKMsOut_NoNascent <- quantifyExpressionsFromTrCounts(txdb=txdb,allcounts=allcounts,experimentalDesign=experimentalDesign,DESeq2=TRUE)
#' 
#' rpkms_NoNascent <- makeRPKMsOut_NoNascent$rpkms
#' counts_NoNascent <- makeRPKMsOut_NoNascent$counts
#' annotations_NoNascent <- makeRPKMsOut_NoNascent$annotation
#' dispersion_parameters_DESeq2_NoNascent <- makeRPKMsOut_NoNascent$dispersion_parameters_DESeq2

quantifyExpressionsFromTrCounts <- function(libsize
							  , exonsWidths
							  , intronsWidths
							  , allcounts
							  , by = c('gene','tx')
							  , DESeq2 = TRUE
							  , experimentalDesign
							  , varSamplingCondition = NULL
							  , plgemFits = NULL
							  , returnPlgemFits = FALSE)
{
	exonsCounts <- allcounts$exonsCounts
	intronsCounts <- allcounts$intronsCounts

	############################################
	### CHECK ARGUMENTS ########################
	############################################
	if( !is.logical(DESeq2) & !any(as.character(experimentalDesign)==varSamplingCondition) & is.null(plgemFits))
		stop('makeExpressions: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
	if(all(table(experimentalDesign)==1))
		stop("makeExpressions: at least one replicate is required.")
	if(length(experimentalDesign)!=ncol(exonsCounts))
		stop('makeExpressions: each counts column must be accounted in the experimentalDesign')
	if(all(table(experimentalDesign)==1))
		stop("makeExpressions: at least one replicate is required.")
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
											   , varSamplingCondition = varSamplingCondition
											   , plgemFits = plgemFits
											   , returnPlgemFits = returnPlgemFits))

	}
}