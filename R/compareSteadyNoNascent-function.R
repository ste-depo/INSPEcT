#' Compare mature RNA steady state data
#' @description
#' This function compare exons and introns expression level matrices, from two up to an arbitrary number of samples, in order
#' to identify genes which are oddly regluated, compared to an expected standard behaviour, from the post transcriptional point of view.
#' @param expressionMatrices A list of two expression matrices, for exons and introns
#' @param expressionThreshold Entries of exons and introns expression matrices are not take into during the analysis except for the evaluation of genes standard behaviour
#' @param log2FCThreshold A parameter which sets how many log2 fold changes of distance from the median behaviour are imputable to noise.
#' @examples
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
#' expDes<-rep(c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16),3)
#'   
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(libsize=totalLS
#'                                             ,exonsWidths=exWdths
#'                                             ,intronsWidths=intWdths
#'                                             ,allcounts=matureCounts
#'                                             ,DESeq2=TRUE
#'                                             ,experimentalDesign=expDes)
#' 
#' regGenes<-compareSteadyNoNascent(expressionMatrices=matExp_DESeq2
#' 								   ,expressionThreshold=0.25
#'								   ,log2FCThreshold=2.)
#' head(regGenes)
#' table(regGenes)


compareSteadyNoNascent <- function(expressionMatrices,expressionThreshold=0.25,log2FCThreshold=2.)
{

	#	var_ex <- err$exons
	#	var_int <- err$introns
	#
	#	var_m <- var_int + var_ex
	#	var_p <- var_int

	eiGenes <- intersect(rownames(expressionMatrices$exonsExpressions),rownames(expressionMatrices$intronsExpressions))
	print(paste0("compareSteadyNoNascent: ",length(eiGenes)," genes with introns and exons under analysis."))
	# Mature, premature and total rpkms
	mature <- expressionMatrices$exonsExpressions[eiGenes,] - expressionMatrices$intronsExpressions[eiGenes,]
	premature <- expressionMatrices$intronsExpressions[eiGenes,]
	total <- expressionMatrices$exonsExpressions[eiGenes,]

	#	var_log2p <- var_p/(median(premature,na.rm=TRUE)*log(2))^2
	#	var_log2m <- var_m/(median(mature,na.rm=TRUE)*log(2))^2
	#	sigma_log2p <- sqrt(var_log2p)
	#	sigma_log2m <- sqrt(var_log2m)
	#	err <- (sigma_log2m*cos(atan(sigma_log2m/sigma_log2p)))*sigmas

	#	print(paste0("compareSteadyNoNascent: acceptance parameter equal ",round(err,2)))

	prematureMedian <- apply(premature,1,function(r)median(r,na.rm=T))
	matureMedian <- apply(mature,1,function(r)median(r,na.rm=T))
	
	suppressWarnings(standardCurveFit <- standardCurveFitFunction(p=prematureMedian,m=matureMedian,err=log2FCThreshold))

	premature[premature<=expressionThreshold] <- NA
	mature[mature<=expressionThreshold] <- NA

	suppressWarnings(classificationTmp <- classificationFunction(p=premature,m=mature,alpha=standardCurveFit,err=log2FCThreshold))

	#	classificationTmp[!is.finite(classificationTmp)] <- FALSE

	return(classificationTmp)
}