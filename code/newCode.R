library(devtools)
load_all("Dropbox/INSPEcT_Cluster")

library(DESeq2)
library(plgem)
library(GenomicAlignments)
library(GenomicFeatures)  

source("Dropbox/INSPEcT_Dev/R/AllInternalFunctions.R")
source("Dropbox/INSPEcT_Dev/R/makeRPKMsFromCounts-functions.R")
source("Dropbox/INSPEcT_Dev/R/makeRPKMsFromBams-functions.R")
source("Dropbox/INSPEcT_Dev/R/quantifyExpressionFromTrAbundance.R")
source("Dropbox/INSPEcT_Dev/R/newINSPEcT-functions.R")

#Test on BAMs
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

paths <- system.file('extdata/', c('bamRep1.bam','bamRep2.bam','bamRep3.bam','bamRep4.bam'), package='INSPEcT')

experimentalDesign <- c(0,0,1,1)

makeExpOutFromBAMs <- makeExpressionsFromBAMs(txdb = txdb
											, BAMfiles = paths
											, DESeq2 = TRUE
											, experimentalDesign = experimentalDesign
											, varSamplingCondition = NULL) 

#Test on Counts
data('allcounts4su', package='INSPEcT')

exonsDB <- reduce(exonsBy(txdb ,'gene'))
exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
intronsDB <- intronsDB[elementNROWS(intronsDB)>0]

experimentalDesign <- rep(c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16),3)

testGenes <- rownames(allcounts$exonCounts)

allcounts <- allcounts4su$total

libsize <- colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
exonsWidths <- sapply(width(exonsDB),sum)
intronsWidths <- sapply(width(intronsDB),sum)

makeExpOutFromCounts_DESeq2_total <- quantifyExpressionsFromTrCounts(libsize = libsize
															 , exonsWidths = exonsWidths[testGenes]
															 , intronsWidths = intronsWidths[testGenes]
															 , allcounts = list(exonsCounts = allcounts$exonCounts, intronsCounts = allcounts$intronCounts)
															 , DESeq2 = TRUE
															 , experimentalDesign = experimentalDesign
															 , varSamplingCondition = NULL)

makeExpOutFromCounts_plgem_total <- quantifyExpressionsFromTrCounts(libsize = libsize
																   , exonsWidths = exonsWidths[testGenes]
																   , intronsWidths = intronsWidths[testGenes]
																   , allcounts = list(exonsCounts = allcounts$exonCounts, intronsCounts = allcounts$intronCounts)
																   , DESeq2 = FALSE
																   , experimentalDesign = experimentalDesign
																   , varSamplingCondition = as.character(experimentalDesign[[1]]))

allcounts <- allcounts4su$foursu

libsize <- colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
exonsWidths <- sapply(width(exonsDB),sum)
intronsWidths <- sapply(width(intronsDB),sum)

makeExpOutFromCounts_DESeq2_nascent <- quantifyExpressionsFromTrCounts(libsize = libsize
															 , exonsWidths = exonsWidths[testGenes]
															 , intronsWidths = intronsWidths[testGenes]
															 , allcounts = list(exonsCounts = allcounts$exonCounts, intronsCounts = allcounts$intronCounts)
															 , DESeq2 = TRUE
															 , experimentalDesign = experimentalDesign
															 , varSamplingCondition = NULL)

makeExpOutFromCounts_plgem_nascent <- quantifyExpressionsFromTrCounts(libsize = libsize
																	, exonsWidths = exonsWidths[testGenes]
																	, intronsWidths = intronsWidths[testGenes]
																	, allcounts = list(exonsCounts = allcounts$exonCounts, intronsCounts = allcounts$intronCounts)
																	, DESeq2 = FALSE
																	, experimentalDesign = experimentalDesign
																	, varSamplingCondition = as.character(experimentalDesign[[1]]))

#New INSPEcT
inspectObject_DESeq2_withNascent <- newINSPEcT(tpts = c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
						  					 , labeling_time = 1/6
						  					 , nascentExpressions = makeExpOutFromCounts_DESeq2_nascent
						  					 , totalExpressions = makeExpOutFromCounts_DESeq2_total)

inspectObject_DESeq2 <- newINSPEcT(tpts = c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
						  		 , labeling_time = NULL
						  		 , nascentExpressions = NULL
						  		 , totalExpressions = makeExpOutFromCounts_DESeq2_total)

inspectObject_plgem_withNascent <- newINSPEcT(tpts = c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
											, labeling_time = 1/6
											, nascentExpressions = makeExpOutFromCounts_plgem_nascent
											, totalExpressions = makeExpOutFromCounts_plgem_total)

inspectObject_plgem <- newINSPEcT(tpts = c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
								, labeling_time = NULL
								, nascentExpressions = NULL
								, totalExpressions = makeExpOutFromCounts_plgem_total)

pdf("plgemVSdeseq2.pdf")
par(mfrow=c(3,2))

smoothScatter(log10(ratesFirstGuess(inspectObject_DESeq2_withNascent,"synthesis")),log10(ratesFirstGuess(inspectObject_plgem_withNascent,"synthesis")),xlab="",ylab="",main="k1 - nascent");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuess(inspectObject_DESeq2,"synthesis")),log10(ratesFirstGuess(inspectObject_plgem,"synthesis")),xlab="",ylab="",main="k1");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuess(inspectObject_DESeq2_withNascent,"processing")),log10(ratesFirstGuess(inspectObject_plgem_withNascent,"processing")),xlab="",ylab="",main="k2 - nascent");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuess(inspectObject_DESeq2,"processing")),log10(ratesFirstGuess(inspectObject_plgem,"processing")),xlab="",ylab="",main="k2");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuess(inspectObject_DESeq2_withNascent,"degradation")),log10(ratesFirstGuess(inspectObject_plgem_withNascent,"degradation")),xlab="",ylab="",main="k3 - nascent");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuess(inspectObject_DESeq2,"degradation")),log10(ratesFirstGuess(inspectObject_plgem,"degradation")),xlab="",ylab="",main="k3");abline(0,1,col="red")

par(mfrow=c(2,2))

smoothScatter(log10(ratesFirstGuessVar(inspectObject_DESeq2_withNascent,"total")),log10(ratesFirstGuessVar(inspectObject_plgem_withNascent,"total")),xlab="DESeq2",ylab="plgem",main="total - nascent");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuessVar(inspectObject_DESeq2,"total")),log10(ratesFirstGuessVar(inspectObject_plgem,"total")),xlab="DESeq2",ylab="plgem",main="total");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuessVar(inspectObject_DESeq2_withNascent,"preMRNA")),log10(ratesFirstGuessVar(inspectObject_plgem_withNascent,"preMRNA")),xlab="DESeq2",ylab="plgem",main="pre - nascent");abline(0,1,col="red")
smoothScatter(log10(ratesFirstGuessVar(inspectObject_DESeq2,"preMRNA")),log10(ratesFirstGuessVar(inspectObject_plgem,"preMRNA")),xlab="DESeq2",ylab="plgem",main="pre");abline(0,1,col="red")

dev.off()


