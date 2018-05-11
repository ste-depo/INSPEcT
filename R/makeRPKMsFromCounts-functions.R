#' Calculate RPKM and count values on introns and exons from bam/sam files
#' @description Given a TranscriptDb object and a matrix of counts for total and eventually RNA 
#' experiments, "makeRPKMsFromCounts" function calculates RPKM on exonic and intronic 
#' features per each gene. Reads that fall where intronic and exonic features overlaps are 
#' univoquely assigned to exons.
#' @param txdb A TranscriptDB object
#' @param allcounts A list object containing intronic and exonic counts from total and eventually 4sU experiments with the associated statistics redarding the assigned reads.
#' @param by A character, either "gene" or "tx", indicating if rpkms and counts should be summarized at the levels of genes or transcripts. "gene" by default
#' @param DESeq2 A logical, if TRUE the RPKMs from exons and introns and associated variances are evaluated through the package DESeq2
#' @param temporalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the columns of the count matrices submitted for the analysis; these labels define conditions and replicates.
#' @return A list containing rpkms, counts and the annotation extracted from TxDB for exons and introns, if DESeq2 = TRUE the output also contains a set of data
#' needed to estimate rpkms variances.
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' data('allcounts4su', package='INSPEcT')
#' 
#' tpts <- c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
#' temporalDesign <- rep(tpts,3)
#' 
#' #4sU analysis with DESeq2
#' makeRPKMsOut_4sU <- makeRPKMsFromCounts(txdb=txdb,allcounts=allcounts4su,temporalDesign=temporalDesign,DESeq2=TRUE)
#' 
#' rpkms_4sU <- makeRPKMsOut_4sU$rpkms
#' counts_4sU <- makeRPKMsOut_4sU$counts
#' annotations_4sU <- makeRPKMsOut_4sU$annotation
#' dispersion_parameters_DESeq2_4sU <- makeRPKMsOut_4sU$dispersion_parameters_DESeq2
#' 
#' #4sU analysis without DESeq2
#' makeRPKMsOut_4sU <- makeRPKMsFromCounts(txdb=txdb,allcounts=allcounts4su,temporalDesign=temporalDesign,DESeq2=FALSE)
#' 
#' rpkms_4sU <- makeRPKMsOut_4sU$rpkms
#' counts_4sU <- makeRPKMsOut_4sU$counts
#' annotations_4sU <- makeRPKMsOut_4sU$annotation
#' 
#' #No4sU analysis with DESeq2
#' allcounts <- allcounts4su$total
#' 
#' makeRPKMsOut_No4sU <- makeRPKMsFromCounts(txdb=txdb,allcounts=allcounts,temporalDesign=temporalDesign,DESeq2=TRUE)
#' 
#' rpkms_No4sU <- makeRPKMsOut_No4sU$rpkms
#' counts_No4sU <- makeRPKMsOut_No4sU$counts
#' annotations_No4sU <- makeRPKMsOut_No4sU$annotation
#' dispersion_parameters_DESeq2_No4sU <- makeRPKMsOut_No4sU$dispersion_parameters_DESeq2

makeRPKMsFromCounts <- function(txdb
							  , allcounts = NULL
							  , by = c('gene','tx')
							  , DESeq2 = TRUE
							  , temporalDesign = NULL) 
{

	if(all(table(temporalDesign)==1))stop("makeRPKMs: at least one replicate is required.")
	counts2rpkms <- function(counts, widths, libsize) counts*(10^9/(widths[rownames(counts)]%o%libsize))

	if(is.null(allcounts$foursu)){
		print("makeRPKMs: No4sU mode.")
		No4sU <- TRUE
	}else{
		No4sU <- FALSE
	}

	message('Generating annotation from txdb...')

	by <- by[1]
	if( by=="gene" ) {

		exonsDB <- reduce(exonsBy(txdb ,'gene'))
		exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
		intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
		intronsDB <- intronsDB[elementNROWS(intronsDB)>0]

	} else if( by=="tx" ){

		exonsDB <- exonsBy(txdb ,'tx', use.names=TRUE)
		intronsDB <- intronsByTranscript(txdb, use.names=TRUE)
		intronsDB <- intronsDB[elementNROWS(intronsDB)>0]

	} else {

		stop("makeRPKMs: argument 'by' not recognized.")

	}

	if(No4sU)
	{
		allrpkms <- list(
			total_exons=counts2rpkms(
				allcounts$exonCounts
				, sapply(width(exonsDB),sum)
				, colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
				)
			,total_introns=counts2rpkms(
				allcounts$intronCounts
				, sapply(width(intronsDB),sum)
				, colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
				)
		)
	
		if(length(table(sapply(allrpkms,ncol)))!=1) stop('makeRPKMs: all rpkm matrices must have the same number of columns')

		outPutTemp <- list(rpkms=allrpkms
						 , counts=allcounts
						 , annotation=list(exon=exonsDB
						 				 , intron=intronsDB))

		if(DESeq2)
		{
			if(is.null(temporalDesign)) stop('makeRPKMs: the DESeq2 analysis requires a temporalDesign')
			if(length(temporalDesign)!=ncol(allrpkms$total_exons)) stop('makeRPKMs: each bam file must be accounted in the temporalDesign')
			sampleTptsNames <- factor(signif(temporalDesign,2))
			colData <- data.frame(tpts=sampleTptsNames)
  
			countsTemp <- list(outPutTemp$counts$exonCounts,outPutTemp$counts$intronCounts)

			ddsList <- lapply(countsTemp,function(countData){ori <- DESeqDataSetFromMatrix(countData = countData
																	,colData = colData
																	,design = ~ tpts)
						return(DESeq(ori))})

			countsExonsDESeq<-counts(ddsList[[1]],normalized=FALSE)
			countsIntronsDESeq<-counts(ddsList[[2]],normalized=FALSE)
	
			muExons<-assays(ddsList[[1]])[['mu']]
			alphaExons<-dispersions(ddsList[[1]])
	
			names(alphaExons) <- rownames(muExons)
	
			muIntrons<-assays(ddsList[[2]])[['mu']]
			alphaIntrons<-dispersions(ddsList[[2]])
	
			names(alphaIntrons) <- rownames(muIntrons)
	
			rpkmExonsDESeq<-counts2rpkms(countsExonsDESeq
										   , sapply(width(exonsDB),sum)
										   , colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE]))
			rpkmIntronsDESeq<-counts2rpkms(countsIntronsDESeq
											 , sapply(width(intronsDB),sum)
											 , colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE]))
	
			exonsWidths <- sapply(width(exonsDB),sum)
			intronsWidths <- sapply(width(intronsDB),sum)
			libsizesTotal <- colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
	
			outPutTemp <- list(rpkms = list(total_exons = rpkmExonsDESeq,total_introns = rpkmIntronsDESeq)
							   ,counts = list(total = list(exonCounts = countsExonsDESeq,intronCounts = countsIntronsDESeq))
							   ,annotation = outPutTemp$annotation
							   ,dispersion_parameters_DESeq2 = list(total = list(muExons = muExons
																			   ,muIntrons = muIntrons
																			   ,alphaExons = alphaExons
																			   ,alphaIntrons = alphaIntrons
																			   ,exonsWidths = exonsWidths[names(alphaExons)]
																			   ,intronsWidths = intronsWidths[names(alphaIntrons)]
																			   ,libsizes = libsizesTotal))
							   )
		}
	
		return(outPutTemp)

	}else{
	
		allrpkms <- list(
			foursu_exons=counts2rpkms(
				allcounts$foursu$exonCounts
				, sapply(width(exonsDB),sum)
				, colSums(allcounts$foursu$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
				)
			,foursu_introns=counts2rpkms(
				allcounts$foursu$intronCounts
				, sapply(width(intronsDB),sum)
				, colSums(allcounts$foursu$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
				)
			,total_exons=counts2rpkms(
				allcounts$total$exonCounts
				, sapply(width(exonsDB),sum)
				, colSums(allcounts$total$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
				)
			,total_introns=counts2rpkms(
				allcounts$total$intronCounts
				, sapply(width(intronsDB),sum)
				, colSums(allcounts$total$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
				)
		)

		if(length(table(sapply(allrpkms,ncol)))!=1) stop('makeRPKMs: all rpkm matrices must have the same number of columns')

		outPutTemp <- list(rpkms=allrpkms, counts=allcounts, annotation=list(exon=exonsDB, intron=intronsDB))

		if(DESeq2)
		{

			if(is.null(temporalDesign)) stop('makeRPKMs: the DESeq2 analysis requires a temporalDesign')
			if(length(temporalDesign)!=ncol(allrpkms$total_exons)) stop('makeRPKMs: each bam file must be accounted in the temporalDesign')

			sampleTptsNames <- factor(signif(temporalDesign,2))
			colData <- data.frame(tpts=sampleTptsNames)

			countsTemp <- list(outPutTemp$counts$foursu$exonCounts
							  ,outPutTemp$counts$foursu$intronCounts
							  ,outPutTemp$counts$total$exonCounts
							  ,outPutTemp$counts$total$intronCounts
							  )

			ddsList <- lapply(countsTemp,function(countData){ori <- DESeqDataSetFromMatrix(countData = countData
																	,colData = colData
																	,design = ~ tpts)
			return(DESeq(ori))})

			countsExons4sUDESeq<-counts(ddsList[[1]],normalized=FALSE)
			countsIntrons4sUDESeq<-counts(ddsList[[2]],normalized=FALSE)
			countsExonsNo4sUDESeq<-counts(ddsList[[3]],normalized=FALSE)
			countsIntronsNo4sUDESeq<-counts(ddsList[[4]],normalized=FALSE)

			muExons4sU<-assays(ddsList[[1]])[['mu']]
			alphaExons4sU<-dispersions(ddsList[[1]])

			muIntrons4sU<-assays(ddsList[[2]])[['mu']]
			alphaIntrons4sU<-dispersions(ddsList[[2]])

			muExonsNo4sU<-assays(ddsList[[3]])[['mu']]
			alphaExonsNo4sU<-dispersions(ddsList[[3]])

			muIntronsNo4sU<-assays(ddsList[[4]])[['mu']]
			alphaIntronsNo4sU<-dispersions(ddsList[[4]])

			names(alphaExons4sU) <- rownames(muExons4sU)
			names(alphaIntrons4sU) <- rownames(muIntrons4sU)
			names(alphaExonsNo4sU) <- rownames(muExonsNo4sU)
			names(alphaIntronsNo4sU) <- rownames(muIntronsNo4sU)

			rpkmExons4sUDESeq<-counts2rpkms(countsExons4sUDESeq
									   , sapply(width(exonsDB),sum)
									   , colSums(allcounts$foursu$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE]))
			rpkmIntrons4sUDESeq<-counts2rpkms(countsIntrons4sUDESeq
										 , sapply(width(intronsDB),sum)
										 , colSums(allcounts$foursu$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE]))

			rpkmExonsNo4sUDESeq<-counts2rpkms(countsExonsNo4sUDESeq
									   , sapply(width(exonsDB),sum)
									   , colSums(allcounts$total$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE]))
			rpkmIntronsNo4sUDESeq<-counts2rpkms(countsIntronsNo4sUDESeq
										 , sapply(width(intronsDB),sum)
										 , colSums(allcounts$total$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE]))

			exonsWidths <- sapply(width(exonsDB),sum)
			intronsWidths <- sapply(width(intronsDB),sum)

			libsizesFoursu <- colSums(allcounts$foursu$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])
			libsizesTotal <- colSums(allcounts$total$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])

			outPutTemp <- list(rpkms = list(foursu_exons = rpkmExons4sUDESeq
										   ,foursu_introns = rpkmIntrons4sUDESeq
										   ,total_exons = rpkmExonsNo4sUDESeq
										   ,total_introns = rpkmIntronsNo4sUDESeq)
							  ,counts = list(foursu = list(exonCounts = countsExons4sUDESeq
														  ,intronCounts = countsIntrons4sUDESeq)
											,total = list(exonCounts = countsExonsNo4sUDESeq
														 ,intronCounts = countsIntronsNo4sUDESeq))
							  ,annotation = outPutTemp$annotation
							  ,dispersion_parameters_DESeq2 = list(foursu = list(muExons = muExons4sU
																				,muIntrons = muIntrons4sU
																				,alphaExons = alphaExons4sU
																				,alphaIntrons = alphaIntrons4sU
																				,exonsWidths = exonsWidths[names(alphaExons4sU)]
																				,intronsWidths = intronsWidths[names(alphaIntrons4sU)]
																				,libsizes = libsizesFoursu)
																  ,total = list(muExons = muExonsNo4sU
																			   ,muIntrons = muIntronsNo4sU
																			   ,alphaExons = alphaExonsNo4sU
																			   ,alphaIntrons = alphaIntronsNo4sU
																			   ,exonsWidths = exonsWidths[names(alphaExonsNo4sU)]
																			   ,intronsWidths = intronsWidths[names(alphaIntronsNo4sU)]
																			   ,libsizes = libsizesTotal)))
		
		}

		return(outPutTemp)
	}
}
