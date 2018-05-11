#' Calculate RPKM and count values on introns and exons from bam/sam files
#' @description Given a TranscriptDb object and a list of bam/sam files for total and eventually RNA 
#' experiments, "makeRPKMsFromBams" function calculates read counts and RPKM on exonic and intronic 
#' features per each gene. Reads that fall where intronic and exonic features overlaps are 
#' univoquely assigned to exons.
#' @param txdb A TranscriptDB object
#' @param paths_foursu A vector of paths of 4sU-seq sam files
#' @param paths_total A vector of paths of RNA-seq sam files
#' @param by A character, either "gene" or "tx", indicating if rpkms and counts should be summarized at the levels of genes or transcripts. "gene" by default
#' @param countMultiMappingReads A logical, if multimapping reads should be counted, FALSE by default. Multimap reads are 
#' identified using the tag "NH" in the bam/sam file.
#' @param allowMultiOverlap A logical, indicating if a read is allowed to be assigned to more than one feature, FALSE by default
#' @param strandSpecific Numeric, 0 if no strand-specific read counting should be performed, 1 stranded, 2 reversely-stranded. 0 by default
#' @param isPairedEnd A logical, if paired-end reads are used, FALSE by default
#' @param DESeq2 A logical, if TRUE the RPKMs from exons and introns and associated variances are evaluated through the package DESeq2
#' @param temporalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the sequence of bam files submitted for the analysis, these labels characterize different files as replicates of a given conditions.
#' @return A list containing rpkms, counts and the annotation extracted from TxDB for exons and introns, if DESeq2 = TRUE the output also contains a set of data
#' needed to estimate rpkms variances.
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' paths_4su <- system.file('extdata', c('bamRep1.bam','bamRep2.bam','bamRep3.bam','bamRep4.bam') , package="INSPEcT")
#' paths_total <- system.file('extdata', c('bamRep1.bam','bamRep2.bam','bamRep3.bam','bamRep4.bam') , package="INSPEcT")
#' 
#' temporalDesign <- c(0,0,1,1)
#' 
#' #4sU analysis with DESeq2
#' makeRPKMsOut_4sU <- makeRPKMsFromBams(txdb=txdb,paths_foursu=paths_4su,paths_total=paths_total,temporalDesign=temporalDesign,DESeq2 = TRUE)
#' 
#' rpkms_4sU <- makeRPKMsOut_4sU$rpkms
#' counts_4sU <- makeRPKMsOut_4sU$counts
#' annotations_4sU <- makeRPKMsOut_4sU$annotations
#' dispersion_parameters_DESeq2_4sU <- makeRPKMsOut_4sU$dispersion_parameters_DESeq2
#' 
#' #4sU analysis without DESeq2
#' makeRPKMsOut_4sU <- makeRPKMsFromBams(txdb=txdb,paths_foursu=paths_4su,paths_total=paths_total,temporalDesign=temporalDesign,DESeq2 = FALSE)
#' 
#' rpkms_4sU <- makeRPKMsOut_4sU$rpkms
#' counts_4sU <- makeRPKMsOut_4sU$counts
#' annotations_4sU <- makeRPKMsOut_4sU$annotations
#' 
#' #No4sU analysis with DESeq2
#' makeRPKMsOut_No4sU <- makeRPKMsFromBams(txdb=txdb,paths_foursu=NULL,paths_total=paths_total,temporalDesign=temporalDesign,DESeq2 = TRUE)
#' 
#' rpkms_No4sU <- makeRPKMsOut_No4sU$rpkms
#' counts_No4sU <- makeRPKMsOut_No4sU$counts
#' annotations_No4sU <- makeRPKMsOut_No4sU$annotations
#' dispersion_parameters_DESeq2_No4sU <- makeRPKMsOut_No4sU$dispersion_parameters_DESeq2

makeRPKMsFromBams <- function(txdb
					, paths_foursu = NULL
					, paths_total
					, by = c('gene','tx')
					, countMultiMappingReads = FALSE
					, allowMultiOverlap = FALSE
					, strandSpecific = 0
					, isPairedEnd = FALSE
					, DESeq2 = TRUE
					, temporalDesign = NULL) 
{

	if(all(table(temporalDesign)==1))stop("makeRPKMs: at least one replicate is required.")
	counts2rpkms <- function(counts, widths, libsize) counts*(10^9/(widths[rownames(counts)]%o%libsize))

	## checks on paths_foursu
	if(is.null(paths_foursu)){
		print("makeRPKMs: No4sU mode.")
		No4sU <- TRUE
	}else{
		No4sU <- FALSE
		if( !is.character(paths_foursu) )
			stop('makeRPKMs: argument "paths_foursu" must be a character')
		# if( any(!file.exists(paths_foursu)) )
		# 	stop('makeRPKMs: at least one file specified in "paths_foursu" argument does not exist.')

		if( is.null(names(paths_foursu)) ) names(paths_foursu) <- paths_foursu
	}
	## checks on paths_total
	if( !is.character(paths_total) )
		stop('makeRPKMs: argument "paths_total" must be a character')
	# if( any(!file.exists(paths_total)) )
	# 	stop('makeRPKMs: at least one file specified in "paths_total" argument does not exist.')

	if( is.null(names(paths_total)) ) names(paths_total) <- paths_total

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
		iecounts <- lapply(paths_total, function(bamfile)
		{
			message(paste('##### - File:',bamfile,'- #####'))
			if( countMultiMappingReads ) {
				message('Importing bamfile...')
				if( isPairedEnd )
					samTab <- readGAlignmentPairs(bamfile)
				else
					samTab <- readGAlignments(bamfile)
			} else { # countMultiMappingReads==FALSE
				message('Importing bamfile...')
				if( isPairedEnd )
					samTab <- readGAlignmentPairs(bamfile, param=ScanBamParam(tagFilter=list('NH'=1)))
				else
					samTab <- readGAlignments(bamfile, param=ScanBamParam(tagFilter=list('NH'=1)))
				if( length(samTab)==0 ) stop('No alignments imported.')
			}
	
			if( strandSpecific == 2 ) samTab <- invertStrand(samTab)

			message('Counting reads on exon features...')
			foOut <- findOverlaps(exonsDB,samTab,ignore.strand=strandSpecific==0)
			onfeature <- unique(subjectHits(foOut))

			if( allowMultiOverlap ) {
				Unassigned_Ambiguity <- 0
				Assigned_Exons <- length(onfeature)
			} else {
				ambiguous_reads <- duplicated(subjectHits(foOut))|duplicated(subjectHits(foOut),fromLast=TRUE)
				Unassigned_Ambiguity <- length(unique(subjectHits(foOut)[ambiguous_reads]))
				Assigned_Exons <- length(which(!ambiguous_reads))
				foOut <- foOut[!ambiguous_reads]
			}
			exonCounts <- table(factor(queryHits(foOut), levels=1:queryLength(foOut)))
			names(exonCounts) <- names(exonsDB)

			message('Counting reads on intron features...')
			if( length(onfeature)>0 ) samTab <- samTab[-onfeature] # remove reads falling on exons
			foOut <- findOverlaps(intronsDB,samTab,ignore.strand=strandSpecific==0)
			onfeature <- unique(subjectHits(foOut))

			if( allowMultiOverlap ) {
				Unassigned_Ambiguity <- 0
				Assigned_Introns <- length(onfeature)
			} else {
				ambiguous_reads <- duplicated(subjectHits(foOut))|duplicated(subjectHits(foOut),fromLast=TRUE)
				Unassigned_Ambiguity <- Unassigned_Ambiguity + length(unique(subjectHits(foOut)[ambiguous_reads]))
				Assigned_Introns <- length(which(!ambiguous_reads))
				foOut <- foOut[!ambiguous_reads]
			}
			intronCounts <- table(factor(queryHits(foOut), levels=1:queryLength(foOut)))
			names(intronCounts) <- names(intronsDB)

			Unassigned_NoFeatures <- length(samTab) - length(onfeature)

			stat <- c(
					Unassigned_Ambiguity=Unassigned_Ambiguity,
					Assigned_Exons=Assigned_Exons,
					Assigned_Introns=Assigned_Introns,
					Unassigned_NoFeatures=Unassigned_NoFeatures
				)

			return(list(exonCounts=exonCounts, intronCounts=intronCounts, stat=stat))

		})

		allcounts <- lapply(c(exonCounts="exonCounts",intronCounts="intronCounts",stat="stat")
				, function(name) sapply(iecounts,'[[',name))

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
	
		allcounts <- lapply(list(foursu=paths_foursu, total=paths_total), function(files) {
	
			iecounts <- lapply(files, function(bamfile) {
	
				message(paste('##### - File:',bamfile,'- #####'))
				if( countMultiMappingReads ) {
					message('Importing bamfile...')
					if( isPairedEnd )
						samTab <- readGAlignmentPairs(bamfile)
					else
						samTab <- readGAlignments(bamfile)
				} else { 
					message('Importing bamfile...')
					if( isPairedEnd )
						samTab <- readGAlignmentPairs(bamfile, param=ScanBamParam(tagFilter=list('NH'=1)))
					else
						samTab <- readGAlignments(bamfile, param=ScanBamParam(tagFilter=list('NH'=1)))
					if( length(samTab)==0 ) stop('No alignments imported.')
				}
	
				if( strandSpecific == 2 ) samTab <- invertStrand(samTab)
	
				message('Counting reads on exon features...')
				foOut <- findOverlaps(exonsDB,samTab,ignore.strand=strandSpecific==0)
				onfeature <- unique(subjectHits(foOut))
	
				if( allowMultiOverlap ) {
					Unassigned_Ambiguity <- 0
					Assigned_Exons <- length(onfeature)
				} else {
					ambiguous_reads <- duplicated(subjectHits(foOut))|duplicated(subjectHits(foOut),fromLast=TRUE)
					Unassigned_Ambiguity <- length(unique(subjectHits(foOut)[ambiguous_reads]))
					Assigned_Exons <- length(which(!ambiguous_reads))
					foOut <- foOut[!ambiguous_reads]
				}
				exonCounts <- table(factor(queryHits(foOut), levels=1:queryLength(foOut)))
				names(exonCounts) <- names(exonsDB)
	
				message('Counting reads on intron features...')
				if( length(onfeature)>0 ) samTab <- samTab[-onfeature] # remove reads falling on exons
				foOut <- findOverlaps(intronsDB,samTab,ignore.strand=strandSpecific==0)
				onfeature <- unique(subjectHits(foOut))
	
				if( allowMultiOverlap ) {
					Unassigned_Ambiguity <- 0
					Assigned_Introns <- length(onfeature)
				} else {
					ambiguous_reads <- duplicated(subjectHits(foOut))|duplicated(subjectHits(foOut),fromLast=TRUE)
					Unassigned_Ambiguity <- Unassigned_Ambiguity + length(unique(subjectHits(foOut)[ambiguous_reads]))
					Assigned_Introns <- length(which(!ambiguous_reads))
					foOut <- foOut[!ambiguous_reads]
				}
				intronCounts <- table(factor(queryHits(foOut), levels=1:queryLength(foOut)))
				names(intronCounts) <- names(intronsDB)
	
				Unassigned_NoFeatures <- length(samTab) - length(onfeature)
	
				stat <- c(
						Unassigned_Ambiguity=Unassigned_Ambiguity,
						Assigned_Exons=Assigned_Exons,
						Assigned_Introns=Assigned_Introns,
						Unassigned_NoFeatures=Unassigned_NoFeatures
					)
	
				return(list(exonCounts=exonCounts, intronCounts=intronCounts, stat=stat))
	
			})
	
			return(lapply(c(exonCounts="exonCounts",intronCounts="intronCounts",stat="stat")
				, function(name) sapply(iecounts,'[[',name)))
	
		})
		
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
