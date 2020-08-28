#' Evaluate introns and exons expressions from BAM or SAM files
#' @description Given a TranscriptDb object and a list of BAM or SAM files 
#' "quantifyExpressionsFormBAMs" evaluates exons and introns expressions
#' and the associated variances per each gene.
#' @param txdb A TranscriptDB object
#' @param BAMfiles A vector of paths
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates.
#' Time points must be ordered according to the sequence of files submitted for the analysis, these labels characterize
#' different files as replicates of a given condition.
#' @param by A character, either "gene" or "tx", indicating if expressions and counts should be summarized at the levels of 
#' genes or transcripts. "gene" by default.
#' In case "tx" is selected, we suggest to set argument "allowMultiOverlap" to TRUE, otherwise the reads mapping to overlapping
#' transcripts of the same gene will remain unassigned.
#' @param countMultiMappingReads A logical, if multimapping reads should be counted, FALSE by default. Multimap reads are 
#' identified using the tag "NH" in the bam/sam file.
#' @param allowMultiOverlap A logical, indicating if a read is allowed to be assigned to more than one feature, FALSE by default
#' @param prioritizeExons A logical, indicating whether reads assigned to exon shold not be accounted for intron counts.
#' If set to FALSE, reads with shared overlap between an exon and the following intron will be assigned also to introns. This could improve
#' intronic quantification in experimental settings (including polyA library preparation) or compact genomes were intronic reads are 
#' sampled at a very low rate compared to exonic reads. By default, TRUE.
#' @param libsize A character, either "assigned" or "all", indicating whether the libsize for expression normalization should include all 
#' mapped reads or only the reads assigned to any of the features. By default, "assigned" is selected.
#' @param strandSpecific Numeric, 0 if no strand-specific read counting should be performed, 1 stranded, 2 reversely-stranded. 0 by default
#' @param isPairedEnd A logical, if paired-end reads are used, FALSE by default
#' @param DESeq2 A logical, if TRUE exons and introns variances are evaluated through the package DESeq2, if FALSE through plgem
#' @param varSamplingCondition A character reporting which experimental condition should be used to sample the variance if DESeq2 = FALSE. 
#' @param BPPARAM Parallelization parameters for bplapply. By default SerialParam()
#' By default, the first element of "experimentalDesign" with replicates.
#' @return A list containing expressions and associated variances for exons and introns.
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#'   require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#'   txdb<-TxDb.Mmusculus.UCSC.mm9.knownGene
#'   expDes<-c(0,0,1,1)
#'   
#'   paths_total<-system.file('extdata/', c('bamRep1.bam'
#'                                         ,'bamRep2.bam'
#'                                         ,'bamRep3.bam'
#'                                         ,'bamRep4.bam')
#'                           ,package='INSPEcT')
#'  
#'   matExp<-quantifyExpressionsFromBAMs(txdb=txdb
#'                                      ,BAMfiles=paths_total
#'                                      ,experimentalDesign=expDes)
#' }
quantifyExpressionsFromBAMs <- function(txdb
					, BAMfiles
					, experimentalDesign
					, by = c('gene','tx')
					, countMultiMappingReads = FALSE
					, allowMultiOverlap = FALSE
					, prioritizeExons = TRUE
					, libsize = c('assigned','all')
					, strandSpecific = 0
					, isPairedEnd = FALSE
					, DESeq2 = TRUE
					, varSamplingCondition = NULL
					, BPPARAM = SerialParam())
{

	############################################
	### CHECK ARGUMENTS ########################
	############################################

	# txdb
	if( class(txdb) != 'TxDb' ) 
		stop('quantifyExpressionsFromBAMs: "txdb" must be an object of TxDb class.')
	# BAMfiles
	if( any(!file.exists(BAMfiles)) )
	 	stop('quantifyExpressionsFromBAMs: at least one file specified in "BAMfiles" argument does not exist.')
	# experimentalDesign
	if(length(experimentalDesign)!=length(BAMfiles))
		stop('quantifyExpressionsFromBAMs: each bam file must be accounted in the experimentalDesign.')
	if(all(table(experimentalDesign)==1))
		stop("quantifyExpressionsFromBAMs: at least one replicate is required.")
	# by
	by <- by[1]
	if( !is.character(by) )
		stop('quantifyExpressionsFromBAMs: "by" must be either "tx" or "gene".')
	if( !( by %in% c('gene','tx') ) )
		stop('quantifyExpressionsFromBAMs: "by" must be either "tx" or "gene".')
	# countMultiMappingReads
	if( !is.logical(countMultiMappingReads) )
		stop('quantifyExpressionsFromBAMs: "countMultiMappingReads" must be a logical.')
	# allowMultiOverlap
	if( !is.logical(allowMultiOverlap) )
		stop('quantifyExpressionsFromBAMs: "allowMultiOverlap" must be a logical.')
	# prioritizeExons
	if( !is.logical(prioritizeExons) )
		stop('quantifyExpressionsFromBAMs: "prioritizeExons" must be a logical.')
	# strandSpecific
	if( !( strandSpecific %in% c(0,1,2) ) )
		stop('quantifyExpressionsFromBAMs: "strandSpecific" must be either a numeric between 0 and 2.')
	# isPairedEnd
	if( !is.logical(isPairedEnd) )
		stop('quantifyExpressionsFromBAMs: "isPairedEnd" must be a logical.')
	# DESeq2
	if( !is.logical(DESeq2) )
		stop('quantifyExpressionsFromBAMs: "DESeq2" must be a logical.')
	# varSamplingCondition
	if( !DESeq2 ) {
		if( is.null(varSamplingCondition) ) {
			varSamplingCondition <- names(which(table(experimentalDesign)>1)[1])
		} else {
			if( length(which(as.character(experimentalDesign) == varSamplingCondition)) < 2 )
				stop('quantifyExpressionsFromBAMs: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
		}
	} 

	############################################
	### MAKE ANNOTATION ########################
	############################################

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

		stop("quantifyExpressionsFromBAMs: argument 'by' not recognized.")

	}
	############################################
	### MAKE COUNTS FROM BAM ###################
	############################################

	if( is.null(names(BAMfiles)) ) {
		replicate_id <- unlist(lapply(split(experimentalDesign, experimentalDesign), seq_along))
		names(BAMfiles) <- paste(experimentalDesign, paste0('rep',replicate_id), sep='_')
	}

	iecounts <- bplapply(BAMfiles, function(bamfile)
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

			# match seqlevels based on the name (must be contained within) and 
			# sequence length (must be equal)
			samTab <- matchSeqnames(match_object=samTab, ref_object=exonsDB)

			if( strandSpecific == 2 ) samTab <- invertStrand(samTab)

			message('Counting reads on exon features...')
			foOut <- findOverlaps(samTab,exonsDB,ignore.strand=strandSpecific==0)
			onfeature <- unique(queryHits(foOut))

			if( allowMultiOverlap ) {
				Unassigned_Ambiguity <- 0
				Assigned_Exons <- length(onfeature)
			} else {
				ambiguous_reads <- duplicated(queryHits(foOut))|duplicated(queryHits(foOut),fromLast=TRUE)
				Unassigned_Ambiguity <- length(unique(queryHits(foOut)[ambiguous_reads]))
				Assigned_Exons <- length(which(!ambiguous_reads))
				foOut <- foOut[!ambiguous_reads]
			}
			exonCounts <- table(factor(subjectHits(foOut), levels=1:subjectLength(foOut)))
			names(exonCounts) <- names(exonsDB)

			message('Counting reads on intron features...')
			if( length(onfeature)>0 & prioritizeExons ) samTab <- samTab[-onfeature] # remove reads falling on exons
			foOut <- findOverlaps(samTab,intronsDB,ignore.strand=strandSpecific==0)
			onfeature <- unique(queryHits(foOut))

			if( allowMultiOverlap ) {
				Unassigned_Ambiguity <- 0
				Assigned_Introns <- length(onfeature)
			} else {
				ambiguous_reads <- duplicated(queryHits(foOut))|duplicated(queryHits(foOut),fromLast=TRUE)
				Unassigned_Ambiguity <- Unassigned_Ambiguity + length(unique(queryHits(foOut)[ambiguous_reads]))
				Assigned_Introns <- length(which(!ambiguous_reads))
				foOut <- foOut[!ambiguous_reads]
			}
			intronCounts <- table(factor(subjectHits(foOut), levels=1:subjectLength(foOut)))
			names(intronCounts) <- names(intronsDB)

			Unassigned_NoFeatures <- length(samTab) - length(onfeature)

			stat <- c(
					Unassigned_Ambiguity=Unassigned_Ambiguity,
					Assigned_Exons=Assigned_Exons,
					Assigned_Introns=Assigned_Introns,
					Unassigned_NoFeatures=Unassigned_NoFeatures
				)

			return(list(exonCounts=exonCounts, intronCounts=intronCounts, countsStats=stat))
	},BPPARAM=BPPARAM)
	allcounts <- lapply(c(exonsCounts="exonCounts",intronsCounts="intronCounts",countsStats="countsStats")
			, function(name) sapply(iecounts,'[[',name))
	libsize <- colSums(allcounts$countsStats[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])

	exonsWidths <- sapply(width(exonsDB),sum)
	intronsWidths <- sapply(width(intronsDB),sum)

	out <- quantifyExpressionsFromTrCounts(allcounts = allcounts
											 , experimentalDesign = experimentalDesign
											 , exonsWidths=exonsWidths
											 , intronsWidths=intronsWidths
											 , libsize=libsize
											 , DESeq2 = DESeq2
											 , varSamplingCondition = varSamplingCondition)
	out <- c( out, list(exonsWidths=exonsWidths, intronsWidths=intronsWidths), allcounts )
	return(out)
}


# look for possible sequence matches between two genomic annotaion objects
# where the name of the match sequence is included within the name 
# of the reference and with the same sequence length
matchSeqnames <- function(match_object, ref_object) {

	match_seqnames <- seqnames(seqinfo(match_object))
	ref_seqnames <- seqnames(seqinfo(ref_object))

	possible_name_matches <- lapply(match_seqnames, function(x) grep(x, ref_seqnames))

	match_seqlengths <- seqlengths(seqinfo(match_object))
	ref_seqlengths <- seqlengths(seqinfo(ref_object))
	possible_length_matches <- lapply(match_seqlengths, function(x) which(ref_seqlengths==x))

	both_matches <- lapply(seq_along(possible_name_matches), function(i) {
		intersect(possible_name_matches[[i]], possible_length_matches[[i]])
		})
	names(both_matches) <- match_seqnames

	good_matches <- sapply(both_matches, length) == 1
	seqlevels(match_object)[good_matches] <- ref_seqnames[unlist(both_matches[good_matches])]

	return(match_object)

}