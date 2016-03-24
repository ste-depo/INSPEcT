#' Calculate RPKM and count values on introns and exons from bam/sam files
#'
#' @description Given a TranscriptDb object and a list of bam/sam files for 4su and total RNA 
#' experiments, "makeRPKMs" function calculates read counts and RPKM on exonic and intronic 
#' features per each gene. Reads that fall where intronic and exonic features overlaps are 
#' univoquely assigned to exons.
#' @param txdb A TranscriptDB object
#' @param paths_foursu A vector of paths of 4sU-seq sam files
#' @param paths_total A vector of paths of RNA-seq sam files
#' @param by A character, either "gene" or "tx", indicating if rpkms and counts should be summarized at the levels of genes or transcripts. "gene" by default
#' @param countMultiMappingReads A logical, if multimapping reads should be counted, FALSE by default. Multimap reads are 
#' identified using the tag "NH" in the bam/sam file.
#' @param allowMultiOverlap A logical, indicating if a read is allowed to be assigned to more than one feature, FALSE by default
#' @param strandSpecific A logical, if strand-specific read counting should be performed, FALSE by default
#' @param isPairedEnd A logical, if paired-end reads are used, FALSE by default
#' @return A list containing rpkms, counts and the annotation extracted from TxDB for exons and introns
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' files4su <- system.file('extdata', '4sURNA_0h.bam', package="INSPEcT")
#' filesTotal <- system.file('extdata', 'totalRNA_0h.bam', package="INSPEcT")
#' makeRPKMsOut <- makeRPKMs(txdb, files4su, filesTotal)
#' rpkms <- makeRPKMsOut$rpkms
#' counts <- makeRPKMsOut$counts
#' annotation <- makeRPKMsOut$annotation
makeRPKMs <- function(txdb, paths_foursu, paths_total, by=c('gene','tx'),
	countMultiMappingReads=FALSE,allowMultiOverlap=FALSE,strandSpecific=FALSE,isPairedEnd=FALSE) 
{

	## checks on paths_foursu
	if( !is.character(paths_foursu) )
		stop('makeRPKMs: argument "paths_foursu" must be a character')
	if( any(!file.exists(paths_foursu)) )
		stop('makeRPKMs: at least one file specified in "paths_foursu" argument does not exist.')
	## checks on paths_total
	if( !is.character(paths_total) )
		stop('makeRPKMs: argument "paths_total" must be a character')
	if( any(!file.exists(paths_total)) )
		stop('makeRPKMs: at least one file specified in "paths_total" argument does not exist.')

	if( is.null(names(paths_foursu)) ) names(paths_foursu) <- paths_foursu
	if( is.null(names(paths_total)) ) names(paths_total) <- paths_total

	message('Generating annotation from txdb...')

	by <- by[1]
	if( by=="gene" ) {

		exonsDB <- reduce(exonsBy(txdb ,'gene'))
		exonsDB <- exonsDB[elementLengths(range(exonsDB))==1]
		intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
		intronsDB <- intronsDB[elementLengths(intronsDB)>0]

	} else if( by=="tx" ){

		exonsDB <- exonsBy(txdb ,'tx', use.names=TRUE)
		intronsDB <- intronsByTranscript(txdb, use.names=TRUE)
		intronsDB <- intronsDB[elementLengths(intronsDB)>0]

	} else {

		stop("makeRPKMs: argument 'by' not recognized.")

	}

	allcounts <- lapply(list(foursu=paths_foursu, total=paths_total), function(files) {

		iecounts <- lapply(files, function(bamfile) {

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

			message('Counting reads on exon features...')
			foOut <- findOverlaps(exonsDB,samTab,ignore.strand=!strandSpecific)
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
			foOut <- findOverlaps(intronsDB,samTab,ignore.strand=!strandSpecific)
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

	counts2rpkms <- function(counts, widths, libsizes) t(t(counts/widths)/libsizes)*10^9

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

	return(list(rpkms=allrpkms, counts=allcounts, annotation=list(exon=exonsDB, intron=intronsDB)))

}
