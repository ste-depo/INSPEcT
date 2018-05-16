#' Calculate Expression and count values on introns and exons from bam/sam files
#' @description Given a TranscriptDb object and a list of bam/sam files for total and eventually RNA 
#' experiments, "makeExpressionsFromBAMs" function calculates read counts and Expression on exonic and intronic 
#' features per each gene. Reads that fall where intronic and exonic features overlaps are 
#' univoquely assigned to exons.
#' @param txdb A TranscriptDB object
#' @param paths_foursu A vector of paths of Nascent-seq sam files
#' @param BAMfiles A vector of paths of RNA-seq sam files
#' @param by A character, either "gene" or "tx", indicating if expressions and counts should be summarized at the levels of genes or transcripts. "gene" by default
#' @param countMultiMappingReads A logical, if multimapping reads should be counted, FALSE by default. Multimap reads are 
#' identified using the tag "NH" in the bam/sam file.
#' @param allowMultiOverlap A logical, indicating if a read is allowed to be assigned to more than one feature, FALSE by default
#' @param strandSpecific Numeric, 0 if no strand-specific read counting should be performed, 1 stranded, 2 reversely-stranded. 0 by default
#' @param isPairedEnd A logical, if paired-end reads are used, FALSE by default
#' @param DESeq2 A logical, if TRUE the Expressions from exons and introns and associated variances are evaluated through the package DESeq2
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the sequence of bam files submitted for the analysis, these labels characterize different files as replicates of a given conditions.
#' @return A list containing expressions, counts and the annotation extracted from TxDB for exons and introns, if DESeq2 = TRUE the output also contains a set of data
#' needed to estimate expressions variances.
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' 
#' paths_Nascent <- system.file('extdata', c('bamRep1.bam','bamRep2.bam','bamRep3.bam','bamRep4.bam') , package="INSPEcT")
#' BAMfiles <- system.file('extdata', c('bamRep1.bam','bamRep2.bam','bamRep3.bam','bamRep4.bam') , package="INSPEcT")
#' 
#' experimentalDesign <- c(0,0,1,1)
#' 
#' #Nascent analysis with DESeq2
#' makeExpressionsOut_Nascent <- makeExpressionsFromBAMs(txdb=txdb,paths_foursu=paths_Nascent,BAMfiles=BAMfiles,experimentalDesign=experimentalDesign,DESeq2 = TRUE)
#' 
#' expressions_Nascent <- makeExpressionsOut_Nascent$expressions
#' counts_Nascent <- makeExpressionsOut_Nascent$counts
#' annotations_Nascent <- makeExpressionsOut_Nascent$annotations
#' dispersion_parameters_DESeq2_Nascent <- makeExpressionsOut_Nascent$dispersion_parameters_DESeq2
#' 
#' #Nascent analysis without DESeq2
#' makeExpressionsOut_Nascent <- makeExpressionsFromBAMs(txdb=txdb,paths_foursu=paths_Nascent,BAMfiles=BAMfiles,experimentalDesign=experimentalDesign,DESeq2 = FALSE)
#' 
#' expressions_Nascent <- makeExpressionsOut_Nascent$expressions
#' counts_Nascent <- makeExpressionsOut_Nascent$counts
#' annotations_Nascent <- makeExpressionsOut_Nascent$annotations
#' 
#' #NoNascent analysis with DESeq2
#' makeExpressionsOut_NoNascent <- makeExpressionsFromBAMs(txdb=txdb,paths_foursu=NULL,BAMfiles=BAMfiles,experimentalDesign=experimentalDesign,DESeq2 = TRUE)
#' 
#' expressions_NoNascent <- makeExpressionsOut_NoNascent$expressions
#' counts_NoNascent <- makeExpressionsOut_NoNascent$counts
#' annotations_NoNascent <- makeExpressionsOut_NoNascent$annotations
#' dispersion_parameters_DESeq2_NoNascent <- makeExpressionsOut_NoNascent$dispersion_parameters_DESeq2

quantifyExpressionsFromBAMs <- function(txdb
					, BAMfiles
					, experimentalDesign
					, by = c('gene','tx')
					, countMultiMappingReads = FALSE
					, allowMultiOverlap = FALSE
					, strandSpecific = 0
					, isPairedEnd = FALSE
					, DESeq2 = TRUE
					, varSamplingCondition = NULL
					, plgemFits = NULL
					, returnPlgemFits = FALSE)
{

	############################################
	### CHECK ARGUMENTS ########################
	############################################
	if( !is.logical(DESeq2) & !any(as.character(experimentalDesign)==varSamplingCondition) & is.null(plgemFits))
		stop('makeExpressions: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
	if(length(experimentalDesign)!=length(BAMfiles))
		stop('makeExpressions: each bam file must be accounted in the experimentalDesign')
	if(all(table(experimentalDesign)==1) & is.null(plgemFits))
		stop("makeExpressions: at least one replicate is required.")
	if( any(!file.exists(BAMfiles)) )
	 	stop('makeExpressions: at least one file specified in "BAMfiles" argument does not exist.')

	if( is.null(names(BAMfiles)) ) names(BAMfiles) <- BAMfiles

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

		stop("makeExpressions: argument 'by' not recognized.")

	}

	############################################
	### MAKE COUNTS FROM BAM ###################
	############################################

	iecounts <- lapply(BAMfiles, function(bamfile)
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
	allcounts <- lapply(c(exonsCounts="exonCounts",intronsCounts="intronCounts",stat="stat")
			, function(name) sapply(iecounts,'[[',name))
	libsize <- colSums(allcounts$stat[c('Assigned_Exons','Assigned_Introns'),,drop=FALSE])

	exonsWidths <- sapply(width(exonsDB),sum)
	intronsWidths <- sapply(width(intronsDB),sum)

	return(quantifyExpressionsFromTrCounts(libsize=libsize
										 , exonsWidths=exonsWidths
										 , intronsWidths=intronsWidths
										 , allcounts = allcounts
										 , by = by
										 , DESeq2 = DESeq2
										 , experimentalDesign = experimentalDesign
										 , varSamplingCondition = varSamplingCondition
										 , plgemFits = plgemFits
										 , returnPlgemFits = returnPlgemFits))

}
