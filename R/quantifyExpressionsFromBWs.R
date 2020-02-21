#' Evaluate introns and exons expressions from BAM or SAM files
#' @description Given a TranscriptDb object and a list of bigWig (BW) files 
#' "quantifyExpressionsFormBWs" evaluates exons and introns expressions
#' and the associated variances per each gene.
#' @param txdb A TranscriptDB object
#' @param BWfiles A vector of paths
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates.
#' Time points must be ordered according to the sequence of files submitted for the analysis, these labels characterize
#' different files as replicates of a given condition.
#' @param readLength A numerical that indicates the read length of the RNA-seq experiment. Used to normalize the coverage. By default, 50.
#' @param by A character, either "gene" or "tx", indicating if expressions and counts should be summarized at the levels of 
#' genes or transcripts. "gene" by default.
#' In case "tx" is selected, we suggest to set argument "allowMultiOverlap" to TRUE, otherwise the reads mapping to overlapping
#' transcripts of the same gene will remain unassigned.
#' @param allowMultiOverlap A logical, indicating if a read is allowed to be assigned to more than one feature, FALSE by default
#' @param prioritizeExons A logical, indicating whether reads assigned to exon shold not be accounted for intron counts.
#' If set to FALSE, reads with shared overlap between an exon and the following intron will be assigned to the intron. This could improve
#' intronic quantification in experimental settings (including polyA library preparation) or compact genomes were intronic reads are 
#' sampled at a very low rate compared to exonic reads. By default, TRUE.
#' @param libsize A character, either "assigned" or "all", indicating whether the libsize for expression normalization should include all 
#' mapped reads or only the reads assigned to any of the features. By default, "assigned" is selected.
#' @param DESeq2 A logical, if TRUE exons and introns variances are evaluated through the package DESeq2, if FALSE through plgem
#' @param varSamplingCondition A character reporting which experimental condition should be used to sample the variance if DESeq2 = FALSE. 
#' @param BPPARAM Parallelization parameters for bplapply. By default SerialParam()
#' By default, the first element of "experimentalDesign" with replicates.
#' @return A list containing expressions and associated variances for exons and introns.
quantifyExpressionsFromBWs <- function(txdb
																			 , BWfiles
																			 , experimentalDesign
																			 , readLength = 50
																			 , by = c('gene','tx')
																			 , libsize = c('assigned','all')
																			 , DESeq2 = TRUE
																			 , varSamplingCondition = NULL
																			 , BPPARAM = SerialParam())
{
	
	############################################
	### CHECK ARGUMENTS ########################
	############################################
	
	# txdb
	if( class(txdb) != 'TxDb' ) 
		stop('quantifyExpressionsFromBWs: "txdb" must be an object of TxDb class.')
	# BWfiles
	if( any(!file.exists(BWfiles)) )
		stop('quantifyExpressionsFromBWs: at least one file specified in "BWfiles" argument does not exist.')
	# experimentalDesign
	if(length(experimentalDesign)!=length(BWfiles))
		stop('quantifyExpressionsFromBWs: each bam file must be accounted in the experimentalDesign.')
	if(all(table(experimentalDesign)==1))
		warning("quantifyExpressionsFromBWs: at least one replicate is required.")
	# by
	by <- by[1]
	if( !is.character(by) )
		stop('quantifyExpressionsFromBWs: "by" must be either "tx" or "gene".')
	if( !( by %in% c('gene','tx') ) )
		stop('quantifyExpressionsFromBWs: "by" must be either "tx" or "gene".')
	# DESeq2
	if( !is.logical(DESeq2) )
		stop('quantifyExpressionsFromBWs: "DESeq2" must be a logical.')
	# varSamplingCondition
	if( !DESeq2 ) {
		if( is.null(varSamplingCondition) ) {
			varSamplingCondition <- names(which(table(experimentalDesign)>1)[1])
		} else {
			if( length(which(as.character(experimentalDesign) == varSamplingCondition)) < 2 )
				stop('quantifyExpressionsFromBWs: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
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
		
		stop("quantifyExpressionsFromBWs: argument 'by' not recognized.")
		
	}
	############################################
	### MAKE COUNTS FROM BAM ###################
	############################################
	
	if( is.null(names(BWfiles)) ) {
		if( any(table(experimentalDesign)>1) ) {
			replicate_id <- unlist(lapply(split(experimentalDesign, experimentalDesign), seq_along))
			names(BWfiles) <- paste(experimentalDesign, paste0('rep',replicate_id), sep='_')
		} else {
			names(BWfiles) <- experimentalDesign
		}
	}
	# BWfiles <- BigWigFileList(BWfiles)
	
	iecounts <- bplapply(BWfiles, function(bwfile)
	{
		message(paste('##### - File:',bwfile,'- #####'))
		coverage <- import(bwfile, as = 'RleList')
		# normalize the coverege by the read length
		coverage <- coverage/readLength
		
		# check compatibility between BigWig and annotation
		bigWigSeqLengths <- sapply(coverage, length)
		annotationSeqLengths <- seqlengths(seqinfo(exonsDB))
		commonSeq <- intersect(names(bigWigSeqLengths), names(annotationSeqLengths))
		if( length(commonSeq) == 0 ) 
			stop('quantifyExpressionsFromBWs: annotation and bigwig file have no common sequences')
		if( !identical(bigWigSeqLengths[commonSeq], annotationSeqLengths[commonSeq]) ) 
			stop('quantifyExpressionsFromBWs: annotation and bigwig file have different sequence lengths')

		message('Counting reads on exon features...')
		exonCounts <- round(unlist(lapply(names(coverage), function(chr) {
			exonsDBchr <- exonsDB[seqnames(exonsDB) == chr]
			counts <- sum(Views(coverage[[chr]], ranges(unlist(exonsDBchr))))
			tapply(counts, names(counts), sum)
		})))
		# names(exonCounts) <- names(exonsDB)
		Assigned_Exons <- sum(exonCounts)
		
		message('Counting reads on intron features...')
		intronCounts <- round(unlist(lapply(names(coverage), function(chr) {
			intronsDBchr <- intronsDB[seqnames(intronsDB) == chr]
			counts <- sum(Views(coverage[[chr]], ranges(unlist(intronsDBchr))))
			tapply(counts, names(counts), sum)
		})))
		# names(intronCounts) <- names(intronsDB)
		Assigned_Introns <- sum(intronCounts)
		
		stat <- c(
			Unassigned_Ambiguity=0,
			Assigned_Exons=Assigned_Exons,
			Assigned_Introns=Assigned_Introns,
			Unassigned_NoFeatures=NA,
			Libsize=round(sum(sapply(coverage, sum)))
		)
		
		return(list(exonCounts=exonCounts, intronCounts=intronCounts, countsStats=stat))
	},BPPARAM=BPPARAM)
	allcounts <- lapply(c(exonsCounts="exonCounts",intronsCounts="intronCounts",countsStats="countsStats")
											, function(name) sapply(iecounts,'[[',name))
	libsize <- allcounts$countsStats['Libsize',]
	
	exonsWidths <- sapply(width(exonsDB),sum)[rownames(allcounts$exonsCounts)]
	intronsWidths <- sapply(width(intronsDB),sum)[rownames(allcounts$intronsCounts)]
	
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

# importFrom("rtracklayer", "BigWigFileList", "import")
# importFrom("IRanges", "Views")