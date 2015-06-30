#' Calculate RPKM values on introns and exons from sam files
#'
#' @description Given a vector of time points, paths to sam files of the 4sU-seq and RNA-seq experiments and paths of GTF files, "makeRPKMs" function calculates read counts and RPKM on exonic and intronic features per each gene. Reads that fall where intronic and exonic features overlaps are assigned to exons.
#' This function is based on HTSeq software (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html) which is expected to pre installed and on the path.
#' @param tpts A vector of time points corresponding in which experiments have been collected (it admits non-unique values in case of replicates)
#' @param sampaths_4su A vector of paths of 4sU-seq sam files of the same length of tpts
#' @param sampaths_total A vector of paths of RNA-seq sam files of the same length of tpts
#' @param exons_gtf_file A path to the GTF file that defines exon features
#' @param introns_gtf_file A path to the GTF file that defines intron features (if NULL only exons reads are quantified)
#' @param out_folder A path to the output folder (by default the current directory)
#' @return A list of intron and exon rpkms for both RNA-seq and 4sU-seq
#' @examples
#' exons_gtf_file <- system.file('extdata', 'UCSC_mm9_genes_exons.txt', package="INSPEcT")
#' introns_gtf_file <- system.file('extdata', 'UCSC_mm9_genes_introns.txt', package="INSPEcT")
#' samfiles4su <- system.file('extdata', '4sURNA_0h.sam', package="INSPEcT")
#' samfilesTotal <- system.file('extdata', 'totalRNA_0h.sam', package="INSPEcT")
#' \dontrun{
#' rpkms <- makeRPKMs(1, samfiles4su, samfilesTotal, exons_gtf_file, introns_gtf_file)
#' save(rpkms, file='rpkms.rd')
#' }
makeRPKMs <- function(tpts, sampaths_4su, sampaths_total
	, exons_gtf_file, introns_gtf_file, out_folder='.'
	)
{

	################
	## check arguments
	######################
	## checks on tpts
	if( !is.numeric(tpts) )
		stop('makeRPKMs: argument "tpts" must be numeric.')
	if( length(tpts) != length(sampaths_4su) )
		stop('makeRPKMs: argument "tpts" and "sampaths_4su" must be of the same length.')
	if( length(tpts) != length(sampaths_total) )
		stop('makeRPKMs: argument "tpts" and "sampaths_total" must be of the same length.')
	## checks on sampaths_4su
	if( !is.character(sampaths_4su) )
		stop('makeRPKMs: argument "sampaths_4su" must be a character')
	if( any(!file.exists(sampaths_4su)) )
		stop('makeRPKMs: at least one file specified in "sampaths_4su" argument does not exist.')
	## checks on sampaths_total
	if( !is.character(sampaths_total) )
		stop('makeRPKMs: argument "sampaths_total" must be a character')
	if( any(!file.exists(sampaths_total)) )
		stop('makeRPKMs: at least one file specified in "sampaths_total" argument does not exist.')
	## checks on exons_gtf_file
	if( !is.character(exons_gtf_file) )
		stop('makeRPKMs: argument "exons_gtf_file" must be a character.')
	if( length(exons_gtf_file) != 1 )
		stop('makeRPKMs: argument "exons_gtf_file" must be of length one.')
	if( !file.exists(exons_gtf_file) )
		stop('makeRPKMs: file specified in "exons_gtf_file" argument does not exist.')
	## checks on introns_gtf_file
	if( !is.character(introns_gtf_file) )
		stop('makeRPKMs: argument "introns_gtf_file" must be a character.')
	if( length(introns_gtf_file) != 1 )
		stop('makeRPKMs: argument "introns_gtf_file" must be of length one.')
	if( !file.exists(introns_gtf_file) )
		stop('makeRPKMs: file specified in "introns_gtf_file" argument does not exist.')

	multiGtfHTSeq <- file.path(system.file(package="INSPEcT"), 'multiGtfHTSeq.py')

	## distinguish between the exons only mode and introns-exons mode
	if( is.null(introns_gtf_file) ) {
		intronsMode <- FALSE
		message('Introns GTF file have not been provided, running in exons only mode.')
	} else intronsMode <- TRUE

	## give a basename for the feature counts according to the gtf they are comeing from
	exons_gtf_file_basename <- 'exons'
	if( intronsMode )
		introns_gtf_file_basename <- 'introns'

	## in case of replicated time points add an attitional lable to make 
	# unique directories where data relative to each time points will be stored
	if( any(duplicated(tpts)) )
	{
		replicate_id <- rep(NA, length(tpts))
		id <- 1
		while( any(is.na(replicate_id)) )
		{
			ix <- is.na(replicate_id)
			replicate_id[ix][!duplicated(tpts[ix])] <- id
			id <- id+1
		}
		replicate_id[!duplicated(tpts)] <- 1
	} else {
		replicate_id <- NULL
	}

	# #################################
	# load introns and exons features
	# ################################

	extractGeneId <- function(id_field)
		sapply(
			lapply(strsplit(id_field, '; |;'), strsplit, split=' '), 
				function(x) 
					x[[which(sapply(x, function(x) x[1]=='gene_id'))]][2]
			)

	# exon features
	message('Retrieving and processing exon features')
	exonsDB <- read.table(exons_gtf_file, sep='\t', as.is=TRUE)#, nrows=nrows)
	colnames(exonsDB) <- c('chr','genome','feature','start','end'
		, 'coord','strand', '?','gene_id')
	exonsDB$gene_id <- extractGeneId(exonsDB$gene_id)
	exonsDBsplit <- split(exonsDB, exonsDB$gene_id)
	exonWidth <- sapply(exonsDBsplit, function(x) sum(x$end-x$start+1))

	# intron features, only if intronsMode is ON
	message('Retrieving and processing intron features')
	if( intronsMode ) {
		intronsDB <- read.table(introns_gtf_file, sep='\t', as.is=TRUE)#, nrows=nrows)
		colnames(intronsDB) <- c('chr','genome','feature','start','end'
			, 'coord','strand', '?','gene_id')
		intronsDB$gene_id <- extractGeneId(intronsDB$gene_id)
		intronsDBsplit <- split(intronsDB, intronsDB$gene_id)
		intronWidth <- sapply(intronsDBsplit, function(x) sum(x$end-x$start+1))
	}

	# #############
	# 4sU samfiles
	# ##############
	
	sampleNames <- paste('4sU', 't', signif(tpts,2), sep='_')
	if( !is.null(replicate_id) ) 
		sampleNames <- paste(sampleNames, replicate_id, sep='_')
	multiGtfHTSeq_4su_out_folders <- file.path(out_folder, sampleNames)
	# count reads on different gtf files giving priorities to
	# different GTF files: if a read falls into a feature in a gtf
	# of group 1 and group 2 is counted only as overlapping the 
	# feature of group 1
	for( i in 1:length(tpts) )
	{
		execute <- paste(
			'python '
			, multiGtfHTSeq
			, ' -o '
			, multiGtfHTSeq_4su_out_folders[i]
			, ' --exons='
			, exons_gtf_file #, ','
			, if( intronsMode ) ' --introns=' else ''
			, if( intronsMode ) introns_gtf_file else ''
			, ' '
			, sampaths_4su[i]
			, sep=''
			)
		message(paste('Processing file', sampaths_4su[i]))
		message(execute)
		system(execute)
	}
	# library sizes
	libsize_files <- file.path(multiGtfHTSeq_4su_out_folders
		, 'general_stats.counts')
	rtOut <- lapply(libsize_files, read.table)
	# nLines <- nrow(rtOut[[1]])
	rowNames <- rtOut[[1]]$V1
	libsize <- data.frame(
		all=sapply(rtOut, function(x) sum(x[rowNames %in% c('feature', 'ambiguous', 'no_feature'),2]))
		, feature=sapply(rtOut, function(x) x[rowNames %in% 'feature',2])
		, ambiguous=sapply(rtOut, function(x) x[rowNames %in% 'ambiguous',2])
		, nofeature=sapply(rtOut, function(x) x[rowNames %in% 'no_feature',2])
		, introns=if( intronsMode ) sapply(rtOut, function(x) x[rowNames %in% introns_gtf_file_basename,2]) else 0
		, exons=sapply(rtOut, function(x) x[rowNames %in% exons_gtf_file_basename,2])
		)
	rownames(libsize) <- sampleNames
	#### exons ####
	# counts
	exons_files <- file.path(multiGtfHTSeq_4su_out_folders
		, paste(exons_gtf_file_basename, '_gr1.counts', sep=''))
		#, '/exons_mm9_gr1.counts', sep='' )
	rtOut <- lapply(exons_files, read.table)
	exon_counts <- sapply(rtOut, function(x) x[,2])
	featureNames <- as.character(rtOut[[1]][,1])
	rownames(exon_counts) <- featureNames
	colnames(exon_counts) <- sampleNames
	# exons rpkms
	exonWidth <- exonWidth[featureNames]
	rpkms_4su_exons <- t(t(exon_counts/exonWidth)/libsize$feature)*10^9
	#### introns ####
	if( intronsMode ) {
		# counts
		introns_files <- file.path(multiGtfHTSeq_4su_out_folders
			, paste(introns_gtf_file_basename, '_gr2.counts', sep=''))
		# introns_files <- paste(multiGtfHTSeq_4su_out_folders
		# 	, '/introns_mm9_gr2.counts', sep='' )
		rtOut <- lapply(introns_files, read.table)
		intron_counts <- sapply(rtOut, function(x) x[,2])
		featureNames <- as.character(rtOut[[1]][,1])
		rownames(intron_counts) <- featureNames
		colnames(intron_counts) <- sampleNames
		# introns rpkms
		intronWidth <- intronWidth[featureNames]
		rpkms_4su_introns <- t(t(intron_counts/intronWidth )/libsize$feature)*10^9		
	}
	# write counts and rpkms on file
	write.table(libsize, file=file.path(out_folder, 'libsize_4su.txt'), sep='\t')
	write.table(exon_counts, file=file.path(out_folder, 'counts_4su_exons.txt'), sep='\t')
	write.table(rpkms_4su_exons, file=file.path(out_folder, 'rpkms_4su_exons.txt'), sep='\t')
	if( intronsMode ) {
		write.table(intron_counts, file=file.path(out_folder, 'counts_4su_introns.txt'), sep='\t')
		write.table(rpkms_4su_introns, file=file.path(out_folder, 'rpkms_4su_introns.txt'), sep='\t')		
	}

	# ##################
	# totalRNA samfiles
	# ##################

	sampleNames <- paste('total', 't', signif(tpts,2), sep='_')
	if( !is.null(replicate_id) ) 
		sampleNames <- paste(sampleNames, replicate_id, sep='_')
	multiGtfHTSeq_total_out_folders <- file.path(out_folder, sampleNames)
	# count reads on different gtf files giving priorities to
	# different GTF files: if a read falls into a feature in a gtf
	# of group 1 and group 2 is counted only as overlapping the 
	# feature of group 1
	for( i in 1:length(tpts) )
	{
		execute <- paste(
			'python '
			, multiGtfHTSeq
			, ' -o '
			, multiGtfHTSeq_total_out_folders[i]
			, ' --exons='
			, exons_gtf_file #, ','
			, if( intronsMode ) ' --introns=' else ''
			, if( intronsMode ) introns_gtf_file else ''
			, ' '
			, sampaths_total[i]
			, sep=''
			)
		message(paste('Processing file', sampaths_total[i]))
		message(execute)
		system(execute)
	}
	# library sizes
	libsize_files <- file.path(multiGtfHTSeq_total_out_folders
		, 'general_stats.counts')
	rtOut <- lapply(libsize_files, read.table)
	libsize <- data.frame(
		all=sapply(rtOut, function(x) sum(x[rowNames %in% c('feature', 'ambiguous', 'no_feature'),2]))
		, feature=sapply(rtOut, function(x) x[rowNames %in% 'feature',2])
		, ambiguous=sapply(rtOut, function(x) x[rowNames %in% 'ambiguous',2])
		, nofeature=sapply(rtOut, function(x) x[rowNames %in% 'no_feature',2])
		, introns=if( intronsMode ) sapply(rtOut, function(x) x[rowNames %in% introns_gtf_file_basename,2]) else 0
		, exons=sapply(rtOut, function(x) x[rowNames %in% exons_gtf_file_basename,2])
		)
	rownames(libsize) <- sampleNames
	#### exons ####
	# counts
	exons_files <- file.path(multiGtfHTSeq_total_out_folders
		, paste(exons_gtf_file_basename, '_gr1.counts', sep=''))
	rtOut <- lapply(exons_files, read.table)
	exon_counts <- sapply(rtOut, function(x) x[,2])
	featureNames <- as.character(rtOut[[1]][,1])
	rownames(exon_counts) <- featureNames
	colnames(exon_counts) <- sampleNames
	# exons rpkms
	exonWidth <- exonWidth[featureNames]
	rpkms_total_exons <- t(t(exon_counts/exonWidth)/libsize$feature)*10^9
	#### introns ####
	if( intronsMode ) {
		# counts
		introns_files <- file.path(multiGtfHTSeq_total_out_folders
			, paste(introns_gtf_file_basename, '_gr2.counts', sep=''))
		rtOut <- lapply(introns_files, read.table)
		intron_counts <- sapply(rtOut, function(x) x[,2])
		featureNames <- as.character(rtOut[[1]][,1])
		rownames(intron_counts) <- featureNames
		colnames(intron_counts) <- sampleNames
		# introns rpkms
		intronWidth <- intronWidth[featureNames]
		rpkms_total_introns <- t(t(intron_counts/intronWidth)/libsize$feature)*10^9
	}
	# write counts and rpkms on file
	write.table(libsize, file=file.path(out_folder, 'libsize_total.txt'), sep='\t')
	write.table(exon_counts, file=file.path(out_folder, 'counts_total_exons.txt'), sep='\t')
	write.table(rpkms_total_exons, file=file.path(out_folder, 'rpkms_total_exons.txt'), sep='\t')
	if( intronsMode ) {
		write.table(intron_counts, file=file.path(out_folder, 'counts_total_introns.txt'), sep='\t')
		write.table(rpkms_total_introns, file=file.path(out_folder, 'rpkms_total_introns.txt'), sep='\t')
	}

	return(list(
		rpkms_4su_exons=rpkms_4su_exons
		, rpkms_4su_introns=if( intronsMode ) rpkms_4su_introns else NULL
		, rpkms_total_exons=rpkms_total_exons
		, rpkms_total_introns=if( intronsMode ) rpkms_total_introns else NULL
		))

}
