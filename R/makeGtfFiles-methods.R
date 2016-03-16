#' Make a GTF file of exonic or intronic features from a TxDb object
#'
#' @name makeGtfFromDb
NULL

#' @rdname makeGtfFromDb
#'
#' @param object An object of class TxDb
#' @param type Either "gene" or "tx", see details for further information
#' @param filename The name of the GTF file that is created
#' @details In case type "gene" is chosen, the union of all exons of transcripts associated to a gene is used.
#' @return None
#' @examples
#' require(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' isActiveSeq(txdb) <- c(TRUE, rep(FALSE, length(isActiveSeq(txdb)) - 1))
#' res <- makeExonsGtfFromDb(txdb, 'gene', 'exons.gtf')
#' isActiveSeq(txdb) <- rep(TRUE, length(isActiveSeq(txdb)))
setMethod( 'makeExonsGtfFromDb' , 'TxDb', function( object , type, filename ) {
	# controls on arguments
	if( missing( type ) )
		stop('makeExonsGtfFromDb: type has to be specified.')
	if( missing( filename ) )
		stop('makeExonsGtfFromDb: filename has to be specified.')
	if( type != 'gene' & type != 'tx' )
		stop('makeExonsGtfFromDb: type has to be either "gene" or "tx".')
	if( ! is.character( filename ) )
		stop('makeExonsGtfFromDb: filename has to be character class.')
	# use method disjointExons to get a GRangeList of disjoint exons 
	message( 'Retrieving data from TxDb...' )
	if( type == 'gene' ) {
		# get a GRangeList of exons belonging to each gene
		exonsDB <- exonsBy( object , 'gene' )
		# reduce each gene ( merge overlapping exons ) and transform to GRanges
		exonsDB <- unlist( reduce( exonsDB ) )
		# feature name to be passed later as annotation
		featureName <- 'gene_id "'
	} else if ( type == 'tx' ) {
		# get a GRangeList of exons belonging to each transcript
		exonsDB <- exonsBy( object , by='tx' , use.names=TRUE )
		# pass from GRangesList to GRanges
		exonsDB <- unlist( exonsDB )
		# feature name to be passed later as annotation
		featureName <- 'tx_name "'
	}
	# make a dataframe from the GRangeList
	exonsDF <- data.frame(
		# chromosome
		chromosome = seqnames( exonsDB )
		# source ( name of the program who generated the feature )
		, source = object$packageName
		# feature 
		, feature = 'exon'
		# start
		, start = start( exonsDB )
		# end
		, end = end( exonsDB )
		# score ( a floating point value )
		, score = '.'
		# strand
		, strand = strand( exonsDB )
		# frame ( One of '0', '1' or '2'. '0' indicates that the first base 
		# of the feature is the first base of a codon, '1' that the second 
		# base is the first base of a codon, and so on.. )
		, frame = '.'
		# attribute ( A semicolon-separated list of tag-value pairs, providing
		# additional information about each feature. )
		, attribute = paste( featureName , names( exonsDB ) , '"' , sep='' )
	)
	# sorting lines according to genomic position
	message('Sorting according to genomic position...')
	idx <- order( exonsDF$chromosome , exonsDF$start )
	exonsDF <- exonsDF[ idx , ]
	# write the GTF file
	write.table(
		exonsDF 
		, file  = filename 
		, quote = FALSE 
		, sep   = '\t'
		, row.names = FALSE 
		, col.names = FALSE
		)
	message('File created.')
	})

#' @rdname makeGtfFromDb
setMethod( 'makeIntronsGtfFromDb' , 'TxDb', function( object , type, filename ) {
	# controls on arguments
	if( missing( type ) )
		stop('makeIntronsGtfFromDb: type has to be specified.')
	if( missing( filename ) )
		stop('makeIntronsGtfFromDb: filename has to be specified.')
	if( type != 'gene' & type != 'tx' )
		stop('makeIntronsGtfFromDb: type has to be either "gene" or "tx".')
	if( ! is.character( filename ) )
		stop('makeIntronsGtfFromDb: filename has to be character class.')
	# use method disjointExons to get a GRangeList of disjoint exons 
	message( 'Retrieving data from TxDb...' )
	if( type == 'gene' ) {
		# get a GRangeList of exons belonging to each gene
		exonsDB <- exonsBy( object , 'gene' )
		# reduce each gene ( merge overlapping exons ) and transform to GRanges
		exonsDB <- reduce( exonsDB )
		validRanges <- elementLengths(range(exonsDB))==1
		exonsDB <- exonsDB[validRanges]
		# find 'holes'
		intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
		intronsDB <- unlist(intronsDB)
		# 1253426 is the longest intron given by intronsByTranscripts 
		# function in TxDb, filtering those longer
		ix <- width(intronsDB) <= 1253426
		intronsDB <- intronsDB[ix]
		# feature name to be passed later as annotation
		featureName <- 'gene_id "'
	} else if ( type == 'tx' ) {
		# get a GRangeList of exons belonging to each transcript
		intronsDB <- intronsByTranscript( object, use.names=TRUE )
		# pass from GRangesList to GRanges
		intronsDB <- unlist(intronsDB)
		# feature name to be passed later as annotation
		featureName <- 'tx_name "'
	}
	# make a dataframe from the GRangeList
	intronsDF <- data.frame(
		# chromosome
		chromosome = seqnames( intronsDB )
		# source ( name of the program who generated the feature )
		, source = object$packageName
		# feature 
		, feature = 'exon'
		# start
		, start = start( intronsDB )
		# end
		, end = end( intronsDB )
		# score ( a floating point value )
		, score = '.'
		# strand
		, strand = strand( intronsDB )
		# frame ( One of '0', '1' or '2'. '0' indicates that the first base 
		# of the feature is the first base of a codon, '1' that the second 
		# base is the first base of a codon, and so on.. )
		, frame = '.'
		# attribute ( A semicolon-separated list of tag-value pairs, providing 
		# additional information about each feature. )
		, attribute = paste( featureName , names( intronsDB ) , '"' , sep='' )
	)
	# sorting lines according to genomic position
	message('Sorting according to genomic position...')
	idx <- order( intronsDF$chromosome , intronsDF$start )
	intronsDF <- intronsDF[ idx , ]
	# write the GTF file
	write.table(
		intronsDF 
		, file  = filename 
		, quote = FALSE 
		, sep   = '\t'
		, row.names = FALSE 
		, col.names = FALSE
		)
	message('File created.')
	})