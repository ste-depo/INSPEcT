#' Create a new INSPEcT object
#'
#' @description
#' The function newINSPEcT creates a new instance of the class INSPEcT provided the experimental time points, expression data (like RPKMs) of mature
#' and eventually nascent RNA. For the nascent analysis, it is also requires a collecting time and the scaling factor to
#' normalize the nascent RNA-seq libraries. This latter parameter can also be calculated by the function itself if both exonic and intronic
#' expression data are provided; otherwise it must be given as an input and it is essential to guarantee the robustness of the analysis.
#' @param tpts A vector of time points, one for each sample
#' @param labeling_time A number, lenght of the Nascent pulse
#' @param nascentExpressions A list which contains exons and introns expression matrices and variances for the nascent RNA
#' @param matureExpressions A list which contains exons and introns expression matrices and variances for the mature RNA
#' @param BPPARAM Configuration for BiocParallel parallelization. By default is set to bpparam()
#' @param totalMedianNorm A logical to perform median normalization over total RNA exons rpkms, it will apply also on introns
#' @param labeledMedianNorm A logical to perform median normalization over Nascent RNA exons rpkms, it will apply also on introns
#' @param totalSF A vector storing user defined normalization scale over Total RNA exons and introns rpkms
#' @param labeledSF A vector storing user defined normalization scale over Nascent RNA exons and introns rpkms
#' @param totalQuantileNorm A logical to perform to perform median normalization over total RNA exons rpkms, it will apply also on introns
#' @param labeledQuantileNorm A logical to perform to perform median normalization over Nascent RNA exons rpkms, it will apply also on introns
#' @param simulatedData A logical, set to TRUE in case the analysis is on simulated data
#' @param degDuringPulse A logical, set to TRUE in case of a long labelling time. Also degradation of newly synthesized transcripts will be taken into account
#' @param Dmin A numerical, it is the lower bound of the degradation rate domain for the prior optimization 
#' @param Dmax A numerical, it is the upper bound of the degradation rate domain for the prior optimization
#' @param genesFilter, A logical which, if TRUE, filters out genes which have no signal in at least 2/3 of the time points in each feature
#' @return An object of class INSPEcT with a first estimation of the rates which can be accessed by the method \code{\link{ratesFirstGuess}}.
#' @examples

#' data('allcounts', package='INSPEcT')
#' data('featureWidths', package='INSPEcT')
#' data('libsizes', package='INSPEcT')
#' 
#' matureCounts<-allcounts$mature
#' expDes<-rep(c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16),3)
#' 
#' matExp_DESeq2<-quantifyExpressionsFromTrCounts(
#'  	allcounts=matureCounts
#'  	,libsize=totalLS
#'  	,exonsWidths=exWdths
#'  	,intronsWidths=intWdths
#'  	,experimentalDesign=expDes)
#' 
#' matureInspObj<-newINSPEcT(tpts=tpts
#'                          ,labeling_time=NULL
#'                          ,nascentExpressions=NULL
#'                          ,matureExpressions=matExp_DESeq2)

newINSPEcT <- function(tpts
					 , labeling_time = NULL
					 , nascentExpressions = NULL
					 , matureExpressions
					 , BPPARAM = bpparam()
					 , totalMedianNorm = FALSE
					 , labeledMedianNorm = FALSE
					 , totalSF = NULL
					 , labeledSF = NULL
					 , totalQuantileNorm = FALSE
					 , labeledQuantileNorm = FALSE
					 , simulatedData = FALSE
					 , degDuringPulse = FALSE
					 , Dmin = 1e-6
					 , Dmax = 10
					 , genesFilter = TRUE)
{

	# tpts
	if( !is.numeric(tpts) )
		stop('newINSPEcT: tpts must be a numeric')

	# labeling_time & nascentExpressions
	if( is.null(labeling_time) & is.null(nascentExpressions) ) {

		NoNascent <- TRUE
		message('No nascent RNA mode.')

		# matureExpressions
		if( !is.list(matureExpressions) )
			stop('newINSPEcT: "matureExpressions" must be a list with elements "exonsExpressions", "exonsVariance", "intronsExpressions" and "intronsVariance"')
		if( !all(c("exonsExpressions", "exonsVariance", "intronsExpressions", "intronsVariance") %in% names(matureExpressions)) )
			stop('newINSPEcT: "matureExpressions" must be a list with elements "exonsExpressions", "exonsVariance", "intronsExpressions" and "intronsVariance"')
		rpkms_total_exons <- matureExpressions$exonsExpressions
		rpkms_total_exons_variances <- matureExpressions$exonsVariance
		rpkms_total_introns <- matureExpressions$intronsExpressions
		rpkms_total_introns_variances <- matureExpressions$intronsVariance
		if(
			! ( is.matrix(rpkms_total_exons) & is.matrix(rpkms_total_exons_variances) & 
				is.matrix(rpkms_total_introns) & is.matrix(rpkms_total_introns_variances)) |
			( ncol(rpkms_total_exons) != length(tpts) ) |
			( ncol(rpkms_total_exons_variances) != length(tpts) ) |
			( ncol(rpkms_total_introns) != length(tpts) ) |
			( ncol(rpkms_total_introns_variances) != length(tpts) )
			)
			stop('newINSPEcT: all the elements of "matureExpressions" must be matrices number of columns equal to the number of time points.')
		if( length(setdiff(rownames(rpkms_total_introns), rownames(rpkms_total_exons)))>0 )
			stop('newINSPEcT: some genes in "matureExpressions" have introns quantifications but not the corresponding exons quantifications.')

	} else if( !is.null( labeling_time ) & !is.null(nascentExpressions) ) {

		NoNascent <- FALSE

		if( ! ( is.list(matureExpressions) & is.list(nascentExpressions) ) )
			stop('newINSPEcT: "matureExpressions" and "nascentExpressions" must be of class "list".')

		rpkms_total_exons <- matureExpressions$exonsExpressions
		rpkms_total_exons_variances <- matureExpressions$exonsVariance
		rpkms_total_introns <- matureExpressions$intronsExpressions
		rpkms_total_introns_variances <- matureExpressions$intronsVariance

		rpkms_Nascent_exons <- nascentExpressions$exonsExpressions
		rpkms_Nascent_exons_variances <- nascentExpressions$exonsVariance
		rpkms_Nascent_introns <- nascentExpressions$intronsExpressions
		rpkms_Nascent_introns_variances <- nascentExpressions$intronsVariance

		if (
			! is.null(rpkms_total_exons) &
			! is.null(rpkms_total_exons_variances) &
			! is.null(rpkms_total_introns) &
			! is.null(rpkms_total_introns_variances) &
			! is.null(rpkms_Nascent_exons) &
			! is.null(rpkms_Nascent_exons_variances) &
			! is.null(rpkms_Nascent_introns) &
			! is.null(rpkms_Nascent_introns_variances) &
			  is.matrix(rpkms_total_exons) &
			  is.matrix(rpkms_total_exons_variances) &
			  is.matrix(rpkms_total_introns) &
			  is.matrix(rpkms_total_introns_variances) &
			  is.matrix(rpkms_Nascent_exons) &
			  is.matrix(rpkms_Nascent_exons_variances) &
			  is.matrix(rpkms_Nascent_introns) &
			  is.matrix(rpkms_Nascent_introns_variances) &
			  is.numeric(rpkms_total_exons) &
			  is.numeric(rpkms_total_exons_variances) &
			  is.numeric(rpkms_total_introns) &
			  is.numeric(rpkms_total_introns_variances) &
			  is.numeric(rpkms_Nascent_exons) &
			  is.numeric(rpkms_Nascent_exons_variances) &
			  is.numeric(rpkms_Nascent_introns) &
			  is.numeric(rpkms_Nascent_introns_variances)

			) {

			NascentSimple <- FALSE
			message('Nascent RNA mode, with exon and intron quantifications.')

			if(	!(
				ncol(rpkms_total_exons) == length(tpts) &
				ncol(rpkms_total_exons_variances) == length(tpts) &
				ncol(rpkms_total_introns) == length(tpts) &
				ncol(rpkms_total_introns_variances) == length(tpts) &
				ncol(rpkms_Nascent_exons) == length(tpts) &
				ncol(rpkms_Nascent_exons_variances) == length(tpts) &
				ncol(rpkms_Nascent_introns) == length(tpts) &
				ncol(rpkms_Nascent_introns_variances) == length(tpts)
				) )
				stop('All elements of "matureExpressions" and "nascentExpressions" must be matrices with number of columns equal to the number of time points.')

			## check that matrices have rownames
			if( 
				is.null(rownames(rpkms_total_exons)) 
				| is.null(rownames(rpkms_total_exons_variances))
				| is.null(rownames(rpkms_Nascent_exons))
				| is.null(rownames(rpkms_Nascent_exons_variances)) 
				| is.null(rownames(rpkms_total_introns)) 
				| is.null(rownames(rpkms_total_introns_variances))
				| is.null(rownames(rpkms_Nascent_introns))
				| is.null(rownames(rpkms_Nascent_introns_variances)) 
				)
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" must be matrices with rownames.')
			## check that the rownames correspond
			if( 
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_exons_variances)) |
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons)) |
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons_variances))
				)
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" relative to exons must be matrices with idenitcal rownames.')
			if( 
				!identical(rownames(rpkms_Nascent_introns), rownames(rpkms_Nascent_introns_variances)) |
				!identical(rownames(rpkms_Nascent_introns), rownames(rpkms_total_introns)) |
				!identical(rownames(rpkms_Nascent_introns), rownames(rpkms_total_introns_variances))
				)
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" relative to introns must be matrices with idenitcal rownames.')

			if( length(setdiff(rownames(rpkms_Nascent_introns), rownames(rpkms_Nascent_exons)))>0 )
				stop('newINSPEcT: some genes in "nascentExpressions" have introns quantifications but not the corresponding exons quantifications.')

		} else if (

			! is.null(rpkms_total_exons) &
			! is.null(rpkms_total_exons_variances) &
			  is.null(rpkms_total_introns) &
			  is.null(rpkms_total_introns_variances) &
			! is.null(rpkms_Nascent_exons) &
			! is.null(rpkms_Nascent_exons_variances) &
			  is.null(rpkms_Nascent_introns) &
			  is.null(rpkms_Nascent_introns_variances) &
			  is.matrix(rpkms_total_exons) &
			  is.matrix(rpkms_total_exons_variances) &
			  is.matrix(rpkms_Nascent_exons) &
			  is.matrix(rpkms_Nascent_exons_variances) &
			  is.numeric(rpkms_total_exons) &
			  is.numeric(rpkms_total_exons_variances) &
			  is.numeric(rpkms_Nascent_exons) &
			  is.numeric(rpkms_Nascent_exons_variances)
			) {

			NascentSimple <- TRUE

			message('Nascent RNA mode, no intron quantifications.')

			if(	!(
				ncol(rpkms_total_exons) == length(tpts) &
				ncol(rpkms_total_exons_variances) == length(tpts) &
				ncol(rpkms_Nascent_exons) == length(tpts) &
				ncol(rpkms_Nascent_exons_variances) == length(tpts)
				) )
				stop('All elements of "matureExpressions" and "nascentExpressions" must be matrices with number of columns equal to the number of time points.')

			## in the simple mode the scaling factor between the two library is
			# strongly reccomended
			if( is.null(labeledSF) )
				stop('newINSPEcT: Without Intronic quantifications, labeledSF must be provided.')
			## check that matrices have rownames
			if( is.null(rownames(rpkms_total_exons)) 
				| is.null(rownames(rpkms_total_exons_variances))
				| is.null(rownames(rpkms_Nascent_exons))
				| is.null(rownames(rpkms_Nascent_exons_variances)) )
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" must be matrices with idenitcal rownames.')
			## check that the rownames correspond
			if( 
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_exons_variances)) |
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons)) |
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons_variances))
				)
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" must be matrices with idenitcal rownames.')

		} else {

			stop('newINSPEcT: "matureExpressions" and "nascentExpressions" must be lists with elements named "exonsExpressions", "exonsVariance", 
				"intronsExpressions" and "intronsVariance" in case of a full analysis, or with elements named "exonsExpressions", "exonsVariance"
				in case of partial analysis. All elements must be numeric matrices with number of columns equal to the number of time points.')

		}

		if( !is.logical(labeledMedianNorm) )
			stop('newINSPEcT: labeledMedianNorm must be a logical.')
		if( !is.logical(labeledQuantileNorm) )
			stop('newINSPEcT: labeledQuantileNorm must be a logical.')
		if( !is.logical(simulatedData) )
			stop('newINSPEcT: simulatedData must be a logical.')
		if( !is.logical(degDuringPulse) )
			stop('newINSPEcT: degDuringPulse must be a logical.')


	} else {

		stop('newINSPEcT: labeling_time and nascentExpressions must be both NULL or defined')
	
	}

	if( !class(BPPARAM) %in% sapply(registered(),class) )
		stop('newINSPEcT: BPPARAM argument not registered.')
	if( !is.logical(totalMedianNorm) )
		stop('newINSPEcT: totalMedianNorm must be a logical.')
	if( !is.logical(totalQuantileNorm) )
		stop('newINSPEcT: totalQuantileNorm must be a logical.')

	if( NoNascent )
	{
		#Only genes with exons and introns can be modeled without nascent RNA 
		eiGenes <- intersect(rownames(rpkms_total_exons),rownames(rpkms_total_introns))

		rpkms_total_exons <- rpkms_total_exons[eiGenes,,drop=FALSE]
		rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes,,drop=FALSE]
		rpkms_total_introns <- rpkms_total_introns[eiGenes,,drop=FALSE]
		rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes,,drop=FALSE]

		negativeMature <- apply(rpkms_total_exons<rpkms_total_introns,1,any)
		if( any(negativeMature) ) {
			message('Removing genes with intronic quantifications greater than the exonic..')
			eiGenes <- eiGenes[!negativeMature]
			rpkms_total_exons <- rpkms_total_exons[eiGenes, ,drop=FALSE]
			rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes, ,drop=FALSE]
			rpkms_total_introns <- rpkms_total_introns[eiGenes, ,drop=FALSE]
			rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes, ,drop=FALSE]
		}

		if(genesFilter)
		{
			########
			### filter out genes which have no signal in at least 2/3 of the time points in each feature
			#######
			ix1 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
			ix2 <- apply(rpkms_total_exons_variances, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons_variances)
			ix3 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
			ix4 <- apply(rpkms_total_introns_variances, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns_variances)

			filteroutGenes <- rownames(rpkms_total_exons)[ix1>2/3 | ix2>2/3 | ix3>2/3 | ix4>2/3 ]

			if( length(filteroutGenes)>0 ) {
				message('Filtering out genes with more than 2/3 of zeros in their exonic or intronic quantifications..')
				eiGenes <- eiGenes[!eiGenes %in% filteroutGenes]
				rpkms_total_exons <- rpkms_total_exons[eiGenes, ,drop=FALSE]
				rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[eiGenes, ,drop=FALSE]
				rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes, ,drop=FALSE]
			}
		}
		
		message(paste('Number of genes with introns and exons: ', length(eiGenes)))

		## subset the "introns and exons genes" from the total
		totRpkmsIntEx <- list(
				exons=rpkms_total_exons
				, exons_var=rpkms_total_exons_variances
				, introns=rpkms_total_introns
				, introns_var=rpkms_total_introns_variances)

		tpts <- sort(unique(tpts))

		## estimate the rates
		out <- ratesAndConcentrationsNoNascent(totRpkms = totRpkmsIntEx
											   , tpts = tpts
											   , BPPARAM = BPPARAM
											   , modellingParameters = list(Dmin = Dmin, Dmax = Dmax)
											   , genesFilter = genesFilter)

		return(createInspectObject(out, NoNascent, FALSE, character(0)))
	
	} else { #Nascent RNA mode

		if( simulatedData ) {

			intExGenes <- rownames(rpkms_Nascent_exons)
			onlyExGenes <- character(0)
			rpkms_Nascent_introns <- matrix(0, nrow(rpkms_Nascent_exons), ncol(rpkms_Nascent_exons))
			rownames(rpkms_Nascent_introns) <- rownames(rpkms_Nascent_exons)
	
		} else {

			## in case either introns of Nascent or introns of total fraction are not
			# provided (both are required) proceed in the simple mode
			if( NascentSimple )
			{
	
				if( genesFilter ) {

					########
					### filter out genes which have no signal in at least 2/3 of the time points in each feature
					#######
		
					ix1 <- apply(rpkms_Nascent_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_exons)
					ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
					filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1>2/3 | ix2>2/3]
		
					if( length(filteroutGenes)>0 ) {
						message(paste('Filtering out', length(filteroutGenes), 'genes with more than 2/3 of zeros in their exonic quantifications.'))
						rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
						rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
					}

				}
	
				## assign all the genes to the "only exons" mode.
				intExGenes <- character(0)
				onlyExGenes <- rownames(rpkms_Nascent_exons)
	
			} else {

				## in case, apply quantile normalization
				if( totalQuantileNorm ) {
					normMat <- normalize.quantiles(rbind(rpkms_total_exons, rpkms_total_introns))
					exonFeatures <- rownames(rpkms_total_exons)
					rpkms_total_exons <- normMat[1:nrow(rpkms_total_exons),]
					rownames(rpkms_total_exons) <- exonFeatures
					intronFeatures <- rownames(rpkms_total_introns)
					rpkms_total_introns <- normMat[(nrow(rpkms_total_exons)+1):nrow(normMat),]
					rownames(rpkms_total_introns) <- intronFeatures
				}
				if( labeledQuantileNorm ) {
					normMat <- normalize.quantiles(rbind(rpkms_Nascent_exons, rpkms_Nascent_introns))
					exonFeatures <- rownames(rpkms_Nascent_exons)
					rpkms_Nascent_exons <- normMat[1:nrow(rpkms_Nascent_exons),]
					rownames(rpkms_Nascent_exons) <- exonFeatures
					intronFeatures <- rownames(rpkms_Nascent_introns)
					rpkms_Nascent_introns <- normMat[(nrow(rpkms_Nascent_exons)+1):nrow(normMat),]
					rownames(rpkms_Nascent_introns) <- intronFeatures
				}

				if( genesFilter ) {

					########
					### filter out genes which have no signal in at least 2/3 of the time points in each feature
					#######
		
					ix1 <- apply(rpkms_Nascent_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_exons)
					ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
					filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1>2/3 | ix2>2/3] # | ix3>2/3 | ix4>2/3]
					if( length(filteroutGenes)>0 ) {
						message(paste('Filtering out', length(filteroutGenes), 'genes with more than 2/3 of zeros in their exonic quantifications.'))
						rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
						rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
						rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
						rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
					}
		
					ix3 <- apply(rpkms_Nascent_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_introns)
					ix4 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
					filteroutGenes <- rownames(rpkms_Nascent_exons)[ix3>2/3 | ix4>2/3] # | ix3>2/3 | ix4>2/3]
					if( length(filteroutGenes)>0 ) {
						message(paste('Filtering out intronic signal of', length(filteroutGenes), 'genes with more than 2/3 of zero quantifications.'))
						message('(for those genes only synthesis and degradation will be evaluated)')
						rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
						rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
					}

				}
	
				## assign genes to "only exons" or "introns and exons" wheter they have introns or not
				if( identical(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_introns)) ) {
					intExGenes <- rownames(rpkms_Nascent_exons)
					onlyExGenes <- character(0)
				} else {
					## exons rpkms must be provided for every genes, if not: stop
					intExGenes <- intersect(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_introns))
					onlyExGenes <- setdiff(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_introns))
				}
			}
		
		}

		totRpkms <- list()
		labeledRpkms <- list()
	
		tpts <- sort(unique(tpts))

		if( length(intExGenes)>0 ) {
			message(paste('Number of genes with introns and exons: ', length(intExGenes)))
			## subset the "introns and exons genes" from the total
			totRpkmsIntEx <- list(
					exons=rpkms_total_exons[intExGenes, , drop=FALSE]
					, exons_var=rpkms_total_exons_variances[intExGenes, , drop=FALSE]
					, introns=rpkms_total_introns[intExGenes, , drop=FALSE]
					, introns_var=rpkms_total_introns_variances[intExGenes, , drop=FALSE]
			)
			labeledRpkmsIntEx <- list(
					exons=rpkms_Nascent_exons[intExGenes, , drop=FALSE]
					, exons_var=rpkms_Nascent_exons_variances[intExGenes, , drop=FALSE]
					, introns=rpkms_Nascent_introns[intExGenes, , drop=FALSE]
					, introns_var=rpkms_Nascent_introns_variances[intExGenes, , drop=FALSE]
				)
			## estimate the rates
			outIntEx <- .getRatesAndConcentrationsFromRpkms(totRpkms=totRpkmsIntEx
														  , labeledRpkms=labeledRpkmsIntEx
														  , tpts=tpts
														  , tL=labeling_time
														  , simulatedData=simulatedData
														  , BPPARAM=BPPARAM
														  , totalMedianNorm=totalMedianNorm
														  , labeledMedianNorm=labeledMedianNorm
														  , totalSF=totalSF
														  , labeledSF=labeledSF
														  , degDuringPulse=degDuringPulse)
			## set the labeledSF and totalSF, they can be used by "only exons genes" (if present)
			labeledSF <- outIntEx$labeledSF
			totalSF <- outIntEx$totalSF
			totalMedianNorm <- FALSE
			labeledMedianNorm <- FALSE
		}

		if( length(onlyExGenes)>0 ) {
			message(paste('Number of genes with only exons: ', length(onlyExGenes)))
			## subset the "only exons genes" from the total
			totRpkmsOnlyEx <- list(
					exons=rpkms_total_exons[onlyExGenes, , drop=FALSE]
					, exons_var=rpkms_total_exons_variances[onlyExGenes, , drop=FALSE]
				)
			labeledRpkmsOnlyEx <- list(
					exons=rpkms_Nascent_exons[onlyExGenes, , drop=FALSE]
					, exons_var=rpkms_Nascent_exons_variances[onlyExGenes, , drop=FALSE]
				)
			## estimate the rates, eventually using "labeledSF" and "totalSF" calculated from "introns and exons genes"

			outOnlyEx <- .getRatesAndConcentrationsFromRpkms(totRpkmsOnlyEx
														   , labeledRpkmsOnlyEx
														   , tpts
				, tL=labeling_time, simulatedData=simulatedData, BPPARAM=BPPARAM
				, totalMedianNorm=totalMedianNorm, labeledMedianNorm=labeledMedianNorm
				, totalSF=totalSF, labeledSF=labeledSF, degDuringPulse=degDuringPulse)
		}

		## merge the output of "introns and exons genes" and "only exons genes" (if both present)
		if( length(intExGenes)>0 & length(onlyExGenes)>0 ) {
			## internal control, these data must be identical in the two cathegories
			## if not: something unpredicted happened
			stopifnot(identical(outIntEx$tpts, outOnlyEx$tpts))
			stopifnot(identical(outIntEx$labeledSF, outOnlyEx$labeledSF))
			stopifnot(identical(outIntEx$totalSF, outOnlyEx$totalSF))
			stopifnot(identical(outIntEx$tL, outOnlyEx$tL))
			## merge the two lists

			out <- list(
				concentrations=list(
					total=rbind(outIntEx$concentrations$total, outOnlyEx$concentrations$total)
					#, total_var=c(outIntEx$concentrations$total_var, outOnlyEx$concentrations$total_var)
					, total_var=rbind(outIntEx$concentrations$total_var, outOnlyEx$concentrations$total_var)
					, preMRNA=rbind(outIntEx$concentrations$preMRNA, outOnlyEx$concentrations$preMRNA)
					#, preMRNA_var=c(outIntEx$concentrations$preMRNA_var, outOnlyEx$concentrations$preMRNA_var)
					, preMRNA_var=rbind(outIntEx$concentrations$preMRNA_var, outOnlyEx$concentrations$preMRNA_var)
					)
				, rates=list(
					alpha=rbind(outIntEx$rates$alpha, outOnlyEx$rates$alpha)
					, alpha_var=rbind(outIntEx$rates$alpha_var, outOnlyEx$rates$alpha_var)
					, beta=rbind(outIntEx$rates$beta, outOnlyEx$rates$beta)
					, gamma=rbind(outIntEx$rates$gamma, outOnlyEx$rates$gamma)
					)
				, ratesEstimPrec=rbind(outIntEx$ratesEstimPrec, outOnlyEx$ratesEstimPrec)
				, geneNames=c(outIntEx$geneNames, outOnlyEx$geneNames)
				## from here both lists have the same information, pick it up from the first one
				, tpts=outIntEx$tpts
				, labeledSF=outIntEx$labeledSF
				, totalSF=outIntEx$totalSF
				, tL=outIntEx$tL
				)
		} else if( length(intExGenes)>0 ) {
			out <- outIntEx
		} else {
			out <- outOnlyEx
		}
	
		return(createInspectObject(out, NoNascent, NascentSimple, onlyExGenes))

	}
}



createInspectObject  <- function(out, NoNascent, NascentSimple, onlyExGenes) {

	# make an "ExpressionSet" object containing all the information
	nTpts <- length(out$tpts)
	exprData <- cbind(out$concentrations$total
					, out$concentrations$preMRNA
					, out$rates$alpha
					, out$rates$beta
					, out$rates$gamma)

	pData <- data.frame(feature=c(rep('total',nTpts)
								, rep('preMRNA',nTpts)
								, rep('synthesis',nTpts)
								, rep('degradation',nTpts)
								, rep('processing',nTpts) )
					  , time=rep(out$tpts, 5))

	colnames(exprData) <- paste(pData$feature, signif(pData$time,2), sep='_')
	rownames(pData) <- colnames(exprData)
	phenoData <- new('AnnotatedDataFrame', data=pData)
	
	fData <- data.frame(total=out$concentrations$total_var
					  , preMRNA=out$concentrations$preMRNA_var
					  , synthesis=out$rates$alpha_var
					  , degradation=1
					  , processing=1
					  )
	rownames(fData) <- out$geneNames

	if( (! NoNascent) & length(onlyExGenes)>0 ) fData[onlyExGenes,]$processing <- NA

	featureData <- new('AnnotatedDataFrame', data=fData)
	ratesFirstGuess <- ExpressionSet(assayData=exprData
								   , phenoData=phenoData
								   , featureData=featureData)

	# Controls necessary to keep the output of newINSPEcT_NoNascent and newINSPEcT equal

	if( is.null(out$totalSF) ) out$totalSF <- numeric(0)
	if( is.null(out$labeledSF) ) out$labeledSF <- numeric(0)
	if( is.null(out$tL) ) out$tL <- numeric(0)

	## update the object and return it
	object <- new('INSPEcT')
	object@tpts <- out$tpts
	object@totalSF <- out$totalSF
	object@labeledSF <- out$labeledSF
	object@tL <- out$tL	
	object@ratesFirstGuess <- ratesFirstGuess
	object@precision <- out$ratesEstimPrec
	object@model@simple <- TRUE
	
	if( NoNascent ) {
		object@ratesFirstGuessP <- out$ratesFirstGuessP
	} else {
		if( NascentSimple )
			object@model@simple <- TRUE
	}

	return(object)

}











