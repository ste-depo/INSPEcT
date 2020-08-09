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
#' @param preexisting A logical, indicating if the mature expression refers to the pre-exising (unlabeled) population. Not implemented yet for the "degDuringPulse" mode.
#' @param BPPARAM Configuration for BiocParallel parallelization. By default is set to SerialParam()
#' @param labeledSF A vector storing user defined normalization scale over Nascent RNA exons and introns quantifications
#' @param simulatedData A logical, set to TRUE in case the analysis is on simulated data
#' @param degDuringPulse A logical, set to TRUE in case of a long labelling time. Also degradation of newly synthesized transcripts will be taken into account
#' @param Dmin A numerical, it is the lower bound of the degradation rate domain for the prior optimization 
#' @param Dmax A numerical, it is the upper bound of the degradation rate domain for the prior optimization
#' @param genesFilter A logical, if TRUE, filters out genes which have no signal in at least a given fraction (2/3 by default) of the observations
#' @param genesFilterThreshold A number, threshold to use for genes filtering (2/3 by default)
#' @param imputeNAs A logical, if TRUE the rates first guess which are not finite are imputed from the neighbours.
#' @return An object of class INSPEcT with a first estimation of the rates which can be accessed by the method \code{\link{ratesFirstGuess}}
#' @examples
#' data('allcounts', package='INSPEcT')
#' data('featureWidths', package='INSPEcT')
#' data('libsizes', package='INSPEcT')
#' 
#' matureCounts<-allcounts$mature
#' tpts <- c(0,1/6,1/3,1/2,1,1.5,2,4,8,12,16)
#' expDes<-rep(tpts,3)
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
					, preexisting = FALSE
					, BPPARAM = SerialParam()
					, labeledSF = NULL
					, simulatedData = FALSE
					, degDuringPulse = FALSE
					, Dmin = 1e-6
					, Dmax = 10
					, genesFilter = TRUE
					, genesFilterThreshold = 2/3
					, imputeNAs = TRUE)
{

	##################################
	### check input arguments:

	# check tpts
	if( !is.vector(tpts) )
		stop('newINSPEcT: tpts must be a vector')

	# check labeling_time & nascentExpressions
	if( is.null(nascentExpressions) & is.null(labeling_time) ) {

		# in case nascent expression is not provided (and coherently also labelling time is not)
		# set the mode to "No nascent RNA"
		NoNascent <- TRUE
		message('No nascent RNA mode.')
		# check matureExpressions format
		if( !is.list(matureExpressions) )
			stop('newINSPEcT: "matureExpressions" must be a list with elements "exonsExpressions", "exonsVariance", "intronsExpressions" and "intronsVariance"')
		if( !all(c("exonsExpressions", "exonsVariance", "intronsExpressions", "intronsVariance") %in% names(matureExpressions)) )
			stop('newINSPEcT: "matureExpressions" must be a list with elements "exonsExpressions", "exonsVariance", "intronsExpressions" and "intronsVariance"')
		# assign elements of the list "matureExpressions" to individual variables
		# and check that they are matrices with the correct dimensions
		rpkms_total_exons <- matureExpressions$exonsExpressions
		rpkms_total_exons_variances <- matureExpressions$exonsVariance
		rpkms_total_introns <- matureExpressions$intronsExpressions
		rpkms_total_introns_variances <- matureExpressions$intronsVariance
		if(
			! ( is.matrix(rpkms_total_exons) & is.matrix(rpkms_total_exons_variances) & 
				is.matrix(rpkms_total_introns) & is.matrix(rpkms_total_introns_variances)) |
			! ( is.numeric(rpkms_total_exons) & is.numeric(rpkms_total_exons_variances) & 
				is.numeric(rpkms_total_introns) & is.numeric(rpkms_total_introns_variances)) |
			( ncol(rpkms_total_exons) != length(tpts) ) |
			( ncol(rpkms_total_exons_variances) != length(tpts) ) |
			( nrow(rpkms_total_exons) != nrow(rpkms_total_exons_variances) ) |
			( ncol(rpkms_total_introns) != length(tpts) ) |
			( ncol(rpkms_total_introns_variances) != length(tpts) ) |
			( nrow(rpkms_total_introns) != nrow(rpkms_total_introns_variances) )
			)
			stop('newINSPEcT: all the elements of "matureExpressions" must be matrices of numerics with number of columns equal to the number of time points and number of rows matching between expressions and variances.')
		## check that matrices have rownames
		if( 
			is.null(rownames(rpkms_total_exons)) 
			| is.null(rownames(rpkms_total_exons_variances))
			| is.null(rownames(rpkms_total_introns)) 
			| is.null(rownames(rpkms_total_introns_variances))
			)
			stop('newINSPEcT: the elements of "matureExpressions" must be matrices with rownames.')
		## check that the rownames correspond between the exons matrices
		if( !identical(rownames(rpkms_total_exons), rownames(rpkms_total_exons_variances)) )
			stop('newINSPEcT: the elements of "matureExpressions" relative to exons must be matrices with idenitcal rownames.')
		## check that the rownames correspond between the introns matrices
		if( !identical(rownames(rpkms_total_introns), rownames(rpkms_total_introns_variances)) )
			stop('newINSPEcT: the elements of "matureExpressions" relative to introns must be matrices with idenitcal rownames.')
		# check Dmin and Dmax arguments
		if( !is.numeric(Dmin) ) stop('newINSPEcT: Dmin must be a numeric')
		if( !is.numeric(Dmax) ) stop('newINSPEcT: Dmax must be a numeric')

	} else if( !is.null(nascentExpressions) & !is.null( labeling_time ) ) {

		# in case both nascent expression and labelling time are provided
		# set the mode to "Nascent RNA"
		NoNascent <- FALSE
		# check matureExpressions and nascentExpressions format
		if( ! ( is.list(matureExpressions) & is.list(nascentExpressions) ) )
			stop('newINSPEcT: "matureExpressions" and "nascentExpressions" must be of class "list".')
		# assign elements of the list "matureExpressions" to individual variables
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
			# In case total and nascent are defined in both exons and introns
			# and they are matrices composed by numerics, set the mode 
			# "Nascent RNA, full"
			NascentSimple <- FALSE
			message('Nascent RNA mode, with exon and intron quantifications.')
			# Then check the dimensions of the matrices
			if(	!(
				ncol(rpkms_total_exons) == length(tpts) &
				ncol(rpkms_total_exons_variances) == length(tpts) &
				nrow(rpkms_total_exons) == nrow(rpkms_total_exons_variances) &
				ncol(rpkms_total_introns) == length(tpts) &
				ncol(rpkms_total_introns_variances) == length(tpts) &
				nrow(rpkms_total_introns) == nrow(rpkms_total_introns_variances) &
				ncol(rpkms_Nascent_exons) == length(tpts) &
				ncol(rpkms_Nascent_exons_variances) == length(tpts) &
				nrow(rpkms_Nascent_exons) == nrow(rpkms_Nascent_exons_variances) &
				ncol(rpkms_Nascent_introns) == length(tpts) &
				ncol(rpkms_Nascent_introns_variances) == length(tpts) &
				nrow(rpkms_Nascent_introns) == nrow(rpkms_Nascent_introns_variances)
				) )
				stop('All elements of "matureExpressions" and "nascentExpressions" must be matrices with number of columns equal to the number of time points and number of rows matching between expressions and variances.')
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
			## check that the rownames correspond between the exons matrices
			if( 
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_Nascent_exons_variances)) |
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons)) |
				!identical(rownames(rpkms_Nascent_exons), rownames(rpkms_total_exons_variances))
				)
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" relative to exons must be matrices with idenitcal rownames.')
			## check that the rownames correspond between the introns matrices
			if( 
				!identical(rownames(rpkms_Nascent_introns), rownames(rpkms_Nascent_introns_variances)) |
				!identical(rownames(rpkms_Nascent_introns), rownames(rpkms_total_introns)) |
				!identical(rownames(rpkms_Nascent_introns), rownames(rpkms_total_introns_variances))
				)
				stop('newINSPEcT: the elements of "matureExpressions" and "nascentExpressions" relative to introns must be matrices with idenitcal rownames.')
			## check that genes with introns defined have also the corresponding exons (but not vice-versa)
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
			# In case total and nascent are defined in only in their exons
			# and they are matrices composed by numerics, set the mode 
			# "Nascent RNA, simple"
			NascentSimple <- TRUE
			message('Nascent RNA mode, no intron quantifications.')
			# Then check the dimensions of the matrices
			if(	!(
				ncol(rpkms_total_exons) == length(tpts) &
				ncol(rpkms_total_exons_variances) == length(tpts) &
				nrow(rpkms_total_exons) == nrow(rpkms_total_exons_variances) &
				ncol(rpkms_Nascent_exons) == length(tpts) &
				ncol(rpkms_Nascent_exons_variances) == length(tpts) &
				nrow(rpkms_Nascent_exons) == nrow(rpkms_Nascent_exons_variances)
				) )
				stop('newINSPEcT: All elements of "matureExpressions" and "nascentExpressions" must be matrices with number of columns equal to the number of time points and number of rows matching between expressions and variances.')
			## in the simple mode the scaling factor between the two library is mandatory
			if( is.null(labeledSF) )
				stop('newINSPEcT: Without Intronic quantifications, labeledSF must be provided.')
			if( length(labeledSF) != length(tpts) )
				stop('newINSPEcT: The length of labeledSF must equal the number of time points.')
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

			stop('newINSPEcT: "matureExpressions" and "nascentExpressions" must be lists with elements named "exonsExpressions", "exonsVariance", "intronsExpressions" and "intronsVariance" in case of a full analysis, or with elements named "exonsExpressions", "exonsVariance" in case of partial analysis. All elements must be numeric matrices with number of columns equal to the number of time points and number of rows matching between expressions and variances.')

		}

		# check simulatedData
		if( !is.logical(simulatedData) )
			stop('newINSPEcT: simulatedData must be a logical.')
		# check degDuringPulse
		if( !is.logical(degDuringPulse) )
			stop('newINSPEcT: degDuringPulse must be a logical.')

	} else {

		# in case nascent expression is provided and not labelling time (or viceversa)
		# throw an error
		stop('newINSPEcT: labeling_time and nascentExpressions must be both NULL (in case of analysis without nascent) or defined (in case of analysis with nascent).')
	
	}

	# check BPPARAM
	if( !class(BPPARAM) %in% sapply(registered(),class) )
		stop('newINSPEcT: BPPARAM argument not registered.')
	# check genesFilter
	if( !is.logical(genesFilter) )
		stop('newINSPEcT: genesFilter must be a logical.')

	### check input arguments, end.
	##################################

	######################################################################################################
	# start the quantification of rates, according to the analysis mode selected based on the input. #######
	######################################################################################################

	##############################
	###### No nascent mode:

	if( NoNascent )
	{
		# Only genes with exons and introns can be modeled without nascent RNA 
		eiGenes <- intersect(rownames(rpkms_total_exons),rownames(rpkms_total_introns))
		if( !identical(eiGenes, rownames(rpkms_total_exons)) ) {
			message(paste('Removing', length(setdiff(rownames(rpkms_total_exons),eiGenes)), 'gene(s) without intronic quantifications.'))
			rpkms_total_exons <- rpkms_total_exons[eiGenes,,drop=FALSE]
			rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes,,drop=FALSE]
			rpkms_total_introns <- rpkms_total_introns[eiGenes,,drop=FALSE]
			rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes,,drop=FALSE]
		}
		# Remove genes with intron quantification greater than the exon quantification, because they give rise to negative mature abundances.
		negativeMature <- apply(rpkms_total_exons<rpkms_total_introns,1,any)
		if( any(negativeMature) ) {
			message(paste('Removing', length(which(negativeMature)), 'gene(s) with intronic quantifications greater than the exonic.'))
			eiGenes <- eiGenes[!negativeMature]
			rpkms_total_exons <- rpkms_total_exons[eiGenes, ,drop=FALSE]
			rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes, ,drop=FALSE]
			rpkms_total_introns <- rpkms_total_introns[eiGenes, ,drop=FALSE]
			rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes, ,drop=FALSE]
		}

		# Filter genes according to expression levels (in case the flag is active)
		if( genesFilter )
		{
			### filter out genes which have no signal in at least genesFilterThreshold of the time points in each feature
			ix1 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
			ix2 <- apply(rpkms_total_exons_variances, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons_variances)
			ix3 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
			ix4 <- apply(rpkms_total_introns_variances, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns_variances)
			filteroutGenes <- rownames(rpkms_total_exons)[ix1>genesFilterThreshold | ix2>genesFilterThreshold | ix3>genesFilterThreshold | ix4>genesFilterThreshold ]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out', length(filteroutGenes),'gene(s) with more than ',round(100*genesFilterThreshold),'% of zeros in their exonic or intronic quantifications..'))
				eiGenes <- eiGenes[!eiGenes %in% filteroutGenes]
				rpkms_total_exons <- rpkms_total_exons[eiGenes, ,drop=FALSE]
				rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[eiGenes, ,drop=FALSE]
				rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes, ,drop=FALSE]
			}
		} else {
			### in any case, filter out genes which have no signal in all observations
			ix1 <- apply(rpkms_total_exons==0, 1, all)
			ix2 <- apply(rpkms_total_exons_variances==0, 1, all)
			ix3 <- apply(rpkms_total_introns==0, 1, all)
			ix4 <- apply(rpkms_total_introns_variances==0, 1, all)
			filteroutGenes <- rownames(rpkms_total_exons)[ix1 | ix2 | ix3 | ix4 ]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out', length(filteroutGenes),'gene(s) with more all zeros in their exonic or intronic quantifications..'))
				eiGenes <- eiGenes[!eiGenes %in% filteroutGenes]
				rpkms_total_exons <- rpkms_total_exons[eiGenes, ,drop=FALSE]
				rpkms_total_exons_variances <- rpkms_total_exons_variances[eiGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[eiGenes, ,drop=FALSE]
				rpkms_total_introns_variances <- rpkms_total_introns_variances[eiGenes, ,drop=FALSE]
			}		
		}
		message(paste('Number of genes with introns and exons: ', length(eiGenes)))
		## sort according to time point ordering
		if( is.numeric(tpts) ) {
			ord <- order(tpts)
			tpts <- tpts[ord]
		} else {
			ord <- seq_along(tpts)
		}
		## estimate RNA dynamics
		totRpkmsIntEx <- list(
				exons=rpkms_total_exons[,ord,drop=FALSE]
				, exons_var=rpkms_total_exons_variances[,ord,drop=FALSE]
				, introns=rpkms_total_introns[,ord,drop=FALSE]
				, introns_var=rpkms_total_introns_variances[,ord,drop=FALSE])
		if( is.numeric(tpts) ) {
			out <- RNAdynamics_NoNascent(totRpkms = totRpkmsIntEx
												   , tpts = tpts
												   , BPPARAM = BPPARAM
												   , modellingParameters = list(Dmin = Dmin, Dmax = Dmax)
												   , genesFilter = genesFilter
												   , imputeNAs = imputeNAs)
			return(createInspectObject(out, NoNascent=TRUE))
		} else {
			####### generate the INSPEcT_steadyNoNascent-class #########
			object <- new('INSPEcT_steadyNoNascent')
			object@sampleNames <- tpts
			object@geneNames <- rownames(rpkms_total_introns)
			object@premature <- rpkms_total_introns
			object@prematureVar <- rpkms_total_introns_variances
			object@mature <- rpkms_total_exons - rpkms_total_introns
			object@matureVar <- rpkms_total_exons_variances + rpkms_total_introns_variances
			return(object)
		}
	} 

	###### Simulated data mode (do not provide nascent introns to the compute of RNA dynamics): 

	if( simulatedData ) {

		# to remove??
		# intExGenes <- rownames(rpkms_Nascent_exons)
		# onlyExGenes <- character(0)
		# rpkms_Nascent_introns <- matrix(0, nrow(rpkms_Nascent_exons), ncol(rpkms_Nascent_exons))
		# rownames(rpkms_Nascent_introns) <- rownames(rpkms_Nascent_exons)
		# end - to remove

		## sort according to time point ordering
		if( is.numeric(tpts) ) {
			ord <- order(tpts)
			tpts <- tpts[ord]
		} else {
			ord <- seq_along(tpts)
		}
		## estimate RNA dynamics
		totRpkmsIntEx <- list(
				exons=rpkms_total_exons[,ord,drop=FALSE]
				, exons_var=rpkms_total_exons_variances[,ord,drop=FALSE]
				, introns=rpkms_total_introns[,ord,drop=FALSE]
				, introns_var=rpkms_total_introns_variances[,ord,drop=FALSE]
				)
		labeledRpkmsIntEx <- list(
				exons=rpkms_Nascent_exons[,ord,drop=FALSE]
				, exons_var=rpkms_Nascent_exons_variances[,ord,drop=FALSE]
			)
		outSim <- RNAdynamics(totRpkms=totRpkmsIntEx
							, labeledRpkms=labeledRpkmsIntEx
							, tpts=tpts
							, tL=labeling_time
							, simulatedData=TRUE
							, BPPARAM=BPPARAM)
		## return the results in the form of an INSPEcT object.
		return(createInspectObject(outSim))

	}

	###### Nascent mode, simple: 

	if( NascentSimple )
	{
		# Filter genes according to expression levels (in case the flag is active)
		if( genesFilter ) {
			### filter out genes which have no signal in at least genesFilterThreshold of the time points in each feature
			ix1 <- apply(rpkms_Nascent_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_exons)
			ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
			filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1>genesFilterThreshold | ix2>genesFilterThreshold]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out', length(filteroutGenes), 'genes with more than ',round(100*genesFilterThreshold),'% of zeros in their exonic quantifications.'))
				rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_Nascent_exons_variances <- rpkms_Nascent_exons_variances[!rownames(rpkms_Nascent_exons_variances) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons_variances <- rpkms_total_exons_variances[!rownames(rpkms_total_exons_variances) %in% filteroutGenes, ,drop=FALSE]
			}
		} else {
			### filter out genes which have no signal in all oservations of a feature
			ix1 <- apply(rpkms_Nascent_exons, 1, all)
			ix2 <- apply(rpkms_total_exons, 1, all)
			filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1 | ix2]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out', length(filteroutGenes), 'genes with more than all zeros in their exonic quantifications.'))
				rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_Nascent_exons_variances <- rpkms_Nascent_exons_variances[!rownames(rpkms_Nascent_exons_variances) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons_variances <- rpkms_total_exons_variances[!rownames(rpkms_total_exons_variances) %in% filteroutGenes, ,drop=FALSE]
			}			
		}
		## sort according to time point ordering
		if( is.numeric(tpts) ) {
			ord <- order(tpts)
			tpts <- tpts[ord]
		} else {
			ord <- seq_along(tpts)
		}
		## estimate RNA dynamics
		totRpkmsIntEx <- list(
			exons=rpkms_total_exons[,ord,drop=FALSE]
			, exons_var=rpkms_total_exons_variances[,ord,drop=FALSE]
			)
		labeledRpkmsIntEx <- list(
			  exons=rpkms_Nascent_exons[,ord,drop=FALSE]
			, exons_var=rpkms_Nascent_exons_variances[,ord,drop=FALSE]
			)
		if( degDuringPulse ) {
			if( preexisting ) {
				stop('Pre-existing mode in still not implemeted when degDuringPulse is active.')
			}
			outSimple <- RNAdynamicsSimpleDDP(totRpkms=totRpkmsIntEx
											, labeledRpkms=labeledRpkmsIntEx
											, labeledSF=labeledSF
											, tpts=tpts
											, tL=labeling_time
											, BPPARAM=BPPARAM)														
		} else {
			outSimple <- RNAdynamicsSimple(totRpkms=totRpkmsIntEx
											, labeledRpkms=labeledRpkmsIntEx
											, labeledSF=labeledSF
											, tpts=tpts
											, tL=labeling_time
											, preexisting=preexisting
											, BPPARAM=BPPARAM)
		}
		## return the results in the form of an INSPEcT object.
		return(createInspectObject(outSimple, degDuringPulse=degDuringPulse))

	} else {

	###### Nascent mode, full:

	###### In case of the full nascent mode analysis, there could be a subset of the genes that have no introns
	###### and therefore must be analysed using the "simple" nascent mode. In this case the dataset must be 
	###### divided into full (genes with both exons and introns) and simple (genes with only introns),
	###### analysed separately and then merged. Genes from the simple mode will have NA in their 
	###### intronic quantification and splicing rate. 


		# Remove genes with intron quantification greater than the exon quantification, because they give rise to negative mature abundances.
		# (it can only apply to genes with both exons and introns)
		eiGenes <- intersect(rownames(rpkms_total_exons), rownames(rpkms_total_introns))
		negativeMature <- apply(rpkms_total_exons[eiGenes,,drop=FALSE]<rpkms_total_introns[eiGenes,,drop=FALSE],1,any)
		negativeNascent <- apply(rpkms_Nascent_exons[eiGenes,,drop=FALSE]<rpkms_Nascent_introns[eiGenes,,drop=FALSE],1,any)
		toRemove <- eiGenes[negativeMature|negativeNascent]
		if( length(toRemove)>0 ) {
			message(paste('Removing', length(toRemove), 'gene(s) with intronic quantifications greater than the exonic..'))
			rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% toRemove, ,drop=FALSE]
			rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% toRemove, ,drop=FALSE]
			rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% toRemove, ,drop=FALSE]
			rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% toRemove, ,drop=FALSE]
		}

		# Filter genes according to expression levels (in case the flag is active), separately for exonic and intronic features
		if( genesFilter ) {
			### filter out genes which have no signal in at least genesFilterThreshold of the time points in exonic features
			### in that case completely remove the gene (both intronic and exonic signal).
			ix1 <- apply(rpkms_Nascent_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_exons)
			ix2 <- apply(rpkms_total_exons, 1, function(x) length(which(x==0)))/ncol(rpkms_total_exons)
			filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1>genesFilterThreshold | ix2>genesFilterThreshold]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out', length(filteroutGenes), 'gene(s) with more than ',round(100*genesFilterThreshold),'% of zeros in their exonic quantifications.'))
				rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
			}
			### filter out genes which have no signal in at least genesFilterThreshold of the time points in intronic features
			### in that case just remove the intronic signal.
			ix3 <- apply(rpkms_Nascent_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_Nascent_introns)
			ix4 <- apply(rpkms_total_introns, 1, function(x) length(which(x==0)))/ncol(rpkms_total_introns)
			filteroutGenes <- rownames(rpkms_Nascent_introns)[ix3>genesFilterThreshold | ix4>genesFilterThreshold]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out intronic signal of', length(filteroutGenes), 'gene(s) with more than genesFilterThreshold of zero quantifications'))
				message('(for those genes only synthesis and degradation will be evaluated).')
				rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
			}
		} else {
			### filter out genes which have no signal in at least genesFilterThreshold of the time points in exonic features
			### in that case completely remove the gene (both intronic and exonic signal).
			ix1 <- apply(rpkms_Nascent_exons, 1, all)
			ix2 <- apply(rpkms_total_exons, 1, all)
			filteroutGenes <- rownames(rpkms_Nascent_exons)[ix1 | ix2]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out', length(filteroutGenes), 'gene(s) with more than ',round(100*genesFilterThreshold),'% of zeros in their exonic quantifications.'))
				rpkms_Nascent_exons <- rpkms_Nascent_exons[!rownames(rpkms_Nascent_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_exons <- rpkms_total_exons[!rownames(rpkms_total_exons) %in% filteroutGenes, ,drop=FALSE]
				rpkms_Nascent_introns <- rpkms_Nascent_introns[!rownames(rpkms_Nascent_introns) %in% filteroutGenes, ,drop=FALSE]
				rpkms_total_introns <- rpkms_total_introns[!rownames(rpkms_total_introns) %in% filteroutGenes, ,drop=FALSE]
			}
			### filter out genes which have no signal in at least genesFilterThreshold of the time points in intronic features
			### in that case just remove the intronic signal.
			ix3 <- apply(rpkms_Nascent_introns, 1, all)
			ix4 <- apply(rpkms_total_introns, 1, all)
			filteroutGenes <- rownames(rpkms_Nascent_introns)[ix3 | ix4]
			if( length(filteroutGenes)>0 ) {
				message(paste('Filtering out intronic signal of', length(filteroutGenes), 'gene(s) with more than genesFilterThreshold of zero quantifications'))
				message('(for those genes only synthesis and degradation will be evaluated).')
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

		## sort according to time point ordering
		if( is.numeric(tpts) ) {
			ord <- order(tpts)
			tpts <- tpts[ord]
		} else {
			ord <- seq_along(tpts)
		}
		# RNA dynamics for genes to be analysed in the full mode:
		message(paste('Number of genes with introns and exons: ', length(intExGenes)))
		totRpkmsIntEx <- list(
				exons=rpkms_total_exons[intExGenes, ord, drop=FALSE]
				, exons_var=rpkms_total_exons_variances[intExGenes, ord, drop=FALSE]
				, introns=rpkms_total_introns[intExGenes, ord, drop=FALSE]
				, introns_var=rpkms_total_introns_variances[intExGenes, ord, drop=FALSE]
			)
		labeledRpkmsIntEx <- list(
				exons=rpkms_Nascent_exons[intExGenes, ord, drop=FALSE]
				, exons_var=rpkms_Nascent_exons_variances[intExGenes, ord, drop=FALSE]
				, introns=rpkms_Nascent_introns[intExGenes, ord, drop=FALSE]
				, introns_var=rpkms_Nascent_introns_variances[intExGenes, ord, drop=FALSE]
			)
		if( degDuringPulse ) {
			if( preexisting ) {
				stop('Pre-existing mode in still not implemeted when degDuringPulse is active.')
			}
			outIntEx <- RNAdynamicsDDP(totRpkms=totRpkmsIntEx
									, labeledRpkms=labeledRpkmsIntEx
									, tpts=tpts
									, tL=labeling_time
									, BPPARAM=BPPARAM
									)				
		} else {
			if( preexisting ) {
				outIntEx <- RNAdynamicsFromPreex(preexRpkms=totRpkmsIntEx
										, labeledRpkms=labeledRpkmsIntEx
										, tpts=tpts
										, tL=labeling_time
										, BPPARAM=BPPARAM
										)
			} else {
				outIntEx <- RNAdynamics(totRpkms=totRpkmsIntEx
										, labeledRpkms=labeledRpkmsIntEx
										, tpts=tpts
										, tL=labeling_time
										, BPPARAM=BPPARAM
										)				
			}
		}
		# RNA dynamics for genes to be analysed in the simple mode (if they exist),
		# scaling the labeled library using the "labeledSF" caluclated from the full analysis:
		if( length(onlyExGenes)>0 ) {
			message(paste('Number of genes with only exons: ', length(onlyExGenes)))
			labeledSF <- outIntEx$labeledSF
			totRpkmsOnlyEx <- list(
					  exons=rpkms_total_exons[onlyExGenes, ord, drop=FALSE]
					, exons_var=rpkms_total_exons_variances[onlyExGenes, ord, drop=FALSE]
				)
			labeledRpkmsOnlyEx <- list(
					  exons=rpkms_Nascent_exons[onlyExGenes, ord, drop=FALSE]
					, exons_var=rpkms_Nascent_exons_variances[onlyExGenes, ord, drop=FALSE]
				)
			if( degDuringPulse ) {
				if( preexisting ) {
					stop('Pre-existing mode in still not implemeted when degDuringPulse is active.')
				}
				outOnlyEx <- RNAdynamicsSimpleDDP(totRpkms=totRpkmsOnlyEx
												, labeledRpkms=labeledRpkmsOnlyEx
												, labeledSF=labeledSF
												, tpts=tpts
												, tL=labeling_time
												, BPPARAM=BPPARAM)
			} else {
				outOnlyEx <- RNAdynamicsSimple(totRpkms=totRpkmsOnlyEx
												, labeledRpkms=labeledRpkmsOnlyEx
												, labeledSF=labeledSF
												, tpts=tpts
												, tL=labeling_time
												, preexisting=preexisting
												, BPPARAM=BPPARAM)			
			}
			# merge RNA dynamics coming from the full and the simples analysis
			outIntEx <- list(
				concentrations=list(
					total=rbind(outIntEx$concentrations$total, outOnlyEx$concentrations$total)
					, total_var=rbind(outIntEx$concentrations$total_var, outOnlyEx$concentrations$total_var)
					, preMRNA=rbind(outIntEx$concentrations$preMRNA, outOnlyEx$concentrations$preMRNA)
					, preMRNA_var=rbind(outIntEx$concentrations$preMRNA_var, outOnlyEx$concentrations$preMRNA_var)
					, labeled_total=rbind(outIntEx$concentrations$labeled_total, outOnlyEx$concentrations$labeled_total)
					, labeled_total_var=rbind(outIntEx$concentrations$labeled_total_var, outOnlyEx$concentrations$labeled_total_var)
					, labeled_preMRNA=rbind(outIntEx$concentrations$labeled_preMRNA, outOnlyEx$concentrations$labeled_preMRNA)
					, labeled_preMRNA_var=rbind(outIntEx$concentrations$labeled_preMRNA_var, outOnlyEx$concentrations$labeled_preMRNA_var)					
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
				, tL=outIntEx$tL
				)
		}
		## return the results in the form of an INSPEcT object.
		return(createInspectObject(outIntEx, NoNascent, NascentSimple, onlyExGenes, degDuringPulse))
	}
}

###############################
##### LOW-LEVEL FUNCTIONS ########
################################

matrixColNorm <- function(mat, norm) t(t(mat)*norm)		

createInspectObject <- function(out, NoNascent=FALSE, 
	NascentSimple=FALSE, onlyExGenes=character(0), degDuringPulse=FALSE) {

	# make an "ExpressionSet" object containing all the information
	nTpts <- length(out$tpts)
	tptsLabels = if( is.numeric(out$tpts) ) signif(out$tpts,9) else out$tpts
	exprData <- cbind(out$concentrations$total
					, out$concentrations$preMRNA
					, out$concentrations$labeled_total
					, out$concentrations$labeled_preMRNA
					, out$rates$alpha
					, out$rates$beta
					, out$rates$gamma)

	pData <- data.frame(feature=c(rep('total',nTpts)
								, rep('preMRNA',nTpts)
								, rep('labeled_total',nTpts)
								, rep('labeled_preMRNA',nTpts)
								, rep('synthesis',nTpts)
								, rep('degradation',nTpts)
								, rep('processing',nTpts) 
								)
					  , time=rep(out$tpts, 7))

	colnames(exprData) <- paste(pData$feature, tptsLabels, sep='_')
	rownames(pData) <- colnames(exprData)
	phenoData <- new('AnnotatedDataFrame', data=pData)

	ratesFirstGuess <- ExpressionSet(assayData=exprData
								   , phenoData=phenoData)
	featureNames(ratesFirstGuess) <- out$geneNames
	
	varData <- cbind(total=out$concentrations$total_var
					  , preMRNA=out$concentrations$preMRNA_var
					  , labeled_total=out$concentrations$labeled_total_var
					  , labeled_preMRNA=out$concentrations$labeled_preMRNA_var
					  , synthesis=out$rates$alpha_var
					  )

	pData <- data.frame(feature=c(rep('total',nTpts)
								, rep('preMRNA',nTpts)
								, rep('labeled_total',nTpts)
								, rep('labeled_preMRNA',nTpts)
								, rep('synthesis',nTpts)
								)
					  , time=rep(out$tpts, 5))

	colnames(varData) <- paste(pData$feature, tptsLabels, sep='_')
	rownames(pData) <- colnames(varData)
	phenoData <- new('AnnotatedDataFrame', data=pData)

	ratesFirstGuessVar <- ExpressionSet(assayData=varData
								   , phenoData=phenoData)
	featureNames(ratesFirstGuessVar) <- out$geneNames

	# Controls necessary to keep the output of newINSPEcT_NoNascent and newINSPEcT equal

	if( is.null(out$labeledSF) ) out$labeledSF <- numeric(0)
	if( is.null(out$tL) ) out$tL <- numeric(0)

	## update the object and return it
	object <- new('INSPEcT')
	object@tpts <- out$tpts
	object@labeledSF <- out$labeledSF
	object@tL <- out$tL	
	object@ratesFirstGuess <- ratesFirstGuess
	object@ratesFirstGuessVar <- ratesFirstGuessVar
	# object@precision <- out$ratesEstimPrec
	object@model@simple <- TRUE
	object@NoNascent <- NoNascent
	object@degDuringPulse <- degDuringPulse
	
	if( NoNascent ) {
		# object@ratesFirstGuessP <- out$ratesFirstGuessP
	} else {
		if( NascentSimple )
			object@model@simple <- TRUE
	}

	return(object)

}

RNAdynamicsSimpleDDP <- function(totRpkms, labeledRpkms, labeledSF, tpts, tL, BPPARAM=SerialParam()) 
{
		
	## retrieve gene names from rownames of exon total rpkms
	geneNames <- rownames(totRpkms$exons)
		
	##### total fraction
	## rename total exons
	Texo  <- totRpkms$exons
	Texo_var <- totRpkms$exons_var

	##### labeled fraction
	## rename labeled exons
	Lexo <- labeledRpkms$exons
	Lexo_var <- labeledRpkms$exons_var

	##### in case tpts is numeric (timecourse analysis)
	##### estimate the derivative
	if( is.numeric(tpts) ) {
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		TexoDer[, 1] <- 0
	##### otherwise put the derivative to zero
	} else {
		TexoDer = matrix(0, nrow=nrow(Texo), ncol=ncol(Texo))
	}
	message('Estimating all rates...')
	## degradation during pulse
	initOneGene <- function(i, j) {
		c(
			Lexo[i,j]/tL
			,(Lexo[i,j]/tL-TexoDer[i,j])/Texo[i,j]
		)
	}
	sysNascentSmall <- function(x) {
		y <- numeric(2)
		y[1] <- TexoDer[i,j] - x[1] + x[2]*Texo[i,j]
		y[2] <- Lexo[i,j]*labeledSF[j] - x[1]/x[2]*(1-exp(-x[2]*tL))
		y
	}
	## get only the rates
	nGenes <- nrow(Lexo)
	nTpts <- ncol(Lexo)
	ratesEstimPrec <- matrix(NA, nrow=nGenes, ncol=nTpts)
	alphaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
	betaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
	capture.output(suppressWarnings({
		for( i in 1:nGenes ) {
			for( j in 1:nTpts ) {
				init <- initOneGene(i,j)
				# if( all(is.finite(init)) ) {
				mrOut <- tryCatch(
					multiroot(
						sysNascentSmall
						, initOneGene(i,j)
						)
					, error=function(e)
					list(
						root=rep(NA, 2)
						, estim.precis=NA
					)
				)
				ratesEstimPrec[i,j] <- mrOut$estim.precis
				alphaTC[i,j] <- mrOut$root[1]
				betaTC[i,j] <- mrOut$root[2]
			}
		}
	}))
		
	ix <- alphaTC<0
	negativePerc <- length(which(ix))/length(ix)*100
	if( negativePerc>20 ) {
		warning(paste(round(negativePerc), '% of the synthesis rates are negative. Putting them to NA.'))
	}
	alphaTC[ix] <- NA
	## degradation
	ix <- betaTC<0
	negativePerc <- length(which(ix))/length(ix)*100
	if( negativePerc>20 ) {
		warning(paste(round(negativePerc), '% of the degradation rates are negative. Putting them to NA.'))
	}
	betaTC[ix] <- NA

	# impute NA values (in case of time course)
	if( is.numeric(tpts) ) {
		alphaTC <- do.call('rbind',bplapply(1:nrow(alphaTC), 
			function(i) impute_na_tc(tpts, alphaTC[i,]), BPPARAM=BPPARAM))
		betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
			function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=BPPARAM))		
	}
		
	## set preMRNA and gamma to NA
	Tint <- matrix(NA, nrow(Texo), ncol(Texo))
	Tint_var <- matrix(NA, nrow(Texo), ncol(Texo))
	Lint <- matrix(NA, nrow(Tint), ncol(Tint))
	Lint_var <- matrix(NA, nrow(Tint), ncol(Tint))
	gammaTC <- matrix(NA, nrow(betaTC), ncol(betaTC))
	## and also alphaTC variance because in DDP is not possible to calculate
	alphaTC_var <- matrix(NA, nrow(alphaTC), ncol(alphaTC))
	
	## remove dimnames from matrices and names from vectors 
	## that will be output
	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(Texo, 'dimnames') <- NULL
	attr(Tint, 'dimnames') <- NULL
	attr(Lexo, 'dimnames') <- NULL
	attr(alphaTC_var, 'dimnames') <- NULL
	attr(Texo_var, 'dimnames') <- NULL
	attr(Tint_var, 'dimnames') <- NULL
	attr(Lexo_var, 'dimnames') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL
	
	## return scaled results and rates (adding NA values 
	## for preMRNA and gamma)

	if(!is.matrix(Texo_var)){Texo_var <- matrix(rep(Texo_var,length(tpts)),nrow=length(Texo_var),ncol=length(tpts))}
	if(!is.matrix(Tint_var)){Tint_var <- matrix(rep(Tint_var,length(tpts)),nrow=length(Tint_var),ncol=length(tpts))}

	return(list(
		concentrations=list(
			total=Texo
			, total_var=Texo_var
			, preMRNA=Tint
			, preMRNA_var=Tint_var
			, labeled_total=Lexo
			, labeled_total_var=Lexo_var
			, labeled_preMRNA=Lint
			, labeled_preMRNA_var=Lint_var
		)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, geneNames=geneNames
		, tpts=tpts
		, tL=tL))
	
}

RNAdynamicsSimple <- function(totRpkms, labeledRpkms, labeledSF, tpts, tL, preexisting, BPPARAM=SerialParam()) 
{
		
	## retrieve gene names from rownames of exon total rpkms
	geneNames <- rownames(totRpkms$exons)

	##### total fraction
	## rename total exons
	if( preexisting ) {
		Texo  <- totRpkms$exons + matrixColNorm(labeledRpkms$exons, labeledSF)
		Texo_var <- totRpkms$exons_var + matrixColNorm(labeledRpkms$exons_var, labeledSF^2)
	} else {
		Texo  <- totRpkms$exons
		Texo_var <- totRpkms$exons_var		
	}

	##### labeled fraction
	## rename labeled exons
	Lexo <- labeledRpkms$exons
	Lexo_var <- labeledRpkms$exons_var

	## calculate alpha and recalculate the variance
	alphaTC <- matrixColNorm(Lexo, labeledSF/tL)
	alphaTC_var <- matrixColNorm(Lexo_var, (labeledSF/tL)^2)

	if( is.numeric(tpts) ) {

		##### estimate the derivative
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		TexoDer[, 1] <- 0
						
		## calculate beta
		message('Estimating degradation rates...')
		inferKBetaFromIntegral <- function(tpts, alpha, totalRNA, maxBeta=75) 
		{
			solveFun <- function(beta, t0, t1, alpha_t0, alpha_t1, X_t0, X_t1 ) 
			{
				m <- (alpha_t0 - alpha_t1 ) / (t0 - t1 )
				q <- alpha_t0 - m * t0
				X_t1 - X_t0*exp(-beta*(t1 - t0)) - (
				(m*t1*beta + q*beta - m ) / (beta^2) - 
				(m*t0*beta + q*beta - m ) * exp(-beta*(t1 - t0)) / (beta^2)
				)
			}
			bplapply(2:length(tpts), function(j)
				lapply(1:nrow(alpha), function(i) {
				tryCatch(
					uniroot(solveFun
						, c(1e-5, maxBeta)
						, t0 = tpts[j-1]
						, t1 = tpts[j]
						, alpha_t0 = alpha[i,j-1]
						, alpha_t1 = alpha[i,j]
						, X_t0 = totalRNA[i,j-1]
						, X_t1 = totalRNA[i,j]
						)
					, error=function(e) return(list(root=NA, estim.prec=NA, error=e))
				)})
			, BPPARAM=BPPARAM)
		}
		betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / Texo[,1]
		betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
		if( length(tpts)>1 ) {
			betaOut <- inferKBetaFromIntegral(tpts, alphaTC, Texo, 
				maxBeta=quantile(betaT0, na.rm=TRUE, probs=.99)*10
				)
			betaTC <- cbind(betaT0, 
				sapply(betaOut, function(x) sapply(x, '[[', 'root'))
				)
			ratesEstimPrec <- cbind(0,
				sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
				)
		} else {
			betaTC <- as.matrix(betaT0)
			ratesEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
		}

		## impute NA values
		betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
			function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=BPPARAM))

	} else {

		betaTC <- alphaTC/Texo
		ratesEstimPrec <- matrix(NA, nrow=nrow(betaTC), ncol=ncol(betaTC))

	}

	## set preMRNA and gamma to NA
	Tint <- matrix(NA, nrow(Texo), ncol(Texo))
	Tint_var <- matrix(NA, nrow(Texo), ncol(Texo))
	Lint <- matrix(NA, nrow(Tint), ncol(Tint))
	Lint_var <- matrix(NA, nrow(Tint), ncol(Tint))
	gammaTC <- matrix(NA, nrow(betaTC), ncol(betaTC))
	
	## remove dimnames from matrices and names from vectors 
	## that will be output
	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(Texo, 'dimnames') <- NULL
	attr(Tint, 'dimnames') <- NULL
	attr(Lexo, 'dimnames') <- NULL
	attr(alphaTC_var, 'dimnames') <- NULL
	attr(Texo_var, 'dimnames') <- NULL
	attr(Tint_var, 'dimnames') <- NULL
	attr(Lexo_var, 'dimnames') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL
	
	## return scaled results and rates (adding NA values 
	## for preMRNA and gamma)

	if(!is.matrix(Texo_var)){Texo_var <- matrix(rep(Texo_var,length(tpts)),nrow=length(Texo_var),ncol=length(tpts))}
	if(!is.matrix(Tint_var)){Tint_var <- matrix(rep(Tint_var,length(tpts)),nrow=length(Tint_var),ncol=length(tpts))}

	return(list(
		concentrations=list(
			total=Texo
			, total_var=Texo_var
			, preMRNA=Tint
			, preMRNA_var=Tint_var
			, labeled_total=Lexo
			, labeled_total_var=Lexo_var
			, labeled_preMRNA=Lint
			, labeled_preMRNA_var=Lint_var
		)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, geneNames=geneNames
		, tpts=tpts
		, tL=tL))

}

#########################
##### local functions ######
########################
fLint <- function(Lexo,gamma,tL) Lexo/(tL*gamma)*(1-exp(-gamma*tL))
# given the number of labeled molecules, gamma and tL gives back the 
# number of processed molecules
fGamma <- function(Lint, Lexo, tL, maxGamma=1e3) {
# given the number of labeled molecules and the number of processed
# molecules gives back the processing rate. It's an inverse function,
# therefore, the precision (step of the interval evaluated) and the
# max value where gamma is possibly evaluated should be provided
	if( is.na(Lint) | is.na(Lexo) ) return(NA)
	if( Lint >= Lexo ) return(0)
	if( Lint == 0 ) return(Inf)
	errorfun <- function(gamma, Lexo, Lint, tL ) (Lint-fLint(Lexo, gamma, tL))^2
	
	optimize(errorfun, c(0,maxGamma), Lexo=Lexo, Lint=Lint, tL=tL)$minimum
}
# calculate the factor which bring the median of the residuals between
# the modeled preMRNA levels and the measured to zero
sq.median.resids <- function(sf, P, dP, alpha, gamma) sapply(sf, function(i) {
	t1 <- dP + gamma*P
	t2 <- i*alpha
	idx <- is.finite(t1) & is.finite(t2) & t1 > 0 & t2 > 0
	resids <- t1[idx] - t2[idx]
	stats::median(resids , na.rm=TRUE)^2
})

RNAdynamicsDDP <- function(totRpkms, labeledRpkms, tpts, tL, simulatedData=FALSE, BPPARAM=SerialParam()) 
{
	
	## the control during the workflow will be based on the 
	## negation of simulatedData (that is realData)
	realData <- !simulatedData
		
	## retrieve gene names from rownames of exon total rpkms
	geneNames <- rownames(totRpkms$exons)
		
	##### total fraction
	## rename total exons
	Texo  <- totRpkms$exons
	Texo_var <- totRpkms$exons_var
	## rename total introns
	Tint  <- totRpkms$introns
	Tint_var <- totRpkms$introns_var

	##### labeled fraction
	## rename labeled exons
	Lexo <- labeledRpkms$exons
	Lexo_var <- labeledRpkms$exons_var
	## rename total introns
	Lint <- labeledRpkms$introns
	Lint_var <- labeledRpkms$introns_var

	# if tpts is numeric (timecourse analysis) estimate
	# derivatives from time course
	if( is.numeric(tpts) ) {
		TintDer <- as.matrix(sapply(1:nrow(Tint), 
			function(i) {
				if( all(is.finite(Tint[i,] ) ) ) {
					spfun <- splinefun(tpts, Tint[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))

		if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
		TintDer[, 1] <- 0 
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		TexoDer[, 1] <- 0
	# otherwise put the derivatives to zero	
	} else {
		TintDer= matrix(0, nrow=nrow(Tint), ncol=ncol(Tint))
		TexoDer= matrix(0, nrow=nrow(Texo), ncol=ncol(Texo))
	}
	
	if( realData ){
	
		#########################
		#### scale Nascent rpkms ###
		#########################
		
		message('Calculating scaling factor between total and Nascent libraries...')
		
		##################
		#### scale data ###
		##################
		# preMRNA derivative as splines 
		# (force the first time point to have derivative zero )
		# estimate of alpha and gamma from Nascent data

		gammaTC <- do.call('cbind',bplapply(1:length(tpts), function(j) 
			sapply(1:nrow(Lint), function(i, Lint, Lexo, tL) 
				fGamma(Lint[i,j] , Lexo[i,j] , tL)
			, Lint=Lint, Lexo=Lexo, tL=tL)
		,BPPARAM=BPPARAM))

		# scale factor 
		labeledSF_prior <- sapply(1:ncol(Tint), function(j)
			optimize(sq.median.resids, c(0.01,100), P=Tint[,j], dP=TintDer[,j], 
			alpha=Lexo[,j]/tL, gamma=gammaTC[,j] )$minimum)

		## re-calculate scaling factor using the DDP framework

		## select the 500 most synthesized genes at time zero
		## to calculate the scaling factor
		if( nrow(Lexo) > 500 ) {
			geneSubset <- order(Texo[,1],decreasing=TRUE)[1:500]
		} else {
			geneSubset <- 1:nrow(Lexo)
		}
		
		## degradation during pulse
		# system with Nascent scaling factor as a 4th variable to be identified
		initOneGene <- function(i, j) {
			c(Lexo[i,j]/tL
			,(Lexo[i,j]/tL-TexoDer[i,j])/(Texo[i,j]-Tint[i,j])
			,(Lexo[i,j]/tL-TintDer[i,j])/Tint[i,j]
			)
		}

		sysNascentScale <- function(x) {
			y <- numeric(4)
			y[1] <- TintDer[i,j] - x[1] + x[3]*Tint[i,j]
			y[2] <- TexoDer[i,j] - x[1] + x[2]*(Texo[i,j]-Tint[i,j])
			y[3] <- x[4]*Lint[i,j] - x[1]/x[3]*(1-exp(-x[3]*tL))
			y[4] <- x[4]*Lexo[i,j] - (x[1]*exp(-x[2]*tL)*(x[3]^2-x[2]^2*exp((x[2]-x[3])*tL)+
			x[2]^2*exp(x[2]*tL)-x[3]^2*exp(x[2]*tL)))/(x[2]*(x[2]-x[3])*x[3])
			y/c(x[1],x[1],x[4]*Lint[i,j],x[4]*Lexo[i,j])
		}
		## get only Nascent scales
		nGenes <- nrow(Lexo)
		nTpts <- ncol(Lexo)

		labeledSfep <- matrix(NA, nrow=nGenes, ncol=nTpts)
		labeledSf <- matrix(NA, nrow=nGenes, ncol=nTpts)

		capture.output(suppressWarnings({
			for( i in geneSubset ) {
				for( j in 1:nTpts ) {
					init <- initOneGene(i,j)*labeledSF_prior[j]
					mrOut <- tryCatch(
						multiroot(sysNascentScale, c(init,labeledSF_prior[j]))
						, error=function(e) list(root=rep(NA, 4), estim.precis=NA)
						)
					labeledSf[i,j] <- mrOut$root[4]
					labeledSfep[i,j] <- mrOut$estim.precis
				}
			}
		}))

		## chose the best resolved genes to estimate 
		## from them the scale factor
		epTsh <- apply(labeledSfep,2,quantile,probs=.75,na.rm=TRUE)
		ix <- t(apply(labeledSfep,1, function(x) x>epTsh))
		labeledSf[ix] <- NA
		labeledSF <- apply(labeledSf, 2, stats::median, na.rm=TRUE)
		
		
		## simulated data
	} else {
	
		# in case of synthetic data the time course is already scaled
		# therefore just rename the variables and compute varince in 
		# case is not provided
		alphaTC <- Lexo
		alphaTC_var <- Lexo_var

		# scaling factor for synthetic dataset is meaningless
		labeledSF <- rep(1, length(tpts))
	
	}
	
	################################
	## estimate degradation rates ##
	################################
			
	message('Estimating all rates...')
	
	## once set the scale factor calculate the rates
	initOneGene <- function(i, j) {
		c(
			Lexo[i,j]*labeledSF[j]/tL
			,(Lexo[i,j]*labeledSF[j]/tL-TexoDer[i,j])/(Texo[i,j]-Tint[i,j])
			,(Lexo[i,j]*labeledSF[j]/tL-TintDer[i,j])/Tint[i,j]
		)
	}

	sysNascent <- function(x, labeledSF) {
		y <- numeric(3)
		y[1] <- TintDer[i,j] - x[1] + x[3]*Tint[i,j]
		y[2] <- TexoDer[i,j] - x[1] + x[2]*(Texo[i,j]-Tint[i,j])
		y[3] <- labeledSF[j]*Lexo[i,j] - (x[1]*exp(-x[2]*tL)*(x[3]^2-x[2]^2*exp((x[2]-x[3])*tL)+
		x[2]^2*exp(x[2]*tL)-x[3]^2*exp(x[2]*tL)))/(x[2]*(x[2]-x[3])*x[3])
		y

	}
	## get only the rates
	nGenes <- nrow(Lexo)
	nTpts <- ncol(Lexo)
	ratesEstimPrec <- matrix(NA, nrow=nGenes, ncol=nTpts)
	alphaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
	betaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
	gammaTC <- matrix(NA, nrow=nGenes, ncol=nTpts)
	capture.output(suppressWarnings({
		for( i in 1:nGenes ) {
			for( j in 1:nTpts ) {
				mrOut <- tryCatch(
					multiroot(sysNascent, initOneGene(i,j), labeledSF=labeledSF)
						, error=function(e) list(root=rep(NA, 4), estim.precis=NA)
					)
				ratesEstimPrec[i,j] <- mrOut$estim.precis
				alphaTC[i,j] <- mrOut$root[1]
				betaTC[i,j] <- mrOut$root[2]
				gammaTC[i,j] <- mrOut$root[3]
			}
		}
	}))
	
	## put negative values to NA and rise a 
	## warning if they are more than 20% of a specific rate
	ix <- alphaTC<0 | betaTC<0 | gammaTC<0
	negativePerc <- length(which(ix))/length(ix)*100
	if( negativePerc>20 ) {
		warning(paste(round(negativePerc), 
		'% of the genes contains negative rates. Putting them to NA.'))
	}
	alphaTC[ix] <- NA
	betaTC[ix] <- NA
	gammaTC[ix] <- NA
	ratesEstimPrec[ix] <- NA

	# impute NA values (in case of time course)
	if( is.numeric(tpts) ) {
		alphaTC <- do.call('rbind',bplapply(1:nrow(alphaTC), 
			function(i) impute_na_tc(tpts, alphaTC[i,]), BPPARAM=BPPARAM))
		betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
			function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=BPPARAM))
		gammaTC <- do.call('rbind',bplapply(1:nrow(gammaTC), 
			function(i) impute_na_tc(tpts, gammaTC[i,]), BPPARAM=BPPARAM))
	}
		
	## return NA variance associated to alphaTC because in this mode 
	## the variance is not directly associated to the rate
	alphaTC_var <- matrix(NA, nrow(alphaTC), ncol(alphaTC))
	
	## remove dimnames from matrices and names from vectors 
	## that will be output
	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(Texo, 'dimnames') <- NULL
	attr(Tint, 'dimnames') <- NULL
	attr(Lexo, 'dimnames') <- NULL
	attr(Lint, 'dimnames') <- NULL
	attr(alphaTC_var, 'dimnames') <- NULL
	attr(Texo_var, 'dimnames') <- NULL
	attr(Tint_var, 'dimnames') <- NULL
	attr(Lexo_var, 'dimnames') <- NULL
	attr(Lint_var, 'dimnames') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL
	
	## return scaled results and rates (adding NA values 
	## for preMRNA and gamma)

	if(!is.matrix(Texo_var)){Texo_var <- matrix(rep(Texo_var,length(tpts)),nrow=length(Texo_var),ncol=length(tpts))}
	if(!is.matrix(Tint_var)){Tint_var <- matrix(rep(Tint_var,length(tpts)),nrow=length(Tint_var),ncol=length(tpts))}

	return(list(
		concentrations=list(
			total=Texo
			, total_var=Texo_var
			, preMRNA=Tint
			, preMRNA_var=Tint_var
			, labeled_total=Lexo
			, labeled_total_var=Lexo_var
			, labeled_preMRNA=Lint
			, labeled_preMRNA=Lint
			, labeled_preMRNA_var=Lint_var
			, labeled_preMRNA_var=Lint_var
		)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, geneNames=geneNames
		, labeledSF=labeledSF
		, tpts=tpts
		, tL=tL))
	
}

RNAdynamics <- function(totRpkms, labeledRpkms, tpts, tL, simulatedData=FALSE, BPPARAM=SerialParam()) 
{
	
	## the control during the workflow will be based on the 
	## negation of simulatedData (that is realData)
	realData <- !simulatedData
		
	## retrieve gene names from rownames of exon total rpkms
	geneNames <- rownames(totRpkms$exons)
		
	##### total fraction
	## rename total exons
	Texo  <- totRpkms$exons
	Texo_var <- totRpkms$exons_var
	## rename total introns
	Tint  <- totRpkms$introns
	Tint_var <- totRpkms$introns_var

	##### labeled fraction
	## rename labeled exons
	Lexo <- labeledRpkms$exons
	Lexo_var <- labeledRpkms$exons_var
	## rename total introns
	Lint <- labeledRpkms$introns
	Lint_var <- labeledRpkms$introns_var

	# if tpts is numeric (timecourse analysis) estimate
	# derivatives from time course
	if( is.numeric(tpts) ) {
		TintDer <- as.matrix(sapply(1:nrow(Tint), 
			function(i) {
				if( all(is.finite(Tint[i,] ) ) ) {
					spfun <- splinefun(tpts, Tint[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		# if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
		TintDer <- t(TintDer)
		TintDer[, 1] <- 0 
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		# if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		TexoDer <- t(TexoDer)
		TexoDer[, 1] <- 0
	# otherwise put the derivatives to zero	
	} else {
		TintDer= matrix(0, nrow=nrow(Tint), ncol=ncol(Tint))
		TexoDer= matrix(0, nrow=nrow(Texo), ncol=ncol(Texo))
	}
	
	if( realData ){
					
		#########################
		#### scale Nascent rpkms ###
		#########################
		
		message('Calculating scaling factor between total and Nascent libraries...')
		
		##################
		#### scale data ###
		##################
		# preMRNA derivative as splines 
		# (force the first time point to have derivative zero )
		# estimate of alpha and gamma from Nascent data

		gammaTC <- do.call('cbind',bplapply(1:length(tpts), function(j) 
			sapply(1:nrow(Lint), function(i, Lint, Lexo, tL) 
				fGamma(Lint[i,j] , Lexo[i,j] , tL)
			, Lint=Lint, Lexo=Lexo, tL=tL)
		,BPPARAM=BPPARAM))

		# scale factor 
		labeledSF <- sapply(1:ncol(Tint), function(j)
			optimize(sq.median.resids, c(0.01,100), P=Tint[,j], dP=TintDer[,j], 
			alpha=Lexo[,j]/tL, gamma=gammaTC[,j] )$minimum)
		
	} else { ## simulated data
	
		# in case of synthetic data the time course is already scaled
		# therefore just rename the variables and compute varince in 
		# case is not provided
		alphaTC <- Lexo
		alphaTC_var <- Lexo_var

		# scaling factor for synthetic dataset is meaningless
		labeledSF <- rep(1, length(tpts))

		# assign NaN values to Lint e Lint_var
		Lint <- matrix(NA, nrow(Tint), ncol(Tint))
		Lint_var <- matrix(NA, nrow(Tint), ncol(Tint))
	
	}
	
	##############################
	## estimate synthesis rates ##
	##############################
			
	## calculate alpha and recalculate the variance
	alphaTC <- matrixColNorm(Lexo, labeledSF/tL)
	alphaTC_var <- matrixColNorm(Lexo_var, (labeledSF/tL)^2)

	if( is.numeric(tpts) ) {

		################################
		## estimate degradation rates ##
		################################
						
		message('Estimating degradation rates...')
		betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / (Texo[,1] - Tint[,1] )
		betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
		if( length(tpts)>1 ) {
			betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, Texo, Tint, 
				maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
				)
			if( nrow(Texo)==1 ) {
				betaTC <- t(c(betaT0, sapply(betaOut, function(x) sapply(x, '[[', 'root'))))
				betaEstimPrec <- t(c(0,sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))))
			} else {
				betaTC <- cbind(betaT0, sapply(betaOut, function(x) sapply(x, '[[', 'root')))
				betaEstimPrec <- cbind(0,sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec')))
			}
		} else {
			betaTC <- as.matrix(betaT0)
			betaEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
		}

		###################################
		## estimate processing rates #########
		#####################################
		
		# calculate gamma (from  total RNA introns and alphas )
		message('Estimating processing rates...')
		gammaT0 <- ( alphaTC[,1] - TintDer[,1] ) / Tint[,1]
		gammaT0[gammaT0 < 0 | !is.finite(gammaT0)] <- NA
		if( length(tpts)>1 ) {
			gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, Tint, 
				maxGamma=quantile(gammaT0,na.rm=TRUE,probs=.99)*10, BPPARAM=BPPARAM
				)
			if( nrow(Texo)==1 ) {
				gammaTC <- t(c(gammaT0, sapply(gammaOut, function(x) sapply(x, '[[', 'root'))))
				gammaEstimPrec <- t(c(0, sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec'))))
			} else {
				gammaTC <- cbind(gammaT0, sapply(gammaOut, function(x) sapply(x, '[[', 'root')))
				gammaEstimPrec <- cbind(0, sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec')))
			}
		} else {
			gammaTC <- as.matrix(gammaT0)
			gammaEstimPrec <- matrix(0, nrow=nrow(gammaTC), ncol=ncol(gammaTC))
		}

		# ## impute NA values
		betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
			function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=BPPARAM))
		gammaTC <- do.call('rbind',bplapply(1:nrow(gammaTC), 
			function(i) impute_na_tc(tpts, gammaTC[i,]), BPPARAM=BPPARAM))

		ratesEstimPrec <- betaEstimPrec + gammaEstimPrec

	} else {

		betaTC= alphaTC/(Texo-Tint)
		gammaTC= alphaTC/Tint
		ratesEstimPrec= matrix(NA, nrow=nrow(betaTC), ncol=ncol(betaTC))

	}
	
	## remove dimnames from matrices and names from vectors 
	## that will be output
	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(Texo, 'dimnames') <- NULL
	attr(Tint, 'dimnames') <- NULL
	attr(Lexo, 'dimnames') <- NULL
	attr(Lint, 'dimnames') <- NULL
	attr(alphaTC_var, 'dimnames') <- NULL
	attr(Texo_var, 'dimnames') <- NULL
	attr(Tint_var, 'dimnames') <- NULL
	attr(Lexo_var, 'dimnames') <- NULL
	attr(Lint_var, 'dimnames') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL
	
	## return scaled results and rates (adding NA values 
	## for preMRNA and gamma)

	if(!is.matrix(Texo_var)){Texo_var <- matrix(rep(Texo_var,length(tpts)),nrow=length(Texo_var),ncol=length(tpts))}
	if(!is.matrix(Tint_var)){Tint_var <- matrix(rep(Tint_var,length(tpts)),nrow=length(Tint_var),ncol=length(tpts))}

	return(list(
		concentrations=list(
			total=Texo
			, total_var=Texo_var
			, preMRNA=Tint
			, preMRNA_var=Tint_var
			, labeled_total=Lexo
			, labeled_total_var=Lexo_var
			, labeled_preMRNA=Lint
			, labeled_preMRNA_var=Lint_var
		)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, geneNames=geneNames
		, labeledSF=labeledSF
		, tpts=tpts
		, tL=tL))
	
}

RNAdynamicsFromPreex <- function(preexRpkms, labeledRpkms, tpts, tL, BPPARAM=SerialParam()) 
{
			
	## retrieve gene names from rownames of exon total rpkms
	geneNames <- rownames(preexRpkms$exons)
		
	##### total fraction
	## rename total exons
	Pexo  <- preexRpkms$exons
	Pexo_var <- preexRpkms$exons_var
	## rename total introns
	Pint  <- preexRpkms$introns
	Pint_var <- preexRpkms$introns_var

	##### labeled fraction
	## rename labeled exons
	Lexo <- labeledRpkms$exons
	Lexo_var <- labeledRpkms$exons_var
	## rename total introns
	Lint <- labeledRpkms$introns
	Lint_var <- labeledRpkms$introns_var
						
	#########################
	#### scale Nascent rpkms ###
	#########################
	
	message('Calculating scaling factor between Pre-existing and Nascent libraries...')

	estimateSFpreexisting <- function(yf, labeledSF, Lint, Lexo, Pint, gammaTC, tL, tpts, j, BPPARAM=SerialParam()) {

		labeledSF[j] <- yf

		Tint <- matrixColNorm(Lint, labeledSF) + Pint
		Lint <- matrixColNorm(Lint, labeledSF)
		Lexo <- matrixColNorm(Lexo, labeledSF)

		if( is.numeric(tpts) & j > 1 ) {
			TintDer <- as.matrix(sapply(1:nrow(Tint), 
				function(i) {
					if( all(is.finite(Tint[i,] ) ) ) {
						spfun <- splinefun(tpts, Tint[i,])
						return(spfun(tpts, deriv=1) )
					} else return(rep(NA, length(tpts)) )
				}
			))
			if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
		} else {
			TintDer= matrix(0, nrow=nrow(Tint), ncol=ncol(Tint))
		}

		# scale factor 
		optimize(sq.median.resids, c(0.01,100), P=Tint[,j], dP=TintDer[,j], 
			alpha=Lexo[,j]/tL, gamma=gammaTC[,j] )$minimum

	}

	# estimate of alpha and gamma from 4sU data 
	# (they are independent of a scaling factor that is common to 
	# exonic and intronic signal)
	gammaTC <- do.call('cbind',bplapply(
		1:ncol(Lexo), function(j) 
			sapply(1:nrow(Lint), 
				function(i, Lint, Lexo, tL) fGamma(Lint[i,j] , Lexo[i,j] , tL)
				, Lint=Lint, Lexo=Lexo, tL=tL)
			,BPPARAM=BPPARAM))

	labeledSF <- rep(1, ncol(Lint))

	for( j in 1:ncol(Lint) ) {

		labeledSF[j] <- optimize(
			f = function(yf) 
				sum((estimateSFpreexisting(yf, labeledSF, Lint, Lexo, 
					Pint, gammaTC, tL, tpts, j, BPPARAM) - 1)^2),
			interval=c(.1,10)
			)$minimum

	}

	##### total fraction
	## rename total exons
	Texo <- Pexo + matrixColNorm(Lexo, labeledSF)
	Texo_var <- Pexo_var + matrixColNorm(Lexo, labeledSF^2)
	## rename total introns
	Tint <- Pint + matrixColNorm(Lint, labeledSF)
	Tint_var <- Pint_var + matrixColNorm(Lint, labeledSF^2)

	# if tpts is numeric (timecourse analysis) estimate
	# derivatives from time course
	if( is.numeric(tpts) ) {
		TintDer <- as.matrix(sapply(1:nrow(Tint), 
			function(i) {
				if( all(is.finite(Tint[i,] ) ) ) {
					spfun <- splinefun(tpts, Tint[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		if( ncol(TintDer)>1 ) TintDer <- t(TintDer)
		TintDer[, 1] <- 0 
		TexoDer <- as.matrix(sapply(1:nrow(Texo), 
			function(i) {
				if( all(is.finite(Texo[i,] ) ) ) {
					spfun <- splinefun(tpts, Texo[i,])
					return(spfun(tpts, deriv=1) )
				} else return(rep(NA, length(tpts)) )
			}
		))
		if( ncol(TexoDer)>1 ) TexoDer <- t(TexoDer)
		TexoDer[, 1] <- 0
	# otherwise put the derivatives to zero	
	} else {
		TintDer= matrix(0, nrow=nrow(Tint), ncol=ncol(Tint))
		TexoDer= matrix(0, nrow=nrow(Texo), ncol=ncol(Texo))
	}

	##############################
	## estimate synthesis rates ##
	##############################
			
	## calculate alpha and recalculate the variance
	alphaTC <- matrixColNorm(Lexo, labeledSF/tL)
	alphaTC_var <- matrixColNorm(Lexo_var, (labeledSF/tL)^2)

	if( is.numeric(tpts) ) {

		################################
		## estimate degradation rates ##
		################################
						
		message('Estimating degradation rates...')
		betaT0 <- ( alphaTC[,1] - TexoDer[,1] ) / (Texo[,1] - Tint[,1] )
		betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
		if( length(tpts)>1 ) {
			betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, Texo, Tint, 
				maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
				)
			betaTC <- cbind(betaT0, 
				sapply(betaOut, function(x) sapply(x, '[[', 'root'))
				)
			betaEstimPrec <- cbind(0,
				sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))
				)
		} else {
			betaTC <- as.matrix(betaT0)
			betaEstimPrec <- matrix(0, nrow=nrow(betaTC), ncol=ncol(betaTC))
		}

		###################################
		## estimate processing rates #########
		#####################################
		
		# calculate gamma (from  total RNA introns and alphas )
		message('Estimating processing rates...')
		gammaT0 <- ( alphaTC[,1] - TintDer[,1] ) / Tint[,1]
		gammaT0[gammaT0 < 0 | !is.finite(gammaT0)] <- NA
		if( length(tpts)>1 ) {
			gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, Tint, 
				maxGamma=quantile(gammaT0,na.rm=TRUE,probs=.99)*10, BPPARAM=BPPARAM
				)
			gammaTC <- cbind(gammaT0, 
				sapply(gammaOut, function(x) sapply(x, '[[', 'root'))
				)
			gammaEstimPrec <- cbind(0, 
				sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec'))
				)
		} else {
			gammaTC <- as.matrix(gammaT0)
			gammaEstimPrec <- matrix(0, nrow=nrow(gammaTC), ncol=ncol(gammaTC))
		}

		# ## impute NA values
		betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
			function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=BPPARAM))
		gammaTC <- do.call('rbind',bplapply(1:nrow(gammaTC), 
			function(i) impute_na_tc(tpts, gammaTC[i,]), BPPARAM=BPPARAM))

		ratesEstimPrec <- betaEstimPrec + gammaEstimPrec

	} else {

		betaTC= alphaTC/(Texo-Tint)
		gammaTC= alphaTC/Tint
		ratesEstimPrec= matrix(NA, nrow=nrow(betaTC), ncol=ncol(betaTC))

	}
	
	## remove dimnames from matrices and names from vectors 
	## that will be output
	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(Texo, 'dimnames') <- NULL
	attr(Tint, 'dimnames') <- NULL
	attr(Lexo, 'dimnames') <- NULL
	attr(Lint, 'dimnames') <- NULL
	attr(alphaTC_var, 'dimnames') <- NULL
	attr(Texo_var, 'dimnames') <- NULL
	attr(Tint_var, 'dimnames') <- NULL
	attr(Lexo_var, 'dimnames') <- NULL
	attr(Lint_var, 'dimnames') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL
	
	## return scaled results and rates (adding NA values 
	## for preMRNA and gamma)

	if(!is.matrix(Texo_var)){Texo_var <- matrix(rep(Texo_var,length(tpts)),nrow=length(Texo_var),ncol=length(tpts))}
	if(!is.matrix(Tint_var)){Tint_var <- matrix(rep(Tint_var,length(tpts)),nrow=length(Tint_var),ncol=length(tpts))}

	return(list(
		concentrations=list(
			total=Texo
			, total_var=Texo_var
			, preMRNA=Tint
			, preMRNA_var=Tint_var
			, labeled_total=Lexo
			, labeled_total_var=Lexo_var
			, labeled_preMRNA=Lint
			, labeled_preMRNA_var=Lint_var
		)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, geneNames=geneNames
		, labeledSF=labeledSF
		, tpts=tpts
		, tL=tL))
	
}

mcsapply <- function( X, FUN, ... ) do.call('cbind', bplapply( X, FUN, ... ))

lineCoefficients_NoNascent <- function(xInitial
								  ,xFinal
								  ,yInitial
								  ,yFinal)
{
  return(c(m = (yFinal-yInitial)/(xFinal-xInitial)
          ,q = (yInitial*xFinal-yFinal*xInitial)/(xFinal-xInitial)))
}

firstStep_NoNascent <- function(tpts
						   ,mature
						   ,premature
						   ,matureVariance
						   ,Dmin
						   ,Dmax) 
{
	fits <- lapply(1:nrow(mature), function(row)
  	{
	    optimize(firstStepError_NoNascent
	            ,c(Dmin, Dmax)
	            ,k2K3Ratio = mature[row,1]/premature[row,1]
	            ,tpts = tpts
	            ,premature = premature[row,]
	            ,mature = mature[row,]
	            ,matureVariance = matureVariance[row,])
   })
   
	out <- t(sapply(fits,unlist))
	out[,2] <- pchisq(out[,2], length(tpts)-1)
	colnames(out) <- c('k3','p')
	return(out)
}

firstStepError_NoNascent <- function(tpts
								,mature
								,premature
								,matureVariance
								,k3
								,k2K3Ratio)
{
	matureEstimated <- numeric(length = length(mature))
	matureEstimated[1] <- mature[1]

	for(t in 2:length(mature))
	{
		tInitial <- tpts[t-1]
		tFinal <- tpts[t]

		prematureInitial <- premature[t-1]
		prematureFinal <- premature[t]

		matureInitial <- matureEstimated[t-1]

		coefficients <- lineCoefficients_NoNascent(xInitial = tInitial
		                                	  ,xFinal = tFinal
		                                	  ,yInitial = prematureInitial
		                                	  ,yFinal = prematureFinal)
		m <- coefficients[1]
		q <- coefficients[2]

    	matureEstimated[t] <- (exp(-k3*(tFinal-tInitial))*matureInitial
                        	+ exp(-k3*tFinal)*k2K3Ratio*(q*(exp(k3*tFinal)-exp(k3*tInitial))
                        	+ m*(exp(k3*tFinal)*(k3*tFinal-1)/k3
                        	- exp(k3*tInitial)*(k3*tInitial-1)/k3)))
  	}

   return(sum((mature[-1] - matureEstimated[-1])^2/matureVariance[-1]))
}

secondStepError_NoNascent <- function(tpts
								 ,mature
								 ,premature
								 ,matureVariance
								 ,k2k3)
{
	k2 <- k2k3[1]
	k3 <- k2k3[2]
	
	matureEstimated <- numeric(length = length(mature))
	matureEstimated[1] <- mature[1]
  
	for(t in 2:length(mature))
	{
		tInitial <- tpts[t-1]
		tFinal <- tpts[t]

		prematureInitial <- premature[t-1]
		prematureFinal <- premature[t]

		matureInitial <- matureEstimated[t-1]

		coefficients <- lineCoefficients_NoNascent(xInitial = tInitial
											  ,xFinal = tFinal
											  ,yInitial = prematureInitial
											  ,yFinal = prematureFinal)
		m <- coefficients[1]
		q <- coefficients[2]

		matureEstimated[t] <- (exp(-k3*(tFinal-tInitial))*matureInitial
							 + exp(-k3*tFinal)*k2/k3*(q*(exp(k3*tFinal)-exp(k3*tInitial))
							 + m*(exp(k3*tFinal)*(k3*tFinal-1)/k3
							 - exp(k3*tInitial)*(k3*tInitial-1)/k3)))
	}

	return(mean((mature[-1] - matureEstimated[-1])^2/matureVariance[-1]))
}

RNAdynamics_NoNascent <- function(totRpkms
								, tpts
								, BPPARAM=SerialParam()
								, modellingParameters=list(Dmin = 1e-6, Dmax = 10)
								, genesFilter
								, imputeNAs
								)
{
	
	Dmin <- modellingParameters$Dmin
	Dmax <- modellingParameters$Dmax

	eiGenes <- rownames(totRpkms$exons)

	# Mature, premature and total rpkms
	mature <- totRpkms$exons - totRpkms$introns
	premature <- totRpkms$introns
	total <- totRpkms$exons

	# Mature, premature and total variance, rpkms
	matureVariance <- totRpkms$introns_var+totRpkms$exons_var
	prematureVariance <- totRpkms$introns_var
	totalVariance <- totRpkms$exons_var

	# Constant post transcriptional rates and fixed post transcriptional ratio 
	suppressWarnings(k3Prior <- firstStep_NoNascent(tpts = tpts
							  ,mature = mature
							  ,premature = premature
							  ,matureVariance = matureVariance
							  ,Dmin = Dmin
							  ,Dmax = Dmax))
	rownames(k3Prior) <- eiGenes

	# Constant post transcriptional rates and variable post transcriptiona ratio
	# (if there is only one gene, constraint for positive rates, i.e. inspectFromPCR)
	if( nrow(mature)==1 ) {
		row <- 1
		fits <- optimPositive(par = c(mature[row,1]/premature[row,1]*k3Prior[row,'k3'], k3Prior[row,'k3'])
					,fn = secondStepError_NoNascent
					,tpts = tpts
					,premature = premature[row,]
					,mature = mature[row,]
					,matureVariance = matureVariance[row,])
		fits <- t(unlist(fits[1:4]))
	} else if(imputeNAs){
		fits <- t(mcsapply(1:nrow(mature), function(row)
		{
			unlist(
				tryCatch(
					optim(par = c(mature[row,1]/premature[row,1]*k3Prior[row,'k3'], k3Prior[row,'k3'])
								,fn = secondStepError_NoNascent
								,tpts = tpts
								,premature = premature[row,]
								,mature = mature[row,]
								,matureVariance = matureVariance[row,])[c('par','value','convergence')]
					,error=function(e)list(par = c(NaN,NaN), value = NaN, convergence = NaN)))
		},BPPARAM = BPPARAM))
	} else {
		fits <- t(mcsapply(1:nrow(mature), function(row)
		{
			unlist(
				tryCatch(
					optimPositive(par = c(mature[row,1]/premature[row,1]*k3Prior[row,'k3'], k3Prior[row,'k3'])
								,fn = secondStepError_NoNascent
								,tpts = tpts
								,premature = premature[row,]
								,mature = mature[row,]
								,matureVariance = matureVariance[row,])[c('par','value','convergence')]
					,error=function(e)list(par = c(NaN,NaN), value = NaN, convergence = NaN)))
		},BPPARAM = BPPARAM))
	}

	fits[,3] <- pchisq(fits[,3], length(tpts)-3)
	colnames(fits) <- c('k2','k3','p','convergence')
	rownames(fits) <- eiGenes

	# Correction of negative priors with the median
	if(genesFilter){
		fits[fits[,'k2']<0,'k2'] <- NaN
		fits[fits[,'k3']<0,'k3'] <- NaN

		notFiniteRates <- !is.finite(fits[,'k2']) | !is.finite(fits[,'k3'])

		fits[notFiniteRates,'k2'] <- median(fits[is.finite(fits[,'k2']),'k2'])
		fits[notFiniteRates,'k3'] <- median(fits[is.finite(fits[,'k3']),'k3'])

		fits[notFiniteRates,'p'] <- NaN
	}

	# Data formatting
	constantModels <- list(models = fits
						 , premature = premature
						 , mature = mature
						 , prematureVariance = prematureVariance
						 , matureVariance = matureVariance)

	ratesConstantPriors <- constantModels$models  

	# betaTC <- matrix(rep(ratesConstantPriors[,'k3'],length(tpts)),ncol=length(tpts))
	gammaTC <- matrix(rep(ratesConstantPriors[,'k2'],length(tpts)),ncol=length(tpts))

	prematureDer <- as.matrix(t(sapply(1:nrow(premature),function(i){
		if(all(is.finite(premature[i,]))){
			spfun <- splinefun(tpts, premature[i,])
			return(spfun(tpts, deriv=1))
		} else return(rep(NA,length(tpts)))
	})))
	prematureDer[,1] <- 0 # force the steady state at time 0

	alphaTC <- prematureDer + gammaTC * premature
	alphaTC[alphaTC<0] <- NaN

	if(!imputeNAs)
	{
		message(paste0(table(apply(alphaTC,1,function(r)any(!is.finite(r))))['TRUE']),' genes were removed because of negative syntesis rate.')
		eiGenes <- eiGenes[which(apply(alphaTC,1,function(r)all(is.finite(r))))]

		total <- total[eiGenes,]
		premature <- premature[eiGenes,]
		mature <- mature[eiGenes,]
		totalVariance <- totalVariance[eiGenes,]
		prematureVariance <- prematureVariance[eiGenes,]
		matureVariance <- matureVariance[eiGenes,]

		alphaTC <- alphaTC[eiGenes,]
	}

	#Evaluate beta as constant between intervals

	betaT0 <- alphaTC[,1] / mature[,1]
	betaT0[betaT0 < 0 | !is.finite(betaT0)] <- NA
	
	betaOut <- inferKBetaFromIntegralWithPre(tpts, alphaTC, total, premature, 
				maxBeta=quantile(betaT0,na.rm=TRUE,probs=.99)*10,BPPARAM=BPPARAM
				)
	if( nrow(mature)==1 ) {
		betaTC <- t(c(betaT0, sapply(betaOut, function(x) sapply(x, '[[', 'root'))))
		betaEstimPrec <- t(c(0,sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec'))))
	} else {
		betaTC <- cbind(betaT0, sapply(betaOut, function(x) sapply(x, '[[', 'root')))
		betaEstimPrec <- cbind(0,sapply(betaOut, function(x) sapply(x, '[[', 'estim.prec')))
	}

	#Evaluate gamma as constant between intervals

	gammaT0 <- alphaTC[,1] / premature[,1]
	gammaT0[gammaT0 < 0 | !is.finite(gammaT0)] <- NA

	gammaOut <- inferKGammaFromIntegral(tpts, alphaTC, premature, 
		maxGamma=quantile(gammaT0,na.rm=TRUE,probs=.99)*10, BPPARAM=BPPARAM
		)
	if( nrow(mature)==1 ) {
		gammaTC <- t(c(gammaT0, sapply(gammaOut, function(x) sapply(x, '[[', 'root'))))
		gammaEstimPrec <- t(c(0, sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec'))))
	} else {
		gammaTC <- cbind(gammaT0, sapply(gammaOut, function(x) sapply(x, '[[', 'root')))
		gammaEstimPrec <- cbind(0, sapply(gammaOut, function(x) sapply(x, '[[', 'estim.prec')))
	}
	# ## impute NA values
	alphaTC <- do.call('rbind',bplapply(1:nrow(alphaTC), 
		function(i) impute_na_tc(tpts, alphaTC[i,]), BPPARAM=BPPARAM))
	betaTC <- do.call('rbind',bplapply(1:nrow(betaTC), 
		function(i) impute_na_tc(tpts, betaTC[i,]), BPPARAM=BPPARAM))
	gammaTC <- do.call('rbind',bplapply(1:nrow(gammaTC), 
		function(i) impute_na_tc(tpts, gammaTC[i,]), BPPARAM=BPPARAM))

	## caluculate error through integration of alphaTC, betaTC, gammaTC?

	pModel <- fits[,"p"]
	pModel[apply(alphaTC,1,function(row)any(!is.finite(row)))] <- NaN

	alphaTC_var <- matrix(NaN, nrow(alphaTC), ncol(alphaTC))
	ratesEstimPrec <- betaEstimPrec + gammaEstimPrec


	attr(alphaTC, 'dimnames') <- NULL
	attr(betaTC, 'dimnames') <- NULL
	attr(gammaTC, 'dimnames') <- NULL
	attr(total, 'dimnames') <- NULL
	attr(premature, 'dimnames') <- NULL
	attr(alphaTC_var, 'dimnames') <- NULL
	attr(totalVariance, 'dimnames') <- NULL
	attr(prematureVariance, 'dimnames') <- NULL
	attr(ratesEstimPrec, 'dimnames') <- NULL

	return(list(
		concentrations=list(
			total=total
			, total_var=totalVariance
			, preMRNA=premature
			, preMRNA_var=prematureVariance
			, labeled_total=matrix(NA, nrow=nrow(total), ncol=ncol(total))
			, labeled_total_var=matrix(NA, nrow=nrow(total), ncol=ncol(total))
			, labeled_preMRNA=matrix(NA, nrow=nrow(total), ncol=ncol(total))
			, labeled_preMRNA_var=matrix(NA, nrow=nrow(total), ncol=ncol(total))
		)
		, rates=list(
			alpha=alphaTC
			, alpha_var=alphaTC_var
			, beta=betaTC
			, gamma=gammaTC)
		, ratesEstimPrec=ratesEstimPrec
		, ratesFirstGuessP = pModel
		, geneNames=eiGenes
		, tpts=tpts
		, tL=NaN))

}
