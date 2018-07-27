#' Given introns and exons abundances (for example RPKMs) this method returns their variances evaluated thorugh plgem.
#' @param exonsAbundances A matrix containing the exons abundances for each experimental condition.
#' @param intronsAbundances A matrix containing the intorns abundances for each experimental condition.
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the columns of the count matrices submitted for the analysis; these labels define conditions and replicates.
#' @param varSamplingCondition A character reporting which experimental condition should be used to sample the variance if DESeq2 = FALSE.
#' @param simulatedData A boolean which is TRUE if the data under analysis are simulated.
#' @return A list containing RPKMs and associated variances for exons and introns.

quantifyExpressionsFromTrAbundance <- function(trAbundaces
											 , experimentalDesign
											 , varSamplingCondition = NULL
											 , simulatedData = FALSE)
{

	if( !is.logical(simulatedData) )
		stop('quantifyExpressionsFromTrAbundance: "simulatedData" must be a logical.')
	if( !simulatedData ) {
		# trAbundaces
		if( ! is.list(trAbundaces) )
			stop('quantifyExpressionsFromTrAbundance: "trAbundaces" must be a list with elements "exonsAbundances" and "intronsAbundances"')
		if( ! all(c('exonsAbundances','intronsAbundances') %in% names(trAbundaces)) )
			stop('quantifyExpressionsFromTrAbundance: "trAbundaces" must be a list with elements "exonsAbundances" and "intronsAbundances"')
		exonsAbundances <- trAbundaces$exonsAbundances
		intronsAbundances <- trAbundaces$intronsAbundances
		if( !( is.matrix(exonsAbundances) & is.matrix(intronsAbundances) ) )
			stop('quantifyExpressionsFromTrAbundance: the elements "exonsAbundances" and "intronsAbundances" of "trAbundaces" must be matrices with the same numebr of columns.')
		if( ncol(exonsAbundances) != ncol(intronsAbundances) )
			stop('quantifyExpressionsFromTrAbundance: the elements "exonsAbundances" and "intronsAbundances" of "trAbundaces" must be matrices with the same numebr of columns.')
		# experimentalDesign
		if(all(table(experimentalDesign)==1))
			stop("quantifyExpressionsFromTrAbundance: at least one condition with replicates is required.")
		if(length(experimentalDesign)!=ncol(exonsAbundances))
			stop('quantifyExpressionsFromTrAbundance: each counts column must be accounted in the experimentalDesign')
		# varSamplingCondition
		if( is.null(varSamplingCondition) ) {
			varSamplingCondition <- names(which(table(experimentalDesign)>1)[1])
		} else {
			if( length(which(as.character(experimentalDesign) == varSamplingCondition)) < 2 )
				stop('quantifyExpressionsFromTrAbundance: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
		}
	} else {
		# trAbundaces
		if( ! is.list(trAbundaces) )
			stop('quantifyExpressionsFromTrAbundance: "trAbundaces" must be a list with element "exonsAbundances".')
		if( ! 'exonsAbundances' %in% names(trAbundaces) )
			stop('quantifyExpressionsFromTrAbundance: "trAbundaces" must be a list with element "exonsAbundances".')
		exonsAbundances <- trAbundaces$exonsAbundances
		if( !is.matrix(exonsAbundances) )
			stop('quantifyExpressionsFromTrAbundance: the element "exonsAbundances" of "trAbundaces" must be matrix.')
		# experimentalDesign
		if(all(table(experimentalDesign)==1))
			stop("quantifyExpressionsFromTrAbundance: at least one condition with replicates is required.")
		if(length(experimentalDesign)!=ncol(exonsAbundances))
			stop('quantifyExpressionsFromTrAbundance: each counts column must be accounted in the experimentalDesign')
		# varSamplingCondition
		if( is.null(varSamplingCondition) ) {
			varSamplingCondition <- names(which(table(experimentalDesign)>1)[1])
		} else {
			if( length(which(as.character(experimentalDesign) == varSamplingCondition)) < 2 )
				stop('quantifyExpressionsFromTrAbundance: if DESeq2 is FALSE varSamplingCondition must be an experimental condition with replicates.')
		}		
	}

	expressionsExons <- exonsAbundances
	
	print('Powerlaw fit for exons')
	
	exprData <- expressionsExons
	pData <- data.frame(feature=experimentalDesign)
	colnames(exprData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
	rownames(pData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
	phenoData <- new('AnnotatedDataFrame', data=pData)
	fData <- data.frame(exprData)
	featureData <- new('AnnotatedDataFrame', data=fData)
	foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)

	exonsPlgemFit <- suppressWarnings(plgem.fit(data=foe
											  , fitCondition=varSamplingCondition
											  , p=10
											  , q=.5
											  , zeroMeanOrSD="trim"
											  , plot.file=FALSE
											  , fittingEval=FALSE
											  , verbose=FALSE))

	exonsPlgemLaw <- function(expression){(exp(exonsPlgemFit$INTERCEPT)*expression^(exonsPlgemFit$SLOPE))^2}

	varianceExpressionsExons <- t(sapply(1:nrow(expressionsExons),function(r)
	{
		sapply(1:ncol(expressionsExons),function(c)
		{
			exonsPlgemLaw(expressionsExons[r,c])
		})
	}))

	expressionsExons <- t(apply(expressionsExons,1,function(x)tapply(x, experimentalDesign, mean)))
	varianceExpressionsExons <- t(apply(varianceExpressionsExons,1,function(x)tapply(x, experimentalDesign, mean)))
	rownames(varianceExpressionsExons) <- rownames(expressionsExons)

	if(!simulatedData)
	{
		expressionsIntrons <- intronsAbundances
	
		print('Powerlaw fit for introns')
	
		exprData <- intronsAbundances
		pData <- data.frame(feature=experimentalDesign)
		colnames(exprData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
		rownames(pData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
		phenoData <- new('AnnotatedDataFrame', data=pData)
		fData <- data.frame(exprData)
		featureData <- new('AnnotatedDataFrame', data=fData)
		foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
		intronsPlgemFit <- suppressWarnings(plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,zeroMeanOrSD="trim",plot.file=FALSE,fittingEval=FALSE, verbose=FALSE))
		intronsPlgemLaw <- function(expression){(exp(intronsPlgemFit$INTERCEPT)*expression^(intronsPlgemFit$SLOPE))^2}		

		varianceExpressionsIntrons <- t(sapply(1:nrow(expressionsIntrons),function(r)
		{
			sapply(1:ncol(expressionsIntrons),function(c)
			{
				intronsPlgemLaw(expressionsIntrons[r,c])
			})
		}))

		expressionsIntrons <- t(apply(expressionsIntrons,1,function(x)tapply(x, experimentalDesign, mean)))
		varianceExpressionsIntrons <- t(apply(varianceExpressionsIntrons,1,function(x)tapply(x, experimentalDesign, mean)))
		rownames(varianceExpressionsIntrons) <- rownames(expressionsIntrons)

		return(list(exonsExpressions = expressionsExons
				  , intronsExpressions = expressionsIntrons
				  , exonsVariance = varianceExpressionsExons
				  , intronsVariance = varianceExpressionsIntrons))
	} else {

		return(list(exonsExpressions = expressionsExons
				  , intronsExpressions = NULL
				  , exonsVariance = varianceExpressionsExons
				  , intronsVariance = NULL))

	}

}