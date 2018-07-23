#' Given introns and exons abundances (for example RPKMs) this method returns their variances evaluated thorugh plgem.
#' @param exonsAbundances A matrix containing the exons abundances for each experimental condition.
#' @param intronsAbundances A matrix containing the intorns abundances for each experimental condition.
#' @param experimentalDesign A numerical which reports the desing of the experiment in terms of time points and replicates. The time points must be ordered according
#' to the columns of the count matrices submitted for the analysis; these labels define conditions and replicates.
#' @param varSamplingCondition A character reporting which experimental condition should be used to sample the variance if DESeq2 = FALSE.
#' @param simulatedData A boolean which is TRUE if the data under analysis are simulated.
#' @return A list containing RPKMs and associated variances for exons and introns.

quantifyExpressionsFromTrAbundance <- function(exonsAbundances
										     , intronsAbundances
											 , experimentalDesign
											 , varSamplingCondition
											 , simulatedData = FALSE)
{

	if(table(experimentalDesign)[varSamplingCondition]==1)
		stop('quantifyExpressions: the varSamplingCondition must have at least one replicate.')
	if(!simulatedData)
	{
		if(ncol(exonsAbundances)!=ncol(intronsAbundances))
			stop('makeExpressions: dimensionality issue in input data.')
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
	}

	return(list(exonsExpressions = expressionsExons
			  , intronsExpressions = NULL
			  , exonsVariance = varianceExpressionsExons
			  , intronsVariance = NULL))

}