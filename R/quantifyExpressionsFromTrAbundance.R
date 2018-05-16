quantifyExpressionsFromTrAbundance <- function(exonsAbundances
										     , intronsAbundances
											 , experimentalDesign
											 , varSamplingCondition
											 , plgemFits = NULL
											 , returnPlgemFits = FALSE)
{

	if(table(experimentalDesign)[varSamplingCondition]==1 & is.null(plgemFits))
		stop('quantifyExpressions: the varSamplingCondition must have at least one replicate.')
	if(ncol(exonsAbundances)!=ncol(intronsAbundances))
		stop('makeExpressions: dimensionality issue in input data.')

	expressionsExons <- exonsAbundances
	expressionsIntrons <- intronsAbundances
	
	if(is.null(plgemFits))
	{
		print('Powerlaw fit for exons')
		
		exprData <- expressionsExons
		pData <- data.frame(feature=experimentalDesign)
		colnames(exprData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
		rownames(pData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
		phenoData <- new('AnnotatedDataFrame', data=pData)
		fData <- data.frame(exprData)
		featureData <- new('AnnotatedDataFrame', data=fData)
		foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
		exonsPlgemFit <- plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,plot.file=FALSE,fittingEval=FALSE, verbose=FALSE)
		exonsPlgemLaw <- function(expression){(exp(exonsPlgemFit$INTERCEPT)*expression^(exonsPlgemFit$SLOPE))^2}
	
		print('Powerlaw fit for introns')
	
		exprData <- intronsAbundances
		pData <- data.frame(feature=experimentalDesign)
		colnames(exprData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
		rownames(pData) <- paste0(pData$feature,"_",seq_along(experimentalDesign))
		phenoData <- new('AnnotatedDataFrame', data=pData)
		fData <- data.frame(exprData)
		featureData <- new('AnnotatedDataFrame', data=fData)
		foe <- ExpressionSet(assayData=exprData, phenoData=phenoData, featureData=featureData)
		intronsPlgemFit <- plgem.fit(data=foe,fitCondition=varSamplingCondition,p=10,q=.5,plot.file=FALSE,fittingEval=FALSE, verbose=FALSE)
		intronsPlgemLaw <- function(expression){(exp(intronsPlgemFit$INTERCEPT)*expression^(intronsPlgemFit$SLOPE))^2}

	}else{
		exonsPlgemFit <- plgemFits$exonsPlgemFit
		exonsPlgemLaw <- function(expression){(exp(exonsPlgemFit$INTERCEPT)*expression^(exonsPlgemFit$SLOPE))^2}

		intronsPlgemFit <- plgemFits$intronsPlgemFit
		intronsPlgemLaw <- function(expression){(exp(intronsPlgemFit$INTERCEPT)*expression^(intronsPlgemFit$SLOPE))^2}		
	}

	varianceExpressionsExons <- t(sapply(1:nrow(expressionsExons),function(r)
	{
		sapply(1:ncol(expressionsExons),function(c)
		{
			exonsPlgemLaw(expressionsExons[r,c])
		})
	}))

	varianceExpressionsIntrons <- t(sapply(1:nrow(expressionsIntrons),function(r)
	{
		sapply(1:ncol(expressionsIntrons),function(c)
		{
			intronsPlgemLaw(expressionsIntrons[r,c])
		})
	}))

	expressionsExons <- t(apply(expressionsExons,1,function(x)tapply(x, experimentalDesign, mean)))
	expressionsIntrons <- t(apply(expressionsIntrons,1,function(x)tapply(x, experimentalDesign, mean)))
	varianceExpressionsExons <- t(apply(varianceExpressionsExons,1,function(x)tapply(x, experimentalDesign, mean)))
	varianceExpressionsIntrons <- t(apply(varianceExpressionsIntrons,1,function(x)tapply(x, experimentalDesign, mean)))

	rownames(varianceExpressionsExons) <- rownames(expressionsExons)
	rownames(varianceExpressionsIntrons) <- rownames(expressionsIntrons)

	if(returnPlgemFits)
	{
		return(list(exonsExpressions = expressionsExons
				  , intronsExpressions = expressionsIntrons
				  , exonsVariance = varianceExpressionsExons
				  , intronsVariance = varianceExpressionsIntrons
				  , exonsPlgemFit = exonsPlgemFit
				  , intronsPlgemFit = intronsPlgemFit))

	}else{
		return(list(exonsExpressions = expressionsExons
				  , intronsExpressions = expressionsIntrons
				  , exonsVariance = varianceExpressionsExons
				  , intronsVariance = varianceExpressionsIntrons))
	}

}