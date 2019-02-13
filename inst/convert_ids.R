convert_ids = function(ids, degDuringPulse=FALSE) {

	## rearrange ratesFirstGuess and ratesFirstGuessVar
	fgVar= fData(ids@ratesFirstGuess)
	fgVar= fgVar[, grep('t0',colnames(fgVar),invert=TRUE)]
	totVar= sapply(seq_along(ids@tpts), function(i) fgVar[,grep('total', colnames(fgVar))])
	preVar= sapply(seq_along(ids@tpts), function(i) fgVar[,grep('preMRNA', colnames(fgVar))])
	synVar= sapply(seq_along(ids@tpts), function(i) fgVar[,grep('synthesis', colnames(fgVar))])
	proVar <- degVar <- matrix(1, nrow(synVar), ncol(synVar))
	varData = cbind(totVar, preVar, synVar, proVar, degVar)
	pData <- phenoData(ids@ratesFirstGuess)
	colnames(varData) = rownames(pData@data)
	rownames(varData) = featureNames(ids@ratesFirstGuess)
	ids@ratesFirstGuess <- ExpressionSet(assayData=exprs(ids@ratesFirstGuess), phenoData=pData)
	ids@ratesFirstGuessVar <- ExpressionSet(assayData=varData, phenoData=pData)

	## add new slots
	ids@NoNascent <- FALSE
	ids@degDuringPulse <- degDuringPulse[1]
	ids@model@params$preferPValue <- FALSE
	ids@model@params$padj <- FALSE

	## re-order brown thresholds
	ids@model@params$thresholds$brown <- 
		ids@model@params$thresholds$brown[c('synthesis','processing','degradation')]

	return(ids)

}