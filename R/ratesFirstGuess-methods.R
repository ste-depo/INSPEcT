#' @rdname ratesFirstGuess
#'
#' @description
#' This method allow to access to the estimated synthesis, degradation, processing rates and pre mRNA and total mRNA 
#' concentrations the way they were calculated by the constructor function \code{\link{newINSPEcT}}.
#' @param object An object of class INSPEcT
#' @param feature A character indicating the feature to retireve, "synthesis", "degradation", "processing" for rates, "total" for total mRNA concentrations or "preMRNA" for premature mRNA concentrations
#' @return A numeric matrix containing the values for the selected feature
#' @seealso \code{\link{newINSPEcT}}, \code{\link{ratesFirstGuessVar}}
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' 
#' ratesFirstGuess(mycerIds10, 'total')
#' ratesFirstGuess(mycerIds10, 'preMRNA')
#' ratesFirstGuess(mycerIds10, 'synthesis')
#' ratesFirstGuess(mycerIds10, 'processing')
#' ratesFirstGuess(mycerIds10, 'degradation')

setMethod('ratesFirstGuess', 'INSPEcT', function(object, feature) {
	ix <- grep(feature,pData(object@ratesFirstGuess)$feature)
	exprs(object@ratesFirstGuess)[,ix, drop=FALSE]
	})

#' @rdname ratesFirstGuessVar
#'
#' @description
#' This method allow to access to the estimated variance of synthesis rates and pre mRNA and total mRNA 
#' concentrations the way they were calculated by the constructor function \code{\link{newINSPEcT}}.
#' @param object An object of class INSPEcT
#' @param feature A character indicating the feature to retireve, "synthesis", "degradation", "processing" for rates, "total" for total mRNA concentrations or "preMRNA" for premature mRNA concentrations
#' @return A numeric vector containing the values for the selected feature
#' @seealso \code{\link{newINSPEcT}}, \code{\link{ratesFirstGuess}}
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' 
#' ratesFirstGuessVar(mycerIds10, 'total')
#' ratesFirstGuessVar(mycerIds10, 'preMRNA')
#' ratesFirstGuessVar(mycerIds10, 'synthesis')
#' ratesFirstGuessVar(mycerIds10, 'processing')
#' ratesFirstGuessVar(mycerIds10, 'degradation')
setMethod('ratesFirstGuessVar', 'INSPEcT', function(object, feature) {
	temp <- object@ratesFirstGuess@featureData@data[,grep("_t0",grep(feature,names(object@ratesFirstGuess@featureData@data),value=T),invert=T,value=T)]
	if(class(temp)=="numeric")
	{
		temp <- data.frame(temp)
		rownames(temp) <- rownames(object@ratesFirstGuess@featureData@data)
		colnames(temp) <- feature
	}
	as.matrix(temp)
	})
