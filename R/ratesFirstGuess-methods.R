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
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' 
#' ratesFirstGuess(nascentInspObj10, 'total')
#' ratesFirstGuess(nascentInspObj10, 'preMRNA')
#' ratesFirstGuess(nascentInspObj10, 'synthesis')

setMethod('ratesFirstGuess', 'INSPEcT', function(object, feature) {
	ix <- pData(object@ratesFirstGuess)$feature == feature
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
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' 
#' ratesFirstGuessVar(nascentInspObj10, 'total')
#' ratesFirstGuessVar(nascentInspObj10, 'preMRNA')
#' ratesFirstGuessVar(nascentInspObj10, 'synthesis')
setMethod('ratesFirstGuessVar', 'INSPEcT', function(object, feature) {
	ix <- pData(object@ratesFirstGuessVar)$feature == feature
	exprs(object@ratesFirstGuessVar)[,ix, drop=FALSE]
	})