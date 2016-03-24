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
#' data('rpkms', package='INSPEcT')
#' tpts <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
#' tL <- 1/6
#' mycerIds <- newINSPEcT(tpts, tL, rpkms$foursu_exons, rpkms$total_exons, 
#' 	rpkms$foursu_introns, rpkms$total_introns)
#' # get estimated synthesis rates
#' ratesFirstGuess(mycerIds, 'synthesis')
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
#' data('rpkms', package='INSPEcT')
#' tpts <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
#' tL <- 1/6
#' mycerIds <- newINSPEcT(tpts, tL, rpkms$foursu_exons, rpkms$total_exons, 
#' 	rpkms$foursu_introns, rpkms$total_introns)
#' ratesFirstGuessVar(mycerIds, 'synthesis')
setMethod('ratesFirstGuessVar', 'INSPEcT', function(object, feature) {
	fData(object@ratesFirstGuess)[[feature]]
	})
