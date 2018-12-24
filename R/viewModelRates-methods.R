#' @rdname viewModelRates
#'
#' @description
#' A method to access the modeld rates via the method \code{\link{modelRates}}
#' @param object An object of class INSPEcT
#' @param feature A character indicating the feature to retireve, "synthesis", "degradation", "processing" for rates, "total" for total mRNA concentrations or "preMRNA" for premature mRNA concentrations
#' @return A numeric matrix containing the values for the selected feature
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' viewModelRates(nascentInspObj10, 'synthesis')
setMethod('viewModelRates', 'INSPEcT', function(object, feature) {

	ix <- grep(feature,pData(object@modelRates)$feature)
	exprs(object@modelRates)[,ix, drop=FALSE]

	})
