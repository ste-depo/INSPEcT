#' @rdname viewConfidenceIntervals
#'
#' @description
#' A method to access the modeld confidence intervals via the method \code{\link{confidenceIntervals}}
#' @param object An object of class INSPEcT
#' @param feature A character indicating the feature to retireve: "synthesis", "degradation", "processing".
#' @return A numeric matrix containing the values for the selected feature
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' viewConfidenceIntervals(nascentInspObj10, 'synthesis')
setMethod('viewConfidenceIntervals', 'INSPEcT', function(object, feature) {

	if(object@NoNascent)stop('viewConfidenceIntervals: no confidence intervals for INSPEcT objects without nascent RNA.')

	ix <- grep(feature,pData(object@confidenceIntervals)$feature)
	exprs(object@confidenceIntervals)[,ix, drop=FALSE]
})