#' @rdname viewConfidenceIntervals
#'
#' @description
#' A method to access the modeld confidence intervals computed via the method \code{\link{computeConfidenceIntervals}}
#' @param object An object of class INSPEcT
#' @param feature A character indicating the feature to retireve: "synthesis", "degradation", "processing".
#' @return A numeric matrix containing the values for the selected feature
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' viewConfidenceIntervals(nascentInspObj10, 'synthesis')
setMethod('viewConfidenceIntervals', 'INSPEcT', function(object, feature) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	ix <- grep(feature,pData(object@confidenceIntervals)$feature)
	out <- exprs(object@confidenceIntervals)[,ix, drop=FALSE]
	if(all(dim(out)==c(0,0)))
	{
		message("Confidence intervals must be computed!")
		return(NaN)
	}else{return(out)}
})