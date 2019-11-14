#' @rdname modelSelection
#'
#' @description
#' Method to visualize the criteria used to assess variability of rates.
#'
#' @param object An object of class INSPEcT or INSPEcT_model
#'
#' @return 
#' \itemize{
#' \item modelSelection 'aic' compares nested models closest to the one with lowest AIC, 'llr' compares all nested models, 
#' 'hib' is a mix between the previous two. (default 'aic')
#' \item preferPValue a logical, if TRUE (default) limit the search for best models among the ones with succeded the goodness of fit test.
#' \item padj a logical, if TRUE (default) correct the p-values for multiple testing
#' \item goodness_of_fit a numeric, the threshold for the goodness-of-fit test (default = .1)
#' \item variability a numeric, a vector with the thresholds for the variability test (one threshold for each rate, default = c('s'=05, 'p'=.05, 'd'=.05))
#' \item limitModelComplexity a logical that limits the complexity of the function used to describe dynamics to the length of the time-course (default = FALSE)
#'     
#' }
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' modelSelection(nascentInspObj10)
setMethod('modelSelection', 'INSPEcT', function(object) {
	return(modelSelection(object@model))
	})
#' @rdname modelSelection
setMethod('modelSelection', 'INSPEcT_model', function(object) {
	return(object@params)
	})
