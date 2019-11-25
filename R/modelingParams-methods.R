#' @rdname modelingParams
#'
#' @description
#' A method to get the parameters used for modeling rates and 
#' concentrations by the method \code{\link{modelRates}}
#' @param object An object of class INSPEcT
#' @return List of parameters and their values
#' \itemize{
#' \item estimateRatesWith Either "int" or "der". With "int" the degradation and processing
#'    rates are estimated integrating the system between one time point and the following. 
#'    With "der" degradation and processing rates are estimated using the derivative of total
#'    and pre mRNA. 
#' \item useSigmoidFun A logical, whether to choose between sigmoid and impulse function 
#'    to fit rates and concentrations. In case not, always impulse function is used. 
#' \item testOnSmooth A logical, wheter models should be tested on smoothed pre-mRNA, 
#'    total mRNA and eventually synthesis rates or not.
#' \item nInit number of optimization to find the best functional representation of each rate
#' \item nIter number of max iteration during optimization
#' \item Dmin lower bondary for degradation rates in the NoNascent mode
#' \item Dmax upper bondary for degradation rates in the NoNascent mode
#' \item seed A numeric, indicatindg the seed set for reproducible results.
#' }
#' @seealso \code{\link{modelRates}}
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' modelingParams(nascentInspObj10)
setMethod('modelingParams', 'INSPEcT', function(object) {
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	return(object@params)
	})