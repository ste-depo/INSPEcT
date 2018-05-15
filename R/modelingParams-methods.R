#' @rdname modelingParams
#'
#' @description
#' A method to get and set the parameters that will be used in the modeling of estimated rates and 
#' concentrations by the method \code{\link{modelRates}}
#' @param object An object of class INSPEcT
#' @return List of parameters and their values
#' \itemize{
#'  \item nInit number of optimization to find the best functional representation of each 
#'    rate (by default 10)
#'  \item nIter number of max iteration during optimization (default is 300)
#'  \item na.rm A logical wheter missing values should be removed from estimated rates 
#'    (default is TRUE)
# #'  \item nCores The number of cores to be used during parallelization (default is 2)
#'  \item verbose A logical wheter to be verbose or not (default is TRUE)
#'  \item estimateRatesWith Either "int" or "der". With "int" the degradation and processing
#'    rates are estimated integrating the system between one time point and the following. 
#'    With "der" degradation and processing rates are estimated using the derivative of total
#'    and pre mRNA. (default is "int")
#'  \item useSigmoidFun A logical, whether to choose between sigmoid and impulse function 
#'    to fit rates and concentrations. In case not, always impulse function is used. 
#'    (default is TRUE)
#'  \item testOnSmooth A logical, wheter models should be tested on smoothed pre-mRNA, 
#'    total mRNA and synthesis rates or not. (default is TRUE)
#' }
#' @seealso \code{\link{modelRates}}
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' modelingParams(mycerIds10)
setMethod('modelingParams', 'INSPEcT', function(object) {
	return(object@params)
	})
#' @rdname modelingParams
#'
#' @param value A list with new parameters
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' mycerIds10 <- removeModel(mycerIds10)
#' modelingParams(mycerIds10)$useSigmoidFun <- FALSE
setReplaceMethod('modelingParams', 'INSPEcT', function(object, value) {

	if(!is.list(value))
		stop('modelingParams: value argument must be a list')
	if( !is.logical(value$NoNascent) )
		stop('modelingParams: NoNascent element of value argument must be a logical')
	if( !is.numeric(value$nInit) )
		stop('modelingParams: nInit element of value argument must be a numeric')
	if( !is.numeric(value$nIter) )
		stop('modelingParams: nIter element of value argument must be a numeric')		
	if( !is.logical(value$'na.rm') )
		stop('modelingParams: na.rm element of value argument must be a logical')		
	if( !is.logical(value$verbose) )
		stop('modelingParams: verbose element of value argument must be a logical')		
	if( !is.character(value$estimateRatesWith) )
		stop('modelingParams: estimateRatesWith element of value argument must be a character')
	if( !is.logical(value$useSigmoidFun) )
		stop('modelingParams: useSigmoidFun element of value argument must be a logical')		
	if( !is.logical(value$testOnSmooth) )
		stop('modelingParams: testOnSmooth element of value argument must be a logical')
	if( !all(sapply(value, length)==1) )
		stop('modelingParams: all elements of value argument must have length 1')
	if( !value$estimateRatesWith %in% c('int', 'der') )
		stop('modelingParams: estimateRatesWith element of value argument must either "int" or "der"')
	if( length(object@model@ratesSpecs) > 0 )
		stop('Remove the model before changing modeling parameters. (See "?removeModel")')
	value$nInit <- as.integer(value$nInit)
	value$nIter <- as.integer(value$nIter)
	# value$nCores <- as.integer(value$nCores)
	object@params <- value
	object
	})
