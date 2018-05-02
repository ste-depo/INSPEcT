#' Get or set parameters for model test and selection
#'
#' @description
#' With this methods the user can personalize the criteria by which INSPEcT
#' selects a rate to be variable or constant. In particular, the model selection
#' criteria can be selected between log-likelihood ratio test and Akaike's information
#' criterion. In case log-likelihood ratio test is selected, the thresholds of 
#' chi-squared and Brown's method can be set (see Details section).
#'
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param value A list or a character that will substitute the set of parameters
#' \itemize{
#'   \item modelSelection: A character, either "llr" to test whether a rate is 
#'     varying using log-likelihood testing framework or "aic" to choose
#'     the best model via Akaike Information Criterion (Default: "llr").
#'   \item thresholds: A named list containing the threshold that is used to 
#'     consider a model as accepted in terms of the chi-squared test and three 
#'     thresholds (one per each rate) that are used to consider a rate as variable
#'     using the Brown's method on the log-likelihood ratio tests
#'   \item llrtests: A list of three elements that represent, for each rate, the pairs 
#'     of models that will be compared via log likelihood ratio test to assess 
#'     whether the rate is variable or not
#' }
#'
#' @return See "value"
#'
#' @details
#' When log-likelihood is chosen as a criterion for model selection, different nested 
#' models can be compared to assess wheter a single rate is varying or constant. 
#' For example, in case we want to establish whether synthesis rate is constant or not
#' we can test the null hypothesis "all the rates are constant" against the alternative 
#' hypothesis "synthesis rate is changing". The null hypothesis is a special case
#' of the alternative hypothesis, therefore the models are nested. We can also assess
#' whether synthesis rate is constant or not by comparing the null hypothesis 
#' "degradation rate is changing" against the alternative hypothesis "degradation and
#' synthesis are changing". The method \code{\link{llrtests}} set the models 
#' that are compared to assess the variability of eache rate. Different comparisons will
#' be combined using Brown's method for combinig p-values.
#' Models are named with a short notation where synthesis is "a", degradation is "b"
#' and processing is "c". "0" is the model where all genes are kept constant
#' and "ab", for example is the model where synthesis rate and degradation rate 
#' are changing.
#' The user can also set the thresholds for Brown's p-value and chi-suqared p-value.
#' While the former set the threshold to assess whether a rate is variable or not over time,
#' the latter set the chi-squared threshold for a pair of model to be used via the
#' log-likelihood ratio test. In order for a pair to be used, at least one model of the 
#' pair should have a chi-squared p-value (goodness of fit) below the threshold.
#' The construction of a synthetic data-set can help in the choice of the correct 
#' parameters for the test (\code{\link{makeSimModel}}, \code{\link{makeSimDataset}}).
#' @seealso \code{\link{makeSimModel}}, \code{\link{makeSimDataset}}
#'
#' @name testingParams
NULL

#' @rdname testingParams
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' modelSelection(mycerIds10)
#' modelSelection(mycerIds10) <- 'aic'
setMethod('modelSelection', 'INSPEcT', function(object) {
	return(object@model@params$modelSelection)
	})
#' @rdname testingParams
setReplaceMethod('modelSelection', 'INSPEcT', function(object, value) {
	if( !is.character(value) )
		stop('modelSelection: value argument must be a character')
	if( !value %in% c('aic', 'llr') )
		stop('modelSelection: value argument must either "llr" or "aic"')
	# if( value == 'aic' ) {
	# 	warning('modelSelection: aic option is still not active')
	# 	return(object)
	# }
	if( !identical(value, object@model@params$modelSelection) ) {
		object@model@params$modelSelection <- value
		# if modeled rates already exist, update them
		if( length(object@model@ratesSpecs)>1 ) {
			message('Updating modeled rates...')
			object <- makeModelRates(object)
		}
	}
	return(object)
	})
#' @rdname testingParams
setMethod('modelSelection', 'INSPEcT_model', function(object) {
	return(object@params$modelSelection)
	})
#' @rdname testingParams
setReplaceMethod('modelSelection', 'INSPEcT_model', function(object, value) {
	if( !is.character(value) )
		stop('modelSelection: value argument must be a character')
	if( !value %in% c('aic', 'llr') )
		stop('modelSelection: value argument must either "llr" or "aic"')
	object@params$modelSelection <- value
	return(object)
	})

#' @rdname testingParams
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' thresholds(mycerIds10)
#' thresholds(mycerIds10)$chisquare <- 1e-3
#' thresholds(mycerIds10)$brown['synthesis'] <- 1e-3
setMethod('thresholds', 'INSPEcT', function(object) {
	return(object@model@params$thresholds)
	})
#' @rdname testingParams
setReplaceMethod('thresholds', 'INSPEcT', function(object, value) {
	if( !is.list(value) )
		stop('thresholds: value argument must be a list')
	if( !identical(names(value),c('chisquare', 'brown')) )
		stop('thresholds: value argument list must be named. Names must be "chisquare" and "brown"')
	if( !is.numeric(value$chisquare) )
		stop('thresholds: chisquare element of value argument must be a numeric')
	if( !is.numeric(value$brown) )
		stop('thresholds: brown element of value argument must be a numeric')
	if( length(value$chisquare) != 1 )
		stop('thresholds: chisquare element of value argument must be of length 1')
	if( length(value$brown) != 3 )
		stop('thresholds: brown element of value argument must be of length 3')
	if( !identical(value, object@model@params$thresholds) ) {
		object@model@params$thresholds <- value
		# if modeled rates already exist, update them
		if( length(object@model@ratesSpecs)>1 ) {
			message('Updating modeled rates...')
			object <- makeModelRates(object)
		}
	}
	return(object)
	})
#' @rdname testingParams
setMethod('thresholds', 'INSPEcT_model', function(object) {
	return(object@params$thresholds)
	})
#' @rdname testingParams
setReplaceMethod('thresholds', 'INSPEcT_model', function(object, value) {
	if( !is.list(value) )
		stop('thresholds: value argument must be a list')
	if( !identical(names(value), c('chisquare', 'brown')) )
		stop('thresholds: value argument must be named list. Names must be "chisquare" and "brown"')
	if( !is.numeric(value$chisquare) )
		stop('thresholds: chisquare element of value argument must be a numeric')
	if( !is.numeric(value$brown) )
		stop('thresholds: brown element of value argument must be a numeric')
	if( length(value$chisquare) != 1 )
		stop('thresholds: chisquare element of value argument must be of length 1')
	if( length(value$brown) != 3 )
		stop('thresholds: brown element of value argument must be of length 3')
	object@params$thresholds <- value
	return(object)
	})

#' @rdname testingParams
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' llrtests(mycerIds10)
#' llrtests(mycerIds10)$synthesis <- list(c('0','a'), c('b','ab'))
setMethod('llrtests', 'INSPEcT', function(object) {
	return(object@model@params$llrtests)
	})
#' @rdname testingParams
setReplaceMethod('llrtests', 'INSPEcT', function(object, value) {
	if( !is.list(value) )
		stop('llrtests: value argument must be a list')
	if( !identical(names(value), c('synthesis', 'degradation', 'processing')) )
		stop('llrtests: value argument must be a named list. Names must be "synthesis", "degradation" and "processing"')
	if( !is.list(value$synthesis) )
		stop('llrtests: synthesis element of value argument must be a list')
	if( !is.list(value$degradation) )
		stop('llrtests: degradation element of value argument must be a list')
	if( !is.list(value$processing) )
		stop('llrtests: processing element of value argument must be a list')
	if( any(sapply(value$synthesis, length)!=2) )
		stop('llrtests: all elements of synthesis element of value must be of length 2')
	if( any(sapply(value$degradation, length)!=2) )
		stop('llrtests: all elements of degradation element of value must be of length 2')
	if( any(sapply(value$processing, length)!=2) )
		stop('llrtests: all elements of processing element of value must be of length 2')
	if( !all(unlist(sapply(value$synthesis, strsplit, split='')) %in% c('0','a','b','c')) )
		stop('llrtests: all elements of synthesis element of value contains not valid characters. Valid characters are: "0", "a", "b", "c"')
	if( !all(unlist(sapply(value$degradation, strsplit, split='')) %in% c('0','a','b','c')) )
		stop('llrtests: all elements of degradation element of value contains not valid characters. Valid characters are: "0", "a", "b", "c"')
	if( !all(unlist(sapply(value$processing, strsplit, split='')) %in% c('0','a','b','c')) )
		stop('llrtests: all elements of processing element of value contains not valid characters. Valid characters are: "0", "a", "b", "c"')
	# update the object
	if( !identical(value, object@model@params$llrtests) ) {
		object@model@params$llrtests <- value
		# if modeled rates already exist, update them
		if( length(object@model@ratesSpecs)>1 ) {
			message('Updating modeled rates...')
			object <- makeModelRates(object)
		}
	}
	return(object)
	})
#' @rdname testingParams
setMethod('llrtests', 'INSPEcT_model', function(object) {
	return(object@params$llrtests)
	})
#' @rdname testingParams
setReplaceMethod('llrtests', 'INSPEcT_model', function(object, value) {
	if( !is.list(value) )
		stop('llrtests: value argument must be a list')
	if( !identical(names(value), c('synthesis', 'degradation', 'processing')) )
		stop('llrtests: value argument must be a named list. Names must be "synthesis", "degradation" and "processing"')
	if( !is.list(value$synthesis) )
		stop('llrtests: synthesis element of value argument must be a list')
	if( !is.list(value$degradation) )
		stop('llrtests: degradation element of value argument must be a list')
	if( !is.list(value$processing) )
		stop('llrtests: processing element of value argument must be a list')
	if( any(sapply(value$synthesis, length)!=2) )
		stop('llrtests: all elements of synthesis element of value must be of length 2')
	if( any(sapply(value$degradation, length)!=2) )
		stop('llrtests: all elements of degradation element of value must be of length 2')
	if( any(sapply(value$processing, length)!=2) )
		stop('llrtests: all elements of processing element of value must be of length 2')
	if( !all(unlist(sapply(value$synthesis, strsplit, split='')) %in% c('0','a','b','c')) )
		stop('llrtests: all elements of synthesis element of value contains not valid characters. Valid characters are: "0", "a", "b", "c"')
	if( !all(unlist(sapply(value$degradation, strsplit, split='')) %in% c('0','a','b','c')) )
		stop('llrtests: all elements of degradation element of value contains not valid characters. Valid characters are: "0", "a", "b", "c"')
	if( !all(unlist(sapply(value$processing, strsplit, split='')) %in% c('0','a','b','c')) )
		stop('llrtests: all elements of processing element of value contains not valid characters. Valid characters are: "0", "a", "b", "c"')
	object@params$llrtests <- value
	return(object)
	})

