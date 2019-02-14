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
#' synthesis are changing". Different comparisons will be combined using Brown's method 
#' for combinig p-values.
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
#' @name modelSelection
NULL

#' @rdname modelSelection
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' modelSelection(nascentInspObj10)
#' modelSelection(nascentInspObj10)$modelSelection <- 'aic'
setMethod('modelSelection', 'INSPEcT', function(object) {
	return(modelSelection(object@model))
	})

#' @rdname modelSelection
setReplaceMethod('modelSelection', 'INSPEcT', function(object, value) {
	pre_val <- modelSelection(object@model)
	modelSelection(object@model) <- value
	if( (!identical(value, pre_val)) & length(object@model@ratesSpecs)>1 ) {
		message('Updating modeled rates...')
		bTsh <- value$thresholds$brown
		if( any(bTsh == 0) ) {
			rates_to_aviod <- names(bTsh)[bTsh == 0]
			message(paste('--',paste(rates_to_aviod, collapse=', and '),'rate(s) excluded from hypothesis of variability --'))
		}
		object <- makeModelRates(object)
	}
	return(object)
	})

#' @rdname modelSelection
setMethod('modelSelection', 'INSPEcT_model', function(object) {
	return(object@params)
	})

#' @rdname modelSelection
setReplaceMethod('modelSelection', 'INSPEcT_model', function(object, value) {
	if(!is.list(value))
		stop('modelSelection: value argument must be a list')
	for( arg in c('modelSelection','preferPValue','padj','thresholds','limitModelComplexity') ) {
		if( !any(grepl(paste0('^',arg), names(value))) )
			stop(paste('modelSelection: value must contain an element named "', arg, '"'))
	}
	# if( !identical(names(value),c('modelSelection','preferPValue','padj','thresholds')) )
	# 	stop('modelSelection: value argument list must be named. Names must be "modelSelection", "preferPValue", "padj" and"thresholds"')
	if( !is.character(value$modelSelection) )
		stop('modelSelection: value argument must be a character')
	if( !value$modelSelection %in% c('aic', 'llr') )
		stop('modelSelection: value argument must either "llr" or "aic"')
	if( !is.logical(value$preferPValue) )
		stop('modelSelection: preferPValue element of value argument must be a logical')
	if( !is.logical(value$padj) )
		stop('modelSelection: padj element of value argument must be a logical')
	if( !is.logical(value$limitModelComplexity) )
		stop('modelSelection: limitModelComplexity element of value argument must be a logical')
	if( !is.list(value$thresholds) )
		stop('modelSelection: thresholds element of value argument must be a list')
	if( !identical(names(value$thresholds),c('chisquare', 'brown')) )
		stop('modelSelection: thresholds element of value argument must be named. Names must be "chisquare" and "brown"')
	if( !is.numeric(value$thresholds$chisquare) )
		stop('modelSelection: chisquare element of value argument must be a numeric')
	if( !is.numeric(value$thresholds$brown) )
		stop('modelSelection: brown element of value argument must be a numeric')
	if( length(value$thresholds$chisquare) != 1 )
		stop('modelSelection: chisquare element of value argument must be of length 1')
	if( length(value$thresholds$brown) != 3 )
		stop('modelSelection: brown element of value argument must be of length 3')
	object@params <- value
	return(object)
	})
