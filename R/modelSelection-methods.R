#' Get or set parameters for model test and selection
#'
#' @description
#' With this methods the user can personalize the criteria by which INSPEcT
#' selects a rate to be variable or constant. In particular, the model selection
#' criteria can be selected between log-likelihood ratio test and Akaike's information
#' criterion (AIC). In case log-likelihood ratio test is selected, the thresholds of 
#' chi-squared and Brown's method can be set (see Details section).
#'
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param value A list or a character that will substitute the set of parameters
#' \itemize{
#'   \item modelSelection: A character, either "llr" to test whether a rate is 
#'     varying using log-likelihood testing framework or "aic" to choose
#'     the best model via Akaike Information Criterion (Default: "llr").
#'   \item thresholds: A named list containing the thresholds for the goodness of fit 
#'     (chisquare) and variability (brown) tests. Both must be comprised from 0 to 1.
#'     The goodness of fit threshold defines which models are considered valid or not
#'     (0 most stringent, 1 most permissive). The variability threshold (one per each rate)
#'     defines the stringency of the call for the variability of each rate (0 the most
#'     stringent, 1 the most permissive). When set to 0, the specific rate is 
#'     excluded from the hypothesis of variability.
#'   \item preferPValue: when model selection is "llr", preferPValue means that 
#'     if the selected model has a goodness of fit below threshold the model with
#'     the best goodness of fit is returned in place of it. When the model selection
#'     is "aic", with preferPValue the best model is tested against the closest 
#'     nested models to test the hypothesis of variability and only when this 
#'     pvalue is below threshold the rate is considered as varible. Otherwise the 
#'     model selection is just based on the lowest AIC (Default: TRUE).
#'   \item padj: whether to correct pvalues with the Benjamini-Hochberg procedure 
#'     or not. (Default: TRUE).
#'   \item limitModelComplexity: limit the complexity of the models associated to
#'     variable rates to the length of the time courses. Particularly helpful
#'     for short time series. (Default:FALSE)
#'     
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
	if( value$limitModelComplexity != pre_val$limitModelComplexity ) {
		message('Updating models accoring to the new complexity...')
		tclength <- length(tpts(object))
		nGenes <- length(object@model@ratesSpecs)
		gene_classes <- names(object@model@ratesSpecs[[1]])

		if( object@model@params$limitModelComplexity ) {
			k <- t(sapply(1:nGenes, function(i) 
				sapply(object@model@ratesSpecs[[1]], function(x) {
					sum(sapply(x[1:3], function(y) min(y$df, tclength)))
				})))
		} else {
			k <- t(sapply(1:nGenes, function(i) 
				sapply(object@model@ratesSpecs[[1]], function(x) {
					sum(sapply(x[1:3], '[[', 'df'))
				})))
		}

		tot_exp = ratesFirstGuess(object, 'total')
		pre_exp = ratesFirstGuess(object, 'preMRNA')
		syn_exp = ratesFirstGuess(object, 'synthesis')

		tot_var = ratesFirstGuessVar(object, 'total')
		pre_var = ratesFirstGuessVar(object, 'preMRNA')
		syn_var = ratesFirstGuessVar(object, 'synthesis')

		D <- sapply(gene_classes, function(gc) {
			objectTmp <- object
			objectTmp@model@ratesSpecs <- lapply(1:nGenes, 
				function(i) objectTmp@model@ratesSpecs[[i]][gc])
			objectTmp <- makeModelRates(objectTmp)
			tot_mod = viewModelRates(objectTmp, 'total')
			pre_mod = viewModelRates(objectTmp, 'preMRNA')
			syn_mod = viewModelRates(objectTmp, 'synthesis')
			apply((tot_mod - tot_exp)^2/tot_var +
				(pre_mod - pre_exp)^2/pre_var +
				if( objectTmp@NoNascent ) 0 else (syn_mod - syn_exp)^2/syn_var
				, 1, sum)
		})
		n <- if( object@NoNascent ) 2*tclength else 3*tclength
		df <- n - k
		chisqtest = sapply(seq_along(gene_classes), function(i) {
			if(df[i]>0) pchisq(D[,i], df[i], lower.tail=TRUE) else rep(1, nGenes)
		})
		colnames(chisqtest) = gene_classes
		AIC <- - 2*logLik(object) + 2*k
		AICc <- AIC + 2*k*(k+1)/(n-k-1)
		## replace in the object
		for( i in 1:nGenes ) {
			for( gc in gene_classes ) {
				object@model@ratesSpecs[[i]][[gc]]$test <- log(chisqtest[i,gc])
				object@model@ratesSpecs[[i]][[gc]]$AIC <- AIC[i,gc]
				object@model@ratesSpecs[[i]][[gc]]$AICc <- AICc[i,gc]
			}
		}
	}
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
	if( !value$modelSelection %in% c('aic', 'llr', 'hib') )
		stop('modelSelection: value argument must either "llr", "aic" or "hib"')
	if( !is.logical(value$preferPValue) )
		stop('modelSelection: preferPValue element of value argument must be a logical')
	if( !is.logical(value$padj) )
		stop('modelSelection: padj element of value argument must be a logical')
	if( !is.logical(value$limitModelComplexity) )
		stop('modelSelection: limitModelComplexity element of value argument must be a logical')
	if( !is.list(value$thresholds) )
		stop('modelSelection: thresholds element of value argument must be a list')
	if( !identical(names(value$thresholds),c('chisquare', 'brown', 'CI')) )
		stop('modelSelection: thresholds element of value argument must be named. Names must be: "chisquare", "brown" or "CI')
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
