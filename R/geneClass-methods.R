#' @rdname geneClass
#'
#' @description
#' This method returns a factor that summarise the gene class (transcriptional regulatory mechanism) that
#' INSPEcT has assigned to each gene. The classification depends on the chi-squared and Brown's method
#' thresholds, that can be both provided as arguments. If the user decides a different thresholding respect to
#' the default, these new values can be permanently set within the object.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @param bTsh A numeric representing the p-value threshold for considering a rate as variable. P-values are calculated through \code{\link{ratePvals}}
#' @param cTsh A numeric representing the threshold for the chi-squared test to consider a model as valid
#' @return A character containing the regulatory class for each gene
#' @seealso \code{\link{ratePvals}}
#' @examples
#' data('nascentInspObj10', package='INSPEcT')
#' geneClass(nascentInspObj10)
#' # see the classification with another threshold for chi-squared test 
#' geneClass(nascentInspObj10, cTsh=.2)
#' # set the new threshold permanently within the object
#' thresholds(nascentInspObj10)$chisquare <- .2
setMethod('geneClass', 'INSPEcT_model', 
	function(object, bTsh=NULL, cTsh=NULL) {
		## get ratesSpec field
		ratesSpecs <- object@ratesSpecs
		## in case some elements of ratesSpecs are longer than one,
		# meaning that a unique choiche for a model has not been done yet,
		# choose one using "bestModel" method
		if( any(sapply(ratesSpecs, length)!=1) )
			ratesSpecs <- .bestModel(object, bTsh, cTsh)@ratesSpecs
		## get a logical matrix with 3 colums per gene stating wheter
		# alpha, beta or gamma are varible or not
		acceptedVarModels <- do.call('rbind', lapply(ratesSpecs, function(geneRates) 
			sapply(geneRates[[1]][c('alpha','beta','gamma')], '[[', 'df')>1))
		## transform the previous information into a string character per gene
		# where the presence of the letter means that the rate is variable
		geneClass <- apply(acceptedVarModels, 1, 
			function(accepted) paste(c('a','b','c')[accepted],collapse=''))
		geneClass[geneClass==''] <- '0'
		names(geneClass) <- names(object@ratesSpecs)
		## return
		return(geneClass)
	})

#' @rdname geneClass
setMethod('geneClass', 'INSPEcT', 
	function(object, bTsh=NULL, cTsh=NULL) {
		return(geneClass(object@model, bTsh=bTsh, cTsh=cTsh))
	})
