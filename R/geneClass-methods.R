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
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' geneClass(nascentInspObj10)
#' # see the classification with another threshold for chi-squared test 
#' geneClass(nascentInspObj10, cTsh=.2)
#' # set the new threshold permanently within the object
#' modelSelection(nascentInspObj10)$thresholds$chisquare <- .2
setMethod('geneClass', 'INSPEcT_model', 
	function(object, bTsh=NULL, cTsh=NULL) {
		## get ratesSpec field
		ratesSpecs <- object@ratesSpecs
		## in case some elements of ratesSpecs are longer than one,
		# meaning that a unique choiche for a model has not been done yet,
		# choose one using "bestModel" method
		if( any(sapply(ratesSpecs, length)!=1) )
			ratesSpecs <- .bestModel(object, bTsh, cTsh)@ratesSpecs
		geneClass <- sapply(ratesSpecs,names)
		if(all(is.null(unlist(sapply(ratesSpecs,names)))))
		{
			# Standard names: alpha, beta and gamma for all
			#
			
			## get a logical matrix with 3 colums per gene stating wheter
			# alpha, beta or gamma are varible or not
	
			acceptedVarModels <- do.call('rbind', lapply(ratesSpecs, function(geneRates) 
				sapply(geneRates[[1]][c('alpha','beta','gamma')], '[[', 'df')>1))
			## transform the previous information into a string character per gene
			# where the presence of the letter means that the rate is variable
			geneClass <- apply(acceptedVarModels, 1, 
				function(accepted) paste(c('a','b','c')[accepted],collapse=''))
			geneClass[geneClass==''] <- '0'
		}
		names(geneClass) <- if( !is.null(names(object@ratesSpecs)) ) names(object@ratesSpecs) else 
			seq_along(object@ratesSpecs)
		## return
		# names(geneClass) <- seq_along(geneClass)
		return(geneClass)
	})

#' @rdname geneClass
setMethod('geneClass', 'INSPEcT', function(object, bTsh=NULL, cTsh=NULL)
{
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	if(object@NoNascent & !object@NF){return(geneClass(object@model, bTsh=bTsh, cTsh=cTsh))}
	else{.bestModel_confidenceIntervals(object, bTsh = modelSelection(object)$thresholds$brown)}

})

.bestModel_confidenceIntervals <- function(object,bTsh=NULL)
{
	if(is.null(bTsh)){bTsh <- modelSelection(object)$thresholds$brown}

	ratePvalsTmp <- ratePvals(object)
	geneClass <- apply(ratePvalsTmp,1,function(r)
	{
		r <- r < unlist(bTsh)
		if(all(is.na(r))) return(NA)
		if(!r[1]&!r[2]&!r[3]) return("0")
		if(r[1]&!r[2]&!r[3]) return("a")
		if(!r[1]&r[2]&!r[3]) return("c")
		if(!r[1]&!r[2]&r[3]) return("b")
		if(r[1]&!r[2]&r[3]) return("ab")
		if(r[1]&r[2]&!r[3]) return("ac")
		if(!r[1]&r[2]&r[3]) return("bc")
		if(r[1]&r[2]&r[3]) return("abc")
	})
	return(geneClass)
}

.bestModel <- function(object, bTsh=NULL, cTsh=NULL, Nascent = FALSE) {
	
	preferPValue <- object@params$preferPValue
	
	## in case bTsh or bTsh are provided set them as
	# permanent for the object
	if( is.null(bTsh) )
		bTsh <- object@params$thresholds$brown
	if( is.null(cTsh) )
		cTsh <- object@params$thresholds$chisquare

	if(Nascent) # It must be always VVV, just for makeModelRates
	{
		ratePvals <- matrix(rep(0,length(object@ratesSpecs)*3),ncol=3)
		rownames(ratePvals) <- names(object@ratesSpecs)
		colnames(ratePvals) <- c("synthesis","processing","degradation")
	}else{
		## calculate ratePvals
		ratePvals <- ratePvals(object, bTsh, cTsh)
		ratePvals <- replace(ratePvals,is.na(ratePvals),1)
	}

	geneClass <- apply(ratePvals,1,function(r)
	{
		r <- r < unlist(bTsh)
		if(all(is.na(r))) return(NA)
		if(!r[1]&!r[2]&!r[3]) return("0")
		if(r[1]&!r[2]&!r[3]) return("a")
		if(!r[1]&r[2]&!r[3]) return("c")
		if(!r[1]&!r[2]&r[3]) return("b")
		if(r[1]&!r[2]&r[3]) return("ab")
		if(r[1]&r[2]&!r[3]) return("ac")
		if(!r[1]&r[2]&r[3]) return("bc")
		if(r[1]&r[2]&r[3]) return("abc")
	})

	ratesSpecs <- object@ratesSpecs
	nGenes <- length(ratesSpecs)

	object@ratesSpecs <- lapply(1:nGenes, function(i) ratesSpecs[[i]][geneClass[i]])
	return(object)
	
}

