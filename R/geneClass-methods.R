#' @rdname geneClass
#'
#' @description
#' This method returns a factor that summarise the gene class (transcriptional regulatory mechanism) that
#' INSPEcT has assigned to each gene. The variability of each rate is indicated with a letter, 's' for
#' synthesis, 'p' for processing and 'd' for degradation. In case more than one rate is variable, the 
#' letters associated to each variable rate are merged, for example 'sd' stands for a gene where synthesis
#' and degradation cotributed to transcriptional changes. 'no-reg' is associated to genes with no
#' change in transcription. The classification depends on the thresholds of the goodness-of-fit and
#' rate variability tests that can be changed via the method \code{\link{calculateRatePvals}}.
#' @param object An object of class INSPEcT or INSPEcT_model
#' @return A character containing the regulatory class for each gene
#' @seealso \code{\link{ratePvals}}
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' geneClass(nascentInspObj10)
#' # see the classification with another threshold for rate variability
#' nascentInspObj10 <- calculateRatePvals(nascentInspObj10, p_variability=rep(1,3))
#' geneClass(nascentInspObj10)
setMethod('geneClass', 'INSPEcT', function(object, ...)
{
	if( !.hasSlot(object, 'version') ) {
		stop("This object is OBSOLETE and cannot work with the current version of INSPEcT.")
	}
	if( !is.numeric(tpts(object)) ) {
		stop("Run 'compareSteady' method on this object to evaluate differential rates.")
	}
	p_variability <- modelSelection(object)$p_variability
	ratePvalsTmp <- ratePvals(object)
	geneClass <- apply(ratePvalsTmp,1,function(r)
	{
		r <- r < unlist(p_variability)
		if(all(is.na(r))) return(NA)
		if(!r[1]&!r[2]&!r[3]) return("no-reg") # 0
		if(r[1]&!r[2]&!r[3]) return("s") # a
		if(!r[1]&r[2]&!r[3]) return("p") # c
		if(!r[1]&!r[2]&r[3]) return("d") # b
		if(r[1]&!r[2]&r[3]) return("sd") # ab
		if(r[1]&r[2]&!r[3]) return("sp") # ac
		if(!r[1]&r[2]&r[3]) return("pd") # bc
		if(r[1]&r[2]&r[3]) return("spd") # abc
	})
	return(geneClass)
	
})

#' @rdname geneClass
setMethod('geneClass', 'INSPEcT_model', function(object, ...)
{
	acceptedVarModels <- do.call('rbind', lapply(object@ratesSpecs, function(geneRates) 
		sapply(geneRates[[1]][c('alpha','beta','gamma')], '[[', 'df')>1))
	## transform the previous information into a string character per gene
	# where the presence of the letter means that the rate is variable
	allResponses <- apply(acceptedVarModels, 1, 
												function(accepted) paste(c('a','b','c')[accepted],collapse=''))
	allResponses[allResponses==''] <- '0'
	return(allResponses)
})

geneClassInternal <- function(object) {
	p_variability <- modelSelection(object)$p_variability
	ratePvalsTmp <- ratePvals(object)
	geneClass <- apply(ratePvalsTmp,1,function(r)
	{
		r <- r < unlist(p_variability)
		if(all(is.na(r))) return(NA)
		if(!r[1]&!r[2]&!r[3]) return("0") # no-reg 
		if(r[1]&!r[2]&!r[3]) return("a") # s 
		if(!r[1]&r[2]&!r[3]) return("c") # p 
		if(!r[1]&!r[2]&r[3]) return("b") # d 
		if(r[1]&!r[2]&r[3]) return("ab") # sd 
		if(r[1]&r[2]&!r[3]) return("ac") # sp 
		if(!r[1]&r[2]&r[3]) return("bc") # pd 
		if(r[1]&r[2]&r[3]) return("abc") # spd 
	})
	return(geneClass)
}
