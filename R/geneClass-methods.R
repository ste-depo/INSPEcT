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
#' thresholds(nascentInspObj10)$chisquare <- .2
setMethod('geneClass', 'INSPEcT_model', 
	function(object, bTsh=NULL, cTsh=NULL) {
		## get ratesSpec field
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

.bestModel <- function(object, bTsh=NULL, cTsh=NULL) {

	preferPValue <- object@params$preferPValue
	
	## in case bTsh or bTsh are provided set them as
	# permanent for the object
	if( is.null(bTsh) )
		bTsh <- object@params$thresholds$brown
	if( is.null(cTsh) )
		cTsh <- object@params$thresholds$chisquare
	## calculate ratePvals
	ratePvals <- ratePvals(object, cTsh)
	ratePvals <- replace(ratePvals,is.na(ratePvals),1)

	if(preferPValue)
	{
		pValues <- chisqtest(object)
		rownames(pValues) <- rownames(ratePvals)

		geneClass <- sapply(rownames(ratePvals),function(i)
		{

			acceptableModelsTemp <- which(pValues[i,] <= cTsh)
			if(length(acceptableModelsTemp)==0){return("0")}
			if(length(acceptableModelsTemp)==1){return(colnames(pValues)[acceptableModelsTemp])}    

			nameTemp <- rep("K",3)
			if(ratePvals[i,"synthesis"] <= bTsh["synthesis"]){nameTemp[[1]] <- "V"}
			if(ratePvals[i,"processing"] <= bTsh["processing"]){nameTemp[[2]] <- "V"}
			if(ratePvals[i,"degradation"] <= bTsh["degradation"]){nameTemp[[3]] <- "V"}

			nameTemp <- paste0(nameTemp[[1]],nameTemp[[2]],nameTemp[[3]])
			nameTemp <- c("KKK" = "0"
						 ,"VKK" = "a"
						 ,"KVK" = "c"
						 ,"KKV" = "b"
						 ,"VVK" = "ac"
						 ,"VKV" = "ab"
						 ,"KVV" = "bc"
						 ,"VVV" = "abc")[nameTemp]

			if(pValues[i,nameTemp]<=cTsh|!is.finite(pValues[i,nameTemp])){return(nameTemp)}else{
				return(names(which.min(pValues[i,])))}
    	})

		ratesSpecs <- object@ratesSpecs
		## select the best model (according to geneClass) per gene
		nGenes <- length(ratesSpecs)
		object@ratesSpecs <- lapply(1:nGenes, 
			function(i) ratesSpecs[[i]][geneClass[i]])
		return(object)

	}else{

		## give a discrete classification per each rate per each gene
		# according to the brown's threshold for the pvalues
		acceptedVarModels <- sapply(1:3, function(i) ratePvals[,i]<bTsh[i])
		if( !is.matrix(acceptedVarModels) )
			acceptedVarModels <- t(as.matrix(acceptedVarModels))
		# nonResolvedGenes <- apply(acceptedVarModels, 1, 
		# 	function(x) all(is.na(x)))
		acceptedVarModels[is.na(acceptedVarModels)] <- FALSE
		rownames(acceptedVarModels) <- rownames(ratePvals)
		colnames(acceptedVarModels) <- colnames(ratePvals)
		geneClass <- apply(acceptedVarModels, 1, 
			function(accepted) paste(c('a','b','c')[accepted],collapse=''))
		geneClass[geneClass==''] <- '0'
		# geneClass[nonResolvedGenes] <- NA
		## retrive all the models
		ratesSpecs <- object@ratesSpecs
		## select the best model (according to geneClass) per gene
		nGenes <- length(ratesSpecs)
		object@ratesSpecs <- lapply(1:nGenes, 
			function(i) ratesSpecs[[i]][geneClass[i]])
		return(object)

	}
}

