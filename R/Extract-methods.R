#' Extract Parts of an INSPEcT or an INSPEcT_model Object
#'
#' @description
#' Operators acting on INSPEcT, INSPEcT_model or INSPEcT_diffsteady objects 
#' to extract parts. INSPEcT_model objects can be subsetted only by gene.
#' INSPEcT objects can be subsetted either by gene id or time point. In case
#' of subsetting an INSPEcT object by time point, the model should be empty.
#' @seealso removeModel
#'
#' @docType methods
#' @name Extract
NULL


#' @rdname Extract
#' @param x An object of class INSPEcT or INSPEcT_model
#' @param i A numeric, a vector of logicals or a vector of names indicating the 
#' features to be extracted
#' @return An Object of class INSPEcT
setMethod('[', 'INSPEcT_model', function(x, i) {
		if( length(x@ratesSpecs)>0 ) {
			x@ratesSpecs <- x@ratesSpecs[i]
			x@params$sim$noiseVar <- lapply(x@params$sim$noiseVar
				, function(x) x[i])
			x@params$sim$foldchange <- lapply(x@params$sim$foldchange
				, function(x) x[i])
		}
		return( x )
	})

#' @rdname Extract
#' @param j A numeric, a vector of logicals indicating the 
#' time points to be extracted
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' nascentInspObj10_5genes <- nascentInspObj10[1:5]
#' \dontrun{
#' ## This will turn out into an error:
#' nascentInspObj10_5genes_5tpts <- nascentInspObj10[1:5, 1:5]
#' }
#' ## Before subsetting time points, the model should be removed:
#' nascentInspObj10_5genes_5tpts <- removeModel(nascentInspObj10)[1:5, 1:5]
setMethod('[', 'INSPEcT', function(x, i, j) {
	# subset the expressionSet slots (if populated)
	if( !missing(i) ) {
		if( is.logical(i) ) i <- which(i)
		if( nrow(x@ratesFirstGuess)>0 ) x@ratesFirstGuess <- x@ratesFirstGuess[i]
		if( nrow(x@ratesFirstGuessVar)>0 ) x@ratesFirstGuessVar <- x@ratesFirstGuessVar[i]
		if( nrow(x@modelRates)>0 ) x@modelRates <- x@modelRates[i]
		# subset the INSPEcT_model slot
		x@model <- x@model[i]		
	}
	if( !missing(j) ) {
		if( length(x@model@ratesSpecs) > 0 ) {
			stop('Remove the model before subsetting time points. (See "?removeModel")')
		} else {
			if( is.logical(j) ) j <- which(j)
			if( ncol(x@ratesFirstGuess)>0 ) {
				x@tpts <- x@tpts[j]
				x@labeledSF <- x@labeledSF[j]
				ix <- pData(x@ratesFirstGuess)$time %in% x@tpts
				x@ratesFirstGuess <- x@ratesFirstGuess[,ix]
				ix <- pData(x@ratesFirstGuessVar)$time %in% x@tpts
				x@ratesFirstGuessVar <- x@ratesFirstGuessVar[,ix]
			}
		}
	}
	# return
	return( x )
	})

#' @rdname Extract
setMethod('[', 'INSPEcT_diffsteady', function(x, i, j) {
	if( !missing(i) ) {
		x@synthesis <- x@synthesis[i,]
		x@degradation <- x@degradation[i,]
		x@processing <- x@processing[i,]
	}
	if( !missing(j) ) {
		x@synthesis <- x@synthesis[,j]
		x@degradation <- x@degradation[,j]
		x@processing <- x@processing[,j]
	}
	return( x )
	})


#' @rdname removeModel
#' 
#' @description
#' Remove the model from an INSPEcT object. It is required when subsetting
#' an INSPEcT object per time points because when removing time points
#' the modeling is not valid anymore.
#' @param object An Object of class INSPEcT
#' @return An Object of class INSPEcT
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' nascentInspObj10_5genes <- nascentInspObj10[1:5]
#'
#' ## This will turn out into an error:
#' \dontrun{nascentInspObj10_5genes_5tpts <- nascentInspObj10[1:5, 1:5]
#' }
#' 
#' ## Before subsetting time points, the model should be removed:
#' nascentInspObj10_5genes_5tpts <- removeModel(nascentInspObj10)[1:5, 1:5]
#'
#' ## Also this will turn out into an error:
#' \dontrun{nascentInspObj10 <- modelRates(nascentInspObj10)}
#' 
#' ## Before running the model again, or changing modeling parameters,
#' ## the previous model should be removed:
#' nascentInspObj10_old <- nascentInspObj10
#' nascentInspObj10_new <- removeModel(nascentInspObj10)
#' modelingParams(nascentInspObj10_new)$useSigmoidFun <- FALSE
#' \dontrun{nascentInspObj10_new <- modelRates(nascentInspObj10_new)}
setMethod('removeModel', 'INSPEcT', function(object) {
	object@model@ratesSpecs <- list()
	return(object)
	})

#' Combine different Objects of Class INSPEcT
#'
#' @description
#' This method combines the information coming from different Objects of INSPEcT class.
#' Requirements for two or more object to be combined together are:
#' \itemize{
#' \item they must be either modeled or either not modeled
#' \item they must have the same time points
#' \item they must have the same modeling parameters
#' }
#'
#' @param x An object of class INSPEcT
#' @param y An object of class INSPEcT
#' @param ... Additional objects of class INSPEcT
#' @return An Object of class INSPEcT
#'
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' nascentInspObj10_2genes <- nascentInspObj10[1:2]
#' nascentInspObj10_5genes <- nascentInspObj10[6:10]
#' nascentInspObj10_7genes <- combine(nascentInspObj10_2genes, nascentInspObj10_5genes)
#'
#' @details
#' In case the same gene is contained in more than one object that the user
#' tries to combine, the information from one object will be used and a 
#' warning will be reported
#'
#' @docType methods
#' @name combine
NULL

#' @rdname combine
setMethod('combine', signature(x='INSPEcT', y='INSPEcT'), function(x, y, ...) {
	dots <- c(x, y, list(...))
	dotsClasses <- sapply(dots, class)
	if( !all(dotsClasses == 'INSPEcT') )
		stop('combine: some of the elemets of provided as argument in not of class INSPEcT')
	modeledObjects <- sapply(dots, function(x) length(x@model@ratesSpecs)>0)
	if( any(modeledObjects) && !all(modeledObjects) )
		stop('combine: either all the object provided should be modeled or not. Model all the objects or use method "removeModel" to remove the models.')
	if( !all(sapply(dots[-1], function(x) identical(x@model@params, dots[[1]]@model@params))) )
		stop('combine: testing parameters are different. Modify them via "modelSelection", "thresholds" and "llrtests"')
	if( !all(sapply(dots[-1], function(x) identical(x@params[!names(x@params) %in% c('verbose')], dots[[1]]@params[!names(x@params) %in% c('verbose')]))) )
		stop('combine: testing parameters are different. Modify them via "modelSelection", "thresholds" and "llrtests"')
	if( !all(sapply(dots[-1], function(x) identical(x@tpts, dots[[1]]@tpts))) )
		stop('combine: trying to merging objects which contains different time points')
	if( any(duplicated(do.call('c', lapply(dots, featureNames)))) )
		warning('combine: there are genes that are contained in more than one object: only one is kept')
	# re-biuld the object
	newObject <- new('INSPEcT')
	newObject@model@params <- dots[[1]]@model@params
	newObject@params[!names(x@params) %in% c('verbose')] <- dots[[1]]@params[!names(x@params) %in% c('verbose')]
	newObject@tpts <- dots[[1]]@tpts
	if( all(sapply(dots[-1], function(x) identical(x@labeledSF, dots[[1]]@labeledSF))) )
		newObject@labeledSF <- dots[[1]]@labeledSF
	else
		newObject@labeledSF <- rep(NA, length(dots[[1]]@tpts))
	if( all(sapply(dots[-1], function(x) identical(x@tL, dots[[1]]@tL))) )
		newObject@tL <- dots[[1]]@tL
	else
		newObject@tL <- NA
	newObject@ratesFirstGuess <- do.call('combine', lapply(dots, function(x) x@ratesFirstGuess))
	if( all(modeledObjects) ) {
		newObject@modelRates <- do.call('combine', lapply(dots, function(x) x@modelRates))
		ratesSpecs <- do.call('c', lapply(dots, function(x) x@model@ratesSpecs))
		newObject@model@ratesSpecs <- ratesSpecs[featureNames(newObject@modelRates)]
	}
	return(newObject)

	})

#' Divide an INSPEcT Object into groups
#'
#' @description
#' Divides the INSPEcT object into the groups defined by 'f',
#' @param x An object of class INSPEcT
#' @param f A vector of length equal to the number of genes in x which defines the groups
#' @param drop A logical belonging to the generic funciton, useless in this context.
#' @param ... Additional arguments to match the generic function
#' @return A list containing objects of class INSPEcT
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' splitIdx <- c(1,1,1,2,2,2,3,3,3,4)
#' nascentInspObj10Split <- split(nascentInspObj10, splitIdx)
#' @docType methods
#' @name split
NULL


#' @rdname split
setMethod('split', 'INSPEcT', function(x, f, drop = FALSE, ...) {
	if( nGenes(x) != length(f) )
		stop('split: length of f must match the number of genes of x')
	f <- as.factor(f)
	return(lapply(levels(f), function(l) x[f==l]))
	})
