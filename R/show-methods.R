#' @rdname INSPEcT_model-class
#' @param object An object of class INSPEcT_model
#' @return Method show for objects of class INSPEcT_model returns the number of
#' the genes that have been modeled
setMethod('show', 'INSPEcT_model', function(object) {
	message(paste('object of class INSPEcT_model of length', 
		length(object@ratesSpecs)))
	})

#' @rdname INSPEcT-class
#' @param object An object of class INSPEcT
#' @return Method show for objects of class INSPEcT displays the main features
#' of the slots ratesFirstGuess, model and modelRates
setMethod('show', 'INSPEcT', function(object) {

	message('\nslot @ratesFirstGuess:')
	print(object@ratesFirstGuess)
	message('\nslot @model:')
	print(object@model)
	message('\nslot @modelRates:')
	print(object@modelRates)

	})

#' @rdname INSPEcT_diffsteady-class
#' @param object An object of class INSPEcT_model
#' @return Method show for objects of class INSPEcT_model returns the number of
#' the genes that have been modeled
setMethod('show', 'INSPEcT_diffsteady', function(object) {
	message('Object of class INSPEcT_diffsteady')

	for(sn in slotNames(object) ) {
		data <- slot(object, sn)
		if(nrow(data)>6) {
			message('Head of slot ',sn,':')
			print(head(data))
			message('... and other ', nrow(data)-6, ' hidden genes.')
		} else {
			message('Slot "',sn,'":')
			print(data)
		}
		message('')
	}

	})

#' @title Gene Names Associated with an Object of Class INSPEcT
#' @description
#' A method to visualize gene names associated with the object of class INSPEcT
#' @param object An object of class INSPEcT
#' @return A character that contains gene names associated with the object of class INSPEcT
setMethod('featureNames', 'INSPEcT', function(object) {
	featureNames(object@ratesFirstGuess)
	})

#' @rdname featureNames-INSPEcT-method
#' @param value A character that will replace the current feature names
setReplaceMethod('featureNames', signature(object='INSPEcT', value='ANY')
	, function(object, value) {
	featureNames(object@ratesFirstGuess) <- value
	if( nrow(object@modelRates) > 0 ) {
		featureNames(object@modelRates) <- value
		names(object@model@ratesSpecs) <- value
	}
	object
	})

#' @rdname nGenes
#' @description
#' A method to obtain the number of the genes associated with the object of class INSPEcT
#' @param object An object of class INSPEcT
#' @return A numeric that indicates the number of genes within the object
#' @examples
#' data('mycerIds10',package="INSPEcT")
#' nGenes(mycerIds10)
setMethod('nGenes', 'INSPEcT', function(object) {
	length(featureNames(object@ratesFirstGuess))
	})

#' @rdname nTpts
#' @description
#' A method to obtain the number of the tpts associated with the object of class INSPEcT
#' @param object An object of class INSPEcT
#' @return A numeric that indicates the number of time points contained the object
#' @examples
#' data('mycerIds10',package="INSPEcT")
#' nTpts(mycerIds10)
setMethod('nTpts', 'INSPEcT', function(object) {
	length(object@tpts)
	})

#' @title Dimensions of an Object of Class INSPEcT
#' @description
#' A method to obtain the dimension of the object of class INSPEcT reported as a vector
#' containing of the genes and the number of time points
#' @param x An object of class INSPEcT
#' @seealso \code{\link{nGenes}}, \code{\link{nTpts}}
#' @return A numeric that indicates the number of genes within the object and 
#' the number of time points contained the object
setMethod('dim', 'INSPEcT', function(x) {
	c(
	length(featureNames(x@ratesFirstGuess))
	, length(x@tpts)
	)
	})

#' @rdname tpts
#' @description
#' Accessor to obtain the tpts associated with the object of class INSPEcT
#' @param object An object of class INSPEcT
#' @return A numeric that indicates time points contained the object
#' @examples
#' data('mycerIds10',package="INSPEcT")
#' tpts(mycerIds10)
setMethod('tpts', 'INSPEcT', function(object) {
	object@tpts
	})

#' @rdname totalSF
#' @description
#' Accessor to obtain the totalSF slot associated with the object of class INSPEcT
#' @param object An object of class INSPEcT
#' @return A numeric that indicates the scaling factors applied between time points
#' of the data coming from RNA-seq library (applies to total-mRNAs
#' and pre-mRNAs concentrations)
#' @examples
#' data('mycerIds10',package="INSPEcT")
#' totalSF(mycerIds10)
setMethod('totalSF', 'INSPEcT', function(object) {
	object@totalSF
	})

#' @rdname labeledSF
#' @description
#' Accessor to obtain the labeledSF slot associated with the object of class INSPEcT
#' @param object An object of class INSPEcT
#' @return A numeric that indicates the scaling factors applied between time points
#' of the data coming from 4sU-seq library (applies directly to synthesis
#' rates and indirectly to degradation rates)
#' @examples
#' data('mycerIds10',package="INSPEcT")
#' labeledSF(mycerIds10)
setMethod('labeledSF', 'INSPEcT', function(object) {
	object@labeledSF
	})
