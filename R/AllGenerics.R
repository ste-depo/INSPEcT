####################################
# generics for class INSPEcT_model ####
####################################

#' @rdname testingParams
setGeneric('modelSelection', function(object) 
	standardGeneric('modelSelection'))
#' @rdname testingParams
setGeneric('modelSelection<-', function(object, value) 
	standardGeneric('modelSelection<-'))
#' @rdname testingParams
setGeneric('thresholds', function(object) 
	standardGeneric('thresholds'))
#' @rdname testingParams
setGeneric('thresholds<-', function(object, value) 
	standardGeneric('thresholds<-'))
#' @rdname testingParams
setGeneric('llrtests', function(object) 
	standardGeneric('llrtests'))
#' @rdname testingParams
setGeneric('llrtests<-', function(object, value) 
	standardGeneric('llrtests<-'))
#' Retrieve all results of chi-squared test
setGeneric('chisqtest', function(object, ...) 
	standardGeneric('chisqtest'))
#' Retrieve results of chi-squared test for the selected models
setGeneric('chisqmodel', function(object, ...) 
	standardGeneric('chisqmodel'))
#' Retrieve results of log likelihood test
setGeneric('logLik', function(object, ...) 
	standardGeneric('logLik'))
#' Retrieve a single p-value for each rate
setGeneric('ratePvals', function(object, cTsh=NULL) 
	standardGeneric('ratePvals'))
#' Retrieve the regulatory class for each gene
setGeneric('geneClass', function(object, bTsh=NULL, cTsh=NULL)
	standardGeneric('geneClass'))
#' Calculate modeled rates and concentrations
setGeneric('makeModelRates', function(object, ...) 
	standardGeneric('makeModelRates'))
#' Generate synthetic rates and concentrations
setGeneric('makeSimDataset', function(object, tpts, nRep, No4sU = FALSE, seed=NULL)
	standardGeneric('makeSimDataset'))
#' Display rate classification performance
setGeneric('rocCurve', function(object, object2, cTsh=NULL, plot=TRUE) 
	standardGeneric('rocCurve'))
#' Display rate classification performance with thresholds visible at x-axis
setGeneric('rocThresholds', function(object, object2, cTsh=NULL, bTsh=NULL, xlim=c(1e-5,1)) 
	standardGeneric('rocThresholds'))

#############################
# generics for class INSPEcT ####
###############################

#' Accessor to the slot tpts of an INSPEcT object
setGeneric('tpts', function(object) 
	standardGeneric('tpts'))
#' Accessor to the slot totalSF of an INSPEcT object
setGeneric('totalSF', function(object) 
	standardGeneric('totalSF'))
#' Accessor to the slot labeledSF of an INSPEcT object
setGeneric('labeledSF', function(object) 
	standardGeneric('labeledSF'))
#' remove modelling information from INSPEcT object
setGeneric('removeModel', function(object) 
	standardGeneric('removeModel'))
#' Get the number of genes within the INSPEcT object
setGeneric('nGenes', function(object) 
	standardGeneric('nGenes'))
#' Get the number of time points within the INSPEcT object
setGeneric('nTpts', function(object) 
	standardGeneric('nTpts'))
#' Get and set number parameters for the modeling
setGeneric('modelingParams', function(object) 
	standardGeneric('modelingParams'))
#' @rdname modelingParams
setGeneric('modelingParams<-', function(object, value) 
	standardGeneric('modelingParams<-'))
#' Get or replace INSPEcT_model object within INSPEcT object
setGeneric('getModel', function(object) 
	standardGeneric('getModel'))
#' @rdname getModel
setGeneric('getModel<-', function(object, value) 
	standardGeneric('getModel<-'))
# #' Get and set number of cores to be used for the modeling
# setGeneric('nCores', function(object, ...) 
# 	standardGeneric('nCores'))
# #' @rdname nCores
# setGeneric('nCores<-', function(object, value) 
# 	standardGeneric('nCores<-'))
#' A nice plot to see scaling factors used for RNA-seq and 4sU-seq libraries
setGeneric('sfPlot', function(object) 
	standardGeneric('sfPlot'))
#' Retrieve pre-modeling rates and concentrations
setGeneric('ratesFirstGuess', function(object, feature) 
	standardGeneric('ratesFirstGuess'))
#' Retrieve pre-modeling rates and concentrations variance
setGeneric('ratesFirstGuessVar', function(object, feature) 
	standardGeneric('ratesFirstGuessVar'))
#' @title Launch the modeling process
#' @description Launch the modeling process with parameters set with \code{\link{modelingParams}}
setGeneric('modelRates', function(object, seed=NULL, BPPARAM=bpparam(), verbose=NULL) 
	standardGeneric('modelRates'))
#' Build the synthetic rates shaped on a dataset
setGeneric('makeSimModel', function(object, nGenes, newTpts=NULL
		, probs=c(constant=.5,sigmoid=.3,impulse=.2), na.rm=TRUE, seed=NULL) 
	standardGeneric('makeSimModel'))
#' Retrieve the modeled rates and concentrations
setGeneric('viewModelRates', function(object, feature) 
	standardGeneric('viewModelRates'))
#' Plot the pre-modeled and modeled profiles for one gene
setGeneric('plotGene', function(object, ix, fix.yaxis=FALSE) 
	standardGeneric('plotGene'))
#' Heatmap that represent the fold changes of all the five features
setGeneric('inHeatmap', function(object, type='pre-model'
	, breaks=seq(-1,1,length.out=51)
	, palette=colorRampPalette(c("green", "black", "firebrick3"))
	, plot_matureRNA=FALSE, absoluteExpression=TRUE
	, rowLabels=NULL, clustering=TRUE, clustIdx=3:5)
	standardGeneric('inHeatmap'))

#' Generate an object of class INSPEcT_diffsteady from two objects of class INSPEcT
setGeneric('compareSteady', function(inspectIds1, inspectIds2) 
	standardGeneric('compareSteady'))

##########################################
# generics for class INSPEcT_diffsteady ####
##########################################

#' @rdname INSPEcT_diffsteady-class
setGeneric('synthesis', function(object) standardGeneric('synthesis'))
#' @rdname INSPEcT_diffsteady-class
setGeneric('processing', function(object) standardGeneric('processing'))
#' @rdname INSPEcT_diffsteady-class
setGeneric('degradation', function(object) standardGeneric('degradation'))

# #############################
# # generics for class TxDb ####
# ###############################

# #' @rdname makeGtfFromDb
# setGeneric( 'makeExonsGtfFromDb' , function( object , type, filename ) standardGeneric( 'makeExonsGtfFromDb' ) )
# #' @rdname makeGtfFromDb
# setGeneric( 'makeIntronsGtfFromDb' , function( object , type, filename ) standardGeneric( 'makeIntronsGtfFromDb' ) )

