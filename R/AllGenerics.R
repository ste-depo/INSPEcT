####################################
# generics for class INSPEcT_model ####
####################################

#' @rdname modelSelection
setGeneric('modelSelection', function(object) 
	standardGeneric('modelSelection'))
#' @rdname modelSelection
setGeneric('modelSelection<-', function(object, value) 
	standardGeneric('modelSelection<-'))
#' Retrieve all results of chi-squared test
setGeneric('chisqtest', function(object, ...) 
	standardGeneric('chisqtest'))
#' Retrieve results of chi-squared test for the selected models
setGeneric('chisqmodel', function(object, gc=NULL, tpts=NULL, ...) 
	standardGeneric('chisqmodel'))
#' Retrieve results of log likelihood test
setGeneric('logLik', function(object, ...) 
	standardGeneric('logLik'))
#' Retrieve a single p-value for each rate
setGeneric('ratePvals', function(object, bTsh=NULL, cTsh=NULL) 
	standardGeneric('ratePvals'))
#' Retrieve the regulatory class for each gene
setGeneric('geneClass', function(object, bTsh=NULL, cTsh=NULL)
	standardGeneric('geneClass'))
#' Calculate modeled rates and concentrations
setGeneric('makeModelRates', function(object, ...) 
	standardGeneric('makeModelRates'))
#' Compute confidence intervals
setGeneric('computeConfidenceIntervals', function(object, BPPARAM=bpparam()) 
	standardGeneric('computeConfidenceIntervals'))
#' Set confidence intervals
setGeneric('setConfidenceIntervals', function(object, confidenceIntervals) 
	standardGeneric('setConfidenceIntervals'))
#' Generate synthetic rates and concentrations
setGeneric('makeSimDataset', function(object, tpts, nRep, NoNascent = FALSE, seed=NULL, b = 0.3, tL = 1/6, noise_sd = 4.0)
	standardGeneric('makeSimDataset'))
#' Display rate classification performance
setGeneric('correlationPlot', function(object, object2, plot=TRUE) 
	standardGeneric('correlationPlot'))
#' Display rate classification performance
setGeneric('rocCurve', function(object, object2, cTsh=NULL, plot=TRUE, comparative=FALSE) 
	standardGeneric('rocCurve'))
#' Display rate classification performance with thresholds visible at x-axis
setGeneric('rocThresholds', function(object, object2, cTsh=NULL, bTsh=NULL, xlim=c(1e-5,1), plot=TRUE) 
	standardGeneric('rocThresholds'))

#############################
# generics for class INSPEcT ####
###############################

#' Accessor to the slot tpts of an INSPEcT object
setGeneric('tpts', function(object) 
	standardGeneric('tpts'))
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
#' A nice plot to see scaling factors used for RNA-seq and Nascent-seq libraries
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
#' @rdname modelRates
setGeneric('modelRatesNF', function(object, BPPARAM=SerialParam()) 
	standardGeneric('modelRatesNF'))
#' Build the synthetic rates shaped on a dataset
setGeneric('makeSimModel', function(object, nGenes, newTpts=NULL
		, probs=c(constant=.5,sigmoid=.3,impulse=.2), na.rm=TRUE, seed=NULL) 
	standardGeneric('makeSimModel'))
#' Build the synthetic rates with oscillatory pattern
setGeneric('makeOscillatorySimModel', function(object, nGenes
		, oscillatoryk3=FALSE, k3delay=NULL, na.rm=TRUE, seed=NULL) 
	standardGeneric('makeOscillatorySimModel'))
#' Retrieve the modeled rates and concentrations
setGeneric('viewModelRates', function(object, feature) 
	standardGeneric('viewModelRates'))
#' Retrieve the modeled Confidence Intervals
setGeneric('viewConfidenceIntervals', function(object, feature) 
	standardGeneric('viewConfidenceIntervals'))
#' Plot the pre-modeled and modeled profiles for one gene
setGeneric('plotGene', function(object, ix, fix.yaxis=FALSE, priors=TRUE) 
	standardGeneric('plotGene'))
#' Heatmap that represent the fold changes of all the five features
setGeneric('inHeatmap', function(object, type='pre-model'
	, breaks=seq(-1,1,length.out=51)
	, palette=colorRampPalette(c("green", "black", "firebrick3"))
	, plot_matureRNA=FALSE, absoluteExpression=TRUE
	, rowLabels=NULL, clustering=TRUE, clustIdx=3:5)
	standardGeneric('inHeatmap'))
#' Generate an object of class INSPEcT_diffsteady from an object of class INSPEcT
setGeneric('compareSteady', function(inspectIds, BPPARAM=bpparam()) 
	standardGeneric('compareSteady'))
#' Compare mature RNA steady state data from an object of class INSPEcT
setGeneric('compareSteadyNoNascent', function(inspectIds,
	expressionThreshold=0.25, log2FCThreshold=2., trivialAngle=NULL, 
	returnNormScores=FALSE, referenceCondition=NULL) 
	standardGeneric('compareSteadyNoNascent'))
#' Classify genes as delayed by the processing using the delta and tau metrics
setGeneric('processingDelay', function(inspectIds, tauThreshold=1.2,deltaThreshold=1.0, silent=TRUE) 
	standardGeneric('processingDelay'))
#' @rdname processingDelay
setGeneric('calculateDelta', function(inspectIds, silent=FALSE) 
	standardGeneric('calculateDelta'))
#' @rdname processingDelay
setGeneric('calculateTau', function(inspectIds, silent=FALSE) 
	standardGeneric('calculateTau'))

##########################################
# generics for class INSPEcT_diffsteady ####
##########################################

#' @rdname INSPEcT_diffsteady-class
setGeneric('synthesis', function(object) standardGeneric('synthesis'))
#' @rdname INSPEcT_diffsteady-class
setGeneric('processing', function(object) standardGeneric('processing'))
#' @rdname INSPEcT_diffsteady-class
setGeneric('degradation', function(object) standardGeneric('degradation'))

