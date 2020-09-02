#' @rdname plotGene
#'
#' @description
#' A method to see the shapes of the estimated synthesis, degradation and processing rates, pre-RNA and total RNA 
#' concentrations (solid thin lines) their variances (dashed lines) and the modeled rates and concentrations 
#' (ticker solid line) of a single gene. 
#' @param object An object of class INSPEcT
#' @param ix Eiher a rowname or a row number to select one single gene
#' @param relative_expression A logical, indicating whether expressions are rates 
#' should be plotted relative to their initial value (Default=FALSE).
#' @param fix.yaxis A logical, indicating whether the limits for y-axis of 
#' degradation and processing rates should be fixed
#' relative to their distributions
#' @param priors A logical, if true the priors of the rates are plotted
#' @param constantModel A logical, if true the constant model for the + nascent modeling are shown
#' @return A list containing total RNA levels and their confidence interval (levels plus and minus
#' one standard deviation), pre-RNA lelevs and their confidence intervals, synthsis rates and 
#' their confidence intervals, degradation rates and processing rates of the selected gene.
#' @examples
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' plotGene(nascentInspObj10, 1)
setMethod('plotGene', 'INSPEcT', function(object, ix, relative_expression=FALSE, fix.yaxis=FALSE, 
                                          priors=TRUE, constantModel=FALSE) {
  checkINSPEcTObjectversion(object)
  
  oldMfrow <- par()$mfrow
  oldMar <- par()$mar
  
  ix <- ix[1]
  tpts <- object@tpts
  oneGene <- object[ix]
  
  if( nrow(oneGene@modelRates) > 0 ) {
    
    modeled_rates <- TRUE
    
    if( oneGene@NoNascent & !oneGene@NF )
      foe <- capture.output(oneGene <- computeConfidenceIntervals(oneGene))
    pValues <- formatC(as.numeric(ratePvals(object)[ix,]),format = "e", digits = 1)
    names(pValues) <- c('synthesis','processing','degradation')
    
  } else {
    
    modeled_rates <- FALSE
    
    pValues <- rep(NA, 3)
    names(pValues) <- c('synthesis','processing','degradation')
    
  }
  
  if( fix.yaxis ) {
    degradationYlim <- quantile(ratesFirstGuess(object, 'degradation'), 
                                probs=c(.02, .98), na.rm=TRUE)
    processingYlim <- quantile(ratesFirstGuess(object, 'processing'), 
                               probs=c(.02, .98), na.rm=TRUE)
  } else {
    processingYlim <- degradationYlim <- NULL
  }
  
  x <- seq_along(tpts)
  par(mfrow=c(1,5), mar=.1+c(5,2.5,4,.5))
  
  plotSingleRNADynamic(dyn_name = 'total-RNA', tag = '', 
                       simtimeplot = x, 
                       simprofile = if(modeled_rates) {
                         c(viewModelRates(oneGene, 'total'))
                       } else rep(NA, length(x)), 
                       ci_left = rep(NA, length(x)), 
                       ci_right = rep(NA, length(x)), 
                       plot_exp = TRUE, exptimeplot = x, 
                       ref_exp = c(ratesFirstGuess(oneGene, 'total')), 
                       sec_exp = rep(NA, length(x)), 
                       ssd_exp = sqrt(c(ratesFirstGuessVar(oneGene, 'total'))), 
                       show_relexpr = relative_expression, ylim = NULL, rate_p = NULL, 
                       command_line = TRUE, col = 1,
                       prior = if(!modeled_rates) {
                         c(ratesFirstGuess(oneGene, 'total'))
                       } else NULL)
  
  axis(1, at=x, labels=signif(tpts, 2), las=3)
  title('total-RNA')
  
  if( any(!is.na(ratesFirstGuess(oneGene, 'preMRNA'))) )
  {

    plotSingleRNADynamic(dyn_name = 'pre-RNA', tag = '', 
                         simtimeplot = x, 
                         simprofile = if(modeled_rates) {
                           c(viewModelRates(oneGene, 'preMRNA'))
                         } else rep(NA, length(x)), 
                         ci_left = rep(NA, length(x)), 
                         ci_right = rep(NA, length(x)), 
                         plot_exp = TRUE, exptimeplot = x, 
                         ref_exp = c(ratesFirstGuess(oneGene, 'preMRNA')), 
                         sec_exp = rep(NA, length(x)), 
                         ssd_exp = sqrt(c(ratesFirstGuessVar(oneGene, 'preMRNA'))), 
                         show_relexpr = relative_expression, ylim = NULL, rate_p = NULL, 
                         command_line = TRUE, col = 2, 
                         prior = if(!modeled_rates) {
                           c(ratesFirstGuess(oneGene, 'preMRNA'))
                         } else NULL)
    
    axis(1, at=x, labels=signif(tpts, 2), las=3)
    title('pre-RNA')
    
  } else {
    
    plot(x, rep(1, length(x)), pch='', main='-', xaxt='n', xlab='time', ylab='')
    axis(1, at=x, labels=signif(tpts, 2), las=3)
    
  }
  
  plotSingleRNADynamic(dyn_name = 'synthesis', tag = 's', 
                       simtimeplot = x, 
                       simprofile = if(modeled_rates) {
                         c(viewModelRates(oneGene, 'synthesis'))
                       } else rep(NA, length(x)), 
                       ci_left = if(modeled_rates) {
                         c(viewConfidenceIntervals(oneGene, 'synthesis_left'))
                       } else rep(NA, length(x)), 
                       ci_right = if(modeled_rates) {
                         c(viewConfidenceIntervals(oneGene, 'synthesis_right'))
                       } else rep(NA, length(x)), 
                       plot_exp = TRUE, exptimeplot = x, 
                       ref_exp = c(ratesFirstGuess(oneGene, 'synthesis')), 
                       sec_exp = rep(NA, length(x)), 
                       ssd_exp = sqrt(c(ratesFirstGuessVar(oneGene, 'synthesis'))), 
                       show_relexpr = relative_expression, ylim = NULL, rate_p = NULL, 
                       command_line = TRUE, col = 3,
                       prior = if(!modeled_rates) {
                         c(ratesFirstGuess(oneGene, 'synthesis'))
                       } else NULL,
                       constant = if( constantModel & !oneGene@NoNascent ) {
                         c(viewConfidenceIntervals(oneGene,"synthesis_constant"))
                       } else NULL)

  axis(1, at=x, labels=signif(tpts, 2), las=3)
  title(paste0('synthesis\n',pValues['synthesis']))
  
  if( any(!is.na(ratesFirstGuess(oneGene, 'processing'))) )
  {
    plotSingleRNADynamic(dyn_name = 'processing', tag = 'p', 
                         simtimeplot = x, 
                         simprofile = if(modeled_rates) {
                           c(viewModelRates(oneGene, 'processing'))
                         } else rep(NA, length(x)), 
                         ci_left = if(modeled_rates) {
                           c(viewConfidenceIntervals(oneGene, 'processing_left'))
                         } else rep(NA, length(x)), 
                         ci_right = if(modeled_rates) {
                           c(viewConfidenceIntervals(oneGene, 'processing_right'))
                         } else rep(NA, length(x)), 
                         plot_exp = FALSE,
                         show_relexpr = relative_expression, ylim = processingYlim, rate_p = NULL, 
                         command_line = TRUE, col = 5, 
                         prior = if(!modeled_rates | priors) {
                           c(ratesFirstGuess(oneGene, 'processing'))
                         } else NULL,
                         constant = if( constantModel & !oneGene@NoNascent ) {
                           c(viewConfidenceIntervals(oneGene,"processing_constant"))
                         } else NULL)
    axis(1, at=x, labels=signif(tpts, 2), las=3)
    title(paste0('processing\n',pValues['processing']))
  } else {
    plot(x, rep(1, length(x)), pch='', main='-', xaxt='n', xlab='time', ylab='')
    axis(1, at=x, labels=signif(tpts, 2), las=3)
  }
  
  plotSingleRNADynamic(dyn_name = 'degradation', tag = 'p', 
                       simtimeplot = x, 
                       simprofile = if(modeled_rates) {
                         c(viewModelRates(oneGene, 'degradation'))
                       } else rep(NA, length(x)), 
                       ci_left = if(modeled_rates) {
                         c(viewConfidenceIntervals(oneGene, 'degradation_left'))
                       } else rep(NA, length(x)), 
                       ci_right = if(modeled_rates) {
                         c(viewConfidenceIntervals(oneGene, 'degradation_right'))
                       } else rep(NA, length(x)), 
                       plot_exp = FALSE,
                       show_relexpr = relative_expression, ylim = degradationYlim, rate_p = NULL, 
                       command_line = TRUE, col = 4, 
                       prior = if(!modeled_rates | priors) {
                         c(ratesFirstGuess(oneGene, 'degradation'))
                       } else NULL,
                       constant = if( constantModel & !oneGene@NoNascent ) {
                         c(viewConfidenceIntervals(oneGene,"degradation_constant"))
                       } else NULL)
  axis(1, at=x, labels=signif(tpts, 2), las=3)
  title(paste0('degradation\n',pValues['degradation']))
  
  par(mfrow=oldMfrow, mar=oldMar)

})
