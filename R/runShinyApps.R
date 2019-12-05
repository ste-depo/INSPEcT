#' Run shiny applications contained in the package INSPEcT
#' @description
#' Two shiny apps are encoded into the package inspect:
#' - runProcessingRateDelay: plots single genes as well as genome wide plots associated to the
#'   processing induced delay loading the data from an INSPEcT object.
#' - runINSPEcTGUI: is a way to visualize and interact with the RNAdynamics at the level of 
#'   a single gene, either loading the data from ad INSPEcT object or from scratch.
#' @docType methods
#' @name INSPEcT-shinyApps
NULL

#' @rdname INSPEcT-shinyApps
runProcessingRateDelay = function() {
	shinyApp(ui = ProcessingRateDelayshinyAppUI, server = ProcessingRateDelayshinyAppServer)
}

#' @rdname INSPEcT-shinyApps
runINSPEcTGUI <- function() {
	shinyApp(ui = INSPEcTGUIshinyAppUI, server = INSPEcTGUIshinyAppServer)
}
