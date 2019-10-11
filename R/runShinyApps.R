#' Run shiny applications contained in the package INSPEcT
#' @description
#' Two shiny apps are encoded into the package inspect:
#' - runProcessingRateDelay: plots single genes as well as genome wide plots associated to the
#'   processing induced delay loading the data from an INSPEcT object.
#' - runINSPEcT_GUI: is a way to visualize and interact with the RNAdynamics at the level of 
#'   a single gene, either loading the data from ad INSPEcT object or from scratch.
#' @docType methods
#' @name INSPEcT-shinyApps
NULL

#' @rdname INSPEcT-shinyApps
#' @examples
#' # runProcessingRateDelay()
runProcessingRateDelay = function() {
	shiny::runApp(system.file(package='INSPEcT', 'ProcessingRateDelay'))
}
#' @rdname INSPEcT-shinyApps
#' @examples
#' # runINSPEcT_GUI()
runINSPEcT_GUI = function() {
	shiny::runApp(system.file(package='INSPEcT', 'RNAdynamics'))
}
