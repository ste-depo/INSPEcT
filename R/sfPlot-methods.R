#' @rdname sfPlot
#'
#' @description
#' This method generates a plot that immediately shows the scaling factors used to scale
#' RNA- and Nascent-seq libraries and the possible relations between them. The ratio between the RNA- and the 
#' Nascent-seq scaling can be in fact considered as a yield of the synthesis within the cells.
#' @param object An object of class INSPEcT
#' @return None
#' @examples
#' data('mycerIds10', package='INSPEcT')
#' sfPlot(mycerIds10)
setMethod('sfPlot', 'INSPEcT', function(object) {
	mat <- cbind(object@totalSF, object@labeledSF)
	matplot(mat, type='l', lty=1, xlab='time index', ylab='scaling factor')
	legend('bottomleft', legend=c('total','Nascent'), col=1:2, lty=1)
	})
