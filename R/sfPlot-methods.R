#' @rdname sfPlot
#'
#' @description
#' This method generates a plot that immediately shows the scaling factors used to scale
#' RNA- and 4sU-seq libraries and the possible relations between them. The ratio between the RNA- and the 
#' 4sU-seq scaling can be in fact considered as a yield of the synthesis within the cells.
#' @param object An object of class INSPEcT
#' @return None
#' @examples
#' data('rpkms', package='INSPEcT')
#' tpts <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
#' tL <- 1/6
#' mycerIds <- newINSPEcT(tpts, tL, rpkms$foursu_exons, rpkms$total_exons, 
#' 	rpkms$foursu_introns, rpkms$total_introns, BPPARAM=SerialParam())
#' sfPlot(mycerIds)
setMethod('sfPlot', 'INSPEcT', function(object) {
	mat <- cbind(object@totalSF, object@labeledSF)
	matplot(mat, type='l', lty=1, xlab='time index', ylab='scaling factor')
	legend('bottomleft', legend=c('total','4sU'), col=1:2, lty=1)
	})
