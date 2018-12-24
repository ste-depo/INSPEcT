#' @rdname sfPlot
#'
#' @description
#' This method generates a plot that immediately shows the scaling factors used to scale
#' RNA- and Nascent-seq libraries and the possible relations between them. The ratio between the RNA- and the 
#' Nascent-seq scaling can be in fact considered as a yield of the synthesis within the cells.
#' @param object An object of class INSPEcT
#' @return None
#' nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
#' sfPlot(nascentInspObj10)
setMethod('sfPlot', 'INSPEcT', function(object) {
	tpts <- tpts(object)
	log_shift <- find_tt_par(tpts)
	x <- time_transf(tpts, log_shift)
	plot(x, labeledSF(object), type='b', lty=1, pch=1, xlab='time', ylab='labeled RNA scaling factor')
	axis(1, at=x, labels=signif(tpts, 2), las=3)
	})
