#' @rdname correlationPlot
#'
#' @description
#' This function plot the rates of a simulated dataset against the modeled ones and compute their correlations. 
#' @param object An object of class INSPEcT_model with simulated rates.
#' @param object2 An object of class INSPEcT.
#' @param plot A logical indicating whether to draw or not the plot. (default=TRUE)
#' @return An list with the correlation values.

setMethod('correlationPlot', signature(object='INSPEcT_model', object2='INSPEcT'), function(object, object2, plot=TRUE)
{
	k1_real <- log10(sapply(object@ratesSpecs,function(g){g[[1]][["alpha"]]$fun$value(0,g[[1]][["alpha"]]$par)}))
	k2_real <- log10(sapply(object@ratesSpecs,function(g){g[[1]][["gamma"]]$fun$value(0,g[[1]][["gamma"]]$par)}))
	k3_real <- log10(sapply(object@ratesSpecs,function(g){g[[1]][["beta"]]$fun$value(0,g[[1]][["beta"]]$par)}))

	k1_modeled <- log10(viewModelRates(object2,"synthesis")[,1])
	k2_modeled <- log10(viewModelRates(object2,"processing")[,1])
	k3_modeled <- log10(viewModelRates(object2,"degradation")[,1])

	k1_cor <- cor(k1_real,k1_modeled,method="s",use="c")
	k2_cor <- cor(k2_real,k2_modeled,method="s",use="c")
	k3_cor <- cor(k3_real,k3_modeled,method="s",use="c")

	if(plot)
	{
		par(mfrow=c(1,3))

		x <- log10(k1_real)
		y <- log10(k1_modeled)
		smoothScatter(x,y,xlab="Log10(Real)",ylab="Log10(Modeled)",main="Synthesis rate",xlim=c(min(c(x,y),na.rm=TRUE),max(c(x,y),na.rm=TRUE)), ylim=c(min(c(x,y),na.rm=TRUE),max(c(x,y),na.rm=TRUE)))
		abline(0,1,col=2,lwd=2,lty=3)
		
		legend("topleft",legend=paste0("Rho coefficient: ",round(k1_cor,2)),border=FALSE,bty = "n")

		x <- log10(k2_real)
		y <- log10(k2_modeled)
		smoothScatter(x,y,xlab="Log10(Real)",ylab="Log10(Modeled)",main="Processing rate",xlim=c(min(c(x,y),na.rm=TRUE),max(c(x,y),na.rm=TRUE)), ylim=c(min(c(x,y),na.rm=TRUE),max(c(x,y),na.rm=TRUE)))
		abline(0,1,col=2,lwd=2,lty=3)
		
		legend("topleft",legend=paste0("Rho coefficient: ",round(k2_cor,2)),border=FALSE,bty = "n")

		x <- log10(k3_real)
		y <- log10(k3_modeled)
		smoothScatter(x,y,xlab="Log10(Real)",ylab="Log10(Modeled)",main="Degradation rate",xlim=c(min(c(x,y),na.rm=TRUE),max(c(x,y),na.rm=TRUE)), ylim=c(min(c(x,y),na.rm=TRUE),max(c(x,y),na.rm=TRUE)))
		abline(0,1,col=2,lwd=2,lty=3)
		
		legend("topleft",legend=paste0("Rho coefficient: ",round(k3_cor,2)),border=FALSE,bty = "n")		
	}

	return(list("synthesis"=k1_cor,"processing"=k2_cor,"degradation"=k3_cor))
})