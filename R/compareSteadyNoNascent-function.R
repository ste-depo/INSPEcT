compareSteadyNoNascent <- function(expressions,err)
{

	eiGenes <- intersect(rownames(expressions$exonsExpressions),rownames(expressions$intronsExpressions))

	# Mature, premature and total rpkms
	mature <- expressions$exonsExpressions[eiGenes,] - expressions$intronsExpressions[eiGenes,]
	premature <- expressions$intronsExpressions[eiGenes,]
	total <- expressions$exonsExpressions[eiGenes,]

	prematureMedian <- apply(premature,1,function(r)median(r,na.rm=T))
	matureMedian <- apply(mature,1,function(r)median(r,na.rm=T))
	
	standardCurveFit <- standardCurveFitFunction(prematureMedian,matureMedian)
	classificationTmp <- classificationFunction(premature,mature,standardCurveFit)


	return(classificationTmp)
}

standardCurveFitFunction <- function(p,m,err)
{
	n_outliers <- function(alpha, x, y, err) {

		#Conversion
		pi_angle <- alpha * pi/180
		#Angular coefficient
		coef_ang <- tan(pi_angle)

		delta_intercept <- err/cos(pi_angle)
		intercept <- median(y,na.rm=TRUE) - coef_ang*median(x,na.rm=TRUE)

		outliers <- y > coef_ang * x + intercept + delta_intercept |
			y < coef_ang * x + intercept - delta_intercept

		length(which(outliers))
	}

	all_alphas <- seq(-89,90)
	all_alphas_outliers <- sapply(all_alphas, n_outliers, x=log2(p), y=log2(m))
	return(seq(-89,90)[which.min(all_alphas_outliers)])
}

classificationFunction <- function(p,m,alpha,err)
{
	classificationTmp <- sapply(rownames(p),function(g)
	{
		x <- log2(p[g,])
		y <- log2(m[g,])

		pi_angle <- standardCurveFit * pi/180
		coef_ang <- tan(pi_angle)
		delta_intercept <- err/cos(pi_angle)

		intercept <- median(y,na.rm=TRUE) - coef_ang*median(x,na.rm=TRUE)

		outliers <- y > coef_ang * x + intercept + delta_intercept |
				y < coef_ang * x + intercept - delta_intercept
		return(outliers)
	})

	return(t(classificationTmp))
}