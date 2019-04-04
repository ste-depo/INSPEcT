convert_gene_classes <- function(gene_classes) {
	diz <- c('0'='KKK', 'a'='VKK', 'b'='KKV', 'c'='KVK',
		'ab'='VKV', 'ac'='VVK', 'bc'='KVV', 'abc'='VVV')
	unname(diz[gene_classes])
}

reconvert_gene_classes <- function(gene_classes) {
	diz <- c('KKK'='0','VKK'='a','KKV'='b','KVK'='c',
		'VKV'='ab','VVK'='ac','KVV'='bc','VVV'='abc')
	unname(diz[gene_classes])
}

###########################
## function for ranges ####
###########################

define_parameter_ranges <- function(ids, logshift, linshift) {

	range_k1_h_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate <- gene[[1]][['alpha']]
			switch(rate$type,
						 "constant" = rate$params[1],
						 "sigmoid" = rate$params[1:2],
						 "impulse" = rate$params[1:3]
						 )
			}))
		, probs=c(.025, .975))
	range_k1_h_pars <- c(
		floor(range_k1_h_pars[1]*100)/100,
		ceiling(range_k1_h_pars[2])
	)
	
	range_k2_h_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate <- gene[[1]][['gamma']]
			switch(rate$type,
						 "constant" = rate$params[1],
						 "sigmoid" = rate$params[1:2],
						 "impulse" = rate$params[1:3]
			)
		}))
	, probs=c(.025, .975))
	range_k2_h_pars <- c(
		floor(range_k2_h_pars[1]*100)/100,
		ceiling(range_k2_h_pars[2])
	)
	
	range_k3_h_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate <- gene[[1]][['beta']]
			switch(rate$type,
						 "constant" = rate$params[1],
						 "sigmoid" = rate$params[1:2],
						 "impulse" = rate$params[1:3]
			)
		}))
	, probs=c(.025, .975))
	range_k3_h_pars <- c(
		floor(range_k3_h_pars[1]*100)/100,
		ceiling(range_k3_h_pars[2])
	)
	
	range_t_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate_k1 <- gene[[1]][['alpha']]
			k1_t <- switch(rate_k1$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k1$params[3],
										 "impulse" = rate_k1$params[4:5]
			)
			rate_k2 <- gene[[1]][['beta']]
			k2_t <- switch(rate_k2$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k2$params[3],
										 "impulse" = rate_k2$params[4:5]
			)
			rate_k3 <- gene[[1]][['gamma']]
			k3_t <- switch(rate_k3$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k3$params[3],
										 "impulse" = rate_k3$params[4:5]
			)
			c(k1_t, k2_t, k3_t)
		}))
	, probs=c(.025, .975))
	range_t_pars <- time_transf_inv(range_t_pars, logshift, linshift)
	range_t_pars <- c(
		floor(range_t_pars[1]*100)/100, # (arrotonda per difetto al secondo decimale)
		ceiling(range_t_pars[2]*100)/100
	)
	
	range_beta_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate_k1 <- gene[[1]][['alpha']]
			k1_t <- switch(rate_k1$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k1$params[4],
										 "impulse" = rate_k1$params[6]
			)
			rate_k2 <- gene[[1]][['beta']]
			k2_t <- switch(rate_k2$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k2$params[4],
										 "impulse" = rate_k2$params[6]
			)
			rate_k3 <- gene[[1]][['gamma']]
			k3_t <- switch(rate_k3$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k3$params[4],
										 "impulse" = rate_k3$params[6]
			)
			c(k1_t, k2_t, k3_t)
		}))
	, probs=c(.025, .975))
	range_beta_pars <- c(
		floor(range_beta_pars[1]*100)/100,
		ceiling(range_beta_pars[2]*100)/100
	)

	return(list(
		k1_h_pars=range_k1_h_pars,
		k2_h_pars=range_k2_h_pars,
		k3_h_pars=range_k3_h_pars,
		t_pars=range_t_pars,
		beta_pars=range_beta_pars
		))

}

#######################
## PLOT FUNCTION ######
#######################

RNAdynamicsAppPlot <- function(data_selection, show_logtime, logshift, linshift,
	time_min, time_max, experiment, k1_function, k2_function, k3_function, 
	k1_params, k2_params, k3_params) {

	# get experimental values

	if( data_selection == 'Experimental data' ) {
		reference_mRNA <- experiment$mRNA
		secondary_mRNA <- experiment$mRNA_smooth
	} else { 
		reference_mRNA <- experiment$mRNA_smooth
		secondary_mRNA <- experiment$mRNA
	}
	if( data_selection == 'Experimental data' ) {
		reference_preMRNA <- experiment$preMRNA 
		secondary_preMRNA <- experiment$preMRNA_smooth
	} else {
		reference_preMRNA <- experiment$preMRNA_smooth
		secondary_preMRNA <- experiment$preMRNA
	}
	if( data_selection == 'Experimental data' ) {
		reference_synthesis <- experiment$synthesis
		secondary_synthesis <- experiment$synthesis_smooth
	} else {
		reference_synthesis <- experiment$synthesis_smooth
		secondary_synthesis <- experiment$synthesis
	}

	experimental_mRNAsd <- experiment$mRNAsd
	experimental_preMRNAsd <- experiment$preMRNAsd
	experimental_synthesissd <- experiment$synthesissd
	if( !experiment$steady_state ) {
		experiment_tpts <- experiment$tpts	
	} else {
		experiment_tpts <- 1:16
	}
	
	# make the simulation
	
	simulation_time <- seq(time_min,time_max,length.out=1000)
	simulation_time <- sort(unique(c(simulation_time, experiment_tpts)))
	transf_simulation_time <- time_transf(simulation_time, logshift)
	transf_tpts <- time_transf(experiment_tpts, logshift)
	
	sim <- deterministic_simulation(
		logshift, linshift, simulation_time, 
		k1_function, k2_function, k3_function, 
		k1_params, k2_params, k3_params)

	# function to define y-limits

	deltaylim <- function( yrange ) {
		deltarange <- yrange[2] * .05
		ylim <- yrange + c(-deltarange, deltarange)
	}

	# start plot routine
	
	par(mfrow=c(5,1))
	# layout(matrix(c(1,1,1,2,2,2,3,4,5), nrow=3, byrow=TRUE))
	par(mar=c(2.5,5,0,1)+.1)
	
	# plot k1

	plot_k1_experiment = ! (data_selection == 'User defined' | experiment$no_nascent)

	if( plot_k1_experiment ) {
		yrange <- range(c(sim[,'k1'], 
			c(secondary_synthesis + experimental_synthesissd, 
				reference_synthesis + experimental_synthesissd) , 
			c(secondary_synthesis - experimental_synthesissd, 
				reference_synthesis - experimental_synthesissd)))
		ylim <- deltaylim(yrange)
	} else {
		ylim <- deltaylim( range(sim[,'k1']) )
	}
	if( show_logtime ) {
		plot(transf_simulation_time, sim[,'k1'], 
			xaxs='i', yaxs='i', xaxt = 'n',
			ylab = 'synthesis', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(transf_simulation_time) 
				+ diff(range(transf_simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
		if( plot_k1_experiment ) {
			points( transf_tpts, secondary_synthesis, pch=1, col='grey')
			points( transf_tpts, reference_synthesis, pch=19)
			segments( transf_tpts , reference_synthesis - experimental_synthesissd 
				, transf_tpts , reference_synthesis + experimental_synthesissd )
		}
	} else {
		plot(simulation_time, sim[,'k1'], 
			xaxs='i', yaxs='i', xaxt = ifelse( plot_k1_experiment , 'n', 's'),
			ylab = 'synthesis', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(simulation_time) 
				+ diff(range(simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
		if( plot_k1_experiment ) {
			if( !experiment$steady_state ) {
				points( experiment_tpts, secondary_synthesis, pch=1, col='grey')
				points( experiment_tpts, reference_synthesis, pch=19)
				segments( experiment_tpts , reference_synthesis - experimental_synthesissd 
					, experiment_tpts , reference_synthesis + experimental_synthesissd )
			} else {
				points( experiment_tpts[1], secondary_synthesis, pch=1, col='grey')
				points( experiment_tpts[1], reference_synthesis, pch=19)
				segments( experiment_tpts[1] , reference_synthesis - experimental_synthesissd 
					, experiment_tpts[1] , reference_synthesis + experimental_synthesissd )				
			}
		}
	}

	# plot pre-RNA dynamics

	if( data_selection != 'User defined' ) {
		yrange <- range(c(sim[,'p'], 
			c(secondary_preMRNA + experimental_preMRNAsd, 
				reference_preMRNA + experimental_preMRNAsd) , 
			c(secondary_preMRNA - experimental_preMRNAsd, 
				reference_preMRNA - experimental_preMRNAsd)))
		ylim <- deltaylim( yrange )
	} else {
		ylim <- deltaylim( range(sim[,'p']) )
	}	
	if( show_logtime ) {
		plot(transf_simulation_time, sim[,'p'], 
				xaxs='i', yaxs='i', xaxt = 'n',
				ylab = 'pre-RNA', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
				xlim = range(transf_simulation_time) 
					+ diff(range(transf_simulation_time)) * c(-.05, .05),
				ylim = ylim
				)
		points( transf_tpts, secondary_preMRNA, pch=1, col='grey')
		points( transf_tpts, reference_preMRNA, pch=19)
		segments( transf_tpts , reference_preMRNA - experimental_preMRNAsd 
			, transf_tpts , reference_preMRNA + experimental_preMRNAsd )
	} else {
		plot(simulation_time, sim[,'p'], 
			xaxs='i', yaxs='i', xaxt = ifelse( data_selection != 'User defined' , 'n', 's'),
			ylab = 'pre-RNA', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(simulation_time) 
				+ diff(range(simulation_time)) * c(-.05, .05),
			ylim = ylim)
		if( data_selection != 'User defined' ) {
			if( !experiment$steady_state ) {
				points( experiment_tpts, secondary_preMRNA, pch=1, col='grey')
				points( experiment_tpts, reference_preMRNA, pch=19)
				segments( experiment_tpts , reference_preMRNA - experimental_preMRNAsd 
					, experiment_tpts , reference_preMRNA + experimental_preMRNAsd )
			} else {
				points( experiment_tpts[1], secondary_preMRNA, pch=1, col='grey')
				points( experiment_tpts[1], reference_preMRNA, pch=19)
				segments( experiment_tpts[1] , reference_preMRNA - experimental_preMRNAsd 
					, experiment_tpts[1] , reference_preMRNA + experimental_preMRNAsd )								
			}
		}
	}

	# plot k2

	ylim <- deltaylim( range(sim[,'k2']) )
	
	if( show_logtime ) {
	plot(transf_simulation_time, sim[,'k2'], 
			xaxs='i', yaxs='i', xaxt = 'n',
			ylab = 'processing', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(transf_simulation_time) 
				+ diff(range(transf_simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
		# axis(1, at=transf_tpts, labels = signif(experiment_tpts,2) )
	} else {
	plot(simulation_time, sim[,'k2'], 
			xaxs='i', yaxs='i', xaxt = ifelse( data_selection != 'User defined' , 'n', 's'),
			ylab = 'processing', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(simulation_time) 
				+ diff(range(simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
	# if( data_selection != 'User defined' )
		# axis(1, at=experiment_tpts, labels = signif(experiment_tpts,2) )
	}
	
	# plot mRNA dynamics 
	
	if( data_selection != 'User defined' ) {
		yrange <- range(c(sim[,'m'], 
			c(secondary_mRNA + experimental_mRNAsd, 
				reference_mRNA + experimental_mRNAsd) , 
			c(secondary_mRNA - experimental_mRNAsd, 
				reference_mRNA - experimental_mRNAsd)))
		ylim <- deltaylim( yrange )
	} else {
		ylim <- deltaylim( range(sim[,'m']) )
	}
	if( show_logtime ) {
		plot(transf_simulation_time, sim[,'m'], 
			xaxs='i', yaxs='i', xaxt = 'n',
			ylab = 'mature RNA', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(transf_simulation_time) 
				+ diff(range(transf_simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
		points( transf_tpts, secondary_mRNA, pch=1, col='grey')
		points( transf_tpts, reference_mRNA, pch=19)
		segments( transf_tpts , reference_mRNA - experimental_mRNAsd 
			, transf_tpts , reference_mRNA + experimental_mRNAsd )
		# axis(1, at=transf_tpts, labels = signif(experiment_tpts,2) )
	} else {
		plot(simulation_time, sim[,'m'], 
			xaxs='i', yaxs='i', xaxt = ifelse( data_selection != 'User defined' , 'n', 's'),
			ylab = 'mature RNA', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(simulation_time) 
				+ diff(range(simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
		if( data_selection != 'User defined' ) {
			if( !experiment$steady_state ) {
				points( experiment_tpts, secondary_mRNA, pch=1, col='grey')
				points( experiment_tpts, reference_mRNA, pch=19)
				segments( experiment_tpts , reference_mRNA - experimental_mRNAsd 
					, experiment_tpts , reference_mRNA + experimental_mRNAsd )
			} else {
				points( experiment_tpts[1], secondary_mRNA, pch=1, col='grey')
				points( experiment_tpts[1], reference_mRNA, pch=19)
				segments( experiment_tpts[1] , reference_mRNA - experimental_mRNAsd 
					, experiment_tpts[1] , reference_mRNA + experimental_mRNAsd )				
			}
			# axis(1, at=experiment_tpts, labels = signif(experiment_tpts,2) )
		}
	}
	
	# plot k3

	ylim <- deltaylim( range(sim[,'k3']) )

	if( show_logtime ) {
	plot(transf_simulation_time, sim[,'k3'], 
			xaxs='i', yaxs='i', xaxt = 'n',
			ylab = 'degradation', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(transf_simulation_time) 
				+ diff(range(transf_simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
		axis(1, at=transf_tpts, labels = signif(experiment_tpts,2) , cex.axis = 1.3)
	} else {
	plot(simulation_time, sim[,'k3'], 
			xaxs='i', yaxs='i', xaxt = ifelse( data_selection != 'User defined' , 'n', 's'),
			ylab = 'degradation', type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			xlim = range(simulation_time) 
				+ diff(range(simulation_time)) * c(-.05, .05),
			ylim = ylim
			)
	if( data_selection != 'User defined' )
		axis(1, at=experiment_tpts, labels = signif(experiment_tpts,2) , cex.axis = 1.3)
	}
			
	# calculate the scores of the modeling and assign to output
	
	if( data_selection != 'User defined' & !experiment$steady_state ) {
		
		scores <- list()
		
		loglik_score <- 
			
			logLikelihoodFunction(reference_preMRNA, 
								sim[simulation_time %in% experiment_tpts,'p'], 
								experimental_preMRNAsd^2) +
			
			logLikelihoodFunction(reference_mRNA, 
								sim[simulation_time %in% experiment_tpts,'m'], 
								experimental_mRNAsd^2) +
			
			ifelse(experiment$no_nascent, 0, logLikelihoodFunction(reference_synthesis, 
								sim[simulation_time %in% experiment_tpts,'k1'], 
								experimental_synthesissd^2))
		
		chisq_score <- 
			
			chisqFunction(reference_preMRNA, 
								sim[simulation_time %in% experiment_tpts,'p'], 
								experimental_preMRNAsd^2) +

			chisqFunction(reference_mRNA, 
								sim[simulation_time %in% experiment_tpts,'m'], 
								experimental_mRNAsd^2) +

			ifelse(experiment$no_nascent, 0, chisqFunction(reference_synthesis, 
								sim[simulation_time %in% experiment_tpts,'k1'], 
								experimental_synthesissd^2))
		
		k <- length(c(ifelse(experiment$no_nascent, 0, k1_params), k2_params, k3_params))

		scores$pchisq <- pchisq( chisq_score, 3*length(experiment_tpts) - k )
		scores$aic <- 2*k - 2*loglik_score
		
		#print(scores)
		
	} else {
		
		scores <- list()
		scores$pchisq <- NA
		scores$aic <- NA
		
	}
	
	return(scores)
	
}

constantModelRNApp <- function(x , par, log_shift, lin_shift) rep(par, length(x))

sigmoidModelRNApp <- function(x, par, log_shift, lin_shift=0) 
	par[1]+(par[2]-par[1])*(1/(1+exp(-par[4]*(time_transf(x,log_shift,lin_shift)-time_transf(par[3],log_shift,lin_shift)))))

impulseModelRNApp <- function(x, par, log_shift, lin_shift=0) 
	1/par[2]*(par[1]+(par[2]-par[1])*(1/(1+exp(-par[6]*(time_transf(x,log_shift,lin_shift)-time_transf(par[4],log_shift,lin_shift))))))*
	(par[3]+(par[2]-par[3])*(1/(1+exp(par[6]*(time_transf(x,log_shift,lin_shift)-time_transf(par[5],log_shift,lin_shift))))))


#############################
## SIMULATION FUNCTION ######
#############################

rxnrateMatureRNA <- function(t,c,parms){
	alpha <- parms$alpha
	beta  <- parms$beta
	gamma <- parms$gamma
	r=rep(0,length(c))
	r[1] <- alpha(t) - gamma(t) * c["p"]
	r[2] <- gamma(t) * c["p"] - beta(t) * c["m"]
	return(list(r))
}

deterministic_simulation <- function(log_shift, lin_shift, simulation_time,
	 k1_function, k2_function, k3_function,
	 k1_params, k2_params, k3_params
	 )
# in this version of the function, 
# all three rates are modeled as impulse models
{
		
	tpts <- simulation_time
	
	params <- list(
		alpha = function(x) 
			switch(k1_function, 
						 "Constant" = constantModelRNApp(x, k1_params),
						 "Sigmoidal" = sigmoidModelRNApp(x, k1_params, log_shift, lin_shift),
						 "Impulsive" = impulseModelRNApp(x, k1_params, log_shift, lin_shift)
						 )
		,
		beta = function(x) 
			switch(k3_function, 
						 "Constant" = constantModelRNApp(x, k3_params),
						 "Sigmoidal" = sigmoidModelRNApp(x, k3_params, log_shift, lin_shift),
						 "Impulsive" = impulseModelRNApp(x, k3_params, log_shift, lin_shift)
			)
		,
		gamma = function(x) 
			switch(k2_function, 
						 "Constant" = constantModelRNApp(x, k2_params),
						 "Sigmoidal" = sigmoidModelRNApp(x, k2_params, log_shift, lin_shift),
						 "Impulsive" = impulseModelRNApp(x, k2_params, log_shift, lin_shift)
			)
	)
	
	cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1])
						 , params$alpha(tpts[1]) / params$beta(tpts[1]))
	names(cinit) <- c('p', 'm')
	model <- ode(y=cinit, times=tpts, func=rxnrateMatureRNA, parms=params)
	#model[,2:3] <- t(t(model[,2:3])/model[1,2:3])
	synthesis <- params$alpha(tpts)
	processing <- params$gamma(tpts)
	degradation <- params$beta(tpts)
	model <- cbind(time=model[,1]
		, k1=synthesis
		, k2=processing
		, k3=degradation
		, p=model[,2]
		, m=model[,3]
		)
	return(model)
}

smoothModel <- function(tpts, experiment, nInit=10, nIter=500, seed=1234)
{
	
	if( !is.null(seed) ) set.seed(seed)

	optimFailOut <- function(e) list(par=NA, value=NA, counts=NA, convergence=1, message=e)

	im.parguess <- function(tpts , values, log_shift ) {
		tpts <- time_transf(tpts, log_shift)
		ntp   <- length(tpts)
		peaks <- which(diff(sign(diff(values)))!=0)+1
		if( length(peaks) == 1 ) peak <- peaks
		if( length(peaks)  > 1 ) peak <- sample(peaks, 1)
		if( length(peaks) == 0 ) peak <- round(length(tpts)/2)
		#
		initial_values <- runif( 1, min=min(values[1:3])
														 , max=max(values[1:3]))
		intermediate_values <- values[peak]
		if( intermediate_values==0 ) intermediate_values <- mean(values[seq(peak-1,peak+1)])
		end_values <- runif( 1, min=min(values[(ntp-2):ntp])
												 , max=max(values[(ntp-2):ntp]))
		time_of_first_response  <- tpts[peak-1]
		time_of_second_response <- tpts[peak+1]
		slope_of_response <- diff(range(tpts)) / 
			(time_of_second_response-time_of_first_response)
		#
		return(c(h0=initial_values, h1=intermediate_values
						 , h2=end_values, t1=time_of_first_response
						 , t2=time_of_second_response, b=slope_of_response))
	}

	im.chisq <- function(par, tpts, experiment, log_shift) 
	{
		model <- impulseModelRNApp(tpts, par, log_shift)
		chisqFunction(experiment, model, 1)
	}
	
	log_shift <- find_tt_par(tpts)
	
	outIM <- sapply(1:nInit, function(x) 
			suppressWarnings(tryCatch(optim(
				par=im.parguess(tpts, experiment, log_shift)
				, fn=im.chisq, tpts=tpts
				, experiment=experiment
				, log_shift=log_shift
				, control=list(maxit=nIter)
			), error=function(e) optimFailOut(e))))
	
	bestIM <- which.min(unlist(outIM[2,]))
	impulseModelRNApp( tpts, outIM[,bestIM]$par, log_shift)
	
}

###############################
## MINIMIZATION FUNCTION ######
###############################

modelChisqMatureRNA <- function(par, tpts, fun, df, alpha_exp, alpha_var #, pval
	, mature_exp, mature_var, preMRNA_exp, preMRNA_var, log_shift, lin_shift, no_nascent)
{
	splitpar <- split(par 
		, c(rep('alpha',df[1]), rep('beta',df[2]), rep('gamma',df[3])) 
		)
	#

	params <- list()
	params$alpha <- function(x) 
		fun$alpha$value(x, splitpar$alpha, log_shift, lin_shift)
	params$beta  <- function(x)
		fun$beta$value(x, splitpar$beta, log_shift, lin_shift)
	params$gamma <- function(x)
		fun$gamma$value(x, splitpar$gamma, log_shift, lin_shift)
	#
	cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1])
		, params$alpha(tpts[1]) / params$beta(tpts[1]))
	names(cinit) <- c('p', 'm')
	model <- ode(y=cinit, times=tpts, func=rxnrateMatureRNA, parms=params)
	#
	alpha_model <- params$alpha(tpts)
	matureodel <- model[,'m']
	preMRNA_model <- model[,'p']
	#
	ifelse(no_nascent, 0, chisqFunction(alpha_exp, alpha_model, alpha_var)) +
		chisqFunction(mature_exp, matureodel, mature_var) +
		chisqFunction(preMRNA_exp, preMRNA_model, preMRNA_var)
}

optimParamsMatureRNA <- function(interpRates, tpts_exp, alpha_exp, alpha_var, mature_exp
	, mature_var, preMRNA_exp, preMRNA_var, maxit=500, method='NM', log_shift, lin_shift, no_nascent)
{
	if( method == 'NM' ) method <- 'Nelder-Mead' 
	tryCatch({
		optOut <- optim(
			par=unlist(sapply(interpRates, '[[', 'params'))
			, fn=modelChisqMatureRNA , tpts=tpts_exp
			, fun=sapply(interpRates, '[[', 'fun')
			, df=unlist(sapply(interpRates, '[[', 'df'))
			, alpha_exp=alpha_exp, mature_exp=mature_exp
			, preMRNA_exp=preMRNA_exp
			, alpha_var=alpha_var, mature_var=mature_var
			, preMRNA_var=preMRNA_var
			, log_shift = log_shift
			, lin_shift = lin_shift
			, no_nascent = no_nascent
			, control = list(maxit = maxit)
			, method = method
			)
		splitpar <- split(optOut$par
			, c(rep('alpha',interpRates$alpha$df)
				, rep('beta',interpRates$beta$df)
				, rep('gamma',interpRates$gamma$df))
			)
		interpRates$alpha$params <- splitpar$alpha
		interpRates$beta$params  <- splitpar$beta
		interpRates$gamma$params <- splitpar$gamma
		return(list(
			alpha=interpRates$alpha
			, beta=interpRates$beta
			, gamma=interpRates$gamma
			, counts=optOut$counts[1]
			, convergence=optOut$convergence
			, message=optOut$message
			))}
		, error=function(e) {
			print(e)
			return(list(
				alpha=interpRates$alpha
				, beta=interpRates$beta
				, gamma=interpRates$gamma
				, counts=0
				, convergence=1
				, message=e
			))
			})
}



