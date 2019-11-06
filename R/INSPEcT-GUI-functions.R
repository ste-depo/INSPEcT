#function for the indicator "Loading...". This function was taken from
#xiaodaigh and dcurrier (https://github.com/AnalytixWare/ShinySky/blob/master/R/busy-indicator.r)
#and corrected. Maybe the real time is inside setInterval function
.busyIndicator <- function(text = "Processing..."
													 , image = "http://i.giphy.com/l3V0EQrPMh1nnfbFe.gif"
													 , wait=1000) {
	tagList(
		singleton(tags$head(
			tags$link(rel = "stylesheet"
								, type = "text/css" 
								,href = file.path("panel","inst","extdata","busyIndicator.css")
			)))
		,div(class = "mybusyindicator",p(text)) #,img(src=image))
		,tags$script(sprintf(
			" setInterval(function(){
			if ($('html').hasClass('shiny-busy')) {
				setTimeout(function() {
					if ($('html').hasClass('shiny-busy')) {
						$('div.mybusyindicator').show()
					}
				}, %d)          
			} else {
				$('div.mybusyindicator').hide()
			}
		},1000)",wait)
		)
	) 
}

##########################

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

define_parameter_ranges <- function(ids) {

	range_k1_h_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate <- gene[[1]][[ 1 ]]
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
			rate <- gene[[1]][[ 2 ]]
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
			rate <- gene[[1]][[ 3 ]]
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
			rate_k1 <- gene[[1]][[ 1 ]]
			k1_t <- switch(rate_k1$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k1$params[3],
										 "impulse" = rate_k1$params[4:5]
			)
			rate_k2 <- gene[[1]][[ 2 ]]
			k2_t <- switch(rate_k2$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k2$params[3],
										 "impulse" = rate_k2$params[4:5]
			)
			rate_k3 <- gene[[1]][[ 3 ]]
			k3_t <- switch(rate_k3$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k3$params[3],
										 "impulse" = rate_k3$params[4:5]
			)
			c(k1_t, k2_t, k3_t)
		}))
	, probs=c(.025, .975))
	# range_t_pars <- timetransf_inv(range_t_pars, logshift, linshift)
	range_t_pars <- c(
		floor(range_t_pars[1]*100)/100, # (arrotonda per difetto al secondo decimale)
		ceiling(range_t_pars[2]*100)/100
	)
	
	range_beta_pars <- quantile(
		unlist(lapply(ids@model@ratesSpecs, function(gene) {
			rate_k1 <- gene[[1]][[ 1 ]]
			k1_t <- switch(rate_k1$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k1$params[4],
										 "impulse" = rate_k1$params[6]
			)
			rate_k2 <- gene[[1]][[ 2 ]]
			k2_t <- switch(rate_k2$type,
										 "constant" = NULL,
										 "sigmoid" = rate_k2$params[4],
										 "impulse" = rate_k2$params[6]
			)
			rate_k3 <- gene[[1]][[ 3 ]]
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

RNAdynamicsAppMake <- function(data_selection,
															 time_min, time_max, experiment,
															 k1_function, k2_function,
															 k3_function, k1_params,
															 k2_params, k3_params, 
															 mod_method) {

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
		simulation_time <- seq(time_min,time_max,length.out=1000)
		simulation_time <- sort(unique(c(simulation_time, experiment_tpts)))
	} else {
		experiment_tpts <- 0
		simulation_time <- seq(0,16,length.out=1000)
	}
	
	if( mod_method == 'int' ) {
		
		sim <- deterministic_simulation(
			simulation_time, 
			k1_function, k2_function, k3_function, 
			k1_params, k2_params, k3_params)

	} else { # mod_method == 'der'

		gene_class <- paste0(
			switch(k1_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V"),
			switch(k2_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V"),
			switch(k3_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V")
		)
		
		sim <- derivative_solution(
			simulation_time, gene_class,
			k1_function, k2_function, k3_function, 
			k1_params, k2_params, k3_params)
			
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
		
		
	} else {
		
		scores <- list()
		scores$pchisq <- NA
		scores$aic <- NA
		
	}
	
	conf_int <- list(
		k1 = cbind(left=rep(NA, length(simulation_time)),right=rep(NA, length(simulation_time))),
		k2 = cbind(left=rep(NA, length(simulation_time)),right=rep(NA, length(simulation_time))),
		k3 = cbind(left=rep(NA, length(simulation_time)),right=rep(NA, length(simulation_time)))
	)
	
	return(list(sim = sim, conf_int = conf_int, scores = scores))
	
}

RNAdynamicsAppMakeConfInt <- function(data_selection, 
																			time_min, time_max, experiment,
																			k1_function, k2_function,
																			k3_function, k1_params,
																			k2_params, k3_params, 
																			mod_method) {
	
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
		simulation_time <- seq(time_min,time_max,length.out=1000)
		simulation_time <- sort(unique(c(simulation_time, experiment_tpts)))
	} else {
		experiment_tpts <- 0
		simulation_time <- seq(0,16,length.out=1000)
	}
	
	if( mod_method == 'int' ) {
		
		gene_class <- paste0(
			switch(k1_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V"),
			switch(k2_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V"),
			switch(k3_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V")
		)
		
		conf_int <- compute_ci_Integrative_Nascent(c(k1_params, k2_params, k3_params),
																							 tpts = experiment_tpts,
																							 model_tpts = simulation_time,
																							 classTmp = gene_class,
																							 experimentalP = reference_preMRNA,
																							 experimentalM = reference_mRNA,
																							 experimentalA = reference_synthesis,
																							 varianceP = experiment$preMRNAsd^2,
																							 varianceM = experiment$mRNAsd^2,
																							 varianceA = experiment$synthesissd^2,
																							 confidenceThreshold = qchisq(.95,1)
		)
		
	} else { # mod_method == 'der'
		
		gene_class <- paste0(
			switch(k1_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V"),
			switch(k2_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V"),
			switch(k3_function, "Constant"="K", "Sigmoidal"="V", "Impulsive"="V")
		)
		
		conf_int <- compute_ci_Derivative_Nascent(c(k1_params, k2_params, k3_params),
																							tpts = experiment_tpts,
																							model_tpts = simulation_time,
																							classTmp = reconvert_gene_classes(gene_class),
																							experimentalP = reference_preMRNA,
																							experimentalM = reference_mRNA,
																							experimentalA = reference_synthesis,
																							varianceP = experiment$preMRNAsd^2,
																							varianceM = experiment$mRNAsd^2,
																							varianceA = if(is.null(experiment$synthesissd)) NULL else experiment$synthesissd^2,
																							confidenceThreshold = qchisq(.95,1)
		)
		
		
	}
	
	# calculate the scores of the modeling and assign to output
	
	p_k1 <- rate_var_p(conf_int$k1[simulation_time %in% experiment_tpts,])
	p_k2 <- rate_var_p(conf_int$k2[simulation_time %in% experiment_tpts,])
	p_k3 <- rate_var_p(conf_int$k3[simulation_time %in% experiment_tpts,])
	rate_p <- c(k1=p_k1, k2=p_k2, k3=p_k3)
	
	return(list(conf_int = conf_int, rate_p = rate_p))
	
}

RNAdynamicsAppPlot <- function(data_selection, 
															 show_logtime, show_relexpr,
															 logshift, linshift, 
															 time_min, time_max, 
															 experiment, 
															 simdata,
															 ylims,
															 rate_p
															 ) {

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
		simulation_time <- seq(time_min,time_max,length.out=1000)
		simulation_time <- sort(unique(c(simulation_time, experiment_tpts)))
	} else {
		experiment_tpts <- seq(0,16,by=4)
		simulation_time <- seq(0,16,length.out=1000)
		reference_mRNA <- c(reference_mRNA, rep(NA, 4))
		secondary_mRNA <- c(secondary_mRNA, rep(NA, 4))
		reference_preMRNA <- c(reference_preMRNA, rep(NA, 4))
		secondary_preMRNA <- c(secondary_preMRNA, rep(NA, 4))
		reference_synthesis <- c(reference_synthesis, rep(NA, 4))
		secondary_synthesis <- c(secondary_synthesis, rep(NA, 4))
		experimental_mRNAsd <- c(experimental_mRNAsd, rep(NA, 4))
		experimental_preMRNAsd <- c(experimental_preMRNAsd, rep(NA, 4))
		experimental_synthesissd <- c(experimental_synthesissd, rep(NA, 4))
	}
	
	# make the simulation

	if( !show_logtime ) {
		simtimeplot <- simulation_time 
		exptimeplot <- experiment_tpts
	} else {
		simtimeplot <- timetransf(simulation_time, logshift)
		exptimeplot <- timetransf(experiment_tpts, logshift)
	}
	
	sim <- simdata$sim
	conf_int <- simdata$conf_int

	# start plot routine
	
	par(mfrow=c(5,1))
	par(mar=c(2.5,8,0,1)+.1)
	
	# plot k1

	plot_k1_experiment = ! (data_selection == 'User defined' | experiment$no_nascent)
	k1_ylim <- plotSingleRNADynamic( 'synthesis', simtimeplot, sim[,'k1'], conf_int$k1[,'left'], conf_int$k1[,'right'], 
												plot_k1_experiment, exptimeplot, reference_synthesis, secondary_synthesis, experimental_synthesissd, show_relexpr, ylims$k1_ylim, rate_p = rate_p['k1'] )
	
	# plot pre-RNA dynamics

	p_ylim <- plotSingleRNADynamic( 'pre-RNA', simtimeplot, sim[,'p'], rep(NA, length(simtimeplot)), rep(NA, length(simtimeplot)), 
												data_selection != 'User defined', exptimeplot, reference_preMRNA, secondary_preMRNA, experimental_preMRNAsd, show_relexpr, ylims$p_ylim )
	
	# plot k2

	k2_ylim <- plotSingleRNADynamic( 'processing', simtimeplot, sim[,'k2'], conf_int$k2[,'left'], conf_int$k2[,'right'], 
												FALSE, show_relexpr = show_relexpr, ylim = ylims$k2_ylim, rate_p = rate_p['k2'] )#, exptimeplot, reference_synthesis, secondary_synthesis, experimental_synthesissd )
	
	# plot mRNA dynamics 

	m_ylim <- plotSingleRNADynamic( 'mature RNA', simtimeplot, sim[,'m'], rep(NA, length(simtimeplot)), rep(NA, length(simtimeplot)), 
												data_selection != 'User defined', exptimeplot, reference_mRNA, secondary_mRNA, experimental_mRNAsd, show_relexpr, ylims$m_ylim )
	
	# plot k3

	k3_ylim <- plotSingleRNADynamic( 'degradation', simtimeplot, sim[,'k3'], conf_int$k3[,'left'], conf_int$k3[,'right'], 
												FALSE, show_relexpr = show_relexpr, ylim = ylims$k3_ylim, rate_p = rate_p['k3'] )#, exptimeplot, reference_synthesis, secondary_synthesis, experimental_synthesissd )
	
	# draw x-axis
	if( show_logtime ) {
		axis(1, at=exptimeplot, labels = signif(experiment_tpts,2) , cex.axis = 1.3)	
	} else {
		axis(1, at=experiment_tpts, labels = signif(experiment_tpts,2) , cex.axis = 1.3)	
	}
	
	# return ylims upon request
	ylims <- list(
		k1_ylim = k1_ylim, 
		k2_ylim = k2_ylim, 
		k3_ylim = k3_ylim, 
		p_ylim = p_ylim, 
		m_ylim = m_ylim
	)
}

plotSingleRNADynamic <- function( dyn_name, simtimeplot, simprofile, ci_left, ci_right, plot_exp, exptimeplot, ref_exp, sec_exp, ssd_exp, show_relexpr = FALSE, ylim, rate_p = NULL ) {
	
	if( !is.null(rate_p) ) {
		dyn_name <- paste(dyn_name, paste0('(p=',signif(rate_p,2),')'), sep = '\n')
	} else {
		dyn_name <- paste(dyn_name, '', sep = '\n')
	}
	deltaylim <- function( yrange ) {
		deltarange <- yrange[2] * .05
		ylim <- yrange + c(-deltarange, deltarange)
	}

	if( plot_exp ) {
		sec_exp_plus_ssd <- sec_exp + ssd_exp
		ref_exp_plus_ssd <- ref_exp + ssd_exp
		sec_exp_minus_ssd <- sec_exp - ssd_exp
		ref_exp_minus_ssd <- ref_exp - ssd_exp
	}
	
	if(show_relexpr) {
		refexpression <- simprofile[1]
		simprofile <- simprofile/refexpression
		ci_left <- ci_left/refexpression
		ci_right <- ci_right/refexpression
		if( plot_exp ) {
			sec_exp <- sec_exp/refexpression
			ref_exp <- ref_exp/refexpression
			sec_exp_plus_ssd <- sec_exp_plus_ssd/refexpression
			ref_exp_plus_ssd <- ref_exp_plus_ssd/refexpression
			sec_exp_minus_ssd <- sec_exp_minus_ssd/refexpression
			ref_exp_minus_ssd <- ref_exp_minus_ssd/refexpression
		}
	}
	
	if( is.null(ylim) ) {
		if( plot_exp ) {
			yrange <- range(c(simprofile, 
												sec_exp_plus_ssd, 
												ref_exp_plus_ssd, 
												sec_exp_minus_ssd, 
												ref_exp_minus_ssd), na.rm=TRUE)
			ylim <- deltaylim(yrange)
		} else {
			ylim <- deltaylim( range(c(simprofile, ci_left, ci_right), na.rm=TRUE) )
		}
	}
	plot(simtimeplot, simprofile, 
			 xaxs='i', yaxs='i', xaxt = 'n',
			 ylab = dyn_name, type='l', xlab='', lwd=2, cex.lab = 1.7, cex.axis=1.3,  
			 xlim = range(simtimeplot) 
			 + diff(range(simtimeplot)) * c(-.05, .05),
			 ylim = ylim
	)
	matlines(simtimeplot, cbind(ci_left, ci_right), lty=2, col=1)	
	if( plot_exp ) {
		points( exptimeplot, sec_exp, pch=1, col='grey')
		points( exptimeplot, ref_exp, pch=19)
		segments( exptimeplot , ref_exp_minus_ssd 
							, exptimeplot , ref_exp_plus_ssd )
	}
	# return ylim upon request
	ylim <- ylim
}


rate_var_p <- function(rate_conf_int) {
	k_start <- mean(rate_conf_int[,2],na.rm=TRUE)
	if(!is.finite(k_start)) return(NaN) #return(list(par=NaN, value=NaN))
	k_scores_out <- optim(k_start, k_score_fun, method='BFGS', rate_conf_int=rate_conf_int)
	pchisq(k_scores_out$value,nrow(rate_conf_int)-1,lower.tail=FALSE)
}

constantModelRNApp <- function(x , par, log_shift, lin_shift) rep(par, length(x))

sigmoidModelRNApp <- function(x, par, log_shift, lin_shift=0)
	par[1]+(par[2]-par[1])*(1/(1+exp(-par[4]*(timetransf(x,log_shift,lin_shift)-timetransf(par[3],log_shift,lin_shift)))))

impulseModelRNApp <- function(x, par, log_shift, lin_shift=0)
	1/par[2]*(par[1]+(par[2]-par[1])*(1/(1+exp(-par[6]*(timetransf(x,log_shift,lin_shift)-timetransf(par[4],log_shift,lin_shift))))))*
	(par[3]+(par[2]-par[3])*(1/(1+exp(par[6]*(timetransf(x,log_shift,lin_shift)-timetransf(par[5],log_shift,lin_shift))))))


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

deterministic_simulation <- function(simulation_time,
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
						 "Constant" = constantModel(x, k1_params),
						 "Sigmoidal" = sigmoidModel(x, k1_params),
						 "Impulsive" = impulseModel(x, k1_params)
			)
		,
		beta = function(x) 
			switch(k3_function, 
						 "Constant" = constantModel(x, k3_params),
						 "Sigmoidal" = sigmoidModel(x, k3_params),
						 "Impulsive" = impulseModel(x, k3_params)
			)
		,
		gamma = function(x) 
			switch(k2_function, 
						 "Constant" = constantModel(x, k2_params),
						 "Sigmoidal" = sigmoidModel(x, k2_params),
						 "Impulsive" = impulseModel(x, k2_params)
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
	model <- data.frame(
		time=model[,1]
		, k1=synthesis
		, k2=processing
		, k3=degradation
		, p=model[,2]
		, m=model[,3]
		)
	return(model)
}

derivative_solution <- function(simulation_time, gene_class,
																k1_function, k2_function, k3_function,
																k1_params, k2_params, k3_params
)
	# in this version of the function, 
	# all three rates are modeled as impulse models
{
	
	tpts <- simulation_time
	model <- data.frame(
		time = tpts,
		k1 =  switch(gene_class,
								 "KKK"=k1KKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKK"=k1VKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVK"=k1KVK_Der_App(tpts, c(k1_params, k2_params, k3_params)),
								 "KKV"=k1KKV_Der_App(tpts, c(k1_params, k2_params, k3_params)),
								 "VVK"=k1VVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKV"=k1VKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVV"=k1KVV_Der_App(tpts, c(k1_params, k2_params, k3_params)),
								 "VVV"=k1VVV_Der(tpts, c(k1_params, k2_params, k3_params))
			)
		,
		k2 = switch(k2_function, 
						 "Constant" = constantModel(tpts, k2_params),
						 "Sigmoidal" = sigmoidModel(tpts, k2_params),
						 "Impulsive" = impulseModel(tpts, k2_params)
			)
		,
		k3 = switch(k3_function, 
						 "Constant" = constantModel(tpts, k3_params),
						 "Sigmoidal" = sigmoidModel(tpts, k3_params),
						 "Impulsive" = impulseModel(tpts, k3_params)
			)
		,
		p = switch(gene_class,
							 "KKK"=prematureKKK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VKK"=prematureVKK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KVK"=prematureKVK_Der_App(tpts, c(k1_params, k2_params, k3_params)),
							 "KKV"=prematureKKV_Der_App(tpts, c(k1_params, k2_params, k3_params)),
							 "VVK"=prematureVVK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VKV"=prematureVKV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KVV"=prematureKVV_Der_App(tpts, c(k1_params, k2_params, k3_params)),
							 "VVV"=prematureVVV_Der(tpts, c(k1_params, k2_params, k3_params))
		)
		,
		m = switch(k1_function, 
							 "Constant" = constantModel(tpts, k1_params),
							 "Sigmoidal" = sigmoidModel(tpts, k1_params),
							 "Impulsive" = impulseModel(tpts, k1_params)
		)
	)
	return(model)
}

derivative_solution_no_nascent <- function(simulation_time, gene_class,
																k1_function, k2_function, k3_function,
																k1_params, k2_params, k3_params
)
	# in this version of the function, 
	# all three rates are modeled as impulse models
{
	
	tpts <- simulation_time
	model <- data.frame(
		time = tpts,
		k1 =  switch(gene_class,
								 "KKK"=k1KKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKK"=k1VKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVK"=k1KVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KKV"=k1KKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VVK"=k1VVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKV"=k1VKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVV"=k1KVV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VVV"=k1VVV_Der(tpts, c(k1_params, k2_params, k3_params))
		)
		,
		k2 =  switch(gene_class,
								 "KKK"=k2KKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKK"=k2VKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVK"=k2KVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KKV"=k2KKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VVK"=k2VVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKV"=k2VKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVV"=k2KVV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VVV"=k2VVV_Der(tpts, c(k1_params, k2_params, k3_params))
		)
		,
		k3 =  switch(gene_class,
								 "KKK"=k3KKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKK"=k3VKK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVK"=k3KVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KKV"=k3KKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VVK"=k3VVK_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VKV"=k3VKV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "KVV"=k3KVV_Der(tpts, c(k1_params, k2_params, k3_params)),
								 "VVV"=k3VVV_Der(tpts, c(k1_params, k2_params, k3_params))
		)
		,
		p = switch(gene_class,
							 "KKK"=prematureKKK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VKK"=prematureVKK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KVK"=prematureKVK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KKV"=prematureKKV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VVK"=prematureVVK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VKV"=prematureVKV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KVV"=prematureKVV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VVV"=prematureVVV_Der(tpts, c(k1_params, k2_params, k3_params))
		)
		,
		m = switch(gene_class,
							 "KKK"=matureKKK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VKK"=matureVKK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KVK"=matureKVK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KKV"=matureKKV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VVK"=matureVVK_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VKV"=matureVKV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "KVV"=matureKVV_Der(tpts, c(k1_params, k2_params, k3_params)),
							 "VVV"=matureVVV_Der(tpts, c(k1_params, k2_params, k3_params))
		)
	)
	return(model)
}

smoothModel <- function(tpts, experiment, nInit=10, nIter=500, seed=1234)
{
	
	if( !is.null(seed) ) set.seed(seed)

	optimFailOut <- function(e) list(par=NA, value=NA, counts=NA, convergence=1, message=e)

	im.parguess <- function(tpts , values, log_shift ) {
		tpts <- timetransf(tpts, log_shift)
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
	
	log_shift <- findttpar(tpts)
	
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

chisqFunction <- function(experiment, model, variance=NULL)
{
	if( is.null(variance)) variance <- stats::var(experiment)
	sum((experiment - model )^2/variance )
}

logLikelihoodFunction <- function(experiment, model, variance=NULL)
{
	if( is.null(variance)) variance <- stats::var(experiment)
	sum(log(2*pnorm(-abs(experiment-model),mean=0,sd=sqrt(variance))))
}

###############################
## MINIMIZATION FUNCTION ######
###############################

modelChisqMatureRNA <- function(par, tpts, fun, df, alpha_exp, alpha_var #, pval
	, mature_exp, mature_var, preMRNA_exp, preMRNA_var #, log_shift, lin_shift
	, no_nascent)
{
	splitpar <- split(par 
		, c(rep('alpha',df['alpha']), rep('beta',df['beta']), rep('gamma',df['gamma'])) 
		)
	#

	params <- list()
	params$alpha <- function(x) 
		fun$alpha$value(x, splitpar$alpha)#, log_shift, lin_shift)
	params$beta  <- function(x)
		fun$beta$value(x, splitpar$beta)#, log_shift, lin_shift)
	params$gamma <- function(x)
		fun$gamma$value(x, splitpar$gamma)#, log_shift, lin_shift)
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
	, mature_var, preMRNA_exp, preMRNA_var, maxit=500, method='NM' #, log_shift, lin_shift
	, no_nascent, mod_method)
{
	if( method == 'NM' ) method <- 'Nelder-Mead'
	
	if( mod_method == 'int') {
		tryCatch({
			optOut <- optim(
				par=c(interpRates$alpha$params, 
							interpRates$beta$params, 
							interpRates$gamma$params)
				, fn=modelChisqMatureRNA , tpts=tpts_exp
				, fun=sapply(interpRates, '[[', 'fun')
				, df=unlist(sapply(interpRates, '[[', 'df'))
				, alpha_exp=alpha_exp, mature_exp=mature_exp
				, preMRNA_exp=preMRNA_exp
				, alpha_var=alpha_var, mature_var=mature_var
				, preMRNA_var=preMRNA_var
				# , log_shift = log_shift
				# , lin_shift = lin_shift
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
	} else { ## derivative modeling
		tryCatch({
			k1_function <- interpRates[['alpha']]$type
			k2_function <- interpRates[['gamma']]$type
			k3_function <- interpRates[['beta']]$type
			gene_class <- paste0(
				switch(k1_function, "constant"="K", "sigmoid"="V", "impulse"="V"),
				switch(k2_function, "constant"="K", "sigmoid"="V", "impulse"="V"),
				switch(k3_function, "constant"="K", "sigmoid"="V", "impulse"="V")
			)
			if(gene_class=='KKK') { ## KKK error function requires less arguments
				optOut <- optim(unname(unlist(lapply(interpRates, '[[', 'params')))
												, errorKKK_Der
												, tpts = tpts_exp
												, premature = preMRNA_exp
												, mature = mature_exp
												, alpha = alpha_exp
												, prematureVariance = preMRNA_var
												, matureVariance = mature_var
												, alphaVariance = alpha_var
												, control = list(maxit = maxit * 1000)
												, method = method
												)				
			} else {

				error_fun_Der <- switch(gene_class,
																"VKK"=errorVKK_Der,
																"KVK"=errorKVK_Der_App,
																"KKV"=errorKKV_Der_App,
																"VVK"=errorVVK_Der,
																"VKV"=errorVKV_Der,
																"KVV"=errorKVV_Der_App,
																"VVV"=errorVVV_Der
				)
				
				if( no_nascent ) {
					## in case of no nascent optimization the native error function (not the _App) should 
					## be used also for K-- models and KKK, initialChisquare and initialDistances should
					## be set. In order to skip them, all their values are set to 1 and initialPenalityRelevance to 0
					
					optOut <- optim(unname(unlist(lapply(interpRates, '[[', 'params')))
													, error_fun_Der
													, tpts = tpts_exp
													, premature = preMRNA_exp
													, mature = mature_exp
													, alpha = alpha_exp
													, prematureVariance = preMRNA_var
													, matureVariance = mature_var
													, alphaVariance = alpha_var
													, KKK = c(1,1,1)
													, initialChisquare = 1
													, initialDistances = 1
													, initialPenalityRelevance = 0
													, derivativePenalityRelevance = 10^-50
													, clean = FALSE
													, control = list(maxit = maxit * 1000))
				} else { # with nascent

					optOut <- optim(unname(unlist(lapply(interpRates, '[[', 'params')))
													, error_fun_Der
													, tpts = tpts_exp
													, premature = preMRNA_exp
													, mature = mature_exp
													, alpha = alpha_exp
													, prematureVariance = preMRNA_var
													, matureVariance = mature_var
													, alphaVariance = alpha_var
													, KKK = NULL
													, initialChisquare = NULL
													, initialDistances = NULL
													, initialPenalityRelevance = 1
													, derivativePenalityRelevance = 10^-50
													, clean = FALSE
													, control = list(maxit = maxit * 1000))
				}
			}
			splitpar <- split(optOut$par
												, c(rep('alpha',interpRates$alpha$df)
														, rep('gamma',interpRates$gamma$df)
														, rep('beta',interpRates$beta$df)
												)
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
}

#####################################
########## confidence intervals ##############
########################################

compute_ci_Integrative_Nascent <- function(parameters,
																					tpts,
																					model_tpts = tpts,
																					classTmp,
																					experimentalP,
																					experimentalM,
																					experimentalA,
																					varianceP,
																					varianceM,
																					varianceA,
																					confidenceThreshold
)
{
	
	if( is.null(names(parameters)) ) {
		names(parameters) <- as.character(seq_along(parameters))
	}
	
	foe <- capture.output({ # Just to capture the output of multiroot function
		suppressWarnings({
			intervals <- sapply(names(parameters),function(parname)
			{
				par <- parameters[parname]
				
				mOut = list(
					left_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e-2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
					left_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1/2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
					center = tryCatch(multiroot(f = logLikelihoodCIerror, start = par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
					right_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1.5*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = FALSE),error=function(e)return(emptyList)),
					right_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = FALSE),error=function(e)return(emptyList))
				)
				precis = sapply(mOut, '[[', 'f.root')
				
				if( length(which(precis<1e-2))>0 )  {
					conf_int = sapply(mOut[which(precis<1e-2)], '[[', 'root')
					low_int = min(conf_int)
					high_int = max(conf_int)
					
					left = ifelse( low_int < par, low_int, NA)
					right = ifelse( high_int > par, high_int, NA)
					
					left = unname(left)
					right = unname(right)
					
				} else {
					left = NA
					right = NA
				}
				return(c(left,right))
			})
			intervals[1,!is.finite(intervals[2,])] <- NaN
			intervals[2,!is.finite(intervals[1,])] <- NaN
		})
	})
	
	perturbedRates <- matrix(rep(NaN,3*length(model_tpts)),ncol=1)
	for(parname in names(parameters))
	{
		for(extremePar in intervals[,parname])
		{
			perturbedParameters <- parameters
			perturbedParameters[parname] <- extremePar
			
			perturbedRates <- cbind(perturbedRates,rates_integrativeModels(tpts=model_tpts, class=classTmp, parameters=perturbedParameters))
		}
	};perturbedRates <- perturbedRates[,-1]
	perturbedRates[perturbedRates<0] <- 0
	
	optTmp <- rates_integrativeModels(tpts=model_tpts, class=classTmp, parameters=parameters)
	
	k1left <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,min,na.rm=TRUE)
	k1TC <- optTmp[grep("alpha",names(optTmp))]
	k1right <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,max,na.rm=TRUE)
	
	k2left <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,min,na.rm=TRUE)
	k2TC <- optTmp[grep("gamma",names(optTmp))]
	k2right <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,max,na.rm=TRUE)
	
	k3left <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,min,na.rm=TRUE)
	k3TC <- optTmp[grep("beta",names(optTmp))]
	k3right <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,max,na.rm=TRUE)
	
	return(list(
		k1 = cbind(left=k1left, opt=k1TC, right=k1right),
		k2 = cbind(left=k2left, opt=k2TC, right=k2right),
		k3 = cbind(left=k3left, opt=k3TC, right=k3right)
	))
	
}

compute_ci_Derivative_Nascent <- function(parameters,
																					tpts,
																					model_tpts = tpts,
																					classTmp,
																					experimentalP,
																					experimentalM,
																					experimentalA,
																					varianceP,
																					varianceM,
																					varianceA,
																					confidenceThreshold
																					)
{
	
	if( is.null(names(parameters)) ) {
		names(parameters) <- as.character(seq_along(parameters))
	}
	
	foe <- capture.output({ # Just to capture the output of multiroot function
		suppressWarnings({
			intervals <- sapply(names(parameters),function(parname)
			{
				par <- parameters[parname]
				mOut = list(
					left_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e-2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = TRUE, app=TRUE),error=function(e)return(emptyList)),
					left_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1/2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = TRUE, app=TRUE),error=function(e)return(emptyList)),
					center = tryCatch(multiroot(f = logLikelihoodCIerror, start = par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = TRUE, app=TRUE),error=function(e)return(emptyList)),
					right_1 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1.5*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = TRUE, app=TRUE),error=function(e)return(emptyList)),
					right_2 = tryCatch(multiroot(f = logLikelihoodCIerror, start = 1e2*par, name = parname, parameters = parameters, class = classTmp, tpts = tpts, experimentalP = experimentalP, experimentalM = experimentalM, experimentalA = experimentalA, varianceP = varianceP, varianceM = varianceM, varianceA = varianceA, confidenceThreshold = confidenceThreshold, derivative = TRUE, app=TRUE),error=function(e)return(emptyList))
				)
				precis = sapply(mOut, '[[', 'f.root')
				if( length(which(precis<1e-2))>0 )  {
					conf_int = sapply(mOut[which(precis<1e-2)], '[[', 'root')
					low_int = min(conf_int)
					high_int = max(conf_int)
					
					left = ifelse( low_int < par, low_int, NA)
					right = ifelse( high_int > par, high_int, NA)
					
					left = unname(left)
					right = unname(right)
					
				} else {
					left = NA
					right = NA
				}
				return(c(left,right))
			})
			intervals[1,!is.finite(intervals[2,])] <- NaN
			intervals[2,!is.finite(intervals[1,])] <- NaN
		})
	})
	
	optTmp <- rates_derivativeModels(tpts=model_tpts, class=classTmp, parameters=parameters, app=TRUE)
	
	perturbedRates <- matrix(rep(NaN,3*length(model_tpts)),ncol=1)
	for(parname in names(parameters))
	{
		for(extremePar in intervals[,parname])
		{
			perturbedParameters <- parameters
			perturbedParameters[parname] <- extremePar
			
			perturbedRates <- cbind(perturbedRates,rates_derivativeModels(tpts=model_tpts, class=classTmp, parameters=perturbedParameters, app=TRUE))
		}
	};perturbedRates <- perturbedRates[,-1]
	perturbedRates[perturbedRates<0] <- 0
	
	k1left <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,min,na.rm=TRUE)
	k1TC <- optTmp[grep("alpha",names(optTmp))]
	k1right <- apply(perturbedRates[grep("alpha",rownames(perturbedRates)),],1,max,na.rm=TRUE)
	
	k2left <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,min,na.rm=TRUE)
	k2TC <- optTmp[grep("gamma",names(optTmp))]
	k2right <- apply(perturbedRates[grep("gamma",rownames(perturbedRates)),],1,max,na.rm=TRUE)
	
	k3left <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,min,na.rm=TRUE)
	k3TC <- optTmp[grep("beta",names(optTmp))]
	k3right <- apply(perturbedRates[grep("beta",rownames(perturbedRates)),],1,max,na.rm=TRUE)

	return(list(
		k1 = cbind(left=k1left, opt=k1TC, right=k1right),
		k2 = cbind(left=k2left, opt=k2TC, right=k2right),
		k3 = cbind(left=k3left, opt=k3TC, right=k3right)
	))	
	
}

## functions for solve the derivative system specific for the app (always centered on mature RNA)

##### KVK

k1KVK_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8]
		return( ( k3Parameters * matureParameters ) * ( 1 + .DimpulseModel(x, k2Parameters) / impulseModel(x, k2Parameters)^2 ) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6]
		return( ( k3Parameters * matureParameters ) * ( 1 + .DsigmoidModel(x, k2Parameters) / sigmoidModel(x, k2Parameters)^2 ) )
	}
}

k2KVK_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8]
		return( impulseModel(x, k2Parameters) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6]
		return( sigmoidModel(x, k2Parameters) )
	}
}

k3KVK_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8]
		return( rep(k3Parameters, length(tpts)) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6]
		return( rep(k3Parameters, length(tpts)) )
	}
}

prematureKVK_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8]
		return( ( k3Parameters * matureParameters ) / impulseModel(x, k2Parameters)	 )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6]
		return( ( k3Parameters * matureParameters ) / sigmoidModel(x, k2Parameters)	 )
	}
	
}

errorKVK_Der_App <- function(parameters, tpts
												 , premature, mature, alpha
												 , prematureVariance, matureVariance, alphaVariance
												 , KKK = NULL
												 , initialChisquare = NULL
												 , initialDistances = NULL
												 , initialPenalityRelevance = 1
												 , derivativePenalityRelevance = 10^-50
												 , clean)
{
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8]
		
		D0_M <- 0
		D0_k2 <- .DimpulseModel(0,k2Parameters)
		D0_k3 <- 0
		D0_P <- k3Parameters * matureParameters * .DimpulseModel(0, k2Parameters) / impulseModel(0, k2Parameters)^2
		
	} else {
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6]
		
		D0_M <- 0
		D0_k2 <- .DsigmoidModel(0,k2Parameters)
		D0_k3 <- 0
		D0_P <- k3Parameters * matureParameters * .DsigmoidModel(0, k2Parameters) / sigmoidModel(0, k2Parameters)^2
	}
	
	matureEstimated <- rep(matureParameters, length(tpts))
	prematureEstimated <- prematureKVK_Der_App(x = tpts, parameters = parameters)
	alphaEstimated <- k1KVK_Der_App(x = tpts, parameters = parameters)
	
	alphaEstimated[alphaEstimated<0] <- NaN
	prematureEstimated[prematureEstimated<0] <- NaN
	matureEstimated[matureEstimated<0] <- NaN
	
	if(any(!is.finite(alphaEstimated)) | 
		 any(!is.finite(prematureEstimated)) | 
		 any(!is.finite(matureEstimated)) | 
		 !is.finite(D0_M) | 
		 !is.finite(D0_k2) | 
		 !is.finite(D0_k3) | 
		 !is.finite(D0_P)
	) return(NaN)
	
	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)
	
	if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
	{
		alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
		initialPenality <- 0
	}else{
		# stop('errorKVK_Der_App: KKK version not implemented')
		if(clean){initialPenality <- 0}else{
			initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1KVK_Der(0,parameters))^2
																																											 + (k2KKK_Der(0,KKK)-k2KVK_Der(0,parameters))^2
																																											 + (k3KKK_Der(0,KKK)-k3KVK_Der(0,parameters))^2)
		}
		alphaChiSquare <- 0
	}
	
	chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
	penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)
	
	if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}
	
	if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
}


###### KKV

k1KKV_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:8]
		return( matureParameters * ( .DimpulseModel(x, k3Parameters) / k2Parameters + impulseModel(x, k3Parameters) ) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:6]
		return( matureParameters * ( .DsigmoidModel(x, k3Parameters) / k2Parameters + sigmoidModel(x, k3Parameters) ) )
	}
}

k2KKV_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:8]
		return( rep(k2Parameters, length(tpts)) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:6]
		return( rep(k2Parameters, length(tpts)) )
	}
}

k3KKV_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:8]
		return( impulseModel(x, k3Parameters) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:6]
		return( sigmoidModel(x, k3Parameters) )
	}
}

prematureKKV_Der_App <- function(x, parameters) {
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:8]
		return( ( impulseModel(x, k3Parameters) * matureParameters ) / k2Parameters	 )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:6]
		return( ( sigmoidModel(x, k3Parameters) * matureParameters ) / k2Parameters	 )
	}
}

errorKKV_Der_App <- function(parameters, tpts
														 , premature, mature, alpha
														 , prematureVariance, matureVariance, alphaVariance
														 , KKK = NULL
														 , initialChisquare = NULL
														 , initialDistances = NULL
														 , initialPenalityRelevance = 1
														 , derivativePenalityRelevance = 10^-50
														 , clean)
{
	if(length(parameters)==8)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:8]
		
		D0_M <- 0
		D0_k2 <- 0
		D0_k3 <- .DimpulseModel(0,k3Parameters)
		D0_P <- .DimpulseModel(0, k3Parameters) * matureParameters / k2Parameters
		
	} else {
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2]
		k3Parameters <- parameters[3:6]
		
		D0_M <- 0
		D0_k2 <- 0
		D0_k3 <- .DsigmoidModel(0,k3Parameters)
		D0_P <- .DsigmoidModel(0, k3Parameters) * matureParameters / k2Parameters
	}
	
	matureEstimated <- rep(matureParameters, length(tpts))
	prematureEstimated <- prematureKKV_Der_App(x = tpts, parameters = parameters)
	alphaEstimated <- k1KKV_Der_App(x = tpts, parameters = parameters)
	
	alphaEstimated[alphaEstimated<0] <- NaN
	prematureEstimated[prematureEstimated<0] <- NaN
	matureEstimated[matureEstimated<0] <- NaN
	
	if(any(!is.finite(alphaEstimated)) | 
		 any(!is.finite(prematureEstimated)) | 
		 any(!is.finite(matureEstimated)) | 
		 !is.finite(D0_M) | 
		 !is.finite(D0_k2) | 
		 !is.finite(D0_k3) | 
		 !is.finite(D0_P)
	) return(NaN)
	
	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)
	
	if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
	{
		alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
		initialPenality <- 0
	}else{
		# stop('errorKKV_Der_App: KKK version not implemented')
		if(clean){initialPenality <- 0}else{
			initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1KKV_Der(0,parameters))^2
																																											 + (k2KKK_Der(0,KKK)-k2KKV_Der(0,parameters))^2
																																											 + (k3KKK_Der(0,KKK)-k3KKV_Der(0,parameters))^2)
		}
		alphaChiSquare <- 0
	}
	
	chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
	penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)
	
	if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}
	
	if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
}

########## KVV

k1KVV_Der_App <- function(x, parameters) {
	if(length(parameters)==13)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8:13]
		return( matureParameters * ( impulseModel(x, k3Parameters) + ( impulseModel(x, k2Parameters) * .DimpulseModel(x, k3Parameters) - 
						impulseModel(x, k3Parameters) * .DimpulseModel(x, k2Parameters) ) / impulseModel(x, k2Parameters)^2  ) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6:9]
		return( matureParameters * ( sigmoidModel(x, k3Parameters) + ( sigmoidModel(x, k2Parameters) * .DsigmoidModel(x, k3Parameters) - 
						sigmoidModel(x, k3Parameters) * .DsigmoidModel(x, k2Parameters) ) / sigmoidModel(x, k2Parameters)^2  ) )
	}
}

k2KVV_Der_App <- function(x, parameters) {
	if(length(parameters)==13)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8:13]
		return( impulseModel(x, k2Parameters) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6:9]
		return( sigmoidModel(x, k2Parameters) )
	}
}

k3KVV_Der_App <- function(x, parameters) {
	if(length(parameters)==13)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8:13]
		return( impulseModel(x, k3Parameters) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6:9]
		return( sigmoidModel(x, k3Parameters) )
	}
}

prematureKVV_Der_App <- function(x, parameters) {
	if(length(parameters)==13)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8:13]
		return( ( impulseModel(x, k3Parameters) * matureParameters ) / impulseModel(x, k2Parameters) )
	}else{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6:9]
		return( ( sigmoidModel(x, k3Parameters) * matureParameters ) / sigmoidModel(x, k2Parameters) )
	}
}

errorKVV_Der_App <- function(parameters, tpts
														 , premature, mature, alpha
														 , prematureVariance, matureVariance, alphaVariance
														 , KKK = NULL
														 , initialChisquare = NULL
														 , initialDistances = NULL
														 , initialPenalityRelevance = 1
														 , derivativePenalityRelevance = 10^-50
														 , clean)
{
	if(length(parameters)==13)
	{
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:7]
		k3Parameters <- parameters[8:13]
		
		D0_M <- 0
		D0_k2 <- .DimpulseModel(0,k2Parameters)
		D0_k3 <- .DimpulseModel(0,k3Parameters)
		D0_P <- matureParameters * ( ( impulseModel(0, k2Parameters) * .DimpulseModel(0, k3Parameters) - 
																	 	impulseModel(0, k3Parameters) * .DimpulseModel(0, k2Parameters) ) / impulseModel(0, k2Parameters)^2  )
		
	} else {
		matureParameters <- parameters[1]
		k2Parameters <- parameters[2:5]
		k3Parameters <- parameters[6:9]
		
		D0_M <- 0
		D0_k2 <- 0
		D0_k3 <- .DsigmoidModel(0,k3Parameters)
		D0_P <- matureParameters * ( ( sigmoidModel(0, k2Parameters) * .DsigmoidModel(0, k3Parameters) - 
																	 	sigmoidModel(0, k3Parameters) * .DsigmoidModel(0, k2Parameters) ) / sigmoidModel(0, k2Parameters)^2  )
	}
	
	matureEstimated <- rep(matureParameters, length(tpts))
	prematureEstimated <- prematureKVV_Der_App(x = tpts, parameters = parameters)
	alphaEstimated <- k1KVV_Der_App(x = tpts, parameters = parameters)
	
	alphaEstimated[alphaEstimated<0] <- NaN
	prematureEstimated[prematureEstimated<0] <- NaN
	matureEstimated[matureEstimated<0] <- NaN
	
	if(any(!is.finite(alphaEstimated)) | 
		 any(!is.finite(prematureEstimated)) | 
		 any(!is.finite(matureEstimated)) | 
		 !is.finite(D0_M) | 
		 !is.finite(D0_k2) | 
		 !is.finite(D0_k3) | 
		 !is.finite(D0_P)
	) return(NaN)
	
	prematureChiSquare <- sum((premature - prematureEstimated )^2/prematureVariance)
	matureChiSquare <- sum((mature - matureEstimated)^2/matureVariance)
	
	if(is.null(KKK)&is.null(initialChisquare)&is.null(initialDistances)&!is.null(alpha)&!is.null(alphaVariance))
	{
		alphaChiSquare <- sum((alpha - alphaEstimated)^2/alphaVariance)
		initialPenality <- 0
	}else{
		# stop('errorKVV_Der_App: KKK version not implemented')
		if(clean){initialPenality <- 0}else{
			initialPenality <- initialPenalityRelevance*(initialChisquare/initialDistances)*((k1KKK_Der(0,KKK)-k1KVV_Der(0,parameters))^2
																																											 + (k2KKK_Der(0,KKK)-k2KVV_Der(0,parameters))^2
																																											 + (k3KKK_Der(0,KKK)-k3KVV_Der(0,parameters))^2)
		}
		alphaChiSquare <- 0
	}
	
	chiSquare <- sum(c(prematureChiSquare,matureChiSquare,alphaChiSquare))
	penalty <- abs(D0_M)+abs(D0_P)+abs(D0_k2)+abs(D0_k3)
	
	if(penalty <= chiSquare*derivativePenalityRelevance){penalty <- 0}
	
	if(clean){return(chiSquare)}else{return(chiSquare+penalty+initialPenality)}
}