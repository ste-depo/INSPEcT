shinyServer(function(input, output, session) {
	
	# global variables
	
	ranges <- reactiveValues()
	values <- reactiveValues()
	experiment <- reactiveValues()
	inspect <- reactiveValues(loaded=FALSE)
	function_types <- reactiveValues()
	modeling <- reactiveValues()
	
	######################################################################
	######################################################################
	### import of the INSPEcT dataset and update of the imported values
	######################################################################
	######################################################################
	
	# load INSPEcT file
	
	contentsrea <- reactive({

		filename <- input$file1$datapath 
		if( is.null(filename) ) filename <- system.file(package='INSPEcT', 'nascentInspObj10.rds')

		## load file
		if( file.exists(filename) ) {
			ids <- readRDS(filename)

			if( class(ids) != 'INSPEcT' ) {
				return(NULL)

			} else { # (the loaded object is of class INSPEcT
				
				## store inspect global values
				experiment$tpts <- tpts(ids)
				experiment$no_nascent <- ids@NoNascent
				experiment$steady_state <- is.character(experiment$tpts)

				# select only genes with exons and introns
				ids <- ids[apply(is.finite(ratesFirstGuess(ids,'preMRNA')),1,all)]

				if( !experiment$steady_state ) {

					inspect$mod_method <- modelingParams(ids)$estimateRatesWith ## either "der" or "int"
					inspect$classes <- geneClass(ids)
					inspect$logshift <- findttpar(experiment$tpts)
					inspect$linshift <- ifelse( experiment$no_nascent,
						abs(min(timetransf(experiment$tpts,inspect$logshift))),0)

					# ... optionally throw a warning and consider as steady-state
					if( ids@NF ) stop("The RNAdynamics app doesn't work with non-functional(NF) INSPEcT models")
					
					## select only the best model
					
					if( ids@NoNascent ){ # (based on the class)
						ids@model@ratesSpecs <-
							lapply(seq_along(inspect$classes), function(i)
								list(ids@model@ratesSpecs[[i]][[inspect$classes[i]]]))
						names(ids@model@ratesSpecs) <- featureNames(ids)
					} else{ # Nascent RNA object (always VVV)
						ids@model@ratesSpecs <-
							lapply(seq_along(inspect$classes), function(i)
								list(ids@model@ratesSpecs[[i]][['abc']]))
						names(ids@model@ratesSpecs) <- featureNames(ids)
					}
					
					## update (converted) gene classes in the select input box
					classes_table <- table(isolate(inspect$classes))
					names(classes_table) <- convert_gene_classes( names(classes_table) )
					classes_table_string <- paste( names(classes_table) , '(', classes_table, ')' )
					updateSelectInput(session, "select_class", 
						choices = classes_table_string, selected = classes_table_string[1])

					## define ranges

					ranges$time_min <- min(experiment$tpts)
					ranges$time_max <- max(experiment$tpts)

				} else { ## steady state

					inspect$mod_method <- 'int'
					inspect$logshift <- 1
					inspect$linshift <- 0

					updateSelectInput(session, "select_condition", 
						choices = experiment$tpts, 
						selected = experiment$tpts[1])

					updateSelectInput(session, "select", 
						choices = sort(featureNames(ids)), 
						selected = sort(featureNames(ids))[1])

					ranges$time_min <- 0
					ranges$time_max <- 16

				}

				## predifined  

				values$logtime <- FALSE
				values$confint <- FALSE
				
				return(ids)
				
			}
		} else { # (if the file name does not exist)
			return(NULL)
		}

	})
	
	################################################
	## update gene names in the select input box ######
	################################################
	
	observe({
		if( !is.null(input$select_class) ) {
			if( !experiment$steady_state ) {
				selected_class <- reconvert_gene_classes(strsplit( input$select_class , ' ')[[1]][1])
				updateSelectInput(session, "select", selected = NULL,
					choices = sort(featureNames(contentsrea()[inspect$classes == selected_class])))
			}
		}
	})
	
	############################################################################
	## select parameters of each rate (synthesis, processing and degradation) #####
	## and experimental values for the selected gene ###########################
	####################################################
	
	observe({

		ids <- contentsrea()
				
		if(
			input$select != "" & !is.null(ids) &
			input$select %in% featureNames(ids) & 
			!(is.null(input$select_condition) & experiment$steady_state ) 
			) {

			## define the modeling strategy
			
			gene_class <- strsplit( input$select_class , ' ')[[1]][1]
			if( inspect$mod_method == 'int' ) {
				modeling_type <- 'K123'
				model_names <- c('alpha','gamma','beta')
			} else { # inspect$mod_method == 'der'
				if( experiment$no_nascent & (gene_class  %in% c('KVK','KKV','KVV') ) ) {
					if( gene_class %in% c('KVK','KVV') ) {
						modeling_type <- 'TK13'
						model_names <- c('total','alpha','beta')
					} else {
						modeling_type <- 'TK12'
						model_names <- c('total','alpha','gamma')
					}
				} else {
					modeling_type <- 'MK23'
					model_names <- c('mature','gamma','beta')
				}
			}
			
			## assign to global variable
			values$modeling_type <- modeling_type
			values$model_names <- model_names
			
			if( !experiment$steady_state ) {

				experiment$synthesis <- if( !experiment$no_nascent ) {
					ratesFirstGuess(ids[input$select], 'synthesis')	
				} else NULL
				experiment$mRNA <- ratesFirstGuess(ids[input$select], 'total') - ratesFirstGuess(ids[input$select], 'preMRNA')
				experiment$preMRNA <- ratesFirstGuess(ids[input$select], 'preMRNA')
				
				## smooth experiment data
				experiment$synthesis_smooth <- if( !experiment$no_nascent ) {
					smoothModel(ids@tpts, ratesFirstGuess(ids[input$select], 'synthesis'))
				} else NULL
				experiment$mRNA_smooth <- smoothModel(ids@tpts, ratesFirstGuess(ids[input$select], 'total') - ratesFirstGuess(ids[input$select], 'preMRNA'))
				experiment$preMRNA_smooth <- smoothModel(ids@tpts, ratesFirstGuess(ids[input$select], 'preMRNA'))

				experiment$synthesissd <- if( !experiment$no_nascent ) {
					sqrt(ratesFirstGuessVar(ids[input$select], 'synthesis'))
				} else NULL
				experiment$mRNAsd <- sqrt(ratesFirstGuessVar(ids[input$select], 'total') + ratesFirstGuessVar(ids[input$select], 'preMRNA'))
				experiment$preMRNAsd <- sqrt(ratesFirstGuessVar(ids[input$select], 'preMRNA'))

				out <- define_parameter_ranges( ids )

				gene_model <- ids@model@ratesSpecs[[input$select]][[1]]
				modeling$counts <- gene_model$counts[1]
				modeling$convergence <- gene_model$convergence

			} else { ## steady state experiment
				
				condition_id <- which(experiment$tpts == input$select_condition)

				experiment$synthesis_smooth <- experiment$synthesis <- 
					ratesFirstGuess(ids[input$select,condition_id], 'synthesis')
				experiment$mRNA_smooth      <- experiment$mRNA <- 
					ratesFirstGuess(ids[input$select,condition_id], 'total') - 
						ratesFirstGuess(ids[input$select,condition_id], 'preMRNA')
				experiment$preMRNA_smooth   <- experiment$preMRNA <-   
					ratesFirstGuess(ids[input$select,condition_id], 'preMRNA')
				
				experiment$synthesissd <- sqrt(ratesFirstGuessVar(ids[input$select,condition_id], 'synthesis'))
				experiment$mRNAsd <- sqrt(ratesFirstGuessVar(ids[input$select,condition_id], 'total') + 
					ratesFirstGuessVar(ids[input$select,condition_id], 'preMRNA'))
				experiment$preMRNAsd <- sqrt(ratesFirstGuessVar(ids[input$select,condition_id], 'preMRNA'))
				
				rateRange <- function(rate) {
					rate_vals = ratesFirstGuess(ids, rate)
					rate_range = range(rate_vals[is.finite(rate_vals)])
					return(c(floor(rate_range[1]), ceiling(rate_range[2])))
				}

				out <- list(
					k1_h_pars=rateRange('synthesis'),
					k2_h_pars=rateRange('processing'),
					k3_h_pars=rateRange('degradation'),
					t_pars=c(0,16),
					beta_pars=c(0,10)
					)

				gene_model <- list(
					alpha=list(type='constant', fun=newPointer(constantModel),
						params=ratesFirstGuess(ids[input$select,condition_id], 'synthesis'), df=1),
					gamma=list(type='constant', fun=newPointer(constantModel),
						params=ratesFirstGuess(ids[input$select,condition_id], 'processing'), df=1),
					beta=list(type='constant', fun=newPointer(constantModel),
										params=ratesFirstGuess(ids[input$select,condition_id], 'degradation'), df=1)
				)
				
				model_names <- names(gene_model)
				
				modeling$counts <- NA
				modeling$convergence <- NA

			}
			
			## set to false the evaluation of CI on load of a new dataset
			updateCheckboxInput(session, "confint_checkbox", value = FALSE)

			try({ 
				## this try is necessary because the class is updated before the gene, then the model_name
				## corresponding to the new class can also not correspond to the currently selected gene
				
				function_types$k1 <- switch(
					gene_model[[ model_names[1] ]][['type']],
					"constant" = "Constant",
					"sigmoid" = "Sigmoidal",
					"impulse" = "Impulsive"
				)
				
				function_types$k2 <- switch(
					gene_model[[ model_names[2] ]][['type']],
					"constant" = "Constant",
					"sigmoid" = "Sigmoidal",
					"impulse" = "Impulsive"
				)
				
				function_types$k3 <- switch(
					gene_model[[ model_names[3] ]][['type']],
					"constant" = "Constant",
					"sigmoid" = "Sigmoidal",
					"impulse" = "Impulsive"
				)
				
				updateRadioButtons(session, 'k1_function',
													 selected = function_types$k1)
				updateRadioButtons(session, 'k2_function',
													 selected = function_types$k2)
				updateRadioButtons(session, 'k3_function',
													 selected = function_types$k3)
				
				if( function_types$k1 == "Constant" ) {
					values$k1_h0 = round(gene_model[[ model_names[1] ]][["params"]][1],2)
					values$k1_h1 = round(gene_model[[ model_names[1] ]][["params"]][1],2)
					values$k1_h2 = round(gene_model[[ model_names[1] ]][["params"]][1],2)
					values$k1_t1 = round(mean(out$t_pars))
					values$k1_t2 = round(mean(out$t_pars))
					values$k1_beta = round(mean(out$beta_pars))
				}
				if( function_types$k1 == "Sigmoidal" ) {
					values$k1_h0 = round(gene_model[[ model_names[1] ]][["params"]][1],2)
					values$k1_h1 = round(gene_model[[ model_names[1] ]][["params"]][2],2)
					values$k1_h2 = round(gene_model[[ model_names[1] ]][["params"]][2],2)
					values$k1_t1 = round(gene_model[[ model_names[1] ]][["params"]][3],2)
					values$k1_t2 = round(gene_model[[ model_names[1] ]][["params"]][3],2)
					values$k1_beta = round(gene_model[[ model_names[1] ]][["params"]][4],2)
				}
				if( function_types$k1 == "Impulsive" ) {
					values$k1_h0 = round(gene_model[[ model_names[1] ]][["params"]][1],2)
					values$k1_h1 = round(gene_model[[ model_names[1] ]][["params"]][2],2)
					values$k1_h2 = round(gene_model[[ model_names[1] ]][["params"]][3],2)
					values$k1_t1 = round(gene_model[[ model_names[1] ]][["params"]][4],2)
					values$k1_t2 = round(gene_model[[ model_names[1] ]][["params"]][5],2)
					values$k1_beta = round(gene_model[[ model_names[1] ]][["params"]][6],2)
				}
				
				if( function_types$k2 == "Constant" ) {
					values$k2_h0 = round(gene_model[[ model_names[2] ]][["params"]][1],2)
					values$k2_h1 = round(gene_model[[ model_names[2] ]][["params"]][1],2)
					values$k2_h2 = round(gene_model[[ model_names[2] ]][["params"]][1],2)
					values$k2_t1 = round(mean(out$t_pars))
					values$k2_t2 = round(mean(out$t_pars))
					values$k2_beta = round(mean(out$beta_pars))
				}
				if( function_types$k2 == "Sigmoidal" ) {
					values$k2_h0 = round(gene_model[[ model_names[2] ]][["params"]][1],2)
					values$k2_h1 = round(gene_model[[ model_names[2] ]][["params"]][2],2)
					values$k2_h2 = round(gene_model[[ model_names[2] ]][["params"]][2],2)
					values$k2_t1 = round(gene_model[[ model_names[2] ]][["params"]][3],2)
					values$k2_t2 = round(gene_model[[ model_names[2] ]][["params"]][3],2)
					values$k2_beta = round(gene_model[[ model_names[2] ]][["params"]][4],2)
				}
				if( function_types$k2 == "Impulsive" ) {
					values$k2_h0 = round(gene_model[[ model_names[2] ]][["params"]][1],2)
					values$k2_h1 = round(gene_model[[ model_names[2] ]][["params"]][2],2)
					values$k2_h2 = round(gene_model[[ model_names[2] ]][["params"]][3],2)
					values$k2_t1 = round(gene_model[[ model_names[2] ]][["params"]][4],2)
					values$k2_t2 = round(gene_model[[ model_names[2] ]][["params"]][5],2)
					values$k2_beta = round(gene_model[[ model_names[2] ]][["params"]][6],2)
				}
				
				if( function_types$k3 == "Constant" ) {
					values$k3_h0 = round(gene_model[[ model_names[3] ]][["params"]][1],2)
					values$k3_h1 = round(gene_model[[ model_names[3] ]][["params"]][1],2)
					values$k3_h2 = round(gene_model[[ model_names[3] ]][["params"]][1],2)
					values$k3_t1 = round(mean(out$t_pars))
					values$k3_t2 = round(mean(out$t_pars))
					values$k3_beta = round(mean(out$beta_pars))
				}
				if( function_types$k3 == "Sigmoidal" ) {
					values$k3_h0 = round(gene_model[[ model_names[3] ]][["params"]][1],2)
					values$k3_h1 = round(gene_model[[ model_names[3] ]][["params"]][2],2)
					values$k3_h2 = round(gene_model[[ model_names[3] ]][["params"]][2],2)
					values$k3_t1 = round(gene_model[[ model_names[3] ]][["params"]][3],2)
					values$k3_t2 = round(gene_model[[ model_names[3] ]][["params"]][3],2)
					values$k3_beta = round(gene_model[[ model_names[3] ]][["params"]][4],2)
				}
				if( function_types$k3 == "Impulsive" ) {
					values$k3_h0 = round(gene_model[[ model_names[3] ]][["params"]][1],2)
					values$k3_h1 = round(gene_model[[ model_names[3] ]][["params"]][2],2)
					values$k3_h2 = round(gene_model[[ model_names[3] ]][["params"]][3],2)
					values$k3_t1 = round(gene_model[[ model_names[3] ]][["params"]][4],2)
					values$k3_t2 = round(gene_model[[ model_names[3] ]][["params"]][5],2)
					values$k3_beta = round(gene_model[[ model_names[3] ]][["params"]][6],2)
				}
				
				## convert the models into MK23
				if( modeling_type == 'TK13' ) {
					if( gene_class == 'KVK') {
						parameters <- unname(unlist(lapply(gene_model[1:3], '[[', 'params')))
						processingfit <- k2KVK_Der(ids@tpts, parameters)
						processingmodel <- .chooseModel(ids@tpts, experiment = processingfit, variance = 1, 
																						sigmoid = TRUE, impulse = TRUE, polynomial = FALSE,
																						nInit = 10, nIter = 300, sigmoidModel=sigmoidModel, impulseModel=impulseModel,
																						sigmoidModelP=sigmoidModelP, impulseModelP=impulseModelP, .polynomialModelP=.polynomialModelP,
																						computeDerivatives = TRUE)
						maturefit <- matureKVK_Der(ids@tpts, parameters)
						maturemodel <- .chooseModel(ids@tpts, experiment = maturefit, variance = 1, 
																				sigmoid = processingmodel$type == 'sigmoid', impulse = processingmodel$type == 'impulse', polynomial = FALSE,
																				nInit = 10, nIter = 300, sigmoidModel=sigmoidModel, impulseModel=impulseModel,
																				sigmoidModelP=sigmoidModelP, impulseModelP=impulseModelP, .polynomialModelP=.polynomialModelP,
																				computeDerivatives = TRUE)
						k1_rate <- switch(processingmodel$type, 
															"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																								 params=maturemodel$params,df=4),
															"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																								 params=maturemodel$params,df=6)
						)
						k2_rate <- switch(processingmodel$type, 
															"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																								 params=processingmodel$params,df=4),
															"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																								 params=processingmodel$params,df=6)
						)
						k3_rate <- list(type='constant',fun=newPointer(constantModel),
														params=gene_model[[3]]$params,df=1)
						params <- list(alpha=k1_rate, gamma=k2_rate, beta=k3_rate)
						gene_model <- optimParamsMatureRNA(params, ids@tpts
																							 , if(input$data_selection == 'Smooth data') isolate(experiment$synthesis_smooth) else isolate(experiment$synthesis)
																							 , isolate(experiment$synthesissd)^2
																							 , if(input$data_selection == 'Smooth data') isolate(experiment$mRNA_smooth) else isolate(experiment$mRNA)
																							 , isolate(experiment$mRNAsd)^2
																							 , if(input$data_selection == 'Smooth data') isolate(experiment$preMRNA_smooth) else isolate(experiment$preMRNA)
																							 , isolate(experiment$preMRNAsd)^2
																							 , maxit=300, method=isolate(input$opt_method)
																							 , isolate(experiment$no_nascent), mod_method = isolate(inspect$mod_method))
						## update mature values
						if( processingmodel$type == 'sigmoid' ) {
							values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],2)
							values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],2)
							values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],2)
							values$k1_t1 = round(gene_model[[ 1 ]][["params"]][3],2)
							values$k1_t2 = round(gene_model[[ 1 ]][["params"]][3],2)
							values$k1_beta = round(gene_model[[ 1 ]][["params"]][4],2)
						} else {
							values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],2)
							values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],2)
							values$k1_h2 = round(gene_model[[ 1 ]][["params"]][3],2)
							values$k1_t1 = round(gene_model[[ 1 ]][["params"]][4],2)
							values$k1_t2 = round(gene_model[[ 1 ]][["params"]][5],2)
							values$k1_beta = round(gene_model[[ 1 ]][["params"]][6],2)
						}
						## update k2 values
						if( processingmodel$type == 'sigmoid' ) {
							values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],2)
							values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],2)
							values$k2_h2 = round(gene_model[[ 3 ]][["params"]][2],2)
							values$k2_t1 = round(gene_model[[ 3 ]][["params"]][3],2)
							values$k2_t2 = round(gene_model[[ 3 ]][["params"]][3],2)
							values$k2_beta = round(gene_model[[ 3 ]][["params"]][4],2)
						} else {
							values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],2)
							values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],2)
							values$k2_h2 = round(gene_model[[ 3 ]][["params"]][3],2)
							values$k2_t1 = round(gene_model[[ 3 ]][["params"]][4],2)
							values$k2_t2 = round(gene_model[[ 3 ]][["params"]][5],2)
							values$k2_beta = round(gene_model[[ 3 ]][["params"]][6],2)
						}
						values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],2)
						values$k3_h1 = round(gene_model[[ 2 ]][["params"]][1],2)
						values$k3_h2 = round(gene_model[[ 2 ]][["params"]][1],2)
						values$k3_t1 = round(mean(out$t_pars))
						values$k3_t2 = round(mean(out$t_pars))
						values$k3_beta = round(mean(out$beta_pars))
						## update functional forms
						updateRadioButtons(session, 'k1_function',
															 selected = if(maturemodel$type == 'impulse') "Impulsive" else "Sigmoidal")
						updateRadioButtons(session, 'k2_function',
															 selected = if(processingmodel$type == 'impulse') "Impulsive" else "Sigmoidal")
						updateRadioButtons(session, 'k3_function',
															 selected = "Constant")
					} else { # KVV
						modeling_fun <- gene_model[[3]]$type
						parameters <- unname(unlist(lapply(gene_model[1:3], '[[', 'params')))
						processingfit <- k2KVV_Der(ids@tpts, parameters)
						processingmodel <- .chooseModel(ids@tpts, experiment = processingfit, variance = 1, 
																						sigmoid = modeling_fun == 'sigmoid', impulse = modeling_fun == 'impulse', polynomial = FALSE,
																						nInit = 10, nIter = 300, sigmoidModel=sigmoidModel, impulseModel=impulseModel,
																						sigmoidModelP=sigmoidModelP, impulseModelP=impulseModelP, .polynomialModelP=.polynomialModelP,
																						computeDerivatives = TRUE)
						maturefit <- matureKVV_Der(ids@tpts, parameters)
						maturemodel <- .chooseModel(ids@tpts, experiment = maturefit, variance = 1, 
																				sigmoid = modeling_fun == 'sigmoid', impulse = modeling_fun == 'impulse', polynomial = FALSE,
																				nInit = 10, nIter = 300, sigmoidModel=sigmoidModel, impulseModel=impulseModel,
																				sigmoidModelP=sigmoidModelP, impulseModelP=impulseModelP, .polynomialModelP=.polynomialModelP,
																				computeDerivatives = TRUE)
						k1_rate <- switch(modeling_fun, 
															"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																							 params=maturemodel$params,df=4),
															"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																							 params=maturemodel$params,df=6)
						)
						k2_rate <- switch(modeling_fun, 
															"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																							 params=processingmodel$params,df=4),
															"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																							 params=processingmodel$params,df=6)
						)
						k3_rate <- switch(modeling_fun, 
															"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																							 params=gene_model[[3]]$params,df=4),
															"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																							 params=gene_model[[3]]$params,df=6)
						)
						params <- list(alpha=k1_rate, gamma=k2_rate, beta=k3_rate)
						gene_model <- optimParamsMatureRNA(params, ids@tpts
																							 , if(input$data_selection == 'Smooth data') isolate(experiment$synthesis_smooth) else isolate(experiment$synthesis)
																							 , isolate(experiment$synthesissd)^2
																							 , if(input$data_selection == 'Smooth data') isolate(experiment$mRNA_smooth) else isolate(experiment$mRNA)
																							 , isolate(experiment$mRNAsd)^2
																							 , if(input$data_selection == 'Smooth data') isolate(experiment$preMRNA_smooth) else isolate(experiment$preMRNA)
																							 , isolate(experiment$preMRNAsd)^2
																							 , maxit=300, method=isolate(input$opt_method)
																							 , isolate(experiment$no_nascent), mod_method = isolate(inspect$mod_method))
						## update mature values
						if( modeling_fun == 'sigmoid' ) {
							values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],2)
							values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],2)
							values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],2)
							values$k1_t1 = round(gene_model[[ 1 ]][["params"]][3],2)
							values$k1_t2 = round(gene_model[[ 1 ]][["params"]][3],2)
							values$k1_beta = round(gene_model[[ 1 ]][["params"]][4],2)
						} else {
							values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],2)
							values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],2)
							values$k1_h2 = round(gene_model[[ 1 ]][["params"]][3],2)
							values$k1_t1 = round(gene_model[[ 1 ]][["params"]][4],2)
							values$k1_t2 = round(gene_model[[ 1 ]][["params"]][5],2)
							values$k1_beta = round(gene_model[[ 1 ]][["params"]][6],2)
						}
						## update k2 values
						if( modeling_fun == 'sigmoid' ) {
							values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],2)
							values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],2)
							values$k2_h2 = round(gene_model[[ 3 ]][["params"]][2],2)
							values$k2_t1 = round(gene_model[[ 3 ]][["params"]][3],2)
							values$k2_t2 = round(gene_model[[ 3 ]][["params"]][3],2)
							values$k2_beta = round(gene_model[[ 3 ]][["params"]][4],2)
						} else {
							values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],2)
							values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],2)
							values$k2_h2 = round(gene_model[[ 3 ]][["params"]][3],2)
							values$k2_t1 = round(gene_model[[ 3 ]][["params"]][4],2)
							values$k2_t2 = round(gene_model[[ 3 ]][["params"]][5],2)
							values$k2_beta = round(gene_model[[ 3 ]][["params"]][6],2)
						}
						## update k3 values
						if( modeling_fun == 'sigmoid' ) {
							values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],2)
							values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],2)
							values$k3_h2 = round(gene_model[[ 2 ]][["params"]][2],2)
							values$k3_t1 = round(gene_model[[ 2 ]][["params"]][3],2)
							values$k3_t2 = round(gene_model[[ 2 ]][["params"]][3],2)
							values$k3_beta = round(gene_model[[ 2 ]][["params"]][4],2)
						} else {
							values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],2)
							values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],2)
							values$k3_h2 = round(gene_model[[ 2 ]][["params"]][3],2)
							values$k3_t1 = round(gene_model[[ 2 ]][["params"]][4],2)
							values$k3_t2 = round(gene_model[[ 2 ]][["params"]][5],2)
							values$k3_beta = round(gene_model[[ 2 ]][["params"]][6],2)
						}
						## update functional forms
						updateRadioButtons(session, 'k1_function',
															 selected = if(modeling_fun == 'impulse') "Impulsive" else "Sigmoidal")
						updateRadioButtons(session, 'k2_function',
															 selected = if(modeling_fun == 'impulse') "Impulsive" else "Sigmoidal")
						updateRadioButtons(session, 'k3_function',
															 selected = if(modeling_fun == 'impulse') "Impulsive" else "Sigmoidal")
					}
				} else if( modeling_type == 'TK12' ) { # KKV
					
					parameters <- unname(unlist(lapply(gene_model[1:3], '[[', 'params')))
					degradationfit <- k3KKV_Der(ids@tpts, parameters)
					degradationmodel <- .chooseModel(ids@tpts, experiment = degradationfit, variance = 1, 
																					sigmoid = TRUE, impulse = TRUE, polynomial = FALSE,
																					nInit = 10, nIter = 300, sigmoidModel=sigmoidModel, impulseModel=impulseModel,
																					sigmoidModelP=sigmoidModelP, impulseModelP=impulseModelP, .polynomialModelP=.polynomialModelP,
																					computeDerivatives = TRUE)
					maturefit <- matureKVK_Der(ids@tpts, parameters)
					maturemodel <- .chooseModel(ids@tpts, experiment = maturefit, variance = 1, 
																			sigmoid = degradationmodel$type == 'sigmoid', impulse = degradationmodel$type == 'impulse', polynomial = FALSE,
																			nInit = 10, nIter = 300, sigmoidModel=sigmoidModel, impulseModel=impulseModel,
																			sigmoidModelP=sigmoidModelP, impulseModelP=impulseModelP, .polynomialModelP=.polynomialModelP,
																			computeDerivatives = TRUE)
					k1_rate <- switch(degradationmodel$type, 
														"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																						 params=maturemodel$params,df=4),
														"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																						 params=maturemodel$params,df=6)
					)
					k2_rate <- list(type='constant',fun=newPointer(constantModel),
													params=gene_model[[3]]$params,df=1)
					k3_rate <- switch(degradationmodel$type, 
														"sigmoid" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																						 params=degradationmodel$params,df=4),
														"impulse" = list(type='impulse',fun=newPointer(impulseModel),
																						 params=degradationmodel$params,df=6)
					)
					params <- list(alpha=k1_rate, gamma=k2_rate, beta=k3_rate)
					gene_model <- optimParamsMatureRNA(params, ids@tpts
																						 , if(input$data_selection == 'Smooth data') isolate(experiment$synthesis_smooth) else isolate(experiment$synthesis)
																						 , isolate(experiment$synthesissd)^2
																						 , if(input$data_selection == 'Smooth data') isolate(experiment$mRNA_smooth) else isolate(experiment$mRNA)
																						 , isolate(experiment$mRNAsd)^2
																						 , if(input$data_selection == 'Smooth data') isolate(experiment$preMRNA_smooth) else isolate(experiment$preMRNA)
																						 , isolate(experiment$preMRNAsd)^2
																						 , maxit=300, method=isolate(input$opt_method)
																						 , isolate(experiment$no_nascent), mod_method = isolate(inspect$mod_method))
					## update mature values
					if( degradationmodel$type == 'sigmoid' ) {
						values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],2)
						values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],2)
						values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],2)
						values$k1_t1 = round(gene_model[[ 1 ]][["params"]][3],2)
						values$k1_t2 = round(gene_model[[ 1 ]][["params"]][3],2)
						values$k1_beta = round(gene_model[[ 1 ]][["params"]][4],2)
					} else {
						values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],2)
						values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],2)
						values$k1_h2 = round(gene_model[[ 1 ]][["params"]][3],2)
						values$k1_t1 = round(gene_model[[ 1 ]][["params"]][4],2)
						values$k1_t2 = round(gene_model[[ 1 ]][["params"]][5],2)
						values$k1_beta = round(gene_model[[ 1 ]][["params"]][6],2)
					}
					## update k2 values
					values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],2)
					values$k2_h1 = round(gene_model[[ 3 ]][["params"]][1],2)
					values$k2_h2 = round(gene_model[[ 3 ]][["params"]][1],2)
					values$k2_t1 = round(mean(out$t_pars))
					values$k2_t2 = round(mean(out$t_pars))
					values$k2_beta = round(mean(out$beta_pars))
					## update k3 values
					if( degradationmodel$type == 'sigmoid' ) {
						values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],2)
						values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],2)
						values$k3_h2 = round(gene_model[[ 2 ]][["params"]][2],2)
						values$k3_t1 = round(gene_model[[ 2 ]][["params"]][3],2)
						values$k3_t2 = round(gene_model[[ 2 ]][["params"]][3],2)
						values$k3_beta = round(gene_model[[ 2 ]][["params"]][4],2)
					} else {
						values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],2)
						values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],2)
						values$k3_h2 = round(gene_model[[ 2 ]][["params"]][3],2)
						values$k3_t1 = round(gene_model[[ 2 ]][["params"]][4],2)
						values$k3_t2 = round(gene_model[[ 2 ]][["params"]][5],2)
						values$k3_beta = round(gene_model[[ 2 ]][["params"]][6],2)
					}
					## update functional forms
					updateRadioButtons(session, 'k1_function',
														 selected = if(maturemodel$type == 'impulse') "Impulsive" else "Sigmoidal")
					updateRadioButtons(session, 'k2_function',
														 selected = "Constant")
					updateRadioButtons(session, 'k3_function',
														 selected = if(degradationmodel$type == 'impulse') "Impulsive" else "Sigmoidal")
				}
				
				# gene_h_vals <- c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2),
				# 	isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2),
				# 	isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2))
				gene_t_vals <- c(isolate(values$k1_t1),isolate(values$k1_t2),isolate(values$k2_t1),
												 isolate(values$k2_t2),isolate(values$k3_t1),isolate(values$k3_t2))
				gene_beta_vals <- c(isolate(values$k1_beta),isolate(values$k2_beta),isolate(values$k3_beta))
				
				ranges$k1_h_min <- min(c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2),out$k1_h_pars[1]))
				ranges$k1_h_max <- max(c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2),out$k1_h_pars[2]))
				ranges$k2_h_min <- min(c(isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2),out$k2_h_pars[1]))
				ranges$k2_h_max <- max(c(isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2),out$k2_h_pars[2]))
				ranges$k3_h_min <- min(c(isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2),out$k3_h_pars[1]))
				ranges$k3_h_max <- max(c(isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2),out$k3_h_pars[2]))
				ranges$t_min    <- min(c(gene_t_vals,out$t_pars[1]))
				ranges$t_max    <- max(c(gene_t_vals,out$t_pars[2]))
				ranges$beta_min <- min(c(gene_beta_vals,out$beta_pars[1]))
				ranges$beta_max <- max(c(gene_beta_vals,out$beta_pars[2]))
				
			}, silent=TRUE)
		}
		
	})
	
	## in case of derivative modeling, when one functional form is changed also the 
	## others must be updated

	## on change of k1
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) & !is.null(input$k1_function) )
			if( inspect$mod_method == 'der' & input$data_selection != 'User defined' ) {
				updateRadioButtons(session, 'k2_function',
													 selected = if( input$k1_function == "Impulsive" & isolate(input$k2_function) == "Sigmoidal" ) "Impulsive" else
													 	if( input$k1_function == "Sigmoidal" & isolate(input$k2_function) == "Impulsive") "Sigmoidal")
				updateRadioButtons(session, 'k3_function',
													 selected = if( input$k1_function == "Impulsive" & isolate(input$k3_function) == "Sigmoidal" ) "Impulsive" else
													 	if( input$k1_function == "Sigmoidal" & isolate(input$k3_function) == "Impulsive") "Sigmoidal")
			}
	})

	## on change of k2
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) & !is.null(input$k2_function) )
			if( inspect$mod_method == 'der' & input$data_selection != 'User defined' ) {
				updateRadioButtons(session, 'k1_function',
													 selected = if( input$k2_function == "Impulsive" & isolate(input$k1_function) == "Sigmoidal" ) "Impulsive" else
													 	if( input$k2_function == "Sigmoidal" & isolate(input$k1_function) == "Impulsive") "Sigmoidal")
				updateRadioButtons(session, 'k3_function',
													 selected = if( input$k2_function == "Impulsive" & isolate(input$k3_function) == "Sigmoidal" ) "Impulsive" else
													 	if( input$k2_function == "Sigmoidal" & isolate(input$k3_function) == "Impulsive") "Sigmoidal")
				output$convergence <- renderPrint({"not converged"})
			}
	})

	## on change of k3
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) & !is.null(input$k3_function) )
			if( inspect$mod_method == 'der' & input$data_selection != 'User defined' ) {
				updateRadioButtons(session, 'k1_function',
													 selected = if( input$k3_function == "Impulsive" & isolate(input$k1_function) == "Sigmoidal" ) "Impulsive" else
													 	if( input$k3_function == "Sigmoidal" & isolate(input$k1_function) == "Impulsive") "Sigmoidal")
				updateRadioButtons(session, 'k2_function',
													 selected = if( input$k3_function == "Impulsive" & isolate(input$k2_function) == "Sigmoidal" ) "Impulsive" else
													 	if( input$k3_function == "Sigmoidal" & isolate(input$k2_function) == "Impulsive") "Sigmoidal")
			}
	})

	## update convergence data

	observe({
		if( !is.null(modeling$counts) & !is.null(modeling$convergence) )
		output$convergence <- renderPrint({
			# paste(modeling$counts,'iters (', switch(as.character(modeling$convergence),
			# 	"0"="converged",
			# 	"1"="not converged",
			# 	"10"="degenerated"), ')')
			switch(as.character(modeling$convergence),
						 "0"="converged",
						 "1"="not converged",
						 "10"="degenerated",
						 "error")
		})
		})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$opt_method)) {
				output$convergence <- renderPrint({"not converged"})
			}
		}
	})

	###########################################
	## observe_parameters_change_by_user #########
	###########################################
	
	# observe({
	# 	if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
	# 		if(!is.null(input$k1_function) & !is.null(isolate(function_types$k1)))
	# 		{
	# 			if(input$k1_function != isolate(function_types$k1)) {
	# 				output$convergence <- renderPrint({"not converged"})
	# 				function_types$k1 <- input$k1_function
	# 			}
	# 		}
	# 	}
	# })
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_h0) & !is.null(isolate(values$k1_h0)))
			{
				if(input$k1_h0 != isolate(values$k1_h0)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k1_h0 <- input$k1_h0
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_h1) & !is.null(isolate(values$k1_h1)))
			{
				if(input$k1_h1 != isolate(values$k1_h1)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k1_h1 <- input$k1_h1
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_h2) & !is.null(isolate(values$k1_h2)))
			{
				if(input$k1_h2 != isolate(values$k1_h2)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k1_h2 <- input$k1_h2
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_t1) & !is.null(isolate(values$k1_t1)))
			{
				if(input$k1_t1 != isolate(values$k1_t1)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k1_t1 <- input$k1_t1
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_t2) & !is.null(isolate(values$k1_t2)))
			{
				if(input$k1_t2 != isolate(values$k1_t2)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k1_t2 <- input$k1_t2
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_beta) & !is.null(isolate(values$k1_beta)))
			{
				if(input$k1_beta != isolate(values$k1_beta)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k1_beta <- input$k1_beta
				}
			}
		}
	})	
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_h0) & !is.null(isolate(values$k2_h0)))
			{
				if(input$k2_h0 != isolate(values$k2_h0)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k2_h0 <- input$k2_h0
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_h1) & !is.null(isolate(values$k2_h1)))
			{
				if(input$k2_h1 != isolate(values$k2_h1)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k2_h1 <- input$k2_h1
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_h2) & !is.null(isolate(values$k2_h2)))
			{
				if(input$k2_h2 != isolate(values$k2_h2)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k2_h2 <- input$k2_h2
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_t1) & !is.null(isolate(values$k2_t1)))
			{
				if(input$k2_t1 != isolate(values$k2_t1)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k2_t1 <- input$k2_t1
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_t2) & !is.null(isolate(values$k2_t2)))
			{
				if(input$k2_t2 != isolate(values$k2_t2)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k2_t2 <- input$k2_t2
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_beta) & !is.null(isolate(values$k2_beta)))
			{
				if(input$k2_beta != isolate(values$k2_beta)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k2_beta <- input$k2_beta
				}
			}
		}
	})	
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_h0) & !is.null(isolate(values$k3_h0)))
			{
				if(input$k3_h0 != isolate(values$k3_h0)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k3_h0 <- input$k3_h0
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_h1) & !is.null(isolate(values$k3_h1)))
			{
				if(input$k3_h1 != isolate(values$k3_h1)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k3_h1 <- input$k3_h1
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_h2) & !is.null(isolate(values$k3_h2)))
			{
				if(input$k3_h2 != isolate(values$k3_h2)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k3_h2 <- input$k3_h2
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_t1) & !is.null(isolate(values$k3_t1)))
			{
				if(input$k3_t1 != isolate(values$k3_t1)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k3_t1 <- input$k3_t1
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_t2) & !is.null(isolate(values$k3_t2)))
			{
				if(input$k3_t2 != isolate(values$k3_t2)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k3_t2 <- input$k3_t2
				}
			}
		}
	})
	
	observe({
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_beta) & !is.null(isolate(values$k3_beta)))
			{
				if(input$k3_beta != isolate(values$k3_beta)) {
					updateCheckboxInput(session, "confint_checkbox", value = FALSE)
					output$convergence <- renderPrint({"not converged"})
					values$k3_beta <- input$k3_beta
				}
			}
		}
	})	

	##############################################
	## observe_parameters_change_by_user (end) #########
	##############################################
	
	## modeling checkbox

	# observe({
	# 	if( input$data_selection == 'User defined' ) {
	# 		updateCheckboxInput(session, "confint_checkbox", value = FALSE)
	# 		values$confint <- FALSE
	# 	}
	# })
	
	output$modeling_box <- renderUI({
		if( input$data_selection != 'User defined' & !experiment$steady_state ) {
			list(
				h5("goodness of fit (p-value):"),
				verbatimTextOutput("pchisq", TRUE),
				h5("Akaike information criterion:"),
				verbatimTextOutput("aic", TRUE),
				h5("minimization status:"),
				verbatimTextOutput("convergence", TRUE),
				fluidRow(
					column(4,radioButtons("opt_method", "method", 
						choices = c('NM','BFGS'), selected = 'NM')),
					column(4,numericInput("nIter", label = h5("iterations"), value = 100)),
					column(4,h5('Optimization'), actionButton("optimize", "Run"))
					)
				)
		}
		})

	output$modeling_type <- renderUI({
		if( inspect$mod_method == 'int' | input$data_selection == 'User defined' ) {
			h4('Integrative framework')
		} else {
			h4('Derivative framework')
		}
	})
	
	output$select_condition <- renderUI({
		if( experiment$steady_state ) {
			selectInput("select_condition", label = "Select condition", 
				choices = NULL, selected = NULL)
		}
		})

	output$select_class <- renderUI({
		if( !experiment$steady_state ) {
			selectInput("select_class", label = "Select class", 
				choices = NULL, selected = NULL)
		}
		})

	######################################################################
	######################################################################
	### set the interactive part of the UI: ranges and values of the  
	### widgets are static at the beginning but can be changed upon the  
	### import of the INSPEcT dataset
	######################################################################
	######################################################################
	
	## rate pvalues 
	
	output$rate_pvals <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(values$confint) ) {
			if( input$data_selection != 'User defined' & values$confint ) {
				# fluidRow(
				# 	column(2,h5("variability of k1 (p):")),
				# 	column(2,verbatimTextOutput("k1_p", TRUE)),
				# 	column(2,h5("k2 (p):")),
				# 	column(2,verbatimTextOutput("k2_p", TRUE)),
				# 	column(2,h5("k3 (p):")),
				# 	column(2,verbatimTextOutput("k3_p", TRUE))
				# )
				fluidRow(
					column(6,h5('Rate variablity (p-value)')),
					column(6,h5(paste(
						'k1:',
						signif(values$rate_p['k1'],2), 
						'k2:',
						signif(values$rate_p['k2'],2),
						'k3:',
						signif(values$rate_p['k3'],2)
						)))
				)
			} else {
				NULL
			}
		}
	})	
	
	## logarithmic time axis
	
	output$logtime_checkbox_ui <- renderUI({
		if( input$data_selection != 'User defined' & !experiment$steady_state ) {
			checkboxInput("logtime_checkbox", 
					label = "Space time logarithmically", 
					value = values$logtime)
		} else {
			NULL
		}
		})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) )
			values$logtime <- input$logtime_checkbox
	})
	
	## view/compute confidence intervals

	output$confint_checkbox_ui <- renderUI({
		if( input$data_selection != 'User defined' & !experiment$steady_state ) {
			checkboxInput("confint_checkbox",
										label = "Evaluate confidence intervals",
										value = values$confint)
		} else {
			NULL
		}
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) )
			values$confint <- input$confint_checkbox
	})

	################################  
	# widgets for for k1
	################################

	convert_model_name <- function(model_name) {
		if(is.null(model_name)) {
			return('')
		} else {
			switch(model_name,
						 'total'='total',
						 'mature'='mature',
						 'alpha'='synthesis',
						 'gamma'='processing',
						 'beta'='degradation'
			)
		}
	}
	
	output$fun1_name <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) ) {
			if( input$data_selection == 'User defined' | inspect$mod_method == 'int' )	{
				h3("synthesis")
			} else {
				h3("mature")
			}
		}
	})
	
	output$fun1_unit <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(values$model_names) ) {
			if( input$data_selection == 'User defined' | inspect$mod_method == 'int' )	{
				h4("(RPKMs/hour)")
			} else {
				h4("(RPKMs)")
			}
		}
	})
	
	output$function_type_k1 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			radioButtons('k1_function', 'Select function', 
				choices = c('Constant','Sigmoidal','Impulsive'), 
				selected = function_types$k1)
	})

	output$min_h_vals_k1 <- renderUI({
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) )
			numericInput("min_h_k1", label = h5("set min"), 
				value = ranges$k1_h_min, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$min_h_k1) )
			ranges$k1_h_min <- input$min_h_k1 
		})
	
	output$max_h_vals_k1 <- renderUI({
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) )
			numericInput("max_h_k1", label = h5("set max"), 
				value = ranges$k1_h_max, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$max_h_k1) )
			ranges$k1_h_max <- input$max_h_k1 
		})
	
	output$ui_k1 <- renderUI({
		
		ids <- contentsrea()		
		if( !is.null(ids) & !is.null(input$k1_function) & !is.null(function_types$k1) )
			switch(input$k1_function,
				"Constant" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = ranges$k1_h_min,
						max = ranges$k1_h_max,
						value = values$k1_h0,
						step = 0.01)
				),
				"Sigmoidal" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = ranges$k1_h_min,
						max = ranges$k1_h_max,
						value = values$k1_h0,
						step = 0.01),
					sliderInput("k1_h1",
						"final levels:",
						min = ranges$k1_h_min,
						max = ranges$k1_h_max,
						value = values$k1_h1,
						step = 0.01),
					sliderInput("k1_t1",
						"response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k1_t1,
						step = 0.01),
					sliderInput("k1_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k1_beta,
						step = 0.01)
				),
				"Impulsive" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = ranges$k1_h_min,
						max = ranges$k1_h_max,
						value = values$k1_h0,
						step = 0.01),
					sliderInput("k1_h1",
						"intermediate levels:",
						min = ranges$k1_h_min,
						max = ranges$k1_h_max,
						value = values$k1_h1,
						step = 0.01),
					sliderInput("k1_h2",
						"end levels:",
						min = ranges$k1_h_min,
						max = ranges$k1_h_max,
						value = values$k1_h2,
						step = 0.01),
					sliderInput("k1_t1",
						"first response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k1_t1,
						step = 0.01),
					sliderInput("k1_t2",
						"second response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k1_t2,
						step = 0.01),
					sliderInput("k1_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k1_beta,
						step = 0.01)
				)
						 
			)
		
	})

	################################  
	# widgets for for k2
	################################
	
	# output$fun2_name <- renderUI({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) ) {
	# 		if( input$data_selection == 'User defined' )	{
	# 			h3("processing")
	# 		} else {
	# 			h3( convert_model_name(values$model_names[2]) )
	# 		}
	# 	}
	# })
	
	output$function_type_k2 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			radioButtons('k2_function', 'Select function', 
				choices = c('Constant','Sigmoidal','Impulsive'), 
				selected = function_types$k2)
	})

	output$min_h_vals_k2 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("min_h_k2", label = h5("set min"), 
				value = ranges$k2_h_min, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$min_h_k2) )
			ranges$k2_h_min <- input$min_h_k2 
		})

	output$max_h_vals_k2 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("max_h_k2", label = h5("set max"), 
				value = ranges$k2_h_max, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$max_h_k2) )
				ranges$k2_h_max <- input$max_h_k2 
		})
	
	output$ui_k2 <- renderUI({
		
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$k2_function) & !is.null(function_types$k2) )
			switch(input$k2_function,
				"Constant" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = ranges$k2_h_min,
						max = ranges$k2_h_max,
						value = values$k2_h0,
						step = 0.01)
				),
				"Sigmoidal" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = ranges$k2_h_min,
						max = ranges$k2_h_max,
						value = values$k2_h0,
						step = 0.01),
					sliderInput("k2_h1",
						"final levels:",
						min = ranges$k2_h_min,
						max = ranges$k2_h_max,
						value = values$k2_h1,
						step = 0.01),
					sliderInput("k2_t1",
						"response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k2_t1,
						step = 0.01),
					sliderInput("k2_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k2_beta,
						step = 0.01)
				),
				"Impulsive" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = ranges$k2_h_min,
						max = ranges$k2_h_max,
						value = values$k2_h0,
						step = 0.01),
					sliderInput("k2_h1",
						"intermediate levels:",
						min = ranges$k2_h_min,
						max = ranges$k2_h_max,
						value = values$k2_h1,
						step = 0.01),
					sliderInput("k2_h2",
						"end levels:",
						min = ranges$k2_h_min,
						max = ranges$k2_h_max,
						value = values$k2_h2,
						step = 0.01),
					sliderInput("k2_t1",
						"first response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k2_t1,
						step = 0.01),
					sliderInput("k2_t2",
						"second response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k2_t2,
						step = 0.01),
					sliderInput("k2_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k2_beta,
						step = 0.01)
				)
						 
			)
		
	})
	
	
	################################  
	# widgets for for k3
	################################
	
	# output$fun3_name <- renderUI({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) ) {
	# 		if( input$data_selection == 'User defined' )	{
	# 			h3("degradation")
	# 		} else {
	# 			h3( convert_model_name(values$model_names[3]) )
	# 		}
	# 	}
	# })
	
	output$function_type_k3 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			radioButtons('k3_function', 'Select function', 
				choices = c('Constant','Sigmoidal','Impulsive'), 
				selected = function_types$k3)
	})

	output$min_h_vals_k3 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("min_h_k3", label = h5("set min"), 
				value = ranges$k3_h_min, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$min_h_k3) )
				ranges$k3_h_min <- input$min_h_k3 
		})
	
	output$max_h_vals_k3 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("max_h_k3", label = h5("set max"), 
				value = ranges$k3_h_max, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$max_h_k3) )
				ranges$k3_h_max <- input$max_h_k3 
		})

	output$ui_k3 <- renderUI({

		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$k3_function) & !is.null(function_types$k3) )
			switch(input$k3_function,
				"Constant" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = ranges$k3_h_min,
						max = ranges$k3_h_max,
						value = values$k3_h0,
						step = 0.01)
				),
				"Sigmoidal" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = ranges$k3_h_min,
						max = ranges$k3_h_max,
						value = values$k3_h0,
						step = 0.01),
					sliderInput("k3_h1",
						"final levels:",
						min = ranges$k3_h_min,
						max = ranges$k3_h_max,
						value = values$k3_h1,
						step = 0.01),
					sliderInput("k3_t1",
						"response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k3_t1,
						step = 0.01),
					sliderInput("k3_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k3_beta,
						step = 0.01)
				),
				"Impulsive" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = ranges$k3_h_min,
						max = ranges$k3_h_max,
						value = values$k3_h0,
						step = 0.01),
					sliderInput("k3_h1",
						"intermediate levels:",
						min = ranges$k3_h_min,
						max = ranges$k3_h_max,
						value = values$k3_h1,
						step = 0.01),
					sliderInput("k3_h2",
						"end levels:",
						min = ranges$k3_h_min,
						max = ranges$k3_h_max,
						value = values$k3_h2,
						step = 0.01),
					sliderInput("k3_t1",
						"first response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k3_t1,
						step = 0.01),
					sliderInput("k3_t2",
						"second response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k3_t2,
						step = 0.01),
					sliderInput("k3_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k3_beta,
						step = 0.01)
				)
						 
			)
		
	})
	
	######################################################################
	######################################################################
	### PLOT of the gene (and update of some parameters if succeded)
	######################################################################
	######################################################################

	## call the plot function when downloading the image
	output$saveRNAdynamicsPlotButton <- downloadHandler(
		filename =  function() {
			"RNAdynamics.pdf"
		},
		# content is a function with argument file. content writes the plot to the device
		content = function(file) {
			pdf(file) # open the pdf device
			RNAdynamicsAppPlot(
				data_selection = input$data_selection,
				show_logtime = input$logtime_checkbox,
				show_confint = input$confint_checkbox,
				logshift = inspect$logshift,
				linshift = inspect$linshift,
				time_min = ranges$time_min,
				time_max = ranges$time_max,
				experiment = experiment,
				simdata = modeling$simdata#,
				#mod_method = if( input$data_selection != 'User defined' ) inspect$mod_method else 'int'
			)
			dev.off()  # turn the device off
		}
	)

	## call the plot function when downloading the image
	output$saveRNAdynamicsDataButton <- downloadHandler(
		filename =  function() {
			"RNAdynamics.tsv"
		},
		# content is a function with argument file. content writes the plot to the device
		content = function(file) {
			write.table(modeling$simdata$sim, file=file, row.names = FALSE, quote = FALSE, sep = '\t')
		}
	)
	
	output$gene <- renderPlot({
		
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) & !is.null(input$k1_function) &
			!is.null(input$k2_function) & !is.null(input$k3_function) )

			suppressWarnings(try({

				k1_params <- switch(input$k1_function,
						"Constant" = input$k1_h0,
						"Sigmoidal" = c(input$k1_h0, input$k1_h1, input$k1_t1, input$k1_beta),
						"Impulsive" = c(input$k1_h0, input$k1_h1, input$k1_h2, input$k1_t1, input$k1_t2, input$k1_beta)
					)

				expected_length_k1_params <- switch(input$k1_function,
						"Constant" = 1, "Sigmoidal" = 4, "Impulsive" = 6)

				k2_params <- switch(input$k2_function,
						"Constant" = input$k2_h0,
						"Sigmoidal" = c(input$k2_h0, input$k2_h1, input$k2_t1, input$k2_beta),
						"Impulsive" = c(input$k2_h0, input$k2_h1, input$k2_h2, input$k2_t1, input$k2_t2, input$k2_beta)
					)

				expected_length_k2_params <- switch(input$k2_function,
						"Constant" = 1, "Sigmoidal" = 4, "Impulsive" = 6)

				k3_params <- switch(input$k3_function,
						"Constant" = input$k3_h0,
						"Sigmoidal" = c(input$k3_h0, input$k3_h1, input$k3_t1, input$k3_beta),
						"Impulsive" = c(input$k3_h0, input$k3_h1, input$k3_h2, input$k3_t1, input$k3_t2, input$k3_beta)
					)

				expected_length_k3_params <- switch(input$k3_function,
						"Constant" = 1, "Sigmoidal" = 4, "Impulsive" = 6)

				if( length(k1_params) == expected_length_k1_params &
					length(k2_params) == expected_length_k2_params &
					length(k3_params) == expected_length_k3_params )
				{

					if( input$data_selection != 'User defined' ) {
						mod_method <- inspect$mod_method
					} else {
						mod_method <- 'int'
					}

					simdata <- RNAdynamicsAppMake(
						data_selection = input$data_selection,
						show_confint = input$confint_checkbox,
						time_min = ranges$time_min,
						time_max = ranges$time_max,
						experiment = experiment,
						k1_function = input$k1_function, 
						k2_function = input$k2_function, 
						k3_function = input$k3_function,
						k1_params = k1_params,
						k2_params = k2_params,
						k3_params = k3_params,
						mod_method = mod_method
						)
					RNAdynamicsAppPlot(
						data_selection = input$data_selection,
						show_logtime = input$logtime_checkbox,
						show_confint = input$confint_checkbox,
						logshift = inspect$logshift,
						linshift = inspect$linshift,
						time_min = ranges$time_min,
						time_max = ranges$time_max,
						experiment = experiment,
						simdata = simdata
					)
					modeling$simdata <- simdata
					output$pchisq <- renderPrint({signif(modeling$simdata$scores$pchisq,3)})
					output$aic <- renderPrint({signif(modeling$simdata$scores$aic,3)})
					values$rate_p <- modeling$simdata$scores$rate_p # put into a global variable
					
				}

			}, silent = TRUE))
		
	})

	#######################
	## OPTIMIZE ###########
	#######################

	observeEvent(input$optimize, {

		# print('optimization started')

		# log_shift <- inspect$logshift
		# lin_shift <- inspect$linshift
		no_nascent <- experiment$no_nascent
		tpts_exp <- experiment$tpts
		alpha_exp <- if( input$data_selection == 'Experimental data' ) 
			experiment$synthesis else experiment$synthesis_smooth
		mature_exp <- if( input$data_selection == 'Experimental data' ) 
			experiment$mRNA else experiment$mRNA_smooth
		preMRNA_exp <- if( input$data_selection == 'Experimental data' ) 
			experiment$preMRNA else experiment$preMRNA_smooth
		alpha_var <- experiment$synthesissd^2
		mature_var <- experiment$mRNAsd^2
		preMRNA_var <- experiment$preMRNAsd^2

		k1_rate <- switch(input$k1_function, 
											"Constant" = list(type='constant',fun=newPointer(constantModel),
																				params=input$k1_h0,df=1),
											"Sigmoidal" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																				 params=c(input$k1_h0, input$k1_h1, input$k1_t1, input$k1_beta),df=4),
											"Impulsive" = list(type='impulse',fun=newPointer(impulseModel),
																				 params=c(input$k1_h0, input$k1_h1, input$k1_h2, input$k1_t1, input$k1_t2, input$k1_beta),df=6)
		)
		
		k2_rate <- switch(input$k2_function, 
											"Constant" = list(type='constant',fun=newPointer(constantModel),
																				params=input$k2_h0,df=1),
											"Sigmoidal" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																				 params=c(input$k2_h0, input$k2_h1, input$k2_t1, input$k2_beta),df=4),
											"Impulsive" = list(type='impulse',fun=newPointer(impulseModel),
																				 params=c(input$k2_h0, input$k2_h1, input$k2_h2, input$k2_t1, input$k2_t2, input$k2_beta),df=6)
		)
		
		k3_rate <- switch(input$k3_function, 
											"Constant" = list(type='constant',fun=newPointer(constantModel),
																				params=input$k3_h0,df=1),
											"Sigmoidal" = list(type='sigmoid',fun=newPointer(sigmoidModel),
																				 params=c(input$k3_h0, input$k3_h1, input$k3_t1, input$k3_beta),df=4),
											"Impulsive" = list(type='impulse',fun=newPointer(impulseModel),
																				 params=c(input$k3_h0, input$k3_h1, input$k3_h2, input$k3_t1, input$k3_t2, input$k3_beta),df=6)
		)

		params <- list(alpha=k1_rate, gamma=k2_rate, beta=k3_rate)
		gene_model <- optimParamsMatureRNA(params, tpts_exp, alpha_exp, alpha_var, mature_exp
		 	, mature_var, preMRNA_exp, preMRNA_var, maxit=input$nIter
		 	, method=input$opt_method#, log_shift, lin_shift
		 	, no_nascent, mod_method = inspect$mod_method)

		modeling$counts <- modeling$counts + gene_model$counts[1]
		modeling$convergence <- gene_model$convergence

		#########################################
		### update the parameters in the GUI #######
		#########################################

		if( input$k1_function == "Constant" ) {
			k1_h <- gene_model[["alpha"]][["params"]][1]
			values$k1_h0 = round(k1_h,2)
			values$k1_h1 = round(k1_h,2)
			values$k1_h2 = round(k1_h,2)
			values$k1_t1 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k1_t2 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k1_beta = round(mean(c(ranges$beta_min, ranges$beta_max)))
		}
		if( input$k1_function == "Sigmoidal" ) {
			values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],2)
			values$k1_h1 = round(gene_model[["alpha"]][["params"]][2],2)
			values$k1_h2 = round(gene_model[["alpha"]][["params"]][2],2)
			values$k1_t1 = round(gene_model[["alpha"]][["params"]][3],2)
			values$k1_t2 = round(gene_model[["alpha"]][["params"]][3],2)
			values$k1_beta = round(gene_model[["alpha"]][["params"]][4],2)
		}
		if( input$k1_function == "Impulsive" ) {
			values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],2)
			values$k1_h1 = round(gene_model[["alpha"]][["params"]][2],2)
			values$k1_h2 = round(gene_model[["alpha"]][["params"]][3],2)
			values$k1_t1 = round(gene_model[["alpha"]][["params"]][4],2)
			values$k1_t2 = round(gene_model[["alpha"]][["params"]][5],2)
			values$k1_beta = round(gene_model[["alpha"]][["params"]][6],2)
		}
		
		if( input$k2_function == "Constant" ) {
			values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],2)
			values$k2_h1 = round(gene_model[["gamma"]][["params"]][1],2)
			values$k2_h2 = round(gene_model[["gamma"]][["params"]][1],2)
			values$k2_t1 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k2_t2 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k2_beta = round(mean(c(ranges$beta_min, ranges$beta_max)))
		}
		if( input$k2_function == "Sigmoidal" ) {
			values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],2)
			values$k2_h1 = round(gene_model[["gamma"]][["params"]][2],2)
			values$k2_h2 = round(gene_model[["gamma"]][["params"]][2],2)
			values$k2_t1 = round(gene_model[["gamma"]][["params"]][3],2)
			values$k2_t2 = round(gene_model[["gamma"]][["params"]][3],2)
			values$k2_beta = round(gene_model[["gamma"]][["params"]][4],2)
		}
		if( input$k2_function == "Impulsive" ) {
			values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],2)
			values$k2_h1 = round(gene_model[["gamma"]][["params"]][2],2)
			values$k2_h2 = round(gene_model[["gamma"]][["params"]][3],2)
			values$k2_t1 = round(gene_model[["gamma"]][["params"]][4],2)
			values$k2_t2 = round(gene_model[["gamma"]][["params"]][5],2)
			values$k2_beta = round(gene_model[["gamma"]][["params"]][6],2)
		}
		
		if( input$k3_function == "Constant" ) {
			values$k3_h0 = round(gene_model[["beta"]][["params"]][1],2)
			values$k3_h1 = round(gene_model[["beta"]][["params"]][1],2)
			values$k3_h2 = round(gene_model[["beta"]][["params"]][1],2)
			values$k3_t1 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k3_t2 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k3_beta = round(mean(c(ranges$beta_min, ranges$beta_max)))
		}
		if( input$k3_function == "Sigmoidal" ) {
			values$k3_h0 = round(gene_model[["beta"]][["params"]][1],2)
			values$k3_h1 = round(gene_model[["beta"]][["params"]][2],2)
			values$k3_h2 = round(gene_model[["beta"]][["params"]][2],2)
			values$k3_t1 = round(gene_model[["beta"]][["params"]][3],2)
			values$k3_t2 = round(gene_model[["beta"]][["params"]][3],2)
			values$k3_beta = round(gene_model[["beta"]][["params"]][4],2)
		}
		if( input$k3_function == "Impulsive" ) {
			values$k3_h0 = round(gene_model[["beta"]][["params"]][1],2)
			values$k3_h1 = round(gene_model[["beta"]][["params"]][2],2)
			values$k3_h2 = round(gene_model[["beta"]][["params"]][3],2)
			values$k3_t1 = round(gene_model[["beta"]][["params"]][4],2)
			values$k3_t2 = round(gene_model[["beta"]][["params"]][5],2)
			values$k3_beta = round(gene_model[["beta"]][["params"]][6],2)
		}
		
		## update the ranges
		
		k1_h_pars <- c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2))
		if( min(k1_h_pars) < isolate(ranges$k1_h_min) ) {
			ranges$k1_h_min <- min(k1_h_pars)
		}
		if( max(k1_h_pars) > isolate(ranges$k1_h_max) ) {
			ranges$k1_h_max <- max(k1_h_pars)
		}
		
		k2_h_pars <- c(isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2))
		if( min(k2_h_pars) < isolate(ranges$k2_h_min) ) {
			ranges$k2_h_min <- min(k2_h_pars)
		}
		if( max(k2_h_pars) > isolate(ranges$k2_h_max) ) {
			ranges$k2_h_max <- max(k2_h_pars)
		}
		
		k3_h_pars <- c(isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2))
		if( min(k3_h_pars) < isolate(ranges$k3_h_min) ) {
			ranges$k3_h_min <- min(k3_h_pars)
		}
		if( max(k3_h_pars) > isolate(ranges$k3_h_max) ) {
			ranges$k3_h_max <- max(k3_h_pars)
		}
		
		t_pars <- c(isolate(values$k1_t1),isolate(values$k1_t2),isolate(values$k2_t1),
								isolate(values$k2_t2),isolate(values$k3_t1),isolate(values$k3_t2))
		if( min(t_pars) < isolate(ranges$t_min) ) {
			ranges$t_min <- min(t_pars)
		}
		if( max(t_pars) > isolate(ranges$t_max) ) {
			ranges$t_max <- max(t_pars)
		}
		
		beta_pars <- c(isolate(values$k1_beta),isolate(values$k2_beta),isolate(values$k3_beta))
		if( min(beta_pars) < isolate(ranges$beta_min) ) {
			ranges$beta_min <- min(beta_pars)
		}
		if( max(beta_pars) > isolate(ranges$beta_max) ) {
			ranges$beta_max <- max(beta_pars)
		}
		
		# print('optimization finished')

		})
	
})