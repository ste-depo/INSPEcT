INSPEcTGUIshinyAppServer <- function(input, output, session) {
	
	# allow loading datasets up to 50MB
	options(shiny.maxRequestSize=50*1024^2)
	
	# global variables
	ranges <- reactiveValues()
	values <- reactiveValues()
	experiment <- reactiveValues()
	inspect <- reactiveValues(loaded=FALSE)
	modeling <- reactiveValues()
	
	######################################################################
	######################################################################
	### import of the INSPEcT dataset and update of the imported values
	######################################################################
	######################################################################
	
	# load INSPEcT file
	
	contentsrea <- reactive({

		filename <- input$file1$datapath 
		if( is.null(filename) ) filename <- system.file(package='INSPEcT', 'INSPEcT_GUI_sample_dataset.rds')

		## load file
		if( file.exists(filename) ) {
			
			ids <- readRDS(filename)
			updateCheckboxInput(session, "fixyaxis_checkbox", value=FALSE)

			if( class(ids) != 'INSPEcT' ) {
				
				error_on_load_file()
				values$loaded_file <- FALSE
				values$loaded_file_error_message <- 'The loaded file is not of class INSPEcT'
				return(new(Class = 'INSPEcT'))

			}  else if(!(.hasSlot(ids, 'version'))) {

				error_on_load_file()
				values$loaded_file <- FALSE
				values$loaded_file_error_message <- 'The loaded file is obsolete and cannot work with the current version of INSPEcT'
				return(new(Class = 'INSPEcT'))
				
			} else { # the loaded object is of class INSPEcT
				
				## store inspect global values
				experiment$tpts <- tpts(ids)
				experiment$no_nascent <- ids@NoNascent
				experiment$steady_state <- is.character(experiment$tpts)
				
				# select only genes with exons and introns
				ids <- ids[apply(is.finite(ratesFirstGuess(ids,'preMRNA')),1,all)]

				if( !experiment$steady_state & nrow(ids@modelRates) > 0 ) {

					inspect$mod_method <- modelingParams(ids)$estimateRatesWith ## either "der" or "int"
					inspect$classes <- geneClass(ids)
					inspect$classes_internal <- geneClassInternal(ids)
					inspect$logshift <- findttpar(experiment$tpts)
					inspect$linshift <- ifelse( experiment$no_nascent,
						abs(min(timetransf(experiment$tpts,inspect$logshift))),0)

					# ... optionally throw a warning and consider as steady-state
					if( ids@NF ) stop("The RNAdynamics app doesn't work with non-functional(NF) INSPEcT models")
					
					# set default values of the checkboxes
					updateRadioButtons(session, inputId = 'data_selection', selected = 'Experimental data')
					
					## select only the best model
					
					if( ids@NoNascent ){ # (based on the class)
						ids@model@ratesSpecs <-
							lapply(seq_along(inspect$classes_internal), function(i)
								list(ids@model@ratesSpecs[[i]][[inspect$classes_internal[i]]]))
						names(ids@model@ratesSpecs) <- featureNames(ids)
					} else{ # Nascent RNA object (always VVV)
						ids@model@ratesSpecs <-
							lapply(seq_along(inspect$classes), function(i)
								list(ids@model@ratesSpecs[[i]][['abc']]))
						names(ids@model@ratesSpecs) <- featureNames(ids)
					}
					
					## update (converted) gene classes in the select input box
					classes_table <- sort(table(isolate(inspect$classes)), decreasing = TRUE)
					# names(classes_table) <- convert_gene_classes( names(classes_table) )
					classes_table_string <- paste( names(classes_table) , '(', classes_table, ')' )
					updateSelectInput(session, "select_class", 
						choices = classes_table_string, selected = classes_table_string[1])

					## define ranges

					ranges$time_min <- min(experiment$tpts)
					ranges$time_max <- max(experiment$tpts)

				} else { ## steady state

					# set default values of the checkboxes
					# updateRadioButtons(session, inputId = 'data_selection', selected = 'Experimental data')
					
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
				values$loaded_file <- TRUE
				
				return(ids)
			
				}
				
		} else { # (if the file name does not exist)
		
			error_on_load_file()
			values$loaded_file <- FALSE
			values$loaded_file_error_message <- 'The selected file does not exist'
			return(new(Class = 'INSPEcT'))
			
		}

	})
	
	#####################
	## load file error ##
	#####################
	
	error_on_load_file <- function() {
		updateRadioButtons(session, 'data_selection', selected = 'User defined')
		updateSelectInput(session, "select_class", choices = character(0))
		updateRadioButtons(session, 'k1_function', selected = 'Constant')
		updateRadioButtons(session, 'k2_function', selected = 'Constant')
		updateRadioButtons(session, 'k3_function', selected = 'Constant')
	}
	
	output$file_error <- renderUI({
		if( !is.null(contentsrea()) & !values$loaded_file ) {
			h6(paste('Error:', values$loaded_file_error_message))
		} else {
			NULL
		}
	})
	
	################################################
	## update gene names in the select input box ######
	################################################
	
	observe({
		if( !is.null(input$select_class) ) {
			if( !experiment$steady_state ) {
				# selected_class <- reconvert_gene_classes(strsplit( input$select_class , ' ')[[1]][1])
				selected_class <- strsplit( input$select_class , ' ')[[1]][1]
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
			
			# in case of time course experiment, proceed only when the gene and the 
			# regulation class conincide
			proceed <- FALSE
			if( experiment$steady_state ) {
				proceed <- TRUE	
			} else {
				gene_class <- strsplit( input$select_class , ' ')[[1]][1]
				if( experiment$steady_state | geneClass(ids)[input$select] == gene_class ) {
					proceed <- TRUE	
				}
			}
			
			if( proceed ) {

				## define the modeling strategy
				
				if( inspect$mod_method == 'int' ) {
					modeling_type <- 'K123'
					model_names <- c('alpha','gamma','beta')
				} else { # inspect$mod_method == 'der'
					if( experiment$no_nascent & (gene_class  %in% c('p','d','pd') ) ) {  ## c('KVK','KKV','KVV')
						if( gene_class %in% c('p','pd') ) {
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
					if( is.na(out$t_pars[1]) ) out$t_pars[1] <- isolate(ranges$time_min)
					if( is.na(out$t_pars[2]) ) out$t_pars[2] <- isolate(ranges$time_max)
					if( is.na(out$beta_pars[1]) ) out$beta_pars[1] <- 0
					if( is.na(out$beta_pars[2]) ) out$beta_pars[2] <- 10
					
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
						rate_range = quantile(rate_vals[is.finite(rate_vals)], probs=c(.025, .975))
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
				
				# try({ 
				# 	## this try is necessary because the class is updated before the gene, then the model_name
				# 	## corresponding to the new class can also not correspond to the currently selected gene
					
					function_types_k1 <- switch(
						gene_model[[ model_names[1] ]][['type']],
						"constant" = "Constant",
						"sigmoid" = "Sigmoidal",
						"impulse" = "Impulsive"
					)
					
					function_types_k2 <- switch(
						gene_model[[ model_names[2] ]][['type']],
						"constant" = "Constant",
						"sigmoid" = "Sigmoidal",
						"impulse" = "Impulsive"
					)
					
					function_types_k3 <- switch(
						gene_model[[ model_names[3] ]][['type']],
						"constant" = "Constant",
						"sigmoid" = "Sigmoidal",
						"impulse" = "Impulsive"
					)
					
					updateRadioButtons(session, 'k1_function',
														 selected = function_types_k1)
					updateRadioButtons(session, 'k2_function',
														 selected = function_types_k2)
					updateRadioButtons(session, 'k3_function',
														 selected = function_types_k3)
					
					if( function_types_k1 == "Constant" ) {
						values$k1_h0 = round(gene_model[[ model_names[1] ]][["params"]][1],6)
						values$k1_h1 = round(gene_model[[ model_names[1] ]][["params"]][1],6)
						values$k1_h2 = round(gene_model[[ model_names[1] ]][["params"]][1],6)
						values$k1_t1 = round(mean(out$t_pars))
						values$k1_t2 = round(mean(out$t_pars))
						values$k1_beta = round(mean(out$beta_pars))
					}
					if( function_types_k1 == "Sigmoidal" ) {
						values$k1_h0 = round(gene_model[[ model_names[1] ]][["params"]][1],6)
						values$k1_h1 = round(gene_model[[ model_names[1] ]][["params"]][2],6)
						values$k1_h2 = round(gene_model[[ model_names[1] ]][["params"]][2],6)
						values$k1_t1 = round(gene_model[[ model_names[1] ]][["params"]][3],6)
						values$k1_t2 = round(gene_model[[ model_names[1] ]][["params"]][3],6)
						values$k1_beta = round(gene_model[[ model_names[1] ]][["params"]][4],6)
					}
					if( function_types_k1 == "Impulsive" ) {
						values$k1_h0 = round(gene_model[[ model_names[1] ]][["params"]][1],6)
						values$k1_h1 = round(gene_model[[ model_names[1] ]][["params"]][2],6)
						values$k1_h2 = round(gene_model[[ model_names[1] ]][["params"]][3],6)
						values$k1_t1 = round(gene_model[[ model_names[1] ]][["params"]][4],6)
						values$k1_t2 = round(gene_model[[ model_names[1] ]][["params"]][5],6)
						values$k1_beta = round(gene_model[[ model_names[1] ]][["params"]][6],6)
					}
					
					if( function_types_k2 == "Constant" ) {
						values$k2_h0 = round(gene_model[[ model_names[2] ]][["params"]][1],6)
						values$k2_h1 = round(gene_model[[ model_names[2] ]][["params"]][1],6)
						values$k2_h2 = round(gene_model[[ model_names[2] ]][["params"]][1],6)
						values$k2_t1 = round(mean(out$t_pars))
						values$k2_t2 = round(mean(out$t_pars))
						values$k2_beta = round(mean(out$beta_pars))
					}
					if( function_types_k2 == "Sigmoidal" ) {
						values$k2_h0 = round(gene_model[[ model_names[2] ]][["params"]][1],6)
						values$k2_h1 = round(gene_model[[ model_names[2] ]][["params"]][2],6)
						values$k2_h2 = round(gene_model[[ model_names[2] ]][["params"]][2],6)
						values$k2_t1 = round(gene_model[[ model_names[2] ]][["params"]][3],6)
						values$k2_t2 = round(gene_model[[ model_names[2] ]][["params"]][3],6)
						values$k2_beta = round(gene_model[[ model_names[2] ]][["params"]][4],6)
					}
					if( function_types_k2 == "Impulsive" ) {
						values$k2_h0 = round(gene_model[[ model_names[2] ]][["params"]][1],6)
						values$k2_h1 = round(gene_model[[ model_names[2] ]][["params"]][2],6)
						values$k2_h2 = round(gene_model[[ model_names[2] ]][["params"]][3],6)
						values$k2_t1 = round(gene_model[[ model_names[2] ]][["params"]][4],6)
						values$k2_t2 = round(gene_model[[ model_names[2] ]][["params"]][5],6)
						values$k2_beta = round(gene_model[[ model_names[2] ]][["params"]][6],6)
					}
					
					if( function_types_k3 == "Constant" ) {
						values$k3_h0 = round(gene_model[[ model_names[3] ]][["params"]][1],6)
						values$k3_h1 = round(gene_model[[ model_names[3] ]][["params"]][1],6)
						values$k3_h2 = round(gene_model[[ model_names[3] ]][["params"]][1],6)
						values$k3_t1 = round(mean(out$t_pars))
						values$k3_t2 = round(mean(out$t_pars))
						values$k3_beta = round(mean(out$beta_pars))
					}
					if( function_types_k3 == "Sigmoidal" ) {
						values$k3_h0 = round(gene_model[[ model_names[3] ]][["params"]][1],6)
						values$k3_h1 = round(gene_model[[ model_names[3] ]][["params"]][2],6)
						values$k3_h2 = round(gene_model[[ model_names[3] ]][["params"]][2],6)
						values$k3_t1 = round(gene_model[[ model_names[3] ]][["params"]][3],6)
						values$k3_t2 = round(gene_model[[ model_names[3] ]][["params"]][3],6)
						values$k3_beta = round(gene_model[[ model_names[3] ]][["params"]][4],6)
					}
					if( function_types_k3 == "Impulsive" ) {
						values$k3_h0 = round(gene_model[[ model_names[3] ]][["params"]][1],6)
						values$k3_h1 = round(gene_model[[ model_names[3] ]][["params"]][2],6)
						values$k3_h2 = round(gene_model[[ model_names[3] ]][["params"]][3],6)
						values$k3_t1 = round(gene_model[[ model_names[3] ]][["params"]][4],6)
						values$k3_t2 = round(gene_model[[ model_names[3] ]][["params"]][5],6)
						values$k3_beta = round(gene_model[[ model_names[3] ]][["params"]][6],6)
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
								values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],6)
								values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_t1 = round(gene_model[[ 1 ]][["params"]][3],6)
								values$k1_t2 = round(gene_model[[ 1 ]][["params"]][3],6)
								values$k1_beta = round(gene_model[[ 1 ]][["params"]][4],6)
							} else {
								values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],6)
								values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_h2 = round(gene_model[[ 1 ]][["params"]][3],6)
								values$k1_t1 = round(gene_model[[ 1 ]][["params"]][4],6)
								values$k1_t2 = round(gene_model[[ 1 ]][["params"]][5],6)
								values$k1_beta = round(gene_model[[ 1 ]][["params"]][6],6)
							}
							## update k2 values
							if( processingmodel$type == 'sigmoid' ) {
								values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],6)
								values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_h2 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_t1 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_t2 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_beta = round(gene_model[[ 3 ]][["params"]][4],6)
							} else {
								values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],6)
								values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_h2 = round(gene_model[[ 3 ]][["params"]][3],6)
								values$k2_t1 = round(gene_model[[ 3 ]][["params"]][4],6)
								values$k2_t2 = round(gene_model[[ 3 ]][["params"]][5],6)
								values$k2_beta = round(gene_model[[ 3 ]][["params"]][6],6)
							}
							values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],6)
							values$k3_h1 = round(gene_model[[ 2 ]][["params"]][1],6)
							values$k3_h2 = round(gene_model[[ 2 ]][["params"]][1],6)
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
						} else 
						{ # KVV
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
								values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],6)
								values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_t1 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_t2 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_beta = round(gene_model[[ 1 ]][["params"]][4],6)
							} else {
								values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],6)
								values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],6)
								values$k1_t1 = round(gene_model[[ 1 ]][["params"]][4],6)
								values$k1_t2 = round(gene_model[[ 1 ]][["params"]][5],6)
								values$k1_beta = round(gene_model[[ 1 ]][["params"]][6],6)
							}
							## update k2 values
							if( modeling_fun == 'sigmoid' ) {
								values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],6)
								values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_h2 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_t1 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_t2 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_beta = round(gene_model[[ 3 ]][["params"]][4],6)
							} else {
								values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],6)
								values$k2_h1 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_h2 = round(gene_model[[ 3 ]][["params"]][2],6)
								values$k2_t1 = round(gene_model[[ 3 ]][["params"]][4],6)
								values$k2_t2 = round(gene_model[[ 3 ]][["params"]][5],6)
								values$k2_beta = round(gene_model[[ 3 ]][["params"]][6],6)
							}
							## update k3 values
							if( modeling_fun == 'sigmoid' ) {
								values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],6)
								values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],6)
								values$k3_h2 = round(gene_model[[ 2 ]][["params"]][2],6)
								values$k3_t1 = round(gene_model[[ 2 ]][["params"]][2],6)
								values$k3_t2 = round(gene_model[[ 2 ]][["params"]][2],6)
								values$k3_beta = round(gene_model[[ 2 ]][["params"]][4],6)
							} else {
								values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],6)
								values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],6)
								values$k3_h2 = round(gene_model[[ 2 ]][["params"]][2],6)
								values$k3_t1 = round(gene_model[[ 2 ]][["params"]][4],6)
								values$k3_t2 = round(gene_model[[ 2 ]][["params"]][5],6)
								values$k3_beta = round(gene_model[[ 2 ]][["params"]][6],6)
							}
							## update functional forms
							updateRadioButtons(session, 'k1_function',
																 selected = if(modeling_fun == 'impulse') "Impulsive" else "Sigmoidal")
							updateRadioButtons(session, 'k2_function',
																 selected = if(modeling_fun == 'impulse') "Impulsive" else "Sigmoidal")
							updateRadioButtons(session, 'k3_function',
																 selected = if(modeling_fun == 'impulse') "Impulsive" else "Sigmoidal")
						}
					} else if( modeling_type == 'TK12' ) 
					{ # KKV
						
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
							values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],6)
							values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],6)
							values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],6)
							values$k1_t1 = round(gene_model[[ 1 ]][["params"]][2],6)
							values$k1_t2 = round(gene_model[[ 1 ]][["params"]][2],6)
							values$k1_beta = round(gene_model[[ 1 ]][["params"]][4],6)
						} else {
							values$k1_h0 = round(gene_model[[ 1 ]][["params"]][1],6)
							values$k1_h1 = round(gene_model[[ 1 ]][["params"]][2],6)
							values$k1_h2 = round(gene_model[[ 1 ]][["params"]][2],6)
							values$k1_t1 = round(gene_model[[ 1 ]][["params"]][4],6)
							values$k1_t2 = round(gene_model[[ 1 ]][["params"]][5],6)
							values$k1_beta = round(gene_model[[ 1 ]][["params"]][6],6)
						}
						## update k2 values
						values$k2_h0 = round(gene_model[[ 3 ]][["params"]][1],6)
						values$k2_h1 = round(gene_model[[ 3 ]][["params"]][1],6)
						values$k2_h2 = round(gene_model[[ 3 ]][["params"]][1],6)
						values$k2_t1 = round(mean(out$t_pars))
						values$k2_t2 = round(mean(out$t_pars))
						values$k2_beta = round(mean(out$beta_pars))
						## update k3 values
						if( degradationmodel$type == 'sigmoid' ) {
							values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],6)
							values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],6)
							values$k3_h2 = round(gene_model[[ 2 ]][["params"]][2],6)
							values$k3_t1 = round(gene_model[[ 2 ]][["params"]][2],6)
							values$k3_t2 = round(gene_model[[ 2 ]][["params"]][2],6)
							values$k3_beta = round(gene_model[[ 2 ]][["params"]][4],6)
						} else {
							values$k3_h0 = round(gene_model[[ 2 ]][["params"]][1],6)
							values$k3_h1 = round(gene_model[[ 2 ]][["params"]][2],6)
							values$k3_h2 = round(gene_model[[ 2 ]][["params"]][2],6)
							values$k3_t1 = round(gene_model[[ 2 ]][["params"]][4],6)
							values$k3_t2 = round(gene_model[[ 2 ]][["params"]][5],6)
							values$k3_beta = round(gene_model[[ 2 ]][["params"]][6],6)
						}
						## update functional forms
						updateRadioButtons(session, 'k1_function',
															 selected = if(maturemodel$type == 'impulse') "Impulsive" else "Sigmoidal")
						updateRadioButtons(session, 'k2_function',
															 selected = "Constant")
						updateRadioButtons(session, 'k3_function',
															 selected = if(degradationmodel$type == 'impulse') "Impulsive" else "Sigmoidal")
					}
					
					gene_t_vals <- c(isolate(values$k1_t1),isolate(values$k1_t2),isolate(values$k2_t1),
													 isolate(values$k2_t2),isolate(values$k3_t1),isolate(values$k3_t2))
					gene_beta_vals <- c(isolate(values$k1_beta),isolate(values$k2_beta),isolate(values$k3_beta))
					
					ranges$k1_h_min <- floor(min(c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2),out$k1_h_pars[1])))
					try(updateNumericInput(session, "min_h_k1", value = isolate(ranges$k1_h_min)))
					ranges$k1_h_max <- ceiling(max(c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2),out$k1_h_pars[2])))
					try(updateNumericInput(session, "max_h_k1", value = isolate(ranges$k1_h_max)))
					ranges$k2_h_min <- floor(min(c(isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2),out$k2_h_pars[1])))
					try(updateNumericInput(session, "min_h_k2", value = isolate(ranges$k2_h_min)))
					ranges$k2_h_max <- ceiling(max(c(isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2),out$k2_h_pars[2])))
					try(updateNumericInput(session, "max_h_k2", value = isolate(ranges$k2_h_max)))
					ranges$k3_h_min <- floor(min(c(isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2),out$k3_h_pars[1])))
					try(updateNumericInput(session, "min_h_k3", value = isolate(ranges$k3_h_min)))
					ranges$k3_h_max <- ceiling(max(c(isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2),out$k3_h_pars[2])))
					try(updateNumericInput(session, "max_h_k3", value = isolate(ranges$k3_h_max)))
					ranges$t_min    <- floor(min(c(gene_t_vals,out$t_pars[1])))
					# try(updateNumericInput(session, "t_min", value = ranges$t_min))
					ranges$t_max    <- ceiling(max(c(gene_t_vals,out$t_pars[2])))
					# try(updateNumericInput(session, "t_max", value = ranges$t_max))
					ranges$beta_min <- floor(min(c(gene_beta_vals,out$beta_pars[1])))
					# try(updateNumericInput(session, "beta_min", value = ranges$beta_min))
					ranges$beta_max <- ceiling(max(c(gene_beta_vals,out$beta_pars[2])))
					# try(updateNumericInput(session, "beta_max", value = ranges$beta_max))
					
				# }, silent = TRUE)
								
			}

		}
		
	})
	
	## in case of derivative modeling, when one functional form is changed also 
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
	
	observeEvent(input$k1_h0, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_h0) & !is.null(isolate(values$k1_h0))) {
				if(input$k1_h0 != isolate(values$k1_h0)) {
					output$convergence <- renderPrint({"not converged"})
					values$k1_h0 <- input$k1_h0
				}
			}
		}
	})
	
	observeEvent(input$k1_h1, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_h1) & !is.null(isolate(values$k1_h1))) {
				if(input$k1_h1 != isolate(values$k1_h1)) {
					output$convergence <- renderPrint({"not converged"})
					values$k1_h1 <- input$k1_h1
				}
			}
		}
	})
	
	observeEvent(input$k1_h2, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_h2) & !is.null(isolate(values$k1_h2))) {
				if(input$k1_h2 != isolate(values$k1_h2)) {
					output$convergence <- renderPrint({"not converged"})
					values$k1_h2 <- input$k1_h2
				}
			}
		}
	})
	
	observeEvent(input$k1_t1, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_t1) & !is.null(isolate(values$k1_t1))) {
				if(input$k1_t1 != isolate(values$k1_t1)) {
					output$convergence <- renderPrint({"not converged"})
					values$k1_t1 <- input$k1_t1
				}
			}
		}
	})
	
	observeEvent(input$k1_t2, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_t2) & !is.null(isolate(values$k1_t2))) {
				if(input$k1_t2 != isolate(values$k1_t2)) {
					output$convergence <- renderPrint({"not converged"})
					values$k1_t2 <- input$k1_t2
				}
			}
		}
	})
	
	observeEvent(input$k1_beta, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k1_beta) & !is.null(isolate(values$k1_beta))) {
				if(input$k1_beta != isolate(values$k1_beta)) {
					output$convergence <- renderPrint({"not converged"})
					values$k1_beta <- input$k1_beta
				}
			}
		}
	})
	
	observeEvent(input$k2_h0, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_h0) & !is.null(isolate(values$k2_h0))) {
				if(input$k2_h0 != isolate(values$k2_h0)) {
					output$convergence <- renderPrint({"not converged"})
					values$k2_h0 <- input$k2_h0
				}
			}
		}
	})
	
	observeEvent(input$k2_h1, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_h1) & !is.null(isolate(values$k2_h1))) {
				if(input$k2_h1 != isolate(values$k2_h1)) {
					output$convergence <- renderPrint({"not converged"})
					values$k2_h1 <- input$k2_h1
				}
			}
		}
	})
	
	observeEvent(input$k2_h2, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_h2) & !is.null(isolate(values$k2_h2))) {
				if(input$k2_h2 != isolate(values$k2_h2)) {
					output$convergence <- renderPrint({"not converged"})
					values$k2_h2 <- input$k2_h2
				}
			}
		}
	})
	
	observeEvent(input$k2_t1, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_t1) & !is.null(isolate(values$k2_t1))) {
				if(input$k2_t1 != isolate(values$k2_t1)) {
					output$convergence <- renderPrint({"not converged"})
					values$k2_t1 <- input$k2_t1
				}
			}
		}
	})
	
	observeEvent(input$k2_t2, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_t2) & !is.null(isolate(values$k2_t2))) {
				if(input$k2_t2 != isolate(values$k2_t2)) {
					output$convergence <- renderPrint({"not converged"})
					values$k2_t2 <- input$k2_t2
				}
			}
		}
	})
	
	observeEvent(input$k2_beta, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k2_beta) & !is.null(isolate(values$k2_beta))) {
				if(input$k2_beta != isolate(values$k2_beta)) {
					output$convergence <- renderPrint({"not converged"})
					values$k2_beta <- input$k2_beta
				}
			}
		}
	})
	
	observeEvent(input$k3_h0, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_h0) & !is.null(isolate(values$k3_h0))) {
				if(input$k3_h0 != isolate(values$k3_h0)) {
					output$convergence <- renderPrint({"not converged"})
					values$k3_h0 <- input$k3_h0
				}
			}
		}
	})
	
	observeEvent(input$k3_h1, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_h1) & !is.null(isolate(values$k3_h1))) {
				if(input$k3_h1 != isolate(values$k3_h1)) {
					output$convergence <- renderPrint({"not converged"})
					values$k3_h1 <- input$k3_h1
				}
			}
		}
	})
	
	observeEvent(input$k3_h2, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_h2) & !is.null(isolate(values$k3_h2))) {
				if(input$k3_h2 != isolate(values$k3_h2)) {
					output$convergence <- renderPrint({"not converged"})
					values$k3_h2 <- input$k3_h2
				}
			}
		}
	})
	
	observeEvent(input$k3_t1, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_t1) & !is.null(isolate(values$k3_t1))) {
				if(input$k3_t1 != isolate(values$k3_t1)) {
					output$convergence <- renderPrint({"not converged"})
					values$k3_t1 <- input$k3_t1
				}
			}
		}
	})
	
	observeEvent(input$k3_t2, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_t2) & !is.null(isolate(values$k3_t2))) {
				if(input$k3_t2 != isolate(values$k3_t2)) {
					output$convergence <- renderPrint({"not converged"})
					values$k3_t2 <- input$k3_t2
				}
			}
		}
	})
	
	observeEvent(input$k3_beta, {
		if( !is.null(inspect$mod_method) & !is.null(input$data_selection) ) {
			if(!is.null(input$k3_beta) & !is.null(isolate(values$k3_beta))) {
				if(input$k3_beta != isolate(values$k3_beta)) {
					output$convergence <- renderPrint({"not converged"})
					values$k3_beta <- input$k3_beta
				}
			}
		}
	})
	
	##############################################
	## observe_parameters_change_by_user (end) #########
	##############################################
	
	## confidence intervals (only for derivative, non steady state, data-driven mode)
	
	output$confint_box <- renderUI({
		if( input$data_selection != 'User defined' & !experiment$steady_state & inspect$mod_method == 'der' ) {
			actionButton("conf_int_button", "Rate variability p-values")
		}
	})
	
	## modeling checkbox
	
	output$modeling_box <- renderUI({
		if( input$data_selection != 'User defined' & !experiment$steady_state ) {
			list(
				h4("Modeling box"),
				h5("goodness of fit (p-value)"),
				verbatimTextOutput("pchisq", TRUE),
				h5("Akaike information criterion"),
				verbatimTextOutput("aic", TRUE),
				h5("minimization status"),
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
	
	## rate pvalues 
	
	output$modeling_type <- renderUI({
		if( !is.null(contentsrea()) & input$data_selection != 'User defined' ) {
			if( experiment$steady_state ) {
				p(paste('Loaded a steady-state INSPEcT object with',nGenes(contentsrea()),'genes and',nTpts(contentsrea()),'conditions'))
			} else {
				if( inspect$mod_method == 'int' ) {
					p(paste('Loaded INSPEcT object with',nGenes(contentsrea()),'genes modeled with the integrative framework'))
				} else {
					p(paste('Loaded INSPEcT object with',nGenes(contentsrea()),'genes modeled with the derivative framework'))
				}
			}
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
	
	## logarithmic time axis
	
	output$logtime_checkbox_ui <- renderUI({
		if( input$data_selection != 'User defined' & !experiment$steady_state ) {
			checkboxInput("logtime_checkbox", 
					label = "Log time", 
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
				selected = NULL)
	})

	output$min_h_vals_k1 <- renderUI({
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) )
			numericInput("min_h_k1", label = h5("set min"), 
				value = ranges$k1_h_min, width='200px')
	})

	# observe({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) & !is.null(input$min_h_k1) )
	# 		ranges$k1_h_min <- input$min_h_k1 
	# 	})
	
	output$max_h_vals_k1 <- renderUI({
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) )
			numericInput("max_h_k1", label = h5("set max"), 
				value = ranges$k1_h_max, width='200px')
	})

	# observe({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) & !is.null(input$max_h_k1) )
	# 		ranges$k1_h_max <- input$max_h_k1 
	# 	})
	
	output$ui_k1 <- renderUI({
		
		ids <- contentsrea()		
		if( !is.null(ids) & !is.null(input$k1_function) ) {
			
			switch(input$k1_function,
				"Constant" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = input$min_h_k1,
						max = input$max_h_k1,
						value = values$k1_h0,
						step = 0.001)
				),
				"Sigmoidal" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = input$min_h_k1,
						max = input$max_h_k1,
						value = values$k1_h0,
						step = 0.001),
					sliderInput("k1_h1",
						"final levels:",
						min = input$min_h_k1,
						max = input$max_h_k1,
						value = values$k1_h1,
						step = 0.001),
					sliderInput("k1_t1",
						"response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k1_t1,
						step = 0.001),
					sliderInput("k1_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k1_beta,
						step = 0.001)
				),
				"Impulsive" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = input$min_h_k1,
						max = input$max_h_k1,
						value = values$k1_h0,
						step = 0.001),
					sliderInput("k1_h1",
						"intermediate levels:",
						min = input$min_h_k1,
						max = input$max_h_k1,
						value = values$k1_h1,
						step = 0.001),
					sliderInput("k1_h2",
						"end levels:",
						min = input$min_h_k1,
						max = input$max_h_k1,
						value = values$k1_h2,
						step = 0.001),
					sliderInput("k1_t1",
						"first response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k1_t1,
						step = 0.001),
					sliderInput("k1_t2",
						"second response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k1_t2,
						step = 0.001),
					sliderInput("k1_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k1_beta,
						step = 0.001)
				)
						 
			)
			
		}
		
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
				selected = NULL)
	})

	output$min_h_vals_k2 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("min_h_k2", label = h5("set min"), 
				value = ranges$k2_h_min, width='200px')
	})

	# observe({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) & !is.null(input$min_h_k2) )
	# 		ranges$k2_h_min <- input$min_h_k2 
	# 	})

	output$max_h_vals_k2 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("max_h_k2", label = h5("set max"), 
				value = ranges$k2_h_max, width='200px')
	})

	# observe({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) & !is.null(input$max_h_k2) )
	# 			ranges$k2_h_max <- input$max_h_k2 
	# 	})
	
	output$ui_k2 <- renderUI({
		
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$k2_function) )
			switch(input$k2_function,
				"Constant" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = input$min_h_k2,
						max = input$max_h_k2,
						value = values$k2_h0,
						step = 0.001)
				),
				"Sigmoidal" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = input$min_h_k2,
						max = input$max_h_k2,
						value = values$k2_h0,
						step = 0.001),
					sliderInput("k2_h1",
						"final levels:",
						min = input$min_h_k2,
						max = input$max_h_k2,
						value = values$k2_h1,
						step = 0.001),
					sliderInput("k2_t1",
						"response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k2_t1,
						step = 0.001),
					sliderInput("k2_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k2_beta,
						step = 0.001)
				),
				"Impulsive" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = input$min_h_k2,
						max = input$max_h_k2,
						value = values$k2_h0,
						step = 0.001),
					sliderInput("k2_h1",
						"intermediate levels:",
						min = input$min_h_k2,
						max = input$max_h_k2,
						value = values$k2_h1,
						step = 0.001),
					sliderInput("k2_h2",
						"end levels:",
						min = input$min_h_k2,
						max = input$max_h_k2,
						value = values$k2_h2,
						step = 0.001),
					sliderInput("k2_t1",
						"first response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k2_t1,
						step = 0.001),
					sliderInput("k2_t2",
						"second response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k2_t2,
						step = 0.001),
					sliderInput("k2_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k2_beta,
						step = 0.001)
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
				selected = NULL)
	})

	output$min_h_vals_k3 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("min_h_k3", label = h5("set min"), 
				value = ranges$k3_h_min, width='200px')
	})

	# observe({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) & !is.null(input$min_h_k3) )
	# 			ranges$k3_h_min <- input$min_h_k3 
	# 	})
	
	output$max_h_vals_k3 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("max_h_k3", label = h5("set max"), 
				value = ranges$k3_h_max, width='200px')
	})

	# observe({
	# 	ids <- contentsrea()
	# 	if( !is.null(ids) & !is.null(input$max_h_k3) )
	# 			ranges$k3_h_max <- input$max_h_k3 
	# 	})

	output$ui_k3 <- renderUI({

		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$k3_function) )
			switch(input$k3_function,
				"Constant" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = input$min_h_k3,
						max = input$max_h_k3,
						value = values$k3_h0,
						step = 0.001)
				),
				"Sigmoidal" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = input$min_h_k3,
						max = input$max_h_k3,
						value = values$k3_h0,
						step = 0.001),
					sliderInput("k3_h1",
						"final levels:",
						min = input$min_h_k3,
						max = input$max_h_k3,
						value = values$k3_h1,
						step = 0.001),
					sliderInput("k3_t1",
						"response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k3_t1,
						step = 0.001),
					sliderInput("k3_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k3_beta,
						step = 0.001)
				),
				"Impulsive" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = input$min_h_k3,
						max = input$max_h_k3,
						value = values$k3_h0,
						step = 0.001),
					sliderInput("k3_h1",
						"intermediate levels:",
						min = input$min_h_k3,
						max = input$max_h_k3,
						value = values$k3_h1,
						step = 0.001),
					sliderInput("k3_h2",
						"end levels:",
						min = input$min_h_k3,
						max = input$max_h_k3,
						value = values$k3_h2,
						step = 0.001),
					sliderInput("k3_t1",
						"first response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k3_t1,
						step = 0.001),
					sliderInput("k3_t2",
						"second response time:",
						min = ranges$t_min,
						max = ranges$t_max,
						value = values$k3_t2,
						step = 0.001),
					sliderInput("k3_beta",
						"slope:",
						min = ranges$beta_min,
						max = ranges$beta_max,
						value = values$k3_beta,
						step = 0.001)
				)
						 
			)
		
	})
	
	###########################
	## confidence intervals ###
	###########################

	observeEvent(input$conf_int_button, {
		
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
					
					ci_res <- RNAdynamicsAppMakeConfInt(
						data_selection = input$data_selection,
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
					
					# put into global variables
					simdata <- modeling$simdata
					simdata$conf_int <- ci_res$conf_int
					modeling$simdata <- simdata
					values$rate_p <- ci_res$rate_p 
					
				}
				
			}, silent = TRUE))
		
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
			suppressWarnings(try({
				RNAdynamicsAppPlot(
					data_selection = input$data_selection,
					show_logtime = input$logtime_checkbox,
					show_relexpr = input$relativexpr_checkbox,
					logshift = inspect$logshift,
					linshift = inspect$linshift,
					time_min = ranges$time_min,
					time_max = ranges$time_max,
					experiment = experiment,
					simdata = modeling$simdata,
					ylims = if(input$fixyaxis_checkbox) isolate(ranges$ylims) else NULL,
					rate_p = values$rate_p
				)
			}, silent = TRUE))
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
	
	observe({
		
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
				
					modeling$simdata <- simdata
					output$pchisq <- renderPrint({signif(isolate(modeling$simdata$scores$pchisq),3)})
					output$aic <- renderPrint({signif(isolate(modeling$simdata$scores$aic),3)})
					values$rate_p <- NULL

				}

			}, silent = TRUE))
		
	})

	output$gene <- renderPlot({
		suppressWarnings(try({
			if( !is.null(modeling$simdata) ) {
				ylims <- RNAdynamicsAppPlot(
					data_selection = input$data_selection,
					show_logtime = input$logtime_checkbox,
					show_relexpr = input$relativexpr_checkbox,
					logshift = inspect$logshift,
					linshift = inspect$linshift,
					time_min = ranges$time_min,
					time_max = ranges$time_max,
					experiment = experiment,
					simdata = modeling$simdata,
					ylims = if(input$fixyaxis_checkbox) isolate(ranges$ylims) else NULL,
					rate_p = values$rate_p
				)
				ranges$ylims <- ylims
			}
		}, silent = TRUE))
	})
	
	#######################
	## OPTIMIZE ###########
	#######################

	observeEvent(input$optimize, {

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
			values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],6)
			values$k1_h1 = round(gene_model[["alpha"]][["params"]][2],6)
			values$k1_h2 = round(gene_model[["alpha"]][["params"]][2],6)
			values$k1_t1 = round(gene_model[["alpha"]][["params"]][2],6)
			values$k1_t2 = round(gene_model[["alpha"]][["params"]][2],6)
			values$k1_beta = round(gene_model[["alpha"]][["params"]][4],6)
		}
		if( input$k1_function == "Impulsive" ) {
			values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],6)
			values$k1_h1 = round(gene_model[["alpha"]][["params"]][2],6)
			values$k1_h2 = round(gene_model[["alpha"]][["params"]][3],6)
			values$k1_t1 = round(gene_model[["alpha"]][["params"]][4],6)
			values$k1_t2 = round(gene_model[["alpha"]][["params"]][5],6)
			values$k1_beta = round(gene_model[["alpha"]][["params"]][6],6)
		}
		
		if( input$k2_function == "Constant" ) {
			values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],6)
			values$k2_h1 = round(gene_model[["gamma"]][["params"]][1],6)
			values$k2_h2 = round(gene_model[["gamma"]][["params"]][1],6)
			values$k2_t1 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k2_t2 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k2_beta = round(mean(c(ranges$beta_min, ranges$beta_max)))
		}
		if( input$k2_function == "Sigmoidal" ) {
			values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],6)
			values$k2_h1 = round(gene_model[["gamma"]][["params"]][2],6)
			values$k2_h2 = round(gene_model[["gamma"]][["params"]][2],6)
			values$k2_t1 = round(gene_model[["gamma"]][["params"]][2],6)
			values$k2_t2 = round(gene_model[["gamma"]][["params"]][2],6)
			values$k2_beta = round(gene_model[["gamma"]][["params"]][4],6)
		}
		if( input$k2_function == "Impulsive" ) {
			values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],6)
			values$k2_h1 = round(gene_model[["gamma"]][["params"]][2],6)
			values$k2_h2 = round(gene_model[["gamma"]][["params"]][3],6)
			values$k2_t1 = round(gene_model[["gamma"]][["params"]][4],6)
			values$k2_t2 = round(gene_model[["gamma"]][["params"]][5],6)
			values$k2_beta = round(gene_model[["gamma"]][["params"]][6],6)
		}
		
		if( input$k3_function == "Constant" ) {
			values$k3_h0 = round(gene_model[["beta"]][["params"]][1],6)
			values$k3_h1 = round(gene_model[["beta"]][["params"]][1],6)
			values$k3_h2 = round(gene_model[["beta"]][["params"]][1],6)
			values$k3_t1 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k3_t2 = round(mean(c(ranges$t_min, ranges$t_max)))
			values$k3_beta = round(mean(c(ranges$beta_min, ranges$beta_max)))
		}
		if( input$k3_function == "Sigmoidal" ) {
			values$k3_h0 = round(gene_model[["beta"]][["params"]][1],6)
			values$k3_h1 = round(gene_model[["beta"]][["params"]][2],6)
			values$k3_h2 = round(gene_model[["beta"]][["params"]][2],6)
			values$k3_t1 = round(gene_model[["beta"]][["params"]][2],6)
			values$k3_t2 = round(gene_model[["beta"]][["params"]][2],6)
			values$k3_beta = round(gene_model[["beta"]][["params"]][4],6)
		}
		if( input$k3_function == "Impulsive" ) {
			values$k3_h0 = round(gene_model[["beta"]][["params"]][1],6)
			values$k3_h1 = round(gene_model[["beta"]][["params"]][2],6)
			values$k3_h2 = round(gene_model[["beta"]][["params"]][3],6)
			values$k3_t1 = round(gene_model[["beta"]][["params"]][4],6)
			values$k3_t2 = round(gene_model[["beta"]][["params"]][5],6)
			values$k3_beta = round(gene_model[["beta"]][["params"]][6],6)
		}
		
		## update the ranges
		
		k1_h_pars <- c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2))
		if( min(k1_h_pars) < isolate(ranges$k1_h_min) ) {
			ranges$k1_h_min <- floor(min(k1_h_pars))
		}
		if( max(k1_h_pars) > isolate(ranges$k1_h_max) ) {
			ranges$k1_h_max <- ceiling(max(k1_h_pars))
		}

		k2_h_pars <- c(isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2))
		if( min(k2_h_pars) < isolate(ranges$k2_h_min) ) {
			ranges$k2_h_min <- floor(min(k2_h_pars))
		}
		if( max(k2_h_pars) > isolate(ranges$k2_h_max) ) {
			ranges$k2_h_max <- ceiling(max(k2_h_pars))
		}

		k3_h_pars <- c(isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2))
		if( min(k3_h_pars) < isolate(ranges$k3_h_min) ) {
			ranges$k3_h_min <- floor(min(k3_h_pars))
		}
		if( max(k3_h_pars) > isolate(ranges$k3_h_max) ) {
			ranges$k3_h_max <- ceiling(max(k3_h_pars))
		}
		
		t_pars <- c(isolate(values$k1_t1),isolate(values$k1_t2),isolate(values$k2_t1),
								isolate(values$k2_t2),isolate(values$k3_t1),isolate(values$k3_t2))
		if( min(t_pars) < isolate(ranges$t_min) ) {
			ranges$t_min <- floor(min(t_pars))
		}
		if( max(t_pars) > isolate(ranges$t_max) ) {
			ranges$t_max <- ceiling(max(t_pars))
		}
		
		beta_pars <- c(isolate(values$k1_beta),isolate(values$k2_beta),isolate(values$k3_beta))
		if( min(beta_pars) < isolate(ranges$beta_min) ) {
			ranges$beta_min <- floor(min(beta_pars))
		}
		if( max(beta_pars) > isolate(ranges$beta_max) ) {
			ranges$beta_max <- ceiling(max(beta_pars))
		}
		
		})
	
}