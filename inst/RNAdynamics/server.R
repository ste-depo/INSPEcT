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
				inspect$classes <- geneClass(ids)
				inspect$logshift <- find_tt_par(experiment$tpts)

				## remove models from other classes
				ids@model@ratesSpecs <- 
					lapply(seq_along(inspect$classes), function(i) 
						list(ids@model@ratesSpecs[[i]][[inspect$classes[i]]]))
				names(ids@model@ratesSpecs) <- featureNames(ids)

				## update (converted) gene classes in the select input box
				classes_table <- table(isolate(inspect$classes))
				names(classes_table) <- convert_gene_classes( names(classes_table) )
				classes_table_string <- paste( names(classes_table) , '(', classes_table, ')' )
				updateSelectInput(session, "select_class", 
					choices = classes_table_string, selected = classes_table_string[1])

				## define ranges

				ranges$time_min <- min(experiment$tpts)
				ranges$time_max <- max(experiment$tpts)

				## predifined logtime 

				values$logtime <- FALSE

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
		selected_class <- reconvert_gene_classes(strsplit( input$select_class , ' ')[[1]][1])
		updateSelectInput(session, "select", selected = NULL,
			choices = sort(featureNames(contentsrea()[inspect$classes == selected_class])))
	})
	
	############################################################################
	## select parameters of each rate (synthesis, processing and degradation) #####
	## and experimental values for the selected gene ###########################
	####################################################
	
	observe({

		ids <- contentsrea()
				
		if( input$select != "" & !is.null(ids) &
			input$select %in% featureNames(ids) ) {

			experiment$synthesis <- ratesFirstGuess(ids[input$select], 'synthesis')
			experiment$mRNA <- ratesFirstGuess(ids[input$select], 'total') - ratesFirstGuess(ids[input$select], 'preMRNA')
			experiment$preMRNA <- ratesFirstGuess(ids[input$select], 'preMRNA')
			
			## test on smooth experiment data???
			experiment$synthesis_smooth <- smoothModel(ids@tpts, ratesFirstGuess(ids[input$select], 'synthesis'))
			experiment$mRNA_smooth <- smoothModel(ids@tpts, ratesFirstGuess(ids[input$select], 'total') - ratesFirstGuess(ids[input$select], 'preMRNA'))
			experiment$preMRNA_smooth <- smoothModel(ids@tpts, ratesFirstGuess(ids[input$select], 'preMRNA'))
			
			experiment$synthesissd <- sqrt(ratesFirstGuessVar(ids[input$select], 'synthesis'))
			experiment$mRNAsd <- sqrt(ratesFirstGuessVar(ids[input$select], 'total') + ratesFirstGuessVar(ids[input$select], 'preMRNA'))
			experiment$preMRNAsd <- sqrt(ratesFirstGuessVar(ids[input$select], 'preMRNA'))

			out <- define_parameter_ranges( ids, isolate(inspect$logshift) )

			gene_model <- ids@model@ratesSpecs[[input$select]][[1]]
			
			modeling$counts <- gene_model$counts[1]
			modeling$convergence <- gene_model$convergence

			function_types$k1 <- switch(
				gene_model[['alpha']][['type']],
				"constant" = "Constant",
				"sigmoid" = "Sigmoidal",
				"impulse" = "Impulsive"
			)

			function_types$k2 <- switch(
				gene_model[['gamma']][['type']],
				"constant" = "Constant",
				"sigmoid" = "Sigmoidal",
				"impulse" = "Impulsive"
			)

			function_types$k3 <- switch(
				gene_model[['beta']][['type']],
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
				values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],2)
				values$k1_h1 = round(gene_model[["alpha"]][["params"]][1],2)
				values$k1_h2 = round(gene_model[["alpha"]][["params"]][1],2)
				values$k1_t1 = round(mean(out$t_pars))
				values$k1_t2 = round(mean(out$t_pars))
				values$k1_beta = round(mean(out$beta_pars))
			}
			if( function_types$k1 == "Sigmoidal" ) {
				values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],2)
				values$k1_h1 = round(gene_model[["alpha"]][["params"]][2],2)
				values$k1_h2 = round(gene_model[["alpha"]][["params"]][2],2)
				values$k1_t1 = round(time_transf_inv(gene_model[["alpha"]][["params"]][3], inspect$logshift),2)
				values$k1_t2 = round(time_transf_inv(gene_model[["alpha"]][["params"]][3], inspect$logshift),2)
				values$k1_beta = round(gene_model[["alpha"]][["params"]][4],2)
			}
			if( function_types$k1 == "Impulsive" ) {
				values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],2)
				values$k1_h1 = round(gene_model[["alpha"]][["params"]][2],2)
				values$k1_h2 = round(gene_model[["alpha"]][["params"]][3],2)
				values$k1_t1 = round(time_transf_inv(gene_model[["alpha"]][["params"]][4], inspect$logshift),2)
				values$k1_t2 = round(time_transf_inv(gene_model[["alpha"]][["params"]][5], inspect$logshift),2)
				values$k1_beta = round(gene_model[["alpha"]][["params"]][6],2)
			}
			
			if( function_types$k2 == "Constant" ) {
				values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],2)
				values$k2_h1 = round(gene_model[["gamma"]][["params"]][1],2)
				values$k2_h2 = round(gene_model[["gamma"]][["params"]][1],2)
				values$k2_t1 = round(mean(out$t_pars))
				values$k2_t2 = round(mean(out$t_pars))
				values$k2_beta = round(mean(out$beta_pars))
			}
			if( function_types$k2 == "Sigmoidal" ) {
				values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],2)
				values$k2_h1 = round(gene_model[["gamma"]][["params"]][2],2)
				values$k2_h2 = round(gene_model[["gamma"]][["params"]][2],2)
				values$k2_t1 = round(time_transf_inv(gene_model[["gamma"]][["params"]][3], inspect$logshift),2)
				values$k2_t2 = round(time_transf_inv(gene_model[["gamma"]][["params"]][3], inspect$logshift),2)
				values$k2_beta = round(gene_model[["gamma"]][["params"]][4],2)
			}
			if( function_types$k2 == "Impulsive" ) {
				values$k2_h0 = round(gene_model[["gamma"]][["params"]][1],2)
				values$k2_h1 = round(gene_model[["gamma"]][["params"]][2],2)
				values$k2_h2 = round(gene_model[["gamma"]][["params"]][3],2)
				values$k2_t1 = round(time_transf_inv(gene_model[["gamma"]][["params"]][4], inspect$logshift),2)
				values$k2_t2 = round(time_transf_inv(gene_model[["gamma"]][["params"]][5], inspect$logshift),2)
				values$k2_beta = round(gene_model[["gamma"]][["params"]][6],2)
			}
			
			if( function_types$k3 == "Constant" ) {
				values$k3_h0 = round(gene_model[["beta"]][["params"]][1],2)
				values$k3_h1 = round(gene_model[["beta"]][["params"]][1],2)
				values$k3_h2 = round(gene_model[["beta"]][["params"]][1],2)
				values$k3_t1 = round(mean(out$t_pars))
				values$k3_t2 = round(mean(out$t_pars))
				values$k3_beta = round(mean(out$beta_pars))
			}
			if( function_types$k3 == "Sigmoidal" ) {
				values$k3_h0 = round(gene_model[["beta"]][["params"]][1],2)
				values$k3_h1 = round(gene_model[["beta"]][["params"]][2],2)
				values$k3_h2 = round(gene_model[["beta"]][["params"]][2],2)
				values$k3_t1 = round(time_transf_inv(gene_model[["beta"]][["params"]][3], inspect$logshift),2)
				values$k3_t2 = round(time_transf_inv(gene_model[["beta"]][["params"]][3], inspect$logshift),2)
				values$k3_beta = round(gene_model[["beta"]][["params"]][4],2)
			}
			if( function_types$k3 == "Impulsive" ) {
				values$k3_h0 = round(gene_model[["beta"]][["params"]][1],2)
				values$k3_h1 = round(gene_model[["beta"]][["params"]][2],2)
				values$k3_h2 = round(gene_model[["beta"]][["params"]][3],2)
				values$k3_t1 = round(time_transf_inv(gene_model[["beta"]][["params"]][4], inspect$logshift),2)
				values$k3_t2 = round(time_transf_inv(gene_model[["beta"]][["params"]][5], inspect$logshift),2)
				values$k3_beta = round(gene_model[["beta"]][["params"]][6],2)
			}

			# out <- define_parameter_ranges( ids, isolate(inspect$logshift) )

			gene_h_vals <- c(isolate(values$k1_h0),isolate(values$k1_h1),isolate(values$k1_h2),
				isolate(values$k2_h0),isolate(values$k2_h1),isolate(values$k2_h2),
				isolate(values$k3_h0),isolate(values$k3_h1),isolate(values$k3_h2))
			gene_t_vals <- c(isolate(values$k1_t1),isolate(values$k1_t2),isolate(values$k2_t1),
				isolate(values$k2_t2),isolate(values$k3_t1),isolate(values$k3_t2))
			gene_beta_vals <- c(isolate(values$k1_beta),isolate(values$k2_beta),isolate(values$k3_beta))
				
			ranges$k1_h0_min <- ranges$k1_h1_min <- ranges$k1_h2_min <- min(c(gene_h_vals,out$k1_h_pars[1]))
			ranges$k1_h0_max <- ranges$k1_h1_max <- ranges$k1_h2_max <- max(c(gene_h_vals,out$k1_h_pars[2]))
			ranges$k2_h0_min <- ranges$k2_h1_min <- ranges$k2_h2_min <- min(c(gene_h_vals,out$k2_h_pars[1]))
			ranges$k2_h0_max <- ranges$k2_h1_max <- ranges$k2_h2_max <- max(c(gene_h_vals,out$k2_h_pars[2]))
			ranges$k3_h0_min <- ranges$k3_h1_min <- ranges$k3_h2_min <- min(c(gene_h_vals,out$k3_h_pars[1]))
			ranges$k3_h0_max <- ranges$k3_h1_max <- ranges$k3_h2_max <- max(c(gene_h_vals,out$k3_h_pars[2]))
			ranges$t_min     <- min(c(gene_t_vals,out$t_pars[1]))
			ranges$t_max     <- max(c(gene_t_vals,out$t_pars[2]))
			ranges$beta_min  <- min(c(gene_beta_vals,out$beta_pars[1]))
			ranges$beta_max  <- max(c(gene_beta_vals,out$beta_pars[2]))
			
		}
		
	})
	
	## update convergence data

	observe({
		if( !is.null(modeling$counts) & !is.null(modeling$convergence) )
		output$convergence <- renderPrint({
			paste(modeling$counts,'iters (', switch(as.character(modeling$convergence),
				"0"="converged",
				"1"="not converged",
				"10"="degenerated"), ')')
			})
		})

	## modeling checkbox

	output$modeling_box <- renderUI({
		if( input$data_selection != 'User defined' ) {
			list(
				h5("pvalue of the chi-squared statistic:"),
				verbatimTextOutput("pchisq", TRUE),
				h5("Akaike information criterion:"),
				verbatimTextOutput("aic", TRUE),
				h5("minimization status:"),
				verbatimTextOutput("convergence", TRUE),
				fluidRow(
					column(4,h5('Refine model'), actionButton("optimize", "Optimize")),
					column(4,numericInput("nIter", label = h5("n iter"), value = 100)),
					column(4,radioButtons("opt_method", "method", 
						choices = c('NM','BFGS'), selected = 'NM'))
					)
				)
		}
		})


	######################################################################
	######################################################################
	### set the interactive part of the UI: ranges and values of the  
	### widgets are static at the beginning but can be changed upon the  
	### import of the INSPEcT dataset
	######################################################################
	######################################################################
	
	output$logtime_checkbox_ui <- renderUI({
		if( input$data_selection != 'User defined' ) {
			checkboxInput("logtime_checkbox", 
					label = "Space time logarithmically", 
					value = values$logtime)
		} else {
			NULL
			# checkboxInput("logtime_checkbox", 
			# 		label = "Space time logarithmically", 
			# 		value = FALSE)
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
				value = ranges$k1_h0_min, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$min_h_k1) )
			ranges$k1_h2_min <- ranges$k1_h1_min <- 
				ranges$k1_h0_min <- input$min_h_k1 
		})
	
	output$max_h_vals_k1 <- renderUI({
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) )
			numericInput("max_h_k1", label = h5("set max"), 
				value = ranges$k1_h0_max, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$max_h_k1) )
			ranges$k1_h2_max <- ranges$k1_h1_max <- 
				ranges$k1_h0_max <- input$max_h_k1 
		})
	
	output$ui_k1 <- renderUI({
		
		ids <- contentsrea()		
		if( !is.null(ids) & !is.null(input$k1_function) & !is.null(function_types$k1) )
			switch(input$k1_function,
				"Constant" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = ranges$k1_h0_min,
						max = ranges$k1_h0_max,
						value = values$k1_h0,
						step = 0.01)
				),
				"Sigmoidal" = list(
					sliderInput("k1_h0",
						"starting levels:",
						min = ranges$k1_h1_min,
						max = ranges$k1_h1_max,
						value = values$k1_h0,
						step = 0.01),
					sliderInput("k1_h1",
						"final levels:",
						min = ranges$k1_h1_min,
						max = ranges$k1_h1_max,
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
						min = ranges$k1_h0_min,
						max = ranges$k1_h0_max,
						value = values$k1_h0,
						step = 0.01),
					sliderInput("k1_h1",
						"intermediate levels:",
						min = ranges$k1_h1_min,
						max = ranges$k1_h1_max,
						value = values$k1_h1,
						step = 0.01),
					sliderInput("k1_h2",
						"end levels:",
						min = ranges$k1_h2_min,
						max = ranges$k1_h2_max,
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
				value = ranges$k2_h0_min, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$min_h_k2) )
			ranges$k2_h2_min <- ranges$k2_h1_min <- 
				ranges$k2_h0_min <- input$min_h_k2 
		})

	output$max_h_vals_k2 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("max_h_k2", label = h5("set max"), 
				value = ranges$k2_h0_max, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$max_h_k2) )
			ranges$k2_h2_max <- ranges$k2_h1_max <- 
				ranges$k2_h0_max <- input$max_h_k2 
		})
	
	output$ui_k2 <- renderUI({
		
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$k2_function) & !is.null(function_types$k2) )
			switch(input$k2_function,
				"Constant" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = ranges$k2_h0_min,
						max = ranges$k2_h0_max,
						value = values$k2_h0,
						step = 0.01)
				),
				"Sigmoidal" = list(
					sliderInput("k2_h0",
						"starting levels:",
						min = ranges$k2_h0_min,
						max = ranges$k2_h0_max,
						value = values$k2_h0,
						step = 0.01),
					sliderInput("k2_h1",
						"final levels:",
						min = ranges$k2_h1_min,
						max = ranges$k2_h1_max,
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
						min = ranges$k2_h0_min,
						max = ranges$k2_h0_max,
						value = values$k2_h0,
						step = 0.01),
					sliderInput("k2_h1",
						"intermediate levels:",
						min = ranges$k2_h1_min,
						max = ranges$k2_h1_max,
						value = values$k2_h1,
						step = 0.01),
					sliderInput("k2_h2",
						"end levels:",
						min = ranges$k2_h2_min,
						max = ranges$k2_h2_max,
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
				value = ranges$k3_h0_min, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$min_h_k3) )
			ranges$k3_h2_min <- ranges$k3_h1_min <- 
				ranges$k3_h0_min <- input$min_h_k3 
		})
	
	output$max_h_vals_k3 <- renderUI({
		ids <- contentsrea()
		if( !is.null(ids) )
			numericInput("max_h_k3", label = h5("set max"), 
				value = ranges$k3_h0_max, width='200px')
	})

	observe({
		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$max_h_k3) )
			ranges$k3_h2_max <- ranges$k3_h1_max <- 
				ranges$k3_h0_max <- input$max_h_k3 
		})

	output$ui_k3 <- renderUI({

		ids <- contentsrea()
		if( !is.null(ids) & !is.null(input$k3_function) & !is.null(function_types$k3) )
			switch(input$k3_function,
				"Constant" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = ranges$k3_h0_min,
						max = ranges$k3_h0_max,
						value = values$k3_h0,
						step = 0.01)
				),
				"Sigmoidal" = list(
					sliderInput("k3_h0",
						"starting levels:",
						min = ranges$k3_h0_min,
						max = ranges$k3_h0_max,
						value = values$k3_h0,
						step = 0.01),
					sliderInput("k3_h1",
						"final levels:",
						min = ranges$k3_h1_min,
						max = ranges$k3_h1_max,
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
						min = ranges$k3_h0_min,
						max = ranges$k3_h0_max,
						value = values$k3_h0,
						step = 0.01),
					sliderInput("k3_h1",
						"intermediate levels:",
						min = ranges$k3_h1_min,
						max = ranges$k3_h1_max,
						value = values$k3_h1,
						step = 0.01),
					sliderInput("k3_h2",
						"end levels:",
						min = ranges$k3_h2_min,
						max = ranges$k3_h2_max,
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
	
	output$gene <- renderPlot({
		
		ids <- contentsrea()
		if( input$select != "" & !is.null(ids) & !is.null(input$k1_function) &
			!is.null(input$k2_function) & !is.null(input$k3_function))

			try({

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
					scores <- RNAdynamicsAppPlot(
						data_selection = input$data_selection,
						show_logtime = input$logtime_checkbox,
						logshift = inspect$logshift,
						time_min = ranges$time_min,
						time_max = ranges$time_max,
						experiment = experiment,
						k1_function = input$k1_function, 
						k2_function = input$k2_function, 
						k3_function = input$k3_function,
						k1_params = k1_params,
						k2_params = k2_params,
						k3_params = k3_params
					)
					
					output$pchisq <- renderPrint({signif(scores$pchisq,3)})
					output$aic <- renderPrint({signif(scores$aic,3)})
				}

			}, silent = FALSE)
		
	})

	#######################
	## OPTIMIZE ###########
	#######################

	observeEvent(input$optimize, {

		# print('optimization started')

		log_shift <- inspect$logshift
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
			"Constant" = list(
				type='constant',
				fun=newPointer(constantModelRNApp),
				params=input$k1_h0,
				df=1
				),
			"Sigmoidal" = list(
				type='sigmoid',
				fun=newPointer(sigmoidModelRNApp),
				params=c(input$k1_h0, input$k1_h1, input$k1_t1, input$k1_beta),
				df=4
				),
			"Impulsive" = list(
				type='constant',
				fun=newPointer(impulseModelRNApp),
				params=c(input$k1_h0, input$k1_h1, 
						input$k1_h2, input$k1_t1, input$k1_t2, input$k1_beta),
				df=6
				)
			)

		k2_rate <- switch(input$k2_function, 
			"Constant" = list(
				type='constant',
				fun=newPointer(constantModelRNApp),
				params=input$k2_h0,
				df=1
				),
			"Sigmoidal" = list(
				type='sigmoid',
				fun=newPointer(sigmoidModelRNApp),
				params=c(input$k2_h0, input$k2_h1, input$k2_t1, input$k2_beta),
				df=4
				),
			"Impulsive" = list(
				type='constant',
				fun=newPointer(impulseModelRNApp),
				params=c(input$k2_h0, input$k2_h1, 
						input$k2_h2, input$k2_t1, input$k2_t2, input$k2_beta),
				df=6
				)
			)

		k3_rate <- switch(input$k3_function, 
			"Constant" = list(
				type='constant',
				fun=newPointer(constantModelRNApp),
				params=input$k3_h0,
				df=1
				),
			"Sigmoidal" = list(
				type='sigmoid',
				fun=newPointer(sigmoidModelRNApp),
				params=c(input$k3_h0, input$k3_h1, input$k3_t1, input$k3_beta),
				df=4
				),
			"Impulsive" = list(
				type='constant',
				fun=newPointer(impulseModelRNApp),
				params=c(input$k3_h0, input$k3_h1, 
						input$k3_h2, input$k3_t1, input$k3_t2, input$k3_beta),
				df=6
				)
			)

		params <- list(alpha=k1_rate, beta=k3_rate, gamma=k2_rate)
		gene_model <- optimParamsMatureRNA(params, tpts_exp, alpha_exp, alpha_var, mature_exp
		 	, mature_var, preMRNA_exp, preMRNA_var, maxit=input$nIter
		 	, method=input$opt_method, log_shift)

		modeling$counts <- modeling$counts + gene_model$counts[1]
		modeling$convergence <- gene_model$convergence

		#########################################
		### update the parameters in the GUI #######
		#########################################

		if( input$k1_function == "Constant" ) {
			values$k1_h0 = round(gene_model[["alpha"]][["params"]][1],2)
			values$k1_h1 = round(gene_model[["alpha"]][["params"]][1],2)
			values$k1_h2 = round(gene_model[["alpha"]][["params"]][1],2)
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

		# print('optimization finished')

		})
	
})