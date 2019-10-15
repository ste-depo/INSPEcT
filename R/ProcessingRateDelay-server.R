ProcessingRateDelayshinyAppServer <- function(input, output, session) {

	metrics_gene1 <- reactiveValues()
	values  <- reactiveValues()
	inspect <- reactiveValues()
	slider_ranges <- reactiveValues()
	filter_selected <- reactiveValues()

	# load INSPEcT file
	
	readInspectObject <- reactive({

		filename <- input$file1$datapath 
		if( is.null(filename) ) filename <- system.file(package='INSPEcT', 'nascentInspObj.rds')

		## load file
		if( file.exists(filename) ) {
			ids <- readRDS(filename)

			if( class(ids) != 'INSPEcT' ) {
			
				return(NULL)

			} else{ # (the loaded object is of class INSPEcT

				# update rates
				rates <- data.frame(
					k1 = ratesFirstGuess(ids, 'synthesis')[,1],
					k2 = ratesFirstGuess(ids, 'processing')[,1],
					k3 = ratesFirstGuess(ids, 'degradation')[,1]
					)

				# select only genes with finite rates
				filter <- is.finite(rates$k1) & is.finite(rates$k2) &
					is.finite(rates$k3)
				rates <- rates[filter,]
				ids <- ids[filter]

				# update gene names in the input box
				updateSelectInput(session, "select", 
					choices = sort(featureNames(ids)), selected = NULL)

				# update metrics
				metrics <- data.frame(
					m1 = sapply(1:nrow(ids), function(i) 
						tau_fun( rates$k2[i] , rates$k3[i] )),
					m2 = sapply(1:nrow(ids), function(i) 
						delta_fun( rates$k1[i], rates$k2[i] , rates$k3[i] ))
					)
				rownames(metrics) <- rownames(rates)

				# define ranges for the plot of rates distribution
				allrates_range <- quantile(unlist(rates), 
					probs=c(0.025,.975), na.rm=TRUE)

				# by default assign all genes
				filter_selected$ix <- featureNames(ids)

				return(list(ids=ids, rates=rates, 
					metrics=metrics , allrates_range=allrates_range))

			}

		} else { # (if the file name does not exist)

			return(NULL)

		}
	})

	############################################################
	##### update slider_ranges to include the rates of the #######
	##### selected gene within the ranges #####################
	########################################

	observe({

		data <- readInspectObject()

		if( !is.null(data) & input$select != "" & 
			input$select %in% featureNames(data$ids)) {

			ids <- data$ids
			values$k1 <- ratesFirstGuess(ids, 'synthesis')[input$select,1]
			values$k2 <- ratesFirstGuess(ids, 'processing')[input$select,1]
			values$k3 <- ratesFirstGuess(ids, 'degradation')[input$select,1]

			# define ranges 
			rates <- data$rates
			range_k1 <- quantile(rates$k1, probs=c(0.025,0.975), na.rm=TRUE)
			range_k2 <- quantile(rates$k2, probs=c(0.025,0.975), na.rm=TRUE)
			range_k3 <- quantile(rates$k3, probs=c(0.025,0.975), na.rm=TRUE)

			# include the selected gene within the range
			range_k1 <- range(range_k1, isolate(values$k1))
			range_k2 <- range(range_k2, isolate(values$k2))
			range_k3 <- range(range_k3, isolate(values$k3))

			# update the reactive value
			slider_ranges$k1_min <- floor(range_k1[1]*100)/100
			slider_ranges$k1_max <- ceiling(range_k1[2])
			slider_ranges$k2_min <- floor(range_k2[1]*100)/100
			slider_ranges$k2_max <- ceiling(range_k2[2])
			slider_ranges$k3_min <- floor(range_k3[1]*100)/100
			slider_ranges$k3_max <- ceiling(range_k3[2])

		}

	})

	####################################################
	##### plot sliders and max levels for each rate #####
	###################################################

	output$ui_gene1 <- renderUI({

		data <- readInspectObject()
		if( !is.null(data) & input$select != "" & 
			input$select %in% featureNames(data$ids)) {
			metrics <- data$metrics
			fluidRow(
				column(8,
					sliderInput("k1_gene1",
								"synthesis (RPKMs/hour):",
								min = slider_ranges$k1_min,
								max = slider_ranges$k1_max,
								value = values$k1,
								step = 0.01),
					sliderInput("k2_gene1",
								"processing (1/hour):",
								min = slider_ranges$k2_min,
								max = slider_ranges$k2_max,
								value = values$k2,
								step = 0.01),
					sliderInput("k3_gene1",
								"degradation (1/hour):",
								min = slider_ranges$k3_min,
								max = slider_ranges$k3_max,
								value = values$k3,
								step = 0.01)
				),
				column(4,
					numericInput("max_k1", label = h5("max levels"), 
						value = slider_ranges$k1_max, width='75px'),
					numericInput("max_k2", label = h5("max levels"), 
						value = slider_ranges$k2_max, width='75px'),
					numericInput("max_k3", label = h5("max levels"), 
						value = slider_ranges$k3_max, width='75px')
				)
			)

		}

	})

	#######################################
	##### observe changes in max levels #####
	#######################################

	observe({

		slider_ranges$k1_max <- input$max_k1
		slider_ranges$k2_max <- input$max_k2
		slider_ranges$k3_max <- input$max_k3

	})

	########################################
	##### observe click on metrics plot ########
	#######################################

	observe({
		data <- readInspectObject()
		if( !is.null(data) ) {
			metrics <- data$metrics
			ids     <- data$ids
  			new_selected <- rownames(nearPoints(metrics, input$metrics_click, 
  				xvar='m2', yvar='m1', threshold = 10, maxpoints = 1))
  			if( length(new_selected)>0 )
				updateSelectInput(session, "select", 
					choices = sort(featureNames(ids)), selected = new_selected)
  		}
  		})

	########################################
	##### observe brush on metrics plot ########
	#######################################

	observe({
		data <- readInspectObject()
		metrics <- data$metrics
		filter_selected_ix <- brushedPoints(metrics, input$metrics_brush, 
			xvar='m2', yvar='m1', allRows = TRUE)$selected
		if( any(filter_selected_ix) ) {
			filter_selected$ix <- filter_selected_ix
			updateSelectInput(session, "select", 
				choices = sort(featureNames(data$ids)[filter_selected_ix]), 
				selected = NULL)
			session$resetBrush("metrics_brush")
		}
	})

	##############################################
	##### observe double click on metrics plot ########
	#############################################

	observe({
		if( !is.null(input$metrics_dblclick) ) {
			data <- readInspectObject()
			filter_selected$ix <- featureNames(data$ids)
			updateSelectInput(session, "select", 
				choices = sort(featureNames(data$ids)), 
				selected = NULL)
		}
	})

  	############################################
  	## plot mature mRNA of the selected gene #######
  	############################################
	
	output$gene1 <- renderPlot({

		if( length(input$k1_gene1)>0 ) {

			out <- shinyProcessingDelayPlot(
				k1=input$k1_gene1,
				k2=input$k2_gene1,
				k3=input$k3_gene1,
				absval=input$absval_gene1, 
				metrics=input$metrics_gene1
			)
			
			metrics_gene1$tau_value <- out$tau_value
			metrics_gene1$delta_value <- out$delta_value

		}
		
	})

  	#####################################
  	## plot mature rates distribution #######
  	#####################################

	output$rates_distribution <- renderPlot({

		data <- readInspectObject()
		if( !is.null(data) & length(input$k1_gene1)>0 ) {

			all_rates <- data$rates
			selected_rates <- all_rates[filter_selected$ix,]
			allrates_range <- data$allrates_range

			par(mar=c(5,4,2,2))

			x_log <- seq(
				log10(allrates_range[1]),
				log10(allrates_range[2]),
				length.out = 512
				)
			y <- cbind(
				density( log10(all_rates$k1) , from = x_log[1], to = x_log[512], n = 512 )$y,
				density( log10(all_rates$k2) , from = x_log[1], to = x_log[512], n = 512 )$y,
				density( log10(all_rates$k3) , from = x_log[1], to = x_log[512], n = 512 )$y,
				density( log10(selected_rates$k1) , from = x_log[1], to = x_log[512], n = 512 )$y,
				density( log10(selected_rates$k2) , from = x_log[1], to = x_log[512], n = 512 )$y,
				density( log10(selected_rates$k3) , from = x_log[1], to = x_log[512], n = 512 )$y
				)
			matplot(10^x_log, y, type='l', lty=c(1,1,1,2,2,2), lwd=3, log='x', col=c(1:3,1:3),
				xlab='rate value', ylab='probability density')
			legend('topright', lty=1, lwd=3, col=1:3, bty='n',
				legend=c('synthesis','processing','degradation'))

			k1_gene1_idx <- which.min(abs(x_log - log10(input$k1_gene1)))
			k2_gene1_idx <- which.min(abs(x_log - log10(input$k2_gene1)))
			k3_gene1_idx <- which.min(abs(x_log - log10(input$k3_gene1)))

			points(10^(x_log[k1_gene1_idx]), 0, pch=19, col=1, cex=1.3)
			points(10^(x_log[k2_gene1_idx]), 0, pch=19, col=2, cex=1.3)
			points(10^(x_log[k3_gene1_idx]), 0, pch=19, col=3, cex=1.3)

			title( expression( 
				'distribution of all (solid) and selected (dashed) genes' ))

		}
	})

  	############################
  	## plot metrics dotplot #######
  	############################

	output$metrics_dotplot <- renderPlot({

		data <- readInspectObject()
		if( !is.null(data) & length(input$k1_gene1)>0 ) {

			metrics <- data$metrics[filter_selected$ix,]

			par(mar=c(5,4,2,2))

			if( input$density_color ) {

				denspalette <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
					"#FCFF00", "#FF9400", "#FF3100"))
				cols <- with(metrics, 
					densCols(log(m1), log(m2), nbin=512, colramp=denspalette))

			} else {

				cols <- rep('black', length(metrics$m1))
				cols[metrics$m1>1.5 & metrics$m2>1] <- 'darkgoldenrod1'
				cols[metrics$m1>1.5 & metrics$m2<1] <- 'blue'
				cols[metrics$m1<1.5 & metrics$m2>1] <- 'red'

			}

			plot(metrics$m2, metrics$m1, log='xy', pch=20, col = cols , 
				xlim = range(c(metrics$m2, metrics_gene1$delta_value)), 
				ylim = range(c(metrics$m1, metrics_gene1$tau_value)),
				xlab = expression(tau), ylab = expression(Delta),
				cex.lab = 1.5
				)
			abline(h=1.5, col='black', lty=3, lwd=1)
			abline(v=1, col='black', lty=3, lwd=1)

			points( metrics_gene1$delta_value , 
				metrics_gene1$tau_value , cex=1.5, lwd=5, col='gray50')

			if( input$density_color ) {

				title(
					expression(paste('processing affecting ', tau, 
						', ', Delta, ' or ', 'both')), 
					line=1)

			} else {

				## build the title with different colors

				title(
					expression(paste('processing affecting ', phantom(tau), 
						', ', phantom(Delta), ' or ', phantom('both'))), 
					line=1)
				title(
					expression(paste(phantom('processing affecting '), tau, 
						phantom(', '), phantom(Delta), phantom(' or '), phantom('both'))), 
					col.main='blue', line=1)
				title(
					expression(paste(phantom('processing affecting '), phantom(tau), 
						phantom(', '), Delta, phantom(' or '), phantom('both'))), 
					col.main='red', line=1)
				title(
					expression(paste(phantom('processing affecting '), phantom(tau), 
						phantom(', '), phantom(Delta), phantom(' or '), 'both')), 
					col.main='darkgoldenrod1', line=1)

			}

		}

	})

	####################
	## show help #########
	####################

	output$help <- renderText({
		# if( !is.null(input$gene_hover) ) return('gene plot')
		# if( !is.null(input$distributions_hover) ) return('distributions plot')
		if( !is.null(input$metrics_hover) ) return(
			'Single click to pick a gene<br/>
			Drag to select a set of genes<br/>
			Double click to visualize all genes again')
		# return('Select a gene from the loaded INSPEcT object to see the corresponding metrics and how it is positioned compared the population')
		return('')
		})

}