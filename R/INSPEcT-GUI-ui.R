# library(shiny)
# options(shiny.maxRequestSize=100*1024^2) 
# 

INSPEcTGUIshinyAppUI <- fluidPage(
	
	# version number
	
	tags$head(
		tags$style(HTML("hr {border-top: 1px solid #000000;}"))
	),	
	h4('INSPEcT-GUI'),
	fluidRow(
		column(2,
					 
					 fileInput("file1", "Choose INSPEcT File", accept = ".rds", width = NULL,
					 					buttonLabel = "Browse...", placeholder = "INSPEcT_GUI_sample_dataset.rds"),
					 uiOutput("file_error"),
					 uiOutput('modeling_type'),
					 uiOutput('select_class'),
					 uiOutput('select_condition'),
					 selectInput("select", label = "Select gene", choices = NULL, selected = NULL),
					 hr(),
					 radioButtons('data_selection', 'Explore RNA dynamics using:',
					 						 choiceValues = c('Experimental data', 'Smooth data', 'User defined'), 
					 						 choiceNames = c('Raw experimental data', 'Smooth experimental data', 'User defined (No input)'), 
					 						 selected = 'Experimental data'),
					 hr(),
					 uiOutput('modeling_box')
		),
		column(4,
					 plotOutput("gene", height = "600px"),#, width="500px"),
					 fluidRow(
					 	column(2, .busyIndicator(text="Loading..." , wait=0 , image='gif.gif')),
					 	column(4, uiOutput('confint_box')),
					 	column(3, list(
					 		checkboxInput("relativexpr_checkbox", label = "Rel. Expr.", value = FALSE),
					 		uiOutput("logtime_checkbox_ui"),
					 		checkboxInput("fixyaxis_checkbox", label = "Fix Y-axis", value = FALSE)
					 	)),
					 	column(3, list(
					 		downloadButton('saveRNAdynamicsPlotButton', 'Get PDF'),
					 		downloadButton('saveRNAdynamicsDataButton', 'Get TSV')
					 		))
					 )
		),
		
		# input panel
		
		# column 1: k1 - synthesis rate
		
		column(2,
					 uiOutput("fun1_name"), # h3("synthesis"), 
					 uiOutput("fun1_unit"), # h4("(RPKMs/hour)"), 
					 uiOutput("function_type_k1"),
					 fluidRow(
					 	column(5,uiOutput("min_h_vals_k1")),
					 	column(1),
					 	column(5,uiOutput("max_h_vals_k1"))
					 ),
					 uiOutput("ui_k1")
		),
		
		# column 2: k2 - processing rate
		
		column(2,
					 h3("processing"), 
					 h4("(1/hour)"), 
					 uiOutput("function_type_k2"),
					 fluidRow(
					 	column(5,uiOutput("min_h_vals_k2")),
					 	column(1),
					 	column(5,uiOutput("max_h_vals_k2"))
					 ),
					 uiOutput("ui_k2")
		),
		
		# column 3: k3 - degradation rate
		
		column(2,
					 h3("degradation"), 
					 h4("(1/hour)"), 
					 uiOutput("function_type_k3"),
					 fluidRow(
					 	column(5,uiOutput("min_h_vals_k3")),
					 	column(1),
					 	column(5,uiOutput("max_h_vals_k3"))
					 ),
					 uiOutput("ui_k3")
		)
		
	)
	
)