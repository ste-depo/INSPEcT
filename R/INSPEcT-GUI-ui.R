# library(shiny)
# options(shiny.maxRequestSize=100*1024^2) 
# 

INSPEcTGUIshinyAppUI <- fluidPage(
	
	# version number
	
	h4('Inspect your INSPEcT data'),
	
	### h3("Select rates from file"),
	
	# plot pre-mRNA and mRNA dynamics
	
	fluidRow(
		column(2,
					 
					 fileInput("file1", "Choose INSPEcT File", accept = ".rds", width = NULL,
					 					buttonLabel = "Browse...", placeholder = "nascentInspObj10.rds"),
					 uiOutput('modeling_type'),
					 uiOutput('select_condition'),
					 uiOutput('select_class'),
					 selectInput("select", label = "Select gene", 
					 						choices = NULL, selected = NULL),
					 radioButtons('data_selection', 'Select input:',
					 						 choiceValues = c('Smooth data', 'Experimental data','User defined'), 
					 						 choiceNames = c('Smooth data', 'Experimental data','User defined (No input)'), 
					 						 selected = 'Smooth data'),
					 uiOutput('modeling_box')
		),
		column(4, 
					 fluidRow(column(2,NULL), column(8,uiOutput("rate_pvals")), column(2,NULL)),
					 plotOutput("gene", height = "600px"),#, width="500px"),
					 fluidRow(
					 	column(4, .busyIndicator(text="Loading..." , wait=1000 , image='gif.gif')), 
					 	# column(4, uiOutput("saveRNAdynamicsPlot")),
					 	column(4, list(
					 		downloadButton('saveRNAdynamicsPlotButton', 'Get PDF'),
					 		downloadButton('saveRNAdynamicsDataButton', 'Get XLS')
					 		)),
					 	column(4, uiOutput("logtime_checkbox_ui"), uiOutput("confint_checkbox_ui"))
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