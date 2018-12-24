library(shiny)
options(shiny.maxRequestSize=100*1024^2) 

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
    ,div(class = "mybusyindicator",p(text),img(src=image))
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
    },1000)
    ",wait)
    )
  ) 
}

shinyUI(fluidPage(
	
	# version number
	
	h4('Inspect your INSPEcT data'),
	
	### h3("Select rates from file"),
		
	# plot pre-mRNA and mRNA dynamics

	fluidRow(
		column(4,

			fileInput("file1", "Choose INSPEcT File", accept = "Rfiles/rds", width = NULL,
				buttonLabel = "Browse...", placeholder = "mycerIds_20genes.rds"),
			selectInput("select_class", label = "Select class", 
				choices = NULL, selected = NULL),
			selectInput("select", label = "Select gene", 
				choices = NULL, selected = NULL),
			radioButtons('data_selection', 'Select reference:',
				choices = c('User defined','Experimental data','Smooth data'), 
				selected = 'Smooth data'),
			uiOutput('modeling_box')
			),
		column(8, 
			plotOutput("gene", height = "600px"),#, width="500px"),
			fluidRow(
				column(8, .busyIndicator(text="Loading..." , wait=1000 , image="gif.gif")), 
				column(4, uiOutput("logtime_checkbox_ui")))
			)
		
	),
	
	# input panel
	
	fluidRow(
		
		# column 1: k1 - synthesis rate
	
		column(4,
			h3("synthesis"), 
			h4("(RPKMs/hour)"), 
			uiOutput("function_type_k1"),
			fluidRow(
				column(5,uiOutput("min_h_vals_k1")),
				column(1),
				column(5,uiOutput("max_h_vals_k1"))
				),
			uiOutput("ui_k1")
		),
		
		# column 2: k2 - processing rate
		
		column(4,
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
		
		column(4,
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
)