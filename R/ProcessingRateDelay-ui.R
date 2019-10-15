# library(shiny)
# options(shiny.maxRequestSize=100*1024^2) 

# Define UI for application that draws a histogram
ProcessingRateDelayshinyAppUI <- fluidPage(
	
	# Application title
	# titlePanel("Hello Shiny! v4"),
	
	fluidRow(
		column(6,
			fileInput("file1", "Choose INSPEcT File", accept = "Rfiles/rds", 
				width = NULL, buttonLabel = "Browse...", 
				placeholder = "nascentInspObj.rds"),
			selectInput("select", label = "Select or type gene name", 
				choices = NULL, selected = NULL),
			uiOutput('ui_gene1')
		),
		column(6,
			plotOutput("gene1", width="400px", hover = 'gene_hover'),
			checkboxInput("metrics_gene1", "Show metrics", TRUE),
			checkboxInput("absval_gene1", "Show absolute expression", TRUE)
		)
	),
	hr(),
	fluidRow(
		column(6,plotOutput("rates_distribution", width="400px",
			hover = 'distributions_hover')),
		column(6,plotOutput("metrics_dotplot", width="400px", 
			click = "metrics_click", brush = "metrics_brush",
			dblclick = "metrics_dblclick", hover = 'metrics_hover'))
		),
	fluidRow(
		column(1),
		column(6, .busyIndicator(text="Loading..." , wait=1000)), # , image="gif.gif")),
		column(5, 
			checkboxInput("density_color", "Color dots by density", FALSE),
			htmlOutput("help")
			)
		)
)