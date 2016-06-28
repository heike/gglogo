library(shiny)
library(shinythemes)
library(seqinr)

fluidPage(theme = shinytheme("cerulean"),

    titlePanel("GGLogo Prototype"),
    
    sidebarLayout(
        sidebarPanel(width = 3,
            h4("Sequence"),
            
            selectizeInput("mychoice", label = "Choose Input Method", choices = c("Upload Seq Data" = "upload", "Type Seq Data" = "type")),
            
            conditionalPanel(condition = "input.mychoice == 'type'", 
                helpText("Input sequencing data"),
                tags$textarea(id="sequence", rows=3, cols=40, "")
            ),
            
            conditionalPanel(condition = "input.mychoice == 'upload'", 
                fileInput("data", "Upload Data (FASTA)")
            ),
            
            #actionButton("confirm", "Build Logo Plot"),
            
            conditionalPanel(condition = "output.plotbuilt == true",
                downloadButton("download", "Download Logo Plot")
            ),
            
            hr(),
            
            checkboxInput("advanced", "Show Advanced Configuration"),
            
            conditionalPanel(condition = "input.advanced == true",
                h4("Advanced Configuration"),
                
                h5("Plot Options"),
                selectizeInput("facetvar", label = "Facet Variable", choices = c("None", "Factor")),
                selectizeInput("colorvar", label = "Color By", choices = c("Polarity", "Water")),
                         
                hr(),
                       
                h5("Plot Labels"),
                textInput("title", "Plot Title", value = ""),
                textInput("xlab", "X Axis Label", value = "position"),
                textInput("ylab", "Y Axis Label", value = "bits"),
                
                hr(),
                
                h5("Download Options"),
                selectizeInput("image_format", "Image Format", choices = c("PNG" = "png", "JPG" = "jpg", "PDF" = "pdf")),
                numericInput("dpi", "DPI", value = 300, min = 1, max = 1000)
            )
        ),
        
        mainPanel(width = 9,
            sliderInput("zoom", "Sequence Region", min = 1, max = 231, value = c(1, 30), width = "100%", dragRange = TRUE, animate = TRUE),
            plotOutput("logoplot", height = "500px")
        )
    )
)
