library(shiny)
library(gglogo)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(seqinr)

data(aacids)

function(input, output, session) {
    
    values <- reactiveValues(seqs = "")
    
    observeEvent(input$confirm, {
        values$seqs <<- input$sequence
    })
    
    ## The initial uploaded dataset
    seqs.initial <- reactive({
        if (is.null(input$data)) return(NULL)
        else return(read.fasta(input$data$datapath))
    })
    
    output$plotbuilt <- reactive({
        return(nrow(my.seqdata()) > 0)
    })
    outputOptions(output, 'plotbuilt', suspendWhenHidden = FALSE)

    my.df <- reactive({
        if (nchar(values$seqs) == 0 && is.null(seqs.initial())) return(NULL)
        
        test <- seqs.initial()
        mydf <- ldply(test, function(my.seq) {
            paste(toupper(as.character(my.seq)), collapse = "")
        })
        names(mydf) <- c("factor", "peptide")
        
        return(mydf)
    })
    
    observe({
        if (!is.null(my.df())) {
            my.df <- my.df()
            updateSliderInput(session, "zoom", max = nchar(as.character(my.df$peptide[1])), value = c(1, min(30, nchar(as.character(my.df$peptide[1])))))
        }
    })
    
    my.seqdata <- reactive({
        if (is.null(my.df())) return(NULL)
        my.df <- my.df()
        my.df$peptide <- sapply(strsplit(my.df$peptide, ""), function(x){paste(x[input$zoom[1]:input$zoom[2]], collapse = "")})
        
        return(my.df)
    })
    
    myplot <- reactive({
        if (is.null(my.seqdata())) return(NULL)
        
        withProgress(message = "Building logo plot, please wait...", expr = {
            test <- my.seqdata()
            
            dm2 <- splitSequence(test, "peptide")
            cols <- c(input$col1, input$col2, input$col3, input$col4)
            
            my.trt <- if (input$facetvar == "Factor") "factor" else NULL
            
            dm3 <- calcInformation(dm2, trt = my.trt, pos="position", elems="element", k=21)
            dm3b <- merge(dm3, aacids, by.x = "element", by.y = "AA", all.x = TRUE)

            dm3b$position <- as.numeric(as.character(dm3b$position))
            dm3b$position <- dm3b$position + input$zoom[1] - 1
            dm3b$position <- factor(dm3b$position, levels = sort(unique(dm3b$position)))
            
            #dm3b$facet_group <- cut(as.numeric(dm3b$position), seq(0, max(as.numeric(dm3b$position)) + 29, by = 30), labels = FALSE)
            
            dm3b$x_var <- if (input$facetvar == "Factor") dm3b$factor else dm3b$position
            my_text <- element_text(angle = ifelse(input$facetvar == "Factor", 90, 0), vjust = ifelse(input$facetvar == "Factor", 0.5, 1))
            
            dm3b$Color <- dm3b[,input$colorvar]
            dm3bb <- subset(dm3b, element != "-")
            
            ggplot(dm3b, aes(x = x_var, y = bits, group = element, label = element, fill = Color), alpha = 0.8) + 
                geom_hline(yintercept=-log(1/21, base=2), colour="grey30", size=0.5) + 
                geom_logo() + 
                scale_fill_discrete(na.value = "white") +
                geom_hline(yintercept=0, colour="white", size=0.5) + 
                geom_hline(yintercept=0, colour="grey30", size=0.125) + 
                theme_bw() + theme(legend.position="top", plot.margin=unit(c(0,0,0,0), "cm"), axis.text.x = my_text) + 
                xlab(input$xlab) +
                ylab(input$ylab) + 
                ggtitle(input$title) +
                scale_y_continuous(breaks=c(-1,0,1,2,4), labels=c(1,0,1,2,4)) +
                if (input$facetvar == "Factor") facet_wrap(~position)
        })
    })
    
    output$logoplot <- renderPlot({
        print(myplot())
    })
    
    output$download <- downloadHandler(
        filename = function() { paste("logoplot", input$image_format, sep = ".") },
        content = function(file) {
            ggsave(file, plot = myplot(), dpi = input$dpi, height = 5, width = 10)
        }
    )
    
}
