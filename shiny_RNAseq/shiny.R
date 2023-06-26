library(shiny)
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(pheatmap))

# Increase maximum file upload size to 50MB
options(shiny.maxRequestSize = 50*1024^2)


ui <- fluidPage(
  titlePanel("Upload CSV Files and Perform DESeq2 Analysis"),
  verticalLayout(
    fileInput("countFile", "Upload Expected Counts File", accept = ".csv"),
    fileInput("colFile", "Upload Coldata File", accept = ".csv"),
    numericInput("feature", "Feature to compare", value = 3),
    textInput("filename", "Output filename", value = "Tcell.csv"),
    actionButton("analyze", "Analyze"),
    verbatimTextOutput("countsStr"),
    verbatimTextOutput("countsHead"),
    verbatimTextOutput("coldataStr"),
    verbatimTextOutput("coldataHead"),
    tableOutput("resultTable"),
    downloadButton('downloadData', 'Download')
  )
)

server <- function(input, output) {
  results <- reactiveVal()
  
  output$countsStr <- renderPrint({
    req(input$countFile)
    counts <- read.csv(input$countFile$datapath)
    print(str(counts))
  })
  
  output$countsHead <- renderPrint({
    req(input$countFile)
    counts <- read.csv(input$countFile$datapath, stringsAsFactors = TRUE)
    print(counts[1:10,1:6])
  })
  
  output$coldataStr <- renderPrint({
    req(input$colFile)
    coldata <- read.csv(input$colFile$datapath, stringsAsFactors = TRUE)
    coldata <- coldata[,-1]
    print(str(coldata))
  })
  
  output$coldataHead <- renderPrint({
    req(input$colFile)
    coldata <- read.csv(input$colFile$datapath, stringsAsFactors = TRUE)
    print(head(coldata))
  })
  
  observeEvent(input$analyze, {
    req(input$countFile, input$colFile)
    
    counts <- read.csv(input$countFile$datapath)
    coldata <- read.csv(input$colFile$datapath, stringsAsFactors = TRUE)
    
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    countmat <- round(apply(counts, 2, as.numeric))
    #countmat = round(as.numeric((counts[-1,])))
    rownames(countmat) = rownames(counts)
    coldata <- coldata[,-1]
    names(coldata) <- c("BioSample", "feature1", "feature2")
    
    
    dds = DESeqDataSetFromMatrix(countData = countmat,
                                 colData   = coldata,
                                 design    = ~ feature1 + feature2 + feature1:feature2)
    
    keep <- which(rowSums(counts(dds)) > 2)
    dds <- dds[keep,]
    
    dds <- DESeq(dds)
    
    res_name <- resultsNames(dds)
    result <- lfcShrink(dds, coef = res_name[input$feature], type = "apeglm")
    result <- result[order(result$padj),]
    
    results(result)
  })
  
  output$resultTable <- renderTable({
    req(results())
    head(results())
  }, rownames = TRUE)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      input$filename
    },
    content = function(file) {
      write.csv(results(), file, row.names = T)
    }
  )
}

shinyApp(ui = ui, server = server)
