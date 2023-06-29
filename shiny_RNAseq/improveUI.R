library(shiny)
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(pheatmap))


# Increase maximum file upload size to 50MB
options(shiny.maxRequestSize = 50*1024^2)


ui <- fluidPage(
  includeCSS("../css/first.css"),
  titlePanel("Upload CSV Files and Perform DESeq2 Analysis"),
  tabsetPanel(
    tabPanel("Data Upload",
             wellPanel(
               fileInput("countFile", "Upload Expected Counts File", accept = ".csv"),
               verbatimTextOutput("countsStr"),
               verbatimTextOutput("countsHead"),
               fileInput("colFile", "Upload Coldata File", accept = ".csv"),
               verbatimTextOutput("coldataStr"),
               verbatimTextOutput("coldataHead"),
               helpText("Upload your CSV files here for further analysis.")
             )
    ),
    tabPanel("DESeq2 Analysis",
             wellPanel(
               numericInput("feature", "Feature to compare", value = 3),
               helpText("Enter the feature you want to compare."),
               actionButton("analyze", "Analyze"),
               tableOutput("resultTable"),
               textInput("filename", "Output filename", value = "Tcell.csv"),
               downloadButton('downloadData', 'Download Full Gene'),
               numericInput("lfc_threshold", "Log2 fold change threshold", value = 1.5),
               numericInput("padj_threshold", "Adjusted p-value threshold", value = 0.05),
               actionButton("filter", "Filter The genes"),
               tableOutput("filteredTable"),
               downloadButton('downloadFiltered', 'Download Filtered Results'),
             )
    ),
    tabPanel("Plots",
             sidebarLayout(
               sidebarPanel(
                 actionButton("plot", "Plot Top Genes"),
                 numericInput("n_plots", "Number of genes to plot (1-10)", value = 10, min = 1, max = 10),
                 actionButton("volcano", "Volcano Plot"),
                 actionButton("pca", "PCA Plot"),
                 actionButton("heatmap", "heatmap Plot")
               ),
               mainPanel(
                 plotOutput("genePlots"),
                 plotOutput("volcanoPlot"),
                 plotOutput("pcaPlot"),
                 plotOutput("heatPlot")
               )
             )
    )
  )
)




server <- function(input, output) {
  counts_df <- reactiveVal()
  coldata_df <- reactiveVal()
  
  analysis_results <- reactiveVal()
  dds_data <- reactiveVal()
  filtered_results <- reactiveVal()
  plots <- reactiveVal()
  volcano_plot <- reactiveVal()
  PCA_plot <- reactiveVal()
  heatmap_plot <- reactiveVal()
  
  output$countsStr <- renderPrint({
    req(input$countFile)
    counts <- read.csv(input$countFile$datapath)
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    #print(str(counts))
    
    counts_df(counts)
  })
  
  output$countsHead <- renderPrint({
    req(counts_df())
    counts <- counts_df()
    print(head(counts))
  })
  
  output$coldataStr <- renderPrint({
    req(input$colFile)
    coldata <- read.csv(input$colFile$datapath, stringsAsFactors = TRUE)
    coldata <- coldata[,-1]
    print(str(coldata))
    
    coldata_df(coldata)
  })
  
  output$coldataHead <- renderPrint({
    req(coldata_df())
    coldata <- coldata_df()
    print(head(coldata))
  })
  
  observeEvent(input$analyze, {
    req(coldata_df(), counts_df())
    
    counts <- counts_df()
    coldata <- coldata_df()
    

    countmat <- round(counts)
    #countmat = round(as.numeric((counts[-1,])))
    rownames(countmat) = rownames(counts)
    
    
    names(coldata) <- c("BioSample", "feature1", "feature2")
    

    dds = DESeqDataSetFromMatrix(countData = countmat,
                                 colData   = coldata,
                                 design    = ~ feature1 + feature2 + feature1:feature2)
    #print(dds)
    #print(counts(dds))
    keep <- which(rowSums(counts(dds)) > 2)
    dds <- dds[keep,]
    
    dds <- DESeq(dds)
    
    res_name <- resultsNames(dds)
    result <- lfcShrink(dds, coef = res_name[input$feature], type = "apeglm")
    result <- result[order(result$padj),]
    
    dds_data(dds)
    analysis_results(result)
  })
  
  output$resultTable <- renderTable({
    req(analysis_results())
    head(analysis_results())
  },rownames = T)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$filename, ".csv")
    },
    content = function(file) {
      write.csv(analysis_results(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$filter, {
    req(analysis_results())
    threshold <- c(input$lfc_threshold,input$padj_threshold)
    de_genes <- subset(analysis_results(), abs(log2FoldChange) > threshold[1] & padj < threshold[2])
    filtered_results(de_genes)
  })
  
  output$filteredTable <- renderTable({
    req(filtered_results())
    head(filtered_results())
  },rownames = T)
  
  output$downloadFiltered <- downloadHandler(
    filename = function() {
      paste0("filtered_", input$filename, ".csv")
    },
    content = function(file) {
      write.csv(filtered_results(), file, row.names = TRUE)
    }
  )
  
  observeEvent(input$plot, {
    req(filtered_results())
    top_genes <- rownames(head(filtered_results(), input$n_plots))
    plots(top_genes)
  })
  
  output$genePlots <- renderPlot({
    req(plots(),dds_data())
    dds <- dds_data() 
    coldata <- colData(dds)
    par(mfrow=c(2,ceiling(input$n_plots/2)))
    for (gene in plots()) {
      plot(factor(coldata$feature1), log(assay(dds)[gene,]), main = gene)
    
    }
  })
  
  observeEvent(input$volcano, {
    req(analysis_results(),  filtered_results())
    par(mfrow=c(1,1))
    results <- analysis_results()
    de_genes <- filtered_results()
    de_color = as.numeric(rownames(results) %in%  rownames(de_genes)) + 3
    df <- data.frame(x = results$log2FoldChange, y = -log10(results$padj))

    vo <- ggplot(df, aes(x,y, color = de_color)) + geom_point()
    volcano_plot(vo)
  })
  
  output$volcanoPlot <- renderPlot({
    req(volcano_plot())
  })
  
  observeEvent(input$pca, {
    req(dds_data() )
    dds <- dds_data()
    vds = vst(dds)
    counts_vst = assay(vds)
    grou <- colnames(colData(dds))
    pca <- DESeq2::plotPCA(vds, intgroup = c(grou[2], grou[3]))
    PCA_plot(pca)
  })
  
  output$pcaPlot <- renderPlot({
    req(PCA_plot())
  })
  
  observeEvent(input$pca, {
    req(dds_data())
    dds <- dds_data()
    vds = vst(dds)

    grou <- colnames(colData(dds))
    pca <- plotPCA(vds, intgroup = c(grou[2], grou[3]))
    PCA_plot(pca)
  })
  
  output$pcaPlot <- renderPlot({
    req(PCA_plot())
  })
  
  # observeEvent(input$heatmap, {
  #   output$heatPlot <- renderPlot({
  #     req(dds_data(), filtered_results())
  #     dds <- dds_data()
  #     de_genes <- filtered_results()
  #     vds = vst(dds)
  #     counts_vst = assay(vds)
  #     
  #     de_counts_vst = counts_vst[rownames(de_genes),]
  #     col_annot = as.data.frame(colData(dds)[,-1])
  #     print(str(de_counts_vst))
  #     pheatmap(de_counts_vst, annotation_col = col_annot, show_rownames = FALSE)
  #   })
  # })
  
  
}



shinyApp(ui = ui, server = server)



