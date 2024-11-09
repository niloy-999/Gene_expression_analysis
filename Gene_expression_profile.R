library(shiny)
library(DESeq2)
library(pheatmap)
library(VennDiagram)
library(plotly)
library(RColorBrewer)


ui <- fluidPage(
  titlePanel("Gene Expression Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload RNA sequence", accept = ".csv"),
      fileInput("file2", "Upload Meta data", accept = ".csv"),
      actionButton("analyze1", "Analyze the Differentially Expressed Gene"),
      actionButton("showTopGenes", "Show Top Expression gene in order")
    ),
    mainPanel(
      uiOutput("st")
    )
  )
)

server <- function(input, output) {
  
  data1 <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, header = TRUE, sep = ',')
  })
  
  
  
  data2 <- reactive({
    req(input$file2)
    read.csv(input$file2$datapath, header = TRUE, sep = ',')
  })
  
  
  output$data_table1 <- renderTable({
    data1()
  })
  
  output$data_table2 <- renderTable({
    data2()
  })
  
  
  
  
  output$analyze1 <- eventReactive(input$analyze1, {
    
    
    dds <- DESeqDataSetFromMatrix( countData = data1(),
                                   colData = data2(),
                                   design = ~ dexamethasone
    )
    
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    
    
    dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
    dds <- DESeq(dds)
    res <- results(dds)
    
    
    
    res  
  })
  
  
  
  output$data3 <- renderTable({
    req(input$analyze1)
    res <- results(dds)
    res
  })
  
  
  
  output$plot1 <- renderPlotly({
    req(input$analyze1)
    res <- results(dds)
    
    plot_ly(
      x = ~res$log2FoldChange, 
      y = ~-log10(res$padj), 
      
      mode = "markers",
      text = ~paste("Gene: ", rownames(res), "<br>log2FC: ", res$log2FoldChange, "<br>p-value: ", res$pvalue),
      hoverinfo = "text"  # Specify that the 'text' attribute is used for tooltips
    ) %>%
      add_markers(color = ~ifelse(res$padj < 0.05, "red", "blue")) %>%
      layout(
        title = "Volcano Plot", 
        xaxis = list(title = "log2FoldChange"), 
        yaxis = list(title = "-log10(padj)"),
        hovermode = "closest"
      )
  })
  
  
  
  
  output$plot5 <- renderPlot({
    req(input$analyze1)
    
    plotMA(res, alpha = 0.05)
  })
  
  
  
  output$plot2 <- renderPlotly({
    req(input$analyze1)
    
    
    disp_est_plot <- plot_ly(
      x = res$baseMean,
      y = res$dispGeneEst,
      mode = "markers",
      type = "scatter",
      marker = list(color = "643843"),
      text = ~paste("Gene: ", rownames(res), "<br>Dispersion: ", res$dispGeneEst),
      hoverinfo = "text"  # Specify that the 'text' attribute is used for tooltips
    ) %>%
      layout(title = "Dispersion Plot", xaxis = list(title = "Mean Expression"), yaxis = list(title = "Dispersion"))
    
    disp_est_plot
    
    
  })
  
  
  
  output$plot4 <- renderPlot({
    req(input$analyze1)
    plotDispEsts(dds)
  })
  
  
  
  
  output$plot3 <- renderPlotly({
    req(input$analyze1)
    top <- head(order(res$padj), 20)
    topc <- counts(dds)[top,]
    
    
    heatmap_plot <- plot_ly(
      z = topc,
      colorscale = "Viridis",
      type = "heatmap",
      colorbar = list(title = "Expression"),
      showscale = TRUE,
      hoverinfo = "analysis"  
    ) %>%
      layout(title = "HeatMap", xaxis = list(title = "Samples"), yaxis = list(title = "Genes"))
    
    heatmap_plot
  })
  
  
  
  
  output$plot6 <- renderPlot({
    req(input$analyze1)
    
    
    top_genes <- head(order(res$padj), 50)
    top_genes_expr <- assay(dds)[top_genes, ]
    
    dist_mat <- dist(t(top_genes_expr))
    hc <- hclust(dist_mat)
    
    
    heatmap(t(top_genes_expr)[hc$order, ], 
            col = colorRampPalette(c("navy", "yellow", "firebrick3"))(100),
            scale = "row", margins = c(10, 5))
  })
  
  
  
  
  
  output$showTopGenesTable <- renderTable({
    req(input$showTopGenes)
    
    
    top_genes <- head(order(res$padj), 10)
    top_genes_data <- data.frame(Gene = rownames(res)[top_genes],
                                 log2FoldChange = res$log2FoldChange[top_genes])
    top_genes_data
  })
  
  
  
  
  output$downloadPlot1 <- downloadHandler(
    filename = function() {
      "volcano_plot.png"
    },
    content = function(file) {
      ggsave(file, plot_ly_to_ggplotly(output$plot1))
    }
  )
  
  output$downloadPlot2 <- downloadHandler(
    filename = function() {
      "Dispersion_Plot.png"
    },
    content = function(file) {
      ggsave(file, plot_ly_to_ggplotly(output$plot2))
    }
  )
  
  output$downloadPlot3 <- downloadHandler(
    filename = function() {
      "HeatMap.png"
    },
    content = function(file) {
      ggsave(file, plot_ly_to_ggplotly(output$plot3))
    }
  )
  
  
  output$downloadPlot4 <- downloadHandler(
    filename = function() {
      "Dispersion_Plot.png"
    },
    content = function(file) {
      ggsave(file, plot_ly_to_ggplotly(output$plot4))
    }
  )
  
  output$downloadPlot5 <- downloadHandler(
    filename = function() {
      "volcano_plot.png"
    },
    content = function(file) {
      ggsave(file, plot_ly_to_ggplotly(output$plot5))
    }
  )
  
  output$downloadPlot6 <- downloadHandler(
    filename = function() {
      "clustered_heatmap.png"
    },
    content = function(file) {
      ggsave(file, plot_ly_to_ggplotly(output$plot6))
    }
  )
  
  
  
  
  output$st <- renderUI({
    tabsetPanel(
      tabPanel("RNA Sequence", tableOutput("data_table1")),
      tabPanel("Meta Data", tableOutput("data_table2")),
      tabPanel("Summary & Statistics", tableOutput("data3")),
      tabPanel("Volcano plot", plotlyOutput("plot1"),downloadButton("downloadPlot1", "Download Volcano Plot"), plotOutput("plot5"),downloadButton("downloadPlot5", "Download Volcano Plot")),
      tabPanel("Dispersion Plot", plotlyOutput("plot2"),downloadButton("downloadPlot2", "Download Dispersion Plot"),plotOutput("plot4"),downloadButton("downloadPlot4", "Download dispersion Plot-2")),
      tabPanel("HeatMap", plotlyOutput("plot3"),downloadButton("downloadPlot3", "Download heatmap Plot")),
      tabPanel("Clustered Heatmap", plotOutput("plot6"),downloadButton("downloadPlot6", "Download clustered heatmap Plot")),
      tabPanel("Top expressed Genes", tableOutput("showTopGenesTable"))
    )
  })
}


shinyApp(ui, server)
