# Title     : Shiny script for RNA DEG analysis
# Objective : DEG, Enrichment, Draw
# Created by: stl
# Created on: 2025/1/13

library(shiny)
library(DESeq2)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(EnhancedVolcano)
library(DT)
library(readxl)
library(org.Hs.eg.db) # Default human, can be changed dynamically

ui <- fluidPage(
  titlePanel("RNA-seq DEG Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("countFile", "Upload Input Matrix (CSV/TXT/XLS)", accept = c(".csv",".txt", ".xls", ".xlsx")),
      fileInput("metaFile", "Upload Metadata File (CSV/TXT/XLS)", accept = c(".csv",".txt", ".xls", ".xlsx")),
      textInput("groupCol", "Metadata Group Column", value = "group"),
      selectInput("method", "Method for DEG Analysis",
                  choices = c("DESeq2(Count)" = "deseq2", "Wilcoxon (TPM/FPKM)" = "wilcoxon")),
      selectInput("filterMethod", "Filter Method",
                  choices = c("padj" = "padj", "pvalue" = "pvalue")),
      numericInput("filterThreshold", "Filter Threshold", value = 0.05, step = 0.01),
      numericInput("foldChange", "log2 Fold Change Cutoff", value = 1, step = 0.1),
      textInput("comparison", "Group Comparison (e.g., control,case)", value = "A,B"),
      selectInput("species", "Species",
                  choices = c("Human" = "org.Hs.eg.db", "Mouse" = "org.Mm.eg.db")),
      checkboxGroupInput("enrichMethods", "Enrichment Methods",
                         choices = c("GO" = "go", "KEGG" = "kegg", "Hallmark" = "hallmark", "Reactome" = "reactome", "Wikipathways" = "wikipathways"),
                         selected = "go"),
      #selectInput("enrichMethods", "Enrichment Methods",
      #      choices = c("GO" = "go", "KEGG" = "kegg","Reactome" = "reactome"
      #                  #"Hallmark" = "hallmark", "Wikipathways" = "wikipathways"
      #      )),
      actionButton("runAnalysis", "Run Analysis")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("DEG Table", dataTableOutput("degTable")),
        tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
        tabPanel("Enrichment Analysis", plotOutput("enrichPlot")),
        tabPanel("Download", plotOutput("download"))
      )
    )
  )
)

server <- function(input, output, session) {

  # Reactive values to store data
  data <- reactiveValues(
    counts = NULL,
    meta = NULL,
    results = NULL,
    enriched = list()
  )

  observeEvent(input$runAnalysis, {
    req(input$countFile, input$metaFile)
    print(paste("Count file path:", input$countFile$datapath))
    print(paste("Metadata file path:", input$metaFile$datapath))

    # Load count data
    tryCatch({
       if (grepl("\\.csv$", input$countFile$name)) {
        data$counts <- read.csv(input$countFile$datapath, row.names = 1)
       }else if (grepl("\\.txt$|\\.xls", input$countFile$name)) {
        data$counts <- as.data.frame(read.table(input$countFile$datapath, row.names = 1,
                                                header=TRUE,sep="\t",check.names=FALSE))
       }else if (grepl("\\.xlsx$", input$countFile$name)) {
        data$counts <- as.data.frame(read_excel(input$countFile$datapath, col_names = TRUE))
        rownames(data$counts) <- data$counts[[1]]
         data$counts <- data$counts[,-1]
  }}, error = function(e) {
        showNotification(paste("Error loading input file:", e$message), type = "error")
         return()
})

    # Load metadata
    tryCatch({
        if (grepl("\\.csv$", input$metaFile$name)) {
        data$meta <- read.csv(input$metaFile$datapath,row.names = 1)
        }else if (grepl("\\.txt$|\\.xls", input$metaFile$name)) {
          data$meta <- as.data.frame(read.table(input$metaFile$datapath,sep="\t",header=TRUE,row.names = 1))
        }else if (grepl("\\.xlsx$", input$metaFile$name)) {
          data$meta <- as.data.frame(read_excel(input$metaFile$datapath, col_names = TRUE))
          rownames(data$meta) <- data$meta[[1]]
          data$meta <- data$meta[,-1]
        }}, error = function(e) {
          showNotification(paste("Error loading metadata file:", e$message), type = "error")
          return()
})

    # Check if data is loaded correctly
    if (is.null(data$counts) || is.null(data$meta)) {
      showNotification("Failed to load input files.", type = "error")
      return()
    }

    # Check group column and split comparison groups
    groupCol <- input$groupCol
    comparison <- strsplit(input$comparison, ",")[[1]]
    #print(groupCol)
    #print(comparison)
    #print(data$meta)
    if (!(groupCol %in% colnames(data$meta)) || length(comparison) != 2) {
      showNotification("Invalid group column or comparison format.", type = "error")
      return()
    }

    # DEG analysis
    method <- input$method
    tryCatch({
      if (method == "deseq2") {

        data$meta <- data$meta[order(match(data$meta[[groupCol]],comparison)),,drop=FALSE]
        conditon <- data$meta[[groupCol]]
        conditon <- factor(conditon,levels=comparison)
        genecolData <- data.frame(row.names=colnames(data$counts)[1:length(colnames(data$counts))], conditon)
        genecountData <- data$counts[, rownames(genecolData)]

        dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(round(genecountData)),
        colData = as.data.frame(genecolData),
        #design = as.formula(paste0("~", condition)
        design = ~ conditon
    )

        dds <- DESeq(dds)
        res <- results(dds, contrast = c(groupCol, comparison[1], comparison[2]))
        #res <- results(dds)

    } else if (method == "wilcoxon") {
      tpm <- log2(data$counts + 1)
      group1 <- which(data$meta[[groupCol]] == comparison[1])
      group2 <- which(data$meta[[groupCol]] == comparison[2])
      pvalues <- apply(tpm, 1, function(x) {wilcox.test(x[group1], x[group2])$p.value})
      res <- data.frame(
        row.names = rownames(tpm),
        log2FoldChange = rowMeans(tpm[, group2]) - rowMeans(tpm[, group1]),
        pvalue = pvalues,
        padj = p.adjust(pvalues, method = "BH")
      )
    }
  res <- as.data.frame(res)
  res <- res[complete.cases(res), ]
}, error = function(e) {
  showNotification(paste("Error in DEG analysis:", e$message), type = "error")
  res <- NULL
})

    # Apply filter
    filterMethod <- input$filterMethod
    threshold <- input$filterThreshold
    data$results <- res[res[[filterMethod]] < threshold, ]

    # Enrichment analysis
    deg_genes <- rownames(data$results)
    #enrich_methods <- input$enrichMethods
    species_db <- input$species

    if (length(deg_genes) == 0) {
       showNotification("No DEGs passed the filter criteria. Enrichment analysis skipped.", type = "warning")
        return()
    }else{
      if (grepl('^ENSG',deg_genes[1])){
    deg_genes_changed <- bitr(deg_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=get(species_db))
    deg_genes_changed <- deg_genes_changed$ENTREZID
    }else if (grepl('^NG_',deg_genes[1])| grepl('^NM_',deg_genes[1]) |grepl('^NP_',deg_genes[1]) ){
    deg_genes_changed <- bitr(deg_genes, fromType="REFSEQ", toType="ENTREZID", OrgDb=get(species_db))
    deg_genes_changed <- deg_genes_changed$ENTREZID
    }else if (grepl('^[0-9]+',deg_genes[1])){
    deg_genes_changed <- deg_genes
    }else {
    deg_genes_changed <- bitr(deg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=get(species_db))
    deg_genes_changed <- deg_genes_changed$ENTREZID
    }
    }

    #data$enriched <- list()
    data$enriched$go <- enrichGO(
        gene = deg_genes_changed,
        OrgDb = species_db,
        #keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable=TRUE)

    data$enriched$kegg <- enrichKEGG(
        gene = deg_genes_changed,
        organism = ifelse(species_db == "org.Hs.eg.db", "hsa", "mmu"),
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        use_internal_data = F
      )


    data$enriched$reactome <- enrichPathway(
        gene = deg_genes_changed,
        organism = ifelse(species_db == "org.Hs.eg.db", "human", "mouse"),
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable=TRUE
      )

    # Add logic for Hallmark and Wikipathways as needed
  }
  )


  output$degTable <- renderDataTable({
    req(data$results)
    datatable(data$results, options = list(pageLength = 10), rownames = TRUE)
  })

  output$volcanoPlot <- renderPlot({
    req(data$results)
    res <- data$results
    EnhancedVolcano(
      res,
      lab = rownames(res),
      x = 'log2FoldChange',
      y = input$filterMethod,
      pCutoff = input$filterThreshold,
      FCcutoff = input$foldChange,
      pointSize = 2.0,
      labSize = 3.0
    )
  })

#  selectedEnrichPlot <- reactive({
#    req(data$enriched, input$enrichMethods)
#    enrich_result <- data$enriched[[input$enrichMethods]]
#    if (!is.null(enrich_result)) {
#        dotplot(enrich_result, showCategory = 20)
#    } else {
#        NULL
#    }
#})

  output$enrichPlot <- renderPlot({
    req(data$enriched)
    selected_method <- input$enrichMethods
    enrich_result <- data$enriched[[selected_method]]
    #print(enrich_result)
    if (!is.null(enrich_result)) {
        dotplot(enrich_result, showCategory = 20) +
          ggtitle(paste("Enrichment Analysis -", toupper(selected_method)))
    } else {
        showNotification(paste("No enrichment results for", selected_method), type = "warning")
    }
})


#  output$enrichPlot <- renderPlot({
#    req(data$enriched)
#    plots <- lapply(data$enriched, function(enrich) {dotplot(enrich, showCategory = 20)})
#    gridExtra::grid.arrange(grobs = plots, ncol = 1)
#  })
}

shinyApp(ui = ui, server = server)

