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
library(gridExtra)
library(grid)
library(ggrepel)
library(ggfun)
library(tidyverse)

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
                  choices = c("pvalue" = "pvalue","padj" = "padj")),
      #numericInput("filterThreshold", "Filter Threshold", value = 0.05, step = 0.01),
      numericInput("filterThreshold", "Filter Threshold", value = 0.05),
      numericInput("foldChange", "log2 Fold Change Cutoff", value = 1),
      textInput("comparison", "Group Comparison (e.g., control,case)", value = "control,case"),
      selectInput("species", "Species",
                  choices = c("Human" = "org.Hs.eg.db", "Mouse" = "org.Mm.eg.db")),
      #checkboxGroupInput("enrichMethods", "Enrichment Methods",
      #                   choices = c("GO" = "go", "KEGG" = "kegg", "Hallmark" = "hallmark", "Reactome" = "reactome", "Wikipathways" = "wikipathways"),
      #                   selected = "go"),
      selectInput("enrichMethods", "Enrichment Methods",
            choices = c("GO" = "go", "KEGG" = "kegg","Reactome" = "reactome"
                        #"Hallmark" = "hallmark", "Wikipathways" = "wikipathways"
            )),
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
    up_enriched = list(),
    down_enriched = list()
  )

  observeEvent(input$runAnalysis, {
    req(input$countFile, input$metaFile)
    print(paste("Input file path:", input$countFile$datapath))
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
    comparison <- reactive(strsplit(input$comparison, ",")[[1]])
    #print(groupCol)
    #print(comparison)
    #print(data$meta)
    if (!(groupCol %in% colnames(data$meta)) || length(comparison()) != 2) {
      showNotification("Invalid group column or comparison format.", type = "error")
      return()
    }

    # DEG analysis
    method <- input$method
    tryCatch({
      if (method == "deseq2") {

        data$meta <- data$meta[order(match(data$meta[[groupCol]],comparison())),,drop=FALSE]
        conditon <- data$meta[[groupCol]]
        conditon <- factor(conditon,levels=comparison())
        genecolData <- data.frame(row.names=colnames(data$counts)[1:length(colnames(data$counts))], conditon)
        genecountData <- data$counts[, rownames(genecolData)]

        dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(round(genecountData)),
        colData = as.data.frame(genecolData),
        #design = as.formula(paste0("~", condition)
        design = ~ conditon
    )

        dds <- DESeq(dds)
        #res <- results(dds, contrast = c(groupCol, comparison[1], comparison[2]))
        res <- results(dds)
        res <- res %>% tibble::rownames_to_column(var = "GeneID")

    } else if (method == "wilcoxon") {
      logtpm <- log2(data$counts + 1)
      group1 <- which(data$meta[[groupCol]] == comparison()[1])
      group2 <- which(data$meta[[groupCol]] == comparison()[2])
      pvalue <- apply(logtpm, 1, function(x) {wilcox.test(x[group1], x[group2])$p.value})
      res <- data.frame(
        row.names = rownames(logtpm),
        GeneID = rownames(logtpm),
        log2FoldChange = rowMeans(logtpm[, group2]) - rowMeans(logtpm[, group1]),
        pvalue = pvalue,
        padj = p.adjust(pvalue, method = "BH")
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
    fcCutoff <- input$foldChange
    #data$results <- res[res[[filterMethod]] < threshold, ]
    data$results <- res
    up_genes <- rownames(data$results[data$results$log2FoldChange > fcCutoff & data$results[[filterMethod]] < threshold, ])
    down_genes <- rownames(data$results[data$results$log2FoldChange < -fcCutoff & data$results[[filterMethod]] < threshold, ])
    data$results$sig <- ifelse(rownames(data$results) %in% up_genes,"Up",
                               ifelse(rownames(data$results) %in% down_genes,"Down","None"))
    # Enrichment analysis
    #deg_genes <- rownames(data$results)
    deg_genes <- c(up_genes,down_genes)
    #enrich_methods <- input$enrichMethods
    species_db <- input$species

    if (length(deg_genes) == 0) {
       showNotification("No DEGs passed the filter criteria. Enrichment analysis skipped.", type = "warning")
        return()
    }else{
      if (grepl('^ENSG',deg_genes[1])){
    up_deg_genes_changed <- bitr(up_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=get(species_db))
    up_deg_genes_changed <- up_deg_genes_changed$ENTREZID
    down_deg_genes_changed <- bitr(down_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=get(species_db))
    down_deg_genes_changed <- down_deg_genes_changed$ENTREZID
    }else if (grepl('^NG_',deg_genes[1])| grepl('^NM_',deg_genes[1]) |grepl('^NP_',deg_genes[1]) ){
    up_deg_genes_changed <- bitr(up_genes, fromType="REFSEQ", toType="ENTREZID", OrgDb=get(species_db))
    up_deg_genes_changed <- up_deg_genes_changed$ENTREZID
    down_deg_genes_changed <- bitr(down_genes, fromType="REFSEQ", toType="ENTREZID", OrgDb=get(species_db))
    down_deg_genes_changed <- down_deg_genes_changed$ENTREZID
    }else if (grepl('^[0-9]+',deg_genes[1])){
    up_deg_genes_changed <- up_genes
    down_deg_genes_changed <- down_genes
    }else {
    up_deg_genes_changed <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=get(species_db))
    up_deg_genes_changed <- up_deg_genes_changed$ENTREZID
    down_deg_genes_changed <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=get(species_db))
    down_deg_genes_changed <- down_deg_genes_changed$ENTREZID
    }
    }

    #data$enriched <- list()
    data$up_enriched$go <- enrichGO(
        gene = up_deg_genes_changed,
        OrgDb = species_db,
        #keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable=TRUE)

    data$up_enriched$kegg <- enrichKEGG(
        gene = up_deg_genes_changed,
        organism = ifelse(species_db == "org.Hs.eg.db", "hsa", "mmu"),
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        use_internal_data = F
      )


    data$up_enriched$reactome <- enrichPathway(
        gene = up_deg_genes_changed,
        organism = ifelse(species_db == "org.Hs.eg.db", "human", "mouse"),
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable=TRUE
      )

    # Add logic for Hallmark and Wikipathways as needed
    
    data$down_enriched$go <- enrichGO(
      gene = down_deg_genes_changed,
      OrgDb = species_db,
      #keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable=TRUE)
    
    data$down_enriched$kegg <- enrichKEGG(
      gene = down_deg_genes_changed,
      organism = ifelse(species_db == "org.Hs.eg.db", "hsa", "mmu"),
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      use_internal_data = F
    )
    
    
    data$down_enriched$reactome <- enrichPathway(
      gene = down_deg_genes_changed,
      organism = ifelse(species_db == "org.Hs.eg.db", "human", "mouse"),
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable=TRUE
    )
    
  }
  )


  output$degTable <- renderDataTable({
    req(data$results)
    datatable(data$results, options = list(pageLength = 10), rownames = TRUE)
  })
  
  # output$volcanoPlot <- renderPlot({
  #   req(data$results)
  #   req(input$comparison)
  #   comparison <- strsplit(input$comparison, ",")[[1]]
  #   res <- data$results
  #   EnhancedVolcano(
  #     res,
  #     lab = rownames(res),
  #    # title = paste0(comparison[2]," vs ",comparison[1]),  ## case vs control
  #     x = 'log2FoldChange',
  #     #y = ifelse(input$filterMethod == "pvalue","-Log10P","-Log10Padj"),
  #     y = input$filterMethod, 
  #     pCutoff = input$filterThreshold,
  #     FCcutoff = input$foldChange,
  #     cutoffLineType = "twodash",
  #     cutoffLineWidth = 0.8,
  #     pointSize = c(ifelse(res$sig != "None",6,4)),
  #     labSize = 6.0,
  #     #col = c("grey",""),
  #     colAlpha = 1,
  #     legendLabels = c("Not Sig","Log2FC","Pvalue","Pvalue&Log2FC"),
  #     legendPosition = "top",
  #     legendLabSize = 12,
  #     legendIconSize = 4.0,
  #     drawConnectors = TRUE,
  #     widthConnectors = 0.75
  #   )
  # })
  output$volcanoPlot <- renderPlot({
    req(data$results)
    req(input$comparison)
    req(input$filterMethod)
    comparison <- strsplit(input$comparison, ",")[[1]]
    print(names(data$results))
    res <- data$results
    logFC.limit <- input$foldChange
    
    filtermethod <- input$filterMethod
    filterThreshold <- input$filterThreshold
    group_name <- paste0(comparison[2]," vs ",comparison[1])
    if(filtermethod == "pvalue"){
    
    ggplot(data = res) + 
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), 
                     color = log2FoldChange,
                     size = -log10(pvalue))) + 
      geom_point(data =  res %>%
                   tidyr::drop_na() %>%
                   dplyr::filter(sig == "Up") %>%
                   dplyr::arrange(desc(-log10(pvalue))) %>%
                   dplyr::slice(1:8),
                 aes(x = log2FoldChange, y = -log10(pvalue),
                     # fill = log2FoldChange,
                     size = -log10(pvalue)),
                 shape = 21, show.legend = F, color = "#000000") +
      geom_text_repel(data =  res %>% 
                        tidyr::drop_na() %>% 
                        dplyr::filter(sig == "Up") %>%
                        dplyr::arrange(desc(-log10(pvalue))) %>%
                        dplyr::slice(1:8),
                      aes(x = log2FoldChange, y = -log10(pvalue), label = GeneID),
                      box.padding = 0.5,
                      nudge_x = 0.5,
                      nudge_y = 0.2,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      # segment.angle = 10,
                      direction = "y", 
                      hjust = "left"
      ) + 
      geom_point(data =  res %>%
                   tidyr::drop_na() %>%
                   dplyr::filter(sig == "Down") %>%
                   dplyr::arrange(desc(-log10(pvalue))) %>%
                   dplyr::slice(1:8),
                 aes(x = log2FoldChange, y = -log10(pvalue),
                     # fill = log2FoldChange,
                     size = -log10(pvalue)),
                 shape = 21, show.legend = F, color = "#000000") +
      geom_text_repel(data =  res %>% 
                        tidyr::drop_na() %>% 
                        dplyr::filter(sig == "Down") %>%
                        dplyr::arrange(desc(-log10(pvalue))) %>%
                        dplyr::slice(1:8),
                      aes(x = log2FoldChange, y = -log10(pvalue), label = GeneID),
                      box.padding = 0.5,
                      nudge_x = -0.2,
                      nudge_y = 0.2,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20,
                      direction = "y", 
                      hjust = "right"
      ) + 
      scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                            values = seq(0, 1, 0.2)) +
      scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                           values = seq(0, 1, 0.2)) +
      geom_vline(xintercept = c(-logFC.limit, logFC.limit), linetype = 2) +
      geom_hline(yintercept = -log10(filterThreshold), linetype = 4) + 
      scale_size(range = c(1,7)) + 
      ggtitle(#label = "Volcano Plot",
        #subtitle = "volcano plot"
        label = group_name
      ) + 
      #xlim(c(-3, 3)) + 
      #ylim(c(-1, 6)) + 
        xlim(min(res$log2FoldChange, na.rm = TRUE) - 0.5, 
             max(res$log2FoldChange, na.rm = TRUE) + 0.5)+
        ylim(0, max(-log10(res$pvalue), na.rm = TRUE) + 1)+  # Auto-adjust the upper limit
      theme_bw() + 
      theme(panel.grid = element_blank(),
            legend.background = element_roundrect(color = "#808080", linetype = 1),
            axis.text = element_text(size = 13, color = "#000000"),
            axis.title = element_text(size = 15),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
      ) + 
      annotate(geom = "text", x = 2.5, y = 0.25, label = "p = 0.05", size = 5) + 
      coord_cartesian(clip = "off") + 
      annotation_custom(
        grob = grid::segmentsGrob(
          y0 = unit(-10, "pt"),
          y1 = unit(-10, "pt"),
          arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
          gp = grid::gpar(lwd = 3, col = "#74add1")
        ), 
        xmin = -3, 
        xmax = -1,
        ymin = 5.5,
        ymax = 5.5
      ) +
      annotation_custom(
        grob = grid::textGrob(
          label = "Down",
          gp = grid::gpar(col = "#74add1")
        ),
        xmin = -3, 
        xmax = -1,
        ymin = 5.5,
        ymax = 5.5
      ) +
      annotation_custom(
        grob = grid::segmentsGrob(
          y0 = unit(-10, "pt"),
          y1 = unit(-10, "pt"),
          arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
          gp = grid::gpar(lwd = 3, col = "#d73027")
        ), 
        xmin = 1, 
        xmax = 3,
        ymin = 5.5,
        ymax = 5.5
      ) +
      annotation_custom(
        grob = grid::textGrob(
          label = "Up",
          gp = grid::gpar(col = "#d73027")
        ),
        xmin = 1, 
        xmax = 3,
        ymin = 5.5,
        ymax = 5.5
      ) 
    
    }else if(filtermethod == "padj"){
      ggplot(data = res) + 
        geom_point(aes(x = log2FoldChange, y = -log10(padj), 
                       color = log2FoldChange,
                       size = -log10(padj))) + 
        geom_point(data =  res %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(sig == "Up") %>%
                     dplyr::arrange(desc(-log10(padj))) %>%
                     dplyr::slice(1:8),
                   aes(x = log2FoldChange, y = -log10(padj),
                       # fill = log2FoldChange,
                       size = -log10(padj)),
                   shape = 21, show.legend = F, color = "#000000") +
        geom_text_repel(data =  res %>% 
                          tidyr::drop_na() %>% 
                          dplyr::filter(sig == "Up") %>%
                          dplyr::arrange(desc(-log10(padj))) %>%
                          dplyr::slice(1:8),
                        aes(x = log2FoldChange, y = -log10(adj), label = GeneID),
                        box.padding = 0.5,
                        nudge_x = 0.5,
                        nudge_y = 0.2,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        # segment.angle = 10,
                        direction = "y", 
                        hjust = "left"
        ) + 
        geom_point(data =  res %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(sig == "Down") %>%
                     dplyr::arrange(desc(-log10(padj))) %>%
                     dplyr::slice(1:8),
                   aes(x = log2FoldChange, y = -log10(padj),
                       # fill = log2FoldChange,
                       size = -log10(padj)),
                   shape = 21, show.legend = F, color = "#000000") +
        geom_text_repel(data =  res %>% 
                          tidyr::drop_na() %>% 
                          dplyr::filter(sig == "Down") %>%
                          dplyr::arrange(desc(-log10(padj))) %>%
                          dplyr::slice(1:8),
                        aes(x = log2FoldChange, y = -log10(padj), label = GeneID),
                        box.padding = 0.5,
                        nudge_x = -0.2,
                        nudge_y = 0.2,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y", 
                        hjust = "right"
        ) + 
        scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                              values = seq(0, 1, 0.2)) +
        scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                             values = seq(0, 1, 0.2)) +
        geom_vline(xintercept = c(-logFC.limit, logFC.limit), linetype = 2) +
        geom_hline(yintercept = -log10(filterThreshold), linetype = 4) + 
        scale_size(range = c(1,7)) + 
        ggtitle(#label = "Volcano Plot",
          #subtitle = "volcano plot"
          label = group_name
        ) + 
        #xlim(c(-3, 3)) + 
        #ylim(c(-1, 6)) + 
        xlim(min(res$log2FoldChange, na.rm = TRUE) - 0.5, 
             max(res$log2FoldChange, na.rm = TRUE) + 0.5)+
        ylim(0, max(-log10(res$pvalue), na.rm = TRUE) + 1)+  # Auto-adjust the upper limit
        theme_bw() + 
        theme(panel.grid = element_blank(),
              legend.background = element_roundrect(color = "#808080", linetype = 1),
              axis.text = element_text(size = 13, color = "#000000"),
              axis.title = element_text(size = 15),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
        ) + 
        annotate(geom = "text", x = 2.5, y = 0.25, label = "p = 0.05", size = 5) + 
        coord_cartesian(clip = "off") + 
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
            gp = grid::gpar(lwd = 3, col = "#74add1")
          ), 
          xmin = -3, 
          xmax = -1,
          ymin = 5.5,
          ymax = 5.5
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label = "Down",
            gp = grid::gpar(col = "#74add1")
          ),
          xmin = -3, 
          xmax = -1,
          ymin = 5.5,
          ymax = 5.5
        ) +
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
            gp = grid::gpar(lwd = 3, col = "#d73027")
          ), 
          xmin = 1, 
          xmax = 3,
          ymin = 5.5,
          ymax = 5.5
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label = "Up",
            gp = grid::gpar(col = "#d73027")
          ),
          xmin = 1, 
          xmax = 3,
          ymin = 5.5,
          ymax = 5.5
        ) 
    }
    
  })
  


  output$enrichPlot <- renderPlot({
    req(data$up_enriched,data$down_enriched)
    #req(comparison())  # Ensure comparison() is ready
    req(input$comparison)
    comparison <- strsplit(input$comparison, ",")[[1]]
    selected_method <- input$enrichMethods
    up_enrich_result <- data$up_enriched[[selected_method]]
    down_enrich_result <- data$down_enriched[[selected_method]]
    #print(enrich_result)
    if (!is.null(up_enrich_result)) {
       p1 <-  dotplot(up_enrich_result, showCategory = 20) +
         theme(plot.margin = margin(10, 10, 10, 10)) + 
          ggtitle(paste(comparison[2],"Up Enrichment Analysis -", toupper(selected_method)))
    } else {
        showNotification(paste("No up enrichment results for", selected_method), type = "warning")
    }
    
    if (!is.null(down_enrich_result)) {
     p2 <- dotplot(down_enrich_result, showCategory = 20) +
       theme(plot.margin = margin(10, 10, 10, 10)) +
        ggtitle(paste(comparison[2],"Down Enrichment Analysis -", toupper(selected_method)))
    } else {
      showNotification(paste("No down enrichment results for", selected_method), type = "warning")
    }
    grid.arrange(p1,p2,ncol=1,heights = c(1, 1)
                 #top = textGrob("enrichment",
                 #                just=c("center"),
                #                gp = gpar(fontsize = 20)
                #                )
                 )
    
})

}


shinyApp(ui = ui, server = server)

