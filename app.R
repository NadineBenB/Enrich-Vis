#############
# Libraries #
#############

if(!requireNamespace("DT", quietly = TRUE)) install.packages("DT")
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if(!requireNamespace("png", quietly = TRUE)) install.packages("png")
if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
if(!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if(!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")
if(!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
if(!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")

library(tidyverse)
library(DT)
library(ggplot2)
library(png)
library(pathview)
library(biomaRt)
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)

##################
# Functions file #
##################

source("functions.R")

###########################
# User Interface of Shiny #
###########################
ui <- fluidPage(
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      # Load input file
      fileInput("upload", NULL, accept = c(".csv", ".tsv")),
      # Clear button
      actionButton("reset", "Clear", icon = icon("refresh")),
      width = 3
    ),
    
    # Main panel
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Home",
                           h3("Enrichement analysis using Gene Ontolgy and KEGG databases"),
                           p("Using this application and the set of deferentially expressed genes (DEGs), we will
                             first generate the GO enrichment analysis to find which GO terms are over-represented. 
                             Then we will perform a KEGG pathway enrichment analysis. It determine which pathway are over-represented given 
                             the same set of DEGs."),
                           br(),
                           p("Load your list of DEGs, it will generate a table with the results of the GO enrichment
                             analysis, then three different figures to represent the results :"),
                           p("- A dot plot of the enrichment term that represent the top enrichment term,"), 
                           p("- A gene-concept network that simply give a network with the enrichment terms and
                             associated genes,"), 
                           p("- A tree plot that groups the enriched terms according to their similarities."),
                           p("Then, for the KEGG pathway enrichment analysis, it generates a table with the results, and we
                             can retrieve the KEGG pathway diagrams with the DEGs highlighted according to their log
                             fold change."),
                           br(),
                           p("After loading the data, the analysis can take some time, wait untill the progress bar is complte! Also, according to the results of the GO and
                             KEGG enrichment analysis, some or all figures wont be generated and a warning will be issued (it is
                             the case for the DEGs issued from the analysis MutantP8 vs MutantP12 and WTP8 vs WTP1")),
                  
                  tabPanel("GO Enrichement table",
                           p("This table contains the results of the GO enrichment analysis using terms related to biological
                             processes (BP), cellular compartment (CC) and molecular functions (MF). You can use the 
                             search tabe to select specific term or gene. You can also select the sub categoy of the ontology (CC, BP, MF).
                             By defaults the terms are organized in order of significativity (adjusted p-value)."),
                           selectInput("ont3",
                                       label = "Gene Ontology",
                                       choices = c("BP","MF","CC")),
                           dataTableOutput("dtable")),
                  
                  tabPanel("Dot plot",
                           p("Dot plot of the top 10 most significant termns (per sub category). You can choose the number
                             of terms to be represented. Each term is ploted agaisnt the gene ratio: number of genes represented 
                             by the term, over the total number DEGs that have an anotation in the selected category. The size of the 
                             dot correspond to the number of genes represented by the term and the color is the associated adjusted p-value."),
                           numericInput("n1", label = "Top significant terms to show", value = 10),
                           plotOutput("dotplot",height = 1000)),
                  
                  tabPanel("Gene-concept network",
                           p("Representation of the top 5 signigicant terms of the BP sub category, and the associated genes. 
                             You can change the number of terms to show and the category. The two plots side by side are the same, one with the terms name,
                             and one with the genes name."),
                           fluidRow(
                             column(3, numericInput("n2", label = "Top significant terms to show", value = 5)),
                             column(3, selectInput("ont1",
                                      label = "Gene Ontology",
                                      choices = c("BP","MF","CC")))),
                           plotOutput("cneplot",height = 1000),
                           p("Here you can choose the term you want to plot"),
                           htmlOutput("terms"),
                           plotOutput("cneplot2",height = 1000)),
                  
                  tabPanel("Tree plot",
                           p("A tree plot that summarizes the similarities between the enriched terms in the BP category. 
                             You can change the category and the number of cluster to create. The size of the dots represent the number of genes per term, 
                             and the colors the p-value."),
                           fluidRow(
                             column(3, numericInput("n3", label = "Number of clusters", value = 5)),
                             column(3, selectInput("ont2",
                                                   label = "Gene Ontology",
                                                   choices = c("BP","MF","CC")))),
                           plotOutput("treeplot",height = 800)),
                  
                  #tabPanel("Enrichment map",
                  #         plotOutput("emapplot")),
                  
                  tabPanel("KEGG pathways enrichement",
                           p("The table contains the result of the KEGG pathway enrichement analysis"),
                           dataTableOutput("resKEGG"),
                           p("Here you can choose the KEGG pathway to visualize, the DEGs are highlited according to their log fold change."),
                           htmlOutput("listSigPathways"),
                           plotOutput("pathKEGG",height = 1000))
      )
    )
  )
)

############################
# Server function of Shiny #
############################

server = function(input, output, session) {
  
  data = reactiveValues()
  
  observeEvent(input$upload, {
    withProgress(message = 'Performing the analysis, please wait ...', value = 0,{
      data$DEG = read.csv(input$upload$datapath, sep = ";", header = TRUE, row.names = 1)
      incProgress(0.25)
      data$EGO = enrich_GO(rownames(data$DEG), "ALL")
      incProgress(0.25)
      data$readEGO = setReadable(data$EGO, 'org.Mm.eg.db')
      incProgress(0.25)
      data$EKEGG = enrich_KEGG(data$DEG)
      incProgress(0.25)
  })})
  
  observeEvent(input$reset, {
    session$reload()
  })
  
  output$dtable = renderDT(as.data.frame(filter(data$readEGO, ONTOLOGY == input$ont3)))
  
  dotplot = reactive(dotPlot(data$EGO, input$n1))
  output$dotplot = renderPlot(dotplot())
  
  cneplot = reactive(cnePlot(data$EGO, data$DEG, input$n2, input$ont1))
  output$cneplot = renderPlot(cneplot())
  
  output$terms = renderUI({
    choices = data$EGO[, 3]
    selectInput("term",
                label = "Select a significant GO term to show",
                choices = choices)
  })
  cneplot2 = reactive(cnePlot(data$EGO, data$DEG, input$term, input$ont1))
  output$cneplot2 = renderPlot(cneplot2())
  
  treeplot = reactive(treePlot(data$EGO, input$ont2, input$n3))
  output$treeplot = renderPlot(treeplot())
  
  #emapplot = reactive(emapPlot(data$EGO))
  #output$emapplot = renderPlot(emapplot())
  

  output$resKEGG = renderDT(filter(data$EKEGG@result[-1], p.adjust < 0.05))
  output$listSigPathways = renderUI({
    data = filter(data$EKEGG@result, p.adjust < 0.05)
    choices = data[, 2]
    selectInput("sigPathways",
                label = "Significant KEGG pathway",
                choices = choices)
  })
  KEGG = reactive(view_KEGG (input$sigPathways, data$EKEGG, data$DEG))
  output$pathKEGG = renderPlot(KEGG())
}

# Run the application
shinyApp(ui = ui, server = server)
