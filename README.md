# Shiny App for Gene Ontology and KEGG Pathway Enrichment Analysis

This repository hosts a Shiny application for performing and visualising enrichment analysis using Gene Ontology (GO) and KEGG pathway databases. The app provides interactive visualisations and detailed results based on a set of differentially expressed genes (DEGs). 


- **GO enrichment analysis**: 
  - Dot plots for the top enrichment terms.
  - Gene-concept networks showing terms and associated genes.
  - Tree plots clustering enriched terms based on similarity.

- **KEGG pathway enrichment analysis**: 
  - Results table with enriched pathways.
  - Visualisation of KEGG pathways with DEGs highlighted by their log fold change.

- **Interactive interface**:
  - Upload your DEG list.
  - Explore results through dynamic tables and plots.
  - Retrieve pathway diagrams directly from the KEGG database.


## Installation

Clone this repository and ensure the required R packages are installed. Install the required R libraries using the following script:

```
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("DT", quietly = TRUE)) install.packages("DT")
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if(!requireNamespace("png", quietly = TRUE)) install.packages("png")
if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if(!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
if(!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if(!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")
if(!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
if(!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
```

## Usage

Place the app.R and functions.R files in the same directory. Run the application using R or RStudio: ```shiny::runApp("app.R")```
Upload your DEG list in .csv or .tsv format. The analysis will begin automatically. Explore the results across various tabs in the app interface.
