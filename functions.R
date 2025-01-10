# GO enrichment #
#################
enrich_GO = function(genes, ont) {
  enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = 'ENSEMBL',
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  )
}

# GO visualization #
####################
dotPlot = function (ego, n_category) {
  dotplot(ego, showCategory = n_category, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free") + scale_y_discrete(guide = guide_axis(check.overlap = TRUE))
}

cnePlot = function (ego, DEG, n_category, ontology) {
  ego_read = setReadable(ego, 'org.Mm.eg.db')
  if (ontology != "ALL"){
    ego_read = filter(ego_read, ONTOLOGY == ontology)
  }
  genes = DEG$log2FoldChange
  names(genes) = rownames(DEG)
  p1 = cnetplot(
    ego_read,
    showCategory = n_category,
    node_label = "gene",
    foldChange = genes
  )
  p2 = cnetplot(
    ego_read,
    showCategory = n_category,
    node_label = "category",
    foldChange = genes
  )
  p1 + p2
}

treePlot = function (ego, ontology,n) {
  ego_read = setReadable(ego, 'org.Mm.eg.db')
  if (ontology != "ALL"){
    ego_read = filter(ego_read, ONTOLOGY == ontology)
  }
  ego_read_pt = pairwise_termsim(ego_read)
  treeplot(ego_read_pt, nCluster = n)
}

emapPlot = function (ego, ontology) {
  if (ontology != "ALL"){
    ego = filter(ego, ONTOLOGY == ontology)
  }
  ego = pairwise_termsim(ego)
  emapplot(ego, node_label = "group"
  )
}


# KEGG enrichment #
###################
enrich_KEGG = function (DEG) {
  genes = bitr(
    rownames(DEG),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db"
  )
  enrichKEGG(
    gene = genes$ENTREZID,
    organism = 'mmu',
    pvalueCutoff = 0.05
  )
}


# KEGG visualization #
######################
view_KEGG = function (pathway, resKEGG, DEG) {
  pathway = rownames(filter(resKEGG@result, Description == pathway))
  DEG$entrez = mapIds(
    org.Mm.eg.db,
    keys = row.names(DEG),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  x = dplyr::select(DEG, log2FoldChange, entrez)
  x = filter(x, is.na(entrez) == FALSE)
  genes = x$log2FoldChange
  names(genes) = x$entrez
  out = pathview(
    gene.data  = genes,
    pathway.id = pathway,
    species    = "mmu",
    same.layer = F, 
    limit = list(gene=max(abs(genes)), cpd=1)
  )

  img = readPNG(paste (pathway, ".pathview.png", sep =""))
  plot.new() 
  rasterImage(img,0,0,1,1)
  file.remove(paste(pathway, ".pathview.png", sep =""))
  file.remove(paste(pathway, ".png", sep =""))
  file.remove(paste(pathway, ".xml", sep =""))
}

#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "mmusculus_gene_ensembl",
#                         host = "http://www.ensembl.org")
#
# DEG$entrez <- getBM(filters = "ensembl_gene_id",
#                attributes = c("ensembl_gene_id","entrezgene_id"),
#                values = row.names(DEG),
#                mart = mart)
