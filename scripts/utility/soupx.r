

SoupCorrect <- function(raw_counts_path, filtered_counts_path, contamination_rate=NULL) {
  
  raw.matrix <- Read10X_h5(raw_counts_path)
  filt.matrix <- Read10X_h5(filtered_counts_path)

  object  <- CreateSeuratObject(counts=filt.matrix[['Gene Expression']])
  soup.channel  <- SoupX::SoupChannel(raw.matrix[['Gene Expression']], filt.matrix[['Gene Expression']])
  
  object <- object %>% 
    NormalizeData() %>%
    CleanVarGenes() %>% 
    ScaleData(verbose=FALSE) %>%
    RunPCA(verbose=FALSE) %>% 
    RunUMAP(dims=1:50, verbose=FALSE) %>% 
    FindNeighbors(dims=1:50, verbose=FALSE) %>% 
    FindClusters(verbose=FALSE, algorithm=4, resolution=0.8)
  
  meta <- object@meta.data
  umap <- object@reductions$umap@cell.embeddings
  soup.channel <- SoupX::setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- SoupX::setDR(soup.channel, umap)
  
  if (is.null(contamination_rate)) {
    soup.channel  <- SoupX::autoEstCont(soup.channel, forceAccept=TRUE, doPlot=FALSE)
  } else {
    soup.channel <- SoupX::setContaminationFraction(soup.channel, contamination_rate, forceAccept=TRUE)
  }

  adj.matrix  <- SoupX::adjustCounts(soup.channel)

  return(adj.matrix)
  
}

  