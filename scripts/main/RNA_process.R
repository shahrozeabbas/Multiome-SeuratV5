args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(args[3]))
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(args[1])


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

    
DefaultAssay(object) <- 'RNA'

all.genes <- rownames(object)

object <- object %>% 

    NormalizeData() %>% 
    CleanVarGenes(nHVG=2000) %>% 
    ScaleData(features=all.genes, verbose=FALSE)
    # CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes, assay='RNA')


saveRDS(object, args[2])