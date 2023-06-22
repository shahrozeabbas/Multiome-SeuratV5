args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(args[3]))
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(args[1])


DefaultAssay(object) <- 'ATAC'

object <- object %>% suppressWarnings(RunTFIDF()) %>% FindTopFeatures(min.cutoff='q25')


saveRDS(object, args[2])

