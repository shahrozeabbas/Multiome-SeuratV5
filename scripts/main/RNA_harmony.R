args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(tail(args, n=1)))
options(future.globals.maxSize=50 * 1000 * 1024^2)

assay <- 'RNA'

input.files <- Filter(Negate(is.null), lapply(args, function(argument) {
    if (argument %like% paste0(assay, '_05.rds')) argument
}))


object.list <- future.apply::future_lapply(input.files, readRDS)

top.genes <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

message('merging samples...')
object <- merge(x=object.list[[1]], y=object.list[-1])
message('Done.')

batch <- 'sample'; noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores')

DefaultAssay(object) <- assay

object <- object %>% 

            ScaleData(vars.to.regress=noise, features=top.genes, verbose=FALSE) %>% 
            RunPCA(features=top.genes, npcs=50, verbose=FALSE) %>% 
            
            harmony::RunHarmony(reduction='pca', group.by.vars=batch, max.iter.harmony=50, reduction.save='harmony_rna')


saveRDS(object[['harmony_rna']], tail(args, n=2)[1])