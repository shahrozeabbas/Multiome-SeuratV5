args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(tail(args, n=1)))
options(future.globals.maxSize=50 * 1000 * 1024^2)

assay <- 'ATAC'

input.files <- Filter(Negate(is.null), lapply(args, function(argument) {
    if (argument %like% paste0(assay, '_05.rds')) argument
}))


object.list <- future.apply::future_lapply(input.files, readRDS)

top.peaks <- Reduce(union, lapply(object.list, VariableFeatures))

message('merging samples...')
object <- merge(x=object.list[[1]], y=object.list[-1])
message('Done.')

batch <- 'sample'

DefaultAssay(object) <- assay

object <- object %>%

            RunSVD(features=top.peaks) %>% 
            harmony::RunHarmony(reduction='lsi', group.by.vars=batch, project.dim=FALSE, dims.use=2:50, reduction.save='harmony_atac')


saveRDS(object[['harmony_atac']], tail(args, n=2)[1])