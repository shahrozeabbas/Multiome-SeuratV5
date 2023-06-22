args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(tail(args, n=1)))
options(future.globals.maxSize=50 * 1000 * 1024^2)


object.list <- future.apply::future_lapply(head(args, -3), readRDS)


m <- rbindlist(lapply(object.list, function(object) {
    m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
}))

cutoff <- NormalMixCutoff(mixtools::normalmixEM(m[, doublet_scores], k=2))

m[, `:=` (
    project=rep(tail(args, n=3)[1], nrow(m)),
    predicted_gmm_doublets=fifelse(doublet_scores < cutoff, 'singlet', 'doublet')
)]

fwrite(x=m, file=tail(args, n=2)[1])

