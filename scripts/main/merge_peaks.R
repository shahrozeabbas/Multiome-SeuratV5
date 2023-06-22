args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(tail(args, 1)))
options(future.globals.maxSize=50 * 1000 * 1024^2)


peaks.list <- future.apply::future_lapply(head(args, -2), readRDS)

peaks.use <- suppressWarnings(reduce(x=unlist(as(peaks.list, 'GRangesList'))))

saveRDS(peaks.use, tail(args, 2)[1])