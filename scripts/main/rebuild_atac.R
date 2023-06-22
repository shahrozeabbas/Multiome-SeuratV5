args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(args[4]))
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(args[2])
peaks.use <- readRDS(args[1])

annotations <- get_annotations()
    
project <- Project(object)
batch <- levels(object$batch)

fragments <- Fragments(object[['ATAC']])[[1]]
cells.use <- GetFragmentData(fragments, slot='cells')               
fragpath <- paste0(data_path, batch, '/Multiome/', project, '/outs/atac_fragments.tsv.gz')

atac_counts <- FeatureMatrix( 
    cells=cells.use,
    features=peaks.use,
    fragments=fragments,
    process_n=5000
)

object <- object %>% subset(cells=cells.use)

object[['ATAC']] <- as(object=CreateChromatinAssay(
    counts=BPCells::convert_matrix_type(atac_counts),
    fragments=fragpath,
    annotation=annotations
), Class='Assay5')


saveRDS(object, args[3])
