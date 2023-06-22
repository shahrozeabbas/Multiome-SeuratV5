args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')


object <- readRDS(args[1])

hub <- AnnotationHub::AnnotationHub()
reference <- hub[['AH75011']]

project <- Project(object)
batch <- levels(object$batch)

frag.file <- paste0(data_path, batch, '/Multiome/', project, '/outs/atac_fragments.tsv.gz')
counts.file <- paste0(data_path, batch, '/Multiome/', project, '/outs/filtered_feature_bc_matrix.h5')

counts <- BPCells::open_matrix_10x_hdf5(counts.file, feature_type='Peaks')


object <- object %>% 

    AddChromiumAssay(atac_counts=counts, enDB=reference, frag.file=frag.file) %>% 
        
    NucleosomeSignal(assay='ATAC') %>% TSSEnrichment(assay='ATAC')


saveRDS(object, args[2])
