

PKGS <- c('Seurat', 'Signac', 'dplyr', 'ggplot2', 'data.table', 'BPCells')

invisible(sapply(PKGS, require, character.only=TRUE))
invisible(lapply(list.files('scripts/utility', full.names=TRUE, pattern='\\.r$'), source))


data_path <- '/data/CARD_singlecell/Brain_atlas/NABEC_multiome/'

setDTthreads(threads=1)


options(Seurat.object.assay.version='v5')
options(future.globals.maxSize=100 * 1000 * 1024^2)