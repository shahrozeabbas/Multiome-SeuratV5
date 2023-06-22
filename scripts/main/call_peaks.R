
args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')


object <- readRDS(args[2])
        
object <- object %>% subset(
    
    subset=
            nCount_RNA > 300 & nCount_ATAC > 300 &
            nFeature_RNA > 300 & nFeature_ATAC > 300 &
            nucleosome_signal < 2 & TSS.enrichment > 2 &
            percent.mt < 10,
    cells=
            fread(args[1])[
                sample %chin% Project(object) & 
                predicted_gmm_doublets %chin% 'singlet', cells
            ]  
    )

saveRDS(object, args[4])

macs_path <- '/data/abbass2/mambaforge/envs/macs/bin/macs3'

peaks <- object %>% 

            CallPeaks(assay='ATAC', macs2.path=macs_path) %>% 
            GenomeInfoDb::keepStandardChromosomes(pruning.mode='coarse') %>% 
            IRanges::subsetByOverlaps(ranges=blacklist_hg38_unified, invert=TRUE)


saveRDS(peaks, args[3])

