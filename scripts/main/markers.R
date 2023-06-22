args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)

object <- readRDS(snakemake@input[['seurat_object']])

Idents(object) <- 'wnn_clusters'

g <- data.table(genes=rownames(object[['RNA']]))
clean.genes <- g[!(genes %like% '^MT' | genes %like% '^RP[LS]'), genes]

markers <- object %>% FindAllMarkers(assay='RNA', min.pct=0.25, features=clean.genes)

fwrite(x=markers, file=snakemake@output[['markers']])


groups <- c('wnn_clusters', 'sample', 'batch')
titles <- c('Clusters', 'Sample', 'Batch')

plot.list <- lapply(seq_along(groups), function(counter) {
    object %>% 
        DimPlot(group.by=groups[counter], label=TRUE, repel=TRUE, reduction='wnn.umap') + 
            ggtitle(titles[counter]) + xlab('') + ylab('') + NoLegend()
})

g <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=2, nrow=2)

ggsave(plot=g, width=13, height=8, filename=snakemake@output[['umap']])
