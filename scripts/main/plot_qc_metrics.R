args <- commandArgs(trailingOnly=TRUE)

working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=as.numeric(args[5]))
options(future.globals.maxSize=50 * 1000 * 1024^2)


m <- fread(args[1])

noise <- c('nCount_ATAC', 'nCount_RNA', 'nFeature_ATAC', 'nFeature_RNA', 'percent.mt', 'percent.rb', 'doublet_scores', 'nucleosome_signal', 'TSS.enrichment')


plot.list <- future.apply::future_lapply(noise, function(feature) {
    
    m %>%

        ggplot(aes(x=args[4], y=get(feature))) + 
        
        geom_violin(fill='steelblue', color='black') + theme_bw() + 
        
        theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + ylab(feature)
        
})

p <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=3)

ggsave(plot=p, width=15, height=8, filename=args[2])


plot1 <- m %>% ggplot(aes(x=get(noise[2]), y=get(noise[4]))) + geom_point(alpha=0.3) 
    theme_bw() + xlab('Number of UMIs - RNA') + ylab('Number of Genes - RNA')

plot2 <- m %>% ggplot(aes(x=get(noise[1]), y=get(noise[3]))) + geom_point(alpha=0.3) 
    theme_bw() + xlab('Number of UMIs - ATAC') + ylab('Number of Genes - ATAC')


ggsave(plot=plot1 + plot2, width=18, height=9, filename=args[3])







