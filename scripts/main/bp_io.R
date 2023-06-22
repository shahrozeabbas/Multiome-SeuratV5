working_dir <- '/data/CARD_singlecell/seurat_v5_multiome_dev/'; setwd(working_dir)

source('scripts/main/load_packages.r')

dataset <- snakemake@params[['sample']]
batch <- data.table::fread(snakemake@input[['samples']])[sample %chin% dataset, batch]


raw.counts <- paste0(data_path, batch, '/Multiome/', dataset, '/outs/raw_feature_bc_matrix.h5')
filtered.counts <- paste0(data_path, batch, '/Multiome/', dataset, '/outs/filtered_feature_bc_matrix.h5')

adj.matrix <- suppressWarnings(SoupCorrect(raw.counts, filtered.counts, contamination_rate=snakemake@params[['soup_rate']]))

# raw.counts <- BPCells::open_matrix_10x_hdf5(raw.counts, feature_type='Gene Expression')
# filtered.counts <- BPCells::open_matrix_10x_hdf5(filtered.counts, feature_type='Gene Expression')

object <- CreateSeuratObject(convert_matrix_type(adj.matrix), min.cells=0, min.features=0, project=dataset)



