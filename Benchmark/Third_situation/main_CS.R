### import functions
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/methods_functions.R')
library(Seurat)
library(SeuratData)
data("pbmcsca")
data("panc8")

run_func <- function(exp_sc_mat, ref.mtx, ref.labels) {
    # pred.tags <- func_scMAGIC(exp_sc_mat, ref.mtx, ref.labels, atlas='HCL')
    # method <- 'scMAGIC'
    
    # pred.tags <- func_scmapcluster(exp_sc_mat, ref.mtx, ref.labels, threshold = 0.2)
    # method <- 'scmapcluster'
    
    # pred.tags <- func_scmapcell(exp_sc_mat, ref.mtx, ref.labels, w = 1, threshold = 0.4)
    # method <- 'scmapcell'
    
    # pred.tags <- func_scPred(exp_sc_mat, ref.mtx, ref.labels, nfeatures = 1000, npcs = 50)
    # method <- 'scPred'
    # pred.tags <- func_singleCellNet(exp_sc_mat, ref.mtx, ref.labels,
    #                                 nTopGenes = 10, nTopGenePairs = 25, nTrees = 500)
    # method <- 'singleCellNet'
    # pred.tags <- func_scClassify(exp_sc_mat, ref.mtx, ref.labels, tree = "HOPACH",
    #                              algorithm = "KNN", k = 10, topN = 50)
    # method <- 'scClassify'
    # pred.tags <- func_SingleR(exp_sc_mat, ref.mtx, ref.labels, gene_method='de', quantile.use=0.8)
    # method <- 'SingleR'
    
    pred.tags <- func_sciBet(exp_sc_mat, ref.mtx, ref.labels, k=600)
    method <- 'sciBet'
    
    # pred.tags <- func_CHETAH(exp_sc_mat, ref.mtx, ref.labels, clust_method = 'complete', n_genes = 0.8)
    # method <- 'CHETAH'
    
    return(list(pred.tags = pred.tags, method = method))
}

path_home <- '/mdshare/node9/zy/'
path_CS <- paste0(path_home, 'MAGIC/Benchmark/cross_species/')

################ Baron -> Segerstolpe #############
#########################################
file_panc_indrop <- '/mdshare/node9/zy/scRef/Benchmark/cross_species/BaronM.Rdata'
list.data <- readRDS(file_panc_indrop)
ref.mtx <- list.data$mat_exp
ref.labels <- list.data$label[,1]

list.data <- readRDS('/mdshare/node9/zy/scRef/Benchmark/cross_species/panc8_smartseq2.Rdata')
exp_sc_mat <- list.data$mat_exp
label_sc <- list.data$label[,1]
exp_sc_mat <- transformHomoloGene(exp_sc_mat)


list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CS, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Baron_Segerstolpe.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Baron_Segerstolpe.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Baron_Segerstolpe.Rdata')
saveRDS(ref.labels, file = file.ref)

#########################################

################ Baron -> Muraro #############
#########################################
file_panc_indrop <- '/mdshare/node9/zy/scRef/Benchmark/cross_species/BaronM.Rdata'
list.data <- readRDS(file_panc_indrop)
ref.mtx <- list.data$mat_exp
ref.labels <- list.data$label[,1]

list.data <- readRDS('/mdshare/node9/zy/scRef/Benchmark/cross_species/panc8_celseq2.Rdata')
exp_sc_mat <- list.data$mat_exp
label_sc <- list.data$label[,1]
exp_sc_mat <- transformHomoloGene(exp_sc_mat)

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CS, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Baron_Muraro.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Baron_Muraro.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Baron_Muraro.Rdata')
saveRDS(ref.labels, file = file.ref)

#########################################

################ MCA -> 10Xv2 #############
#########################################
list.data <- readRDS('/mdshare/node9/zy/scRef/Benchmark/cross_species/MCA1.Rdata')
ref.mtx <- list.data$mat_exp
ref.labels <- list.data$label[,1]

exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
exp_sc_mat <- transformHomoloGene(exp_sc_mat)

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CS, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_10Xv2.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_10Xv2.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_10Xv2.Rdata')
saveRDS(ref.labels, file = file.ref)

#########################################

################ MCA -> pbmc3k #############
#########################################
list.data <- readRDS('/mdshare/node9/zy/scRef/Benchmark/cross_species/MCA1.Rdata')
ref.mtx <- list.data$mat_exp
ref.labels <- list.data$label[,1]

list.data <- readRDS('/mdshare/node9/zy/scRef/Benchmark/cross_species/pbmc3k.Rdata')
exp_sc_mat <- list.data$mat_exp
label_sc <- list.data$label[,1]
exp_sc_mat <- transformHomoloGene(exp_sc_mat)

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_pbmc3k.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_pbmc3k.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_pbmc3k.Rdata')
saveRDS(ref.labels, file = file.ref)

#########################################




