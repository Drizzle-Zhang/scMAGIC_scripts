### import functions
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/methods_functions.R')
library(Seurat)
library(SeuratData)
data("pbmcsca")
data("panc8")

run_func <- function(exp_sc_mat, ref.mtx, ref.labels) {
    # pred.tags <- func_scMAGIC(exp_sc_mat, ref.mtx, ref.labels, altas = 'HCL',
    #                           n_gene=2000, conf_cutoff=5, num_threads = 8)
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
    
    # pred.tags <- func_sciBet(exp_sc_mat, ref.mtx, ref.labels, k=600)
    # method <- 'sciBet'
    
    # pred.tags <- func_CHETAH(exp_sc_mat, ref.mtx, ref.labels, clust_method = 'complete', n_genes = 0.8)
    # method <- 'CHETAH'
    
    return(list(pred.tags = pred.tags, method = method))
}

path_home <- '/mdshare/node9/zy/'
path_CP <- paste0(path_home, 'MAGIC/Benchmark/cross_platform/')

################ Baron -> Muraro #############
#########################################
file_panc_indrop <- '/mdshare/node9/zy/MAGIC/sc_data/CV/panc_indrop.Rdata'
list.data <- readRDS(file_panc_indrop)
ref.mtx <- list.data$mat_exp
ref.labels <- list.data$label[,1]

exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
label_sc <- as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')]


list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Baron_Muraro.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Baron_Muraro.Rdata')
saveRDS(true.tags, file = file.true)

#########################################

################ Baron -> Xin #############
#########################################
file_panc_indrop <- '/mdshare/node9/zy/MAGIC/sc_data/CV/panc_indrop.Rdata'
list.data <- readRDS(file_panc_indrop)
ref.mtx <- list.data$mat_exp
ref.labels <- list.data$label[,1]

list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/Xin.Rdata')
exp_sc_mat <- list.data$mat_exp
label_sc <- list.data$label[,1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Baron_Xin.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Baron_Xin.Rdata')
saveRDS(true.tags, file = file.true)

#########################################

################ Muraro -> Baron #############
#########################################
ref.mtx <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
ref.labels <- as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')]

file_panc_indrop <- '/mdshare/node9/zy/MAGIC/sc_data/CV/panc_indrop.Rdata'
list.data <- readRDS(file_panc_indrop)
exp_sc_mat <- list.data$mat_exp
label_sc <- list.data$label[,1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Muraro_Baron.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Muraro_Baron.Rdata')
saveRDS(true.tags, file = file.true)

#########################################

################ Muraro -> Xin #############
#########################################
ref.mtx <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
ref.labels <- as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')]

list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/Xin.Rdata')
exp_sc_mat <- list.data$mat_exp
label_sc <- list.data$label[,1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Muraro_Xin.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Muraro_Xin.Rdata')
saveRDS(true.tags, file = file.true)

#########################################

################ Segerstolpe -> Muraro #############
#########################################
ref.mtx <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('smartseq2')])
ref.labels <- as.character(panc8$celltype)[panc8$dataset %in% c('smartseq2')]

exp_sc_mat <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
label_sc <- as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Segerstolpe_Muraro.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Segerstolpe_Muraro.Rdata')
saveRDS(true.tags, file = file.true)

#########################################

################ Campbell -> Tasic #############
#########################################
OUT <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/mouse_brain/Campbell.Rdata')
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label$CellType

OUT <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/mouse_brain/Tasic.Rdata')
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[,1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Campbell_Tasic.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Campbell_Tasic.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Campbell_Tasic.Rdata')
saveRDS(ref.labels, file = file.ref)

#########################################

################ Campbell -> Zeisel #############
#########################################
OUT <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/mouse_brain/Campbell.Rdata')
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label$CellType

OUT <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/mouse_brain/Zeisel.Rdata')
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[,1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_CP, method, '/other_CP/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Campbell_Zeisel.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Campbell_Zeisel.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Campbell_Zeisel.Rdata')
saveRDS(ref.labels, file = file.ref)

#########################################




