### import functions
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/methods_functions.R')

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
    # pred.tags <- func_sciBet(exp_sc_mat, ref.mtx, ref.labels, k=600)
    # method <- 'sciBet'
    pred.tags <- func_CHETAH(exp_sc_mat, ref.mtx, ref.labels, clust_method = 'complete', n_genes = 0.8)
    method <- 'CHETAH'
    
    return(list(pred.tags = pred.tags, method = method))
}

path_home <- '/mdshare/node9/zy/'
path_MB <- paste0(path_home, 'MAGIC/Benchmark/mouse_brain/')

######### Tasic -> Tasic2018VISp #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Tasic2.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[, 1]

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/VISp_Tasic.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Tasic_Tasic2018VISp.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Tasic_Tasic2018VISp.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Tasic_Tasic2018VISp.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### Tasic -> Tasic2018ALM #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Tasic2.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[, 1]

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/ALM_Tasic.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Tasic_Tasic2018ALM.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Tasic_Tasic2018ALM.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Tasic_Tasic2018ALM.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### MCA -> Tasic2018VISp #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/MCA.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/VISp_Tasic.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_Tasic2018VISp.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_Tasic2018VISp.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_Tasic2018VISp.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### MCA -> Tasic2018ALM #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/MCA.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/ALM_Tasic.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_Tasic2018ALM.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_Tasic2018ALM.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_Tasic2018ALM.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################


######### Tasic -> Campbell #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Tasic.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[, 1]

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Campbell.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Tasic_Campbell.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Tasic_Campbell.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Tasic_Campbell.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### Tasic -> Hochgerner #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Tasic.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[, 1]

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/HochgernerA.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 2]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Tasic_Hochgerner.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Tasic_Hochgerner.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Tasic_Hochgerner.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### Tasic -> Mizrak #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Tasic.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label[, 1]

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Mizrak.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_Tasic_Mizrak.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_Tasic_Mizrak.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_Tasic_Mizrak.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### MCA -> Campbell #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/MCA.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Campbell.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_Campbell.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_Campbell.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_Campbell.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### MCA -> Hochgerner #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/MCA.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/HochgernerA.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 2]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_Hochgerner.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_Hochgerner.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_Hochgerner.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################

######### MCA -> Mizrak #############
###############################################
OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/MCA.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <- OUT$label

OUT <- readRDS(paste0(path_home, 'scMAGIC/Benchmark/mouse_brain/Mizrak.Rdata'))
exp_sc_mat <- OUT$mat_exp
label_sc <- OUT$label[, 1]

list_func <- run_func(exp_sc_mat, ref.mtx, ref.labels)
pred.tags <- list_func$pred.tags
method <- list_func$method

path_data <- paste0(path_MB, method, '/')
if (!file.exists(path_data)) {
    dir.create(path_data)
}
saveRDS(pred.tags, file = paste0(path_data, 'pred_MCA_Mizrak.Rdata'))
true.tags <- label_sc
file.true <- paste0(path_data, 'true_MCA_Mizrak.Rdata')
saveRDS(true.tags, file = file.true)
file.ref <- paste0(path_data, 'ref_MCA_Mizrak.Rdata')
saveRDS(ref.labels, file = file.ref)

###############################################
