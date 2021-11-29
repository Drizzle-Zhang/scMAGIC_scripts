### import functions
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/Cross_validation/Cross_Validation.R')
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/Cross_validation/method_functions.R')
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/methods_functions.R')

library(Seurat)
library(SeuratData)
data("pbmcsca")
data("panc8")

run_all_methods <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir) {
    # run_scMAGIC(DataPath,LabelsPath,CV_RDataPath,OutputDir,
    #             atlas='HCL',NumGenes=2000,conf_cutoff=5)
    # run_sciBet(DataPath,LabelsPath,CV_RDataPath,OutputDir,k=600)
    # run_singleCellNet(DataPath,LabelsPath,CV_RDataPath,OutputDir,
    #                   nTopGenes = 10, nTopGenePairs = 25, nTrees = 500)
    # run_scPred(DataPath,LabelsPath,CV_RDataPath,OutputDir,
    #            nfeatures = 1000, npcs = 50)
    # run_scmapcluster(DataPath,LabelsPath,CV_RDataPath,OutputDir,
    #                  threshold = 0.2)
    # run_SingleR(DataPath,LabelsPath,CV_RDataPath,OutputDir,
    #             genes = "de", quantile.use = 0.8)
    # run_scClassify(DataPath,LabelsPath,CV_RDataPath,OutputDir, tree = "HOPACH", 
    #                algorithm = "KNN", k = 10, topN = 50)
    # run_CHETAH(DataPath,LabelsPath,CV_RDataPath,OutputDir,  
    #            clust_method = 'complete', n_genes = 1000)
    # run_scmapcell(DataPath,LabelsPath,CV_RDataPath,OutputDir,
    #               w = 1, threshold = 0.4)
    
    
}

pathout_CV <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_validation/'


############ Baron_mouse #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/cross_species/BaronM.Rdata')
mat_panc <- list.data$mat_exp
label_panc <- list.data$label[,1]

OutputDir <- paste0(pathout_CV, 'Baron_mouse/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)
###########################################

############ Baron_human #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/panc_indrop.Rdata')
mat_panc <- list.data$mat_exp
label_panc <- list.data$label[,1]

OutputDir <- paste0(pathout_CV, 'Baron_human/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Muraro #############
###########################################
mat_panc <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('celseq2')])
label_panc <- data.frame(
    annotations = as.character(panc8$celltype)[panc8$dataset %in% c('celseq2')],
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'Muraro/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Segerstolpe #############
###########################################
mat_panc <- as.matrix(panc8@assays$RNA@counts[, panc8$dataset %in% c('smartseq2')])
label_panc <- data.frame(
    annotations = as.character(panc8$celltype)[panc8$dataset %in% c('smartseq2')],
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'Segerstolpe/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Xin #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/Xin_human_pancreatic_SMARTer.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata$cell.type,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'Xin/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ CellBench_10X #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/CellBench_human_lung_10Xchromium.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata$cell_line_demuxlet,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'CellBench_10X/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ CellBench_CELSeq2 #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/CellBench_human_lung_CELseq2.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata$cell_line_demuxlet,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'CellBench_CELSeq2/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ TM #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/TM_mouse_20tissues_SMARTseq2.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata,
    row.names = colnames(list.data$mat))

OutputDir <- paste0(pathout_CV, 'TM/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ AMB_Level1 #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/Tasic2018_AMB_1.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata$cell_class,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'AMB_Level1/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ AMB_Level2 #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/Tasic2018_AMB_2.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata$cell_class,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'AMB_Level2/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ AMB_Level3 #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/Tasic2018_AMB_3.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata$cell_class,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'AMB_Level3/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Zhang_sorted #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/zheng_sorted_pbmc_10Xchromium.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'Zhang_sorted/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Zhang_68K #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/MAGIC/sc_data/CV/zheng_68K_pbmc_10Xchromium.rds')
mat_panc <- list.data$mat
label_panc <- data.frame(
    annotations = list.data$metadata,
    row.names = colnames(mat_panc))

OutputDir <- paste0(pathout_CV, 'Zhang_68K/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Tasic #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/mouse_brain/Tasic.Rdata')
mat_panc <- list.data$mat_exp
label_panc <- list.data$label[,1]

OutputDir <- paste0(pathout_CV, 'Tasic/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################

############ Campbell #############
###########################################
list.data <- readRDS('/mdshare/node9/zy/scMAGIC/Benchmark/mouse_brain/Campbell.Rdata')
mat_panc <- list.data$mat_exp
label_panc <- list.data$label[,1]

OutputDir <- paste0(pathout_CV, 'Campbell/')
if (!file.exists(OutputDir)) {
    dir.create(OutputDir)
}
LabelsPath <- paste0(OutputDir, 'label.txt')
write.table(label_panc, file = LabelsPath, sep = '\t', quote = F)
Cross_Validation(LabelsPath, OutputDir)
CV_RDataPath <- paste0(OutputDir, 'CV_folds.RData')
OUT <- list()
OUT$data.filter <- mat_panc
OUT$label.filter <- data.frame(annotations = label_panc,
                               row.names = colnames(OUT$data.filter))
DataPath <- paste0(OutputDir, 'Data.Rdata')
saveRDS(OUT, file = DataPath)

run_all_methods(DataPath,LabelsPath,CV_RDataPath,OutputDir)
print(OutputDir)

###########################################


