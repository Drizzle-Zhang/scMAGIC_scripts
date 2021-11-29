library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')
library(Seurat)
library(SeuratData)
data("pbmcsca")

# evaluation
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/Cross_validation/evaluate.R')
simple.evaluation <- function(true.tag, scRef.tag, df.ref.names, df.sc.names) {
    # uniform tags
    for (j in 1:dim(df.ref.names)[1]) {
        scRef.tag[scRef.tag == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
    }
    scRef.tag[!(scRef.tag %in% df.ref.names$name)] <- 'Unassigned'
    for (j in 1:dim(df.sc.names)[1]) {
        true.tag[true.tag == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
    }
    
    percent.unassigned <- sum(scRef.tag == 'Unassigned')/sum(true.tag == 'Unassigned')
    
    
    # true.labels <- setdiff(unique(true.tag), 'Unassigned')
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    balanced_acc <- metrics$balanced_accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    macro_f1.rm.unassigned <- metrics$f1_score(true.tag.rm, our.tag.rm, average = 'macro', labels = unique(true.tag.rm))
    balanced.accuracy.rm.unassigned <-
        metrics$balanced_accuracy_score(true.tag.rm, our.tag.rm)
    
    f1 <- c()
    for (label in true.labels) {
        tmp.true.tag <- true.tag
        tmp.our.tag <- our.tag
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.f1 <- metrics$f1_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        f1 <- c(f1, sub.f1)
    }
    names(f1) <- true.labels
    
    our.labels <- setdiff(unique(our.tag), 'Unassigned')
    precision <- c()
    for (label in our.labels) {
        tmp.true.tag <- true.tag.rm
        tmp.our.tag <- our.tag.rm
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.precision <- metrics$precision_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        precision <- c(precision, sub.precision)
        
    }
    names(precision) <- our.labels
    mean.precision <- mean(precision)
    
    out <- list()
    out$percent.unassigned <- percent.unassigned
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$balanced.accuracy <- balanced_acc
    out$f1 <- f1
    out$med_f1 <- median(f1)
    out$accuracy.rm.unassigned <- accuracy.rm.unassigned
    out$macro_f1.rm.unassigned <- macro_f1.rm.unassigned
    out$precision.rm.unassigned <- precision
    out$balanced.accuracy.rm.unassigned <- balanced.accuracy.rm.unassigned
    out$mean.precision.rm.unassigned <- mean.precision
    out$conf <- table(true.tag, our.tag)
    
    return(out)
    
}

uniform_labels <- function(true.tag, scRef.tag, df.ref.names, df.sc.names) {
    # uniform tags
    for (j in 1:dim(df.ref.names)[1]) {
        scRef.tag[scRef.tag == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
    }
    scRef.tag[!(scRef.tag %in% df.ref.names$name)] <- 'Unassigned'
    for (j in 1:dim(df.sc.names)[1]) {
        true.tag[true.tag == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
    }
    return(list(ref_labels = scRef.tag, label_sc = true.tag))
}

cluster_tags <- function(tags, clusters) {
    tag_table <- table(tags, clusters)
    cluster_query <- c()
    if (dim(tag_table)[1] == 1) {
        cluster_query <- rep(rownames(tag_table), dim(tag_table)[2])
    } else {
        for (clust in colnames(tag_table)) {
            sub_col <- tag_table[, clust]
            cluster_query <- c(cluster_query, names(sub_col)[sub_col == max(sub_col)][1])
        }
    }
    
    return(cluster_query)
}

func_eval_cluster <- function(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster) {
    list_labels <- uniform_labels(true.tags, pred.tags, df.ref.names, df.sc.names)
    pred_labels <- list_labels$ref_labels
    query_label <- list_labels$label_sc
    clusters <- df_cluster$cluster
    cluster_query <- cluster_tags(query_label, clusters)
    cluster_pred <- cluster_tags(pred_labels, clusters)
    num = 0
    for (i in 1:length(cluster_pred)) {
        if (cluster_query[i] == cluster_pred[i]) {
            num = num + 1
        }
    }
    return(num)
}

# path
path_CV <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_validation/'
path_supp <- '/mdshare/node9/zy/MAGIC/Benchmark/supplementary/'
path_CP <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/'
path_eval_dataset <- '/mdshare/node9/zy/MAGIC/evaluation/datasets/'
path_cluster <- '/mdshare/node9/zy/MAGIC/MAGIC_cluster/'

# methods
vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')

# dataset
vec.folders <- c('Baron_mouse', 'Baron_human', 'Muraro', 'Segerstolpe', 'Xin',
                 'CellBench_10X', 'CellBench_CELSeq2', 'Tasic', 'Campbell',
                 'TM', 'AMB_Level1', 'AMB_Level2', 'AMB_Level3', 'Zhang_sorted', 'Zhang_68K', 'Ding')
vec.datasets <- c('Baron, Mouse pancreas', 'Baron, Human pancreas', 'Muraro, Human pancreas',
                  'Segerstolpe, Human pancreas', 'Xin, Human pancreas', 'Tian, CellBench 10X',
                  'Tian, CellBench CEL-Seq2', 'Tasic2015, Mouse brain', 'Campbell, Mouse brain', 
                  'The Tabula Muris Consortium', 'Tasic2018, Mouse brain(Level 1)',
                  'Tasic2018, Mouse brain(Level 2)', 'Tasic2018, Mouse brain(Level 3)',
                  'Zhang, Human PBMC(sorted)', 'Zhang, Human PBMC(68K)', 'Ding, Human PBMC')

# CV
df_CV <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    for (i in 1:length(vec.folders)) {
        folder <- vec.folders[i]
        dataset <- vec.datasets[i]
        pathout.dataset <- paste0(path_CV, folder, '/')
        TrueLabelsPath <- paste0(pathout.dataset, method, '_True_Labels.csv')
        PredLabelsPath <- paste0(pathout.dataset, method, '_Pred_Labels.csv')
        res.out <- evaluate(TrueLabelsPath, PredLabelsPath)
        df.sub <- data.frame(term = 'Accuracy', method = method,
                             dataset = dataset,
                             value = res.out$accuracy, stringsAsFactors = F)
        df.sub <- rbind(df.sub, 
                          data.frame(term = 'Balanced accuracy', method = method,
                                     dataset = dataset,
                                     value = res.out$balanced.accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Median F1', method = method,
                                   dataset = dataset,
                                   value = res.out$MedF1, stringsAsFactors = F))
        # cluster
        df_cluster <- readRDS(paste0(path_cluster, folder, '.Rdata'))
        df_cluster <- df_cluster[colnames(exp_sc_mat), ]
        cell_id <- colnames(exp_sc_mat)[true.tags_ori%in%overlap_cells & pred.tags_ori%in%overlap_cells]
        df_cluster <- df_cluster[intersect(df_cluster$id, cell_id),]
        all_cluster <- length(table(df_cluster$cluster))
        correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Number of correctly labeled clusters', method = method,
                                   dataset = dataset,
                                   value = correct_cluster, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Total number of clusters', method = method,
                                   dataset = dataset,
                                   value = all_cluster, stringsAsFactors = F))
        
        # dataset
        list_dataset <- readRDS(paste0(path_eval_dataset, folder, '.Rdata'))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'F1_asw', method = method,
                                   dataset = dataset,
                                   value = list_dataset$F1_asw, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Depth of reference', method = method,
                                   dataset = dataset,
                                   value = list_dataset$depth_ref, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Cell proportion covered by reference', method = method,
                                   dataset = dataset,
                                   value = list_dataset$prop_cover, stringsAsFactors = F))

        df_CV <- rbind(df_CV, df.sub)
        print(c(method, folder))
    }
}
df_CV$type <- rep('Cross-validation', nrow(df_CV))
df_CV$type[df_CV$dataset == 'Ding, Human PBMC'] <- 'Cross-sample'

df_crosssample <- data.frame(stringsAsFactors = F)
vec.platform <- c('CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
for (method in vec.methods) {
    path.supp <- paste0(path_supp, method, '/pbmcsca/')
    for (platform in vec.platform) {
        platform_name <- platform
        ref <- 'pbmc1'
        query <- 'pbmc2'
        if (platform == '10Xv2') {
            platform <- '10x Chromium (v2)'
            ref.labels1 <- as.character(pbmcsca$CellType)[(pbmcsca$Method == '10x Chromium (v2) A') & (pbmcsca$Experiment == 'pbmc1')]
            ref.labels2 <- as.character(pbmcsca$CellType)[(pbmcsca$Method == '10x Chromium (v2) B') & (pbmcsca$Experiment == 'pbmc1')]
            ref.labels = c(ref.labels1, ref.labels2)
        } else {
            if (platform == '10Xv3') {
                platform <- '10x Chromium (v3)'
            } else {
                if (platform == 'inDrops') {
                    ref <- 'pbmc2'
                    query <- 'pbmc1'
                }
            }
            ref.labels <- as.character(pbmcsca$CellType)[(pbmcsca$Method == platform) & (pbmcsca$Experiment == ref)]
        }
        
        
        exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, (pbmcsca$Method == platform) & (pbmcsca$Experiment == query)])
        label_sc <- as.character(pbmcsca$CellType)[(pbmcsca$Method == platform) & (pbmcsca$Experiment == query)]
        
        overlap_cells <- intersect(unique(ref.labels), unique(label_sc))
        ref.labels <- ref.labels[ref.labels %in% overlap_cells]
        label_sc <- label_sc[label_sc %in% overlap_cells]
        ref.names <- names(table(ref.labels))
        df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
        all.cell <- names(table(label_sc))
        df.sc.names <- data.frame(sc.name = all.cell, name = all.cell)
        
        true.tags_ori <- readRDS(paste0(path.supp, 'true_', platform_name, '_', platform_name, '.Rdata'))
        if (method == 'CALLR') {
            pred.tags_list <- readRDS(paste0(path.supp, 'pred_', platform_name, '_', platform_name, '.Rdata'))
            true.tags_ori <- true.tags_ori[pred.tags_list[[2]]]
            pred.tags_ori <- pred.tags_list[[1]]
        } else {
            pred.tags_ori <- readRDS(paste0(path.supp, 'pred_', platform_name, '_', platform_name, '.Rdata'))
        }
        pred.tags <- pred.tags_ori[true.tags_ori%in%overlap_cells & pred.tags_ori%in%overlap_cells]
        true.tags <- true.tags_ori[true.tags_ori%in%overlap_cells & pred.tags_ori%in%overlap_cells]
        res.out <- simple.evaluation(true.tags, pred.tags, df.ref.names, df.sc.names)
        df.sub <- data.frame(term = 'Weighted macro F1', method = method,
                             value = res.out$weighted_macro_f1, stringsAsFactors = F)
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Macro F1', method = method,
                                   value = res.out$macro_f1, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Accuracy', method = method,
                                   value = res.out$accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Accuracy', method = method,
                                   value = res.out$accuracy.rm.unassigned, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
        # cluster
        df_cluster <- readRDS(paste0(path_cluster, platform_name, '_', platform_name, '.Rdata'))
        df_cluster <- df_cluster[colnames(exp_sc_mat), ]
        cell_id <- colnames(exp_sc_mat)[true.tags_ori%in%overlap_cells & pred.tags_ori%in%overlap_cells]
        df_cluster <- df_cluster[intersect(df_cluster$id, cell_id),]
        all_cluster <- length(table(df_cluster$cluster))
        correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Number of correctly labeled clusters', method = method,
                                   value = correct_cluster, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Total number of clusters', method = method,
                                   value = all_cluster, stringsAsFactors = F))
        
        # dataset
        list_dataset <- readRDS(paste0(path_eval_dataset, platform_name, '_', platform_name, '.Rdata'))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'F1_asw', method = method,
                                   value = list_dataset$F1_asw, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Depth of reference', method = method,
                                   value = list_dataset$depth_ref, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Cell proportion covered by reference', method = method,
                                   value = list_dataset$prop_cover, stringsAsFactors = F))
        
        if (platform_name == '10Xv2') {
            platform_name <- '10x Chromium (v2)'
        } else {
            if (platform_name == '10Xv3') {
                platform_name <- '10x Chromium (v3)'
            } else {
                platform_name <- platform_name
            }
        }
        df.sub$dataset <- rep(paste0(platform_name, ' -> ', platform_name, ', Human PBMC'), nrow(df.sub))
        df_crosssample <- rbind(df_crosssample, df.sub)
        print(c(method, platform_name))
    }
    
}

vec.ref <- c('Tasic2018VISp', 'Tasic2018ALM')
vec.query <- c('Tasic2018ALM', 'Tasic2018VISp')
vec.dataset <- c('Tasic2018, VISp -> ALM', 'Tasic2018, ALM -> VISp')
for (method in vec.methods) {
    path.supp <- paste0(path_supp, method, '/mouse_brain/')
    for (i in 1:length(vec.ref)) {
        ref <- vec.ref[i]
        query <- vec.query[i]
        true.tags <- readRDS(paste0(path.supp, 'true_', ref, '_', query, '.Rdata'))
        pred.tags <- readRDS(paste0(path.supp, 'pred_', ref, '_', query, '.Rdata'))
        ref.names <- names(table(pred.tags))
        all.cell <- names(table(true.tags))
        df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
        df.sc.names <- data.frame(sc.name = all.cell, name = all.cell)
        res.out <- simple.evaluation(true.tags, pred.tags, df.ref.names, df.sc.names)
        df.sub <- data.frame(term = 'Weighted macro F1', method = method,
                             value = res.out$weighted_macro_f1, stringsAsFactors = F)
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Macro F1', method = method,
                                   value = res.out$macro_f1, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Accuracy', method = method,
                                   value = res.out$accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Accuracy', method = method,
                                   value = res.out$accuracy.rm.unassigned, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
        # cluster
        df_cluster <- readRDS(paste0(path_cluster, platform_name, '.Rdata'))
        all_cluster <- length(table(df_cluster$cluster))
        correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Number of correctly labeled clusters', method = method,
                                   value = correct_cluster, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Total number of clusters', method = method,
                                   value = all_cluster, stringsAsFactors = F))
        # dataset
        list_dataset <- readRDS(paste0(path_eval_dataset, ref, '_', query, '.Rdata'))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'F1_asw', method = method,
                                   value = list_dataset$F1_asw, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Depth of reference', method = method,
                                   value = list_dataset$depth_ref, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Cell proportion covered by reference', method = method,
                                   value = list_dataset$prop_cover, stringsAsFactors = F))
        dataset <- vec.dataset[i]
        df.sub$dataset <- rep(dataset, nrow(df.sub))
        df_crosssample <- rbind(df_crosssample, df.sub)
        print(c(method, ref, query))
    }
}
df_crosssample$type <- rep('Cross-sample', nrow(df_crosssample))


# other CP
vec.ref <- c('Baron', 'Baron', 'Muraro', 'Muraro', 'Segerstolpe')
vec.query <- c('Muraro', 'Xin', 'Baron', 'Xin', 'Muraro')
vec.dataset <- c('inDrops -> CEL-Seq2, Human Pancreas', 
                 'inDrops -> SMARTer, Human Pancreas', 'CEL-Seq2 -> inDrops, Human Pancreas',
                 'CEL-Seq2 -> SMARTer, Human Pancreas', 'Smart-seq2 -> CEL-Seq2, Human Pancreas')

df_CP <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path.CP <- paste0(path_CP, method, '/other_CP/')
    for (i in 1:length(vec.ref)) {
        ref <- vec.ref[i]
        query <- vec.query[i]
        true.tags <- readRDS(paste0(path.CP, 'true_', ref, '_', query, '.Rdata'))
        pred.tags <- readRDS(paste0(path.CP, 'pred_', ref, '_', query, '.Rdata'))
        ref.names <- names(table(pred.tags))
        all.cell <- names(table(true.tags))
        df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
        df.sc.names <- data.frame(sc.name = all.cell, name = all.cell)
        res.out <- simple.evaluation(true.tags, pred.tags, df.ref.names, df.sc.names)
        df.sub <- data.frame(term = 'Weighted macro F1', method = method,
                             value = res.out$weighted_macro_f1, stringsAsFactors = F)
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Macro F1', method = method,
                                   value = res.out$macro_f1, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Accuracy', method = method,
                                   value = res.out$accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Accuracy', method = method,
                                   value = res.out$accuracy.rm.unassigned, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Median F1', method = method,
                                   value = res.out$med_f1, stringsAsFactors = F))
        # cluster
        df_cluster <- readRDS(paste0(path_cluster, ref, '_', query, '.Rdata'))
        all_cluster <- length(table(df_cluster$cluster))
        correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Number of correctly labeled clusters', method = method,
                                   value = correct_cluster, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Total number of clusters', method = method,
                                   value = all_cluster, stringsAsFactors = F))
        
        # dataset
        list_dataset <- readRDS(paste0(path_eval_dataset, ref, '_', query, '.Rdata'))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'F1_asw', method = method,
                                   value = list_dataset$F1_asw, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Depth of reference', method = method,
                                   value = list_dataset$depth_ref, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Cell proportion covered by reference', method = method,
                                   value = list_dataset$prop_cover, stringsAsFactors = F))
        dataset <- vec.dataset[i]
        df.sub$dataset <- rep(dataset, nrow(df.sub))
        df_CP <- rbind(df_CP, df.sub)
        print(c(method, ref, query))
    }
    
    # Campbell -> Tasic
    file.pred <- paste0(path.CP, 'pred_Campbell_Tasic.Rdata')
    file.true <- paste0(path.CP, 'true_Campbell_Tasic.Rdata')
    file.ref <- paste0(path.CP, 'ref_Campbell_Tasic.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    ref.names <- names(table(ref.labels))
    uniform.names <- c("Astrocytes", "Endothelial cells", "Ependymocytes", "Mural cells", "Neurons",
                       "Oligodendrocytes", "OPC", "Pars tuberalis", "PVMs & Microglia", "Tanycytes", "VLMCs")
    df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocytes", "Endothelial cells",
                       "PVMs & Microglia", "Neurons", "Oligodendrocytes", "OPC")
    df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)
    res.out <- simple.evaluation(true.tags, pred.tags, df.ref.names, df.sc.names)
    df.sub <- data.frame(term = 'Weighted macro F1', method = method,
                         value = res.out$weighted_macro_f1, stringsAsFactors = F)
    df.sub <- rbind(df.sub,
                    data.frame(term = 'Macro F1', method = method,
                               value = res.out$macro_f1, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'Accuracy', method = method,
                               value = res.out$accuracy, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'Balanced accuracy', method = method,
                               value = res.out$balanced.accuracy, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'labeled-Accuracy', method = method,
                               value = res.out$accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'labeled-Balanced accuracy', method = method,
                               value = res.out$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Median F1', method = method,
                               value = res.out$med_f1, stringsAsFactors = F))
    # cluster
    df_cluster <- readRDS(paste0(path_cluster, 'Campbell_Tasic.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Campbell_Tasic.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Drop-seq -> SMARTer, Mouse brain', nrow(df.sub))
    df_CP <- rbind(df_CP, df.sub)
    
    # Campbell -> Zeisel
    file.pred <- paste0(path.CP, 'pred_Campbell_Zeisel.Rdata')
    file.true <- paste0(path.CP, 'true_Campbell_Zeisel.Rdata')
    file.ref <- paste0(path.CP, 'ref_Campbell_Zeisel.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    ref.names <- names(table(ref.labels))
    uniform.names <- c("astrocytes_ependymal", "endothelial-mural", "astrocytes_ependymal", "endothelial-mural", 
                       "neurons", "oligodendrocytes", "OPC", "Pars tuberalis", "microglia", "Tanycytes", "VLMCs")
    df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
    all.cell <- names(table(true.tags))
    df.sc.names <- data.frame(sc.name = all.cell, name = all.cell)
    res.out <- simple.evaluation(true.tags, pred.tags, df.ref.names, df.sc.names)
    df.sub <- data.frame(term = 'Weighted macro F1', method = method,
                         value = res.out$weighted_macro_f1, stringsAsFactors = F)
    df.sub <- rbind(df.sub,
                    data.frame(term = 'Macro F1', method = method,
                               value = res.out$macro_f1, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'Accuracy', method = method,
                               value = res.out$accuracy, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'Balanced accuracy', method = method,
                               value = res.out$balanced.accuracy, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'labeled-Accuracy', method = method,
                               value = res.out$accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub <- rbind(df.sub,
                    data.frame(term = 'labeled-Balanced accuracy', method = method,
                               value = res.out$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Median F1', method = method,
                               value = res.out$med_f1, stringsAsFactors = F))
    # cluster
    df_cluster <- readRDS(paste0(path_cluster, 'Campbell_Zeisel.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Campbell_Zeisel.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Drop-seq -> Fluidigm C1, Mouse brain', nrow(df.sub))
    df_CP <- rbind(df_CP, df.sub)
    
}
df_CP$type <- rep('Cross-platform', nrow(df_CP))

# PBMC
vec.ref <- c('10Xv2', '10Xv2', '10Xv2', '10Xv2', '10Xv2', '10Xv2', 
             '10Xv3', '10Xv3', '10Xv3', '10Xv3', '10Xv3', '10Xv3', 
             'Drop-seq', 'Drop-seq', 'Drop-seq', 'Drop-seq', 'Drop-seq', 'Drop-seq', 
             'inDrops', 'inDrops', 'inDrops', 'inDrops', 'inDrops', 'inDrops', 
             'CEL-Seq2', 'Smart-seq2')
vec.query <- c('10Xv3', 'Drop-seq', 'inDrops', 'Seq-Well', 'CEL-Seq2', 'Smart-seq2',
               '10Xv2', 'Drop-seq', 'inDrops', 'Seq-Well', 'CEL-Seq2', 'Smart-seq2',
               '10Xv2', '10Xv3', 'inDrops', 'Seq-Well', 'CEL-Seq2', 'Smart-seq2',
               '10Xv2', '10Xv3', 'Drop-seq', 'Seq-Well', 'CEL-Seq2', 'Smart-seq2',
               'Smart-seq2', 'CEL-Seq2')
vec.dataset <- c('10x Chromium (v2) -> 10x Chromium (v3), Human PBMC', 
                 '10x Chromium (v2) -> Drop-seq, Human PBMC', '10x Chromium (v2) -> inDrops, Human PBMC',
                 '10x Chromium (v2) -> Seq-Well, Human PBMC', '10x Chromium (v2) -> CEL-Seq2, Human PBMC',
                 '10x Chromium (v2) -> Smart-seq2, Human PBMC',
                 '10x Chromium (v3) -> 10x Chromium (v2), Human PBMC', 
                 '10x Chromium (v3) -> Drop-seq, Human PBMC', '10x Chromium (v3) -> inDrops, Human PBMC',
                 '10x Chromium (v3) -> Seq-Well, Human PBMC', '10x Chromium (v3) -> CEL-Seq2, Human PBMC',
                 '10x Chromium (v3) -> Smart-seq2, Human PBMC',
                 'Drop-seq -> 10x Chromium (v2), Human PBMC', 'Drop-seq -> 10x Chromium (v3), Human PBMC',
                 'Drop-seq -> inDrops, Human PBMC', 'Drop-seq -> Seq-Well, Human PBMC', 
                 'Drop-seq -> CEL-Seq2, Human PBMC', 'Drop-seq -> Smart-seq2, Human PBMC',
                 'inDrops -> 10x Chromium (v2), Human PBMC', 'inDrops -> 10x Chromium (v3), Human PBMC',
                 'inDrops -> Drop-seq, Human PBMC', 'inDrops -> Seq-Well, Human PBMC', 
                 'inDrops -> CEL-Seq2, Human PBMC', 'inDrops -> Smart-seq2, Human PBMC',
                 'CEL-Seq2 -> Smart-seq2, Human PBMC', 'Smart-seq2 -> CEL-Seq2, Human PBMC')

df.pbmc <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path.pbmc <- paste0(path_CP, method, '/PBMC/')
    for (i in 1:length(vec.ref)) {
        ref <- vec.ref[i]
        query <- vec.query[i]
        if (ref == '10Xv2') {
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (ref == '10Xv3') {
                ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
            }
        }
        ref.names <- names(table(ref.labels))
        df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
        if (query == '10Xv2') {
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        all.cell <- names(table(label_sc))
        uniform.names <- all.cell
        uniform.names[!(uniform.names %in% setdiff(ref.names, 'Unassigned'))] <- 'Unassigned'
        df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)
        pred.tags <- readRDS(paste0(path.pbmc, ref, '_', query, '.Rdata'))
        res.out <- simple.evaluation(label_sc, pred.tags, df.ref.names, df.sc.names)
        df.sub <- data.frame(term = 'Weighted macro F1', method = method,
                             value = res.out$weighted_macro_f1, stringsAsFactors = F)
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Macro F1', method = method,
                                   value = res.out$macro_f1, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Accuracy', method = method,
                                   value = res.out$accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Accuracy', method = method,
                                   value = res.out$accuracy.rm.unassigned, stringsAsFactors = F))
        df.sub <- rbind(df.sub,
                        data.frame(term = 'labeled-Balanced accuracy', method = method,
                                   value = res.out$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Median F1', method = method,
                                   value = res.out$med_f1, stringsAsFactors = F))
        # cluster
        df_cluster <- readRDS(paste0(path_cluster, ref, '_', query, '.Rdata'))
        all_cluster <- length(table(df_cluster$cluster))
        correct_cluster <- func_eval_cluster(label_sc, pred.tags, df.ref.names, df.sc.names, df_cluster)
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Number of correctly labeled clusters', method = method,
                                   value = correct_cluster, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Total number of clusters', method = method,
                                   value = all_cluster, stringsAsFactors = F))
        # dataset
        list_dataset <- readRDS(paste0(path_eval_dataset, ref, '_', query, '.Rdata'))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'F1_asw', method = method,
                                   value = list_dataset$F1_asw, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Depth of reference', method = method,
                                   value = list_dataset$depth_ref, stringsAsFactors = F))
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Cell proportion covered by reference', method = method,
                                   value = list_dataset$prop_cover, stringsAsFactors = F))
        dataset <- vec.dataset[i]
        df.sub$dataset <- rep(dataset, nrow(df.sub))
        df.pbmc <- rbind(df.pbmc, df.sub)
        print(c(method, ref, query))
    }
}

df.pbmc$type <- rep('Cross-platform', nrow(df.pbmc))


df_1situ <- rbind(df_CV, df_crosssample, df_CP, df.pbmc)

file_1situ <- '/mdshare/node9/zy/MAGIC/evaluation/first_situation.txt'
write.table(df_1situ, file = file_1situ, row.names = F, sep = '\t')


df.acc <- df_1situ[df_1situ$term=='Accuracy',]


