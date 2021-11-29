library(reticulate)
use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')
library(Seurat)
library(SeuratData)
data("pbmcsca")

# evaluation
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
    accuracy.rm.unassigned <- metrics$accuracy_score(our.tag.rm, true.tag.rm)
    balanced.accuracy.rm.unassigned <- 
        metrics$balanced_accuracy_score(our.tag.rm, true.tag.rm)
    macro_f1.rm.unassigned <- metrics$f1_score(true.tag.rm, our.tag.rm, average = 'macro', labels = unique(true.tag.rm))
    
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
path_CS <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_species/'
path_eval_dataset <- '/mdshare/node9/zy/MAGIC/evaluation/datasets/'
path_cluster <- '/mdshare/node9/zy/MAGIC//MAGIC_cluster/'

vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')
vec.ref <- c('Baron', 'Baron', 'MCA', 'MCA')
vec.query <- c('Segerstolpe', 'Muraro', '10Xv2', 'pbmc3k')
vec.dataset <- c('Baron -> Muraro, Pancreas', 'Baron -> Segerstolpe, Pancreas', 
                 'Han -> Butler, PBMC', 'Han -> Ding, PBMC', 'Mean')

df_CS <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path_data <- paste0(path_CS, method, '/')
    
    ######### Baron -> Segerstolpe #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Baron_Segerstolpe.Rdata')
    file.true <- paste0(path_data, 'true_Baron_Segerstolpe.Rdata')
    file.ref <- paste0(path_data, 'ref_Baron_Segerstolpe.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    all.cell <- names(table(true.tags))
    df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
    uniform.names <- c("Unassigned", "activated_stellate", "alpha", "beta",
                       "delta", "ductal", "endothelial", "Unassigned",
                       "gamma", "macrophage", "Unassigned", "quiescent_stellate",
                       "schwann")
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
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Percent of unassigned', method = method,
                               value = res.out$percent.unassigned, stringsAsFactors = F))
    # cluster
    df_cluster <- readRDS(paste0(path_cluster, 'BaronM_Segerstolpe.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'BaronM_Segerstolpe.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Baron -> Segerstolpe, Pancreas', nrow(df.sub))
    df_CS <- rbind(df_CS, df.sub)
    ###############################################
    
    ######### Baron -> Muraro #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Baron_Muraro.Rdata')
    file.true <- paste0(path_data, 'true_Baron_Muraro.Rdata')
    file.ref <- paste0(path_data, 'ref_Baron_Muraro.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    all.cell <- names(table(true.tags))
    df.ref.names <- data.frame(ref.name = ref.names, name = ref.names)
    uniform.names <- c("Unassigned", "activated_stellate", "alpha", "beta",
                       "delta", "ductal", "endothelial", "Unassigned",
                       "gamma", "macrophage", "Unassigned", "quiescent_stellate",
                       "schwann")
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
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Percent of unassigned', method = method,
                               value = res.out$percent.unassigned, stringsAsFactors = F))
    # cluster
    df_cluster <- readRDS(paste0(path_cluster, 'BaronM_Muraro.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'BaronM_Muraro.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Baron -> Muraro, Pancreas', nrow(df.sub))
    df_CS <- rbind(df_CS, df.sub)
    ###############################################
    
    ################ MCA -> 10Xv2 #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_10Xv2.Rdata')
    file.true <- paste0(path_data, 'true_MCA_10Xv2.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_10Xv2.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    all.cell <- names(table(true.tags))
    uniform.names <- c("B cell", "Basophil", "Dendritic cell", "Erythroid cell",
                       "Monocyte", "Neutrophil", "NK cell", "T cell")
    df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
    uniform.names <- c("B cell", "Monocyte", "Monocyte", 
                       'T cell', 'T cell', "Dendritic cell", 
                       "Unassigned", "NK cell", "Dendritic cell")
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
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Percent of unassigned', method = method,
                               value = res.out$percent.unassigned, stringsAsFactors = F))
    # cluster
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_10Xv2.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_10Xv2.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Han -> Ding, PBMC', nrow(df.sub))
    df_CS <- rbind(df_CS, df.sub)
    ###############################################
    
    ################ MCA -> pbmc3k #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_pbmc3k.Rdata')
    file.true <- paste0(path_data, 'true_MCA_pbmc3k.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_pbmc3k.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    all.cell <- names(table(true.tags))
    uniform.names <- c("B cell", "Basophil", "Dendritic cell", "Erythroblast",
                       "Monocyte", "Neutrophil", "NK cell", "T cell")
    df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
    uniform.names <- c("B cell", "Monocyte", "T cell", "Dendritic cell",
                       "Monocyte", "T cell", "T cell", "NK cell",
                       "Unassigned")
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
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Percent of unassigned', method = method,
                               value = res.out$percent.unassigned, stringsAsFactors = F))
    # cluster
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_pbmc3k.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_pbmc3k.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Han -> Butler, PBMC', nrow(df.sub))
    df_CS <- rbind(df_CS, df.sub)
    ###############################################
    print(method)
}

df_CS$type <- rep('Cross-species', nrow(df_CS))

file_3situ <- '/mdshare/node9/zy/MAGIC/evaluation/third_situation.txt'
write.table(df_CS, file = file_3situ, row.names = F, sep = '\t')


