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
path_CP <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/'
path_MB <- '/mdshare/node9/zy/MAGIC/Benchmark/mouse_brain/'
path_eval_dataset <- '/mdshare/node9/zy/MAGIC/evaluation/datasets/'
path_cluster <- '/mdshare/node9/zy/MAGIC/MAGIC_cluster/'

# methods
vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')

# high-depth
# PBMC
vec.ref <- c('CEL-Seq2', 'CEL-Seq2', 'CEL-Seq2', 'CEL-Seq2', 'CEL-Seq2', 
             'Smart-seq2', 'Smart-seq2', 'Smart-seq2', 'Smart-seq2', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'Drop-seq', 'inDrops', 'Seq-Well',
               '10Xv2', '10Xv3', 'Drop-seq', 'inDrops', 'Seq-Well')
vec.dataset <- c('CEL-Seq2 -> 10x Chromium (v2), Human PBMC', 'CEL-Seq2 -> 10x Chromium (v3), Human PBMC',
                 'CEL-Seq2 -> Drop-seq, Human PBMC', 'CEL-Seq2 -> inDrops, Human PBMC',
                 'CEL-Seq2 -> Seq-Well, Human PBMC', 
                 'Smart-seq2 -> 10x Chromium (v2), Human PBMC', 'Smart-seq2 -> 10x Chromium (v3), Human PBMC',
                 'Smart-seq2 -> Drop-seq, Human PBMC', 'Smart-seq2 -> inDrops, Human PBMC',
                 'Smart-seq2 -> Seq-Well, Human PBMC')

df.plot.pbmc.HD <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path.pbmc <- paste0(path_CP, method, '/PBMC/')
    for (i in 1:length(vec.ref)) {
        ref <- vec.ref[i]
        query <- vec.query[i]
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
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
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Percent of unassigned', method = method,
                                   value = res.out$percent.unassigned, stringsAsFactors = F))
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
        df.plot.pbmc.HD <- rbind(df.plot.pbmc.HD, df.sub)
        print(c(method, ref, query))
    }
}

# mouse brain
df.plot.mouse.HD <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path_data <- paste0(path_MB, method, '/')
    
    ######### Tasic -> Tasic2018VISp #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Tasic_Tasic2018VISp.Rdata')
    file.true <- paste0(path_data, 'true_Tasic_Tasic2018VISp.Rdata')
    file.ref <- paste0(path_data, 'ref_Tasic_Tasic2018VISp.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Endothelial Cell", "GABA-ergic Neuron", "Glutamatergic Neuron",
                 "PVM & Microglia", "Oligodendrocyte & OPC", "Oligodendrocyte & OPC")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocyte", "Endothelial Cell", "GABA-ergic Neuron", "Glutamatergic Neuron",
                       "Oligodendrocyte & OPC", "Unassigned", "PVMs & Microglia", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'Tasic_Tasic2018VISp.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Tasic_Tasic2018VISp.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('SMARTer -> SMART-Seqv4(VISp), Mouse brain', nrow(df.sub))
    df.plot.mouse.HD <- rbind(df.plot.mouse.HD, df.sub)
    ###############################################
    
    ######### Tasic -> Tasic2018ALM #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Tasic_Tasic2018ALM.Rdata')
    file.true <- paste0(path_data, 'true_Tasic_Tasic2018ALM.Rdata')
    file.ref <- paste0(path_data, 'ref_Tasic_Tasic2018ALM.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Endothelial Cell", "GABA-ergic Neuron", "Glutamatergic Neuron",
                 "PVM & Microglia", "Oligodendrocyte & OPC", "Oligodendrocyte & OPC")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocyte", "Endothelial Cell", "GABA-ergic Neuron", "Glutamatergic Neuron",
                       "Oligodendrocyte & OPC", "Unassigned", "PVMs & Microglia", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'Tasic_Tasic2018ALM.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Tasic_Tasic2018ALM.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('SMARTer -> SMART-Seqv4(ALM), Mouse brain', nrow(df.sub))
    df.plot.mouse.HD <- rbind(df.plot.mouse.HD, df.sub)
    ###############################################
    
    ######### Tasic -> Campbell #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Tasic_Campbell.Rdata')
    file.true <- paste0(path_data, 'true_Tasic_Campbell.Rdata')
    file.ref <- paste0(path_data, 'ref_Tasic_Campbell.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocytes", "Endothelial cells",
                 "PVMs & Microglia", "Neurons", "Oligodendrocytes", "OPC")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocytes", "Endothelial cells", "Unassigned", "Unassigned", "Neurons",
                       "Oligodendrocytes", "OPC", "Unassigned", "PVMs & Microglia", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'Tasic_Campbell.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Tasic_Campbell.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('SMARTer -> Drop-seq(Campbell), Mouse brain', nrow(df.sub))
    df.plot.mouse.HD <- rbind(df.plot.mouse.HD, df.sub)
    ###############################################
    
    ######### Tasic -> Hochgerner #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Tasic_Hochgerner.Rdata')
    file.true <- paste0(path_data, 'true_Tasic_Hochgerner.Rdata')
    file.ref <- paste0(path_data, 'ref_Tasic_Hochgerner.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Endothelial Cell",
                 "Microglia", "Neuron", "Oligodendrocyte", "OPC")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocyte", "Unassigned", "Endothelial Cell", "Neuron", "Neuron",
                       "Neuron", "Neuron", "Microglia", "Neuron", "Neuron", "Neuron",
                       "Unassigned", "Unassigned", "Unassigned", "Unassigned", "Oligodendrocyte", 
                       "OPC", "Unassigned", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'Tasic_Hochgerner.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Tasic_Hochgerner.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('SMARTer -> 10x Chromium (v2), Mouse brain', nrow(df.sub))
    df.plot.mouse.HD <- rbind(df.plot.mouse.HD, df.sub)
    ###############################################
    
    
    ######### Tasic -> Mizrak #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_Tasic_Mizrak.Rdata')
    file.true <- paste0(path_data, 'true_Tasic_Mizrak.Rdata')
    file.ref <- paste0(path_data, 'ref_Tasic_Mizrak.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Endothelial Cell",
                 "Microglia", "Neuron", "Oligodendrocyte", "OPC")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Unassigned", "Astrocyte", "Unassigned", "Endothelial Cell", "Unassigned",
                       "Microglia", "Unassigned", "Neuron", "Oligodendrocyte", "OPC", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'Tasic_Mizrak.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'Tasic_Mizrak.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('SMARTer -> Drop-seq(Mizrak), Mouse brain', nrow(df.sub))
    df.plot.mouse.HD <- rbind(df.plot.mouse.HD, df.sub)
    ###############################################
    print(c(method))
}

df_HD <- rbind(df.plot.pbmc.HD, df.plot.mouse.HD)
df_HD$type <- rep('High-depth reference', nrow(df_HD))

# low-depth
# PBMC
vec.ref <- c('Seq-Well', 'Seq-Well', 'Seq-Well', 'Seq-Well', 'Seq-Well', 'Seq-Well')
vec.query <- c('10Xv2', '10Xv3', 'Drop-seq', 'inDrops', 'CEL-Seq2', 'Smart-seq2')
vec.dataset <- c('Seq-Well -> 10x Chromium (v2), Human PBMC', 'Seq-Well -> 10x Chromium (v3), Human PBMC',
                 'Seq-Well -> Drop-seq, Human PBMC', 'Seq-Well -> inDrops, Human PBMC',
                 'Seq-Well -> CEL-Seq2, Human PBMC', 'Seq-Well -> Smart-seq2, Human PBMC')

df.plot.pbmc.LD <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path.pbmc <- paste0(path_CP, method, '/PBMC/')
    for (i in 1:length(vec.ref)) {
        ref <- vec.ref[i]
        query <- vec.query[i]
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
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
        df.sub <- rbind(df.sub, 
                        data.frame(term = 'Percent of unassigned', method = method,
                                   value = res.out$percent.unassigned, stringsAsFactors = F))
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
        df.plot.pbmc.LD <- rbind(df.plot.pbmc.LD, df.sub)
        print(c(method, ref, query))
    }
}

# mouse brain
df.plot.mouse.LD <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    path_data <- paste0(path_MB, method, '/')
    
    ######### MCA -> Tasic2018VISp #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_Tasic2018VISp.Rdata')
    file.true <- paste0(path_data, 'true_MCA_Tasic2018VISp.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_Tasic2018VISp.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Astroglial cell", "Granulocyte", "Ependymocytes",
                 "PVMs & Microglia", "PVMs & Microglia", "Oligodendrocyte & OPC", "Neuron", 
                 "Oligodendrocyte & OPC", "Schwann cell")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocyte", "Unassigned", "Neuron", "Neuron", "Oligodendrocyte & OPC",
                       "Unassigned", "PVMs & Microglia", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_Tasic2018VISp.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_Tasic2018VISp.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Microwell-Seq -> SMART-Seqv4(VISp), Mouse brain', nrow(df.sub))
    df.plot.mouse.LD <- rbind(df.plot.mouse.LD, df.sub)
    ###############################################
    
    ######### MCA -> Tasic2018ALM #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_Tasic2018ALM.Rdata')
    file.true <- paste0(path_data, 'true_MCA_Tasic2018ALM.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_Tasic2018ALM.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Astroglial cell", "Granulocyte", "Ependymocytes",
                 "PVMs & Microglia", "PVMs & Microglia", "Oligodendrocyte & OPC", "Neuron", 
                 "Oligodendrocyte & OPC", "Schwann cell")
    if (method == "scSemiCluster") {
        sc.name <- c("Astrocyte", "Astroglial cell", "Granulocyte", "Ependymocytes",
                     "PVMs & Microglia", "PVMs & Microglia", "Oligodendrocyte & OPC", "Neuron", 
                     "Oligodendrocyte & OPC")
    }
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocyte", "Unassigned", "Neuron", "Neuron", "Oligodendrocyte & OPC",
                       "Unassigned", "PVMs & Microglia", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_Tasic2018ALM.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_Tasic2018ALM.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Microwell-Seq -> SMART-Seqv4(ALM), Mouse brain', nrow(df.sub))
    df.plot.mouse.LD <- rbind(df.plot.mouse.LD, df.sub)
    ###############################################
    
    ######### MCA -> Campbell #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_Campbell.Rdata')
    file.true <- paste0(path_data, 'true_MCA_Campbell.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_Campbell.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Astroglial cell", "Granulocyte", "Ependymocytes",
                 "PVMs & Microglia", "PVMs & Microglia", "Oligodendrocyte", "Neurons", 
                 "OPC", "Schwann cell")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocytes", "Unassigned", "Ependymocytes", "Unassigned", "Neurons",
                       "Oligodendrocyte", "OPC", "Unassigned", "PVMs & Microglia", "Unassigned", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_Campbell.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_Campbell.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Microwell-Seq -> Drop-seq(Campbell), Mouse brain', nrow(df.sub))
    df.plot.mouse.LD <- rbind(df.plot.mouse.LD, df.sub)
    ###############################################
    
    ######### MCA -> Hochgerner #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_Hochgerner.Rdata')
    file.true <- paste0(path_data, 'true_MCA_Hochgerner.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_Hochgerner.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- names(table(ref.labels))
    sc.name <- c("Astrocyte", "Astroglial cell", "Granulocyte", "Ependymocytes",
                 "PVM", "Microglia", "Oligodendrocyte", "Neuron", 
                 "OPC", "Schwann cell")
    df.ref.names <- data.frame(ref.name = ref.names, name = sc.name)
    all.cell <- names(table(true.tags))
    uniform.names <- c("Astrocyte", "Unassigned", "Unassigned", "Neuron", "Neuron",
                       "Neuron", "Neuron", "Microglia", "Neuron", "Neuron", 
                       "Neuron", "Unassigned", "Unassigned", "Unassigned", "Unassigned", 
                       "Oligodendrocyte", "OPC", "Unassigned", "PVM", "Unassigned")
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
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_Hochgerner.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_Hochgerner.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Microwell-Seq -> 10x Chromium (v2), Mouse brain', nrow(df.sub))
    df.plot.mouse.LD <- rbind(df.plot.mouse.LD, df.sub)
    ###############################################
    
    
    ######### MCA -> Mizrak #############
    ###############################################
    file.pred <- paste0(path_data, 'pred_MCA_Mizrak.Rdata')
    file.true <- paste0(path_data, 'true_MCA_Mizrak.Rdata')
    file.ref <- paste0(path_data, 'ref_MCA_Mizrak.Rdata')
    true.tags <- readRDS(file.true)
    ref.labels <- readRDS(file.ref)
    pred.tags <- readRDS(file.pred)
    
    ref.names <- unique(ref.labels)
    all.cell <- unique(true.tags)
    uniform.names <- c("Oligodendrocyte", "Microglia", "Astrocyte", "Neuron",
                       "Macrophage", "Granulocyte", "OPC",
                       "Schwann cell", "Astrocyte", "Ependymal")
    df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
    uniform.names <- c("Neuron", "Oligodendrocyte", "Unassigned", "Unassigned", "Ependymal", 
                       "Microglia", "Unassigned", "Astrocyte", "OPC", "Unassigned", 
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
    df_cluster <- readRDS(paste0(path_cluster, 'MCA_Mizrak.Rdata'))
    all_cluster <- length(table(df_cluster$cluster))
    correct_cluster <- func_eval_cluster(true.tags, pred.tags, df.ref.names, df.sc.names, df_cluster)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Number of correctly labeled clusters', method = method,
                               value = correct_cluster, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Total number of clusters', method = method,
                               value = all_cluster, stringsAsFactors = F))
    # dataset
    list_dataset <- readRDS(paste0(path_eval_dataset, 'MCA_Mizrak.Rdata'))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'F1_asw', method = method,
                               value = list_dataset$F1_asw, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Depth of reference', method = method,
                               value = list_dataset$depth_ref, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Cell proportion covered by reference', method = method,
                               value = list_dataset$prop_cover, stringsAsFactors = F))
    df.sub$dataset <- rep('Microwell-Seq -> Drop-seq(Mizrak), Mouse brain', nrow(df.sub))
    df.plot.mouse.LD <- rbind(df.plot.mouse.LD, df.sub)
    ###############################################
    print(c(method))
}

df_LD <- rbind(df.plot.pbmc.LD, df.plot.mouse.LD)
df_LD$type <- rep('Low-depth reference', nrow(df_LD))

df_2situ <- rbind(df_HD, df_LD)

file_2situ <- '/mdshare/node9/zy/MAGIC/evaluation/second_situation.txt'
write.table(df_2situ, file = file_2situ, row.names = F, sep = '\t')


df.acc <- df_2situ[df_2situ$term=='Accuracy',]

