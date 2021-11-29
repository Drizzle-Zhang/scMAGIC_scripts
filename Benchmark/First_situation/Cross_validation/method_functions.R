run_scMAGIC<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, atlas,NumGenes=2000,conf_cutoff=5,
                    GeneOrderPath = NULL){
  "
  run scMAGIC
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  library(stringr)
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[, Cells_to_Keep]
  colnames(Data) <- str_replace_all(colnames(Data), '_', '.')
  colnames(Data) <- str_replace_all(colnames(Data), '-', '.')
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               scMAGIC                                     #
  #############################################################################
  library(scMAGIC)
  
  True_Labels_scMAGIC <- list()
  Pred_Labels_scMAGIC <- list()
  Total_Time_scMAGIC <- list()
  
  for (i in c(1:n_folds)) {
    train_data <- Data[, Train_Idx[[i]]]
    train_label <- Labels[Train_Idx[[i]]]
    # train_set <- .generate_ref(train_data, train_label)
    test_set <- Data[, Test_Idx[[i]]]
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
      start_time <- Sys.time()
      # scMAGIC = scMAGIC(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
      #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
      #                   Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    } else {
      start_time <- Sys.time()
      output.scMAGIC <- scMAGIC(test_set, train_data, train_label,
                                atlas = atlas,
                                corr_use_HVGene1 = n_gene, corr_use_HVGene2 = n_gene,
                                threshold = conf_cutoff, num_threads = 4)
      pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
      label.scMAGIC <- as.character(pred.scMAGIC)
      # print(table(Labels[Test_Idx[[i]]], label.scMAGIC))
      end_time <- Sys.time()
    }
    Total_Time_scMAGIC[i] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_scMAGIC[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scMAGIC[i] <- list(label.scMAGIC)
  }
  True_Labels_scMAGIC <- as.vector(unlist(True_Labels_scMAGIC))
  Pred_Labels_scMAGIC <- as.vector(unlist(Pred_Labels_scMAGIC))
  Total_Time_scMAGIC <- as.vector(unlist(Total_Time_scMAGIC))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    write.csv(
      True_Labels_scMAGIC,
      paste('scMAGIC_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_scMAGIC,
      paste('scMAGIC_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Total_Time_scMAGIC,
      paste('scMAGIC_', NumGenes, '_Total_Time.csv', sep = ''),
      row.names = FALSE
    )
  } else{
    write.csv(True_Labels_scMAGIC, 'scMAGIC_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_scMAGIC, 'scMAGIC_Pred_Labels.csv', row.names = FALSE)
    write.csv(Total_Time_scMAGIC, 'scMAGIC_Total_Time.csv', row.names = FALSE)
  }
}


run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,
                      genes = "de", quantile.use = 0.8,
                      GeneOrderPath = NULL,NumGenes = NULL){
    "
  run SingleR
  Wrapper script to run SingleR on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
    
    # Data <- read.delim(DataPath,row.names = 1)
    # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    INPUT <- readRDS(DataPath)
    Data <- INPUT$data.filter
    Labels <- INPUT$label.filter
    Labels <- Labels[,1]
    load(CV_RDataPath)
    # Labels <- as.vector(Labels[,col_Index])
    Data <- Data[,Cells_to_Keep]
    Labels <- Labels[Cells_to_Keep]
    if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                               SingleR                                     #
    #############################################################################
    library(SingleR)
    library(Seurat)
    True_Labels_SingleR <- list()
    Pred_Labels_SingleR <- list()
    Total_Time_SingleR <- list()

    for (i in c(1:n_folds)){
      train_data <- as.matrix(Data[,Train_Idx[[i]]])
      train_label <- Labels[Train_Idx[[i]]]
      test_set <- as.matrix(Data[,Test_Idx[[i]]])
      if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
            start_time <- Sys.time()
            # singler = SingleR(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]], 
            #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]], 
            #                   Labels[Train_Idx[[i]]], numCores = 1)
            end_time <- Sys.time()
        }
        else{
            start_time <- Sys.time()
            ref.mtx <- train_data
            ref.labels <- train_label
            exp_sc_mat <- test_set
            label_sc <- Labels[Test_Idx[[i]]]
            set.seed(123)
            if (ncol(ref.mtx) > 6000) {
              sample.ref <- sample(colnames(ref.mtx), 6000)
              ref.labels <- ref.labels[colnames(ref.mtx) %in% sample.ref]
              train_set <- as.matrix(ref.mtx)[,sample.ref]
            } else {
              train_set <- as.matrix(ref.mtx)
            }
            if (ncol(exp_sc_mat) > 3000) {
              sample.cell <- sample(colnames(exp_sc_mat), 3000)
              label_sc <- label_sc[colnames(exp_sc_mat) %in% sample.cell]
              test_set <- as.matrix(exp_sc_mat)[,sample.cell]
            } else {
              test_set <- as.matrix(exp_sc_mat)
            }
            singler = SingleR(method = "single", sc_data = test_set, 
                              ref_data = train_data, 
                              types = train_label, 
                              genes = genes, quantile.use = quantile.use,
                              numCores = 2)
            end_time <- Sys.time()
        }
        Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
        
        True_Labels_SingleR[i] <- list(as.vector(label_sc))
        Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
    }
    True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
    Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
    Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))
    
    setwd(OutputDir)
    
    if(!is.null(GeneOrderPath) & !is.null(NumGenes)){
        write.csv(True_Labels_SingleR,paste('SingleR_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
        write.csv(Pred_Labels_SingleR,paste('SingleR_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
        write.csv(Total_Time_SingleR,paste('SingleR_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
    }
    else{
        write.csv(True_Labels_SingleR,'SingleR_True_Labels.csv',row.names = FALSE)
        write.csv(Pred_Labels_SingleR,'SingleR_Pred_Labels.csv',row.names = FALSE)
        write.csv(Total_Time_SingleR,'SingleR_Total_Time.csv',row.names = FALSE)
    }
}


run_scmapcluster <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                             threshold = 0.2,
                      GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scmap
  Wrapper script to run scmap on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  load(CV_RDataPath)
  # Labels <- as.vector(Labels[, col_Index])
  Data <- Data[,Cells_to_Keep]
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                 scmap                                     #
  #############################################################################
  library(scmap)
  library(SingleCellExperiment)
  True_Labels_scmapcluster <- list()
  Pred_Labels_scmapcluster <- list()
  True_Labels_scmapcell <- list()
  Pred_Labels_scmapcell <- list()
  Training_Time_scmapcluster <- list()
  Testing_Time_scmapcluster <- list()
  Training_Time_scmapcell <- list()
  Testing_Time_scmapcell <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)) {
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
      sce <-
        SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                      1, Train_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <-
        selectFeatures(sce, n_features = NumGenes, suppress_plot = TRUE)
      
      sce_test <-
        SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                      1, Test_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    }
    else{
      sce <-
        SingleCellExperiment(list(normcounts = Data[, Train_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <- selectFeatures(sce, suppress_plot = TRUE)
      
      sce_test <-
        SingleCellExperiment(list(normcounts = Data[, Test_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    }
    
    # scmap-cluster
    start_time <- Sys.time()
    sce <- indexCluster(sce)
    end_time <- Sys.time()
    Training_Time_scmapcluster[i] <-
      as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    start_time <- Sys.time()
    scmapCluster_results <-
      scmapCluster(projection = sce_test,
                   index_list = list(metadata(sce)$scmap_cluster_index), 
                   threshold = threshold)
    end_time <- Sys.time()
    Testing_Time_scmapcluster[i] <-
      as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_scmapcluster[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scmapcluster[i] <-
      list(scmapCluster_results$combined_labs)
    
  }
  
  True_Labels_scmapcluster <-
    as.vector(unlist(True_Labels_scmapcluster))
  Pred_Labels_scmapcluster <-
    as.vector(unlist(Pred_Labels_scmapcluster))
  Training_Time_scmapcluster <-
    as.vector(unlist(Training_Time_scmapcluster))
  Testing_Time_scmapcluster <-
    as.vector(unlist(Testing_Time_scmapcluster))

  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    write.csv(
      True_Labels_scmapcluster,
      paste('scmapcluster_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_scmapcluster,
      paste('scmapcluster_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Training_Time_scmapcluster,
      paste('scmapcluster_', NumGenes, '_Training_Time.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Testing_Time_scmapcluster,
      paste('scmapcluster_', NumGenes, '_Testing_Time.csv', sep = ''),
      row.names = FALSE
    )
  }
  else{
    write.csv(True_Labels_scmapcluster,
              'scmapcluster_True_Labels.csv',
              row.names = FALSE)
    write.csv(Pred_Labels_scmapcluster,
              'scmapcluster_Pred_Labels.csv',
              row.names = FALSE)
    write.csv(Training_Time_scmapcluster,
              'scmapcluster_Training_Time.csv',
              row.names = FALSE)
    write.csv(Testing_Time_scmapcluster,
              'scmapcluster_Testing_Time.csv',
              row.names = FALSE)
  }
}


run_scmapcell <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                          w = 1, threshold = 0.4,
                      GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scmap
  Wrapper script to run scmap on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  load(CV_RDataPath)
  # Labels <- as.vector(Labels[, col_Index])
  Data <- Data[,Cells_to_Keep]
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                 scmap                                     #
  #############################################################################
  library(scmap)
  library(SingleCellExperiment)
  True_Labels_scmapcluster <- list()
  Pred_Labels_scmapcluster <- list()
  True_Labels_scmapcell <- list()
  Pred_Labels_scmapcell <- list()
  Training_Time_scmapcluster <- list()
  Testing_Time_scmapcluster <- list()
  Training_Time_scmapcell <- list()
  Testing_Time_scmapcell <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)) {
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
      sce <-
        SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                      1, Train_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <-
        selectFeatures(sce, n_features = NumGenes, suppress_plot = TRUE)
      
      sce_test <-
        SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                      1, Test_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    }
    else{
      sce <-
        SingleCellExperiment(list(normcounts = Data[, Train_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <- selectFeatures(sce, suppress_plot = TRUE)
      
      sce_test <-
        SingleCellExperiment(list(normcounts = Data[, Test_Idx[[i]]]),
                             colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    }

    # scmap-cell
    start_time <- Sys.time()
    set.seed(1)
    sce <- indexCell(sce)
    end_time <- Sys.time()
    Training_Time_scmapcell[i] <-
      as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    start_time <- Sys.time()
    scmapCell_results <-
      scmapCell(sce_test, list(metadata(sce)$scmap_cell_index), w = w)
    w_2 <- ceiling(w * 3/10)
    scmapCell_clusters <-
      scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$cell_type1)),
                        w=w_2, threshold = threshold)
    end_time <- Sys.time()
    Testing_Time_scmapcell[i] <-
      as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_scmapcell[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scmapcell[i] <-
      list(scmapCell_clusters$combined_labs)
  }
  
  True_Labels_scmapcell <- as.vector(unlist(True_Labels_scmapcell))
  Pred_Labels_scmapcell <- as.vector(unlist(Pred_Labels_scmapcell))
  Training_Time_scmapcell <-
    as.vector(unlist(Training_Time_scmapcell))
  Testing_Time_scmapcell <-
    as.vector(unlist(Testing_Time_scmapcell))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    write.csv(
      True_Labels_scmapcell,
      paste('scmapcell_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_scmapcell,
      paste('scmapcell_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Training_Time_scmapcell,
      paste('scmapcell_', NumGenes, '_Training_Time.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Testing_Time_scmapcell,
      paste('scmapcell_', NumGenes, '_Testing_Time.csv', sep = ''),
      row.names = FALSE
    )
  }
  else{
    write.csv(True_Labels_scmapcell,
              'scmapcell_True_Labels.csv',
              row.names = FALSE)
    write.csv(Pred_Labels_scmapcell,
              'scmapcell_Pred_Labels.csv',
              row.names = FALSE)
    write.csv(Training_Time_scmapcell,
              'scmapcell_Training_Time.csv',
              row.names = FALSE)
    write.csv(Testing_Time_scmapcell,
              'scmapcell_Testing_Time.csv',
              row.names = FALSE)
  }
}


run_CHETAH<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                     clust_method = 'complete', n_genes = 0.8,
                     GeneOrderPath = NULL,NumGenes = NULL){
  "
  run CHETAH
  Wrapper script to run CHETAH on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    INPUT <- readRDS(DataPath)
    Data <- INPUT$data.filter
    Labels <- INPUT$label.filter
    Labels <- Labels[, 1]
    load(CV_RDataPath)
    # Labels <- as.vector(Labels[, col_Index])
    Data <- Data[, Cells_to_Keep]
    Labels <- Labels[Cells_to_Keep]
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                                CHETAH                                     #
    #############################################################################
    library(CHETAH)
    library(SingleCellExperiment)
    True_Labels_CHETAH <- list()
    Pred_Labels_CHETAH <- list()
    Total_Time_CHETAH <- list()
    Data = as.matrix(Data)
    
    for (i in c(1:n_folds)) {
        if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
            sce <-
                SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                                     1, Train_Idx[[i]]]),
                                     colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
            
            sce_test <-
                SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                                     1, Test_Idx[[i]]]),
                                     colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
            start_time <- Sys.time()
            depth_ref <- round(median(colSums(ref.mtx != 0)))
            n_genes <- n_genes * depth_ref
            sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce,
                                         clust_method = clust_method, n_genes = n_genes)
            end_time <- Sys.time()
        }
        else{
            sce <-
                SingleCellExperiment(assays = list(counts = Data[, Train_Idx[[i]]]),
                                     colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
            
            sce_test <-
                SingleCellExperiment(assays = list(counts = Data[, Test_Idx[[i]]]),
                                     colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
            start_time <- Sys.time()
            sce_test <-
              CHETAHclassifier(input = sce_test, ref_cells = sce, 
                               clust_method = clust_method, n_genes = n_genes)
            end_time <- Sys.time()
        }
        
        Total_Time_CHETAH[i] <-
            as.numeric(difftime(end_time, start_time, units = 'secs'))
        
        True_Labels_CHETAH[i] <- list(Labels[Test_Idx[[i]]])
        Pred_Labels_CHETAH[i] <- list(sce_test$celltype_CHETAH)
    }
    True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
    Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
    Total_Time_CHETAH <- as.vector(unlist(Total_Time_CHETAH))
    
    setwd(OutputDir)
    
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
        write.csv(
            True_Labels_CHETAH,
            paste('CHETAH_', NumGenes, '_True_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Pred_Labels_CHETAH,
            paste('CHETAH_', NumGenes, '_Pred_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Total_Time_CHETAH,
            paste('CHETAH_', NumGenes, '_Total_Time.csv', sep = ''),
            row.names = FALSE
        )
    }
    else{
        write.csv(True_Labels_CHETAH, 'CHETAH_True_Labels.csv', row.names = FALSE)
        write.csv(Pred_Labels_CHETAH, 'CHETAH_Pred_Labels.csv', row.names = FALSE)
        write.csv(Total_Time_CHETAH, 'CHETAH_Total_Time.csv', row.names = FALSE)
    }
}


run_scPred<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                     nfeatures = 1000, npcs = 50,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scPred
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    INPUT <- readRDS(DataPath)
    Data <- INPUT$data.filter
    Labels <- INPUT$label.filter
    Labels <- Labels[, 1]
    load(CV_RDataPath)
    # Labels <- as.vector(Labels[, col_Index])
    Data <- Data[, Cells_to_Keep]
    Labels <- Labels[Cells_to_Keep]
    # print(dim(Data))
    # print(length(Labels))
    
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                                scPred                                     #
    #############################################################################
    library("scPred")
    library("Seurat")
    library("magrittr")
    True_Labels_scPred <- list()
    Pred_Labels_scPred <- list()
    Training_Time_scPred <- list()
    Testing_Time_scPred <- list()
    Data = as.matrix(Data)
    
    for (i in c(1:n_folds)) {
        if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
            sce <-
                SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                                1, Train_Idx[[i]]]),
                                     colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
            sce_counts <- normcounts(sce)
            sce_cpm <-
                apply(sce_counts, 2, function(x)
                    (x / sum(x)) * 1000000)
            sce_metadata <- as.data.frame(colData(sce))
            
            sce_test <-
                SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
                                                                1, Test_Idx[[i]]]),
                                     colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
            sce_counts_test <- normcounts(sce_test)
            sce_cpm_test <-
                apply(sce_counts_test, 2, function(x)
                    (x / sum(x)) * 1000000)
            sce_metadata_test <- as.data.frame(colData(sce_test))
        } else{
            reference <- CreateSeuratObject(counts = Data[, Train_Idx[[i]]])
            reference@meta.data$cell_type <- Labels[Train_Idx[[i]]]
            query <- CreateSeuratObject(counts = Data[, Test_Idx[[i]]])
            # sce <- SingleCellExperiment(list(normcounts = Data[,Train_Idx[[i]]]),
            #                             colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
            # sce_counts <- normcounts(sce)
            # sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
            # sce_metadata <- as.data.frame(colData(sce))
            #
            # sce_test <- SingleCellExperiment(list(normcounts = Data[,Test_Idx[[i]]]),
            #                                  colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
            # sce_counts_test <- normcounts(sce_test)
            # sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
            # sce_metadata_test <- as.data.frame(colData(sce_test))
        }
        
        
        # scPred Training
        start_time <- Sys.time()
        ref.mtx <- Data[, Train_Idx[[i]]]
        ref.labels <- Labels[Train_Idx[[i]]]
        exp_sc_mat <- Data[, Test_Idx[[i]]]
        # scPred Training    
        reference <- CreateSeuratObject(counts = ref.mtx)
        reference@meta.data$cell_type <- ref.labels
        reference <- reference %>%
          NormalizeData() %>%
          FindVariableFeatures(nfeatures = nfeatures) %>%
          ScaleData() %>%
          RunPCA(npcs = npcs)
        reference <- getFeatureSpace(reference, "cell_type")
        reference <- trainModel(reference)
        # scPred Prediction
        query <- CreateSeuratObject(counts = exp_sc_mat)
        query <- NormalizeData(query)
        query <- scPredict(query, reference)
        pred.scPred <- query@meta.data$scpred_prediction
        end_time <- Sys.time()
        Testing_Time_scPred[i] <-
            as.numeric(difftime(end_time, start_time, units = 'secs'))
        
        True_Labels_scPred[i] <- list(Labels[Test_Idx[[i]]])
        Pred_Labels_scPred[i] <- list(pred.scPred)
        # Pred_Labels_scPred[i] <- list(getPredictions(scp)$predClass)
    }
    True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
    Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
    Training_Time_scPred <- as.vector(unlist(Training_Time_scPred))
    Testing_Time_scPred <- as.vector(unlist(Testing_Time_scPred))
    Pred_Labels_scPred <-
        gsub('pyramidal.SS', 'pyramidal SS', Pred_Labels_scPred)
    Pred_Labels_scPred <-
        gsub('pyramidal.CA1', 'pyramidal CA1', Pred_Labels_scPred)
    Pred_Labels_scPred <-
        gsub('endothelial.mural', 'endothelial-mural', Pred_Labels_scPred)
    
    setwd(OutputDir)
    
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
        write.csv(
            True_Labels_scPred,
            paste('scPred_', NumGenes, '_True_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Pred_Labels_scPred,
            paste('scPred_', NumGenes, '_Pred_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Training_Time_scPred,
            paste('scPred_', NumGenes, '_Training_Time.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Testing_Time_scPred,
            paste('scPred_', NumGenes, '_Testing_Time.csv', sep = ''),
            row.names = FALSE
        )
    }
    else{
        write.csv(True_Labels_scPred, 'scPred_True_Labels.csv', row.names = FALSE)
        write.csv(Pred_Labels_scPred, 'scPred_Pred_Labels.csv', row.names = FALSE)
        write.csv(Training_Time_scPred,
                  'scPred_Training_Time.csv',
                  row.names = FALSE)
        write.csv(Testing_Time_scPred,
                  'scPred_Testing_Time.csv',
                  row.names = FALSE)
    }
}


run_sciBet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,k=1000,
                     GeneOrderPath = NULL,NumGenes = NULL){
  "
  run sciBet
  Wrapper script to run sciBet on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    INPUT <- readRDS(DataPath)
    Data <- INPUT$data.filter
    Labels <- INPUT$label.filter
    Labels <- Labels[, 1]
    load(CV_RDataPath)
    Data <- Data[, Cells_to_Keep]
    Labels <- Labels[Cells_to_Keep]
    df_Labels <- data.frame(Labels = Labels, row.names = colnames(Data))
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                               sciBet                                     #
    #############################################################################
    suppressMessages(library(tidyverse))
    suppressMessages(library(scibet))
    suppressMessages(library(viridis))
    suppressMessages(library(ggsci))
    True_Labels_sciBet <- list()
    Pred_Labels_sciBet <- list()
    Total_Time_sciBet <- list()
    # Data = t(as.matrix(Data))
    # names(Data) <- 1:dim(Data)[2]
    Data <- Data / 1.0
    
    for (i in c(1:n_folds)) {
        train_cells <- colnames(Data[, Train_Idx[[i]]])
        test_cells <- colnames(Data[, Test_Idx[[i]]])
        train_labels <- df_Labels[train_cells, 'Labels']
        test_labels <- df_Labels[test_cells, 'Labels']
        train_set <- as.data.frame(t(Data[, Train_Idx[[i]]]))
        train_cells <- rownames(train_set)
        train_set$label <- train_labels
        test_set <- as.data.frame(t(Data[, Test_Idx[[i]]]))
        if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
            start_time <- Sys.time()
            # sciBet = sciBet(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
            #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
            #                   Labels[Train_Idx[[i]]], numCores = 1)
            end_time <- Sys.time()
        } else {
            start_time <- Sys.time()
            sciBet = SciBet(train_set, test_set, k=k)
            end_time <- Sys.time()
        }
        Total_Time_sciBet[i] <-
            as.numeric(difftime(end_time, start_time, units = 'secs'))
        
        True_Labels_sciBet[i] <- list(test_labels)
        Pred_Labels_sciBet[i] <- list(sciBet)
    }
    True_Labels_sciBet <- as.vector(unlist(True_Labels_sciBet))
    Pred_Labels_sciBet <- as.vector(unlist(Pred_Labels_sciBet))
    Total_Time_sciBet <- as.vector(unlist(Total_Time_sciBet))
    
    setwd(OutputDir)
    
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
        write.csv(
            True_Labels_sciBet,
            paste('sciBet_', NumGenes, '_True_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Pred_Labels_sciBet,
            paste('sciBet_', NumGenes, '_Pred_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Total_Time_sciBet,
            paste('sciBet_', NumGenes, '_Total_Time.csv', sep = ''),
            row.names = FALSE
        )
    }
    else{
        write.csv(True_Labels_sciBet, 'sciBet_True_Labels.csv', row.names = FALSE)
        write.csv(Pred_Labels_sciBet, 'sciBet_Pred_Labels.csv', row.names = FALSE)
        write.csv(Total_Time_sciBet, 'sciBet_Total_Time.csv', row.names = FALSE)
    }
}


run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                            nTopGenes = 10, nTopGenePairs = 25, nTrees = 1000,
                            GeneOrderPath = NULL,NumGenes = NULL){
  "
  run singleCellNet
  Wrapper script to run singleCellNet on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
    INPUT <- readRDS(DataPath)
    Data <- INPUT$data.filter
    Labels <- INPUT$label.filter
    Labels <- Labels[, 1]
    # Data <- read.delim(DataPath, row.names = 1)
    genes <- rownames(Data)
    genes <- gsub('_', '.', genes)
    rownames(Data) <- genes
    # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    load(CV_RDataPath)
    # Labels <- as.vector(Labels[,col_Index])
    Data <- Data[, Cells_to_Keep]
    Labels <- Labels[Cells_to_Keep]
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                              singleCellNet                                #
    #############################################################################
    library(singleCellNet)
    library(dplyr)
    True_Labels_singleCellNet <- list()
    Pred_Labels_singleCellNet <- list()
    Training_Time_singleCellNet <- list()
    Testing_Time_singleCellNet <- list()
    # Data = as.matrix(Data)            # deals also with sparse matrix
    
    for (i in c(1:n_folds)) {
        if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
            DataTrain <-
                Data[as.vector(GenesOrder[c(1:NumGenes), i]) + 1, Train_Idx[[i]]]
            DataTest <-
                Data[as.vector(GenesOrder[c(1:NumGenes), i]) + 1, Test_Idx[[i]]]
        } else{
            DataTrain <- as.matrix(Data[, Train_Idx[[i]]])
            LabelsTrain <-
                data.frame(Annotation = Labels[Train_Idx[[i]]],
                           row.names = colnames(DataTrain))
            DataTest <- as.matrix(Data[, Test_Idx[[i]]])
        }
        
        start_time <- Sys.time()
        class_info <-
            scn_train(stTrain = LabelsTrain,
                      expTrain = DataTrain,
                      dLevel = "Annotation", 
                      nTopGenes = nTopGenes, nTopGenePairs = nTopGenePairs, nTrees = nTrees)
        end_time <- Sys.time()
        Training_Time_singleCellNet[i] <-
            as.numeric(difftime(end_time, start_time, units = 'secs'))
        
        start_time <- Sys.time()
        classRes <-
            scn_predict(cnProc = class_info[['cnProc']],
                        expDat = DataTest,
                        nrand = 50)
        end_time <- Sys.time()
        Testing_Time_singleCellNet[i] <-
            as.numeric(difftime(end_time, start_time, units = 'secs'))
        
        True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
        Pred_Labels_singleCellNet[i] <-
            list((rownames(classRes)[apply(classRes, 2, which.max)])[1:length(Test_Idx[[i]])])
    }
    True_Labels_singleCellNet <-
        as.vector(unlist(True_Labels_singleCellNet))
    Pred_Labels_singleCellNet <-
        as.vector(unlist(Pred_Labels_singleCellNet))
    Training_Time_singleCellNet <-
        as.vector(unlist(Training_Time_singleCellNet))
    Testing_Time_singleCellNet <-
        as.vector(unlist(Testing_Time_singleCellNet))
    
    setwd(OutputDir)
    
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
        write.csv(
            True_Labels_singleCellNet,
            paste('singleCellNet_', NumGenes, '_True_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Pred_Labels_singleCellNet,
            paste('singleCellNet_', NumGenes, '_Pred_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Training_Time_singleCellNet,
            paste('singleCellNet_', NumGenes, '_Training_Time.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Testing_Time_singleCellNet,
            paste('singleCellNet_', NumGenes, '_Testing_Time.csv', sep = ''),
            row.names = FALSE
        )
    }
    else{
        write.csv(True_Labels_singleCellNet,
                  'singleCellNet_True_Labels.csv',
                  row.names = FALSE)
        write.csv(Pred_Labels_singleCellNet,
                  'singleCellNet_Pred_Labels.csv',
                  row.names = FALSE)
        write.csv(Training_Time_singleCellNet,
                  'singleCellNet_Training_Time.csv',
                  row.names = FALSE)
        write.csv(Testing_Time_singleCellNet,
                  'singleCellNet_Testing_Time.csv',
                  row.names = FALSE)
    }
}


run_CaSTLe<-function(DataPath,LabelsPath,CV_RDataPath, OutputDir, GeneOrderPath = NULL, NumGenes = NULL){
  "
  run CaSTLe
  Wrapper script to run CaSTLe on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.delim(DataPath,row.names = 1)
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                CaSTLe                                     #
  #############################################################################
  library(igraph)
  library(xgboost)
  True_Labels_Castle <- list()
  Pred_Labels_Castle <- list()
  Training_Time_Castle <- list()
  Testing_Time_Castle <- list()
  Data <- t(as.matrix(Data))
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  for(i in c(1:n_folds)){
    # 1. Load datasets
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      ds1 = Data[Train_Idx[[i]],as.vector(GenesOrder[c(1:NumGenes),i])+1]
      ds2 = Data[Test_Idx[[i]],as.vector(GenesOrder[c(1:NumGenes),i])+1]
    }
    else{
      ds1 = Data[Train_Idx[[i]],]
      ds2 = Data[Test_Idx[[i]],]
    }
    
    sourceCellTypes = as.factor(Labels[Train_Idx[[i]]])
    targetCellTypes = as.factor(Labels[Test_Idx[[i]]])
    
    start_time <- Sys.time()
    # 2. Unify sets, excluding low expressed genes
    source_n_cells_counts = apply(ds1, 2, function(x) { sum(x > 0) } )
    target_n_cells_counts = apply(ds2, 2, function(x) { sum(x > 0) } )
    common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                              colnames(ds2)[target_n_cells_counts>10])
    remove(source_n_cells_counts, target_n_cells_counts)
    ds1 = ds1[, colnames(ds1) %in% common_genes]
    ds2 = ds2[, colnames(ds2) %in% common_genes]
    ds = rbind(ds1[,common_genes], ds2[,common_genes])
    isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
    remove(ds1, ds2)
    
    # 3. Highest mean in both source and target
    topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
    end_time <- Sys.time()
    Training_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    # for each cell - what is the most probable classification?
    L = length(levels(sourceCellTypes))
    targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))
    
    for (cellType in levels(sourceCellTypes)) {
      
      inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
      
      # 4. Highest mutual information in source
      topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
      
      # 5. Top n genes that appear in both mi and avg
      selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
      
      # 6. remove correlated features
      tmp = cor(ds[,selectedFeatures], method = "pearson")
      tmp[!lower.tri(tmp)] = 0
      selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
      remove(tmp)
      
      # 7,8. Convert data from continous to binned dummy vars
      # break datasets to bins
      dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
      # use only bins with more than one value
      nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
      # convert to dummy vars
      ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
      remove(dsBins, nUniq)
      
      cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
      
      inTypeSource = sourceCellTypes == cellType
      # 9. Classify
      xg=xgboost(data=ds0[isSource,] , 
                 label=inTypeSource,
                 objective="binary:logistic", 
                 eta=0.7 , nthread=1, nround=20, verbose=0,
                 gamma=0.001, max_depth=5, min_child_weight=10)
      
      # 10. Predict
      inTypeProb = predict(xg, ds0[!isSource, ])
      
      targetClassification[cellType,] = inTypeProb
    }
    end_time <- Sys.time()
    Testing_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_Castle[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_Castle[i] <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  }
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  Training_Time_Castle <- as.vector(unlist(Training_Time_Castle))
  Testing_Time_Castle <- as.vector(unlist(Testing_Time_Castle))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_Castle,paste('True_Labels_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_Castle,paste('Pred_Labels_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_Castle,paste('Training_Time_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_Castle,paste('Testing_Time_Castle_',NumGenes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_Castle,'CaSTLe_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_Castle,'CaSTLe_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_Castle,'CaSTLe_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_Castle,'CaSTLe_Testing_Time.csv',row.names = FALSE)
  }
  
}


run_scID<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                   logFC = 0.5,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scID
  Wrapper script to run scID on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  # Labels <- as.vector(Labels[, col_Index])
  Data <- Data[, Cells_to_Keep]
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                 scID                                      #
  #############################################################################
  library(scID)
  library(Seurat)
  True_Labels_scID <- list()
  Pred_Labels_scID <- list()
  Total_Time_scID <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)) {
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
      # Train_Labels <- list(Labels[Train_Idx[[i]]])
      # names(Train_Labels[[1]]) <-
      #     colnames(Data[as.vector(GenesOrder[c(1:NumGenes), i]) + 1, Train_Idx[[i]]])
      # start_time <- Sys.time()
      # scID_output <-
      #     scid_multiclass(Data[as.vector(GenesOrder[c(1:NumGenes), i]) + 1, Test_Idx[[i]]],
      #                     Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
      #                              1, Train_Idx[[i]]],
      #                     Train_Labels[[1]])
      # end_time <- Sys.time()
    }
    else{
      Train_Labels <- list(Labels[Train_Idx[[i]]])
      names(Train_Labels[[1]]) <- colnames(Data[, Train_Idx[[i]]])
      start_time <- Sys.time()
      scID_output <-
        scid_multiclass(Data[, Test_Idx[[i]]], Data[, Train_Idx[[i]]], Train_Labels[[1]],
                        logFC = logFC)
      end_time <- Sys.time()
    }
    Total_Time_scID[i] <-
      as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_scID[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scID[i] <- list(as.vector(scID_output$labels))
  }
  True_Labels_scID <- as.vector(unlist(True_Labels_scID))
  Pred_Labels_scID <- as.vector(unlist(Pred_Labels_scID))
  Total_Time_scID <- as.vector(unlist(Total_Time_scID))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    write.csv(
      True_Labels_scID,
      paste('scID_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_scID,
      paste('scID_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Total_Time_scID,
      paste('scID_', NumGenes, '_Total_Time.csv', sep = ''),
      row.names = FALSE
    )
  }
  else{
    write.csv(True_Labels_scID, 'scID_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_scID, 'scID_Pred_Labels.csv', row.names = FALSE)
    write.csv(Total_Time_scID, 'scID_Total_Time.csv', row.names = FALSE)
  }
}


run_scClassify<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, tree = "HOPACH", 
                         algorithm = "KNN", k = 10, topN = 50,
                    GeneOrderPath = NULL,NumGenes = NULL){
    "
  run scClassify
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
    
    library(stringr)
    INPUT <- readRDS(DataPath)
    Data <- INPUT$data.filter
    Labels <- INPUT$label.filter
    Labels <- Labels[, 1]
    # Data <- read.delim(DataPath,row.names = 1)
    # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
    load(CV_RDataPath)
    Data <- Data[, Cells_to_Keep]
    colnames(Data) <- str_replace_all(colnames(Data), '_', '.')
    colnames(Data) <- str_replace_all(colnames(Data), '-', '.')
    Labels <- Labels[Cells_to_Keep]
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
        GenesOrder = read.csv(GeneOrderPath)
    }
    
    #############################################################################
    #                               scClassify                                     #
    #############################################################################
    library("scClassify")
    library(Matrix)
    True_Labels_scClassify <- list()
    Pred_Labels_scClassify <- list()
    Total_Time_scClassify <- list()
    
    for (i in c(1:n_folds)) {
        train_data <- Data[, Train_Idx[[i]]]
        train_label <- Labels[Train_Idx[[i]]]
        # train_set <- .generate_ref(train_data, train_label)
        test_set <- Data[, Test_Idx[[i]]]
        if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
            start_time <- Sys.time()
            # scMAGIC = scMAGIC(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
            #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
            #                   Labels[Train_Idx[[i]]], numCores = 1)
            end_time <- Sys.time()
        } else {
            start_time <- Sys.time()
            exprsMat_train <- as(as.matrix(log1p(train_data)), "dgCMatrix")
            exp_sc_mat <- as(as.matrix(log1p(test_set)), "dgCMatrix")
            scClassify_res <- scClassify(exprsMat_train = exprsMat_train,
                                         cellTypes_train = train_label,
                                         exprsMat_test = list(one = exp_sc_mat),
                                         tree = tree,
                                         algorithm = algorithm,
                                         k = k,
                                         topN = topN,
                                         returnList = FALSE,
                                         verbose = FALSE)
            pred.scClassify <- scClassify_res$testRes$one[[1]]$predRes
            end_time <- Sys.time()
        }
        Total_Time_scClassify[i] <-
            as.numeric(difftime(end_time, start_time, units = 'secs'))
        
        True_Labels_scClassify[i] <- list(Labels[Test_Idx[[i]]])
        Pred_Labels_scClassify[i] <- list(pred.scClassify)
    }
    True_Labels_scClassify <- as.vector(unlist(True_Labels_scClassify))
    Pred_Labels_scClassify <- as.vector(unlist(Pred_Labels_scClassify))
    Total_Time_scClassify <- as.vector(unlist(Total_Time_scClassify))
    
    setwd(OutputDir)
    
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
        write.csv(
            True_Labels_scClassify,
            paste('scClassify_', NumGenes, '_True_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Pred_Labels_scClassify,
            paste('scClassify_', NumGenes, '_Pred_Labels.csv', sep = ''),
            row.names = FALSE
        )
        write.csv(
            Total_Time_scClassify,
            paste('scClassify_', NumGenes, '_Total_Time.csv', sep = ''),
            row.names = FALSE
        )
    } else{
        write.csv(True_Labels_scClassify, 'scClassify_True_Labels.csv', row.names = FALSE)
        write.csv(Pred_Labels_scClassify, 'scClassify_Pred_Labels.csv', row.names = FALSE)
        write.csv(Total_Time_scClassify, 'scClassify_Total_Time.csv', row.names = FALSE)
    }
}


run_CALLR<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,k=0.9,
                    GeneOrderPath = NULL){
  "
  run CALLR
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  NumGenes = NULL
  library(stringr)
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[, Cells_to_Keep]
  colnames(Data) <- str_replace_all(colnames(Data), '_', '.')
  colnames(Data) <- str_replace_all(colnames(Data), '-', '.')
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               CALLR                                     #
  #############################################################################
  source('/local/wzk/data/scmagic/code/source/scMAGIC_script/methods_functions.R')
  
  True_Labels_CALLR <- list()
  Pred_Labels_CALLR <- list()
  Total_Time_CALLR <- list()
  
  for (i in c(1:n_folds)) {
    print(paste0('Running Fold ', as.character(i)))
    train_data <- Data[, Train_Idx[[i]]]
    train_label <- Labels[Train_Idx[[i]]]
    # train_set <- .generate_ref(train_data, train_label)
    test_set <- Data[, Test_Idx[[i]]]
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
      start_time <- Sys.time()
      # CALLR = CALLR(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
      #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
      #                   Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    } else {
      start_time <- Sys.time()
      pred.list <- func_CALLR(test_set, train_data, train_label, k=k)
      pred.CALLR = pred.list[[1]]
      pred.index = pred.list[[2]]
      label.CALLR <- as.character(pred.CALLR)
      end_time <- Sys.time()
    }
    Total_Time_CALLR[i] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_CALLR[i] <- list(Labels[Test_Idx[[i]]][pred.index])
    Pred_Labels_CALLR[i] <- list(label.CALLR)
  }
  True_Labels_CALLR <- as.vector(unlist(True_Labels_CALLR))
  Pred_Labels_CALLR <- as.vector(unlist(Pred_Labels_CALLR))
  Total_Time_CALLR <- as.vector(unlist(Total_Time_CALLR))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    write.csv(
      True_Labels_CALLR,
      paste('CALLR_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_CALLR,
      paste('CALLR_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Total_Time_CALLR,
      paste('CALLR_', NumGenes, '_Total_Time.csv', sep = ''),
      row.names = FALSE
    )
  } else{
    write.csv(True_Labels_CALLR, 'CALLR_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_CALLR, 'CALLR_Pred_Labels.csv', row.names = FALSE)
    write.csv(Total_Time_CALLR, 'CALLR_Total_Time.csv', row.names = FALSE)
  }
}



run_SVM<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir, 
                     nfeatures = 2000, threshold = 0.1,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run SVM
  Wrapper script to run SVM on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  load(CV_RDataPath)
  # Labels <- as.vector(Labels[, col_Index])
  Data <- Data[, Cells_to_Keep]
  Labels <- Labels[Cells_to_Keep]
  # print(dim(Data))
  # print(length(Labels))
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                SVM                                     #
  #############################################################################
  source('/local/zy/scMAGIC_script/methods_functions.R')
  True_Labels_SVM <- list()
  Pred_Labels_SVM <- list()
  Training_Time_SVM <- list()
  Testing_Time_SVM <- list()
  Data = as.matrix(Data)
  
  for (i in c(1:n_folds)) {
    if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
      # sce <-
      #   SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
      #                                                 1, Train_Idx[[i]]]),
      #                        colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      # sce_counts <- normcounts(sce)
      # sce_cpm <-
      #   apply(sce_counts, 2, function(x)
      #     (x / sum(x)) * 1000000)
      # sce_metadata <- as.data.frame(colData(sce))
      # 
      # sce_test <-
      #   SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes), i]) +
      #                                                 1, Test_Idx[[i]]]),
      #                        colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      # sce_counts_test <- normcounts(sce_test)
      # sce_cpm_test <-
      #   apply(sce_counts_test, 2, function(x)
      #     (x / sum(x)) * 1000000)
      # sce_metadata_test <- as.data.frame(colData(sce_test))
    } else{
      # reference <- CreateSeuratObject(counts = Data[, Train_Idx[[i]]])
      # reference@meta.data$cell_type <- Labels[Train_Idx[[i]]]
      # query <- CreateSeuratObject(counts = Data[, Test_Idx[[i]]])
      # sce <- SingleCellExperiment(list(normcounts = Data[,Train_Idx[[i]]]),
      #                             colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      # sce_counts <- normcounts(sce)
      # sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
      # sce_metadata <- as.data.frame(colData(sce))
      #
      # sce_test <- SingleCellExperiment(list(normcounts = Data[,Test_Idx[[i]]]),
      #                                  colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      # sce_counts_test <- normcounts(sce_test)
      # sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
      # sce_metadata_test <- as.data.frame(colData(sce_test))
    }
    
    
    #  Training
    start_time <- Sys.time()
    ref.mtx <- Data[, Train_Idx[[i]]]
    ref.labels <- Labels[Train_Idx[[i]]]
    exp_sc_mat <- Data[, Test_Idx[[i]]]
    pred.SVM <- func_SVM(exp_sc_mat, ref.mtx, ref.labels, nfeatures, threshold)
    end_time <- Sys.time()
    Testing_Time_SVM[i] <-
      as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_SVM[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SVM[i] <- list(pred.SVM)
    # Pred_Labels_SVM[i] <- list(getPredictions(scp)$predClass)
  }
  True_Labels_SVM <- as.vector(unlist(True_Labels_SVM))
  Pred_Labels_SVM <- as.vector(unlist(Pred_Labels_SVM))
  Training_Time_SVM <- as.vector(unlist(Training_Time_SVM))
  Testing_Time_SVM <- as.vector(unlist(Testing_Time_SVM))
  # Pred_Labels_SVM <-
  #   gsub('pyramidal.SS', 'pyramidal SS', Pred_Labels_SVM)
  # Pred_Labels_SVM <-
  #   gsub('pyramidal.CA1', 'pyramidal CA1', Pred_Labels_SVM)
  # Pred_Labels_SVM <-
  #   gsub('endothelial.mural', 'endothelial-mural', Pred_Labels_SVM)
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null (NumGenes)) {
    write.csv(
      True_Labels_SVM,
      paste('SVM_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_SVM,
      paste('SVM_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Training_Time_SVM,
      paste('SVM_', NumGenes, '_Training_Time.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Testing_Time_SVM,
      paste('SVM_', NumGenes, '_Testing_Time.csv', sep = ''),
      row.names = FALSE
    )
  }
  else{
    write.csv(True_Labels_SVM, 'SVM_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_SVM, 'SVM_Pred_Labels.csv', row.names = FALSE)
    write.csv(Training_Time_SVM,
              'SVM_Training_Time.csv',
              row.names = FALSE)
    write.csv(Testing_Time_SVM,
              'SVM_Testing_Time.csv',
              row.names = FALSE)
  }
}


run_scSemiCluster<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,k=c(1000,128,64,32),
                            GeneOrderPath = NULL){
  "
  run scSemiCluster
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  NumGenes = NULL
  library(stringr)
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[, Cells_to_Keep]
  colnames(Data) <- str_replace_all(colnames(Data), '_', '.')
  colnames(Data) <- str_replace_all(colnames(Data), '-', '.')
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               scSemiCluster                                     #
  #############################################################################
  source('/local/wzk/data/scmagic/code/source/scMAGIC_script/methods_functions.R')
  
  True_Labels_scSemiCluster <- list()
  Pred_Labels_scSemiCluster <- list()
  Total_Time_scSemiCluster <- list()
  
  for (i in c(1:n_folds)) {
    train_data <- Data[, Train_Idx[[i]]]
    train_label <- Labels[Train_Idx[[i]]]
    # train_set <- .generate_ref(train_data, train_label)
    test_set <- Data[, Test_Idx[[i]]]
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
      start_time <- Sys.time()
      # scSemiCluster = scSemiCluster(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
      #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
      #                   Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    } else {
      start_time <- Sys.time()
      pred.scSemiCluster <- func_scSemiCluster(test_set, train_data, train_label, k=k)
      label.scSemiCluster <- as.character(pred.scSemiCluster)
      end_time <- Sys.time()
    }
    Total_Time_scSemiCluster[i] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_scSemiCluster[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scSemiCluster[i] <- list(label.scSemiCluster)
  }
  True_Labels_scSemiCluster <- as.vector(unlist(True_Labels_scSemiCluster))
  Pred_Labels_scSemiCluster <- as.vector(unlist(Pred_Labels_scSemiCluster))
  Total_Time_scSemiCluster <- as.vector(unlist(Total_Time_scSemiCluster))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    write.csv(
      True_Labels_scSemiCluster,
      paste('scSemiCluster_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_scSemiCluster,
      paste('scSemiCluster_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Total_Time_scSemiCluster,
      paste('scSemiCluster_', NumGenes, '_Total_Time.csv', sep = ''),
      row.names = FALSE
    )
  } else{
    write.csv(True_Labels_scSemiCluster, 'scSemiCluster_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_scSemiCluster, 'scSemiCluster_Pred_Labels.csv', row.names = FALSE)
    write.csv(Total_Time_scSemiCluster, 'scSemiCluster_Total_Time.csv', row.names = FALSE)
  }
}


run_GSVAanno<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,pred_cluster_id,k=0,
                       GeneOrderPath = NULL){
  "
  run GSVAanno
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  NumGenes = NULL
  library(stringr)
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[, Cells_to_Keep]
  colnames(Data) <- str_replace_all(colnames(Data), '_', '.')
  colnames(Data) <- str_replace_all(colnames(Data), '-', '.')
  Labels <- Labels[Cells_to_Keep]
  pred_cluster_id <- pred_cluster_id[Cells_to_Keep, ]
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               GSVAanno                                     #
  #############################################################################
  source('/local/wzk/data/scmagic/code/source/scMAGIC_script/methods_functions.R')
  
  True_Labels_GSVAanno <- list()
  Pred_Labels_GSVAanno <- list()
  Total_Time_GSVAanno <- list()
  
  Train_Idx = sample(1:length(colnames(Data)), round(length(colnames(Data))/2)+1)
  Test_Idx = (1:length(colnames(Data)))[!((1:length(colnames(Data))) %in% Train_Idx)]
  
  
  train_data <- Data[, Train_Idx]
  train_label <- Labels[Train_Idx]
  # train_set <- .generate_ref(train_data, train_label)
  test_set <- Data[, Test_Idx]
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    start_time <- Sys.time()
    # GSVAanno = GSVAanno(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
    #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
    #                   Labels[Train_Idx[[i]]], numCores = 1)
    end_time <- Sys.time()
  } else {
    start_time <- Sys.time()
    sub_pred_cluster_id <- pred_cluster_id[Test_Idx, ]
    pred.GSVAanno <- func_GSVAanno(test_set, train_data, train_label, sub_pred_cluster_id, Labels[Test_Idx])
    label.GSVAanno <- as.character(pred.GSVAanno[[1]])
    end_time <- Sys.time()
  }
  Total_Time_GSVAanno[1] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  
  True_Labels_GSVAanno[1] <- list(pred.GSVAanno[[2]])
  Pred_Labels_GSVAanno[1] <- list(label.GSVAanno)
  
  True_Labels_GSVAanno <- as.vector(unlist(True_Labels_GSVAanno))
  Pred_Labels_GSVAanno <- as.vector(unlist(Pred_Labels_GSVAanno))
  Total_Time_GSVAanno <- as.vector(unlist(Total_Time_GSVAanno))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    write.csv(
      True_Labels_GSVAanno,
      paste('GSVAanno_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_GSVAanno,
      paste('GSVAanno_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Total_Time_GSVAanno,
      paste('GSVAanno_', NumGenes, '_Total_Time.csv', sep = ''),
      row.names = FALSE
    )
  } else{
    write.csv(True_Labels_GSVAanno, 'GSVAanno_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_GSVAanno, 'GSVAanno_Pred_Labels.csv', row.names = FALSE)
    write.csv(Total_Time_GSVAanno, 'GSVAanno_Total_Time.csv', row.names = FALSE)
  }
}


run_Seurat<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,k=c(30, 5, 200, 30),
                     GeneOrderPath = NULL){
  "
  run Seurat
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.tsv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.tsv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  NumGenes = NULL
  library(stringr)
  INPUT <- readRDS(DataPath)
  Data <- INPUT$data.filter
  Labels <- INPUT$label.filter
  Labels <- Labels[, 1]
  # Data <- read.delim(DataPath,row.names = 1)
  # Labels <- as.matrix(read.delim(LabelsPath, row.names = 1))
  load(CV_RDataPath)
  Data <- Data[, Cells_to_Keep]
  colnames(Data) <- str_replace_all(colnames(Data), '_', '.')
  colnames(Data) <- str_replace_all(colnames(Data), '-', '.')
  Labels <- Labels[Cells_to_Keep]
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               Seurat                                     #
  #############################################################################
  source('/local/wsy/method/methods_functions.R')
  
  True_Labels_seurat <- list()
  Pred_Labels_seurat <- list()
  Total_Time_seurat <- list()
  
  for (i in c(1:n_folds)) {
    train_data <- Data[, Train_Idx[[i]]]
    train_label <- Labels[Train_Idx[[i]]]
    # train_set <- .generate_ref(train_data, train_label)
    test_set <- Data[, Test_Idx[[i]]]
    if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
      start_time <- Sys.time()
      # scSemiCluster = scSemiCluster(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
      #                   Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
      #                   Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    } else {
      start_time <- Sys.time()
      pred.seurat <- func_seurat(test_set, train_data, train_label, k=k)
      label.seurat <- as.character(pred.seurat)
      end_time <- Sys.time()
    }
    Total_Time_seurat[i] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
    
    True_Labels_seurat[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_seurat[i] <- list(label.seurat)
  }
  True_Labels_seurat <- as.vector(unlist(True_Labels_seurat))
  Pred_Labels_seurat <- as.vector(unlist(Pred_Labels_seurat))
  Total_Time_seurat <- as.vector(unlist(Total_Time_seurat))
  
  setwd(OutputDir)
  
  if (!is.null(GeneOrderPath) & !is.null(NumGenes)) {
    write.csv(
      True_Labels_seura,
      paste('seura_', NumGenes, '_True_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Pred_Labels_seura,
      paste('seura_', NumGenes, '_Pred_Labels.csv', sep = ''),
      row.names = FALSE
    )
    write.csv(
      Total_Time_seurat,
      paste('seura_', NumGenes, '_Total_Time.csv', sep = ''),
      row.names = FALSE
    )
  } else{
    write.csv(True_Labels_seurat, 'seurat_True_Labels.csv', row.names = FALSE)
    write.csv(Pred_Labels_seurat, 'seurat_Pred_Labels.csv', row.names = FALSE)
    write.csv(Total_Time_seurat, 'seurat_Total_Time.csv', row.names = FALSE)
  }
}


