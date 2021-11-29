
get_overlap_genes <- function(exp_sc_mat, exp_ref_mat) {
    exp_ref_mat <- as.data.frame(exp_ref_mat)
    exp_sc_mat <- as.data.frame(exp_sc_mat)
    # get overlap genes
    exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc <- rownames(exp_sc_mat)
    gene_ref <- rownames(exp_ref_mat)
    gene_over <- gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat <- exp_sc_mat[gene_over,]
    exp_ref_mat <- exp_ref_mat[gene_over,]

    out.overlap <- list()
    out.overlap$exp_sc_mat <- exp_sc_mat
    out.overlap$exp_ref_mat <- exp_ref_mat
    out.overlap$gene_over <- gene_over
    return(out.overlap)
}


transformHomoloGene <- function(exp_sc_mat, inTaxID = 9606, outTaxID = 10090) {
    library(homologene)
    genes.in <- rownames(exp_sc_mat)
    res.home <- homologene(genes.in, inTax = inTaxID, outTax = outTaxID)
    res.home <- res.home[!duplicated(res.home[, 1]),]
    res.home <- res.home[!duplicated(res.home[, 2]),]
    genes.out <- res.home[, 1]
    genes.homo <- res.home[, 2]
    exp.out <- exp_sc_mat[genes.out,]
    rownames(exp.out) <- genes.homo
    return(exp.out)
}



### scMAGIC
func_scMAGIC <- function(exp_sc_mat, ref.mtx, ref.labels, atlas, num_threads = 8) {
    library(scMAGIC)
    output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels,
                              atlas = atlas, num_threads = num_threads)
    pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
}

### sciBet
func_sciBet <- function(exp_sc_mat, ref.mtx, ref.labels, k=600) {
    suppressMessages(library(tidyverse))
    suppressMessages(library(scibet))
    suppressMessages(library(viridis))
    suppressMessages(library(ggsci))
    train_set <- as.data.frame(t(ref.mtx))/1.0
    train_set$label <- ref.labels
    test_set <- as.data.frame(t(exp_sc_mat))/1.0
    sciBet <- SciBet(train_set, test_set, k=k)
    return(sciBet)
}

### singleCellNet
func_singleCellNet <- function(exp_sc_mat, ref.mtx, ref.labels,
                               nTopGenes = 10, nTopGenePairs = 25, nTrees = 1000) {
    library(singleCellNet)
    library(dplyr)
    out <- get_overlap_genes(exp_sc_mat, ref.mtx)
    train_set <- as.matrix(out$exp_ref_mat)
    LabelsTrain <- data.frame(Annotation = ref.labels, row.names = colnames(ref.mtx))
    test_set <- as.matrix(out$exp_sc_mat)
    class_info <- scn_train(stTrain = LabelsTrain, expTrain = train_set, dLevel = "Annotation",
                            nTopGenes = nTopGenes, nTopGenePairs = nTopGenePairs, nTrees = nTrees)
    classRes <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_set, nrand = 50)
    classRes <- classRes[, colnames(test_set)]
    tags.singleCellNet <- rownames(classRes)[apply(classRes,2,which.max)]
    tags.singleCellNet[tags.singleCellNet == 'rand'] <- 'Unassigned'
    return(tags.singleCellNet)
}

### SingleR
func_SingleR <- function(exp_sc_mat, ref.mtx, ref.labels, gene_method='de', quantile.use=0.8) {
    library(SingleR)
    train_set <- as.matrix(ref.mtx)
    test_set <- as.matrix(exp_sc_mat)
    singler = SingleR(sc_data = test_set, ref_data = train_set, types = ref.labels,
                      genes = gene_method, quantile.use = quantile.use,
                      numCores = 2)
    pred.singleR <- singler$labels
    return(pred.singleR)
}

### CALLR
##########################
"dist2" = function( x, c = NA ) {

    # set the parameters for x
    if(is.na(c)) {
        c = x
    }

    # compute the dimension
    n1 = nrow(x)
    d1 = ncol(x)
    n2 = nrow(c)
    d2 = ncol(c)
    if(d1!=d2) {
        stop("Data dimension does not match dimension of centres.")
    }

    # compute the distance
    dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) +
        (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) -
        2 * (x%*%t(c))

    return(dist)

}

#the function to project a vector y to simplex
"ps" = function(y){


    z=sort(y)
    n=length(y)
    t=seq(n)
    for(i in 1:n-1){
        t[i]=(sum(z[(i+1):n])-1)/(n-i)
    }
    t[n]=z[n]-1
    tt=n-length(which(t<z))

    if(tt==0){
        tt=(sum(y)-1)/n
    }else{
        tt=t[tt]
    }

    x=y-tt
    x[x<0]<-0
    return(x)
}

#this is the function to get a training set matrix T from a known true label vector x
#r is the ratio of training set, can be set around 5%

"rand" = function(x,r){

    aaa=c()
    for (i in 1:length(unique(x))){
        aaa=c(aaa,sample(which(x==sort(unique(x))[i]),
                         max(ceiling(length(which(x==sort(unique(x))[i]))*r),2),replace=FALSE))
    }
    F=matrix(0,nrow=2,ncol=length(aaa))
    F[1,]=aaa
    F[2,]=x[aaa]
    return(F)
}



#input: expression matrix X, parameter u, training set matrix T

"callr"=function(X,u,T){

    N<-dim(X)[1]
    K=max(T[2,])

    Xt<-X[T[1,],]

    #do logistic regression on all samples to obtain the initial label matrix U
    logistic<-glmnet::glmnet(Xt,T[2,],family=c("multinomial"))

    b=matrix(0,nrow=K,ncol=dim(X)[2])
    for(i in 1:K){
        b[i,]=logistic[["beta"]][[i]][,length(logistic[["df"]])]
    }
    a=logistic[["a0"]][,length(logistic[["df"]])]

    U_hat=exp(t(t(X%*%t(b))+a))
    U_hat=U_hat/apply(U_hat,1,sum)
    U_hat=log(U_hat)

    for( i in 1:dim(T)[2]){
        U_hat[T[1,i],] = t(rep(-10000,K))
        U_hat[T[1,i],T[2,i]] = 0
    }
    U_hat[is.nan(U_hat)]=-10000
    U_hat[is.infinite(U_hat)]=-10000


    #Compute Laplacian Matrixx L



    D=dist2(X)

    D_sort = t(apply(D,MARGIN=2,FUN=sort));
    n = 17;


    mu=apply(D_sort[,2:n+1],MARGIN=1,FUN=mean)
    e=(matrix(rep(mu,N),N,N)+t(matrix(rep(mu,N),N,N)))/2

    W=exp(-D/(2*e^2))/(sqrt(2*pi)*e);


    for(i in 1:N){
        for(j in 1:N){
            if(D[i,j]>D_sort[i,n+1]&&D[i,j]>D_sort[j,n+1]){
                W[i,j] = 0;
            }else{
            }
        }
    }

    Degree = diag((rowSums(W))^(-0.5));
    L=diag(seq(1,1,length=N))-(Degree)%*%W%*%(Degree)


    #Some defaults and starting value

    dt = 0.005
    NS = 3


    U = matrix(data = runif(N*K,min=0,max=1),nrow=N,ncol=K)
    for(i in 1:N){
        U[i,]=ps(U[i,])
    }


    U[T[1,],] = U_hat[T[1,],]



    U0 = matrix(data = 0,nrow=N,ncol=K)
    U_old = matrix(data = 0,nrow=N,ncol=K)

    #The main iteration

    while (max(abs(U_old-U))>0.000001){
        U_old=U

        #MBO algorithm
        while(max(abs(U0-U))>0.000001){
            U0=U

            for(i in 1:NS){
                U = U+(dt/NS)*(-L%*%U+u*U_hat)
            }

            for(i in 1:N){
                U[i,] = ps(U[i,])

            }
        }
        #end of MBO algorithm

        for(i in 1:N){
            U[i,which.max(U[i,])]=1
            U[i,][U[i,]<1]=0
        }

        for( i in 1:dim(T)[2]){
            U[T[1,i],] = t(seq(0,0,length=K))
            U[T[1,i],T[2,i]] = 1
        }


        lab=seq(N)
        for(i in 1:N){
            lab[i]=which(U[i,]==1)
        }

        logistic<-glmnet::glmnet(X,lab,family=c("multinomial"))
        b=matrix(0,nrow=K,ncol=dim(X)[2])
        for(i in 1:K){
            b[i,]=logistic[["beta"]][[i]][,length(logistic[["df"]])]
        }
        a=logistic[["a0"]][,length(logistic[["df"]])]

        U_hat=exp(t(t(X%*%t(b))+a))
        U_hat=U_hat/apply(U_hat,1,sum)

        U_hat=log(U_hat)

        for( i in 1:dim(T)[2]){
            U_hat[T[1,i],] = t(seq(-10000,-10000,length=K))
            U_hat[T[1,i],T[2,i]] = 0
        }
        U_hat[is.nan(U_hat)]=-10000
        U_hat[is.infinite(U_hat)]=-10000

    }
    #end of the main iteration


    return(lab)
    #return the estimated labels
}


preprocess=function(X){
    m=dim(X)[2]
    for(i in 1:m){
        a=X[,i]
        a=as.matrix(a)[-which(a==0),]
        X[,i]=X[,i]/exp(mean(log(a)))
    }
    X=as.matrix(X)
    return(X)
}


generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM',
                         refnames = FALSE ){
    M <- M
    # print(M)
    min_cell <- min_cell
    refnames <- refnames
    exp_sc_mat <- exp_sc_mat
    TAG <- TAG
    NewRef <- c()
    TAG[, 2] <- as.character(TAG[, 2])
    if (refnames == FALSE) {
        refnames <- names(table(TAG[, 2]))
    }
    else{
        refnames <- refnames
    }
    outnames <- c()
    for (one in refnames) {
        this_col <- which(TAG[, 2] == one)
        if (length(this_col) >= min_cell) {
            outnames <- c(outnames, one)
            if (length(this_col) > 1) {
                if (M == 'SUM') {
                    this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
                } else{
                    this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
                }
            }
            else{
                this_new_ref <- exp_sc_mat[, this_col]
            }
            NewRef <- cbind(NewRef, this_new_ref)
        }
    }
    if (is.null(dim(NewRef))) {
        return(NULL)
    }
    rownames(NewRef) <- rownames(exp_sc_mat)
    colnames(NewRef) <- outnames
    return(NewRef)
}


.get_high_variance_genes <- function(exp_ref_mat, num.genes = 2000, type_ref = 'sum-counts') {
    if (type_ref %in% c('sum-counts', 'sc-counts')) {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref@assays$RNA@data <- exp_ref_mat
    }
    seurat.Ref <- FindVariableFeatures(
        seurat.Ref,
        selection.method = "vst",
        nfeatures = num.genes,
        verbose = F
    )
    return(VariableFeatures(seurat.Ref))

}

#########################
func_CALLR <- function(exp_sc_mat, ref.mtx, label_panc, u = 0.3) {
    library(Seurat)
    library(glmnet)
    library(Matrix)
    # preprocessing
    cell_num = 1000
    cells.sample = colnames(ref.mtx)
    if (length(cells.sample) > cell_num){
        cells.sample = sample(cells.sample, cell_num)
        ref.mtx = ref.mtx[, cells.sample]
        ref.labels = label_panc[cells.sample,]
    }
    if (length(colnames(exp_sc_mat)) > cell_num){
        test.sample = sample(colnames(exp_sc_mat), cell_num)
        exp_sc_mat = exp_sc_mat[, test.sample]
        label_sc = label_panc[test.sample, ]
    }

    TAG = label_panc
    TAG$id = rownames(TAG)
    TAG = TAG[colnames(ref.mtx), c(2, 1)]
    ref_gene_df = generate_ref(ref.mtx, TAG)
    hvg = .get_high_variance_genes(ref_gene_df, num.genes = 2000)
    mix_mat = cbind(ref.mtx[hvg, ], exp_sc_mat[hvg, ])

    type2num = 1:length(unique(ref.labels))
    names(type2num) = unique(ref.labels)
    ref_T = t(matrix(c(1:length(ref.labels),type2num[ref.labels]), ncol = 2))

    mix_mat = preprocess(mix_mat)
    output=callr(t(mix_mat), p, ref_T)
    num2type = names(type2num)
    names(num2type) = type2num
    output_type = num2type[output[1001:2000]]
    return(output_type)
}

### scClassify
func_scClassify <- function(exp_sc_mat, ref.mtx, ref.labels, tree = "HOPACH",
                            algorithm = "KNN", k = 10, topN = 50) {
    library("scClassify")
    library(Matrix)
    exprsMat_train <- as(as.matrix(log1p(ref.mtx)), "dgCMatrix")
    exp_sc_mat <- as(as.matrix(log1p(exp_sc_mat)), "dgCMatrix")
    scClassify_res <- scClassify(exprsMat_train = exprsMat_train,
                                 cellTypes_train = ref.labels,
                                 exprsMat_test = list(one = exp_sc_mat),
                                 tree = tree,
                                 algorithm = algorithm,
                                 k = k,
                                 topN = topN,
                                 returnList = FALSE,
                                 verbose = FALSE)
    pred.scClassify <- scClassify_res$testRes$one[[1]]$predRes
    return(pred.scClassify)
}

### scPred
func_scPred <- function(exp_sc_mat, ref.mtx, ref.labels, nfeatures = 1000, npcs = 50) {
    library("scPred")
    library("Seurat")
    library("magrittr")
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
    return(pred.scPred)
}

### CHETAH
func_CHETAH <- function(exp_sc_mat, ref.mtx, ref.labels,
                        clust_method = 'complete', n_genes = 0.8) {
    library(CHETAH)
    library(SingleCellExperiment)
    ref.labels[ref.labels == 'Unassigned'] <- 'unassigned'
    sce <- SingleCellExperiment(assays = list(counts = ref.mtx),
                                colData = data.frame(celltypes = ref.labels))
    sce_test <- SingleCellExperiment(list(counts = exp_sc_mat))

    # depth_query <- round(median(colSums(exp_sc_mat != 0)))
    depth_ref <- round(median(colSums(ref.mtx != 0)))
    n_genes <- n_genes * depth_ref
    sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce,
                                 clust_method = clust_method, n_genes = n_genes)
    pred.CHETAH <- sce_test$celltype_CHETAH
    pred.CHETAH[!(pred.CHETAH %in% setdiff(unique(ref.labels), 'unassigned'))] <- 'Unassigned'
    return(pred.CHETAH)
}

### scmap-cluster
func_scmapcluster <- function(exp_sc_mat, ref.mtx, ref.labels, threshold = 0.2) {
    library(scmap)
    library(SingleCellExperiment)
    out <- get_overlap_genes(exp_sc_mat, ref.mtx)
    train_set <- as.matrix(out$exp_ref_mat)
    test_set <- as.matrix(out$exp_sc_mat)
    sce <- SingleCellExperiment(list(normcounts = train_set),
                                colData = data.frame(cell_type1 = ref.labels))
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- selectFeatures(sce, suppress_plot = TRUE)

    sce_test <- SingleCellExperiment(list(normcounts = test_set))
    logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
    rowData(sce_test)$feature_symbol <- rownames(sce_test)
    sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData

    # scmap-cluster
    sce <- indexCluster(sce)
    scmapCluster_results <-
        scmapCluster(projection = sce_test,
                     index_list = list(metadata(sce)$scmap_cluster_index),
                     threshold = threshold)
    pred.scmap.cluster <- scmapCluster_results$combined_labs
    return(pred.scmap.cluster)
}

### scmap-cell
func_scmapcell <- function(exp_sc_mat, ref.mtx, ref.labels, w = 1, threshold = 0.4) {
    library(scmap)
    library(SingleCellExperiment)
    out <- get_overlap_genes(exp_sc_mat, ref.mtx)
    train_set <- as.matrix(out$exp_ref_mat)
    test_set <- as.matrix(out$exp_sc_mat)
    sce <- SingleCellExperiment(list(normcounts = train_set),
                                colData = data.frame(cell_type1 = ref.labels))
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- selectFeatures(sce, suppress_plot = TRUE)

    sce_test <- SingleCellExperiment(list(normcounts = test_set))
    logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
    rowData(sce_test)$feature_symbol <- rownames(sce_test)
    sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData

    # scmap-cluster
    set.seed(1)
    sce <- indexCell(sce)
    scmapCell_results <- scmapCell(sce_test,list(metadata(sce)$scmap_cell_index))
    scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)),
                                            w=w, threshold = threshold)
    pred.scmap.cell <- scmapCell_clusters$combined_labs
    return(pred.scmap.cell)
}

### scID
func_scID <- function(exp_sc_mat, ref.mtx, ref.labels, logFC = 0.4) {
    library(scID)
    library(Seurat)
    Train_Labels <- list(ref.labels)
    names(Train_Labels[[1]]) <- colnames(ref.mtx)
    scID_output <- scid_multiclass(exp_sc_mat, ref.mtx, Train_Labels[[1]],
                                   logFC = logFC)
    pred.scID <- scID_output$labels
    return(pred.scID)
}

### CaSTLe
func_CaSTLe <- function(exp_sc_mat, ref.mtx, ref.labels, nFeatures) {
    # 1. Load datasets
    ds1 = exp_sc_mat
    ds2 = ref.mtx
    sourceCellTypes = as.factor(ref.labels)

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
    pred_tags <- rownames(targetClassification)[apply(targetClassification,2,which.max)]
    return(pred_tags)
}


########Seurat
###################
func_seurat <- function(exp_sc_mat, ref.mtx, ref.labels, k = c(50, 8, 320, 48)){
    library(Seurat)
    get_overlap_genes <- function(exp_sc_mat, exp_ref_mat) {
        exp_ref_mat <- as.data.frame(exp_ref_mat)
        exp_sc_mat <- as.data.frame(exp_sc_mat)
        # get overlap genes
        exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)),]
        exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)),]
        gene_sc <- rownames(exp_sc_mat)
        gene_ref <- rownames(exp_ref_mat)
        gene_over <- gene_sc[which(gene_sc %in% gene_ref)]
        exp_sc_mat <- exp_sc_mat[gene_over,]
        exp_ref_mat <- exp_ref_mat[gene_over,]

        out.overlap <- list()
        out.overlap$exp_sc_mat <- exp_sc_mat
        out.overlap$exp_ref_mat <- exp_ref_mat
        out.overlap$gene_over <- gene_over
        return(out.overlap)

    }
    o_list = get_overlap_genes(exp_sc_mat, ref.mtx)
    exp_sc_mat = o_list[[1]]
    ref.mtx = o_list[[2]]
    ref<-CreateSeuratObject(ref.mtx)
    ref<-NormalizeData(ref,verbose = FALSE)
    ref<-FindVariableFeatures(ref,verbose = FALSE)
    ref<-ScaleData(ref)

    sample<-CreateSeuratObject(exp_sc_mat,verbose = FALSE)
    sample<-NormalizeData(sample,verbose = FALSE)
    sample<-FindVariableFeatures(sample,verbose = FALSE)
    sample<-ScaleData(sample)
    a = k[1]
    b = k[2]
    c = k[3]
    d = k[4]
    anchors <- FindTransferAnchors(reference = ref, query = sample, npcs = a, k.anchor = b, k.filter = c, k.score = d, verbose = FALSE)
    ref$annotation=ref.labels
    predictions <- TransferData(
        anchorset = anchors,
        refdata = ref$annotation,verbose = FALSE)
    PredictResult<-predictions[,1]
    return(PredictResult)
}
##########################################


########SVM
####################
func_SVM <-function(exp_sc_mat, ref.mtx, ref.labels, num_genes=2000, threshold=0.1){
    library(Seurat)
    library(stringr)
    colnames(ref.mtx) <- str_replace_all(colnames(ref.mtx), '_', '.')
    colnames(ref.mtx) <- str_replace_all(colnames(ref.mtx), '-', '.')
    colnames(exp_sc_mat) <- str_replace_all(colnames(exp_sc_mat), '_', '.')
    colnames(exp_sc_mat) <- str_replace_all(colnames(exp_sc_mat), '-', '.')
    rownames(ref.mtx) <- str_replace_all(rownames(ref.mtx), '_', '.')
    rownames(ref.mtx) <- str_replace_all(rownames(ref.mtx), '-', '.')
    rownames(exp_sc_mat) <- str_replace_all(rownames(exp_sc_mat), '_', '.')
    rownames(exp_sc_mat) <- str_replace_all(rownames(exp_sc_mat), '-', '.')

    get_overlap_genes <- function(exp_sc_mat, exp_ref_mat) {
        exp_ref_mat <- as.data.frame(exp_ref_mat)
        exp_sc_mat <- as.data.frame(exp_sc_mat)
        # get overlap genes
        exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)),]
        exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)),]
        gene_sc <- rownames(exp_sc_mat)
        gene_ref <- rownames(exp_ref_mat)
        gene_over <- gene_sc[which(gene_sc %in% gene_ref)]
        exp_sc_mat <- exp_sc_mat[gene_over,]
        exp_ref_mat <- exp_ref_mat[gene_over,]

        out.overlap <- list()
        out.overlap$exp_sc_mat <- exp_sc_mat
        out.overlap$exp_ref_mat <- exp_ref_mat
        out.overlap$gene_over <- gene_over
        return(out.overlap)

    }

    generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM',
                             refnames = FALSE ){
        M <- M
        # print(M)
        min_cell <- min_cell
        refnames <- refnames
        exp_sc_mat <- exp_sc_mat
        TAG <- TAG
        NewRef <- c()
        TAG[, 2] <- as.character(TAG[, 2])
        if (refnames == FALSE) {
            refnames <- names(table(TAG[, 2]))
        }
        else{
            refnames <- refnames
        }
        outnames <- c()
        for (one in refnames) {
            this_col <- which(TAG[, 2] == one)
            if (length(this_col) >= min_cell) {
                outnames <- c(outnames, one)
                if (length(this_col) > 1) {
                    if (M == 'SUM') {
                        this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
                    } else{
                        this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
                    }
                }
                else{
                    this_new_ref <- exp_sc_mat[, this_col]
                }
                NewRef <- cbind(NewRef, this_new_ref)
            }
        }
        if (is.null(dim(NewRef))) {
            return(NULL)
        }
        rownames(NewRef) <- rownames(exp_sc_mat)
        colnames(NewRef) <- outnames
        return(NewRef)
    }

    .get_high_variance_genes <- function(exp_ref_mat, num.genes = 2000, type_ref = 'sum-counts') {
        if (type_ref %in% c('sum-counts', 'sc-counts')) {
            seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
            seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                        scale.factor = 1e6, verbose = F)
        }
        if (type_ref %in% c('fpkm', 'tpm', 'rpkm')) {
            exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
            seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
            seurat.Ref@assays$RNA@data <- exp_ref_mat
        }
        seurat.Ref <- FindVariableFeatures(
            seurat.Ref,
            selection.method = "vst",
            nfeatures = num.genes,
            verbose = F
        )
        return(VariableFeatures(seurat.Ref))

    }

    o_list = get_overlap_genes(exp_sc_mat, ref.mtx)
    exp_sc_mat = o_list[[1]]
    ref.mtx = o_list[[2]]
    if (num_genes > 0) {
        TAG = data.frame(id = colnames(ref.mtx), label = ref.labels, stringsAsFactors = F)
        TAG = as.matrix(TAG)
        ref_gene_df = generate_ref(ref.mtx, TAG)
        hvg = .get_high_variance_genes(ref_gene_df, num.genes = num_genes)
        ref<-t(ref.mtx[hvg, ])
        test<-t(exp_sc_mat[hvg, ])
    } else {
        ref<-t(ref.mtx)
        test<-t(exp_sc_mat)
    }

    tmp_root <- "/mdshare/node9/zy/MAGIC/tools/SVM_tmp/"
    for (i in 1:10) {
        tmp_id <- sample(1:100000, 1)
        tmp_folder <- paste0(tmp_root, 'tmp_', tmp_id, '/')
        if (!dir.exists(tmp_folder)) {
            dir.create(tmp_folder)
            break
        }
    }
    # print(tmp_folder)

    setwd(tmp_folder)
    write.table(ref, file='ref.mtx.txt',sep = "\t",eol = "\n",
                na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
    write.table(ref.labels, file='ref.labels.txt',sep = "\t",eol = "\n",
                na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
    write.table(test, file='exp_sc_mat.txt',sep = "\t",eol = "\n",
                na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
    write.table(threshold, file='threshold.txt',sep = "\t",eol = "\n",
                na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
    system(command = '/local/zy/tools/anaconda3/bin/python3 /mdshare/node9/zy/MAGIC/tools/SVM/predict_svm.py', intern=F)
    a<-read.table('output.txt',header=F, sep="\t")
    PredictResult<-unlist(a)
    return(PredictResult)
}

### scSemiCluster
############################
func_scSemiCluster <- function(exp_sc_mat, ref.mtx, ref.labels, k=c(1000,128,64,32)) {
    o_list = get_overlap_genes(exp_sc_mat, ref.mtx)
    exp_sc_mat = o_list[[1]]
    ref.mtx = o_list[[2]]
    mix_mat = cbind(ref.mtx, exp_sc_mat)
    
    type2num = 1:length(unique(ref.labels))
    names(type2num) = unique(ref.labels)
    input_labels = c(type2num[ref.labels], rep(1, length(colnames(exp_sc_mat))))
    batch_labels = c(rep(0, length(ref.labels)), rep(1, length(colnames(exp_sc_mat))))
    tmp_dir = '~/data/scmagic/test_scSemiCluster/tmp/'
    suffix = as.character(unclass(Sys.time()))
    mat_file = paste0(tmp_dir, '/test_mat_', suffix, '.csv')
    labels_file = paste0(tmp_dir, '/test_labels_', suffix, '.csv')
    batch_file = paste0(tmp_dir, '/test_batch_', suffix, '.csv')
    para_file = paste0(tmp_dir, '/test_para_', suffix, '.csv')
    output_file = paste0(tmp_dir, '/test_pred_', suffix, '.csv')
    
    fwrite(t(mix_mat), file = mat_file, col.names = F, row.names = F)
    # write.table(t(mix_mat), file = mat_file, 
    #             row.names = F, col.names = F, sep = ',')
    write.table(input_labels, file = labels_file, 
                row.names = F, col.names = F, sep = ',')
    write.table(batch_labels, file = batch_file, 
                row.names = F, col.names = F, sep = ',')
    write.table(as.numeric(k), file = para_file, 
                row.names = F, col.names = F, sep = ',')
    system(paste0('python /local/wzk/data/scmagic/code/rscript/scSemiCluster_model.py ', mat_file, ' ', labels_file, ' ', batch_file, ' ', para_file, ' ', output_file),
           intern=FALSE,ignore.stdout=FALSE,ignore.stderr=FALSE,
           wait=TRUE,input=NULL,show.output.on.console=TRUE,minimized=FALSE,invisible=TRUE)
    pred_res = read.csv(output_file, header = F)
    pred_vec = as.character(pred_res[1,])
    num2type = names(type2num)
    names(num2type) = type2num
    num2type = c(num2type, 'unassigned')
    names(num2type)[length(num2type)] = '0'
    return(as.character(num2type[pred_vec]))
}

### GSVAanno
############################
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


func_GSVAanno <- function(exp_sc_mat, ref.mtx, ref.labels, pred_cluster_id, label_sc){
    library(Seurat)
    o_list = get_overlap_genes(exp_sc_mat, ref.mtx)
    exp_sc_mat = o_list[[1]]
    ref.mtx = o_list[[2]]
    
    ref.labels = gsub(' ', '.', ref.labels)
    label_sc = gsub(' ', '.', label_sc)
    ref.obj <- CreateSeuratObject(counts = ref.mtx)
    exp.obj <- CreateSeuratObject(counts = exp_sc_mat)
    ref.obj$label <- ref.labels
    exp.obj$cluster <- pred_cluster_id$cluster
    Idents(ref.obj) = ref.obj$label
    ref.obj.markers = FindAllMarkers(ref.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.4, verbose = FALSE)
    
    table_ref = table(ref.obj.markers$cluster)
    ref.id = unique(ref.labels)
    ref.gmt = data.frame()
    ref.gmt = rbind(ref.gmt, ref.id)
    ref.gmt = rbind(ref.gmt, ref.id)
    for (id in ref.id){
        genes = as.character(ref.obj.markers$gene[ref.obj.markers$cluster == id])
        genes_str = paste(genes, collapse = '\t')
        index = which(ref.id == id)
        ref.gmt[3, index] = genes_str
    }
    ref.gmt = t(ref.gmt)
    
    exp.obj = NormalizeData(exp.obj, verbose = F)
    exp.id = unique(levels(exp.obj$cluster))
    exp.data = as.data.frame(exp.obj@assays$RNA@data)
    exp.mat = data.frame()
    for (id in exp.id){
        if (length(which(pred_cluster_id$cluster == id)) == 1){
            exp.mat = rbind(exp.mat, exp.data[, (1:length(rownames(pred_cluster_id)))[pred_cluster_id$cluster == id]])
        }else{
            exp.mat = rbind(exp.mat, rowMeans(exp.data[, (1:length(rownames(pred_cluster_id)))[pred_cluster_id$cluster == id]]))
        }
    }
    exp.mat = t(exp.mat)
    exp.mat[is.na(exp.mat)] = 0
    rownames(exp.mat) = rownames(exp.data)
    colnames(exp.mat) = exp.id
    exp.mat = as.data.frame(exp.mat)
    exp.mat$genes = rownames(exp.mat)
    exp.mat = exp.mat[, c('genes', exp.id)]
    
    tmp_dir = '/mdshare/node8/wzk/scmagic/test_GSVA_anno/tmp/'
    suffix = as.character(unclass(Sys.time()))
    ref_file = paste0(tmp_dir, '/', suffix, '.gmt')
    mat_file = paste0(tmp_dir, '/', suffix, '.txt')
    fwrite(exp.mat, file = mat_file, sep = '\t', row.names = F, col.names = T, quote = F)
    fwrite(ref.gmt, file = ref_file, sep = '\t', row.names = F, col.names = F, quote = F)
    system(paste0('/usr/local/R/bin/Rscript /mdshare/node8/wzk/scmagic/code/scRNAseq_cell_cluster_labeling-master/bin/r_programs/obtains_GSVA_for_MatrixColumns.R -i ', mat_file, 
                  ' -t DGE -c ', ref_file, ' -o ', tmp_dir, ' -p ', suffix),
           intern=FALSE,ignore.stdout=FALSE,ignore.stderr=FALSE,
           wait=TRUE,input=NULL,show.output.on.console=TRUE,minimized=FALSE,invisible=TRUE)
    pred_df = read.table(paste0(tmp_dir, '/GSVA/', suffix, '.GSVA_final_label.tsv'), sep = '\t', header = F)
    pred_df$cluster = NA
    for (pos in 1:length(pred_df[,1])){
        pred_df[pos, 3] = substr(pred_df[pos, 1], 2, nchar(pred_df[pos, 1]))
    }
    pred_df[,3] = as.numeric(pred_df[,3])-2
    pred_df = pred_df[order(pred_df[,3]),]
    true_cluster_tag = cluster_tags(label_sc, pred_cluster_id$cluster)
    # table(true_cluster_tag, pred_df$V2)
    pred.tags = pred_df[,2]
    return(list(pred = pred.tags, true = true_cluster_tag))
}



