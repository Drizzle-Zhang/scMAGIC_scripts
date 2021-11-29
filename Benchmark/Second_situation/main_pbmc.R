library(Seurat)
library(SeuratData)
data("pbmcsca")
source('/local/zy/my_git/scMAGIC_scripts/Benchmark/methods_functions.R')


### scMAGIC
########################################
method <- 'scMAGIC'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_scMAGIC(exp_sc_mat, ref.mtx, ref.labels, atlas='HCL')
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### sciBet
########################################
method <- 'sciBet'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_sciBet(exp_sc_mat, ref.mtx, ref.labels, k=600)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### singleCellNet
########################################
method <- 'singleCellNet'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
# pathout <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/singleCellNet_1/PBMC/'
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_singleCellNet(exp_sc_mat, ref.mtx, ref.labels,
                                        nTopGenes = 10, nTopGenePairs = 25, nTrees = 500)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### scPred
########################################
method <- 'scPred'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
# pathout <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/singleCellNet/PBMC/'
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_scPred(exp_sc_mat, ref.mtx, ref.labels,
                                        nfeatures = 1000, npcs = 50)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### CHETAH
########################################
method <- 'CHETAH'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
# pathout <- '/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/CHETAH/PBMC/'

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }

        # depth_query <- round(median(colSums(exp_sc_mat != 0)))
        depth_ref <- round(median(colSums(ref.mtx != 0)))
        pred.tags <- func_CHETAH(exp_sc_mat, ref.mtx, ref.labels,
                                 clust_method = 'complete', n_genes = depth_ref)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### scClassify
########################################
method <- 'scClassify'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_scClassify(exp_sc_mat, ref.mtx, ref.labels, tree = "HOPACH",
                                     algorithm = "KNN", k = 10, topN = 50)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### scmapcluster
########################################
method <- 'scmapcluster'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_scmapcluster(exp_sc_mat, ref.mtx, ref.labels, threshold = 0.2)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### scmapcell
########################################
method <- 'scmapcell'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_scmapcell(exp_sc_mat, ref.mtx, ref.labels, w = 1, threshold = 0.4)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


### singleR
########################################
method <- 'SingleR'
pathout <- paste0('/mdshare/node9/zy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_SingleR(exp_sc_mat, ref.mtx, ref.labels,
                                  gene_method='de', quantile.use=0.8)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################

### CALLR
########################################
method <- 'CALLR'
pathout <- paste0('/mdshare/node8/wzk/scmagic/code/source/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        print(paste0('===================  ', ref, ' -->> ', query, '  ==================='))
        t1 = proc.time()
        pred.tags <- func_CALLR(exp_sc_mat, ref.mtx, ref.labels, k=0.9)
        t2 = proc.time()
        print(t2-t1)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################

### scSemiCluster
########################################
method <- 'scSemiCluster'
pathout <- paste0('/mdshare/node8/wzk/scmagic/code/source/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        print(paste0('===================  ', ref, ' -->> ', query, '  ==================='))
        if (paste0(ref, '_', query) %in% c('10Xv2_10Xv3', '10Xv2_CEL-Seq2', '10Xv2_Drop-seq', '10Xv2_inDrops', '10Xv2_Seq-Well',
                                           '10Xv2_Smart-seq2', '10Xv3_10Xv2', '10Xv3_CEL-Seq2', '10Xv3_Drop-seq')){
            next
        }else if(paste0(ref, '_', query) == '10Xv3_inDrops'){
            t1 = proc.time()
            pred.tags <- func_scSemiCluster(exp_sc_mat, ref.mtx, ref.labels, k=c(1000,128,64,32))
            t2 = proc.time()
            print(t2-t1)
            saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
        }else{
            next
        }
        
    }
}
########################################

### ### scID
########################################
method <- 'scID'
pathout <- paste0('/mdshare/node8/wzk/scmagic/code/source/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_scID(exp_sc_mat, ref.mtx, ref.labels)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################

### GSVAanno
########################################
method <- 'GSVAanno'
pathout <- paste0('/mdshare/node8/wzk/scmagic/code/source/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')
vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        print(paste0('===================  ', ref, ' -->> ', query, '  ==================='))
        pred_cluster_id = readRDS(paste0('/mdshare/node8/wzk/scmagic/code/source/MAGIC_cluster/', ref, '_', query, '.Rdata'))
        t1 = proc.time()
        pred.tags <- func_GSVAanno(exp_sc_mat, ref.mtx, ref.labels, pred_cluster_id, label_sc)
        t2 = proc.time()
        print(t2-t1)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
        
    }
}
########################################

### SeuratR
########################################
method <- 'SeuratR'
pathout <- paste0('/mdshare/node9/wsy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_seurat(exp_sc_mat, ref.mtx, ref.labels, k = c(50, 8, 320, 48))
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################

### SVMR
########################################
method <- 'SVMR'
pathout <- paste0('/mdshare/node9/wsy/MAGIC/Benchmark/cross_platform/', method, '/PBMC/')

vec.ref <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'Seq-Well', 'Smart-seq2', 'inDrops')
vec.query <- c('10Xv2', '10Xv3', 'CEL-Seq2', 'Drop-seq', 'inDrops', 'Seq-Well', 'Smart-seq2')

for (ref in vec.ref) {
    if (ref == '10Xv2') {
        ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
        ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
    } else {
        if (ref == '10Xv3') {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
        } else {
            ref.mtx <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == ref])
            ref.labels <- as.character(pbmcsca$CellType)[pbmcsca$Method == ref]
        }
    }
    for (query in vec.query) {
        if (query == ref) {next()}
        if (query == '10Xv2') {
            exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v2)'])
            label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v2)']
        } else {
            if (query == '10Xv3') {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == '10x Chromium (v3)'])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == '10x Chromium (v3)']
            } else {
                exp_sc_mat <- as.matrix(pbmcsca@assays$RNA@counts[, pbmcsca$Method == query])
                label_sc <- as.character(pbmcsca$CellType)[pbmcsca$Method == query]
            }
        }
        pred.tags <- func_SVM(exp_sc_mat, ref.mtx, ref.labels)
        saveRDS(pred.tags, file = paste0(pathout, ref, '_', query, '.Rdata'))
    }
}
########################################


