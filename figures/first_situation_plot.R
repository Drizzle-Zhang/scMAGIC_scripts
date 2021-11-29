library(ggplot2)
library(RColorBrewer)

file_1situ <- '/mdshare/node9/zy/MAGIC/evaluation/First_situation.txt'
df_1situ <- read.delim(file_1situ, sep = '\t', stringsAsFactors = F)
path.fig <- '/mdshare/node9/zy/MAGIC/fig/'

# mean
vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')
df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_1situ[df_1situ$method == method,]
    df_sub_CV <- df_sub[df_sub$type == 'Cross-validation',]
    df_sub_CS <- df_sub[df_sub$type == 'Cross-sample',]
    df_sub_CP <- df_sub[df_sub$type == 'Cross-platform',]
    df_sub_mean <- data.frame(stringsAsFactors = F)
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub$value[df_sub$term == 'Accuracy']),
                                    dataset = 'Mean', type = 'Average'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub$value[df_sub$term == 'Balanced accuracy']),
                                    dataset = 'Mean', type = 'Average'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub_CV$value[df_sub_CV$term == 'Accuracy']),
                                    dataset = 'Mean', type = 'Cross-validation'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub_CV$value[df_sub_CV$term == 'Balanced accuracy']),
                                    dataset = 'Mean', type = 'Cross-validation'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub_CS$value[df_sub_CS$term == 'Accuracy']),
                                    dataset = 'Mean', type = 'Cross-sample'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub_CS$value[df_sub_CS$term == 'Balanced accuracy']),
                                    dataset = 'Mean', type = 'Cross-sample'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub_CP$value[df_sub_CP$term == 'Accuracy']),
                                    dataset = 'Mean', type = 'Cross-platform'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub_CP$value[df_sub_CP$term == 'Balanced accuracy']),
                                    dataset = 'Mean', type = 'Cross-platform'))
    df_mean <- rbind(df_mean, df_sub_mean)
}

df_plot <- rbind(df_1situ, df_mean)
color_type <- rep('0', nrow(df_plot))
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy', 'Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value < 0.1)] <- '1'
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy','Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value > 0.85)] <- '1'
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy','Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value > 0.1) & (df_plot$value < 0.85)] <- '2'
df_plot$Color_Type <- color_type


# Accuracy
########################################
sort.datasets <- c(unique(df_1situ$dataset), 'Mean')
labels.datasets <- c("Baron, MP", "Baron, HP", "Muraro, HP", "Segerstolpe, HP", "Xin, HP", 
                     "Tian, CellBench 10X", "Tian, CellBench CL", "Tasic2015, MB",
                     "Campbell, MB", "TM, Mouse", 
                     "Tasic2018, MB(L1)", "Tasic2018, MB(L2)", "Tasic2018, MB(L3)",
                     "Zheng, HPBMC(sorted)", "Zhang, HPBMC(68K)", 
                     "Ding, HPBMC(10XV2)", "Ding, HPBMC(CL)", "Ding, HPBMC(DR)", 
                     "Ding, HPBMC(iD)", "Ding, HPBMC(SW)", "Ding, HPBMC(SM2)", 
                     "Tasic2018, MB(ALM->VISp)", "Tasic2018, MB(VISp->ALM)", 
                     "iD -> CL, HP", "iD -> SM, HP", "CL -> iD, HP", "CL -> SM, HP", "SM2 -> CL, HP",
                     "DR -> SM, MB", "DR -> FC, MB",
                     "10XV2 -> 10XV3, HPBMC", "10XV2 -> DR, HPBMC", "10XV2 -> iD, HPBMC",
                     "10XV2 -> SW, HPBMC", "10XV2 -> CL, HPBMC", "10XV2 -> SM2, HPBMC",
                     "10XV3 -> 10XV2, HPBMC", "10XV3 -> DR, HPBMC", "10XV3 -> iD, HPBMC",
                     "10XV3 -> SW, HPBMC", "10XV3 -> CL, HPBMC", "10XV3 -> SM2, HPBMC",
                     "DR -> 10XV2, HPBMC", "DR -> 10XV3, HPBMC", "DR -> iD, HPBMC",
                     "DR -> SW, HPBMC", "DR -> CL, HPBMC", "DR -> SM2, HPBMC",
                     "iD -> 10XV2, HPBMC", "iD -> 10XV3, HPBMC", "iD -> DR, HPBMC",
                     "iD -> SW, HPBMC", "iD -> CL, HPBMC", "iD -> SM2, HPBMC",
                     "CL -> SM2, HPBMC", "SM2 -> CL, HPBMC", "Mean")
df_plot.acc <- df_plot[df_plot$term %in% c('Accuracy', 'Balanced accuracy'),]
df.mean.acc <- df_plot.acc[df_plot.acc$type == 'Average' & df_plot.acc$term == 'Accuracy',]
sort.methods <- df.mean.acc$method[order(df.mean.acc$value, decreasing = F)]
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
df_plot.acc$method <- factor(df_plot.acc$method, levels = rev(sort.methods), 
                             labels = rev(label.methods))
df_plot.acc$dataset <- factor(df_plot.acc$dataset, levels = rev(sort.datasets), labels = rev(labels.datasets))
df_plot.acc$type <- factor(df_plot.acc$type, 
                           levels = c('Cross-validation', 'Cross-sample', 'Cross-platform', 'Average'))

plot.heatmap.acc <-
    ggplot(data = df_plot.acc, aes(x = method, y = dataset)) +
    geom_tile(aes(fill = value)) +
    facet_grid(type ~ term, scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(c(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:8], 'firebrick'))(100)) +
    scale_y_discrete(position = 'right') +
    labs(x = '', y = '', fill = '') +
    geom_text(aes(label = round(value, 2), color = Color_Type),
              family = "Arial", size = 3) +
    scale_color_manual(breaks = c('0', '1', '2'),
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(
              size = 11, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.x = element_text(
              size = 16, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.y = element_blank(),
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial'),
          legend.position = 'none') +
    guides(color = F)

ggsave(filename = 'heatmap_CV_CP.png', 
       path = path.fig, plot = plot.heatmap.acc,
       units = 'cm', height = 32.5, width = 32)

# boxplot
df.box.acc <- df_plot.acc[df_plot.acc$Color_Type == '0' & df_plot.acc$term == 'Accuracy',]

plot.box.acc <- 
    ggplot(df.box.acc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Cross-validation', 'Cross-sample', 'Cross-platform'),
                       values = c("#1F78B4", "#EE7700", "#B22222")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(ylim = c(0, 1)) + 
    # coord_flip() + 
    labs(title = "", y = '', x = '', color = 'Metrics', fill = '') + 
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'right',
          axis.text.y = element_text(size = 11, color = 'black', family = 'Arial'),
          axis.text.x = element_text(size = 14.5, color = 'black', family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.title = element_text(size = 13, family = 'Arial'),
          legend.title = element_text(
              size = 14, color = "black", family = 'Arial'),
          legend.text = element_text(size = 13, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_CV_CP_accuracy.png', 
       path = path.fig, plot = plot.box.acc,
       units = 'cm', height = 9, width = 19)

df.box.bacc <- df_plot.acc[df_plot.acc$Color_Type == '0' & df_plot.acc$term == 'Balanced accuracy',]

plot.box.bacc <- 
    ggplot(df.box.bacc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Cross-validation', 'Cross-sample', 'Cross-platform'),
                       values = c("#1F78B4", "#EE7700", "#B22222")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(ylim = c(0, 1)) + 
    # coord_flip() + 
    labs(title = "", y = '', x = '', color = 'Metrics', fill = '') + 
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'right',
          axis.text.y = element_text(size = 11, color = 'black', family = 'Arial'),
          axis.text.x = element_text(size = 14.5, color = 'black', family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.title = element_text(size = 13, family = 'Arial'),
          legend.title = element_text(
              size = 14, color = "black", family = 'Arial'),
          legend.text = element_text(size = 13, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_CV_CP_balanceaccuracy.png', 
       path = path.fig, plot = plot.box.bacc,
       units = 'cm', height = 9, width = 19)
########################################


# cluster
########################################
# mean
vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')
vec.datasets <- unique(df_1situ$dataset)

df_cluster_acc <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    for (dataset in vec.datasets) {
        num_correct <- df_1situ[df_1situ$method == method & df_1situ$dataset == dataset & 
                                    df_1situ$term == 'Number of correctly labeled clusters', 'value']
        num_total <- df_1situ[df_1situ$method == method & df_1situ$dataset == dataset & 
                                  df_1situ$term == 'Total number of clusters', 'value']
        type <- df_1situ[df_1situ$method == method & df_1situ$dataset == dataset & 
                             df_1situ$term == 'Total number of clusters', 'type']
        df_cluster_acc <- rbind(df_cluster_acc, data.frame(method = method, dataset = dataset,
                                                           value = num_correct/num_total,
                                                           type = type))
    }
}
df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_cluster_acc[df_cluster_acc$method == method,]
    df_sub_CV <- df_sub[df_sub$type == 'Cross-validation',]
    df_sub_CS <- df_sub[df_sub$type == 'Cross-sample',]
    df_sub_CP <- df_sub[df_sub$type == 'Cross-platform',]
    df_sub_mean <- data.frame(stringsAsFactors = F)
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub$value),
                                    dataset = 'Mean', type = 'Average'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub_CV$value),
                                    dataset = 'Mean', type = 'Cross-validation'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub_CS$value),
                                    dataset = 'Mean', type = 'Cross-sample'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub_CP$value),
                                    dataset = 'Mean', type = 'Cross-platform'))
    df_mean <- rbind(df_mean, df_sub_mean)
}

df_plot <- rbind(df_cluster_acc, df_mean)
color_type <- rep('0', nrow(df_plot))
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy', 'Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value < 0.1)] <- '1'
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy','Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value > 0.85)] <- '1'
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy','Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value > 0.1) & (df_plot$value < 0.85)] <- '2'
df_plot$Color_Type <- color_type

sort.datasets <- c(unique(df_1situ$dataset), 'Mean')
df_plot.acc <- df_plot
df.mean.acc <- df_plot.acc[df_plot.acc$type == 'Average',]
sort.methods <- df.mean.acc$method[order(df.mean.acc$value, decreasing = F)]
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
df_plot.acc$method <- factor(df_plot.acc$method, levels = rev(sort.methods), 
                             labels = rev(label.methods))
df_plot.acc$dataset <- factor(df_plot.acc$dataset, levels = rev(sort.datasets), labels = rev(labels.datasets))
df_plot.acc$type <- factor(df_plot.acc$type, 
                           levels = c('Cross-validation', 'Cross-sample', 'Cross-platform', 'Average'))

plot.heatmap.acc <-
    ggplot(data = df_plot.acc, aes(x = method, y = dataset)) +
    geom_tile(aes(fill = value)) +
    facet_grid(type ~ ., scales = 'free', space = 'free') +
    scale_fill_gradientn(
        colors = colorRampPalette(c(rev(brewer.pal(n = 9, name =  "RdYlBu"))[1:8], 'firebrick'))(100)) +
    scale_y_discrete(position = 'right') +
    labs(x = '', y = '', fill = '') +
    geom_text(aes(label = round(value, 2), color = Color_Type),
              family = "Arial", size = 3) +
    scale_color_manual(breaks = c('0', '1', '2'),
                       values = c('transparent', 'white', 'black')) +
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(
              size = 11, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.x = element_text(
              size = 16, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          strip.text.y = element_blank(),
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial'),
          legend.position = 'none') +
    guides(color = F)

ggsave(filename = 'heatmap_CV_CP_cluster.png', 
       path = path.fig, plot = plot.heatmap.acc,
       units = 'cm', height = 32.5, width = 20)

# boxplot
df.box.acc <- df_plot.acc[df_plot.acc$Color_Type == '0',]

plot.box.acc <- 
    ggplot(df.box.acc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('Cross-validation', 'Cross-sample', 'Cross-platform'),
                       values = c("#1F78B4", "#EE7700", "#B22222")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(ylim = c(0, 1)) + 
    # coord_flip() + 
    labs(title = "", y = '', x = '', color = 'Metrics', fill = '') + 
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'right',
          axis.text.y = element_text(size = 11, color = 'black', family = 'Arial'),
          axis.text.x = element_text(size = 14.5, color = 'black', family = 'Arial',
                                     angle = 315, vjust = 1, hjust = 0),
          axis.title = element_text(size = 13, family = 'Arial'),
          legend.title = element_text(
              size = 14, color = "black", family = 'Arial'),
          legend.text = element_text(size = 13, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_CV_CP_cluster.png', 
       path = path.fig, plot = plot.box.acc,
       units = 'cm', height = 9, width = 20)

########################################


