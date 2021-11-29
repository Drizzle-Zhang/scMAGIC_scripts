library(ggplot2)
library(RColorBrewer)

file_2situ <- '/mdshare/node9/zy/MAGIC/evaluation/Second_situation.txt'
df_2situ <- read.delim(file_2situ, sep = '\t', stringsAsFactors = F)

path.fig <- '/mdshare/node9/zy/MAGIC/fig/'


# mean
vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')
df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_2situ[df_2situ$method == method,]
    df_sub_HD <- df_sub[df_sub$type == 'High-depth reference',]
    df_sub_LD <- df_sub[df_sub$type == 'Low-depth reference',]
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
                         data.frame(term = 'labeled-Accuracy', method = method, 
                                    value = mean(df_sub$value[df_sub$term == 'labeled-Accuracy']),
                                    dataset = 'Mean', type = 'Average'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'labeled-Balanced accuracy', method = method, 
                                    value = mean(df_sub$value[df_sub$term == 'labeled-Balanced accuracy']),
                                    dataset = 'Mean', type = 'Average'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub_HD$value[df_sub_HD$term == 'Accuracy']),
                                    dataset = 'Mean', type = 'High-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub_HD$value[df_sub_HD$term == 'Balanced accuracy']),
                                    dataset = 'Mean', type = 'High-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'labeled-Accuracy', method = method, 
                                    value = mean(df_sub_HD$value[df_sub_HD$term == 'labeled-Accuracy']),
                                    dataset = 'Mean', type = 'High-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'labeled-Balanced accuracy', method = method, 
                                    value = mean(df_sub_HD$value[df_sub_HD$term == 'labeled-Balanced accuracy']),
                                    dataset = 'Mean', type = 'High-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub_LD$value[df_sub_LD$term == 'Accuracy']),
                                    dataset = 'Mean', type = 'Low-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub_LD$value[df_sub_LD$term == 'Balanced accuracy']),
                                    dataset = 'Mean', type = 'Low-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'labeled-Accuracy', method = method, 
                                    value = mean(df_sub_LD$value[df_sub_LD$term == 'labeled-Accuracy']),
                                    dataset = 'Mean', type = 'Low-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'labeled-Balanced accuracy', method = method, 
                                    value = mean(df_sub_LD$value[df_sub_LD$term == 'labeled-Balanced accuracy']),
                                    dataset = 'Mean', type = 'Low-depth reference'))
    df_mean <- rbind(df_mean, df_sub_mean)
}

df_plot <- rbind(df_2situ, df_mean)
color_type <- rep('0', nrow(df_plot))
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy', 'Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value < 0.1)] <- '1'
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy', 'Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value > 0.85)] <- '1'
color_type[(df_plot$dataset %in% c('Mean Accuracy', 'Mean Balanced Accuracy',
                                   'Mean Labeled-Accuracy', 'Mean',
                                   'Mean Labeled-Balanced Accuracy')) & 
               (df_plot$value > 0.1) & (df_plot$value < 0.85)] <- '2'
df_plot$Color_Type <- color_type


# Accuracy
########################################
sort.datasets <- c(unique(df_2situ$dataset), 'Mean')
labels.datasets <- c("CL -> 10XV2, HPBMC", "CL -> 10XV3, HPBMC", "CL -> DR, HPBMC", 
                     "CL -> iD, HPBMC", "CL -> SW, HPBMC", "SM2 -> 10XV2, HPBMC",
                     "SM2 -> 10XV3, HPBMC", "SM2 -> DR, HPBMC", "SM2 -> iD, HPBMC",
                     "SM2 -> SW, HPBMC", "SM -> SM4 (VISp), MB", "SM -> SM4 (ALM), MB", 
                     "SM -> DR (Campbell), MB", "SM -> 10XV2, MB", "SM -> DR (Mizrak), MB",
                     "SW -> 10XV2, HPBMC", "SW -> 10XV3, HPBMC", "SW -> DR, HPBMC",
                     "SW -> iD, HPBMC", "SW -> CL, HPBMC", "SW -> SM2, HPBMC", 
                     "MW -> SM4 (VISp), MB", "MW -> SM4 (ALM), MB", "MW -> DR (Campbell), MB",
                     "MW -> 10XV2, MB", "MW -> DR (Mizrak), MB", "Mean")
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
                           levels = c('High-depth reference', 'Low-depth reference', 'Average'))

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

ggsave(filename = 'heatmap_2situ.png', 
       path = path.fig, plot = plot.heatmap.acc,
       units = 'cm', height = 20, width = 32)


# boxplot
df.box.acc <- df_plot.acc[df_plot.acc$Color_Type == '0' & df_plot.acc$term == 'Accuracy',]

plot.box.acc <- 
    ggplot(df.box.acc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('High-depth reference', 'Low-depth reference'),
                       values = c("#1F78B4", "#EE7700")) +
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
ggsave(filename = 'boxplot_2situ_accuracy.png', 
       path = path.fig, plot = plot.box.acc,
       units = 'cm', height = 9, width = 20)

df.box.bacc <- df_plot.acc[df_plot.acc$Color_Type == '0' & df_plot.acc$term == 'Balanced accuracy',]

plot.box.bacc <- 
    ggplot(df.box.bacc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('High-depth reference', 'Low-depth reference'),
                       values = c("#1F78B4", "#EE7700")) +
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
ggsave(filename = 'boxplot_2situ_balanceaccuracy.png', 
       path = path.fig, plot = plot.box.bacc,
       units = 'cm', height = 9, width = 20)
########################################



# Labeled Accuracy
########################################
df_plot.acc <- df_plot[df_plot$term %in% c('labeled-Accuracy', 'labeled-Balanced accuracy'),]
df.mean.acc <- df_plot.acc[df_plot.acc$type == 'Average' & df_plot.acc$term == 'labeled-Accuracy',]
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
                           levels = c('High-depth reference', 'Low-depth reference', 'Average'))

plot.heatmap.label <-
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

ggsave(filename = 'heatmap_2situ_label.png', 
       path = path.fig, plot = plot.heatmap.label,
       units = 'cm', height = 22, width = 32)

# boxplot
df.box.acc <- df_plot.acc[df_plot.acc$Color_Type == '0' & df_plot.acc$term == 'labeled-Accuracy',]

plot.box.acc <- 
    ggplot(df.box.acc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('High-depth reference', 'Low-depth reference'),
                       values = c("#1F78B4", "#EE7700")) +
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
ggsave(filename = 'boxplot_2situ_accuracy_label.png', 
       path = path.fig, plot = plot.box.acc,
       units = 'cm', height = 9, width = 20)

df.box.bacc <- df_plot.acc[df_plot.acc$Color_Type == '0' & df_plot.acc$term == 'labeled-Balanced accuracy',]

plot.box.bacc <- 
    ggplot(df.box.bacc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('High-depth reference', 'Low-depth reference'),
                       values = c("#1F78B4", "#EE7700")) +
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
ggsave(filename = 'boxplot_2situ_balanceaccuracy_label.png', 
       path = path.fig, plot = plot.box.bacc,
       units = 'cm', height = 9, width = 20)


df.box <- df_plot[df_plot$term %in% c('Percent of unassigned'),]
df.box$value_limit <- df.box$value
df.box$value_limit[df.box$value < exp(-3)] <- exp(-3)
df.box$value_limit[df.box$value > exp(3)] <- exp(3)
df.box$log_ratio <- log(df.box$value_limit)
df.box$log_ratio[df.box$method %in% c("CALLR", "SVM", "CaSTLe", "SingleR", 
                                      "sciBet", "seurat", "scSemiCluster")] <- -3
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
df.box$method <- factor(df.box$method, levels = (sort.methods), 
                             labels = (label.methods))

plot.box <- 
    ggplot(df.box, aes(x = method, y = log_ratio)) +
    geom_boxplot(width = 0.7, outlier.size = 1, color = "#1F78B4") +
    geom_hline(yintercept = 0, color = 'dimgray', size = 0.5, linetype = 2) +
    # scale_color_manual(breaks = c('Accuracy', 'Balanced Accuracy'),
    #                    values = c("#1F78B4", "#FF7F00")) +
    coord_flip() + 
    labs(title = "", y = 'log(Ratio of unassigned cells)', x = '') + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          axis.text.y = element_text(
              size = 11, color = "black", family = 'Arial',
              vjust = 0.5, hjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 9, family = 'Arial'),
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1))
ggsave(filename = 'boxplot_2situ_labeled.png', 
       path = path.fig, plot = plot.box,
       units = 'cm', height = 10, width = 13)


########################################

# cluster
########################################
# mean
vec.methods <- c('scMAGIC', 'scmapcluster', 'scmapcell', 'SingleR', 'scClassify', 'scPred',
                 'sciBet', 'singleCellNet', 'CHETAH', 'CALLR', 'scSemiCluster', 'scID', 
                 'CaSTLe', 'SVM', 'seurat')
vec.datasets <- unique(df_2situ$dataset)
df_cluster_acc <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    for (dataset in vec.datasets) {
        num_correct <- df_2situ[df_2situ$method == method & df_2situ$dataset == dataset & 
                                    df_2situ$term == 'Number of correctly labeled clusters', 'value']
        num_total <- df_2situ[df_2situ$method == method & df_2situ$dataset == dataset & 
                                  df_2situ$term == 'Total number of clusters', 'value']
        type <- df_2situ[df_2situ$method == method & df_2situ$dataset == dataset & 
                             df_2situ$term == 'Total number of clusters', 'type']
        df_cluster_acc <- rbind(df_cluster_acc, data.frame(method = method, dataset = dataset,
                                                           value = num_correct/num_total,
                                                           type = type))
    }
}
df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_cluster_acc[df_cluster_acc$method == method,]
    df_sub_HD <- df_sub[df_sub$type == 'High-depth reference',]
    df_sub_LD <- df_sub[df_sub$type == 'Low-depth reference',]
    df_sub_mean <- data.frame(stringsAsFactors = F)
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub$value),
                                    dataset = 'Mean', type = 'Average'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub_HD$value),
                                    dataset = 'Mean', type = 'High-depth reference'))
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(method = method, value = mean(df_sub_LD$value),
                                    dataset = 'Mean', type = 'Low-depth reference'))
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
                           levels = c('High-depth reference', 'Low-depth reference', 'Average'))

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

ggsave(filename = 'heatmap_2situ_cluster.png', 
       path = path.fig, plot = plot.heatmap.acc,
       units = 'cm', height = 20, width = 20)

# boxplot
df.box.acc <- df_plot.acc[df_plot.acc$Color_Type == '0',]

plot.box.acc <- 
    ggplot(df.box.acc, aes(x = method, y = value, color = type)) +
    geom_boxplot(width = 0.7, outlier.size = 1) +
    scale_color_manual(breaks = c('High-depth reference', 'Low-depth reference'),
                       values = c("#1F78B4", "#EE7700")) +
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
ggsave(filename = 'boxplot_2situ_cluster.png', 
       path = path.fig, plot = plot.box.acc,
       units = 'cm', height = 9, width = 20)

########################################

