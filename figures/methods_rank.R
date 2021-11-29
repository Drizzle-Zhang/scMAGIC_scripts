library(ggplot2)
library(RColorBrewer)

path.fig <- '/mdshare/node9/zy/MAGIC/fig/'

file_1situ <- '/mdshare/node9/zy/MAGIC/evaluation/First_situation.txt'
df_1situ <- read.delim(file_1situ, sep = '\t', stringsAsFactors = F)

file_2situ <- '/mdshare/node9/zy/MAGIC/evaluation/Second_situation.txt'
df_2situ <- read.delim(file_2situ, sep = '\t', stringsAsFactors = F)

file_3situ <- '/mdshare/node9/zy/MAGIC/evaluation/Third_situation.txt'
df_3situ <- read.delim(file_3situ, sep = '\t', stringsAsFactors = F)

df_plot <- rbind(df_1situ, df_2situ, df_3situ)

sort.datasets <- unique(df_plot$dataset)
df.plot.rank <- data.frame(stringsAsFactors = F)
for (i in 1:length(sort.datasets)) {
    dataset <- sort.datasets[i]
    sub.res <- df_plot[df_plot$dataset == dataset, ]
    sub.acc <- sub.res[sub.res$term == 'Accuracy',]
    sub.acc$rank <- nrow(sub.acc) - ceiling(rank(sub.acc$value, ties.method = 'max')) + 1
    sub.bacc <- sub.res[sub.res$term == 'Balanced accuracy',]
    sub.bacc$rank <- nrow(sub.bacc) - ceiling(rank(sub.bacc$value, ties.method = 'max')) + 1
    sub.cluster <- sub.res[sub.res$term == 'Number of correctly labeled clusters',]
    sub.cluster$rank <- nrow(sub.cluster) - ceiling(rank(sub.cluster$value, ties.method = 'max')) + 1
    sub.plot <- rbind(sub.acc, sub.bacc, sub.cluster)
    sub.plot$data_idx <- rep(i, nrow(sub.plot))
    df.plot.rank <- rbind(df.plot.rank, sub.plot)
}

# accuracy
#########################################
df_plot_acc <- df.plot.rank[df.plot.rank$term == 'Accuracy',]
df_rank_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), mean)
colnames(df_rank_acc) <- c('method', 'mean')
df_rank_acc <- as.data.frame(df_rank_acc)

df_median_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), median)
colnames(df_median_acc) <- c('method', 'median')
df_median_acc <- as.data.frame(df_median_acc)

df_rank_acc <- merge(df_rank_acc, df_median_acc, by = 'method')
df_rank_acc$sum <- df_rank_acc$mean + df_rank_acc$median

sort.methods <- df_rank_acc$method[order(df_rank_acc$sum, decreasing = T)]
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
df_plot_acc$method <- factor(df_plot_acc$method, levels = (sort.methods), 
                             labels = (label.methods))
vec_color <- c(brewer.pal(12, 'Paired')[1:5], brewer.pal(12, 'Paired')[7:10], 
                   brewer.pal(12, 'Paired')[12], brewer.pal(12, 'Paired')[6])
getPalette <- colorRampPalette(vec_color)
vec.color <- getPalette(length(sort.methods))

plot.rank <- 
    ggplot(data = df_plot_acc, aes(x = data_idx, y = rank, 
                                   color = method, fill = method, 
                                   group = 1)) + 
    geom_area(alpha = 0.6) + 
    scale_color_manual(values = vec.color) + 
    scale_fill_manual(values = vec.color) + 
    facet_grid(method ~ ., scales = 'free', space = 'free') + 
    # ylim(c(0, 10)) + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          strip.text.y = element_text(
              size = 13, color = "black", family = 'Arial', angle = 0), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position = 'none')
ggsave(filename = 'Line_accuracy.png', 
       path = path.fig, plot = plot.rank,
       units = 'cm', height = 12, width = 20)


df_median_acc$method <- factor(df_median_acc$method, levels = rev(sort.methods))

plot.bar <- 
    ggplot(df_median_acc, aes(x = method, y = median, color = method, fill = method)) +
    geom_bar(stat = 'identity') +
    scale_color_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    scale_fill_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    labs(title = "", y = 'Median of ranks', x = '') + 
    coord_flip() + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'none',
          axis.text.x = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9, color = 'black', family = 'Arial'))
ggsave(filename = 'barplot_median_acc.png', 
       path = path.fig, plot = plot.bar,
       units = 'cm', height = 11, width = 3)

df_plot_acc_1situ <- df_plot_acc[df_plot_acc$type %in% 
                                     c("Cross-validation", "Cross-sample", "Cross-platform"),]
df_median_acc_1situ <- aggregate(df_plot_acc_1situ$rank, by = list(df_plot_acc_1situ$method), median)
colnames(df_median_acc_1situ) <- c('method', 'median')
df_median_acc_1situ <- as.data.frame(df_median_acc_1situ)


df_plot_acc_2situ <- df_plot_acc[df_plot_acc$type %in% 
                                     c("High-depth reference", "Low-depth reference", "Cross-species"),]
df_median_acc_2situ <- aggregate(df_plot_acc_2situ$rank, by = list(df_plot_acc_2situ$method), median)
colnames(df_median_acc_2situ) <- c('method', 'median')
df_median_acc_2situ <- as.data.frame(df_median_acc_2situ)
#########################################


# Balanced accuracy
#########################################
df_plot_acc <- df.plot.rank[df.plot.rank$term == 'Balanced accuracy',]
df_rank_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), mean)
colnames(df_rank_acc) <- c('method', 'mean')
df_rank_acc <- as.data.frame(df_rank_acc)

df_median_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), median)
colnames(df_median_acc) <- c('method', 'median')
df_median_acc <- as.data.frame(df_median_acc)

df_rank_acc <- merge(df_rank_acc, df_median_acc, by = 'method')
df_rank_acc$sum <- df_rank_acc$mean + df_rank_acc$median

sort.methods <- df_rank_acc$method[order(df_rank_acc$sum, decreasing = T)]
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
df_plot_acc$method <- factor(df_plot_acc$method, levels = (sort.methods), 
                             labels = (label.methods))
vec_color <- c(brewer.pal(12, 'Paired')[1:5], brewer.pal(12, 'Paired')[7:10], 
               brewer.pal(12, 'Paired')[12], brewer.pal(12, 'Paired')[6])
getPalette <- colorRampPalette(vec_color)
vec.color <- getPalette(length(sort.methods))

plot.rank.bacc <- 
    ggplot(data = df_plot_acc, aes(x = data_idx, y = rank, 
                                   color = method, fill = method, 
                                   group = 1)) + 
    geom_area(alpha = 0.6) + 
    scale_color_manual(values = vec.color) + 
    scale_fill_manual(values = vec.color) + 
    facet_grid(method ~ ., scales = 'free', space = 'free') + 
    # ylim(c(0, 10)) + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          strip.text.y = element_text(
              size = 13, color = "black", family = 'Arial', angle = 0), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position = 'none')
ggsave(filename = 'Line_balanced_accuracy.png', 
       path = path.fig, plot = plot.rank.bacc,
       units = 'cm', height = 12, width = 20)


df_median_acc$method <- factor(df_median_acc$method, levels = rev(sort.methods))
plot.bar <- 
    ggplot(df_median_acc, aes(x = method, y = median, color = method, fill = method)) +
    geom_bar(stat = 'identity') +
    scale_color_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    scale_fill_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    labs(title = "", y = 'Median of ranks', x = '') + 
    coord_flip() + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'none',
          axis.text.x = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9, color = 'black', family = 'Arial'))
ggsave(filename = 'barplot_median_bacc.png', 
       path = path.fig, plot = plot.bar,
       units = 'cm', height = 11, width = 3)
#########################################


# cluster
#########################################
df_plot_acc <- df.plot.rank[df.plot.rank$term == 'Number of correctly labeled clusters',]
df_rank_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), mean)
colnames(df_rank_acc) <- c('method', 'mean')
df_rank_acc <- as.data.frame(df_rank_acc)

df_median_acc <- aggregate(df_plot_acc$rank, by = list(df_plot_acc$method), median)
colnames(df_median_acc) <- c('method', 'median')
df_median_acc <- as.data.frame(df_median_acc)

df_rank_acc <- merge(df_rank_acc, df_median_acc, by = 'method')
df_rank_acc$sum <- df_rank_acc$mean + df_rank_acc$median

sort.methods <- df_rank_acc$method[order(df_rank_acc$median, decreasing = T)]
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
df_plot_acc$method <- factor(df_plot_acc$method, levels = (sort.methods), 
                             labels = (label.methods))
vec_color <- c(brewer.pal(12, 'Paired')[1:5], brewer.pal(12, 'Paired')[7:10], 
               brewer.pal(12, 'Paired')[12], brewer.pal(12, 'Paired')[6])
getPalette <- colorRampPalette(vec_color)
vec.color <- getPalette(length(sort.methods))

plot.rank.cluster <- 
    ggplot(data = df_plot_acc, aes(x = data_idx, y = rank, 
                                   color = method, fill = method, 
                                   group = 1)) + 
    geom_area(alpha = 0.6) + 
    scale_color_manual(values = vec.color) + 
    scale_fill_manual(values = vec.color) + 
    facet_grid(method ~ ., scales = 'free', space = 'free') + 
    # ylim(c(0, 10)) + 
    theme(panel.background = element_rect(color = 'transparent',
                                          fill = 'transparent'), 
          strip.text.y = element_text(
              size = 13, color = "black", family = 'Arial', angle = 0), 
          strip.background = element_rect(
              color = "transparent", fill = "transparent"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position = 'none')
ggsave(filename = 'Line_cluster_accuracy.png', 
       path = path.fig, plot = plot.rank.cluster,
       units = 'cm', height = 12, width = 20)


df_median_acc$method <- factor(df_median_acc$method, levels = rev(sort.methods))
plot.bar <- 
    ggplot(df_median_acc, aes(x = method, y = median, color = method, fill = method)) +
    geom_bar(stat = 'identity') +
    scale_color_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    scale_fill_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    labs(title = "", y = 'Median of ranks', x = '') + 
    coord_flip() + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = 'lightgray',
                                            size = 0.15, linetype = 'dashed'),
          legend.position = 'none',
          axis.text.x = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9, color = 'black', family = 'Arial'))
ggsave(filename = 'barplot_median_cacc.png', 
       path = path.fig, plot = plot.bar,
       units = 'cm', height = 11, width = 3)
#########################################

