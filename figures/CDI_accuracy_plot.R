library(ggplot2)
library(RColorBrewer)

file_1situ <- '/mdshare/node9/zy/MAGIC/evaluation/First_situation.txt'
df_1situ <- read.delim(file_1situ, sep = '\t', stringsAsFactors = F)

file_2situ <- '/mdshare/node9/zy/MAGIC/evaluation/Second_situation.txt'
df_2situ <- read.delim(file_2situ, sep = '\t', stringsAsFactors = F)

file_3situ <- '/mdshare/node9/zy/MAGIC/evaluation/Third_situation.txt'
df_3situ <- read.delim(file_3situ, sep = '\t', stringsAsFactors = F)

df_plot <- rbind(df_1situ, df_2situ, df_3situ)
df_dataset <- 
    df_plot[(df_plot$term %in% c('F1_asw', 'Depth of reference', 'Cell proportion covered by reference')) & 
                (df_plot$method == "scMAGIC" ),]
df_dataset_index <- data.frame(stringsAsFactors = F)
for (dataset in unique(df_dataset$dataset)) {
    depth_ref <- df_dataset[(df_dataset$dataset == dataset)&
                                (df_dataset$term == 'Depth of reference'), 'value']
    if (depth_ref > 5000) {depth_ref <- 5000}
    if (depth_ref < 500) {depth_ref <- 500}
    df_dataset_index <- 
        rbind(df_dataset_index,
              data.frame(dataset = dataset, 
                         index_1 = df_dataset[(df_dataset$dataset == dataset)&
                                                  (df_dataset$term == 'F1_asw'), 'value'],
                         index_2 = depth_ref,
                         index_3 = df_dataset[(df_dataset$dataset == dataset)&
                                                  (df_dataset$term == 'Cell proportion covered by reference'), 'value']))
}

df_dataset_index$index_1_norm <- (df_dataset_index$index_1 - min(df_dataset_index$index_1)) / (max(df_dataset_index$index_1) - min(df_dataset_index$index_1))
df_dataset_index$index_2_norm <- (df_dataset_index$index_2 - min(df_dataset_index$index_2)) / (max(df_dataset_index$index_2) - min(df_dataset_index$index_2))
df_dataset_index$index_3_norm <- (df_dataset_index$index_3 - min(df_dataset_index$index_3)) / (max(df_dataset_index$index_3) - min(df_dataset_index$index_3))
df_dataset_index$combine_index <- 
    apply(df_dataset_index[, c('index_1', 'index_3')], 1, 
          function(x) {(((1 - x[1])^2 + (x[2])^2)/3)})
df_dataset_index$rank_combine_index <- rank(df_dataset_index$combine_index)

df_plot_index <- merge(df_plot, df_dataset_index[, c('dataset', 'combine_index', 
                                                     'rank_combine_index')], by = 'dataset')

## Accuracy
#########################################
df_plot_index_acc <- df_plot_index[df_plot_index$term == 'Accuracy',]

vec.methods <- unique(df_plot_index_acc$method)
df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_plot_index_acc[df_plot_index_acc$method == method,]
    df_sub_mean <- data.frame(stringsAsFactors = F)
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Accuracy', method = method, 
                                    value = mean(df_sub$value[df_sub$term == 'Accuracy']),
                                    dataset = 'Mean Accuracy', type = 'Average'))
    df_mean <- rbind(df_mean, df_sub_mean)
}
sort.methods <- rev(df_mean$method[order(df_mean$value, decreasing = F)])
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
vec_color <- rev(c(brewer.pal(12, 'Paired')[1:5], brewer.pal(12, 'Paired')[7:10], 
                   brewer.pal(12, 'Paired')[12], brewer.pal(12, 'Paired')[6]))
getPalette <- colorRampPalette(vec_color)
vec.color <- getPalette(length(sort.methods))
vec.shape <- c(19, 0, 2, 3, 1, 4, 23, 6, 10, 8, 12, 17, 15, 13, 11)

plot_index <- 
    ggplot(df_plot_index_acc, aes(x = rank_combine_index, y = value, color = method, shape = method)) + 
    geom_point(size = 1.2) + 
    geom_smooth(method = lm, formula = y ~ poly(x, 1), se = F, size = 0.8) + 
    labs(x = 'Ranks of Classification Difficulty Index (CDI)', y = 'Accuracy', color = 'Method', shape = 'Method') + 
    scale_color_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    scale_shape_manual(breaks = sort.methods, values = vec.shape, labels = label.methods) +
    theme_bw() + 
    theme(panel.background = element_rect(colour = 'black', fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12, color = "black", family = 'Arial'),
          axis.text.x = element_text(size = 10, color = "black", family = 'Arial'),
          axis.text.y = element_text(
              size = 10, color = "black", family = 'Arial'),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial')) + 
    guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

path.fig <- '/mdshare/node9/zy/MAGIC/fig/'
ggsave(filename = 'Difficulty_accuracy.png', 
       path = path.fig, plot = plot_index,
       units = 'cm', height = 11, width = 20)

vec_slope <- c()
for (method in sort.methods) {
    sub_df <- df_plot_index_acc[df_plot_index_acc$method == method, c("value", "rank_combine_index")]
    fit_summary <- summary(lm(value ~ rank_combine_index, data = sub_df))
    print(method)
    print(fit_summary)
    vec_slope <- c(vec_slope, fit_summary$coefficients['rank_combine_index', 'Estimate'])
}
names(vec_slope) <- label.methods
z_slope <- as.data.frame(scale(vec_slope))
df_z_slope <- data.frame(method = rownames(z_slope), slope = z_slope$V1)
df_z_slope$method <- factor(df_z_slope$method, levels = rev(label.methods))
# sort(vec_slope)
plot.bar <- 
    ggplot(df_z_slope, aes(x = method, y = slope, color = method, fill = method)) +
    geom_bar(stat = 'identity') +
    scale_color_manual(breaks = label.methods, values = vec.color, labels = label.methods) +
    scale_fill_manual(breaks = label.methods, values = vec.color, labels = label.methods) +
    labs(title = "", y = "Z-score of\n the fitted line's slope", x = '') + 
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
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 11, color = 'black', family = 'Arial'))
ggsave(filename = 'barplot_slope_acc.png', 
       path = path.fig, plot = plot.bar,
       units = 'cm', height = 11, width = 4)

#########################################

## balanced accuracy
#########################################
df_plot_index_acc <- df_plot_index[df_plot_index$term == 'Balanced accuracy',]

vec.methods <- unique(df_plot_index_acc$method)
df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_plot_index_acc[df_plot_index_acc$method == method,]
    df_sub_mean <- data.frame(stringsAsFactors = F)
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Balanced accuracy', method = method, 
                                    value = mean(df_sub$value[df_sub$term == 'Balanced accuracy']),
                                    dataset = 'Mean Accuracy', type = 'Average'))
    df_mean <- rbind(df_mean, df_sub_mean)
}
sort.methods <- rev(df_mean$method[order(df_mean$value, decreasing = F)])
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
vec_color <- rev(c(brewer.pal(12, 'Paired')[1:5], brewer.pal(12, 'Paired')[7:10], 
                   brewer.pal(12, 'Paired')[12], brewer.pal(12, 'Paired')[6]))
getPalette <- colorRampPalette(vec_color)
vec.color <- getPalette(length(sort.methods))
vec.shape <- c(19, 0, 2, 3, 1, 4, 23, 6, 10, 8, 12, 17, 15, 13, 11)

plot_index <- 
    ggplot(df_plot_index_acc, aes(x = rank_combine_index, y = value, color = method, shape = method)) + 
    geom_point(size = 1.2) + 
    geom_smooth(method = lm, formula = y ~ poly(x, 1), se = F, size = 0.8) + 
    labs(x = 'Ranks of Classification Difficulty Index (CDI)', y = 'Balanced accuracy', color = 'Method', shape = 'Method') + 
    scale_color_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    scale_shape_manual(breaks = sort.methods, values = vec.shape, labels = label.methods) +
    theme_bw() + 
    theme(panel.background = element_rect(colour = 'black', fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12, color = "black", family = 'Arial'),
          axis.text.x = element_text(size = 10, color = "black", family = 'Arial'),
          axis.text.y = element_text(
              size = 10, color = "black", family = 'Arial'),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial')) + 
    guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

path.fig <- '/mdshare/node9/zy/MAGIC/fig/'
ggsave(filename = 'Difficulty_balancedaccuracy.png', 
       path = path.fig, plot = plot_index,
       units = 'cm', height = 11, width = 20)

vec_slope <- c()
for (method in sort.methods) {
    sub_df <- df_plot_index_acc[df_plot_index_acc$method == method, c("value", "rank_combine_index")]
    # print(method)
    fit_summary <- summary(lm(value ~ rank_combine_index, data = sub_df))
    vec_slope <- c(vec_slope, fit_summary$coefficients['rank_combine_index', 'Estimate'])
}
names(vec_slope) <- label.methods
z_slope <- as.data.frame(scale(vec_slope))
df_z_slope <- data.frame(method = rownames(z_slope), slope = z_slope$V1)
df_z_slope$method <- factor(df_z_slope$method, levels = rev(label.methods))
# sort(vec_slope)
plot.bar <- 
    ggplot(df_z_slope, aes(x = method, y = slope, color = method, fill = method)) +
    geom_bar(stat = 'identity') +
    scale_color_manual(breaks = label.methods, values = vec.color, labels = label.methods) +
    scale_fill_manual(breaks = label.methods, values = vec.color, labels = label.methods) +
    labs(title = "", y = "Z-score of\n the fitted line's slope", x = '') + 
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
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 11, color = 'black', family = 'Arial'))
ggsave(filename = 'barplot_slope_bacc.png', 
       path = path.fig, plot = plot.bar,
       units = 'cm', height = 11, width = 4)
#########################################

## cluster accuracy
#########################################
df_plot_index_acc <- df_plot_index[df_plot_index$term %in% 
                                       c('Number of correctly labeled clusters', 
                                         'Total number of clusters'),]
vec.methods <- unique(df_plot_index_acc$method)
vec.datasets <- unique(df_plot_index_acc$dataset)
df_cluster_acc <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    for (dataset in vec.datasets) {
        num_correct <- df_plot_index_acc[df_plot_index_acc$method == method & df_plot_index_acc$dataset == dataset & 
                                    df_plot_index_acc$term == 'Number of correctly labeled clusters', 'value']
        num_total <- df_plot_index_acc[df_plot_index_acc$method == method & df_plot_index_acc$dataset == dataset & 
                                  df_plot_index_acc$term == 'Total number of clusters', 'value']
        type <- df_plot_index_acc[df_plot_index_acc$method == method & df_plot_index_acc$dataset == dataset & 
                             df_plot_index_acc$term == 'Total number of clusters', 'type']
        index_rank <- df_plot_index_acc[df_plot_index_acc$method == method & df_plot_index_acc$dataset == dataset & 
                                            df_plot_index_acc$term == 'Total number of clusters', 'rank_combine_index']
        df_cluster_acc <- rbind(df_cluster_acc, data.frame(method = method, dataset = dataset,
                                                           value = num_correct/num_total,
                                                           index_rank = index_rank,
                                                           type = type))
    }
}

df_mean <- data.frame(stringsAsFactors = F)
for (method in vec.methods) {
    df_sub <- df_cluster_acc[df_cluster_acc$method == method,]
    df_sub_mean <- data.frame(stringsAsFactors = F)
    df_sub_mean <- rbind(df_sub_mean, 
                         data.frame(term = 'Cluster accuracy', method = method, 
                                    value = mean(df_sub$value),
                                    dataset = 'Mean Accuracy', type = 'Average'))
    df_mean <- rbind(df_mean, df_sub_mean)
}
sort.methods <- rev(df_mean$method[order(df_mean$value, decreasing = F)])
label.methods <- sort.methods
label.methods <- gsub('scmapcell', 'scmap-cell', label.methods)
label.methods <- gsub('scmapcluster', 'scmap-cluster', label.methods)
label.methods <- gsub('seurat', 'Seurat v4', label.methods)
label.methods <- gsub('SVM', 'SVM_rejection', label.methods)
vec_color <- rev(c(brewer.pal(12, 'Paired')[1:5], brewer.pal(12, 'Paired')[7:10], 
                   brewer.pal(12, 'Paired')[12], brewer.pal(12, 'Paired')[6]))
getPalette <- colorRampPalette(vec_color)
vec.color <- getPalette(length(sort.methods))
vec.shape <- c(19, 0, 2, 3, 1, 4, 23, 6, 10, 8, 12, 17, 15, 13, 11)

plot_index <- 
    ggplot(df_cluster_acc, aes(x = index_rank, y = value, color = method, shape = method)) + 
    geom_point(size = 1.2) + 
    geom_smooth(method = lm, formula = y ~ poly(x, 1), se = F, size = 0.8) + 
    labs(x = 'Ranks of Classification Difficulty Index (CDI)', y = 'Cluster accuracy', color = 'Method', shape = 'Method') + 
    scale_color_manual(breaks = sort.methods, values = vec.color, labels = label.methods) +
    scale_shape_manual(breaks = sort.methods, values = vec.shape, labels = label.methods) +
    theme_bw() + 
    theme(panel.background = element_rect(colour = 'black', fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12, color = "black", family = 'Arial'),
          axis.text.x = element_text(size = 10, color = "black", family = 'Arial'),
          axis.text.y = element_text(
              size = 10, color = "black", family = 'Arial'),
          legend.title = element_text(
              size = 12, color = "black", family = 'Arial'),
          legend.text = element_text(size = 11, family = 'Arial')) + 
    guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

path.fig <- '/mdshare/node9/zy/MAGIC/fig/'
ggsave(filename = 'Difficulty_clusteraccuracy.png', 
       path = path.fig, plot = plot_index,
       units = 'cm', height = 11, width = 20)

vec_slope <- c()
for (method in sort.methods) {
    sub_df <- df_plot_index_acc[df_plot_index_acc$method == method, c("value", "rank_combine_index")]
    # print(method)
    fit_summary <- summary(lm(value ~ rank_combine_index, data = sub_df))
    vec_slope <- c(vec_slope, fit_summary$coefficients['rank_combine_index', 'Estimate'])
}
names(vec_slope) <- label.methods
z_slope <- as.data.frame(scale(vec_slope))
df_z_slope <- data.frame(method = rownames(z_slope), slope = z_slope$V1)
df_z_slope$method <- factor(df_z_slope$method, levels = rev(label.methods))
# sort(vec_slope)
plot.bar <- 
    ggplot(df_z_slope, aes(x = method, y = slope, color = method, fill = method)) +
    geom_bar(stat = 'identity') +
    scale_color_manual(breaks = label.methods, values = vec.color, labels = label.methods) +
    scale_fill_manual(breaks = label.methods, values = vec.color, labels = label.methods) +
    labs(title = "", y = "Z-score of\n the fitted line's slope", x = '') + 
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
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 11, color = 'black', family = 'Arial'))
ggsave(filename = 'barplot_slope_cacc.png', 
       path = path.fig, plot = plot.bar,
       units = 'cm', height = 11, width = 4)
#########################################


# boxplot CDI
df_CDI <- unique(df_plot_index_acc[, c("dataset", "type", "rank_combine_index")])
level_type <- c("Cross-validation", "Cross-sample", "Cross-platform", 
                "High-depth reference", "Low-depth reference", "Cross-species")
label_type <- c("Cross-validation, 1st situation", "Cross-sample, 1st situation", 
                "Cross-platform, 1st situation", "High-depth reference, 2nd situation", 
                "Low-depth reference, 2nd situation", "Cross-species, 3rd situation")
df_CDI$type <- factor(df_CDI$type, levels = level_type, labels = label_type)

plot.box_CDI <- 
    ggplot(df_CDI, aes(x = type, y = rank_combine_index)) +
    geom_boxplot(width = 0.7, outlier.size = 1, color = "#1F78B4") +
    labs(y = 'Ranks of Classification Difficulty Index (CDI)', x = 'Situations and scenarios of benchmark tests') + 
    coord_flip() + 
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          axis.text.y = element_text(
              size = 10, color = "black", family = 'Arial'),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(size = 9, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'))
ggsave(filename = 'boxplot_CDI.png', 
       path = path.fig, plot = plot.box_CDI,
       units = 'cm', height = 10, width = 15)


