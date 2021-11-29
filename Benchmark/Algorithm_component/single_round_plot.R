library(ggplot2)
library(RColorBrewer)

file_1situ <- '/mdshare/node9/zy/MAGIC/evaluation/First_situation.txt'
df_1situ <- read.delim(file_1situ, sep = '\t', stringsAsFactors = F)
# 
file_2situ <- '/mdshare/node9/zy/MAGIC/evaluation/Second_situation.txt'
df_2situ <- read.delim(file_2situ, sep = '\t', stringsAsFactors = F)

file_3situ <- '/mdshare/node9/zy/MAGIC/evaluation/third_situation.txt'
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
df_plot_index_acc <- df_plot_index[df_plot_index$term == 'Accuracy',]
df_plot_index_acc_scMAGIC <- 
    df_plot_index_acc[df_plot_index_acc$method == 'scMAGIC', 
                      c("dataset", "value", "type", "combine_index", "rank_combine_index")]
names(df_plot_index_acc_scMAGIC) <- c("dataset", "accuracy", "type", "combine_index", "rank_combine_index")

file_single <- '/mdshare/node9/zy/MAGIC/evaluation/scMGAIC_Single.txt'
df_single <- read.delim(file_single, sep = '\t', stringsAsFactors = F)
df_single_acc <- df_single[df_single$term == 'Accuracy', c('dataset', 'value')]
names(df_single_acc) <- c('dataset', 'accuracy_single')
df_single_prop <- df_single[df_single$term == 'Proportion of confidently validated cells', c('dataset', 'value')]
names(df_single_prop) <- c('dataset', 'prop')

df_single_plot <- merge(df_single_acc, df_single_prop, by = 'dataset')
df_single_plot <- merge(df_single_plot, df_plot_index_acc_scMAGIC, by = 'dataset')

plot_single_acc <- 
    ggplot(df_single_plot, aes(x = rank_combine_index, y = accuracy_single)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 1), se = F) + 
    labs(x = 'Ranks of Classification Difficulty Index (CDI)', 
         y = 'Accuracy of the\nconfidently validated query cells', 
         color = 'Method', shape = 'Method') + 
    theme_bw() + 
    theme(panel.background = element_rect(colour = 'black', fill = 'transparent'),
          panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          axis.text = element_text(size = 7, color = "black", family = 'Arial'),
          axis.title = element_text(size = 10, color = "black", family = 'Arial'))
path.fig <- '/mdshare/node9/zy/MAGIC/fig/'
ggsave(filename = 'single_acc.png', 
       path = path.fig, plot = plot_single_acc,
       units = 'cm', height = 6.3, width = 10)

######### lm
summary(lm(rank_combine_index ~ accuracy_single, df_single_plot))

########

plot_single_prop <- 
    ggplot(df_single_plot, aes(x = rank_combine_index, y = prop)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 1), se = F) + 
    labs(x = 'Ranks of Classification Difficulty Index (CDI)', 
         y = 'Percentage of the\nconfidently validated query cells', 
         color = 'Method', shape = 'Method') + 
    theme_bw() + 
    theme(panel.background = element_rect(colour = 'black', fill = 'transparent'),
          panel.border = element_rect(size = 1.2),
          panel.grid = element_blank(),
          axis.text = element_text(size = 7, color = "black", family = 'Arial'),
          axis.title = element_text(size = 10, color = "black", family = 'Arial'))
path.fig <- '/mdshare/node9/zy/MAGIC/fig/'
ggsave(filename = 'single_prop.png', 
       path = path.fig, plot = plot_single_prop,
       units = 'cm', height = 6.3, width = 10)

######### lm
summary(lm(rank_combine_index ~ prop, df_single_plot))

########


plot_single_twoacc <- 
    ggplot(df_single_plot, aes(x = accuracy_single, y = accuracy)) + 
    geom_point() + 
    geom_smooth(method = lm, formula = y ~ poly(x, 1), se = F) + 
    labs(x = 'Accuracy of the confidently validated query cells', 
         y = 'Accuracy of scMAGIC', color = 'Method', shape = 'Method') + 
    theme_bw() + 
    theme(panel.background = element_rect(colour = 'black', fill = 'transparent'),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 1.2),
          axis.text = element_text(size = 7, color = "black", family = 'Arial'),
          axis.title = element_text(size = 10, color = "black", family = 'Arial'))
path.fig <- '/mdshare/node9/zy/MAGIC/fig/'
ggsave(filename = 'single_twoacc.png', 
       path = path.fig, plot = plot_single_twoacc,
       units = 'cm', height = 6.3, width = 9.5)

######### lm
summary(lm(accuracy_single ~ accuracy, df_single_plot))

########


