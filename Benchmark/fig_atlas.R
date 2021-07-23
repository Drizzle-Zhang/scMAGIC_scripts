library(ggplot2)

path.res <- '/mdshare/node9/zy/scRef/figure/atlas_anno/'
vec.ref <- c('MCA', 'MCA', 'HCL', 'HCL')
vec.dataset <- c('Tasic2018', 'Haber_Duodenum', 'panc8_indrop', 'hPBMC_10Xv2')
vec.title <- c('MCA -> Mouse Neocortex', 'MCA -> Mouse Duodenum', 
               'HCL -> Human Pancreas', 'HCL -> Human PBMC')

df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(vec.dataset)) {
    reference <- vec.ref[i]
    dataset <- vec.dataset[i]
    title <- vec.title[i]
    file.res.scMAGIC <- paste0(path.res, 'RES_', reference,'_', dataset, '_scMAGIC.Rdata')
    res.scMAGIC <- readRDS(file.res.scMAGIC)
    df.sub <- data.frame(term = 'Accuracy', method = 'scMAGIC',
                         value = res.scMAGIC$accuracy, stringsAsFactors = F)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'Balanced accuracy', method = 'scMAGIC',
                               value = res.scMAGIC$balanced.accuracy, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'labeled Accuracy', method = 'scMAGIC',
                               value = res.scMAGIC$accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'labeled Balanced accuracy', method = 'scMAGIC',
                               value = res.scMAGIC$balanced.accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub$dataset <- rep(title, nrow(df.sub))
    df.plot <- rbind(df.plot, df.sub)
}
df.plot$dataset <- factor(df.plot$dataset, levels = rev(vec.title))
df.plot$term <- factor(df.plot$term, 
                       levels = c('labeled Balanced accuracy', 'labeled Accuracy',
                                  'Balanced accuracy', 'Accuracy'))
# library(RColorBrewer)
# brewer.pal(11,"Paired")
plot.bar <-
    ggplot(df.plot, aes(x = dataset, y = value, color = term, fill = term)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    # geom_hline(yintercept  = 0.9, color = 'grey', linetype = 2) +
    scale_color_manual(breaks = c('Accuracy', 'Balanced accuracy', 'labeled Accuracy',
                                  'labeled Balanced accuracy'),
                       values = c("#80B1D3", "#5AB4AC", "#BC80BD", "#FFA07A")) +
    scale_fill_manual(breaks = c('Accuracy', 'Balanced accuracy', 'labeled Accuracy',
                                 'labeled Balanced accuracy'),
                      values = c("#80B1D3", "#5AB4AC", "#BC80BD", "#FFA07A")) +
    labs(title = "", y = '', x = '', color = '', fill = '') +
    coord_flip() +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 10, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 10, family = 'Arial'),
          legend.position = 'right',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(filename = 'Balanced_accuracy_atlas.png',
       path = path.res, plot = plot.bar,
       units = 'cm', height = 8, width = 17)

# delete balanced accuracy
df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(vec.dataset)) {
    reference <- vec.ref[i]
    dataset <- vec.dataset[i]
    title <- vec.title[i]
    file.res.scMAGIC <- paste0(path.res, 'RES_', reference,'_', dataset, '_scMAGIC.Rdata')
    res.scMAGIC <- readRDS(file.res.scMAGIC)
    df.sub <- data.frame(term = 'Accuracy', method = 'scMAGIC',
                         value = res.scMAGIC$accuracy, stringsAsFactors = F)
    df.sub <- rbind(df.sub, 
                    data.frame(term = 'labeled Accuracy', method = 'scMAGIC',
                               value = res.scMAGIC$accuracy.rm.unassigned, stringsAsFactors = F))
    df.sub$dataset <- rep(title, nrow(df.sub))
    df.plot <- rbind(df.plot, df.sub)
}
df.plot$dataset <- factor(df.plot$dataset, levels = rev(vec.title))
# library(RColorBrewer)
# brewer.pal(11,"Paired")
plot.bar <-
    ggplot(df.plot, aes(x = dataset, y = value, color = term, fill = term)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    scale_color_manual(breaks = c('Accuracy', 'labeled Accuracy'),
                       values = c("#80B1D3", "#FFA07A")) +
    scale_fill_manual(breaks = c('Accuracy', 'labeled Accuracy'),
                      values = c("#80B1D3", "#FFA07A")) +
    labs(title = "", y = '', x = '', color = '', fill = '') +
    coord_flip() +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 10, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 10, family = 'Arial'),
          legend.text = element_text(size = 10, family = 'Arial'),
          legend.position = 'right',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(filename = 'accuracy_atlas.png',
       path = path.res, plot = plot.bar,
       units = 'cm', height = 6, width = 16)



