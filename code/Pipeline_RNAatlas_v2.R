library(circRNAprofiler)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(gridExtra)
library(data.table)
library(cowplot)
library(dplyr)
library(ggsci)
library(scales)
library(scCustomize)
library(stringr)
library(tibble)
library(tidyr)

################################################################################
# 11. Deconvolution using CIBERSORTX
################################################################################

source('code/Utils.R')

# Results & Visualization
results = data.frame(fread('CA_deconv_perm100.csv'))

rownames(results) = results$Mixture
results$Mixture = NULL
results = results[,!(colnames(results) %in% c('P.value','Correlation','RMSE'))]

# Proportional Bar Plots
results$ID = rownames(results)
results_melted = reshape2::melt(results)

results_melted$ID = factor(results_melted$ID,
                           levels = rev(c('Brain cerebellum','Brain parietal cortex','Brain stem','Brain striatum',
                                          'Fetal Brain','Brain occipital cortex',
                                          'Brain frontal cortex','brain')))

results_melted$variable = factor(results_melted$variable, 
                                 levels = rev(c(
                                   'Inhibitory',
                                   'Excitatory',
                                   'Oligodendrocytes',
                                   'Astrocytes',
                                   'Microglia',
                                   'OPCs'
                                   )))

ggplot(results_melted, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'fill') +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = 'Cell Type Deconvolution', x = '', y = 'Proportion')+
  ggsci::scale_fill_npg()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

################################################################################
# 12. Cell Type Proportion <-> Total SCI
################################################################################
# Load Count
cutoff_prop_samples = 0.2
type = 'tissue'

# SCI Plots
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Except Fetal Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']; results$OPCs = NULL

# Filter out circRNAs with less than 20% of samples
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]
cts_sci_filtered = cts_sci_filtered[,c('circRNA_id',list_brain)]
cts_sci_filtered[1:4,1:4]

Scaled_SCI = scale(colSums(cts_sci_filtered[-1], na.rm = T))%>%data.frame
Sum_SCI = colSums(cts_sci_filtered[-1], na.rm = T)%>%data.frame
SCI = cbind(Scaled_SCI, Sum_SCI)
names(SCI) = c('Scaled_SCI', 'Sum_SCI')
rm(Scaled_SCI) ; rm(Sum_SCI)
p_list = lapply(colnames(results), function(celltype){
  p1 = ggplot(merged, aes(x = merged[,celltype], y = Sum_SCI, label = Row.names, fill = Sum_SCI)) +
    geom_point(alpha = 0.5, size = 3, shape= 21) +
    geom_smooth(method = 'lm', se = T, color = 'black') +
    geom_text_repel() +
    labs(title = celltype, x = celltype, y = 'Sum_SCI') +
    theme_classic() +
    scale_fill_gradient(low = 'blue', high = 'red')
  cor_df[celltype,]$coef %>% round(2) %>% as.character -> coef
  cor_df[celltype,]$p.value %>% round(2) %>% as.character -> p_value
  
  p2 = p1 + annotate('text', x = min(merged[,celltype], na.rm = T), y = max(merged$Sum_SCI, na.rm = T), 
                     label = paste('r =', coef, 'p =', p_value), hjust = 0, vjust = 1)
  
  return(p2)
})
names(p_list) = colnames(results)


################################################################################
# 13. ADAR3-binding circRNA expression vs. Non-binding circRNA expression
################################################################################
rm(list=ls())
gc()
source('scripts/manuscript/Utils.R')

type = 'tissue'
cutoff_prop_samples = 0.2

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Only Use Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
cts_sci = cts_sci[,c('circRNA_id', list_brain)]
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples


# irCLASH binding circRNA
list_ADARbound = readRDS('results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')

core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
core_enriched_adar3 = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
core_enriched_adar = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR_One_Union',]$core_enrichment, pattern = '/'))
list_ADARbound$ADAR3_One_Core_Enriched = core_enriched_adar3
list_ADARbound$ADAR_One_Union_Core_Enriched = core_enriched_adar

# Grouping
cts_sci_filtered$Group = ifelse(cts_sci_filtered$circRNA_id %in% list_ADARbound$ADAR3_One, 'Yes', 'No')
cts_sci_filtered[-1] %>% na.omit %>% group_by(Group) %>% summarise_all(mean, na.rm = T) -> x
log2(x[2,-1] + 1) - log2(x[1,-1] + 1)
cts_sci_filtered[-1] %>% na.omit %>% group_by(Group) %>% summarise_all(median, na.rm = T)


# Box Plot
exp_sum = (rowSums(na.omit(cts_sci_filtered[-c(1,ncol(cts_sci_filtered))])))
exp_bin = ntile(exp_sum, 10)
exp_bin%>%unique%>%length

df_to_box = reshape2::melt(na.omit(cts_sci_filtered)[exp_bin %in% c(1,2,3,4,5,6,7,8,9,10),])
df_to_box$Group = factor(df_to_box$Group, levels = c('Yes','No'))


df_to_box %>%
  group_by(variable) %>%
  mutate(rank = rank(-value)) %>%
  mutate(delabel = ifelse(rank <= n() * 0.01, circRNA_id, NA)) %>%
  ungroup() %>%
  dplyr::select(-rank) -> df_to_box_filtered

lapply(list_brain, function(x){
  p1 = ggplot(df_to_box[df_to_box$variable == x,], aes(x = Group, y = value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(position=position_jitterdodge()) +
    # geom_violin(scale="width")+
    # geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 0.1) +
    labs(
      # title = 'circRNA-Host Correlation', 
      x = 'ADARB2-binding',
      y = 'SCI') +
    theme(legend.position = 'bottom') +
    # stat_compare_means() +
    scale_fill_npg() +
    theme_bw() + # white background
    theme(legend.position = 'right',
          # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
          # axis.line.x = element_line(color="black", size = 0.5),
          # axis.line.y = element_line(color="black", size = 0.5)
    ) +
    ggtitle(x) +
    geom_signif(comparisons = list(c('Yes','No')), 
                map_signif_level = TRUE, 
                test = "wilcox.test",
                y_position = 20000,
                tip_length = 0.01
    ) +
    coord_cartesian(ylim = c(0,1e5))
  
  return(p1)
}) -> list_p

p_combined = ggarrange(plotlist = list_p, ncol = 3, nrow =3, common.legend = T, legend = 'right')

# Gene-level
cts_agg_filtered = readRDS(file = 'results/cts_agg_filtered.rds')
cts_agg_filtered = cts_agg_filtered[,c('hgnc_symbol',list_brain)]

# RIP-seq
df_gene_bind = read.csv(file = 'results/ADAR3_Binding_RIP.txt', sep = '\t')
gene_bind = df_gene_bind$Gene

# Universe Genesets
results_to_save = readRDS(file.path('results/', 'DEG_SCI_Aggregated_Sum.rds'))

cutoff_padj = 0.05
cutoff_LFC = 1
up_host_gene = results_to_save[results_to_save$fold_change > cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
universe_host_gene = na.omit(results_to_save$hgnc_symbol) ; rm(results_to_save)

# Only Genes included in Host Gene Universe
gene_bind_subset = gene_bind[gene_bind %in% universe_host_gene]

df_gene_bind_subset = df_gene_bind[df_gene_bind$Gene %in% gene_bind_subset,]
df_gene_bind_subset = df_gene_bind

apply(df_gene_bind_subset[-c(1,2)], 2, function(x){
  gs_symbol = df_gene_bind_subset$Gene[unlist(x) == 'yes' | unlist(x) == 'Yes']
  return(gs_symbol)
}) -> df_geneset

# Process Genesets and generate curated geneset
names(df_geneset) = str_replace_all(names(df_geneset), pattern = '\\.', replacement = '')
df_geneset$RIPseq = df_gene_bind_subset$Gene
df_geneset$RIPseq_inhibitory = intersect(df_geneset$BoundbyADAR2inHEKcellsinSongetal2020, df_geneset$BoundbyADAR1inBahnetal2015)
df_geneset$RIPseq_inhibitory_sum = union(df_geneset$BoundbyADAR1inBahnetal2015, df_geneset$BoundbyADAR2inHEKcellsinSongetal2020)
df_geneset$RIPseq_inhibitory2 = intersect(df_geneset$RIPseq_inhibitory, df_geneset$BoundbyADAR3inHEKcellsinSongetal2020)

# ADARB2-binding vs. Non-binding
cts_agg_filtered$Group = ifelse(cts_agg_filtered$hgnc_symbol %in% intersect(up_host_gene,df_geneset$RIPseq_inhibitory), 'Yes', 'No')
cts_agg_filtered[-1] %>% na.omit %>% group_by(Group) %>% summarise_all(median, na.rm = T)

exp_sum = (rowSums(na.omit(cts_agg_filtered[-c(1,ncol(cts_agg_filtered))])))
exp_bin = ntile(exp_sum, 10)

df_to_box = reshape2::melt(na.omit(cts_agg_filtered)[exp_bin %in% c(1,2,3,4,5,6,7,8,9,10),])
df_to_box$Group = factor(df_to_box$Group, levels = c('Yes','No'))


df_to_box %>%
  group_by(variable) %>%
  mutate(rank = rank(-value)) %>%
  mutate(delabel = ifelse(rank <= n() * 0.01, hgnc_symbol, NA)) %>%
  ungroup() %>%
  dplyr::select(-rank) -> df_to_box_filtered

df_to_box_filtered$delabel%>%
  table %>%
  data.frame %>%
  arrange(desc(Freq))

lapply(list_brain, function(x){
  p1 = ggplot(df_to_box[df_to_box$variable == x,], aes(x = Group, y = value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(position=position_jitterdodge()) +
    # geom_violin(scale="width")+
    # geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 0.1) +
    labs(
      # title = 'circRNA-Host Correlation', 
      x = 'ADARB2-binding',
      y = 'SCI') +
    theme(legend.position = 'bottom') +
    # stat_compare_means() +
    scale_fill_npg() +
    theme_bw() + # white background
    theme(legend.position = 'right',
          # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
          # axis.line.x = element_line(color="black", size = 0.5),
          # axis.line.y = element_line(color="black", size = 0.5)
    ) +
    ggtitle(x) +
    geom_signif(comparisons = list(c('Yes','No')), 
                map_signif_level = TRUE, 
                test = "wilcox.test",
                y_position = 20000,
                tip_length = 0.01
    ) +
    coord_cartesian(ylim = c(0,1e5))
  return(p1)
}) -> list_p

p_combined = ggarrange(plotlist = list_p, ncol = 3, nrow =3, common.legend = T, legend = 'right')
