getwd()
setwd('/home/local/suwanyu_971004/circRNA/')

rm(list=ls())
gc()

library(circRNAprofiler)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(gridExtra)
library(data.table)
library(cowplot)
library(dplyr)
library(stringr)


################################################################################
# 1. Calculate Linear-to-Circular Ratio
################################################################################
dir.create('results/RESTART/1.calculateSCI')

# Load circRNA only Annotated
annot = data.frame(fread(file = 'results/circRNAprofiler/annotatedBSJs_withAlu_IntronPair.tsv',sep = '\t'))[c('circRNA_id','gene')]
cts_rpm = read.csv(file = 'circRNA_GSE138734/GSE138734_RNAAtlas_circRNA_backspliced_rpm_filtered.txt', sep = '\t')

cts_rpm = cts_rpm[cts_rpm$circRNA_id %in% annot$circRNA_id,]

# Load metadata
meta = read.csv(file = 'circRNA_GSE138734/GSE138734_metadata_v2.txt', sep = '\t')
names(meta)[names(meta) == 'RNA_Atlas_Number'] = 'sample'

# Filter Meta data
samples = names(cts_rpm)[!(names(cts_rpm) %in% c('circRNA_id'))]
meta = meta[meta$sample %in% samples,]
meta = meta %>% arrange(sample)

# Merge circRNA cts with Gene Symbol
cts_rpm = left_join(cts_rpm, annot, by = 'circRNA_id')

# Load linear RNA cts
rna_tpm = data.frame(fread(file = 'circRNA_GSE138734/GSE138734_RNAAtlas_totalRNA_tpm_filtered.txt', sep= '\t'))
rna_tpm$genome_coordinates = NULL
names(rna_tpm)[1] = 'gene'

# Sanity Check: Check whether the order of the columns are the same
if (!all(names(rna_tpm[-1]) == names(cts_rpm[!(names(cts_rpm) %in% c('circRNA_id','gene'))]))){
  stop('The order of the columns are not the same')
}

# Subset circRNA cts and linear RNA cts
list_genes_its = intersect(cts_rpm$gene, rna_tpm$gene)
cts_rpm = cts_rpm[cts_rpm$gene %in% list_genes_its,]
rna_tpm = rna_tpm[rna_tpm$gene %in% list_genes_its,]


# Check Average TPM values by Tissue, Cell Line, Cell Types
lapply(c('cancer cell line','cell line','cell type','tissue','all'), function(type){
  samples_to_subset = meta$sample[meta$type == type]
  if (type == 'all'){
    samples_to_subset = meta$sample
  }
  return(rowMeans(rna_tpm[samples_to_subset]))
}) -> meta_tpm
meta_tpm = data.frame(meta_tpm)
names(meta_tpm) = c('cancer cell line','cell line','cell type','tissue','all')
meta_tpm$gene = rna_tpm$gene
write.csv(meta_tpm, file = 'results/RESTART/1.calculateSCI/AverageTPM_linearRNA.csv', row.names = F)

# Check Average RPM values by Tissue, Cell Line, Cell Types
lapply(c('cancer cell line','cell line','cell type','tissue','all'), function(type){
  samples_to_subset = meta$sample[meta$type == type]
  if (type == 'all'){
    samples_to_subset = meta$sample
  }
  return(rowMeans(cts_rpm[samples_to_subset]))
}) -> meta_rpm
meta_rpm = data.frame(meta_rpm)
names(meta_rpm) = c('cancer cell line','cell line','cell type','tissue','all')
meta_rpm$circRNA = cts_rpm$circRNA
meta_rpm$gene = cts_rpm$gene
write.csv(meta_rpm, file = 'results/RESTART/1.calculateSCI/AverageRPM_circRNA.csv', row.names = F)


################################################################################
# Gene Filtering or Generate NA values
################################################################################
cutoff_tpm_elements = 1

# NA value if linear RNA expression is under 1 TPM
rna_tpm[-1][rna_tpm[-1] <= cutoff_tpm_elements] = NA

# # Gene Filtering
# gghistogram(meta_tpm[meta_tpm$all <= 100,], 
#             x = 'all', 
#             bins = 100, 
#             title = 'Histogram of Average TPM values of Linear RNA',
#             )
# 
# # Gene Filtering
# gghistogram(meta_tpm[meta_tpm$tissue <= 100,],
#             x = 'tissue', 
#             bins = 100, 
#             title = 'Histogram of Average TPM values of Linear RNA',
# )

cutoff_tpm = 0
target_type = 'tissue'

list_genes_filtered = meta_tpm[meta_tpm[[target_type]] > cutoff_tpm,]$gene

# Subset circRNA cts and linear RNA cts
cts_rpm = cts_rpm[cts_rpm$gene %in% list_genes_filtered,]
rna_tpm = rna_tpm[rna_tpm$gene %in% list_genes_filtered,]

# Merge circRNA cts and linear RNA cts
names(rna_tpm)[-1] = paste0(names(rna_tpm)[-1],'_total')
merged = left_join(cts_rpm, rna_tpm, by = 'gene')

# Sanity Check
if (!all(paste0(names((merged[2:297])), '_total') == names(merged[299:ncol(merged)]))){
  stop('The order of the columns are not the same')
}

# SCI with pseudocount
pseudocount = 1e-12
LCratio = (merged[2:297] + pseudocount)/ (merged[299:ncol(merged)] + pseudocount)
LCratio$circRNA_id = merged$circRNA_id

# Check whether Any Value is above 1
table(LCratio == 1)
table(LCratio > 1)
is.na(LCratio) %>% table

# Save Results
fwrite(LCratio, file = 'results/RESTART/1.calculateSCI/LCratio_withNAandPseudocount.txt', sep = '\t')

################################################################################
# 1. Calculate Scaled-Circularity-Index (Including NA)
################################################################################
LC_ratio = fread('results/RESTART/1.calculateSCI/LCratio_withNAandPseudocount.txt')
LC_ratio = data.frame(LC_ratio)
LC_ratio[-ncol(LC_ratio)] = LC_ratio[-ncol(LC_ratio)] * 1e6
LC_ratio = LC_ratio[c(ncol(LC_ratio), 1:(ncol(LC_ratio)-1))]

if (!(names(LC_ratio)[1] == 'circRNA_id')){
  stop('The first column should be circRNA_id')
}

fwrite(LC_ratio, 'results/RESTART/1.calculateSCI/LCratio_withNAandPseudocount_Milion.txt', sep = '\t')


################################################################################
# 2. EDA analysis (background_2.EDA.R)
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/2.EDA', showWarnings = F)
source('scripts/manuscript/Utils.R')

require(factoextra)
require(missMDA)
require(FactoMineR)

cutoff_prop_samples = 0.2
type = 'tissues'
title = 'circRNA SCI values, imputed'

# SCI Plots
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Filter out circRNAs with less than 20% of samples
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]


log_rna_matrix = log2(cts_sci_filtered[-1] + 1)
# log_rna_matrix = log_rna_matrix[,colnames(log_rna_matrix) != 'testicle']
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]

nFeatures = nrow(log_rna_matrix)

# Transpose matrix to have samples as rows and features as columns
log_rna_matrix_t <- data.frame(t(log_rna_matrix))

# Imputation
imputed_data <- missMDA::imputePCA(log_rna_matrix_t, ncp = 2)

# Perform PCA on the imputed data
pca_res = FactoMineR::PCA(imputed_data$completeObs, graph = FALSE, ncp = 2, scale.unit = TRUE)

dir.create('results/RESTART/6.FeaturePlot', showWarnings = F)
saveRDS(pca_res, file = str_glue('results/RESTART/6.FeaturePlot/pca_res_{type}.rds'))

# Prepare data for ggplot
pca_df = as.data.frame(pca_res$ind$coord)
pca_df$Sample = rownames(pca_df)
pca_df = left_join(pca_df, meta, by = c("Sample" = "name"))
saveRDS(pca_df, file = str_glue('results/RESTART/6.FeaturePlot/pca_df_{type}.rds'))

# Extract variance explained for plotting
var_explained <- pca_res$eig
percent_variance_explained <- var_explained[, 2]  # Percentage of variance

# Plotting the first two principal components
ggplot(pca_df, aes(x = Dim.1, y = Dim.2, label = Sample, color = organ_system)) +
  geom_point() +
  xlab(paste("Principal Component 1 (", sprintf("%.2f", percent_variance_explained[1]), "%)", sep="")) +
  ylab(paste("Principal Component 2 (", sprintf("%.2f", percent_variance_explained[2]), "%)", sep="")) +
  ggtitle(paste0(title, ',', nFeatures)) +
  theme_minimal() +
  geom_text_repel() +
  theme(legend.position = "none")


# Feature Plot (SCI)
meta_to_add = data.frame(Sample = colnames(cts_sci_filtered[-1]), Average_SCI = colMeans(cts_sci_filtered[-1], na.rm = T), Sum_SCI = colSums(cts_sci_filtered[-1], na.rm = T))
pca_df_merged = left_join(pca_df, meta_to_add, by = c('Sample' = 'Sample'))
ggplot(pca_df_merged, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(color='black', aes(fill = Average_SCI), size = 3, shape = 21) +
  xlab(paste("Principal Component 1", sep="")) +
  ylab(paste("Principal Component 2", sep="")) +
  ggtitle("Average SCI by tissue") +
  theme_minimal() +
  geom_text_repel() +
  scale_fill_gradient(low = "blue", high = "red")

ggplot(pca_df_merged, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(color='black', aes(fill = Sum_SCI), size = 3, shape = 21) +
  xlab(paste("Principal Component 1", sep="")) +
  ylab(paste("Principal Component 2", sep="")) +
  ggtitle("Sum SCI by tissue") +
  theme_minimal() +
  geom_text_repel() +
  scale_fill_gradient(low = "blue", high = "red")


meta_to_add = data.frame(Sample = colnames(cts_sci_filtered[-1]), Scaled_SCI = scale(colSums(cts_sci_filtered[-1], na.rm = T)))
pca_df_merged = left_join(pca_df, meta_to_add, by = c('Sample' = 'Sample'))

ggplot(pca_df_merged, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(color='black', aes(fill = Scaled_SCI), size = 3, shape = 21) +
  xlab(paste("Principal Component 1", sep="")) +
  ylab(paste("Principal Component 2", sep="")) +
  ggtitle("Scaled SCI by tissue ()") +
  theme_minimal() +
  geom_text_repel() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

################################################################################
# 2. Tissue Specificity (Silhouette Score by Tissues)
################################################################################
rm(list=ls())
gc()

library(cluster)
source('scripts/manuscript/Utils.R')
Boxplot_Silhouette = function(log_rna_matrix, meta, title){
  library(scales)
  library(ggsignif)
  
  # Calculate silhouette  
  rownames(meta) = meta$name
  dist_for_sil = dist(data.frame(t(na.omit(log_rna_matrix)))) 
  
  chr_name = meta[colnames(log_rna_matrix),]$name
  factor_organ = factor(meta[colnames(log_rna_matrix),]$organ_system)
  df_organ = data.frame(organ = factor_organ, cluster = as.numeric(factor_organ), name = chr_name)
  
  if (!all(df_organ$name == colnames(log_rna_matrix))){
    stop()
  }
  
  res_sil = data.frame(silhouette(df_organ$cluster, dist_for_sil))
  df_organ$cluster = NULL
  res_sil = cbind(res_sil, df_organ)
  res_sil$sil_width_scaled = scale(res_sil$sil_width)
  
  # Group samples  
  res_sil$Group = ifelse(grepl(res_sil$name, pattern = 'Brain', ignore.case = T), 'Brain', 'Others')
  res_sil$Group = ifelse(grepl('Heart|ventricle|atruim|atrium', res_sil$name, ignore.case = T), 'Heart',
                         ifelse(res_sil$Group == 'Brain', 'Brain','Others')) %>%
    factor(levels = c('Brain','Heart', 'Others'))
  
  # Plots
  hex <- hue_pal()(15)
  color_map = c(hex[2:3], 'grey')
  names(color_map) = c('Heart','Brain','Others')

  p = ggplot(res_sil, aes(y = sil_width, x = Group, fill = Group)) +
    geom_boxplot() +
    labs(
      title = title, 
      x = 'Sample',
      y = 'Silhouette Coefficient'
      ) +
    theme(legend.position = 'bottom') +
    # stat_compare_means() +
    scale_fill_manual(values = color_map) +
    theme_bw() + # white background
    theme(legend.position = 'right',
          # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12),
          # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
          # axis.line.x = element_line(color="black", size = 0.5),
          # axis.line.y = element_line(color="black", size = 0.5)
    )+ geom_signif(comparisons = list(c("Brain", "Heart"), c("Heart", "Others"), c("Brain", "Others")),
                   map_signif_level = TRUE, test = "wilcox.test", paired = FALSE)
  return(p)
}

cutoff_prop_samples=1
type='tissue'

# SCI 
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Filter out circRNAs with less than 20% of samples
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]
log_rna_matrix = log2(cts_sci_filtered[-1] + 1)
log_rna_matrix = log_rna_matrix[,colnames(log_rna_matrix) != 'testicle']
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]
p_sci = Boxplot_Silhouette(log_rna_matrix, meta, 'SCI by tissue')

# Linear TPM
annot = read.csv('results/circRNAprofiler/annotatedBSJs.csv', row.names = 1) %>% na.omit()
annot$circRNA_id = lapply(str_split(annot$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
# circRNA_id_Tissue = na.omit(get_data(slot = 'rpm', annot = NULL, type = 'tissue', LCratio = T)[[1]])$circRNA_id
circRNA_id_Tissue = cts_sci_filtered$circRNA_id

annot = annot[annot$circRNA_id %in% circRNA_id_Tissue,]

input_list = get_data_linear(slot = 'data', annot = annot, type = type)
cts = input_list[[1]] ;  meta = input_list[[2]]; rm()

log_rna_matrix = log2(cts[-1] + 1)
log_rna_matrix = na.omit(log_rna_matrix)
log_rna_matrix = log_rna_matrix[,colnames(log_rna_matrix) != 'testicle']
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]
p_linear = Boxplot_Silhouette(log_rna_matrix, meta, 'Linear TPM by tissue')

# Circular RPM
list_input = get_data(slot = 'rpm', type = type, LCratio = F, annot = annot)
cts = list_input[[1]]; meta = list_input[[2]]; rownames(meta) = meta$name; rm(list_input)

log_rna_matrix = log2(cts[-1] + 1)
log_rna_matrix = na.omit(log_rna_matrix)
log_rna_matrix = log_rna_matrix[,colnames(log_rna_matrix) != 'testicle']
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]
log_rna_matrix = na.omit(log_rna_matrix)

p_circle = Boxplot_Silhouette(log_rna_matrix, meta, 'Circular RPM by tissue')

combined = ggarrange(plotlist = list(p_sci, p_linear, p_circle), ncol = 3, nrow = 1, common.legend = T)
dir.create(str_glue('results/RESTART/2.EDA/Silhoette_{type}'), showWarnings = F)
ggsave(file.path(str_glue('results/RESTART/2.EDA/Silhoette_{type}'), str_glue('boxplot_{type}.pdf')), width = 6, height = 6)


################################################################################
# 3. DEG analysis from SCI values (Gene-level Aggregation)
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/3.DEG_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

# Load SCI
input_list = get_data(slot = 'rpm', annot = NULL, type = 'tissue', LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]


# Brain Tissues vs. Other Tissues
meta$group = ifelse(grepl('Brain', meta$name, ignore.case = T), 'Brain', 'Others') %>%
  factor(levels = c('Others', 'Brain'))
meta$group%>%table

# Brain Cortex vs. Other Brain vs. Other Tissues
meta$group_detail = ifelse(grepl('Brain', meta$name, ignore.case = T), 
                           ifelse(meta$name %in% c('Fetal Brain','Brain stem','Brain cerebellum'), 'Fetal Brain, Brain stem, Brain cerebellum',
                                  'Brain and Cortex and Striatum'),
                           'Others')

# Preprocess the Data
# cts = na.omit(cts_sci)
cts = cts_sci


# Aggregate circRNA by gene-level
annotatedBSJs = data.frame(fread(file = 'results/circRNAprofiler/annotatedBSJs.csv'))[-1]
annotatedBSJs$circRNA_id = lapply(str_split(annotatedBSJs$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
names(annotatedBSJs)[names(annotatedBSJs) == 'gene'] = 'hgnc_symbol'
cts = left_join(cts, annotatedBSJs[c('circRNA_id','hgnc_symbol')], by = 'circRNA_id')
cts$hgnc_symbol%>%unique%>%length

# Aggregate values row-wise by hgnc_symbol
custom_sum <- function(x) {
  if (all(is.na(x))) {
    NA
  } else {
    sum(x, na.rm = TRUE)
  }
}

# cts_agg = cts %>%
#   group_by(hgnc_symbol) %>%
#   summarise(across(where(is.numeric), custom_sum), .groups = "drop")

cts_agg = cts %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = T), .groups = "drop")



# Filter out circRNAs with less than 20% of samples
# apply(cts_sci[-1], 1, function(x){sum(is.na(x))})%>%table
cutoff_prop_samples = 0.2
cts_agg_filtered= cts_agg[apply(cts_agg[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]

cts_agg_filtered%>%names

list_symbol_id = cts_agg_filtered$hgnc_symbol
cts_agg_filtered$hgnc_symbol = NULL ; cts_agg_filtered[1:4,1:4]
saveRDS(cts_agg_filtered, file = 'results/RESTART/3.DEG_Brain/cts_agg_filtered.rds')


cts_agg_filtered = data.frame(t(cts_agg_filtered))
cts_agg_filtered$name = rownames(cts_agg_filtered)
cts_agg_filtered = left_join(cts_agg_filtered, meta[c('name','group')], by = 'name')
cts_agg_filtered%>%rownames()
colnames(cts_agg_filtered)[length(colnames(cts_agg_filtered))]

matching_tb = data.frame(name = names(cts_agg_filtered)[!(names(cts_agg_filtered) %in% c('name','group'))], hgnc_symbol = list_symbol_id)
rownames(matching_tb) = matching_tb$name

# Wilcox.test
Wilcox_SCI = function(cts){
  cts_to_test = cts[!(colnames(cts) %in% c('name','group'))]
  cts_to_test = cts_to_test[!(colSums(cts_to_test, na.rm = T) == 0)]
  
  results <- apply(cts_to_test, 2, function(x) {
    case_values = x[cts$group != "Others"]
    control_values = x[cts$group == "Others"]
    
    if (length(na.omit(case_values)) <= 3 | length(na.omit(control_values)) <= 3){
      return(c(NA, NA, NA, NA, NA))
    }
    
    test_result <- wilcox.test(x ~ cts$group)
    p_value <- test_result$p.value
    
    mean_case <- log2(mean(na.omit(case_values)) + 1)
    mean_control <- log2(mean(na.omit(control_values)) + 1)
    
    fold_change <- mean_case - mean_control # Log2 fold change
    
    median_case <- log2(median(na.omit(case_values)) + 1)
    median_control <- log2(median(na.omit(control_values)) + 1)
    median_fold_change <- median_case - median_control
    
    
    baseMean = mean(x, na.rm = T)
    
    is_na = any(is.na(x))
    
    return(c(p_value, fold_change, baseMean, is_na, median_fold_change))
  })
  
  results_df <- as.data.frame(t(results))
  colnames(results_df) <- c("p_value", "fold_change", "baseMean", "is_na","median_fold_change")
  
  results_df$p_value <- as.numeric(as.character(results_df$p_value))
  results_df$fold_change <- as.numeric(as.character(results_df$fold_change))
  results_df$baseMean <- as.numeric(as.character(results_df$baseMean))
  results_df$is_na <- as.logical(results_df$is_na)
  results_df$median_fold_change <- as.numeric(as.character(results_df$median_fold_change))
  
  return(results_df)
}

results_df = Wilcox_SCI(cts_agg_filtered)
results_df$name = rownames(results_df)
results_df = left_join(results_df, matching_tb, by = 'name')

# Omit NA
results_df = results_df[!is.na(results_df$p_value),]
results_df$p_adj = p.adjust(results_df$p_value,method = 'BH')
rownames(results_df) = results_df$hgnc_symbol

################################################################################
# 3. Functional Annotation
################################################################################
cutoff_padj = 0.05
cutoff_LFC = 1

Nuclear_RBP = data.frame(fread(file = 'RBP/NuclearRBP.txt', sep = '\t', header = F))
names(Nuclear_RBP) = 'hgnc_symbol'
Nuclear_RBP$Nuclear_RBP = T

results_df$hgnc_symbol = rownames(results_df)

results_to_save = left_join(data.frame(results_df), Nuclear_RBP, by = 'hgnc_symbol')
results_to_save$Nuclear_RBP[is.na(results_to_save$Nuclear_RBP)] = F
results_to_save$Significant = ifelse(results_to_save$p_adj < cutoff_padj & abs(results_to_save$fold_change) > cutoff_LFC, T, F)
results_to_save$Regulator = ifelse(results_to_save$p_adj < cutoff_padj,
                                   ifelse(results_to_save$fold_change >= cutoff_LFC, 'UP', 
                                          ifelse(results_to_save$fold_change <= -cutoff_LFC, 'DOWN', 'None')),
                                   'None')

results_to_save = results_to_save %>% 
  dplyr::select(hgnc_symbol, baseMean, fold_change, median_fold_change, p_value, p_adj, Significant, Regulator,
                is_na, Nuclear_RBP) %>%
  arrange(p_adj)

results_to_save$p_adj%>%is.na%>%table
results_to_save$Regulator%>%table
results_to_save$fold_change%>%is.na%>%table
# fwrite(results_to_save, file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.csv'), row.names = F)
# saveRDS(results_to_save, file = file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))

fwrite(results_to_save, file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Average.csv'), row.names = F)
saveRDS(results_to_save, file = file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Average.rds'))

# df_biomaRt = data.frame(fread(file = 'results/RESTART/3.DEG_Brain/biomaRt_240508.csv'))
# df_biomaRt = df_biomaRt[c('hgnc_symbol','entrezgene_id','description')]
# df_biomaRt = df_biomaRt[!duplicated(df_biomaRt),]
# 
# results_to_excel = left_join(results_to_save, df_biomaRt, by = 'hgnc_symbol')
# # openxlsx::write.xlsx(results_to_save, file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.xlsx'), rowNames = F)
# openxlsx::write.xlsx(results_to_save, file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Average.xlsx'), rowNames = F)

################################################################################
# 3. Visualization (Volcano Plot)
################################################################################
cutoff_padj = 0.05
cutoff_LFC = 1

# devtools::install_github("BioSenior/ggVolcano")
require(ggVolcano)
require(ggVennDiagram)
require(ggplot2)

# results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Average.rds'))
results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))
df_to_plot = results_to_save[!is.na(results_to_save$p_adj),]
df_to_plot <- add_regulate(df_to_plot,
                           log2FC_name = "fold_change",
                           fdr_name = "p_adj",
                           log2FC = cutoff_LFC, fdr = cutoff_padj)

df_to_plot$regulate%>%table
df_to_plot$Regulator%>%table

# df_to_plot = df_to_plot[df_to_plot$baseMean > as.numeric(quantile(df_to_plot$baseMean)['25%']),]


# plot
ggvolcano(df_to_plot, x = "log2FoldChange", y = "padj",
          label = "hgnc_symbol", 
          label_number = 0, 
          output = FALSE)
ggsave(file.path('results/RESTART/3.DEG_Brain/', 'Volcano_SCI_Aggregated_Sum.pdf'), width = 8, height = 8)
# ggsave(file.path('results/RESTART/3.DEG_Brain/', 'Volcano_SCI_Aggregated_Average.pdf'), width = 8, height = 8)




################################################################################
# 4. SCI-Pathway Analysis
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/3.DEG_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

require(DOSE)
require(clusterProfiler)

cutoff_padj = 0.05
cutoff_LFC = 1

results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))

up_host_gene = results_to_save[results_to_save$fold_change > cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
down_host_gene = results_to_save[results_to_save$fold_change < -cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
null_host_gene = results_to_save[results_to_save$p_adj > cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
universe_host_gene = na.omit(results_to_save$hgnc_symbol)


# GSA analysis
res_go = enrichGO(gene = up_host_gene, 
                  # universe = universe_host_gene,
                  keyType = 'SYMBOL', 
                  OrgDb = org.Hs.eg.db, 
                  ont = 'BP', 
                  pvalueCutoff = 0.25, 
                  qvalueCutoff = 1, 
                  readable = T,
                  pAdjustMethod="BH"
)
data.frame(res_go) %>% View
barplot(res_go, showCategory = 20)

# GSEA analysis
geneList = results_to_save$fold_change
names(geneList) = results_to_save$hgnc_symbol
geneList = sort(geneList, decreasing = T)
gse_go = gseGO(geneList = geneList, 
               keyType = 'SYMBOL',
               OrgDb = org.Hs.eg.db,
               ont = "BP", 
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               verbose = T)
saveRDS(gse_go, file = file.path('results/RESTART/3.DEG_Brain/', 'GSEA_BP_SCI_Aggregated_Sum.rds'))
dotplot(gse_go, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)
ggsave(file.path('results/RESTART/3.DEG_Brain/', 'GSEA_BP_SCI_Aggregated_Sum.pdf'), width = 12, height = 8)

?barplot

gse_go = gseGO(geneList = geneList, 
               keyType = 'SYMBOL',
               OrgDb = org.Hs.eg.db,
               ont = "CC", 
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               verbose = T)
saveRDS(gse_go, file = file.path('results/RESTART/3.DEG_Brain/', 'GSEA_CC_SCI_Aggregated_Sum.rds'))
dotplot(gse_go, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)
ggsave(file.path('results/RESTART/3.DEG_Brain/', 'GSEA_CC_SCI_Aggregated_Sum.pdf'), width = 12, height = 8)

gse_go = gseGO(geneList = geneList, 
               keyType = 'SYMBOL',
               OrgDb = org.Hs.eg.db,
               ont = "MF", 
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               verbose = T)
saveRDS(gse_go, file = file.path('results/RESTART/3.DEG_Brain/', 'GSEA_MF_SCI_Aggregated_Sum.rds'))
dotplot(gse_go, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)
ggsave(file.path('results/RESTART/3.DEG_Brain/', 'GSEA_MF_SCI_Aggregated_Sum.pdf'), width = 12, height = 8)

# KEGG
df_biomaRt = data.frame(fread(file = 'results/RESTART/3.DEG_Brain/biomaRt_240508.csv'))
df_biomaRt = df_biomaRt[c('hgnc_symbol','entrezgene_id','description')]
df_biomaRt = df_biomaRt[!duplicated(df_biomaRt),]

results_to_kegg = left_join(results_to_save, df_biomaRt, by = 'hgnc_symbol')
results_to_kegg = results_to_kegg[!is.na(results_to_kegg$entrezgene_id),]

geneList = results_to_kegg$fold_change
names(geneList) = results_to_kegg$entrezgene_id
geneList = sort(geneList, decreasing = T)

gse_kegg = gseKEGG(geneList = geneList, 
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   verbose = T)
saveRDS(gse_kegg, file = file.path('results/RESTART/3.DEG_Brain/', 'GSEA_KEGG_SCI_Aggregated_Sum.rds'))
dotplot(gse_kegg, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)
ggsave(file.path('results/RESTART/3.DEG_Brain/', 'GSEA_KEGG_SCI_Aggregated_Sum.pdf'), width = 12, height = 8)


################################################################################
# 4. RIP-seq
################################################################################
rm(list = ls())
gc()

dir.create('results/RESTART/4.RIP_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

require(clusterProfiler)
require(DOSE)
require(org.Hs.eg.db)

cutoff_padj = 0.05
cutoff_LFC = 1

# Load RIP-seq results
df_gene_bind = read.csv(file = 'results/4.ADARB2_RIP/ADAR3_Binding_RIP.txt', sep = '\t')
df_gene_bind%>%head
gene_bind = df_gene_bind$Gene


# DEG results
results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))

up_host_gene = results_to_save[results_to_save$fold_change > cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
down_host_gene = results_to_save[results_to_save$fold_change < -cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
null_host_gene = results_to_save[results_to_save$p_adj > cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
universe_host_gene = na.omit(results_to_save$hgnc_symbol)

# Only Genes included in Host Gene Universe
gene_bind_subset = gene_bind[gene_bind %in% universe_host_gene]

df_gene_bind_subset = df_gene_bind[df_gene_bind$Gene %in% gene_bind_subset,]
df_gene_bind_subset = df_gene_bind
df_gene_bind_subset[-c(1,2)]%>%head

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

melting_geneset = reshape2::melt(df_geneset)
colnames(melting_geneset) = c('hgnc_symbol','gs_name')
melting_geneset = melting_geneset[c(2,1)]
melting_geneset$gs_name%>%table
melting_geneset = data.frame(melting_geneset)

melting_geneset$gs_name%>%unique

# GSA
res_rip_up = enricher(gene = up_host_gene,
                   TERM2GENE = melting_geneset,
                   universe = universe_host_gene,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1,
                   pAdjustMethod = 'none',
                   minGSSize = 1,
                   maxGSSize = 5000)

df_to_plot = data.frame(res_rip_up)
df_to_plot$Description = factor(df_to_plot$Description, levels = df_to_plot$Description[order(df_to_plot$pvalue,decreasing = T)])
df_to_plot%>%View
saveRDS(df_to_plot, file = file.path('results/RESTART/4.RIP_Brain/', 'GSA_RIP.rds'))

# Create the barplot using ggplot2 with a border around the plot panel
ggplot(df_to_plot, aes(x = Count, y = Description, fill = pvalue)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene Count", y = "Description", fill="p-value") +
  scale_fill_gradient(low = "red", high = "blue") +
  theme_bw() + # white background
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5),  # Add border around the plot panel
        axis.text.y = element_text(family = "sans", color = "black", size = 12),  # Set font and size for y-axis 
        axis.text.x = element_text(family = "sans", color = "black", size = 12), # Set font and size for x-axis 
        axis.title = element_text(family = "sans", color = "black", size = 12),  # Set font and size for axis titles
        plot.title = element_text(family = "sans", color = "black", size = 12), # Set font and size for title
  )
ggsave(file.path('results/RESTART/4.RIP_Brain/GSA_RIP.pdf'), width = 8, height = 8)

# Venn Diagram
gene_list = list(Up = up_host_gene, RIP_inhibi = df_geneset$RIPseq_inhibitory)
ggVennDiagram(gene_list,
              category.names = c('Up-DEG',"Bound by ADAR3")) +
  scale_fill_gradient(low="grey90",high = "red") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  coord_flip()
ggsave(file.path('results/RESTART/4.RIP_Brain/Venn_DEG_ADARB2_RIP.pdf'), width = 8, height = 6)


# # GSEA
# row_filter = results_to_save$baseMean > as.numeric(quantile(results_to_save$baseMean, 0.95))
# 
# geneList = (results_to_save[row_filter,]$fold_change)
# names(geneList) = results_to_save[row_filter,]$hgnc_symbol
# geneList = sort(geneList, decreasing = T)
# 
# res_rip_gsea = GSEA(gene = geneList,
#                       TERM2GENE = melting_geneset,
#                       pvalueCutoff = 1,
#                       pAdjustMethod = 'none',
#                       minGSSize = 1,
#                       maxGSSize = 5000)
# View(data.frame(res_rip_gsea))


# Proportional Bar Plot
prop_up = prop.table(table((up_host_gene %in% df_geneset$RIPseq_inhibitory)))
prop_null = prop.table(table((universe_host_gene[!(universe_host_gene %in% up_host_gene)] %in% df_geneset$RIPseq_inhibitory)))


names(prop_up) = c('Not', 'Bound by ADAR3')
names(prop_null) = c('Not', 'Bound by ADAR3')
prop_up$Group = 'Up-DEG'
prop_null$Group = 'NULL'

prop_all = rbind(data.frame(prop_up), data.frame(prop_null))
prop_all = reshape2::melt(prop_all)

prop_all$Group = factor(prop_all$Group, levels = c('Up-DEG','NULL'))
prop_all$percent = round(prop_all$value * 100)
prop_all$`RIP-seq` = factor(str_replace_all(prop_all$variable,pattern = '\\.',replacement = ' '), levels = c('Bound by ADAR3', 'Not'))
prop_all$proportion = (prop_all$percent/100)

res_to_plot = data.frame(res_rip_up)
ggplot(prop_all, aes(x = Group, y = value, fill = `RIP-seq`)) + 
  geom_col(position = 'stack', colour="black")+
  geom_text(aes(label = proportion), 
            position = position_stack(vjust = .5), 
            color='black', 
            size=4) +    
  xlab('Differentialy Expressed Genes') +
  ylab('Proportion') +
  theme_minimal() +
  theme(
    # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
    plot.title = element_text(family = "sans", color = "black", size = 15, hjust = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
        ) +
  scale_fill_manual(values = c('red', 'lightgrey')) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ggtitle(paste0('p-value ',format.pval(res_to_plot$pvalue[res_to_plot$Description == 'RIPseq_inhibitory'], eps = 0.0001)))

################################################################################
# 5. irCLASH
################################################################################
rm(list = ls())
gc()

dir.create('results/RESTART/5.irCLASH_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

# Annotate whether circRNA are bound by ADAR
lapply(c('ADAR','ADARB1','ADARB2'), function(x){
  annotatedBSJs_withADARfam = data.frame(fread(str_glue('results/4.ADARB2_irCLASH/annotatedBSJs_with{x}.tsv')))
  
  annotatedBSJs_withADARfam$Group = apply(annotatedBSJs_withADARfam, MARGIN = 1, function(row_){
    #row_ = annotatedBSJs_withADARfam[4,]
    upIntronlist = row_[['Substrate.id_UpIntron']]
    downIntronlist = row_[['Substrate.id_DownIntron']]
    
    upIntronlist = strsplit(upIntronlist, ',')[[1]]
    downIntronlist = strsplit(downIntronlist, ',')[[1]]
    
    if (length(intersect(upIntronlist, downIntronlist)) != 0){
      return('Target')
    } else if (any(duplicated(upIntronlist))){
      return('UpIntron')
    } else if (any(duplicated(downIntronlist))){
      return('DownIntron')
    } else if ((length(upIntronlist) >= 1) | (length(downIntronlist) >= 1)){
      return('Moderate') 
    } else {
      return('Non_target')
    }
  }) %>% unlist
  annotatedBSJs_withADARfam$circRNA_id = lapply(str_split(annotatedBSJs_withADARfam$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
  return(annotatedBSJs_withADARfam)
}) -> list_df

# Subset only circRNA which was not filtered by previous
input_list = get_data(slot = 'rpm', annot = NULL, type = 'tissue', LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Filter out circRNAs with less than 20% of samples
cutoff_prop_samples = 0.2
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]

list_circRNA_id = cts_sci_filtered$circRNA_id

list_df[[1]]$Group%>%table
list_df[[2]]$Group%>%table
list_df[[3]]$Group%>%table

# Only Target both flanking introns
list_ADARbound_both = lapply(list_df, function(x){
  x = x[x$circRNA_id %in% list_circRNA_id,]
  # target = x$circRNA_id[!(x$Group %in% c('Non_target','Moderate'))]
  target = x$circRNA_id[x$Group == 'Target']
  return(target)
})

# ggsci::pal_npg()(9)%>%show_col()
# ggsci::pal_jama()(9)%>%show_col()
# ggsci::pal_jco()(9)%>%show_col()

# VennDiagram
ggVennDiagram(list_ADARbound_both,
              category.names = c('Bound by ADAR1',"Bound by ADAR2","Bound by ADAR3"),
              label_alpha = 0
              ) +
  scale_fill_gradient(low="grey90",high = "steelblue") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(file.path('results/RESTART/5.irCLASH_Brain/Venn_ADAR1_ADAR2_ADAR3_both.pdf'), width = 6, height = 6)


# Target at least one flanking intron
list_ADARbound_one = lapply(list_df, function(x){
  x = x[x$circRNA_id %in% list_circRNA_id,]
  target = x$circRNA_id[!(x$Group %in% c('Non_target','Moderate'))]
  # target = x$circRNA_id[x$Group == 'Target']
  return(target)
})

# VennDiagram
ggVennDiagram(list_ADARbound_one,
              category.names = c('Bound by ADAR1',"Bound by ADAR2","Bound by ADAR3"),
              label_alpha = 0
) +
  scale_fill_gradient(low="grey90",high = "steelblue") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(file.path('results/RESTART/5.irCLASH_Brain/Venn_ADAR1_ADAR2_ADAR3_one.pdf'), width = 6, height = 6)


names(list_ADARbound_both) = c('ADAR1_Both','ADAR2_Both','ADAR3_Both')
names(list_ADARbound_one) = c('ADAR1_One','ADAR2_One','ADAR3_One')

list_ADARbound = c(list_ADARbound_both, list_ADARbound_one)
list_ADARbound$ADAR_Both_Int = intersect(intersect(list_ADARbound$ADAR1_Both, list_ADARbound$ADAR2_Both), list_ADARbound$ADAR3_Both)
list_ADARbound$ADAR_Both_Union = union(union(list_ADARbound$ADAR1_Both, list_ADARbound$ADAR2_Both), list_ADARbound$ADAR3_Both)
list_ADARbound$ADAR_One_Int = intersect(intersect(list_ADARbound$ADAR1_One, list_ADARbound$ADAR2_One), list_ADARbound$ADAR3_One)
list_ADARbound$ADAR_One_Union = union(union(list_ADARbound$ADAR1_One, list_ADARbound$ADAR2_One), list_ADARbound$ADAR3_One)
list_ADARbound$ADAR_One_2and3 = intersect(list_ADARbound$ADAR2_One, list_ADARbound$ADAR3_One)
list_ADARbound$ADAR_One_1and2 = intersect(list_ADARbound$ADAR1_One, list_ADARbound$ADAR2_One)
saveRDS(list_ADARbound, file = 'results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')

# # Version 2
# list_ADARbound$ADAR3nADAR1uADAR2 = intersect(list_ADARbound$ADAR3_One, union(list_ADARbound$ADAR2_One, list_ADARbound$ADAR1_One))

# Load SCI DEG
cutoff_padj = 0.05
cutoff_LFC = 1

results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI.rds'))
results_to_save = results_to_save[!is.na(results_to_save$p_adj),]

# deg_sci[deg_sci$circRNA_id %in% list_ADARbound[[3]],]%>%View

up_circ_gene = results_to_save[results_to_save$fold_change > cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$circRNA_id %>% unique() %>% na.omit
down_circ_gene = results_to_save[results_to_save$fold_change < -cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$circRNA_id %>% unique() %>% na.omit
null_circ_gene = results_to_save[results_to_save$p_adj > cutoff_padj,]$circRNA_id %>% unique() %>% na.omit
universe_circ_gene = na.omit(results_to_save$circRNA_id)

melting_geneset = reshape2::melt(list_ADARbound)
colnames(melting_geneset) = c('hgnc_symbol','gs_name')
melting_geneset = melting_geneset[c(2,1)]
melting_geneset$gs_name%>%table
melting_geneset = data.frame(melting_geneset)


# GSA
res_clash_up = enricher(gene = up_circ_gene,
                        TERM2GENE = melting_geneset,
                        universe = universe_circ_gene,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = 'none',
                        minGSSize = 1,
                        maxGSSize = 5000)
print(data.frame(res_clash_up)[c('ID','Description','pvalue','qvalue')])

# GSEA
#row_filter = results_to_save$baseMean > as.numeric(quantile(results_to_save$baseMean, i)) # Gene Filtering
row_filter = rep(TRUE, nrow(results_to_save))

geneList = (results_to_save[row_filter,]$fold_change)*(-log10(results_to_save[row_filter,]$p_value))
# geneList = results_to_save[row_filter,]$fold_change
names(geneList) = results_to_save[row_filter,]$circRNA_id
geneList = sort(geneList, decreasing = T)

set.seed(42)
res_clash_gsea = GSEA(gene = geneList,
                      TERM2GENE = melting_geneset,
                      pvalueCutoff = 1,
                      pAdjustMethod = 'none',
                      minGSSize = 1,
                      maxGSSize = 5000)
print(data.frame(res_clash_gsea)[c("Description","ID","NES","pvalue")]%>%arrange(pvalue))
saveRDS(res_clash_gsea, file.path('results/RESTART/5.irCLASH_Brain/', 'GSEA_results.rds'))


res_clash_gsea = readRDS(file.path('results/RESTART/5.irCLASH_Brain/', 'GSEA_results.rds'))
data.frame(res_clash_gsea)[c("Description","ID","NES","pvalue")]%>%View
dotplot(res_clash_gsea, showCategory = 10, title = 'GSEA')

# melting_geneset_subset = melting_geneset[melting_geneset$gs_name %in% c('ADAR1_One','ADAR2_One','ADAR3_One'),]
# melting_geneset_subset$gs_name = str_replace_all(melting_geneset_subset$gs_name, pattern = '_One', replacement = '')
# 
# res_clash_gsea = GSEA(gene = geneList,
#                       TERM2GENE = melting_geneset_subset,
#                       pvalueCutoff = 1,
#                       pAdjustMethod = 'BH',
#                       minGSSize = 1,
#                       maxGSSize = 5000)
# print(data.frame(res_clash_gsea)[c("Description","ID","NES","pvalue")]%>%head)
# data.frame(res_clash_gsea)[c("Description","ID","NES","pvalue")]%>%View

res_clash_gsea@result%>%data.frame%>%View

gseaplot2(res_clash_gsea, 
          geneSetID = c('ADAR_One_Union','ADAR3_One', 'ADAR2_One','ADAR1_One'),
          #pvalue_table = TRUE,
          color = "black",
          subplots = 1, 
          base_size = 15) + 
  labs(title = "", 
       x = "Rank in Ordered Dataset", 
       y = "Running Enrichment Score") +
  scale_color_manual(name = 'Genesets', 
                    labels = c('All','ADAR1','ADAR2','ADAR3'),
                    values = 
                      # alpha(RColorBrewer::brewer.pal(n = 4, 'Dark2'), .7)
                      ggsci::pal_npg(alpha = 1)(4)
                    ) +
  theme(legend.position = "bottom",  # Adjusting legend position
        # legend.title = element_text(size = 12, face = "bold"),  # Adjusting legend title
        # legend.text = element_text(size = 12),  # Adjusting legend text size
        # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centering title
        # axis.title = element_text(size = 14),
        # axis.text = element_text(size = 12),
        plot.background = element_rect(fill = "white", color = NA),  # Setting plot background to white
        # text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12)
        ) +
  annotate("text", x = Inf, y = Inf, label = sprintf("\tAll, p-value=%.2f, NES=%.2f\nADAR3, p-value=%.2f, NES=%.2f\nADAR2, p-value=%.2f, NES=%.2f\nADAR1, p-value=%.2f, NES=%.2f",
                                                     data.frame(res_clash_gsea)['ADAR_One_Union','pvalue'],
                                                     data.frame(res_clash_gsea)['ADAR_One_Union','NES'],
                                                     data.frame(res_clash_gsea)['ADAR3_One','pvalue'],
                                                     data.frame(res_clash_gsea)['ADAR3_One','NES'],
                                                     data.frame(res_clash_gsea)['ADAR2_One','pvalue'],
                                                     data.frame(res_clash_gsea)['ADAR2_One','NES'],
                                                     data.frame(res_clash_gsea)['ADAR1_One','pvalue'],
                                                     data.frame(res_clash_gsea)['ADAR1_One','NES']
                                                     ),
         hjust = 1.1, vjust = 1.1, size = 3.4, color = "black")

# +
#   annotate("text", x = Inf, y = Inf, label = sprintf("ADAR1, p-value=%.2f, NES=%.2f\nADAR2, p-value=%.2f, NES=%.2f\nADAR3, p-value=%.2f, NES=%.2f\n\tAll, p-value=%.2f, NES=%.2f",
#                                                      data.frame(res_clash_gsea)['ADAR1_One','pvalue'],
#                                                      data.frame(res_clash_gsea)['ADAR1_One','NES'],
#                                                      data.frame(res_clash_gsea)['ADAR2_One','pvalue'],
#                                                      data.frame(res_clash_gsea)['ADAR2_One','NES'],
#                                                      data.frame(res_clash_gsea)['ADAR3_One','pvalue'],
#                                                      data.frame(res_clash_gsea)['ADAR3_One','NES'],
#                                                      data.frame(res_clash_gsea)['ADAR_One_Union','pvalue'],
#                                                      data.frame(res_clash_gsea)['ADAR_One_Union','NES']
#                                                      ),
#          hjust = 1.1, vjust = 1.1, size = 3.4, color = "black")

# ggsave(file.path('results/RESTART/5.irCLASH_Brain//GSEA_ADAR.pdf'), width = 6, height = 6)

################################################################################
# 6. FeaturePlot (SCI)
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/2.EDA', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'cell type'

# Load PCA coordinates
pca_res = readRDS(file = str_glue('results/RESTART/6.FeaturePlot/pca_res_{type}.rds'))
pca_df = readRDS(file = str_glue('results/RESTART/6.FeaturePlot/pca_df_{type}.rds'))

# SCI Plots
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]
# Filter out circRNAs with less than 20% of samples
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]

# SCI values
meta_to_add = data.frame(Sample = colnames(cts_sci_filtered[-1]), Scaled_SCI = scale(colSums(cts_sci_filtered[-1], na.rm = T)))
pca_df_merged = left_join(pca_df, meta_to_add, by = c('Sample' = 'Sample'))

ggplot(pca_df_merged, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(color='black', aes(fill = Scaled_SCI), size = 3, shape = 21) +
  xlab(paste("Principal Component 1", sep="")) +
  ylab(paste("Principal Component 2", sep="")) +
  ggtitle("Scaled SCI by tissue") +
  theme_minimal() +
  geom_text_repel() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(file.path('results/RESTART/6.FeaturePlot/', str_glue('SCI_SumScaled_{type}.pdf')), width = 6, height = 6)


################################################################################
# 6. FeaturePlot (ADAR Expression)
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/2.EDA', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'cell type'

# Load PCA coordinates
pca_res = readRDS(file = str_glue('results/RESTART/6.FeaturePlot/pca_res_{type}.rds'))
pca_df = readRDS(file = str_glue('results/RESTART/6.FeaturePlot/pca_df_{type}.rds'))

# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)
cts = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts = cts[!(rowSums(cts[-1]) == 0),]
rownames(cts) = cts$gene
cts = cts[,-1]


# # Load Raw Linear Counts from Poly-A
# list_input = get_data_linear(slot = 'count', type = type, polyA = TRUE) 
# cts = list_input[[1]]; meta = list_input[[2]]; rownames(meta) = meta$name; rm(list_input)
# meta$group = ifelse(grepl('Brain', meta$name, ignore.case = T), 'Brain', 'Others') %>%
#   factor(levels = c('Others', 'Brain'))
# 
# # DESeq2 Normalization
# require(DESeq2)
# dds = DESeqDataSetFromCountsMatrix(sampleTable = meta, cts = cts, filtering = T, design = ~ group) # Only Zero Filtering
# dds = DESeq(dds)
# 
# dir.create('results/RESTART/7.PolyA', showWarnings = F)
# saveRDS(dds, file.path('results/RESTART/7.PolyA/', str_glue('dds_{type}_polyA.rds')))
# 
# # Extract Normalized Counts
# cts = counts(dds, normalized = T)
# 
# cts = data.frame(cts)
# colnames(cts) = str_replace_all(colnames(cts), pattern = '\\.', replacement = ' ')

list_features = rownames(cts)[str_detect(rownames(cts), 'ADAR')]
meta_to_add = data.frame(t(cts[list_features,]))
meta_to_add = log2(meta_to_add + 1) # Log Transform

colnames(meta_to_add) = str_replace_all(colnames(meta_to_add), pattern = '\\.', replacement = '-')
# rownames(meta_to_add) = colnames(dds) # In case of cell type data
meta_to_add$Sample = rownames(meta_to_add)

# Normalized Counts
pca_df_merged = left_join(pca_df, meta_to_add, by = c('Sample' = 'Sample'))

list_p = lapply(list_features[-length(list_features)], function(feature){
  ggplot(pca_df_merged, aes(x = Dim.1, y = Dim.2, label = Sample)) +
    geom_point(color='black', aes(fill = .data[[feature]]), size = 3, shape = 21) +
    xlab(paste("Principal Component 1", sep="")) +
    ylab(paste("Principal Component 2", sep="")) +
    theme_minimal() +
    geom_text_repel() +
    scale_fill_gradient(low = "white", high = "blue") +
    ggtitle(label = feature) +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(fill = '')
})
p_all = ggarrange(plotlist = list_p, ncol = 3, legend = 'right', common.legend = TRUE)
p_all
ggsave(filename = file.path('results/RESTART/6.FeaturePlot/', str_glue('ADAR_{type}_polyA_count.pdf')), plot = p_all, width = 18, height = 6)


# Scaled Counts
meta_to_add_scaled = data.frame(scale(meta_to_add[-ncol(meta_to_add)]))
colnames(meta_to_add_scaled) = str_replace_all(colnames(meta_to_add_scaled), pattern = '\\.', replacement = '-')
meta_to_add_scaled$Sample = meta_to_add$Sample
pca_df_merged = left_join(pca_df, meta_to_add_scaled, by = c('Sample' = 'Sample'))

list_p = lapply(list_features[-length(list_features)], function(feature){
  ggplot(pca_df_merged, aes(x = Dim.1, y = Dim.2, label = Sample)) +
    geom_point(color='black', aes(fill = .data[[feature]]), size = 3, shape = 21) +
    xlab(paste("Principal Component 1", sep="")) +
    ylab(paste("Principal Component 2", sep="")) +
    theme_minimal() +
    geom_text_repel() +
    scale_fill_gradient2(low = "blue", mid = 'white',high = "red", midpoint = 0) +
    ggtitle(label = feature) +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(fill = '')
})
p_all = ggarrange(plotlist = list_p, ncol = 3, legend = 'right', common.legend = TRUE)
p_all
ggsave(filename = file.path('results/RESTART/6.FeaturePlot/', str_glue('ADAR_{type}_polyA_scaled.pdf')), plot = p_all, width = 18, height = 6)


################################################################################
# 8. Correlation analysis (SCI and RNA binding Protein)
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/8.Cor_RBP', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples

# Sum of SCI
sci_sum = data.frame(Sample = colnames(cts_sci_filtered[-1]), Sum_SCI = colSums(cts_sci_filtered[-1], na.rm = T))
sci_sum$Sum_SCI = log2(sci_sum$Sum_SCI + 1)

# # Load Linear Data (Poly-A) -> DESeq2 Normalized
# dds = readRDS(file.path('results/RESTART/7.PolyA/', str_glue('dds_{type}_polyA.rds')))
# cts_linear = counts(dds, normalized = T) # Extract Normalized Counts
# cts_linear = data.frame(cts_linear)
# colnames(cts_linear) = str_replace_all(colnames(cts_linear), pattern = '\\.', replacement = ' ')
# 
# # Only RBP lists
# genelist_RBP = fread(file = 'RBP/NuclearRBP.txt', header = F) %>% unlist
# cts_linear_subset = cts_linear[rownames(cts_linear) %in% genelist_RBP,]
# cts_linear_subset = log2(cts_linear_subset + 1)


# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)
cts_linear = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts_linear = cts_linear[!(rowSums(cts_linear[-1]) == 0),]
rownames(cts_linear) = cts_linear$gene
cts_linear = cts_linear[,-1]

# Only RBP lists
genelist_RBP = fread(file = 'RBP/NuclearRBP.txt', header = F) %>% unlist
cts_linear_subset = cts_linear[rownames(cts_linear) %in% genelist_RBP,]
cts_linear_subset = log2(cts_linear_subset + 1)


# Calculate Correlation with RBP & Sum of SCI
input_x = as.matrix(t(cts_linear_subset)) ; input_x[1:4,1:4]
input_y = as.matrix(data.frame(sci_sum[-1])) ; input_y%>%head

# Sanity Check
require(Hmisc)
corr_matrix = rcorr(x = input_x, 
                    y = input_y,
                    type = 'pearson')

corr_matrix%>%names
coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

res_cor = data.frame(coeffs, pvalues)
res_cor$padj = p.adjust(res_cor$pvalues, method = 'BH')
saveRDS(res_cor, file = file.path('results/RESTART/8.Cor_RBP', str_glue('res_cor_{type}_polyA.rds')))


# Volcano Plot
log2FC_cutoff = min(abs(res_cor$coeffs[(res_cor$padj < 0.05)]))
df_to_plot <- add_regulate(res_cor, log2FC_name = "coeffs",fdr_name = "padj",
                           log2FC = log2FC_cutoff,
                           fdr = 0.05)
df_to_plot$label = rownames(df_to_plot)
df_to_plot$label = ifelse(df_to_plot$label %in% c('ADAR', 'ADARB1','ADARB2'), df_to_plot$label, NA)


# Generate the plot
df_to_plot$label = rownames(df_to_plot)
df_to_plot$label = ifelse(df_to_plot$label %in% c('ADAR', 'ADARB1','ADARB2'), df_to_plot$label, NA)
df_to_plot$label = factor(df_to_plot$label, levels = c('ADAR', 'ADARB1', 'ADARB2'))

# Define colors
color_map <- c("Up" = "#E41A1C", "Down" = "blue", "Normal" = "grey",
               "ADAR" = "#E64B3599", "ADARB1" =  "#4DBBD599" , "ADARB2" = "#00A08799"
               )

ggplot(df_to_plot, aes(x = log2FoldChange, y = -log10(padj), label = label)) +
  geom_point(aes(color = regulate),
             # position=position_jitter(h=0.01, w=0.01)
              ) +
  
  scale_color_manual(values = color_map) +
  theme_classic() +
  labs(x = 'Correlation Coefficient', y = '-Log10 FDR')

show_col(pal_npg(alpha = 0.6)(3))
pal_npg(alpha = 0.6)(3)
# # Volcano Plot
# ggplot(df_to_plot, aes(x = log2FoldChange, y = -log10(padj), label = label)) +
#   geom_point(aes(color = regulate), alpha = ifelse(df_to_plot$label %in% c("ADAR", "ADARB1", "ADARB2"), 1, 0.1),
#              size = ifelse(df_to_plot$label %in% c("ADAR", "ADARB1", "ADARB2"), 3, 1.5)) +
#   geom_point(data = subset(df_to_plot, label %in% c("ADAR", "ADARB1", "ADARB2")),
#              aes(color = label), size = 3) +
#   labs(x = 'Correlation Coefficient', y = '-Log10 FDR') +
#   theme_classic() +
#   theme(legend.position = "none",
# #        plot.title = element_text(hjust = 0.5, size = 15),
#         text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12)
#         ) +
#   scale_color_manual(values = color_map) +
#   geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.5) +
#   geom_text_repel(data = subset(df_to_plot, label %in% c("ADAR", "ADARB1", "ADARB2")),
#                   aes(color = label),
#                   nudge_x = 0.3,
#                   nudge_y = -0.1,
#                   size = 5,
#                   max.overlaps = 10) 

# Scatter Plot
input_x[1:4,1:4]
input_y%>%head


df_to_scatter = cbind(data.frame(input_x), data.frame(input_y))
df_to_scatter$Sample = rownames(df_to_scatter)
df_to_scatter = left_join(df_to_scatter, meta, by = c('Sample' = 'name'))

# Scatterplot
feature = 'ADAR'
rownames(df_to_scatter) = df_to_scatter$Sample

df_to_scatter$organ_system

ggplot(df_to_scatter, aes(x = Sum_SCI, y = .data[[feature]])) +
  geom_smooth(method = 'lm', se = F, color = 'black') +
  geom_point(alpha = 0.5, size = 1.5) +
  labs(x = 'SCI by Tissues', y = str_glue('{feature} Expression')) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
        text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12)
        ) +
  geom_text_repel(data = subset(df_to_scatter, (ADAR > median(df_to_scatter[[feature]])) & (Sum_SCI > median(df_to_scatter[['Sum_SCI']]))),
                  aes(label = rownames(subset(df_to_scatter, (ADAR > median(df_to_scatter[[feature]])) & (Sum_SCI > median(df_to_scatter[['Sum_SCI']]))))),
                  # nudge_x = 0.5,
                  nudge_y = 0.2,
                  size = 4,
                  max.overlaps = 10) +
  ggtitle(feature)

lapply(c('ADAR', 'ADARB1', 'ADARB2'), function(feature) {
  # Filtering data for text labels before passing to ggplot to simplify the plot code
  # filtered_data <- df_to_scatter[(df_to_scatter[[feature]] > quantile(df_to_scatter[[feature]], 0.75) & df_to_scatter$Sum_SCI > quantile(df_to_scatter$Sum_SCI,  0.75)) | (str_detect(df_to_scatter$Sample, 'Brain|brain')), ]

  filtered_data <- df_to_scatter[df_to_scatter$Sample %in% c('Neuron','Oligodendrocyte Precursor Cell'),] # CEll type
  ggplot(df_to_scatter, aes_string(x = "Sum_SCI", y = feature)) +
    geom_smooth(method = 'lm', se = F, color = 'black') +
    geom_point(alpha = 0.5, size = 2) +
    labs(x = 'SCI by Tissues', y = paste(feature, "Expression")) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
          text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12)
    ) +
    geom_text_repel(data = filtered_data,
                    aes_string(label = "rownames(filtered_data)"),
                    nudge_y = 0.1,
                    #size = 4,
                    max.overlaps = 5
                    ) +
    ggtitle(feature) -> p
  
  return(p)
}) -> list_p

p_all = ggarrange(plotlist = list_p, ncol = 3)
p_all
ggsave(file.path(str_glue('results/RESTART/8.Cor_RBP/Scatter_ADAR_{type}.pdf')), p_all, width = 18, height = 6)


################################################################################
# 9. circRNA-host RNA correlation
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/9.Cor_Host', showWarnings = F)
source('scripts/manuscript/Utils.R')

type = 'tissue'

# Load circRNA Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = F)
cts = input_list[[1]] ;  meta = input_list[[2]]
cts[-1] = log2(cts[-1] + 1)

# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = F)
cts_linear = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts_linear = cts_linear[!(rowSums(cts_linear[-1]) == 0),]
cts_linear[-1] = log2(cts_linear[-1] + 1)


# Only Use Ovelapped Samples
list_overlap = intersect(colnames(cts)[-1], colnames(cts_linear)[-1])

cts = cbind(cts[1],cts[-1][list_overlap])
cts_linear = cbind(cts_linear[1],cts_linear[-1][list_overlap])

# Sanity Check
if (!all(colnames(cts_linear)[-1]==colnames(cts)[-1])){
  stop()
}

# Ordering
annotatedBSJs_withAlu = data.frame(fread(file = 'results/circRNAprofiler/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
list_circ = cts$circRNA_id
list_Host = sapply(cts$circRNA_id, function(x){
  annotatedBSJs_withAlu$gene[annotatedBSJs_withAlu$circRNA_id == x]
})

# Input for calculating correlation
cts_linear_ordered = lapply(list_Host, function(x){
  cts_linear[cts_linear$gene == x, ]
}) %>% do.call(rbind, .)

circRNA_withGene = lapply(list_Host, function(x){
  x %in% cts_linear$gene
}) %>% unlist

cts_ordered = cts[circRNA_withGene,]

# Calculate Correlation
list_samples = colnames(cts[-1])
res_cor = list()

# filter_linear = !(as.numeric(cts_linear_ordered[,sample]) == 0)
# filter_circular = !(as.numeric(cts_ordered[,sample]) == 0)
# filter_both = filter_linear & filter_circular

for (i in 1:length(list_samples)){
  sample = list_samples[i]
  res_cor[[sample]] = cor.test(x = as.numeric(cts_linear_ordered[,sample]), #[filter_both], 
                               y = as.numeric(cts_ordered[,sample]), #[filter_both], 
                               method = 'pearson', 
                               use="complete.obs")
}

df_cor = lapply(res_cor, function(x){
  data.frame(coef = x$estimate, p.value = x$p.value)
}) %>% do.call(rbind, .)


df_cor$tissue = rownames(df_cor)
df_cor$name = df_cor$tissue
df_cor = left_join(df_cor, meta, by = 'name')
fwrite(df_cor, file = str_glue('results/RESTART/9.Cor_Host/Cor_Host_{type}.tsv'), sep = '\t', quote = F, row.names = F)


ggplot(df_cor, aes(x = reorder(name, coef), y = coef)) +
  geom_bar(stat = 'identity', position = 'dodge', colour="black", width=1) +
  labs(title = 'circRNA-Host Correlation', x = 'Sample', y = 'Correlation Coefficient') +
  # theme_classic2()+
  theme_bw() + # white background
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)
        # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12),
        # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
        # axis.line.x = element_line(color="black", size = 0.5),
        # axis.line.y = element_line(color="black", size = 0.5)
  ) +
  coord_flip()
################################################################################
# 9. SCI-host RNA correlation
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/9.Cor_Host_SCI', showWarnings = F)
source('scripts/manuscript/Utils.R')

type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts = input_list[[1]] ;  meta = input_list[[2]]
cts = cts[apply(cts[-1], 1, function(x){sum(!is.na(x))/length(x) >= 1}),] 
cts[-1] = log2(cts[-1] + 1)

# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = F)
cts_linear = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts_linear = cts_linear[!(rowSums(cts_linear[-1]) == 0),]
cts_linear[-1] = log2(cts_linear[-1] + 1)


# Only Use Ovelapped Samples
list_overlap = intersect(colnames(cts)[-1], colnames(cts_linear)[-1])

cts = cbind(cts[1],cts[-1][list_overlap])
cts_linear = cbind(cts_linear[1],cts_linear[-1][list_overlap])

# Sanity Check
if (!all(colnames(cts_linear)[-1]==colnames(cts)[-1])){
  stop()
}

# Ordering
annotatedBSJs_withAlu = data.frame(fread(file = 'results/circRNAprofiler/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
list_circ = cts$circRNA_id
list_Host = sapply(cts$circRNA_id, function(x){
  annotatedBSJs_withAlu$gene[annotatedBSJs_withAlu$circRNA_id == x]
})

# Input for calculating correlation
cts_linear_ordered = lapply(list_Host, function(x){
  cts_linear[cts_linear$gene == x, ]
}) %>% do.call(rbind, .)

circRNA_withGene = lapply(list_Host, function(x){
  x %in% cts_linear$gene
}) %>% unlist

cts_ordered = cts[circRNA_withGene,]

# Calculate Correlation
list_samples = colnames(cts[-1])
res_cor = list()

for (i in 1:length(list_samples)){
  sample = list_samples[i]
  res_cor[[sample]] = cor.test(x = as.numeric(cts_linear_ordered[,sample]), 
                               y = as.numeric(cts_ordered[,sample]), 
                               method = 'pearson', 
                               use="complete.obs")
}

df_cor = lapply(res_cor, function(x){
  data.frame(coef = x$estimate, p.value = x$p.value)
}) %>% do.call(rbind, .)


df_cor$tissue = rownames(df_cor)
df_cor$name = df_cor$tissue
df_cor = left_join(df_cor, meta, by = 'name')
fwrite(df_cor, file = 'results/RESTART/9.Cor_Host_SCI/Cor_Host_tissue.tsv', sep = '\t', quote = F, row.names = F)

df_cor = data.frame(fread(file = 'results/RESTART/9.Cor_Host_SCI/Cor_Host_tissue.tsv', sep = '\t'))
ggplot(df_cor, aes(x = reorder(name, coef), y = coef, fill = organ_system)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'circRNA-Host Correlation', x = 'Sample', y = 'Correlation Coefficient') +
  # theme_classic2()+
  theme_bw() + # white background
  scale_fill_manual(values = color_map) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)
        # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12),
        # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
        # axis.line.x = element_line(color="black", size = 0.5),
        # axis.line.y = element_line(color="black", size = 0.5)
  ) +
  scale_y_continuous(expand = c(NA, 0), limits = c(NA, 0.25)) +
  coord_flip()


# df_cor = df_cor[!(grepl(df_cor$name, pattern = 'Brain',ignore.case = T)),]
df_cor$Group = ifelse(grepl(df_cor$name, pattern = 'Brain', ignore.case = T), 'Brain', 'Others')
df_cor$Group = ifelse(grepl('Heart|ventricle|atruim|atrium', df_cor$name, ignore.case = T), 'Heart',
                      ifelse(df_cor$Group == 'Brain', 'Brain','Others')) %>%
  factor(levels = c('Brain','Heart', 'Others'))

library(scales)
hex <- hue_pal()(15)
color_map = c(hex[2:3], 'grey')
names(color_map) = c('Heart','Brain','Others')

library(ggsignif)
ggplot(df_cor, aes(y = coef, x = Group, fill = Group)) +
  geom_boxplot() +
  labs(
    # title = 'circRNA-Host Correlation', 
    x = 'Sample',
    y = 'Correlation Coefficient') +
  theme(legend.position = 'bottom') +
  # stat_compare_means() +
  scale_fill_manual(values = color_map) +
  theme_bw() + # white background
  theme(legend.position = 'right',
        # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12),
        # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
        # axis.line.x = element_line(color="black", size = 0.5),
        # axis.line.y = element_line(color="black", size = 0.5)
  )+ geom_signif(comparisons = list(c("Brain", "Heart"), c("Heart", "Others"), c("Brain", "Others")),
                 map_signif_level = TRUE, test = "wilcox.test", paired = TRUE)


################################################################################
# 9. SCI-host RNA Regression
################################################################################
rm(list=ls())
gc()
dir.create('results/RESTART/9.Cor_Host_SCI', showWarnings = F)
source('scripts/manuscript/Utils.R')

type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts = input_list[[1]] ;  meta = input_list[[2]]
cts = cts[apply(cts[-1], 1, function(x){sum(!is.na(x))/length(x) >= 1}),] 
cts[-1] = log2(cts[-1] + 1)

# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = F)
cts_linear = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts_linear = cts_linear[!(rowSums(cts_linear[-1]) == 0),]
cts_linear[-1] = log2(cts_linear[-1] + 1)


# Only Use Ovelapped Samples
list_overlap = intersect(colnames(cts)[-1], colnames(cts_linear)[-1])

cts = cbind(cts[1],cts[-1][list_overlap])
cts_linear = cbind(cts_linear[1],cts_linear[-1][list_overlap])

# Sanity Check
if (!all(colnames(cts_linear)[-1]==colnames(cts)[-1])){
  stop()
}

# Ordering
annotatedBSJs_withAlu = data.frame(fread(file = 'results/circRNAprofiler/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
list_circ = cts$circRNA_id
list_Host = sapply(cts$circRNA_id, function(x){
  annotatedBSJs_withAlu$gene[annotatedBSJs_withAlu$circRNA_id == x]
})

# Input for calculating correlation
cts_linear_ordered = lapply(list_Host, function(x){
  cts_linear[cts_linear$gene == x, ]
}) %>% do.call(rbind, .)

circRNA_withGene = lapply(list_Host, function(x){
  x %in% cts_linear$gene
}) %>% unlist

cts_ordered = cts[circRNA_withGene,]

# Calculate Regression
list_samples = colnames(cts[-1])
res_reg = list()

for (i in 1:length(list_samples)){
  sample = list_samples[i]
  # Univariate Regression
  dat = data.frame(linear = as.numeric(cts_linear_ordered[,sample]),
                   circ = as.numeric(cts_ordered[,sample]))
  res_reg[[sample]] <- lm(circ ~ linear, data=dat)
}


df_reg = lapply(res_reg, function(lm_fit){
  return(data.frame(adj.r.squared = summary(lm_fit)$adj.r.squared,
                    coef = summary(lm_fit)$coefficients['linear','Estimate'],
                    coef.p_val = summary(lm_fit)$coefficients['linear','Pr(>|t|)'],
                    model.p_val = as.numeric(broom::glance(lm_fit)['p.value']))
  )
}) %>% do.call(rbind, .)


df_reg$tissue = rownames(df_reg)
df_reg$name = df_reg$tissue

df_reg = left_join(df_reg, meta, by = 'name')
fwrite(df_reg, file = 'results/RESTART/9.Cor_Host_SCI/Reg_Host_tissue.tsv', sep = '\t', quote = F, row.names = F)

df_reg = fread(file = 'results/RESTART/9.Cor_Host_SCI/Reg_Host_tissue.tsv', sep = '\t', header = T)
color_map = readRDS(file = str_glue('results/RESTART/X.Figure/Robject/color_map_scale.v2.rds'))

ggplot(df_reg, aes(x = reorder(name, coef), y = coef, fill = organ_system)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'circRNA-Host Regression', x = 'Sample', y = 'Beta Coefficient') +
  # theme_classic2()+
  theme_bw() + # white background
  scale_fill_manual(values = color_map) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)
  ) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 0.25)) +
  coord_flip()


df_reg$Group = ifelse(grepl(df_reg$name, pattern = 'Brain', ignore.case = T), 'Brain', 'Others')
df_reg$Group = ifelse(grepl('Heart|ventricle|atruim|atrium', df_reg$name, ignore.case = T), 'Heart', 
                      ifelse(df_reg$Group == 'Brain', 'Brain','Others')) %>%
  factor(levels = c('Brain','Heart', 'Others'))

color_map = c(color_map[c('Cardiac System', 'Central Nervous System')], 'grey')
names(color_map) = c('Heart','Brain','Others')

library(ggsignif)
ggplot(df_reg, aes(y = coef, x = Group, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    # title = 'circRNA-Host Correlation', 
    x = 'Sample',
    y = 'Beta Coefficient') +
  theme(legend.position = 'bottom') +
  # stat_compare_means() +
  scale_fill_manual(values = color_map) +
  theme_bw() + # white background
  theme(legend.position = 'right',
        # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
        # axis.line.x = element_line(color="black", size = 0.5),
        # axis.line.y = element_line(color="black", size = 0.5)
  )+ geom_signif(comparisons = list(c("Brain", "Heart"), c("Heart", "Others"), c("Brain", "Others")),
                 map_signif_level = TRUE, test = "wilcox.test")

ggplot(df_reg, aes(y = adj.r.squared, x = Group, fill = Group)) +
  geom_boxplot() +
  labs( x = 'Sample', y = 'adjusted R squared') +
  theme(legend.position = 'bottom') +
  scale_fill_manual(values = color_map) +
  theme_bw() + # white background
  theme(legend.position = 'right',
        # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # panel.border = element_rect(color = "black", fill = NA, linewidth  = 0.5)
        # axis.line.x = element_line(color="black", size = 0.5),
        # axis.line.y = element_line(color="black", size = 0.5)
  )+ geom_signif(comparisons = list(c("Brain", "Heart"), c("Heart", "Others"), c("Brain", "Others")),
                 map_signif_level = TRUE, test = "wilcox.test")


sample = 'esophagus'
dat = data.frame(x = as.numeric(cts_linear_ordered[,sample]),
                 y = as.numeric(cts_ordered[,sample]))

# with correlation coefficient
p = ggplot(dat[ntile(dat$x,n = 4) == 4,], aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = sample, x = 'linearRNA', y = 'circRNA') +
  theme_classic()

library(ggpmisc)
p
p + stat_poly_eq(use_label(c("eq", "R2")), label.x = "right", label.y = "top")