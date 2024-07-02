# Load library
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

# Load circRNA only Annotated
annot = data.frame(fread(file = 'results/annotatedBSJs_withAlu_IntronPair.tsv',sep = '\t'))[c('circRNA_id','gene')]
cts_rpm = read.csv(file = 'data/GSE138734_RNAAtlas_circRNA_backspliced_rpm_filtered.txt', sep = '\t')

cts_rpm = cts_rpm[cts_rpm$circRNA_id %in% annot$circRNA_id,]

# Load metadata
meta = read.csv(file = 'data/GSE138734_metadata_v2.txt', sep = '\t')
names(meta)[names(meta) == 'RNA_Atlas_Number'] = 'sample'

# Filter Meta data
samples = names(cts_rpm)[!(names(cts_rpm) %in% c('circRNA_id'))]
meta = meta[meta$sample %in% samples,]
meta = meta %>% arrange(sample)

# Merge circRNA cts with Gene Symbol
cts_rpm = left_join(cts_rpm, annot, by = 'circRNA_id')

# Load linear RNA cts
rna_tpm = data.frame(fread(file = 'data/GSE138734_RNAAtlas_totalRNA_tpm_filtered.txt', sep= '\t'))
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
write.csv(meta_tpm, file = 'results/AverageTPM_linearRNA.csv', row.names = F)

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
write.csv(meta_rpm, file = 'results/AverageRPM_circRNA.csv', row.names = F)


################################################################################
# 1. Gene Filtering or Generate NA values
################################################################################
cutoff_tpm_elements = 1

# NA value if linear RNA expression is under 1 TPM
rna_tpm[-1][rna_tpm[-1] <= cutoff_tpm_elements] = NA

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

# Save Results
fwrite(LCratio, file = 'results/LCratio_withNAandPseudocount.txt', sep = '\t')

LC_ratio = data.frame(LC_ratio)
LC_ratio[-ncol(LC_ratio)] = LC_ratio[-ncol(LC_ratio)] * 1e6
LC_ratio = LC_ratio[c(ncol(LC_ratio), 1:(ncol(LC_ratio)-1))]

if (!(names(LC_ratio)[1] == 'circRNA_id')){
  stop('The first column should be circRNA_id')
}

fwrite(LC_ratio, 'results/LCratio_withNAandPseudocount_Milion.txt', sep = '\t')

################################################################################
# 2. PCA analysis
################################################################################
source('scripts/manuscript/Utils.R')

require(factoextra)
require(missMDA)
require(FactoMineR)

cutoff_prop_samples = 0.2
type = 'tissues'
title = 'circRNA SCI values'

# SCI Plots
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Filter out circRNAs with less than 20% of samples
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]

log_rna_matrix = log2(cts_sci_filtered[-1] + 1)
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]

nFeatures = nrow(log_rna_matrix)

# Transpose matrix to have samples as rows and features as columns
log_rna_matrix_t <- data.frame(t(log_rna_matrix))

# Imputation
imputed_data <- missMDA::imputePCA(log_rna_matrix_t, ncp = 2)

# Perform PCA on the imputed data
pca_res = FactoMineR::PCA(imputed_data$completeObs, graph = FALSE, ncp = 2, scale.unit = TRUE)

# Prepare data for ggplot
pca_df = as.data.frame(pca_res$ind$coord)
pca_df$Sample = rownames(pca_df)
pca_df = left_join(pca_df, meta, by = c("Sample" = "name"))

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

################################################################################
# 2. Tissue Specificity (Silhouette Score by Tissues)
################################################################################
library(cluster)

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
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]
p_sci = Boxplot_Silhouette(log_rna_matrix, meta, 'SCI by tissue')

# Linear TPM
annot = read.csv('data/annotatedBSJs.csv', row.names = 1) %>% na.omit()
annot$circRNA_id = lapply(str_split(annot$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
circRNA_id_Tissue = cts_sci_filtered$circRNA_id

annot = annot[annot$circRNA_id %in% circRNA_id_Tissue,]

input_list = get_data_linear(slot = 'data', annot = annot, type = type)
cts = input_list[[1]] ;  meta = input_list[[2]]; rm()

log_rna_matrix = log2(cts[-1] + 1)
log_rna_matrix = na.omit(log_rna_matrix)
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]
p_linear = Boxplot_Silhouette(log_rna_matrix, meta, 'Linear TPM by tissue')

# Circular RPM
list_input = get_data(slot = 'rpm', type = type, LCratio = F, annot = annot)
cts = list_input[[1]]; meta = list_input[[2]]; rownames(meta) = meta$name; rm(list_input)

log_rna_matrix = log2(cts[-1] + 1)
log_rna_matrix = na.omit(log_rna_matrix)
log_rna_matrix = log_rna_matrix[rowSums(log_rna_matrix, na.rm = T) > 0,]
log_rna_matrix = na.omit(log_rna_matrix)

p_circle = Boxplot_Silhouette(log_rna_matrix, meta, 'Circular RPM by tissue')

combined = ggarrange(plotlist = list(p_sci, p_linear, p_circle), ncol = 3, nrow = 1, common.legend = T)


################################################################################
# 9. circRNA-host RNA correlation
################################################################################
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
annotatedBSJs_withAlu = data.frame(fread(file = 'data/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
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

################################################################################
# 9. SCI-host RNA correlation
################################################################################
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
annotatedBSJs_withAlu = data.frame(fread(file = 'data/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
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

################################################################################
# SCI-host RNA Regression
################################################################################
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
annotatedBSJs_withAlu = data.frame(fread(file = 'data/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
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


################################################################################
# DEG analysis from SCI values (Gene-level Aggregation)
################################################################################
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

cts = cts_sci

# Aggregate circRNA by gene-level
annotatedBSJs = data.frame(fread(file = 'data/annotatedBSJs.csv'))[-1]
annotatedBSJs$circRNA_id = lapply(str_split(annotatedBSJs$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
names(annotatedBSJs)[names(annotatedBSJs) == 'gene'] = 'hgnc_symbol'
cts = left_join(cts, annotatedBSJs[c('circRNA_id','hgnc_symbol')], by = 'circRNA_id')

# Aggregate values row-wise by hgnc_symbol
custom_sum <- function(x) {
  if (all(is.na(x))) {
    NA
  } else {
    sum(x, na.rm = TRUE)
  }
}

cts_agg = cts %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = T), .groups = "drop")


# Filter out circRNAs with less than 20% of samples
cutoff_prop_samples = 0.2
cts_agg_filtered= cts_agg[apply(cts_agg[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),]

list_symbol_id = cts_agg_filtered$hgnc_symbol
cts_agg_filtered$hgnc_symbol = NULL ; cts_agg_filtered[1:4,1:4]

cts_agg_filtered = data.frame(t(cts_agg_filtered))
cts_agg_filtered$name = rownames(cts_agg_filtered)
cts_agg_filtered = left_join(cts_agg_filtered, meta[c('name','group')], by = 'name')

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

# Functional Annotation
cutoff_padj = 0.05
cutoff_LFC = 1

Nuclear_RBP = data.frame(fread(file = 'data/NuclearRBP.txt', sep = '\t', header = F))
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
saveRDS(results_to_save, file.path('results/', 'DEG_SCI_Aggregated_Sum.rds'))

################################################################################
# SCI-Pathway Analysis
################################################################################

require(DOSE)
require(clusterProfiler)

cutoff_padj = 0.05
cutoff_LFC = 1

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
dotplot(gse_go, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)

gse_go = gseGO(geneList = geneList, 
               keyType = 'SYMBOL',
               OrgDb = org.Hs.eg.db,
               ont = "CC", 
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               verbose = T)
dotplot(gse_go, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)

gse_go = gseGO(geneList = geneList, 
               keyType = 'SYMBOL',
               OrgDb = org.Hs.eg.db,
               ont = "MF", 
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               verbose = T)
dotplot(gse_go, showCategory=10, split=".sign", label_format = 70) + facet_grid(.~.sign)