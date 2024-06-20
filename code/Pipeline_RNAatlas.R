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
# Gene Filtering or Generate NA values
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

# Save 
fwrite(LCratio, file = 'results/LCratio_withNAandPseudocount.txt', sep = '\t')

################################################################################
# 1. Calculate Scaled-Circularity-Index (Including NA)
################################################################################
LC_ratio = fread('/results/LCratio_withNAandPseudocount.txt')
LC_ratio = data.frame(LC_ratio)
LC_ratio[-ncol(LC_ratio)] = LC_ratio[-ncol(LC_ratio)] * 1e6
LC_ratio = LC_ratio[c(ncol(LC_ratio), 1:(ncol(LC_ratio)-1))]

if (!(names(LC_ratio)[1] == 'circRNA_id')){
  stop('The first column should be circRNA_id')
}

fwrite(LC_ratio, 'results/LCratio_withNAandPseudocount_Milion.txt', sep = '\t')

################################################################################
# 2. Tissue Specificity (Silhouette Score by Tissues)
################################################################################
library(cluster)
source('code/Utils.R')
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
annot = read.csv('results/annotatedBSJs.csv', row.names = 1) %>% na.omit()
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

################################################################################
# 3. DEG analysis from SCI values (Gene-level Aggregation)
################################################################################
source('code/Utils.R')

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
annotatedBSJs = data.frame(fread(file = 'results/annotatedBSJs.csv'))[-1]
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

Nuclear_RBP = data.frame(fread(file = 'results/NuclearRBP.txt', sep = '\t', header = F))
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


################################################################################
# 3. Visualization (Volcano Plot)
################################################################################
cutoff_padj = 0.05
cutoff_LFC = 1

require(ggVolcano)
require(ggVennDiagram)
require(ggplot2)

results_to_save = readRDS(file.path('results', 'DEG_SCI_Aggregated_Sum.rds'))
df_to_plot = results_to_save[!is.na(results_to_save$p_adj),]
df_to_plot <- add_regulate(df_to_plot,
                           log2FC_name = "fold_change",
                           fdr_name = "p_adj",
                           log2FC = cutoff_LFC, fdr = cutoff_padj)

# plot
ggvolcano(df_to_plot, x = "log2FoldChange", y = "padj",
          label = "hgnc_symbol", 
          label_number = 0, 
          output = FALSE)

################################################################################
# 4. SCI-Pathway Analysis
################################################################################
source('code/Utils.R')

require(DOSE)
require(clusterProfiler)

cutoff_padj = 0.05
cutoff_LFC = 1

results_to_save = readRDS(file.path('results/', 'DEG_SCI_Aggregated_Sum.rds'))

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

################################################################################
# 4. RIP-seq
################################################################################
source('code/Utils.R')

require(clusterProfiler)
require(DOSE)
require(org.Hs.eg.db)

cutoff_padj = 0.05
cutoff_LFC = 1

# Load RIP-seq results
df_gene_bind = read.csv(file = 'results/ADAR3_Binding_RIP.txt', sep = '\t')
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
df_geneset$RIPseq_inhibitory = intersect(df_geneset$BoundbyADAR2inHEKcellsinSongetal2020, df_geneset$BoundbyADAR1inBahnetal2015)
df_geneset$RIPseq_inhibitory2 = intersect(df_geneset$RIPseq_inhibitory, df_geneset$BoundbyADAR3inHEKcellsinSongetal2020)

melting_geneset = reshape2::melt(df_geneset)
colnames(melting_geneset) = c('hgnc_symbol','gs_name')
melting_geneset = melting_geneset[c(2,1)]
melting_geneset$gs_name%>%table
melting_geneset = data.frame(melting_geneset)

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

# Venn Diagram
gene_list = list(Up = up_host_gene, RIP_inhibi = df_geneset$RIPseq_inhibitory)
ggVennDiagram(gene_list,
              category.names = c('Up-DEG',"Bound by ADAR3")) +
  scale_fill_gradient(low="grey90",high = "red") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  coord_flip()

################################################################################
# 5. irCLASH
################################################################################

source('code/Utils.R')

# Annotate whether circRNA are bound by ADAR
lapply(c('ADAR','ADARB1','ADARB2'), function(x){
  annotatedBSJs_withADARfam = data.frame(fread(str_glue('results/annotatedBSJs_with{x}.tsv')))
  
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

# VennDiagram
ggVennDiagram(list_ADARbound_both,
              category.names = c('Bound by ADAR1',"Bound by ADAR2","Bound by ADAR3"),
              label_alpha = 0
              ) +
  scale_fill_gradient(low="grey90",high = "steelblue") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  scale_x_continuous(expand = expansion(mult = .2))

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

names(list_ADARbound_both) = c('ADAR1_Both','ADAR2_Both','ADAR3_Both')
names(list_ADARbound_one) = c('ADAR1_One','ADAR2_One','ADAR3_One')

# Load SCI DEG
cutoff_padj = 0.05
cutoff_LFC = 1

results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI.rds'))
results_to_save = results_to_save[!is.na(results_to_save$p_adj),]

melting_geneset = reshape2::melt(list_ADARbound)
colnames(melting_geneset) = c('hgnc_symbol','gs_name')
melting_geneset = melting_geneset[c(2,1)]
melting_geneset$gs_name%>%table
melting_geneset = data.frame(melting_geneset)


# GSEA
row_filter = rep(TRUE, nrow(results_to_save))

geneList = (results_to_save[row_filter,]$fold_change)*(-log10(results_to_save[row_filter,]$p_value))
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

################################################################################
# 8. Correlation analysis (SCI and RNA binding Protein)
################################################################################

cutoff_prop_samples = 0.2
type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples

# Sum of SCI
sci_sum = data.frame(Sample = colnames(cts_sci_filtered[-1]), Sum_SCI = colSums(cts_sci_filtered[-1], na.rm = T))
sci_sum$Sum_SCI = log2(sci_sum$Sum_SCI + 1)

# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)
cts_linear = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts_linear = cts_linear[!(rowSums(cts_linear[-1]) == 0),]
rownames(cts_linear) = cts_linear$gene
cts_linear = cts_linear[,-1]

# Only RBP lists
genelist_RBP = fread(file = 'results/NuclearRBP.txt', header = F) %>% unlist
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

coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

res_cor = data.frame(coeffs, pvalues)
res_cor$padj = p.adjust(res_cor$pvalues, method = 'BH')


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
annotatedBSJs_withAlu = data.frame(fread(file = 'results/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
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
source('code/Utils.R')

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
annotatedBSJs_withAlu = data.frame(fread(file = 'results/annotatedBSJs_withAlu_NewName.tsv',sep = '\t'))
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
