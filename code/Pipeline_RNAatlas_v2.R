################################################################################
# 11. Deconvolution using CIBERSORTX
################################################################################
rm(list=ls())
gc()

# devtools::install_github('xuranw/MuSiC')
# install.packages("dtangle")


dir.create('results/RESTART/10.Deconvolution', showWarnings = F)
source('scripts/manuscript/Utils.R')

type = 'tissue'

# Load Linear Data
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)
cts = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
cts = cts[!(rowSums(cts[-1]) == 0),]
rownames(cts) = cts$gene
cts = cts[,-1]


cts$GeneSymbol = rownames(cts)
cts = cts[,c('GeneSymbol', meta$name[meta$organ_system == 'Central Nervous System'])]
duplicated(cts$GeneSymbol)%>%table
cts[1:4,1:4]
rownames(cts) = NULL
# fwrite(cts, 'results/RESTART/11.Deconvolution/mixture_matrix_tpm.tsv', sep = '\t', row.names = F)


# # Load Signature Matrix - VL
# require(openxlsx) 
# signatures = read.xlsx('results/RESTART/11.Deconvolution/signature_matrix.xlsx', sheet = 7, startRow = 1, colNames = T, rowNames = T)
# signatures$Neurons = NULL


# Load Signature Matrix - CA
require(openxlsx)
signatures = read.xlsx('results/RESTART/11.Deconvolution/signature_matrix.xlsx', sheet = 1, startRow = 1, colNames = T, rowNames = T)
signatures$Neurons = NULL
signatures$Endothelia = NULL
signatures$OPCs = NULL
signatures%>%colnames

# # Load Signature Matrix - LK
# require(openxlsx)
# signatures = read.xlsx('results/RESTART/11.Deconvolution/signature_matrix.xlsx', sheet = 8, startRow = 1, colNames = T, rowNames = T)
# signatures%>%colnames
# signatures$Neurons = NULL
# signatures$Endothelia = NULL


# # ID Mapping
# require(biomaRt)
# mart = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
# id_map = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = mart, filter = 'ensembl_gene_id', values = rownames(signatures))
# fwrite(id_map, 'results/RESTART/11.Deconvolution/id_map_LK.csv')

id_map = data.frame(fread('results/RESTART/11.Deconvolution/id_map_CA.csv'))
id_map = id_map[!(id_map$hgnc_symbol == ''),]

signatures = signatures[rownames(signatures) %in% id_map$ensembl_gene_id,]
signatures$GeneSymbol = sapply(rownames(signatures), function(x) id_map$hgnc_symbol[id_map$ensembl_gene_id == x])

is.na(signatures$GeneSymbol) %>% table
table(signatures$GeneSymbol == '')

rownames(signatures) = NULL
signatures = signatures[,c('GeneSymbol', colnames(signatures)[-ncol(signatures)])]
table(rowSums(signatures[-1]) == 0)

# fwrite(signatures, 'results/RESTART/11.Deconvolution/signature_matrix_VL.tsv', sep = '\t')

cts$GeneSymbol = make.names(cts$GeneSymbol, unique=TRUE)
signatures$GeneSymbol = make.names(signatures$GeneSymbol, unique=TRUE)

intersected = intersect(signatures$GeneSymbol, cts$GeneSymbol)

cts = cts[cts$GeneSymbol %in% intersected,]
signatures = signatures[signatures$GeneSymbol %in% intersected,]

fwrite(cts, 'results/RESTART/11.Deconvolution/mixture_matrix_tpm_CA_int_OutOPCs.txt', sep = '\t', row.names = F, quote = F)
fwrite(signatures, 'results/RESTART/11.Deconvolution/signature_matrix_CA_int_OutOPCs.txt', sep = '\t', row.names = F, quote = F)


# Results & Visualization
results = data.frame(fread('results/RESTART/11.Deconvolution/results/CA_deconv_perm100.csv'))

rownames(results) = results$Mixture
results$Mixture = NULL
results = results[,!(colnames(results) %in% c('P.value','Correlation','RMSE'))]


# Proportional Bar Plots
results$ID = rownames(results)
results_melted = reshape2::melt(results)

# merged %>% arrange(desc(Sum_SCI)) %>%
#   .$Row.names -> ID_order
# 
# results_melted$ID = factor(results_melted$ID,
#                            levels = rev(ID_order))

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

saveRDS(results_melted, 'results/RESTART/11.Deconvolution/results/df_to_plot_CA.rds')

ggplot(results_melted, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'fill') +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = 'Cell Type Deconvolution', x = '', y = 'Proportion')+
  ggsci::scale_fill_npg()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

ggsave('results/RESTART/11.Deconvolution/results/Prop_Bar_CA.pdf', width = 6, height = 5)

################################################################################
# 12. Cell Type Proportion <-> Total SCI, ADAR-binding SCI Sum
################################################################################
rm(list = ls())
gc()
source('scripts/manuscript/Utils.R')

# Results & Visualization
results = data.frame(fread('results/RESTART/11.Deconvolution/results/CA_deconv_perm100.csv'))

rownames(results) = results$Mixture
results$Mixture = NULL
results = results[,!(colnames(results) %in% c('P.value','Correlation','RMSE'))]

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

# Summation Expression of irCLASH targeted circRNA
list_ADARbound = readRDS('results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')

core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
core_enriched_adar3 = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
core_enriched_adar = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR_One_Union',]$core_enrichment, pattern = '/'))
list_ADARbound$ADAR3_One_Core_Enriched = core_enriched_adar3
list_ADARbound$ADAR_One_Union_Core_Enriched = core_enriched_adar
lapply(list_ADARbound, function(x){length(x)})

# cts_sci_filtered_omitNA = na.omit(cts_sci_filtered)
lapply(list_ADARbound, function(list_circ){
  df_to_add = cts_sci_filtered[cts_sci_filtered$circRNA_id %in% list_circ,][-1] %>% colSums(na.rm = T) %>% data.frame
  return(df_to_add)
}) %>% do.call(cbind, .) -> list_ADARbound_mean

# # Relative Ratio
# lapply(list_ADARbound, function(list_circ){
#   df_to_add = cts_sci_filtered[cts_sci_filtered$circRNA_id %in% list_circ,][-1] %>% colMeans(na.rm = T) %>% data.frame
#   df_to_div = cts_sci_filtered[!(cts_sci_filtered$circRNA_id %in% list_circ),][-1] %>% colMeans(na.rm = T) %>% data.frame
#   
#   df_to_add = df_to_add/df_to_div
#   
#   return(df_to_add)
# }) %>% do.call(cbind, .) -> list_ADARbound_mean


names(list_ADARbound_mean) = names(list_ADARbound)


# Summation Expression of RIP-seq targeted circRNA
cts_agg_filtered = readRDS(file = 'results/RESTART/3.DEG_Brain/cts_agg_filtered.rds')
df_gene_bind = read.csv(file = 'results/4.ADARB2_RIP/ADAR3_Binding_RIP.txt', sep = '\t')
df_gene_bind%>%head
gene_bind = df_gene_bind$Gene

# DEG results
results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))

cutoff_LFC = 1
cutoff_padj = 0.05
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


lapply(df_geneset, function(x){
  df_to_add = cts_agg_filtered[cts_agg_filtered$hgnc_symbol %in% x,][-1] %>% colSums(na.rm = T) %>% data.frame
  return(df_to_add)
}) %>% do.call(cbind, .) -> list_RIPseq_mean
names(list_RIPseq_mean) = names(df_geneset)

# # Relative Sum of SCI
# lapply(df_geneset, function(x){
#   df_to_add = cts_agg_filtered[cts_agg_filtered$hgnc_symbol %in% x,][-1] %>% colMeans(na.rm = T) %>% data.frame
#   df_to_div = cts_agg_filtered[!(cts_agg_filtered$hgnc_symbol %in% x),][-1] %>% colMeans(na.rm = T) %>% data.frame
#   
#   df_to_add = df_to_add/df_to_div
#   
#   return(df_to_add)
# }) %>% do.call(cbind, .) -> list_RIPseq_mean
# names(list_RIPseq_mean) = names(df_geneset)


# Merge Cell Type Proportion & ADAR irCLASH & RIP-seq results
merged = merge(results, SCI, by.x = 'row.names', by.y = 'row.names')
merged = merge(merged, list_ADARbound_mean, by.x = 'Row.names', by.y = 'row.names')
merged = merge(merged, list_RIPseq_mean, by.x = 'Row.names', by.y = 'row.names')

saveRDS(merged, file = 'results/RESTART/12.Cor_Brain/merged.rds')
# saveRDS(merged, file = 'results/RESTART/12.Cor_Brain/merged_relative.rds')
# ScatterPlot
# ggplot(merged, aes(x = Scaled_SCI, y = Inhibitory)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   labs(title = 'Inhibitory', x = 'Scaled_SCI', y = 'Inhibitory') +
#   theme_classic()
cor_df = lapply(colnames(results), function(celltype){
  cor = cor.test(merged$Sum_SCI, merged[,celltype], method = 'pearson')
  return(data.frame(coef = cor$estimate, p.value = cor$p.value))
}) %>% do.call(rbind,.)

cor_df_irCLASH = lapply(colnames(list_ADARbound_mean), function(SCI_){
  lapply(colnames(results), function(celltype){
    cor = cor.test(merged[[SCI_]], merged[,celltype], method = 'pearson')
    return(data.frame(coef = cor$estimate, p.value = cor$p.value))
  }) %>% do.call(rbind,.) -> cor_df_to_add
  cor_df_to_add$celltype = colnames(results)
  cor_df_to_add$SCI = SCI_
  
  return(cor_df_to_add)
}) %>% do.call(rbind,.)

cor_df_RIPseq = lapply(colnames(list_RIPseq_mean), function(SCI_){
  lapply(colnames(results), function(celltype){
    cor = cor.test(merged[[SCI_]], merged[,celltype], method = 'pearson')
    return(data.frame(coef = cor$estimate, p.value = cor$p.value))
  }) %>% do.call(rbind,.) -> cor_df_to_add
  cor_df_to_add$celltype = colnames(results)
  cor_df_to_add$SCI = SCI_
  
  return(cor_df_to_add)
}) %>% do.call(rbind,.)


rownames(cor_df) = colnames(results)

cor_df_irCLASH[cor_df_irCLASH$celltype == 'Inhibitory',] %>% arrange(p.value) %>% View()
cor_df_irCLASH[cor_df_irCLASH$p.value < 0.05,] %>% arrange(p.value) %>% View()

cor_df_RIPseq%>%arrange(p.value)%>%View
cor_df_RIPseq[cor_df_RIPseq$celltype == 'Inhibitory',] %>% arrange(p.value) %>% View()
cor_df_RIPseq[cor_df_RIPseq$p.value < 0.05,] %>% arrange(p.value) %>% View()

# write.csv(cor_df, 'results/RESTART/11.Deconvolution/results/correlation_SCI_CA.csv', row.names = T)
View(cor_df)

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
p_list

merged = readRDS(file = 'results/RESTART/12.Cor_Brain/merged.rds')
ggplot(merged[c('Row.names', 'Sum_SCI','Inhibitory','ADAR3_One_Core_Enriched','RIPseq_inhibitory')], aes(x=Inhibitory, y = ADAR3_One_Core_Enriched)) +
  geom_point(aes(fill = Sum_SCI), size = 3, shape= 21) +
  geom_text_repel(aes(label = Row.names)) +
  geom_smooth(method = 'lm', se = T, color = 'black') +
  theme_classic() +
  scale_fill_gradient(low = 'blue', high = 'red')



################################################################################
# 13. SCI Sample Correlation Heatmap
################################################################################
rm(list=ls())
gc()

dir.create('results/RESTART/12.Cor_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

require(ComplexHeatmap)

type = 'tissue'
cutoff_prop_samples = 0.2

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]


# Only Use Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']
cts_sci = cts_sci[,c('circRNA_id', list_brain)]
meta = meta[meta$name %in% list_brain,]

cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples

# Calculate Correlation
corr_matrix <- cor(cts_sci_filtered[!(names(cts_sci_filtered) %in% c('circRNA_id'))], use="complete.obs")
dim(corr_matrix)


# Draw Heatmap
df_to_heatmap = corr_matrix ; rm(corr_matrix)

# Check Orders
all(colnames(df_to_heatmap) == meta$name)

# Row Annotation
colourCount = length(unique(meta$organ_system))
getPalette = readRDS('results/RESTART/X.Figure/Robject/color_map_ggsci.rds')


# Add circRNA SCI
sci_sum = data.frame(Sample = colnames(cts_sci_filtered[-1]), Sum_SCI = colSums(cts_sci_filtered[-1], na.rm = T))

meta_merged = merge(meta, data.frame(name = sci_sum$Sample,
                                     Sum_SCI = sci_sum$Sum_SCI),
                    by = 'name')

merged = readRDS('results/RESTART/12.Cor_Brain/merged.rds')
names(merged)[1] = 'name'

meta_merged = merge(meta_merged, merged, by = 'name')
# meta_merged = meta_merged[match(colnames(df_to_heatmap), meta_merged$name),]
meta_merged = meta_merged[match(meta_merged$name, colnames(df_to_heatmap)),]

# Check Orders
all(colnames(df_to_heatmap) == meta_merged$name)


meta_subset = meta_merged[c('name','Inhibitory','Oligodendrocytes')]
rownames(meta_subset) = meta_subset$name
meta_subset$name = NULL

meta_subset = meta_subset[colnames(df_to_heatmap),]

col_ha =  columnAnnotation(
  `Total SCI` = anno_barplot(meta_merged$Sum_SCI.x, gp = gpar(fill = "#D50605E5")),
  # `Inhibitory Neuron` = anno_barplot(meta_merged$Inhibitory, gp = gpar(fill = getPalette[3])),
  # `OLG` = anno_barplot(meta_merged$Oligodendrocytes, gp = gpar(fill = "#05A18E5C")),
  `Cell Type Proportion` = anno_barplot(meta_subset,
                                        gp = gpar(fill = c("#05A18E5C", "#B1776EE5")),
                                        beside = F,
                                        width = 1),
  col = list(`Cell Type Proportion` = c(Inhibitory = "#05A18E5C", OLG = "#B1776EE5"))
)

lgd = Legend(labels = c('Inhibitory Neuron',"OLG"),
             legend_gp = gpar(fill = c("#05A18E5C","#B1776EE5")),
             title = "Cell Type")


input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)
cts = input_list[[1]] ;  meta = input_list[[2]]

cts = cts[,c('gene', list_brain)]
meta = meta[meta$name %in% list_brain,]

cts_rbp = cts[cts$gene %in% c('ADAR','ADARB1','ADARB2'),]
cts_rbp = log2(cts_rbp[-1]+1)

# scale row
cts_rbp = data.frame(t(scale(t(cts_rbp))))
colnames(cts_rbp) = colnames(cts)[-1]

# Reorder
cts_rbp = cts_rbp[list_brain]
rownames(cts_rbp) = c('ADAR','ADARB1','ADARB2')

# Heatmap
Heatmap(df_to_heatmap,
        column_title = 'Correlation Heatmap',
        top_annotation = col_ha,
        # bottom_annotation = down_ha,
        show_column_names = F,
        heatmap_legend_param = list(
          title = 'r'
        ),
        show_column_dend = FALSE,
        show_row_dend = TRUE,
        col = c("#327eba",  "white", "#e06663")
) +
  Heatmap(t(cts_rbp),
          show_column_dend = F,
          cluster_rows = F,
          cluster_columns = F,
          heatmap_legend_param = list(
            title = 'Scaled Expression'
          ),
          col = c("blue",  "white", "red")
  ) -> p_heatmap


# Heatmap with lgd
draw(p_heatmap, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right", 
     annotation_legend_list = list(lgd))


################################################################################
# 14. Correlation in CNS system : ADAR expression -> Gene Set SCI & Cell Type Proportion
################################################################################
rm(list=ls())
gc()

dir.create('results/RESTART/12.Cor_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'tissue'


# Load Cell Type Proportion and SCI by Target Genes
merged = readRDS('results/RESTART/12.Cor_Brain/merged.rds')
rownames(merged) = merged$Row.names
merged$Row.names = NULL

# Load Linear Data (Poly-A) -> TPM
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)

# Only Use Brain
cts_linear = input_list[[1]] ;  meta = input_list[[2]] ; rm(input_list)
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']
cts_linear = cts_linear[,c('gene', list_brain)]

cts_linear = cts_linear[!(rowSums(cts_linear[-1]) == 0),]
rownames(cts_linear) = cts_linear$gene
cts_linear = cts_linear[,-1]

# Only Use Brain
cts_linear = cts_linear[,list_brain]


# Only RBP lists
genelist_RBP = fread(file = 'RBP/NuclearRBP.txt', header = F) %>% unlist
cts_linear_subset = cts_linear[rownames(cts_linear) %in% genelist_RBP,]
cts_linear_subset = log2(cts_linear_subset + 1)

cts_linear_subset = cts_linear_subset[c('ADAR','ADARB1','ADARB2'),]


df_to_add = cts_linear_subset['ADARB2',]/cts_linear_subset['ADAR',]
rownames(df_to_add) = 'ADARB2-ADAR'
cts_linear_subset = rbind(cts_linear_subset, df_to_add)

df_to_add = cts_linear_subset['ADARB2',]/ cts_linear_subset['ADARB1',]
rownames(df_to_add) = 'ADARB2-ADARB1'
cts_linear_subset = rbind(cts_linear_subset, df_to_add)

df_to_add = cts_linear_subset['ADARB2',]/ (cts_linear_subset['ADARB1',] + cts_linear_subset['ADAR',])
rownames(df_to_add) = 'ADARB2-ADARB1-ADAR'
cts_linear_subset = rbind(cts_linear_subset, df_to_add)


# Correlation of ADAR, ADAR-binding circRNA SCI
input_x = as.matrix(t(cts_linear_subset)) ; input_x
input_y = as.matrix(merged) ; input_y

# Sanity Check
if (any(rownames(input_x) != rownames(input_y))){
  input_y = input_y[rownames(input_x),]
}

library(Hmisc)
corr_matrix = rcorr(input_x, input_y)
# coeffs = corr_matrix[[1]]
# pvalues = corr_matrix[[3]]

coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

coeffs_melted = reshape2::melt(coeffs)
pvalues_melted = reshape2::melt(pvalues)

# Sanity Check
if (!all(coeffs_melted$Var1 == pvalues_melted$Var1)){
  stop('Error')
}
if (!all(coeffs_melted$Var2 == pvalues_melted$Var2)){
  stop('Error')
}

res_cor = coeffs_melted
res_cor$pvalue = pvalues_melted$value
colnames(res_cor)[3] = 'coef'

res_cor[res_cor$pvalue < 0.05,] %>% View

res_cor[res_cor$Var2 == 'Inhibitory',]%>%View
res_cor[res_cor$Var1 == 'ADARB2',]%>%View



# Correlation between Inhibitory and ADAR
res_deconv = data.frame(fread('results/RESTART/11.Deconvolution/results/CA_deconv_perm100.csv.csv'))
rownames(res_deconv) = res_deconv$Mixture
res_deconv$Mixture = NULL

# Sanity Check
input_x = as.matrix(t(cts_linear_subset[c('ADAR','ADARB1','ADARB2'),])) ; input_x
input_x = apply(input_x, 2, scale); input_x
input_y = as.matrix(res_deconv[,1:(ncol(res_deconv)-3)]) ; input_y

# Calculate Correlation with RBP & Sum of SCI
require(Hmisc)
corr_matrix_2 = rcorr(x = input_x, 
                      y = input_y,
                      type = 'pearson')

coeffs = corr_matrix_2[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix_2[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

pvalues_signed = (-log10(pvalues)) * (sign(coeffs))
pheatmap::pheatmap(t(pvalues_signed), 
                   main = 'signed p-value between RBP and Cell Type',
                   color = colorRampPalette(c("blue", "white", "red"))(50)
)



################################################################################
# 15. Correlation in CNS system : ADAR expression & circRNA expression individually
################################################################################
rm(list=ls())
gc()
require(trqwe)

dir.create('results/RESTART/12.Cor_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Only Use Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']
cts_sci = cts_sci[,c('circRNA_id', list_brain)]

# Load Cell Type Proportion and SCI by Target Genes
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples
rownames(cts_sci_filtered) = cts_sci_filtered$circRNA_id
cts_sci_filtered = cts_sci_filtered[,-1]
cts_sci_filtered = log2(cts_sci_filtered + 1)


# ADAR expression
input_list = get_data_linear(slot = 'data', annot = NULL, type = type, polyA = T)
cts = input_list[[1]] ;  meta = input_list[[2]]

cts = cts[,c('gene', list_brain)]
meta = meta[meta$name %in% list_brain,]

cts_rbp = cts[cts$gene %in% c('ADAR','ADARB1','ADARB2'),]
cts_rbp = log2(cts_rbp[-1]+1)


# scale row
cts_rbp = data.frame(t(scale(t(cts_rbp))))
colnames(cts_rbp) = colnames(cts)[-1]


if (!all(rownames(merged) == colnames(cts_sci_filtered))){
  merged = merged[colnames(cts_sci_filtered),]
}

# Correlation of ADAR, ADAR-binding circRNA SCI
input_x = as.matrix(t(na.omit(cts_sci_filtered))) ; input_x
input_y = as.matrix(t(cts_rbp)) ; input_y%>%.[1:4,1:4]
dim(input_y)

# Sanity Check
if (any(rownames(input_x) != rownames(input_y))){
  input_y = input_y[rownames(input_x),]
}

# corr_matrix = rcorr(input_x, input_y)
# mcsaveRDS(corr_matrix, file = "results/RESTART/12.Cor_Brain/Cor_ADAR_circ_naomit.rds")

corr_matrix = mcreadRDS("results/RESTART/12.Cor_Brain/Cor_ADAR_circ_naomit.rds")
coeffs = corr_matrix[[1]]
pvalues = corr_matrix[[3]]

coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

coeffs_melted = reshape2::melt(coeffs)
pvalues_melted = reshape2::melt(pvalues)

# Sanity Check
if (!all(coeffs_melted$Var1 == pvalues_melted$Var1)){
  stop('Error')
}
if (!all(coeffs_melted$Var2 == pvalues_melted$Var2)){
  stop('Error')
}

res_cor = coeffs_melted
res_cor$pvalue = pvalues_melted$value
colnames(res_cor)[3] = 'coef'

# res_cor[res_cor$Var1 %in% c('ADARB2'),] %>% View
# res_cor[res_cor$pvalue < 0.05,] %>% View


# Plotting Correlation
names(res_cor) = c('RBP','circRNA_id','coef','pvalue')
res_cor$padj = p.adjust(res_cor$pvalue, method = 'BH')

res_cor$RBP[res_cor$pvalue < 0.05 & res_cor$coef > 0] %>% table
res_cor$RBP[res_cor$pvalue < 0.05 & res_cor$coef < 0] %>% table

res_cor$RBP = factor(res_cor$RBP, levels = c('ADAR','ADARB1','ADARB2'))

# ADAR3-binding circRNA Correlation with ADAR
# circRNA Annotation with  irCLASH binding
list_ADARbound = readRDS('results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')
core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
core_enriched_adar3 = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
core_enriched_adar = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR_One_Union',]$core_enrichment, pattern = '/'))
list_ADARbound$ADAR3_One_Core_Enriched = core_enriched_adar3
list_ADARbound$ADAR_One_Union_Core_Enriched = core_enriched_adar
lapply(list_ADARbound, function(x){length(x)})
df_ADARbound = reshape2::melt(list_ADARbound)
df_ADARbound$col = 1
df_ADARbound = reshape2::dcast(df_ADARbound, value ~ L1, value.var = 'col')
names(df_ADARbound)[1] = 'circRNA_id'


res_cor_merged = left_join(res_cor[res_cor$RBP == 'ADARB2',], df_ADARbound[c('circRNA_id','ADAR3_One')], by = 'circRNA_id')
res_cor_merged$ADAR3_One = ifelse(is.na(res_cor_merged$ADAR3_One), 'No', 'Yes')

# Box Plot
res_cor_merged$ADAR3_One = factor(res_cor_merged$ADAR3_One, levels = c('Yes','No'))
ggplot(res_cor_merged, aes(x = ADAR3_One, y = abs(coef), fill = ADAR3_One)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  ggsci::scale_fill_npg() +
  labs(x = 'ADAR3-binding', y = 'Correlation Coefficient (abs)') +
  theme(legend.position = 'none') +
  ggsignif::geom_signif(
    comparisons = list(c('Yes', 'No'))
  )



# Subset Correlation Matrix
# res_cor_subset = res_cor[res_cor$pvalue < 0.05,]
# res_cor_subset$significance = case_when(
#   res_cor_subset$coef > 0  ~ 1,
#   res_cor_subset$coef < 0 ~ -1
# )

res_cor_subset = res_cor[res_cor$circRNA_id %in% unique(res_cor$circRNA_id[res_cor$pvalue < 0.05]),]
res_cor_subset$significance = res_cor_subset$coef

res_cor_subset[c('RBP','circRNA_id','coef')]%>%nrow
res_cor_subset[c('RBP','circRNA_id','coef')]%>%na.omit%>%nrow

res_cor_subset[c('RBP','circRNA_id','significance')] %>%
  reshape2::dcast(RBP ~ circRNA_id, value.var = 'significance') %>%
  t() %>%
  data.frame -> df_to_heatmap


colnames(df_to_heatmap) = df_to_heatmap[1,]%>%unlist
df_to_heatmap = df_to_heatmap[-1,]

df_to_heatmap[is.na(df_to_heatmap)] = 0
annot = data.frame(circRNA_id = rownames(df_to_heatmap))


# circRNA Annotation with RIP-seq & irCLASH binding
list_ADARbound = readRDS('results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')

core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
core_enriched_adar3 = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
core_enriched_adar = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR_One_Union',]$core_enrichment, pattern = '/'))
list_ADARbound$ADAR3_One_Core_Enriched = core_enriched_adar3
list_ADARbound$ADAR_One_Union_Core_Enriched = core_enriched_adar
lapply(list_ADARbound, function(x){length(x)})
df_ADARbound = reshape2::melt(list_ADARbound)

df_ADARbound$col = 1
df_ADARbound = reshape2::dcast(df_ADARbound, value ~ L1, value.var = 'col')
names(df_ADARbound)[1] = 'circRNA_id'


annot = left_join(annot, df_ADARbound, by = 'circRNA_id')

rowannot = annot
rownames(rowannot) = rowannot$circRNA_id

rowannot = rowannot[c('ADAR1_One','ADAR2_One','ADAR3_One')]
rowannot[rowannot == 1] = 'Yes'
rowannot[is.na(rowannot)] = 'No'


row_filter = !((rowannot$ADAR1_One == 'No') & (rowannot$ADAR2_One == 'No') & (rowannot$ADAR3_One == 'No'))


heatmap_annot_row = rowAnnotation(
  
  ADAR_binding = rowannot$ADAR1_One[row_filter],
  ADARB1_binding = rowannot$ADAR2_One[row_filter],
  ADARB2_binding = rowannot$ADAR3_One[row_filter],
  
  col = list(ADAR_binding = c('Yes' = 'darkred', 'No' = 'white'),
             ADARB1_binding = c('Yes' = 'darkgreen', 'No' = 'white'),
             ADARB2_binding = c('Yes' = 'navy', 'No' = 'white'))
  
)

# Heatmap
Heatmap(df_to_heatmap%>%apply(.,2,as.numeric)%>%.[row_filter,],
        name = 'Correlation',
        cluster_rows = T,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = T,
        show_row_dend = F,
        col = colorRamp2::colorRamp2(c(-1, 0, 1), c("#327eba",  "white", "#e06663" )),
        use_raster = T,
        right_annotation = heatmap_annot_row
)

################################################################################
# 16. Correlation - Inhibitory Proportion & Individual circRNA expression (irCLASH)
################################################################################
rm(list=ls())
gc()
require(trqwe)

getwd()
setwd('/home/local/suwanyu_971004/circRNA/')

dir.create('results/RESTART/12.Cor_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Only Use Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']
cts_sci = cts_sci[,c('circRNA_id', list_brain)]

# Load Cell Type Proportion and SCI by Target Genes
cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples
rownames(cts_sci_filtered) = cts_sci_filtered$circRNA_id
cts_sci_filtered = cts_sci_filtered[,-1]
cts_sci_filtered = log2(cts_sci_filtered + 1)


# Load Cell Type Proportion and SCI by Target Genes
merged = readRDS('results/RESTART/12.Cor_Brain/merged.rds')
rownames(merged) = merged$Row.names
merged$Row.names = NULL
merged = merged[c('Astrocytes', 'Excitatory', 'Inhibitory',   'Microglia', 'Oligodendrocytes')]

if (!all(rownames(merged) == colnames(cts_sci_filtered))){
  merged = merged[colnames(cts_sci_filtered),]
}

# Correlation of ADAR, ADAR-binding circRNA SCI
input_x = as.matrix(t(na.omit(cts_sci_filtered))) ; input_x
input_y = as.matrix(merged) ; input_y%>%.[1:4,1:4]
dim(input_y)

# Sanity Check
if (any(rownames(input_x) != rownames(input_y))){
  input_y = input_y[rownames(input_x),]
}

corr_matrix = mcreadRDS(file = "results/RESTART/12.Cor_Brain/Cor_Proportion_circ_naomit_exceptFetal.rds")
coeffs = corr_matrix[[1]]
pvalues = corr_matrix[[3]]

coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

coeffs_melted = reshape2::melt(coeffs)
pvalues_melted = reshape2::melt(pvalues)

# Sanity Check
if (!all(coeffs_melted$Var1 == pvalues_melted$Var1)){
  stop('Error')
}
if (!all(coeffs_melted$Var2 == pvalues_melted$Var2)){
  stop('Error')
}

res_cor = coeffs_melted
res_cor$pvalue = pvalues_melted$value
colnames(res_cor)[3] = 'coef'

# res_cor[res_cor$Var1 %in% c('ADARB2'),] %>% View
# res_cor[res_cor$pvalue < 0.05,] %>% View


# Plotting Correlation
names(res_cor) = c('circRNA_id','CellType','coef','pvalue')
res_cor$padj = p.adjust(res_cor$pvalue, method = 'BH')

res_cor$CellType[res_cor$pvalue < 0.05 & res_cor$coef > 0] %>% table
res_cor$CellType[res_cor$pvalue < 0.05 & res_cor$coef < 0] %>% table

# ADAR3-binding circRNA Correlation with ADAR
# circRNA Annotation with  irCLASH binding
list_ADARbound = readRDS('results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')
core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
core_enriched_adar3 = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
core_enriched_adar = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR_One_Union',]$core_enrichment, pattern = '/'))
list_ADARbound$ADAR3_One_Core_Enriched = core_enriched_adar3
list_ADARbound$ADAR_One_Union_Core_Enriched = core_enriched_adar
lapply(list_ADARbound, function(x){length(x)})
df_ADARbound = reshape2::melt(list_ADARbound)
df_ADARbound$col = 1
df_ADARbound = reshape2::dcast(df_ADARbound, value ~ L1, value.var = 'col')
names(df_ADARbound)[1] = 'circRNA_id'
list_ADARbound%>%names


res_cor_merged = left_join(res_cor[res_cor$CellType == 'Inhibitory',], df_ADARbound[c('circRNA_id','ADAR_One_Union')], by = 'circRNA_id')
res_cor_merged$ADAR_One_Union = ifelse(is.na(res_cor_merged$ADAR_One_Union), 'No', 'Yes')

# Box Plot
res_cor_merged$ADAR_One_Union = factor(res_cor_merged$ADAR_One_Union, levels = c('Yes','No'))
ggplot(res_cor_merged, aes(x = ADAR_One_Union, y = abs(coef), fill = ADAR_One_Union)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggsci::scale_fill_npg() +
  labs(x = 'ADAR_One_Union', y = 'Correlation Coefficient (abs)') +
  theme(legend.position = 'none') +
  ggsignif::geom_signif(
    comparisons = list(c('Yes', 'No'))
  )
ggplot(res_cor_merged, aes(x = ADAR_One_Union, y = coef, fill = ADAR_One_Union)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggsci::scale_fill_npg() +
  labs(x = 'ADAR_One_Union', y = 'Correlation Coefficient ') +
  theme(legend.position = 'none') +
  ggsignif::geom_signif(
    comparisons = list(c('Yes', 'No'))
  )
# 
# res_cor_subset = res_cor[res_cor$circRNA_id %in% unique(res_cor$circRNA_id[res_cor$pvalue < 0.05]),]
# res_cor_subset$significance = res_cor_subset$coef
# 
# res_cor_subset[c('CellType','circRNA_id','coef')]%>%nrow
# res_cor_subset[c('CellType','circRNA_id','coef')]%>%na.omit%>%nrow
# 
# res_cor_subset[c('CellType','circRNA_id','significance')] %>%
#   reshape2::dcast(CellType ~ circRNA_id, value.var = 'significance') %>%
#   t() %>%
#   data.frame -> df_to_heatmap
# 
# 
# colnames(df_to_heatmap) = df_to_heatmap[1,]%>%unlist
# df_to_heatmap = df_to_heatmap[-1,]
# 
# df_to_heatmap[is.na(df_to_heatmap)] = 0
# annot = data.frame(circRNA_id = rownames(df_to_heatmap))
# 
# 
# # circRNA Annotation with RIP-seq & irCLASH binding
# list_ADARbound = readRDS('results/RESTART/5.irCLASH_Brain/list_ADARbound.rds')
# 
# core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
# core_enriched_adar3 = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
# core_enriched_adar = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR_One_Union',]$core_enrichment, pattern = '/'))
# list_ADARbound$ADAR3_One_Core_Enriched = core_enriched_adar3
# list_ADARbound$ADAR_One_Union_Core_Enriched = core_enriched_adar
# lapply(list_ADARbound, function(x){length(x)})
# df_ADARbound = reshape2::melt(list_ADARbound)
# 
# df_ADARbound$col = 1
# df_ADARbound = reshape2::dcast(df_ADARbound, value ~ L1, value.var = 'col')
# names(df_ADARbound)[1] = 'circRNA_id'
# 
# 
# annot = left_join(annot, df_ADARbound, by = 'circRNA_id')
# 
# rowannot = annot
# rownames(rowannot) = rowannot$circRNA_id
# 
# rowannot = rowannot[c('ADAR1_One','ADAR2_One','ADAR3_One')]
# rowannot[rowannot == 1] = 'Yes'
# rowannot[is.na(rowannot)] = 'No'
# 
# 
# row_filter = !((rowannot$ADAR1_One == 'No') & (rowannot$ADAR2_One == 'No') & (rowannot$ADAR3_One == 'No'))
# 
# 
# heatmap_annot_row = rowAnnotation(
#   
#   ADAR_binding = rowannot$ADAR1_One[row_filter],
#   ADARB1_binding = rowannot$ADAR2_One[row_filter],
#   ADARB2_binding = rowannot$ADAR3_One[row_filter],
#   
#   col = list(ADAR_binding = c('Yes' = 'darkred', 'No' = 'white'),
#              ADARB1_binding = c('Yes' = 'darkgreen', 'No' = 'white'),
#              ADARB2_binding = c('Yes' = 'navy', 'No' = 'white'))
#   
# )
# 
# # Heatmap
# Heatmap(df_to_heatmap%>%apply(.,2,as.numeric)%>%.[row_filter, c('Inhibitory','Oligodendrocytes','Astrocytes','Excitatory')],
#         name = 'Correlation',
#         cluster_rows = T,
#         cluster_columns = T,
#         show_row_names = F,
#         show_column_names = T,
#         show_row_dend = F,
#         col = colorRamp2::colorRamp2(c(-1, 0, 1), c("#327eba",  "white", "#e06663" )),
#         use_raster = T,
#         right_annotation = heatmap_annot_row
# )


################################################################################
# 17. Correlation - Inhibitory Proportion & Individual circRNA expression (RIP-seq)
################################################################################
rm(list = ls())

dir.create('results/RESTART/12.Cor_Brain', showWarnings = F)
source('scripts/manuscript/Utils.R')

cutoff_prop_samples = 0.2
type = 'tissue'

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
meta = input_list[[2]]

# Only Use Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']

# Gene-level
cts_agg_filtered = readRDS(file = 'results/RESTART/3.DEG_Brain/cts_agg_filtered.rds')
cts_agg_filtered = cts_agg_filtered[,c('hgnc_symbol',list_brain)]
cts_agg_filtered = data.frame(cts_agg_filtered)
rownames(cts_agg_filtered) = cts_agg_filtered$hgnc_symbol
cts_agg_filtered$hgnc_symbol = NULL
colnames(cts_agg_filtered) = str_replace_all(colnames(cts_agg_filtered), '\\.', ' ')

# Load Cell Type Proportion and SCI by Target Genes
merged = readRDS('results/RESTART/12.Cor_Brain/merged.rds')
rownames(merged) = merged$Row.names
merged$Row.names = NULL
merged = merged[c('Astrocytes', 'Excitatory', 'Inhibitory',   'Microglia', 'Oligodendrocytes')]

if (!all(rownames(merged) == colnames(cts_agg_filtered))){
  merged = merged[colnames(cts_agg_filtered),]
}

# Correlation of ADAR, ADAR-binding circRNA SCI
input_x = as.matrix(t(na.omit(cts_agg_filtered))) ; input_x
input_y = as.matrix(merged) ; input_y%>%.[1:4,1:4]
dim(input_x)
dim(input_y)

# Sanity Check
if (any(rownames(input_x) != rownames(input_y))){
  input_y = input_y[rownames(input_x),]
}

# corr_matrix = rcorr(input_x, input_y)
# mcsaveRDS(corr_matrix, file = "results/RESTART/12.Cor_Brain/Cor_Proportion_genelevel_naomit_exceptFetal.rds")

corr_matrix = mcreadRDS(file = "results/RESTART/12.Cor_Brain/Cor_Proportion_genelevel_naomit_exceptFetal.rds")
coeffs = corr_matrix[[1]]; coeffs[1:4,1:4]
pvalues = corr_matrix[[3]]; pvalues[1:4,1:4]

coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

coeffs_melted = reshape2::melt(coeffs)
pvalues_melted = reshape2::melt(pvalues)

# Sanity Check
if (!all(coeffs_melted$Var1 == pvalues_melted$Var1)){
  stop('Error')
}
if (!all(coeffs_melted$Var2 == pvalues_melted$Var2)){
  stop('Error')
}

res_cor = coeffs_melted
res_cor$pvalue = pvalues_melted$value
colnames(res_cor)[3] = 'coef'

# res_cor[res_cor$Var1 %in% c('ADARB2'),] %>% View
# res_cor[res_cor$pvalue < 0.05,] %>% View


# Plotting Correlation
names(res_cor) = c('hgnc_symbol','CellType','coef','pvalue')
res_cor$padj = p.adjust(res_cor$pvalue, method = 'BH')

res_cor$CellType[res_cor$pvalue < 0.05 & res_cor$coef > 0] %>% table
res_cor$CellType[res_cor$pvalue < 0.05 & res_cor$coef < 0] %>% table


# RIP-seq
df_gene_bind = read.csv(file = 'results/4.ADARB2_RIP/ADAR3_Binding_RIP.txt', sep = '\t')
gene_bind = df_gene_bind$Gene

# Universe Genesets
results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))

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
df_geneset = reshape2::melt(df_geneset)
df_geneset$col = 1
df_geneset = reshape2::dcast(df_geneset, value ~ L1, value.var = 'col')
names(df_geneset)[1] = 'hgnc_symbol'


res_cor_merged = left_join(res_cor[res_cor$CellType == 'Inhibitory',], df_geneset[c('hgnc_symbol','RIPseq_inhibitory','RIPseq_inhibitory2')], by = 'hgnc_symbol')
res_cor_merged$RIPseq_inhibitory2 = ifelse(is.na(res_cor_merged$RIPseq_inhibitory2), 'No', 'Yes')

# Box Plot
res_cor_merged$RIPseq_inhibitory2 = factor(res_cor_merged$RIPseq_inhibitory2, levels = c('Yes','No'))
ggplot(res_cor_merged, aes(x = RIPseq_inhibitory2, y = coef, fill = RIPseq_inhibitory2)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggsci::scale_fill_npg() +
  labs(x = 'RIPseq_inhibitory2', y = 'Correlation Coefficient ') +
  theme(legend.position = 'none') +
  ggsignif::geom_signif(
    comparisons = list(c('Yes', 'No'))
  )

################################################################################
# Correlation (RBP & Sum SCI)
################################################################################

# Load SCI Data
input_list = get_data(slot = 'rpm', annot = NULL, type = type, LCratio = T)
cts_sci = input_list[[1]] ;  meta = input_list[[2]]

# Only Use Brain
list_brain = meta$name[meta$organ_system == 'Central Nervous System']
list_brain = list_brain[list_brain != 'Fetal Brain']
cts_sci = cts_sci[,c('circRNA_id', list_brain)]

cts_sci_filtered = cts_sci[apply(cts_sci[-1], 1, function(x){sum(!is.na(x))/length(x) >= cutoff_prop_samples}),] # Filter out circRNAs with less than 20% of samples

# Sum of SCI
sci_sum = data.frame(Sample = colnames(cts_sci_filtered[-1]), Sum_SCI = colSums(cts_sci_filtered[-1], na.rm = T))
sci_sum$Sum_SCI = log2(sci_sum$Sum_SCI + 1)

input_x = as.matrix(t(cts_linear_subset)) ; input_x[1:4,1:4]
input_y = as.matrix(data.frame(sci_sum[-1])) ; input_y%>%head

# Calculate Correlation with RBP & Sum of SCI
require(Hmisc)
corr_matrix = rcorr(x = input_x, 
                    y = input_y,
                    type = 'pearson')

corr_matrix%>%names
coeffs = corr_matrix[[1]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]
pvalues = corr_matrix[[3]][1:ncol(input_x), (ncol(input_x)+1):(ncol(input_x)+ncol(input_y))]

res_cor = data.frame(coeffs, pvalues)
res_cor$padj = p.adjust(res_cor$pvalues, method = 'BH')
res_cor%>%head
res_cor%>%is.na%>%table
res_cor = na.omit(res_cor)
saveRDS(res_cor, file = file.path('results/RESTART/12.Cor_Brain/', str_glue('res_cor_{type}_polyA.rds')))


# Volcano Plot
res_cor = readRDS(file = file.path('results/RESTART/12.Cor_Brain/', str_glue('res_cor_{type}_polyA.rds')))
table(res_cor$padj < 0.05)
# log2FC_cutoff = min(abs(res_cor$coeffs[(res_cor$padj < 0.05)]))
# df_to_plot <- add_regulate(res_cor, log2FC_name = "coeffs",fdr_name = "padj",
#                            log2FC = log2FC_cutoff,
#                            fdr = 0.05)
log2FC_cutoff = min(abs(res_cor$coeffs[(res_cor$pvalues < 0.05)]))
# table(df_to_plot$pvalues < 0.05)

# df_to_plot <- add_regulate(res_cor, log2FC_name = "coeffs",fdr_name = "pvalues",
#                            log2FC = log2FC_cutoff,
#                            fdr = 0.05)

df_to_plot$regulate = case_when(df_to_plot$coeffs > 0 & df_to_plot$pvalues < 0.05 ~ 'Up',
                                df_to_plot$coeffs < 0 & df_to_plot$pvalues < 0.05 ~ 'Down',
                                TRUE ~ 'Normal')

# Generate the plot
df_to_plot$label = rownames(df_to_plot)
df_to_plot$label = ifelse(df_to_plot$label %in% c('ADAR', 'ADARB1','ADARB2'), df_to_plot$label, NA)
df_to_plot$label = factor(df_to_plot$label, levels = c('ADAR', 'ADARB1', 'ADARB2'))

# Define colors
color_map <- c("Up" = "#E41A1C", "Down" = "blue", "Normal" = "grey",
               "ADAR" = "#E64B3599", "ADARB1" =  "#4DBBD599" , "ADARB2" = "#00A08799"
)

# Volcano Plot
ggplot(df_to_plot, aes(x = coeffs, 
                       y = -log10(pvalues),
                       # y = -log10(p.adj)
                       label = label)) +
  geom_point(aes(color = regulate), alpha = ifelse(df_to_plot$label %in% c("ADAR", "ADARB1", "ADARB2"), 1, 0.1),
             size = ifelse(df_to_plot$label %in% c("ADAR", "ADARB1", "ADARB2"), 3, 1.5)) +
  geom_point(data = subset(df_to_plot, label %in% c("ADAR", "ADARB1", "ADARB2")),
             aes(color = label), size = 3) +
  labs(x = 'Correlation Coefficient', y = '-Log10 P-value') +
  theme_classic() +
  theme(legend.position = "none",
        #        plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12)
  ) +
  scale_color_manual(values = color_map) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.5) +
  geom_text_repel(data = subset(df_to_plot, label %in% c("ADAR", "ADARB1", "ADARB2")),
                  aes(color = label),
                  nudge_x = 0.3,
                  nudge_y = -0.1,
                  size = 5,
                  max.overlaps = 10)

# Scatter Plot
df_to_scatter = cbind(data.frame(input_x), data.frame(input_y))
df_to_scatter$Sample = rownames(df_to_scatter)
df_to_scatter = left_join(df_to_scatter, meta, by = c('Sample' = 'name'))

# Scatterplot
lapply(c('ADAR', 'ADARB1', 'ADARB2'), function(feature) {
  # Filtering data for text labels before passing to ggplot to simplify the plot code
  # filtered_data <- df_to_scatter[(df_to_scatter[[feature]] > quantile(df_to_scatter[[feature]], 0.75) & df_to_scatter$Sum_SCI > quantile(df_to_scatter$Sum_SCI,  0.75)) | (str_detect(df_to_scatter$Sample, 'Brain|brain')), ]
  
  ggplot(df_to_scatter, aes_string(x = feature, y = "Sum_SCI")) +
    geom_smooth(method = 'lm', se = T, color = 'black') +
    geom_point(alpha = 0.5, size = 2) +
    labs(x = paste(feature, "Expression"), y = 'SCI by Tissues') +
    theme_classic() +
    geom_text_repel(aes(label=Sample)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
          # text = element_text(family = "Helvetica", face = "plain", color = "black", size = 12)
    ) +
    ggtitle(feature) -> p
  
  return(p)
}) -> list_p

p_all = ggarrange(plotlist = list_p, ncol = 3)
p_all
ggsave(file.path(str_glue('results/RESTART/12.Cor_Brain//Scatter_ADAR_{type}.pdf')), p_all, width = 18, height = 6)


################################################################################
# 11. ADAR3-binding circRNA expression vs. Non-binding circRNA expression
################################################################################
rm(list=ls())
gc()

dir.create('results/RESTART/12.Cor_Brain', showWarnings = F)
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
x
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

df_to_box_filtered$delabel%>%
  table %>%
  data.frame %>%
  arrange(desc(Freq)) %>%
  head

# ggplot(df_to_box, aes(x= variable, y = value, color = variable)) +
#   geom_jitter() +
#   geom_text_repel(data = df_to_box_filtered, aes(x= variable, y = value, label = delabel), size = 2) +
#   theme_bw() +
#   scale_color_npg() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = 'none') +
#   labs(x = 'CNS system', y = 'SCI')


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
p_combined


# Gene-level
cts_agg_filtered = readRDS(file = 'results/RESTART/3.DEG_Brain/cts_agg_filtered.rds')
cts_agg_filtered = cts_agg_filtered[,c('hgnc_symbol',list_brain)]


# RIP-seq
df_gene_bind = read.csv(file = 'results/4.ADARB2_RIP/ADAR3_Binding_RIP.txt', sep = '\t')
gene_bind = df_gene_bind$Gene

# Universe Genesets
results_to_save = readRDS(file.path('results/RESTART/3.DEG_Brain/', 'DEG_SCI_Aggregated_Sum.rds'))

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

# ggplot(df_to_box, aes(x= variable, y = value, color = variable)) +
#   geom_jitter() +
#   geom_text_repel(data = df_to_box_filtered, aes(x= variable, y = value, label = delabel), size = 2) +
#   theme_bw() +
#   scale_color_npg() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = 'none') +
#   labs(x = 'CNS system', y = 'SCI')

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
