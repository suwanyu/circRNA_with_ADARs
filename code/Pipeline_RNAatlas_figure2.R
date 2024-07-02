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
# Correlation analysis (SCI and RNA binding Protein - ADAR, ADARB1, ADARB2)
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
genelist_RBP = fread(file = 'data/NuclearRBP.txt', header = F) %>% unlist
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


################################################################################
# RIP-seq
################################################################################
require(clusterProfiler)
require(DOSE)
require(org.Hs.eg.db)

cutoff_padj = 0.05
cutoff_LFC = 1

# Load RIP-seq results
df_gene_bind = read.csv(file = 'data/ADAR3_Binding_RIP.txt', sep = '\t')
df_gene_bind%>%head
gene_bind = df_gene_bind$Gene

# DEG results
results_to_save = readRDS(file.path('results/', 'DEG_SCI_Aggregated_Sum.rds'))

up_host_gene = results_to_save[results_to_save$fold_change > cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
down_host_gene = results_to_save[results_to_save$fold_change < -cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
null_host_gene = results_to_save[results_to_save$p_adj > cutoff_padj,]$hgnc_symbol %>% unique() %>% na.omit
universe_host_gene = na.omit(results_to_save$hgnc_symbol)

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
df_geneset$RIPseq_inhibitory2 = intersect(intersect(df_geneset$BoundbyADAR2inHEKcellsinSongetal2020, df_geneset$BoundbyADAR1inBahnetal2015), df_geneset$BoundbyADAR3inHEKcellsinSongetal2020)

melting_geneset = reshape2::melt(df_geneset)
colnames(melting_geneset) = c('hgnc_symbol','gs_name')
melting_geneset = melting_geneset[c(2,1)]
melting_geneset$gs_name%>%table
melting_geneset = data.frame(melting_geneset)

# GSA analysis
res_rip_up = enricher(gene = up_host_gene,
                      TERM2GENE = melting_geneset,
                      universe = universe_host_gene,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 1,
                      pAdjustMethod = 'none',
                      minGSSize = 1,
                      maxGSSize = 5000)
df_to_plot = data.frame(res_rip_up)

################################################################################
# irCLASH
################################################################################

# Annotate whether circRNA are bound by ADAR
lapply(c('ADAR','ADARB1','ADARB2'), function(x){
  annotatedBSJs_withADARfam = data.frame(fread(str_glue('data/ADARB2_irCLASH/annotatedBSJs_with{x}.tsv')))
  
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

# Only Target both flanking introns
list_ADARbound_both = lapply(list_df, function(x){
  x = x[x$circRNA_id %in% list_circRNA_id,]
  target = x$circRNA_id[x$Group == 'Target']
  return(target)
})

# Target at least one flanking intron
list_ADARbound_one = lapply(list_df, function(x){
  x = x[x$circRNA_id %in% list_circRNA_id,]
  target = x$circRNA_id[!(x$Group %in% c('Non_target','Moderate'))]
  return(target)
})

names(list_ADARbound_both) = c('ADAR1_Both','ADAR2_Both','ADAR3_Both')
names(list_ADARbound_one) = c('ADAR1_One','ADAR2_One','ADAR3_One')


# Load SCI DEG
cutoff_padj = 0.05
cutoff_LFC = 1

results_to_save = readRDS(file.path('results/', 'DEG_SCI.rds'))
results_to_save = results_to_save[!is.na(results_to_save$p_adj),]

up_circ_gene = results_to_save[results_to_save$fold_change > cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$circRNA_id %>% unique() %>% na.omit
down_circ_gene = results_to_save[results_to_save$fold_change < -cutoff_LFC & results_to_save$p_adj < cutoff_padj,]$circRNA_id %>% unique() %>% na.omit
null_circ_gene = results_to_save[results_to_save$p_adj > cutoff_padj,]$circRNA_id %>% unique() %>% na.omit
universe_circ_gene = na.omit(results_to_save$circRNA_id)

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
