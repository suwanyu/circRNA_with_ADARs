library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(patchwork)
library(SingleR)
library(scater)
library(scran)
library(dplyr)
library(SingleCellExperiment)
library(celldex)
library(BiocParallel)
library(scRNAseq)
library(data.table)
library(trqwe)

################################################################################
# Preprocess Raw Datasets
################################################################################

# Load Datasets
meta = data.frame(fread('data/sample.csv', sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('data/',
                  meta$Accession.number,
                  "_gene_count.tsv")
names(list_dir) = meta$Accession.number

list_GSE.Number = c("GSE125288", "GSE67835", "GSE71315", "GSE75140")
list_dir = list_dir[list_GSE.Number]

featuredata = data.frame(fread('data/biomart_240307.csv'))
featuredata = featuredata[c('ensembl_gene_id','hgnc_symbol')]
names(featuredata) = c('GeneID','GeneSymbol')
featuredata = featuredata[featuredata$GeneSymbol!='',]
featuredata = featuredata[!duplicated(featuredata),]

# Convert Gene Name 
dir.create('cts.linear', showWarnings = F)

for (i in 1:length(list_dir)) {
  cts = data.frame(fread(list_dir[i]))
  cts = left_join(cts, featuredata, by = 'GeneID')
  cts = cts[order(rowSums(cts[sapply(cts, is.numeric)]),decreasing = T),]
  cts$GeneID = NULL
  cts = cts[!duplicated(cts$GeneSymbol),]
  
  if (any(duplicated(cts$GeneSymbol))){
    stop()
  }
  cts = cts[!is.na(cts$GeneSymbol),]
  
  cts = cbind(data.frame(GeneSymbol = cts[,ncol(cts)]),cts[,-ncol(cts)])
  fwrite(cts, file = file.path('cts.linear', paste0(names(list_dir)[i], '.csv')), row.names = F)
}

################################################################################
# Integration by 4 Methods
################################################################################

# Load Datasets
meta = data.frame(fread('data/sample.csv', sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('data/',
                  meta$Accession.number,
                  "_gene_count.tsv")
names(list_dir) = meta$Accession.number

list_GSE.Number = c("GSE125288","GSE67835", "GSE71315", "GSE75140")
list_dir = list_dir[list_GSE.Number]

dir_list = file.path('cts.linear', paste0(list_GSE.Number, '.csv'))
names(dir_list) = list_GSE.Number

# Only Intersected Gene Lists
for(i in 1:length(dir_list)){
  # Load Data
  counts <- data.frame(fread(dir_list[i]))
  rownames(counts) = counts$GeneSymbol
  counts$GeneSymbol = NULL
  
  counts = counts[!(rowSums(counts) == 0),]
  
  if (i == 1){
    gene_list = rownames(counts)
  } else {
    gene_list = intersect(gene_list, rownames(counts))
  }
}  

featuredata = data.frame(fread('data/biomaRt.csv'))
gene_list = unique(featuredata[featuredata$gene_biotype == 'protein_coding',]$hgnc_symbol)


# Some of the code below is from https://github.com/bioinfo-biols/Code_for_circSC
# Filtering Cutoffs
minGene=500
maxGene=20000

pbmclist = list()
summary = list()

for(i in 1:length(dir_list)){
  # Load Data
  counts <- data.frame(fread(dir_list[i]))
  rownames(counts) = counts$GeneSymbol
  counts$GeneSymbol = NULL
  counts = counts[!(rowSums(counts) == 0),]
  
  # Load Objects
  pbmclist[[i]] <- CreateSeuratObject(as.sparse(counts), min.cells=1, project  = names(dir_list[i]))
  
  # Remove Outlier Cells
  pbmclist[[i]][["percent.mt"]] <- PercentageFeatureSet(pbmclist[[i]], pattern="^MT-")
  pbmclist[[i]][['percent.ribo']] <- PercentageFeatureSet(pbmclist[[i]], pattern = "^RPS|^RPL")
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(pbmclist[[i]]@assays$RNA))
  HB.genes <- rownames(pbmclist[[i]]@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  pbmclist[[i]][["percent.HB"]]<-PercentageFeatureSet(pbmclist[[i]], features=HB.genes)
  
  p1 = VlnPlot(pbmclist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.HB"))
  pbmclist[[i]] <- subset(pbmclist[[i]], subset=nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < 15)

  # Only Intersected Gene Lists (Protein-coding genes)
  pbmclist[[i]] = pbmclist[[i]][gene_list,]
  
  summary$GEO[i] = names(dir_list)[i]
  summary$nCell_beforeQC[i] = dim(counts)[2]
  summary$nCell_afterQC[i] = dim(pbmclist[[i]])[2]
  summary$nFeature_RNA[i] = mean(pbmclist[[i]]@meta.data$nFeature_RNA)
  summary$percent.mt[i] = mean(pbmclist[[i]]@meta.data$percent.mt, na.rm = T)
  summary$percent.ribo[i] = mean(pbmclist[[i]]@meta.data$percent.ribo, na.rm = T)
  summary$percent.HB[i] = mean(pbmclist[[i]]@meta.data$percent.HB, na.rm = T)
}

summary = data.frame(summary)


# Merge raw samples
obj <- merge(x = pbmclist[[1]],
             y = pbmclist[2:length(pbmclist)],
             merge.data = TRUE)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Integration
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE, 
  k.weight = 90
)

obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)

saveRDS(obj, file = 'data/brain_reintegration.rds')

################################################################################
# Annotation - scType Scoring by Cell
################################################################################

# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R") # load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R") # load cell type annotation function

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets from scType
gs_list = gene_sets_prepare(db_, tissue)
gs_list$gs_positive%>%names

# list for cell types
list_target = c("Astrocytes","Cholinergic neurons","Dopaminergic neurons","Endothelial cells",
                "GABAergic neurons","Glutamatergic neurons","Immature neurons","Immune system cells","Mature neurons",
                "Microglial cells","Myelinating Schwann cells","Neural Progenitor cells",
                "Neural stem cells","Neuroblasts","Neuroepithelial cells","Non myelinating Schwann cells",
                "Oligodendrocyte precursor cells","Oligodendrocytes","Radial glial cells","Schwann precursor cells",
                "Serotonergic neurons","Tanycytes"              
)

list_celltype = names(gs_list$gs_positive)
gs_list$gs_positive = gs_list$gs_positive[(names(gs_list$gs_positive)) %in% list_target]
gs_list$gs_negative = gs_list$gs_negative[(names(gs_list$gs_negative)) %in% list_target]
saveRDS(gs_list, file = 'results/scType/gs_list.rds')

# get cell-type by cell matrix
DefaultAssay(obj) = 'RNA'
obj = ScaleData(obj, features = rownames(obj))

es.max = sctype_score(obj[["RNA"]]['scale.data'], scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

write.csv(es.max, file = 'results/scType/es_max.csv')


################################################################################
# Clustering
################################################################################

reduction_method = 'harmony'
res = 0.14
npc = 14

obj = readRDS(file = 'data/brain_reintegration.rds')

?FindNeighbors()

obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:npc)
obj <- FindClusters(obj, resolution = res, cluster.name = str_glue("harmony_clusters_{res}"))
obj <- RunUMAP(obj, reduction = "integrated.harmony", dims = 1:npc, reduction.name = "umap.harmony")

DimPlot(obj, reduction = "umap.harmony", group.by = str_glue("harmony_clusters_{res}"))

es.max = read.csv(file = 'results/scType/es_max.csv', row.names = 1)
obj = AddMetaData(obj, metadata = t(es.max))
obj$seurat_clusters = obj@meta.data[[str_glue('harmony_clusters_{res}')]]

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

obj@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]
  obj@meta.data$scType[obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Annotation By Cluster
obj$seurat_clusters = obj@meta.data[[str_glue('{reduction_method}_clusters_{res}')]]
saveRDS(obj, file = 'results/clustered.rds')

obj$scType = factor(obj$scType, levels = c("GABAergic neurons", "Oligodendrocytes", "Oligodendrocyte precursor cells",
                                           "Immature neurons", "Neural Progenitor cells","Radial glial cells", "Neuroblasts", 
                                           "Endothelial cells", "Microglial cells"
))

DimPlot(obj,
        group.by = 'scType',
        pt.size = 1,
        repel = T,
        reduction = str_glue('umap.harmony')) +
  ggtitle('Cell Type') 

obj$CellType = obj$scType
levels(obj$CellType) = c('GABA','OLG','OPC','IMN','NPC','RGC','NB','EC','MG')
saveRDS(obj, 'data/brain_annotated.rds')

# Expression By Cluster
obj = readRDS(file = str_glue('data/brain_annotated.rds'))

DefaultAssay(obj)= 'RNA'
obj = JoinLayers(obj)

list_markers = readRDS('results/scType//gs_list.rds')$gs_positive


list_markers = list_markers[levels(obj$scType)]
names(list_markers) = levels(obj$CellType)
list_markers_order = as.character(unlist(list_markers))%>%unique

order = sort(rowSums(obj@assays$RNA), decreasing = T)

obj_agg = AverageExpression(obj, features = list_markers_order, group.by = 'CellType', assays = 'RNA')$RNA
obj_agg = data.frame(obj_agg)
obj_sd = sort(apply(obj_agg, 1, sd), decreasing = T)

list_markers_sd = lapply(list_markers, function(x){
  x = x[x %in% list_markers_order]
  x_sd = sort(obj_sd[x], decreasing = T)
  x_subset = names(x_sd)[x_sd > 0.4]
  return(head(x_subset, 5))
})

list_duplicated = unlist(list_markers_sd)[duplicated(unlist(list_markers_sd))]
list_markers_sd = lapply(list_markers_sd, function(x){
  x = x[!x %in% list_duplicated]
  return(x)
})

DotPlot(obj,
        features = list_markers_sd, 
        group.by = 'CellType',
        assay = 'RNA') +
  Seurat::RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = 'white', high = "red")

################################################################################
# Add Modality (CIRC) to object
################################################################################
reduction_method = 'harmony'
res = 0.14

# Load the circRNA Annotation Data
annot = data.frame(fread(file = 'data/annotatedBSJs_withAlu_NewName_CIRC.tsv',
                         sep = '\t'))

# Load circRNA datasets
meta = data.frame(fread('data/sample.csv', 
                        sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('data/',
                  meta$Accession.number, 
                  "_bsj_circRNAs.csv")
names(list_dir) = meta$Accession.number
list_dir = list_dir[c("GSE125288", "GSE67835", "GSE71315", "GSE75140")]

# Check circRNA annotation id
list_rownames = lapply(list_dir, function(x){
  return(data.frame(fread(x))['V1'])
})
unlist(list_rownames)%>%unique%>%length

# Check Rownames
unique(unlist(list_rownames)) %in% 
  unlist(lapply(str_split(annot$id, pattern = ':'), function(str_list){
    paste0(str_list[4], ':', str_list[5], '|', str_list[6])
  })) %>% table


# Add circRNA_id to Annotation Data
circRNA_ids = unique(unlist(list_rownames))
annot_ids = unlist(lapply(str_split(annot$id, pattern = ':'), function(str_list){
  paste0(str_list[4], ':', str_list[5], '|', str_list[6])
}))

annot$circRNA_id_SC = annot_ids
annot$is_SC = annot$circRNA_id_SC %in% circRNA_ids

saveRDS(annot, file = 'data/annot.rds')

################################################################################
# Add circRNA to Seurat Object
################################################################################
obj = readRDS(file = str_glue('data/brain_annotated.rds'))
annot = readRDS('data/annot.rds')

# Load list dir
meta = data.frame(fread('data/sample.csv', 
                        sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('data/',
                  meta$Accession.number, 
                  "_bsj_circRNAs.csv")
names(list_dir) = meta$Accession.number
list_dir = list_dir[c("GSE125288", "GSE67835", "GSE71315", "GSE75140")]

# Load circRNA Count Matrix
list_cts = lapply(list_dir, function(x){
  df = data.frame(fread(x))
  rownames(df) = df$V1
  df$V1 = NULL
  return(df)
})

# Merge list of Count Matrix
for (i in 1:length(list_cts)){
  if (i == 1){
    cts = list_cts[[i]]
  } else {
    cts = merge(cts, list_cts[[i]], by = 'row.names', all = T)
    rownames(cts) = cts$Row.names
    cts$Row.names = NULL
  }
}
cts[is.na(cts)] = 0


# Filter circRNA matrix
cts = cts[,colnames(cts) %in% colnames(obj)]
cts = cts[rowSums(cts) != 0,]
cts = cts[,colSums(cts) != 0]

# Filter Seurat Object
obj$is_circ = colnames(obj) %in% colnames(cts)
obj$is_circ%>%table
obj_subset = subset(obj, subset = is_circ)


# Add Modality to obj
circ_assay = CreateAssay5Object(counts = as.sparse(cts))
obj_subset[['CIRC']] = circ_assay
saveRDS(obj_subset, file = str_glue('data/brain_annotated.rds'))
