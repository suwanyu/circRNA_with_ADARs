rm(list=ls())
suppressMessages({
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
})
dir.create('/home/local/suwanyu_971004/circRNA/results/RESTART/10.SingleCell/', showWarnings = F)
setwd('/home/local/suwanyu_971004/circRNA/results/RESTART/10.SingleCell/')

################################################################################
# Preprocess Raw Datasets
################################################################################

# Load Datasets
meta = data.frame(fread('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/sample.csv', sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/',
                  meta$Accession.number,
                  "_gene_count.tsv")
names(list_dir) = meta$Accession.number

# list_dir = list_dir[c("GSE125288", "GSE67835", "GSE71315", "GSE75140", "SRP041736")]
list_GSE.Number = c("GSE125288", "GSE67835", "GSE71315", "GSE75140")
list_dir = list_dir[list_GSE.Number] #SRP041736


featuredata = data.frame(fread('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/biomart_240307.csv'))
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
rm(list=ls())
gc()

dir.create('1.Integration_v2')

# Load Datasets
meta = data.frame(fread('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/sample.csv', sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/',
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

# # Select Only Features (Protein-coding genes)
# require(biomaRt)
# 
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# featuredata = getBM(attributes = c('hgnc_symbol', 'gene_biotype', "start_position", "end_position"),
#                     values = gene_list,
#                     filters = 'hgnc_symbol',
#                     mart = ensembl)
# 
# featuredata$gene_biotype%>%table
# featuredata%>%is.na%>%table
# featuredata[featuredata$gene_biotype == 'protein_coding',]$hgnc_symbol%>%length
# fwrite(featuredata, file = '1.Integration_v2/biomaRt.csv', row.names = F)

featuredata = data.frame(fread('1.Integration_v2/biomaRt.csv'))
gene_list = unique(featuredata[featuredata$gene_biotype == 'protein_coding',]$hgnc_symbol)


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
  
  dir.create('1.Integration_v2/QC', showWarnings = F)
  p1 = VlnPlot(pbmclist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.HB"))
  ggsave(p1, file = file.path('1.Integration_v2/QC', paste0(names(dir_list)[i], '_beforeQC.pdf')), width = 8, height = 8)

  pbmclist[[i]] <- subset(pbmclist[[i]], subset=nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < 15)
  p1 = VlnPlot(pbmclist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.HB"))
  ggsave(p1, file = file.path('1.Integration_v2/QC', paste0(names(dir_list)[i], '_afterQC.pdf')), width = 8, height = 8)
  
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
fwrite(summary, file = '1.Integration_v2/QC/summary.csv', row.names = F)



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

saveRDS(obj, file = 'Robject_v2/brain_reintegration.rds')

################################################################################
# Clustering Analysis
################################################################################
set.seed(123)
dir.create('2.Clustering_v2', showWarnings = F)

obj = readRDS(file = 'Robject_v2/brain_reintegration.rds')

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:10)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:10)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:10)
obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:10)

for (resolution in c(0.1,0.15,0.2,0.3,0.5,0.6,0.8)){
  obj <- FindClusters(obj, resolution = resolution, cluster.name = str_glue("cca_clusters_{resolution}"))
  obj <- FindClusters(obj, resolution = resolution, cluster.name = str_glue("mnn_clusters_{resolution}"))
  obj <- FindClusters(obj, resolution = resolution, cluster.name = str_glue("rpca_clusters_{resolution}"))
  obj <- FindClusters(obj, resolution = resolution, cluster.name = str_glue("harmony_clusters_{resolution}"))
}


obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:10, reduction.name = "umap.cca")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:10, reduction.name = "umap.rpca")
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:10, reduction.name = "umap.mnn")
obj <- RunUMAP(obj, reduction = "integrated.harmony", dims = 1:10, reduction.name = "umap.harmony")

dir.create(path = '2.Clustering_v2/cca/', showWarnings = FALSE)
dir.create(path = '2.Clustering_v2/rpca/', showWarnings = FALSE)
dir.create(path = '2.Clustering_v2/mnn/', showWarnings = FALSE)
dir.create(path = '2.Clustering_v2/harmony/', showWarnings = FALSE)

for (resolution in c(0.1,0.15,0.2,0.3,0.5,0.6,0.8)){
  DimPlot(obj, reduction = "umap.cca", group.by = str_glue("cca_clusters_{resolution}"))
  ggsave(str_glue('2.Clustering_v2/cca//cca_clusters_{resolution}.pdf'), width = 8, height = 8)
  
  DimPlot(obj, reduction = "umap.rpca", group.by = str_glue("rpca_clusters_{resolution}"))
  ggsave(str_glue('2.Clustering_v2/rpca//rpca_clusters_{resolution}.pdf'), width = 8, height = 8)
  
  DimPlot(obj, reduction = "umap.mnn", group.by = str_glue("mnn_clusters_{resolution}"))
  ggsave(str_glue('2.Clustering_v2/mnn//mnn_clusters_{resolution}.pdf'), width = 8, height = 8)
  
  DimPlot(obj, reduction = "umap.harmony", group.by = str_glue("harmony_clusters_{resolution}"))
  ggsave(str_glue('2.Clustering_v2/harmony//harmony_clusters_{resolution}.pdf'), width = 8, height = 8)
}

saveRDS(obj, file = 'Robject_v2/brain_clustered.rds')


################################################################################
# Annotation - scType Scoring by Cell
################################################################################
rm(list=ls())
gc()

dir.create('3.Annotation_v2/scType', showWarnings = F, recursive = T)

obj = readRDS(file = 'Robject_v2/brain_clustered.rds')

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
saveRDS(gs_list, file = '3.Annotation_v2/scType/gs_list.rds')

# get cell-type by cell matrix
DefaultAssay(obj) = 'RNA'
obj = ScaleData(obj, features = rownames(obj))

es.max = sctype_score(obj[["RNA"]]['scale.data'], scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

write.csv(es.max, file = '3.Annotation_v2/scType/es_max.csv')

################################################################################
# Annotation by Cluster
################################################################################
rm(list=ls())
gc()

# Cell Type Annotation
obj = readRDS(file = 'Robject_v2/brain_clustered.rds')
es.max = read.csv(file = '3.Annotation_v2/scType/es_max.csv', row.names = 1)
obj = AddMetaData(obj, metadata = t(es.max))


for (reduction_method in c('cca', 'rpca', 'mnn','harmony')){
  
  dir.create(path = str_glue('3.Annotation_v2/{reduction_method}'), showWarnings = FALSE)
  
  for (res in c(0.1, 0.15, 0.2,0.3,0.5,0.6,0.8)){

    # Annotation By Cluster
    obj$seurat_clusters = obj@meta.data[[str_glue('{reduction_method}_clusters_{res}')]]
    
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
    
    DimPlot(obj, group.by = 'scType', label=T, pt.size = 1.5, reduction = str_glue('umap.{reduction_method}')) +
      ggtitle(str_glue('Resolution: {res}'))
    ggsave(str_glue('3.Annotation_v2/{reduction_method}/scType_{reduction_method}_{res}.pdf'), width = 9, height = 6)
  }
}


################################################################################
# Check Cell Type Annotation using Previous Article (Cell Type Markers)
################################################################################
rm(list=ls())
gc()

reduction_method = 'harmony'
res = 0.3

dir.create(str_glue('3.Annotation_v2/{reduction_method}_{res}'), showWarnings = F, recursive = T)

obj = readRDS(file = 'Robject_v2/brain_clustered.rds')
es.max = read.csv(file = '3.Annotation_v2/scType/es_max.csv', row.names = 1)
obj = AddMetaData(obj, metadata = t(es.max))

obj@meta.data%>%colnames

FeaturePlot(obj, 
            features = c("GABAergic neurons", "Glutamatergic neurons", "Oligodendrocyte precursor cells", "Oligodendrocytes",
                         "Immature neurons", "Neural Progenitor cells", "Radial glial cells", "Neuroblasts", "Microglial cells",
                         "Astrocytes", "Endothelial cells"
                         ), 
            pt.size = 0.1,
            reduction = str_glue('umap.{reduction_method}'))
ggsave(str_glue('3.Annotation_v2/{reduction_method}_{res}/FeaturePlot_scType.pdf'), width = 15, height = 14)


# Annotation By Cluster
obj$seurat_clusters = obj@meta.data[[str_glue('{reduction_method}_clusters_{res}')]]
DimPlot(obj, 
        group.by = 'seurat_clusters',
        label=T, 
        pt.size = 1,
        reduction = str_glue('umap.{reduction_method}'))
ggsave(str_glue('3.Annotation_v2/{reduction_method}_{res}/ClusterPlot.pdf'), width = 7, height = 6)


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


obj$scType%>%table%>%sort(., decreasing = T)%>%print()

# obj$scType = factor(obj$scType, levels = c("GABAergic neurons", "Glutamatergic neurons", "Oligodendrocyte precursor cells", "Oligodendrocytes",
#                                            "Immature neurons", "Neural Progenitor cells", "Radial glial cells", "Neuroblasts", "Microglial cells",
#                                            "Astrocytes", "Endothelial cells", "Mature neurons"
#                                            )
#                     )

obj$scType = factor(obj$scType, levels = c("GABAergic neurons", "Oligodendrocytes",
                                           "Immature neurons", "Neural Progenitor cells", "Radial glial cells", "Neuroblasts",
                                           "Astrocytes", "Endothelial cells", "Mature neurons"
)
)


DimPlot(obj,
        label = TRUE, 
        group.by = 'scType',
        pt.size = 1,
        repel = T,
        reduction = str_glue('umap.{reduction_method}')) +
  ggtitle('Cell Type')
ggsave(str_glue('3.Annotation_v2/{reduction_method}_{res}/CellTypePlot.pdf'), width = 7, height = 6)

saveRDS(obj, file = str_glue('Robject_v2/brain_annotated.rds'))

obj$CellType = obj$scType
# levels(obj$CellType) = c('GABA','GLUT','OPC','OLG','IMN','NPC','RGC','NB','MG','ASC','EC')
levels(obj$CellType) = c('GABA','OLG','IMN','NPC','RGC','NB','ASC','EC','MN')


DimPlot(obj,
        label = TRUE, 
        group.by = 'CellType',
        pt.size = 1,
        repel = T,
        reduction = str_glue('umap.{reduction_method}')) +
  ggtitle('Cell Type')
ggsave(str_glue('3.Annotation_v2/{reduction_method}_{res}/CellTypePlot_Abbre.pdf'), width = 7, height = 6)

# Expression By Cluster
DefaultAssay(obj)= 'RNA'
obj = JoinLayers(obj)

list_markers = readRDS('3.Annotation_v2/scType//gs_list.rds')$gs_positive

lapply(list_markers, function(x){sum(!(x %in% rownames(obj)))}) %>% data.frame
lapply(list_markers, function(x){length(x)}) %>% data.frame

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
  x_subset = names(x_sd)[x_sd > 0.5]

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

ggsave(str_glue('3.Annotation_v2/{reduction_method}_{res}/DotPlot.pdf'), width = 10, height = 8)
Idents(obj) = obj$CellType
saveRDS(obj, file = str_glue('Robject_v2/brain_annotated_v2.rds'))

################################################################################
# Feature Plot ADAR Expression
################################################################################

dir.create('4.Result_v2', showWarnings = F)
DefaultAssay(obj) = 'RNA'
FeaturePlot(obj, 
            features = c('ADAR','ADARB1','ADARB2'), 
            ncol = 3,
            order = T,
            reduction = str_glue('umap.{reduction_method}'))
ggsave('4.Result_v2/ADAR_FeaturePlot.pdf', width = 14, height = 4)

VlnPlot(obj, features = c('ADAR','ADARB1','ADARB2'), 
        ncol = 3, 
        group.by = 'CellType',
        pt.size = 0)
ggsave('4.Result_v2/ADAR_VlnPlot.pdf', width = 14, height = 4)



################################################################################
# Overlap Features of RNAatlas & single cell
################################################################################
rm(list=ls())
gc()

reduction_method = 'harmony'
res = 0.3

# Load the circRNA Annotation Data
annot = data.frame(fread(file = '/home/local/suwanyu_971004/circRNA/results/circRNAprofiler/annotatedBSJs_withAlu_NewName.tsv',
                         sep = '\t'))

# Load circRNA datasets
meta = data.frame(fread('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/sample.csv', 
                        sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/',
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
# annot = annot[annot$circRNA_id_SC %in% circRNA_ids,]
annot$is_SC = annot$circRNA_id_SC %in% circRNA_ids


################################################################################
# Add Modality (CIRC) to object
################################################################################
rm(list=ls())
gc()

reduction_method = 'harmony'
res = 0.3

# Load the circRNA Annotation Data
annot = data.frame(fread(file = '/home/local/suwanyu_971004/circRNA/results/circRNAprofiler/annotatedBSJs_withAlu_NewName.tsv',
                         sep = '\t'))

# Load circRNA datasets
meta = data.frame(fread('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/sample.csv', 
                        sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/',
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
# annot = annot[annot$circRNA_id_SC %in% circRNA_ids,]
annot$is_SC = annot$circRNA_id_SC %in% circRNA_ids

saveRDS(annot, file = 'Robject_v2/annot.rds')

################################################################################
# Add circRNA to Seurat Object
################################################################################
rm(list=ls())
gc()

obj = readRDS(file = str_glue('Robject_v2/brain_annotated_v2.rds'))
annot = readRDS('Robject_v2/annot.rds')

# Load list dir
meta = data.frame(fread('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/sample.csv', 
                        sep = '\t'))
meta = meta[meta$Tissue.Type %in% c('Brain','neural stem cell'),]

list_dir = paste0('/home/local/suwanyu_971004/circRNA/results/6.SingleCell/data/',
                  meta$Accession.number, 
                  "_bsj_circRNAs.csv")
names(list_dir) = meta$Accession.number
list_dir = list_dir[c("GSE125288", "GSE67835", "GSE71315", "GSE75140")]

# Load circRNA Count Matrix
list_cts = lapply(list_dir, function(x){
  df = data.frame(fread(x))
  # df = df[df$V1 %in% annot$circRNA_id_SC,]
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
saveRDS(obj_subset, file = str_glue('Robject_v2/brain_annotated_v2_circ.rds'))


# Venn Diagram Single cell And circAtlas
require(ggVennDiagram)
ggVennDiagram(list(circSC = rownames(cts), 
                   RNAatlas = unlist(lapply(str_split(annot$id, pattern = ':'), function(str_list){
                     paste0(str_list[4], ':', str_list[5], '|', str_list[6])
                   }))), 
              list_names = c('circSC', 'RNAatlas'))+
  scale_fill_gradient(low="grey90",high = "steelblue") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  coord_flip()


################################################################################
# Generate circRNA FeatureData -> SingleCell_GenomicAnnotation.R
################################################################################
# rm(list=ls())
# gc()
# 
# # Annotate whether circRNA are bound by ADAR
# lapply(c('ADAR','ADARB1','ADARB2'), function(x){
#   annotatedBSJs_withADARfam = data.frame(fread(str_glue('Robject_v2/annotatedBSJs_with{x}.tsv')))
#   
#   annotatedBSJs_withADARfam$Group = apply(annotatedBSJs_withADARfam, MARGIN = 1, function(row_){
#     #row_ = annotatedBSJs_withADARfam[4,]
#     upIntronlist = row_[['Substrate.id_UpIntron']]
#     downIntronlist = row_[['Substrate.id_DownIntron']]
#     
#     upIntronlist = strsplit(upIntronlist, ',')[[1]]
#     downIntronlist = strsplit(downIntronlist, ',')[[1]]
#     
#     if (length(intersect(upIntronlist, downIntronlist)) != 0){
#       return('Target')
#     } else if (any(duplicated(upIntronlist))){
#       return('UpIntron')
#     } else if (any(duplicated(downIntronlist))){
#       return('DownIntron')
#     } else if ((length(upIntronlist) >= 1) | (length(downIntronlist) >= 1)){
#       return('Moderate') 
#     } else {
#       return('Non_target')
#     }
#   }) %>% unlist
#   annotatedBSJs_withADARfam$circRNA_id = lapply(str_split(annotatedBSJs_withADARfam$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
#   return(annotatedBSJs_withADARfam)
# }) -> list_df
# 
# 
# # Only Target both flanking introns
# list_ADARbound_both = lapply(list_df, function(x){
#   # x = x[x$circRNA_id %in% list_circRNA_id,]
#   # target = x$circRNA_id[!(x$Group %in% c('Non_target','Moderate'))]
#   target = x$circRNA_id[x$Group == 'Target']
#   return(target)
# })
# 
# # ggsci::pal_npg()(9)%>%show_col()
# # ggsci::pal_jama()(9)%>%show_col()
# # ggsci::pal_jco()(9)%>%show_col()
# 
# # VennDiagram
# ggVennDiagram(list_ADARbound_both,
#               category.names = c('Bound by ADAR1',"Bound by ADAR2","Bound by ADAR3"),
#               label_alpha = 0
# ) +
#   scale_fill_gradient(low="grey90",high = "steelblue") +
#   scale_colour_manual(values = c("black", "black", "black","black")) +
#   scale_x_continuous(expand = expansion(mult = .2))
# 
# 
# # Target at least one flanking intron
# list_ADARbound_one = lapply(list_df, function(x){
#   # x = x[x$circRNA_id %in% list_circRNA_id,]
#   target = x$circRNA_id[!(x$Group %in% c('Non_target','Moderate'))]
#   # target = x$circRNA_id[x$Group == 'Target']
#   return(target)
# })
# 
# # VennDiagram
# ggVennDiagram(list_ADARbound_one,
#               category.names = c('Bound by ADAR1',"Bound by ADAR2","Bound by ADAR3"),
#               label_alpha = 0
# ) +
#   scale_fill_gradient(low="grey90",high = "steelblue") +
#   scale_colour_manual(values = c("black", "black", "black","black")) +
#   scale_x_continuous(expand = expansion(mult = .2))
# 
# names(list_ADARbound_both) = c('ADAR1_Both','ADAR2_Both','ADAR3_Both')
# names(list_ADARbound_one) = c('ADAR1_One','ADAR2_One','ADAR3_One')
# 
# list_ADARbound = c(list_ADARbound_both, list_ADARbound_one)
# list_ADARbound$ADAR_Both_Int = intersect(intersect(list_ADARbound$ADAR1_Both, list_ADARbound$ADAR2_Both), list_ADARbound$ADAR3_Both)
# list_ADARbound$ADAR_Both_Union = union(union(list_ADARbound$ADAR1_Both, list_ADARbound$ADAR2_Both), list_ADARbound$ADAR3_Both)
# list_ADARbound$ADAR_One_Int = intersect(intersect(list_ADARbound$ADAR1_One, list_ADARbound$ADAR2_One), list_ADARbound$ADAR3_One)
# list_ADARbound$ADAR_One_Union = union(union(list_ADARbound$ADAR1_One, list_ADARbound$ADAR2_One), list_ADARbound$ADAR3_One)
# list_ADARbound$ADAR_One_2and3 = intersect(list_ADARbound$ADAR2_One, list_ADARbound$ADAR3_One)
# list_ADARbound$ADAR_One_1and2 = intersect(list_ADARbound$ADAR1_One, list_ADARbound$ADAR2_One)
# saveRDS(list_ADARbound, file = 'Robject_v2/list_ADARbound.rds')
# 
# rownames(obj_subset@assays$CIRC)
# list_ADARbound = lapply(list_ADARbound, function(x){
#   unlist(lapply(str_split(x, pattern = "_"), function(y){
#     paste0(y[1], ':',y[2], '-',y[3])
#   }))
# })
# list_ADARbound
# saveRDS(list_ADARbound, file = 'Robject_v2/list_ADARbound.v2.rds')


################################################################################
# Visualize circRNA Expression
################################################################################
rm(list = ls())
gc()


reduction_method = 'harmony'
obj_subset = readRDS(file = str_glue('Robject_v2/brain_annotated_v2_circ.rds'))
DimPlot(obj_subset, group.by = 'CellType', reduction = str_glue('umap.{reduction_method}'))


# Color by celltype
DefaultAssay(obj_subset) = 'CIRC'
obj_subset
list_celltype = levels(obj_subset$CellType)
list_cols = scales::hue_pal()(length(list_celltype))
names(list_cols) = list_celltype

# Visualize nCount_CIRC, nFeature_CIRC
FeaturePlot(obj_subset,
            features = c('nCount_CIRC', 'nFeature_CIRC'),
            ncol = 2,
            order = T,
            reduction = str_glue('umap.{reduction_method}'),
            pt.size = 1)

obj_subset$nFeature_CIRC_scaled = scale(obj_subset$nFeature_CIRC)
obj_subset$nCount_CIRC_scaled = scale(obj_subset$nCount_CIRC)


p1 = FeaturePlot(obj_subset,
                 features = 'nCount_CIRC_scaled',
                 order = T,
                 reduction = str_glue('umap.{reduction_method}'),
                 pt.size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'ivory', high = 'red', 
                        midpoint = 0
                        )
  # scale_color_gradient(low = 'blue', high = 'red')
p2 = FeaturePlot(obj_subset,
            features = 'nFeature_CIRC_scaled',
            order = T,
            reduction = str_glue('umap.{reduction_method}'),
            pt.size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'ivory', high = 'red', 
                        midpoint = 0)

ggpubr::ggarrange(p1, p2, ncol = 2)


VlnPlot(obj_subset, 
        features = c('nCount_CIRC','nFeature_CIRC'),
        group.by = 'CellType', 
        pt.size = 0,
        sort = T, 
        cols = list_cols)

# Summarize nCount, nFeature
obj_subset@meta.data %>%
  group_by(CellType) %>%
  summarise(nCell = n(), nCount = round(median(nCount_CIRC)), nFeature = round(median(nFeature_CIRC))) %>%
  arrange(desc(nCount)) %>%
  View

# Aggregate circRNA expression
DefaultAssay(obj_subset) = 'CIRC'
obj_subset_agg = AggregateExpression(obj_subset, 
                                     assay = 'CIRC',
                                     group.by = 'CellType',
                                     return.seurat = T)
obj_subset_agg = ScaleData(obj_subset_agg, assay = 'CIRC', features = rownames(obj_subset_agg))

obj_subset_agg$nCount = colSums(obj_subset_agg@assays$CIRC$counts)
obj_subset_agg$nFeature = colSums(obj_subset_agg@assays$CIRC$counts > 0)

obj_subset_agg@meta.data[c('nCount','nFeature')]%>%t()%>%data.frame%>%View

list_celltype = unique(obj_subset_agg$CellType)
feature_byCelltype = lapply(list_celltype, function(celltype){
  obj_subset_agg@assays$CIRC$counts[,obj_subset_agg$CellType == celltype] %>%
    as.data.frame %>%
    rownames_to_column('feature') -> df_count
  return(df_count$feature[df_count$. > 0])
})

names(feature_byCelltype) = list_celltype
saveRDS(feature_byCelltype, file = 'Robject_v2/feature_byCelltype.rds')



# list_ADARbound = readRDS(file = 'Robject_v2/list_ADARbound.v2.rds')
# 
# ggVennDiagram(
#   list(
#     GABA = feature_byCelltype$GABA,
#     ADAR1 = list_ADARbound$ADAR1_One,
#     ADAR2 = list_ADARbound$ADAR2_One,
#     ADAR3 = list_ADARbound$ADAR3_One
#   ),
#   list_names = c('GABA', 'ADAR1', 'ADAR2', 'ADAR3')
#   ) +
#   scale_fill_gradient(low="grey90",high = "steelblue") +
#   scale_colour_manual(values = c("black", "black", "black","black"))
# 
# ggVennDiagram(
#   list(
#     NPC = feature_byCelltype$NPC,
#     ADAR1 = list_ADARbound$ADAR1_One,
#     ADAR2 = list_ADARbound$ADAR2_One,
#     ADAR3 = list_ADARbound$ADAR3_One
#   ),
#   list_names = c('NPC', 'ADAR1', 'ADAR2', 'ADAR3')
# ) +
#   scale_fill_gradient(low="grey90",high = "steelblue") +
#   scale_colour_manual(values = c("black", "black", "black","black"))


################################################################################
# Venn Diagram by Cell Type Features
################################################################################
library(UpSetR)
feature_byCelltype = readRDS(file = 'Robject_v2/feature_byCelltype.rds')
feature_byCelltype%>%lapply(.,length)%>%unlist%>%sort
upset(fromList(feature_byCelltype))


################################################################################
# Check circRNA Expressions in Single-Cell Data
################################################################################
rm(list=ls())
gc()
reduction_method = 'harmony'
obj_subset = readRDS(file = str_glue('Robject_v2/brain_annotated_v2_circ.rds'))
annot = readRDS(file = 'Robject_v2/annot.rds')

# Color by celltype
list_celltype = levels(obj_subset$CellType)
list_cols = scales::hue_pal()(length(list_celltype))
names(list_cols) = list_celltype

# Load Core Enriched circRNA
DefaultAssay(obj_subset) = 'CIRC'
core_enriched = data.frame(readRDS('/home/local/suwanyu_971004/circRNA/results/RESTART/5.irCLASH_Brain/GSEA_results.rds'))
core_enriched = unlist(str_split(core_enriched[core_enriched$Description == 'ADAR3_One',]$core_enrichment, pattern = '/'))
core_enriched_sc = annot[annot$circRNA_id %in% core_enriched,]$circRNA_id_SC
core_enriched_sc = str_replace_all(core_enriched_sc, pattern = '\\|', replacement = '-')
core_enriched_sc_filtered = core_enriched_sc[core_enriched_sc %in% rownames(obj_subset)]


# Annotate whether circRNA are bound by ADAR
lapply(c('ADAR','ADARB1','ADARB2'), function(x){
  annotatedBSJs_withADARfam = data.frame(fread(str_glue('/home/local/suwanyu_971004/circRNA/results/4.ADARB2_irCLASH/annotatedBSJs_with{x}.tsv')))
  
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
list_circRNA_id = rownames(GetAssayData(obj_subset, 'CIRC'))
# Target at least one flanking intron
genesets_AUC = lapply(list_df, function(x){
  x$circRNA_id_SC = unlist(lapply(str_split(x$id, pattern = ':'), function(str_list){
    paste0(str_list[4], ':', str_list[5], '-', str_list[6])
  }))
  x = x[x$circRNA_id_SC %in% list_circRNA_id,]
  target = x$circRNA_id_SC[!(x$Group %in% c('Non_target','Moderate'))]
  # target = x$circRNA_id[x$Group == 'Target']
  return(target)
})
names(genesets_AUC) = c('ADAR1','ADAR2','ADAR3')

# Curated Gene Sets
genesets_AUC[['ADAR3_Enriched']] = core_enriched_sc_filtered


# Analyze circRNA obj
obj_subset = NormalizeData(obj_subset, assay = 'CIRC')
# Aggregate Expression
obj_subset_agg = AggregateExpression(obj_subset, 
                                      assay = 'CIRC',
                                      group.by = 'CellType',
                                      return.seurat = T)
obj_subset_agg = ScaleData(obj_subset_agg, assay = 'CIRC', features = rownames(obj_subset_agg))

if (!all(core_enriched_sc_filtered %in% rownames(obj_subset_agg@assays$CIRC$scale.data))){stop()} 


# Compare Expression
pseudobulk = GetAssay(obj_subset_agg, 'CIRC')$data%>%data.frame
pseudobulk$circRNA_id = rownames(pseudobulk) 
pseudobulk = pseudobulk[ncol(pseudobulk):1]

pseudobulk$Group = ifelse(pseudobulk$circRNA_id %in% genesets_AUC$ADAR3, 'Yes', 'No')
pseudobulk[-1] %>% na.omit %>% group_by(Group) %>% summarise_all(mean, na.rm = T) 



# # Expression Heatmap
# pheatmap::pheatmap(obj_subset_agg@assays$CIRC$data[core_enriched_sc_filtered,], 
#                   cluster_rows = T, 
#                   cluster_cols = T, 
#                   show_rownames = F, 
#                   show_colnames = T, 
#                   color = colorRampPalette(c("#327eba","white","#e06663"), bias = 3)(100))
# 
# # Scale Expression Heatmap
# pheatmap::pheatmap(obj_subset_agg@assays$CIRC$scale.data[core_enriched_sc_filtered,], 
#                    cluster_rows = T, 
#                    cluster_cols = T, 
#                    show_rownames = F, 
#                    show_colnames = T)
# 
# # nFeature & nCount of core-enriched circRNAs
# obj_subset_agg@assays$CIRC$data[core_enriched_sc_filtered,] %>%
#   data.frame() %>%
#   colSums() %>%
#   sort(., decreasing = T)
# 
# obj_subset_agg@assays$CIRC$data[core_enriched_sc_filtered,] %>%
#   data.frame() %>%
#   colMeans()%>%
#   sort(., decreasing = T)
# 
# obj_subset_agg@assays$CIRC$data[genesets_AUC[['ADAR3']],] %>%
#   data.frame() %>%
#   colMeans()%>%
#   sort(., decreasing = T) 
# 
# 
# obj_subset_agg@assays$CIRC$data[core_enriched_sc_filtered,] %>%
#   data.frame() %>% 
#   melt() %>%
#   ggplot(aes(x = reorder(variable, -value), y = value, fill = variable)) +
#   geom_violin() +
#   theme_classic2()
# 
# 
# obj_subset_agg@assays$CIRC$data[core_enriched_sc_filtered,] %>%
#   data.frame() %>%
#   apply(2, function(x) (ifelse(x != 0, 1, 0))) %>%
#   colSums()%>%
#   sort(., decreasing = T)
# 
# obj_subset_agg@assays$CIRC$data[core_enriched_sc_filtered,] %>%
#   data.frame() %>%
#   apply(2, function(x) (ifelse(x != 0, 1, 0))) %>%
#   colSums()
# 
# obj_subset@assays$CIRC$data[core_enriched_sc_filtered,] %>%
#   data.frame() %>%
#   apply(2, function(x) (ifelse(x != 0, 1, 0))) %>%
#   colSums() %>%
#   data.frame() -> ExpByCell
# 
# names(ExpByCell)= 'expression'
# ExpByCell$CellType = obj_subset$CellType
# 
# ExpByCell %>%
#   ggplot(aes(x = reorder(CellType, expression), y = expression, fill = expression)) +
#   geom_violin() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ExpByCell%>%
#   group_by(CellType)%>%
#   summarise(nCells = sum(expression != 0), nCells_Total = n(), prop = nCells/nCells_Total) %>%
#   arrange(desc(prop))
# 
# 
# # Calculate Module Scores
# obj_subset = AddModuleScore(obj_subset, 
#                              list(core_enriched_sc_filtered), 
#                              assay = 'CIRC',
#                              name = 'ADAR3_Enriched',
#                              seed = 42,
#                              slot = 'data',
#                              nbin = 64)
# obj_subset$ADAR3_Enriched = obj_subset$ADAR3_Enriched1
# obj_subset$ADAR3_Enriched1 = NULL
# 
# obj_subset = UCell::AddModuleScore_UCell(obj_subset, 
#                                           genesets_AUC,
#                                           slot = 'data')
# 
# p1 = VlnPlot(obj_subset,
#         features = 'ADAR1_UCell', 
#         group.by = 'CellType', 
#         sort = T, 
#         cols = list_cols)
# p2 = VlnPlot(obj_subset,
#         features = 'ADAR',
#         group.by = 'CellType', 
#         sort = T, 
#         cols = list_cols,
#         assay = 'RNA'
#         )
# ggpubr::ggarrange(p1,p2, common.legend = T, legend = 'right')
# 
# p1 = VlnPlot(obj_subset,
#         features = 'ADAR2_UCell', 
#         group.by = 'CellType', 
#         sort = T, 
#         cols = list_cols)
# p2 = VlnPlot(obj_subset,
#              features = 'ADARB1',
#              group.by = 'CellType', 
#              sort = T, 
#              cols = list_cols,
#              assay = 'RNA'
# )
# ggpubr::ggarrange(p1,p2, common.legend = T, legend = 'right')
# 
# 
# 
# p1 = VlnPlot(obj_subset,
#              features = 'ADARB2', 
#              group.by = 'CellType', 
#              sort = T, 
#              cols = list_cols,
#              assay = 'RNA')
# p2 = VlnPlot(obj_subset,
#              features = 'ADAR3_UCell', 
#              group.by = 'CellType',
#              sort = T, 
#              cols = list_cols)
# p3 = VlnPlot(obj_subset,
#              features = 'ADAR3_Enriched_UCell', 
#              group.by = 'CellType',
#              sort = T, 
#              cols = list_cols)
# ggpubr::ggarrange(p1,p2,p3, common.legend = T, legend = 'right', ncol = 3)




# # Summarize Ucell Module Scores
# obj_subset@meta.data%>%
#   data.frame%>%
#   group_by(scType)%>%
#   summarise(score_mean = mean(ADAR3_Enriched_UCell), score_median = median(ADAR3_Enriched_UCell)) %>%
#   arrange(desc(score_mean))
# 
# VlnPlot(pbmc_subset,
#         features = 'ADAR3_Enriched_UCell', 
#         group.by = 'scType', 
#         sort = T, 
#         cols = list_cols)+
#   ggsignif::geom_signif(comparisons = list(
#     c("Inhibitory", "Oligodendrocyte"),
#     c("Inhibitory", "NPC"),
#     c("Oligodendrocyte", "NPC")
#   ),
#   map_signif_level = TRUE, 
#   test = "wilcox.test",
#   y_position = 0.08)
# 
# 
# # Summarize Module Scores
# pbmc_subset@meta.data%>%
#   data.frame%>%
#   group_by(scType)%>%
#   summarise(score_mean = mean(ADAR3_Enriched1), score_median = median(ADAR3_Enriched1)) %>%
#   arrange(desc(score_median))
# 
# # Violin Plots
# VlnPlot(pbmc_subset, 
#         features = 'ADAR3_Enriched1', 
#         group.by = 'scType',
#         sort = T, cols = list_cols) +
#   ggsignif::geom_signif(comparisons = list(
#     c("Inhibitory", "Oligodendrocyte"),
#     c("Inhibitory", "NPC"),
#     c("Oligodendrocyte", "NPC")
#     ),
#     map_signif_level = TRUE, 
#     test = "wilcox.test",
#     y_position = 0.15)
# 
# 
# 
# # Calculate Module Scores and Statistical Test
# pbmc_subset$Group = ifelse(pbmc_subset$scType %in% c('Inhibitory', 'Oligodendrocyte', 'NPC'),
#                            pbmc_subset$scType, 'Others')
# 
# list_cols_2 = list_cols[names(list_cols)%in% c('Inhibitory', 'Oligodendrocyte','NPC')]
# list_cols_2[['Others']] = 'grey'
# VlnPlot(pbmc_subset, 
#         features = 'ADAR3_Enriched1', 
#         group.by = 'Group',
#         sort = T, cols = list_cols_2) +
#   ggsignif::geom_signif(comparisons = list(
#     # c("Inhibitory", "Oligodendrocyte"), 
#     c("Inhibitory", "Others"),
#     c("Oligodendrocyte", "Others"),
#     c("NPC", "Others")),
#               map_signif_level = TRUE, 
#               test = "t.test",
#               y_position = 0.15)
# 
# 
# VlnPlot(pbmc_subset,
#         features = 'nCount_CIRC',
#         group.by = 'Group',
#         sort = T, cols = list_cols_2) +
#   ggsignif::geom_signif(comparisons = list(c("Inhibitory", "Oligodendrocyte"), c("Inhibitory", "Others"), c("Oligodendrocyte", "Others")),
#                         map_signif_level = TRUE, 
#                         test = "wilcox.test",
#                         y_position = 2300)
# 
# 
# pbmc_subset_filtered = subset(pbmc_subset, subset = (scType != 'NPC'))
# pbmc_subset_filtered = AddModuleScore(pbmc_subset_filtered, 
#                              list(core_enriched_sc_filtered), 
#                              assay = 'CIRC',
#                              name = 'ADAR3_Enriched',
#                              seed = 42,
#                              slot = 'data')
# 
# 
# VlnPlot(pbmc_subset_filtered, 
#         features = 'ADAR3_Enriched1', 
#         group.by = 'Group',
#         sort = T, cols = list_cols_2) +
#   ggsignif::geom_signif(comparisons = list(
#     c("Inhibitory", "Oligodendrocyte"), 
#     c("Inhibitory", "Others"), 
#     c("Oligodendrocyte", "Others")
#     ),
#               map_signif_level = TRUE, 
#               test = "t.test",
#               y_position = 0.12)



################################################################################
# Overlap with RNA features
################################################################################
feature_byCelltype = readRDS(file = 'Robject_v2/feature_byCelltype.rds')

list_genes = list(
  GABA = feature_byCelltype$GABA,
  ADAR1 = genesets_AUC$ADAR1,
  ADAR2 = genesets_AUC$ADAR2,
  ADAR3 = genesets_AUC$ADAR3
)


# Venn Diagram Single cell And circAtlas
require(ggVennDiagram)
ggVennDiagram(list_genes)+
  scale_fill_gradient(low="grey90",high = "steelblue") +
  scale_colour_manual(values = c("black", "black", "black","black"))

annot = readRDS(file = 'Robject_v2/annot.rds')
ggVennDiagram(list(circSC = rownames(obj_subset_agg), 
                   RNAatlas = unlist(lapply(str_split(annot$id, pattern = ':'), function(str_list){
                     paste0(str_list[4], ':', str_list[5], '-', str_list[6])
                   }))), 
              list_names = c('circSC', 'RNAatlas'))+
  scale_fill_gradient(low="grey90",high = "steelblue") +
  scale_colour_manual(values = c("black", "black", "black","black")) +
  coord_flip()


# Universe = Intersection of circSC & RNAatlas
universe = intersect(rownames(obj_subset), unlist(lapply(str_split(annot$id, pattern = ':'), function(str_list){
  paste0(str_list[4], ':', str_list[5], '-', str_list[6])
})))

lapply(colnames(obj_subset_agg), function(celltype){
  celltype_id = feature_byCelltype[[celltype]]
  celltype_id = celltype_id[celltype_id%in%universe]
  
  # Hypergeometric test
  melting_geneset = reshape2::melt(genesets_AUC)
  colnames(melting_geneset) = c('hgnc_symbol','gs_name')
  melting_geneset = melting_geneset[c(2,1)]
  melting_geneset$gs_name%>%table
  melting_geneset = data.frame(melting_geneset)
  
  
  # GSA
  require(clusterProfiler)
  enrich = enricher(celltype_id, 
                    TERM2GENE = melting_geneset,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    pAdjustMethod = 'none',
                    minGSSize = 1,
                    maxGSSize = 10000,
                    universe = universe)
  enrich = as.data.frame(enrich)
  enrich$celltype = celltype
  return(enrich)
}) %>% do.call(rbind, .) -> result_gsa

result_gsa[result_gsa$pvalue < 0.05,!(colnames(result_gsa) == 'geneID')] %>% View
