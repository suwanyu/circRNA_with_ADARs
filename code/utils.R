get_data_linear = function(slot = NULL, annot = NULL, type = NULL, polyA = F){
  if (polyA == T){
    if (slot == 'data'){
      cts = read.csv(file = 'data/GSE138734_RNAAtlas_polyA_tpm_filtered.txt', sep= '\t')
      dim(cts)
    } else if (slot == 'count'){
      cts = read.csv(file = 'data/GSE138734_RNAAtlas_polyA_counts_filtered.txt', sep= '\t')
      dim(cts)
    } else {
      stop('error: slot name')
    }
  } else {
    if (slot == 'data'){
      cts = read.csv(file = 'data/GSE138734_RNAAtlas_totalRNA_tpm_filtered.txt', sep = '\t')
      dim(cts)
    } else if (slot == 'count'){
      cts = read.csv(file = 'data/GSE138734_RNAAtlas_totalRNA_counts_filtered.txt', sep = '\t')
      dim(cts)
    } else {
      stop('error: slot name')
    }
  }
  
  # Remove unnecessay column
  cts$genome_coordinates = NULL
  
  # Filter genes not annotated
  if (!is.null(annot)){
    cts = cts[cts$gene_name %in% unique(annot$gene),]
  }
  
  # Sanity Check
  if (names(cts)[1] != 'gene_name'){stop()}
  
  # Load metadata
  meta = read.csv(file = 'data/GSE138734_metadata_v2.txt', sep = '\t')
  meta$type %>% unique
  
  # order metadata by linear count data
  meta = meta[meta$RNA_Atlas_Number %in% names(cts)[-1],]
  
  # Convert Column names to Names
  names(cts)[-1] = sapply(names(cts)[-1], function(RNA_Atlas_Num){
    return(meta$name[meta$RNA_Atlas_Number == RNA_Atlas_Num])
  })
  
  # Only Filter Tissue Samples
  if (!is.null(type)){
    meta = meta[meta$type %in% type,]
    cts = cts[c('gene_name',meta$name)]
  }
  
  # Order cts & meta
  cts = cts[c('gene_name',meta$name)]
  
  names(cts)[1] = 'gene'
  
  return(list(cts = cts, meta = meta))
}


get_data = function(slot = NULL, annot = NULL, type = NULL, LCratio = NULL){
  # Load Count Matrix
  if (slot == 'rpm'){
    if (LCratio){
      cts = read.csv(file = 'data/LCratio_withNAandPseudocount_Milion.txt', sep = '\t')
    } else {
      cts = read.csv(file = 'data/GSE138734_RNAAtlas_circRNA_backspliced_rpm_filtered.txt', sep = '\t')
    }
  } else if (slot == 'count'){
    cts = read.csv(file = 'data/GSE138734_RNAAtlas_circRNA_backspliced_counts_filtered.txt', sep = '\t')
  } else {
    stop('error: slot name')
  }
  
  if (!is.null(annot) & !('circRNA_id' %in% names(annot))){
    annot$circRNA_id = lapply(str_split(annot$id, pattern = ':'), function(x){paste(x[4], x[5], x[6], x[2], sep = '_')}) %>% unlist
  }
  
  # Filter circRNA not annotated
  if (!is.null(annot)){
    cts = cts[cts$circRNA_id %in% annot$circRNA_id,]
  }
  
  # Sanity Check
  if (names(cts)[1] != 'circRNA_id'){stop()}
  
  # Load metadata
  meta = read.csv(file = 'data/GSE138734_metadata_v2.txt', sep = '\t')
  meta$type %>% unique
  # Filter metadata having circRNA ID
  meta = meta[meta$RNA_Atlas_Number %in% names(cts)[-1],]
  
  # Convert Column names to Names
  names(cts)[-1] = sapply(names(cts)[-1], function(RNA_Atlas_Num){
    return(meta$name[meta$RNA_Atlas_Number == RNA_Atlas_Num])
  })
  
  # Only Filter Tissue Samples
  if (!is.null(type)){
    meta = meta[meta$type %in% type,]
    cts = cts[c('circRNA_id',meta$name)]
  }
  
  list_circRNA_id_typespecific = readRDS('data/list_circRNA_id_typespecific.rds')
  cts = cts[cts$circRNA_id %in% list_circRNA_id_typespecific[[type]],]
  
  # Order cts & meta
  cts = cts[c('circRNA_id',meta$name)]
  
  return(list(cts = cts, meta = meta))
}