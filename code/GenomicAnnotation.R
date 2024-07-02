################################################################################
# Genomic Annotation with circRNAprofiler & GenomicRanges
################################################################################

library(circRNAprofiler)
library(gridExtra)
library(data.table)

# Load input files and formatting for circRNAprofiler
cts = read.csv(file = 'data/GSE138734_RNAAtlas_circRNA_backspliced_counts_filtered.txt', sep = '\t')
featuredata = read.csv(file = 'data/GSE138734_circRNA.txt')
gtf = formatGTF(pathToGTF = 'data/gencode.v38.annotation.gtf')

featuredata = featuredata[featuredata$circRNA_id %in% cts$circRNA_id,]
featuredata = featuredata[c('consensus_gene', 'strand', 'chr', 'start', 'end')]
featuredata[featuredata==''] = NA
featuredata = na.omit(featuredata)

# Input formatting
input = data.frame(
  id = paste(featuredata$consensus_gene, featuredata$strand, featuredata$chr, 
             featuredata$chr, featuredata$start, featuredata$end, sep = ':'),
  gene = featuredata$consensus_gene,
  strand = featuredata$strand,
  chrom = featuredata$chr,
  startUpBSE = ifelse(featuredata$strand =='-', as.numeric(featuredata$end), as.numeric(featuredata$start)),
  endDownBSE = ifelse(featuredata$strand =='-', as.numeric(featuredata$start), as.numeric(featuredata$end))
)

annotatedBSJs = annotateBSJs(input, gtf, isRandom = F)
write.csv(annotatedBSJs, file = 'data/annotatedBSJs.csv')

# Functions for Annotating Alu elements
intron2bed = function(annotatedBSJs){
  # Input: Output from circRNAprofiler
  ## startUpIntron: is the 5’ coordinate of the intron immediately upstream the acceptor site in the selected transcript.
  ## endUpIntron: is the 3’ coordinate of the intron immediately upstream the acceptor site in the selected transcript.
  ## startDownIntron: is the 5’ coordinate of the intron immediately downstream the donor site in the selected transcript.
  ## endDownIntron: is the 3’ coordinate of the intron immediately downstream the donor site in the selected transcript.
  
  annotatedBSJs_filtered = na.omit(annotatedBSJs)
  
  UpIntron = data.frame(
    chr = annotatedBSJs_filtered$chrom,
    start = ifelse(annotatedBSJs_filtered$strand =='-', 
                   as.numeric(annotatedBSJs_filtered$endUpIntron),
                   as.numeric(annotatedBSJs_filtered$startUpIntron)),
    end = ifelse(annotatedBSJs_filtered$strand =='-', 
                 as.numeric(annotatedBSJs_filtered$startUpIntron),
                 as.numeric(annotatedBSJs_filtered$endUpIntron)),
    Intron_id = paste0(annotatedBSJs_filtered$id, '_UpIntron'),
    score = 0,
    strand = annotatedBSJs_filtered$strand
  )
  
  DownIntron = data.frame(
    chr = annotatedBSJs_filtered$chrom,
    start = ifelse(annotatedBSJs_filtered$strand =='-', 
                   as.numeric(annotatedBSJs_filtered$endDownIntron),
                   as.numeric(annotatedBSJs_filtered$startDownIntron)),
    end = ifelse(annotatedBSJs_filtered$strand =='-', 
                 as.numeric(annotatedBSJs_filtered$startDownIntron),
                 as.numeric(annotatedBSJs_filtered$endDownIntron)),
    Intron_id = paste0(annotatedBSJs_filtered$id, '_DownIntron'),
    score = 0,
    strand = annotatedBSJs_filtered$strand
  )
  
  Intron = rbind(UpIntron, DownIntron)
  
  return(Intron)
}
FindAluWithinIntron = function(fn_intron, fn_alu){
  require(GenomicRanges)
  
  intron <- read.table(fn_intron, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand"))
  alu <- read.table(fn_alu, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand"))
  
  gr_intron <- GRanges(seqnames = intron$chrom, ranges = IRanges(start = intron$start, end = intron$end), strand = intron$strand,
                       intron_id = intron$name)
  gr_alu <- GRanges(seqnames = alu$chrom, ranges = IRanges(start = alu$start, end = alu$end), strand = alu$strand,
                    alu_id = alu$name)
  gr_merged = mergeByOverlaps(gr_alu, gr_intron, type = 'within', ignore.strand = TRUE)
  merged = as.data.frame(gr_merged)
  
  
  # Sanity Check
  if (!all(as.character(merged$gr_alu.seqnames) == as.character(merged$gr_intron.seqnames))){
    stop('something wrong')
  }
  
  merged$gr_alu.seqnames = NULL
  merged$gr_alu.alu_id = NULL
  merged$gr_intron.seqnames = NULL
  merged$gr_intron.start = NULL
  merged$gr_intron.end = NULL
  merged$gr_intron.width = NULL
  merged$gr_intron.intron_id = NULL
  
  merged$id = str_split_i(merged$intron_id, pattern = '_', i = 1)
  merged$Intron = str_split_i(merged$intron_id, pattern = '_', i = 2)
  merged$intron_id = NULL
  
  names(merged) = c('Alu_start','Alu_end', 'Alu_width', 'Alu_strand', 'Alu_id', 'Intron_strand', 'id','Intron_position')
  
  return(merged)
}
Annotate_Alu = function(annotatedBSJs, AluWithIntron){
  DownIntron = AluWithIntron[AluWithIntron$Intron_position == 'DownIntron',]
  UpIntron = AluWithIntron[AluWithIntron$Intron_position == 'UpIntron',]
  
  names(DownIntron)[1:5] = paste(names(DownIntron)[1:5], 'DownIntron', sep = '_')
  names(UpIntron)[1:5] = paste(names(UpIntron)[1:5], 'UpIntron', sep = '_')
  
  Flatten_DF = function(DF, Intron_position = NULL){
    DF = apply(DF, MARGIN = 2, function(x){
      x = as.character(x)
      x = str_replace_all(x, pattern = ' ', replacement = '')
    }) %>% data.frame
    
    str_flatten_comma = function(x){
      str_flatten(x,collapse = ",")
    }
    
    DF = DF %>%
      reframe(across(str_replace_all(c("Alu_start_DF",
                                       "Alu_end_DF",
                                       "Alu_width_DF",
                                       "Alu_strand_DF",
                                       "Alu_id_DF",
                                       "Intron_strand",
                                       "Intron_position"), 
                                     pattern = 'DF', 
                                     replacement = Intron_position)
                     , str_flatten_comma),
              .by = id) %>%
      distinct()
    
    # Filter columns before merge
    DF$Intron_position = NULL
    DF$Intron_strand = NULL
    
    return(DF)
  }
  
  UpIntron = Flatten_DF(UpIntron, 'UpIntron')
  DownIntron = Flatten_DF(DownIntron, 'DownIntron')
  
  # Merge with annotation file
  annotatedBSJs_withAlu = left_join(annotatedBSJs, UpIntron, by = 'id')
  annotatedBSJs_withAlu = left_join(annotatedBSJs_withAlu, DownIntron, by = 'id')
  return(annotatedBSJs_withAlu)
}

# bed formatting
Intron = intron2bed(annotatedBSJs)
write.table(Intron, file = 'data/FlankingIntron.bed6', sep = '\t', quote = F,row.names = F,col.names = F)

# Find Alu elements in flanking introns
fn_intron = 'data/FlankingIntron.bed6'
fn_alu = "data/Alu_only.bed"
AluWithIntron = FindAluWithinIntron(fn_intron, fn_alu)

# Merge annotation file with Alu element information
annotatedBSJs_withAlu = Annotate_Alu(na.omit(annotatedBSJs), AluWithIntron)
fwrite(annotatedBSJs_withAlu, 
       file = 'data/annotatedBSJs_withAlu.tsv',
       sep = '\t',
       row.names = F)

################################################################################
# Genomic Overlap with ADARB2 Binding Duplex and Flanking Intron
################################################################################

# Load ADAR family binding table
lapply(1:3, function(sheet){
  read.xlsx('data/ADARB2_irCLASH/ADARfam_Binding.xlsx', sheet = sheet, startRow = 2)
}) -> list_df

names(list_df) = c('ADAR','ADARB1','ADARB2')

# Gene of Interest
gene_of_interest = 'ADARB2'
ADARB2_binding = list_df[[gene_of_interest]]

# Subset row with NA in start, end, chrom
ADARB2_binding = ADARB2_binding[!is.na(ADARB2_binding$Start) & !is.na(ADARB2_binding$End) & !is.na(ADARB2_binding$Chromosome),]
ADARB2_binding = ADARB2_binding[!str_detect(ADARB2_binding$Start, pattern = ';') & !str_detect(ADARB2_binding$End, pattern = ';'),]

ADARB2_binding$Start = as.numeric(ADARB2_binding$Start)
ADARB2_binding$End = as.numeric(ADARB2_binding$End)

# LiftOver
library(liftOver)
library(rtracklayer)
require(GenomicRanges)
require(stringr)

binding2bed = function(df){
  # Hybrid.id Arm Chromosome Strand     Start       End Genic.location Gene.name Repeat.annotation
  bed = data.frame(chr = df$Chromosome,
                   start = df$Start,
                   end = df$End,
                   name = paste(
                     df$Subtrate.id,
                     # df$Hybrid.id,
                     df$Arm,
                     df$Genic.location,
                     df$Gene.name,
                     df$Repeat.annotation,
                     df$Affinity,
                     sep = ';'),
                   score = 0,
                   strand = df$Strand)
  return(bed)
}
ADARB2_binding_bed = binding2bed(ADARB2_binding)
gr_adarb2 = GRanges(seqnames = ADARB2_binding_bed$chr, ranges = IRanges(start = ADARB2_binding_bed$start, end = ADARB2_binding_bed$end), strand = ADARB2_binding_bed$strand, name = ADARB2_binding_bed$name)

ch = import.chain('data/hg19ToHg38.over.chain')
seqlevelsStyle(gr_adarb2) = "UCSC"
gr_adarb2_hg38 = liftOver(gr_adarb2, ch)
gr_adarb2_hg38 = unlist(gr_adarb2_hg38)


# Overlap with flanking intron
intron = read.table('data/FlankingIntron.bed6', header = FALSE,col.names = c("chrom", "start", "end", "name", "score", "strand"))
gr_intron = GRanges(seqnames = intron$chrom, ranges = IRanges(start = intron$start, end = intron$end), strand = intron$strand,intron_id = intron$name)

gr_merged = mergeByOverlaps(gr_adarb2, gr_intron, type = 'any', ignore.strand = TRUE)
print(nrow(gr_merged))
rm(gr_merged) # Just for check

gr_merged_hg38 = mergeByOverlaps(gr_adarb2_hg38, gr_intron, type = "any", ignore.strand = TRUE)
print(nrow(gr_merged_hg38))

merged = as.data.frame(gr_merged_hg38)

# Sanity Check
if (!all(as.character(merged$seqnames) == as.character(gr_adarb2_hg38$seqnames))) {
  stop("seqnames are not matched")
}

merged$gr_adarb2_hg38.seqnames = NULL
merged$gr_adarb2_hg38.name = NULL
merged$gr_intron.seqnames = NULL
merged$gr_intron.start = NULL
merged$gr_intron.end = NULL
merged$gr_intron.width = NULL
merged$gr_intron.intron_id = NULL

# merged$Hybrid.id = sapply(strsplit(merged$name, ";"), function(x) x[1])
merged$Substrate.id = sapply(strsplit(merged$name, ";"), function(x) x[1])
merged$Arm = sapply(strsplit(merged$name, ";"), function(x) x[2])
merged$Genic.location = sapply(strsplit(merged$name, ";"), function(x) x[3])
merged$Gene.name = sapply(strsplit(merged$name, ";"), function(x) x[4])
merged$Repeat.annotation = sapply(strsplit(merged$name, ";"), function(x) x[5])
merged$Affinity = sapply(strsplit(merged$name, ";"), function(x) x[6])
merged$name = NULL

merged$id = str_split_i(merged$intron_id, pattern = '_', i = 1)
merged$Intron = str_split_i(merged$intron_id, pattern = '_', i = 2)
merged$intron_id = NULL

names(merged) = c('ADARN_start','ADARN_end','ADARN_width','ADARN_strand','Intron_strand',
                  'Substrate.id','Arm','Genic.location','Gene.name','Repeat.annotation',
                  'Affinity','id','Intron_position')

merged = merged[c('ADARN_start','ADARN_end','ADARN_width','ADARN_strand','Substrate.id',
                  'Arm','Genic.location','Gene.name','Repeat.annotation','Affinity',
                  'Intron_strand','id','Intron_position')
]

UpIntron = merged[merged$Intron_position == 'UpIntron',]
DownIntron = merged[merged$Intron_position == 'DownIntron',]

names(UpIntron)[1:10] = paste(names(UpIntron)[1:10], 'UpIntron', sep = '_')
names(DownIntron)[1:10] = paste(names(DownIntron)[1:10], 'DownIntron', sep = '_')


Flatten_DF = function(DF, Intron_position = NULL){
  DF = apply(DF, MARGIN = 2, function(x){
    x = as.character(x)
    x = str_replace_all(x, pattern = ' ', replacement = '')
  }) %>% data.frame
  
  str_flatten_comma = function(x){
    str_flatten(x,collapse = ",")
  }
  
  DF = DF %>%
    reframe(across(str_replace_all(c("ADARN_start_DF",
                                     "ADARN_end_DF",
                                     "ADARN_width_DF",
                                     "ADARN_strand_DF",
                                     "Substrate.id_DF",
                                     "Arm_DF",
                                     "Genic.location_DF",
                                     "Gene.name_DF",
                                     "Repeat.annotation_DF",
                                     "Affinity_DF",
                                     "Intron_strand",
                                     "Intron_position"), 
                                   pattern = 'DF', 
                                   replacement = Intron_position)
                   , str_flatten_comma),
            .by = id) %>%
    distinct()
  
  # Filter columns before merge
  DF$Intron_position = NULL
  DF$Intron_strand = NULL
  
  return(DF)
}
UpIntron = Flatten_DF(UpIntron, 'UpIntron')
DownIntron = Flatten_DF(DownIntron, 'DownIntron')


# Annotation circRNAs
annotatedBSJs = read.csv(file = 'data/annotatedBSJs.csv', row.names = 1)
annotatedBSJs = na.omit(annotatedBSJs)

# Merge with Annotation file
annotatedBSJs_withADARB2 = left_join(annotatedBSJs, UpIntron, by = 'id')
annotatedBSJs_withADARB2 = left_join(annotatedBSJs_withADARB2, DownIntron, by = 'id')


fwrite(annotatedBSJs_withADARB2,
       str_glue('data/ADARB2_irCLASH/annotatedBSJs_with{gene_of_interest}.tsv'),
       sep = '\t',
       row.names = FALSE)
