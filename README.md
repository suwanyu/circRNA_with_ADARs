# circRNA analysis pipeline

This repository contains code and data used for the analysis of circRNA genesis and relationship with ADARs in RNA atlas dataset and single cell datasets.

**Datasets**
GEO RNA Atlas (GSE138734)
Description: Total RNA sequencing data processed using hg19. circRNAs identified and quantified using find_circ2 and CIRCexplorer2, then converted to hg38 using UCSC liftOver.
Samples: 45 tissues, 158 cell types
Data Used: Raw and RPM-normalized circRNA counts, raw and TPM-normalized linear RNA counts, TPM-normalized RNA counts from Poly-A RNA-seq

**Single-Cell Datasets**
Datasets: GSE67835, GSE71315, GSE75140, GSE125288
Description: Pre-processed raw count matrices from the CircSC database. circRNA quantification using CIRI2 and CIRIquant.
Selection: Healthy human brain tissues; GSE125288 (RamDA-seq), others (SMARTer)
Data Used: Quantified and pre-processed linear and circular RNA raw count matrices


