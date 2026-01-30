# RNA-seq Quantification and Exploratory Cross-Species Comparison of Placental Transporter Genes

## Overview
This repository provides an RNA-seq analysis pipeline for quantifying transcript and gene expression levels from mouse, rat, and human placental samples, and performing exploratory comparison of transporter gene expression (e.g., SLC and ABC families).

The pipeline covers:
- SRA data retrieval
- Quality control
- Adapter trimming
- Transcript-level quantification using Salmon
- RefSeq ID to gene symbol annotation
- Gene-level summarization using tximport
- TMM normalization and CPM calculation using edgeR
- Exploratory visualization (MA plots)

This workflow is intended for reproducible and exploratory analysis rather than formal differential expression testing.

---

## Data

### Mouse
- Paired-end RNA-seq  
- Labyrinth zone (mouse_LZ)

### Rat
- Paired-end RNA-seq  
- Labyrinth zone (rat_LZ)

### Human
- Syncytiotrophoblast RNA-seq (GSE182381)
  - SRR15514355 (single-end)
  - SRR15514356 (single-end)
  - SRR15514387 (paired-end)

Raw sequencing data are obtained from NCBI SRA.

---

## Requirements

### OS
- macOS Sequoia 15.3.2

### Terminal tools
- sra-tools 3.2.1  
- FastQC 0.11.9  
- Trimmomatic 0.39  
- Salmon 1.9.0  

### Reference transcripts (NCBI RefSeq)
- Mouse: GRCm39  
- Rat: GRCr8  
- Human: GRCh38.p14  

### R environment
- R 4.4.2  
- RStudio 2024.09.1+394  

Main R packages:
- tximport (1.34.0)  
- edgeR (4.4.2)  
- limma (3.62.2)  
- clusterProfiler (4.14.6)  
- org.Mm.eg.db (3.20.0)  
- org.Rn.eg.db (3.20.0)  
- org.Hs.eg.db (3.20.0)  
- dplyr (1.1.4)  
- data.table (1.17.0)  
- ggplot2 (3.5.2)  
- cowplot (1.1.3)

---

```

## Pipeline

### 1. Download SRA data
prefetch SRR15514355 SRR15514356 SRR15514387
fasterq-dump -p SRRxxxxx.sra --split-files
gzip *.fq

### 2. Quality control
fastqc *.fq.gz

### 3. Adapter trimming
Paired-end: TruSeq3-PE
Single-end: TruSeq3-SE
trimmomatic PE ...
trimmomatic SE ...

### 4. Post-trimming QC
fastqc trim_*.fq.gz

### 5. Build Salmon index
salmon index -t transcripts.fa -i transcripts_index

### 6. Transcript quantification
salmon quant -i transcripts_index -l A ...

### 7. Import to R
- Convert RefSeq IDs to gene symbols
- Remove RefSeq version suffix
- Subset NM_ transcripts
- Summarize to gene level using tximport

### 8. Normalization
- Create DGEList
- Filter low-expression genes
- TMM normalization
- CPM calculation

### 9. Visualization
MA plots of SLC and ABC family genes
```
---

## Output
- Gene-level count matrix
- CPM matrix
- MAplot figures (PNG)

---

## MA plots comparing:
- Human STB vs Rat LZ
- Mouse LZ vs Rat LZ

---

## Limitations
- This analysis is intended for exploratory comparison of transporter gene expression across species.
- Gene-level merging is performed using gene symbols without explicit ortholog mapping; therefore, strict cross-species quantitative comparison is not guaranteed.
- CPM values are used for visualization purposes only and are not intended for formal differential expression testing.
- Differences between RNA expression levels and protein abundance or transporter activity are not addressed in this pipeline.
- Human sample size is limited, which restricts assessment of intra-species variability.
- An ortholog-mapped version using Ensembl BioMart or HomoloGene is under development.

---

## Notes
- Example directory structure:

raw_fastq/
trimmed_fastq/
qc_reports/
salmon_quant/
scripts/
results/

- Thread numbers in terminal commands should be adjusted according to available computational resources.
- Raw sequencing data are not included in this repository (deposited in the NCBI Gene Expression Omnibus (GEO) repository under the accession number GSE295294).
