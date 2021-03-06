Differential Expression Analysis of *A. thaliana* hsm13
================

Yi initially ran the data in [Salmon](https://combine-lab.github.io/salmon/), and was able to get very close count values between his biological replicates except in one sample. After talking with him, I decided to run [RSEM](https://github.com/deweylab/RSEM) instead. The reads were mapped against the Arabidopsis thaliana TAIR10 reference [transcriptome](ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/arabidopsis_thaliana/). I decided to load each of the counts files into R and use histograms to visualize counts.

### Loading required packages

``` r
library(tidyverse)
library(tximport)
library(DESeq2)
library(ashr)
```

Validation of TPM data
----------------------

One of Yi's issues was that there was a single sample with all FPKM values ~5 times lower than it's biological replicate. I wrote a small function called `counts_hist` that plots the distribution of log transformed counts in two `sample.genes.results` files.

``` r
hsm_L1 <- read_tsv("data/hsm-L1.isoforms.results")
hsm_L2 <- read_tsv("data/hsm-L2.isoforms.results")


counts_hist <- function(file1, file2){ #files 1 and 2 are counts matrices that contain a column called "TPM"
  tpm1 <- file1$TPM
  tpm2 <- file2$TPM
  f1_lc <- density(log(tpm1)) #contains log transformed counts values for file1
  f2_lc <- density(log(tpm2)) #contains log transformed counts values for file2
  plot(f1_lc$x, f1_lc$y*length(tpm1),
       type="l",
       main="TPM Distribution",
       xlab="log10(TPM)",
       ylab="Number of Transcripts")
  lines(f2_lc$x, f2_lc$y*length(tpm2), col="red")
}


counts_hist(hsm_L1, hsm_L2)
```

All of the TPM distributions were very similar between samples, so I decided to proceed with DEG analysis using the DESeq2 package.

Loading data into DESeq2
------------------------

I used tximport to load the TPM data into R. First, I generated a file called `samples.txt` that includes basic information about the samples. I will need to use this later for DESeq2 as well.

``` r
samples <- read_tsv("samples.txt")
tx2gene <- read_csv("A_thal.tx2gene.csv")
dir <- "/Users/andrew/compbio/2019/220119_HSM_RNASeq/data/"
files <- file.path(dir, paste0(samples$sample_id, ".isoforms.results"))
names(files) <- samples$sample_id
all(file.exists(files))
txi.rsem <- tximport(files, type = "rsem", tx2gene = tx2gene, txIn = TRUE, txOut = FALSE)
```

``` r
dds <- DESeqDataSetFromTximport(txi.rsem, colData = samples, design = ~ Genotype + Treatment + Genotype:Treatment)
keep <- rowSums(counts(dds)) >= 10 #pre-filtering of large object to speed up DE analysis
dds <- dds[keep,]
dds$Genotype <- relevel(dds$Genotype, "WT")
dds <- DESeq(dds)
resultsNames(dds)
```

Data QC
-------

Before exporting the data and sending it to Yi for GO enrichment, I ran some QC on the data to make sure that everything was working well. In these analyses, I expected to see that the leaf samples clustered apart from the root samples.

### Transform the data for visualization

``` r
rld <- rlog(dds, blind = FALSE)
ntd <- normTransform(dds)

library(vsn)
meanSdPlot(assay(rld))
meanSdPlot(assay(ntd))
```

### Heatmap

This code produces a heatmap of the top 20 genes by count in all of the samples. It can be useful to see if the replicates are similar to one another. The next two code blocks provide other ways to visualize how samples cluster together.

``` r
library(pheatmap)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Genotype", "Treatment")])
pheatmap(assay(rld)[select,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation = df)
```

### Distance matrix

``` r
library(RColorBrewer)
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$Genotype, ntd$Treatment, sep = "_")
colnames(sampleDistMatrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, col=colours)
```

### PCA

``` r
plotPCA(ntd, intgroup=c("Genotype", "Treatment"))
```

Results and Analysis
--------------------

Since the experimental design is quite complicated we want to pull out DEGs for 4 different comparisons:

-   hsm\_L vs. WT\_L
-   hsm\_R365 vs. hsm\_RCT
-   hsm\_RCT vs. WT\_RCT
-   WT\_R365 vs. WT\_RCT

### Selecting the data for results tables

``` r
#hsm13 leaf vs wt leaf
res_1 <- lfcShrink(dds, contrast=list(c("Genotype_hsm13_vs_WT","Genotypehsm13.TreatmentLeaf")), type="ashr")
#hsm13 roots with WCS365 vs buffer
res_2 <- lfcShrink(dds, contrast=list(c("Treatment_WCS365_vs_Buffer", "Genotypehsm13.TreatmentWCS365")), type = "ashr")
#hsm13 roots with buffer vs wt roots with buffer
res_3 <- lfcShrink(dds, "Genotype_hsm13_vs_WT", type = "ashr")
#wt roots with WCS365 vs bufffer
res_4 <- lfcShrink(dds, "Treatment_WCS365_vs_Buffer", type = "ashr")
#hsm13 roots with WCS365 vs wt roots with WCS365
res_5 <- lfcShrink(dds, contrast=list("Genotype_hsm13_vs_WT", "Treatment_WCS365_vs_Buffer"), type="ashr")
```

### MA plots

MA plots are different from volcano plots because instead of comparing p-value vs. fold change, they compare fold change and mean counts. Volcano plots can be made in ggplot.

``` r
p1 <- plotMA(res_1, ylim=c(-2,2), main="hsm13-leaf vs. WT-leaf")
p2 <- plotMA(res_2, ylim=c(-2,2), main="hsm13 roots (WCS365 vs buffer)")
p3 <- plotMA(res_3, ylim=c(-2,2), main="hsm13 vs WT (roots in buffer)")
p4 <- plotMA(res_4, ylim=c(-2,2), main="WT roots (WCS365 vs buffer)")
p5 <- plotMA(res_5, ylim=c(-2,2), main="hsm13 vs WT (roots with WCS365")
```

### Filtering and export of results

``` r
res_1F <- subset(res_1, padj < 0.1)
res_2F <- subset(res_2, padj < 0.1)
res_3F <- subset(res_3, padj < 0.1)
res_4F <- subset(res_4, padj < 0.1)
res_5F <- subset(res_5, padj < 0.1)

write.csv(as.data.frame(res_1F), "DE_tables/hsm_leaf_vs_wt_leaf.csv")
write.csv(as.data.frame(res_2F), "DE_tables/hsm_wcs365_vs_hsm_buffer.csv")
write.csv(as.data.frame(res_3F), "DE_tables/hsm_buffer_vs_wt_buffer.csv")
write.csv(as.data.frame(res_4F), "DE_tables/wt_wcs365_vs_wt_buffer.csv")
write.csv(as.data.frame(res_5F), "DE_tables/hsm_wcs365_vs_wt_wcs365.csv")
```
