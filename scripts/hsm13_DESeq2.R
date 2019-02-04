## ----Loading packages, message=FALSE, warning=FALSE----------------------
library(tidyverse)
library(tximport)
library(DESeq2)
library(tximport)
library(ashr)

## ----counts_hist definition, eval=TRUE, message=FALSE, echo=FALSE--------
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

## ----TPM distribution plots, eval=TRUE, message=FALSE, warning=FALSE-----
counts_hist(hsm_L1, hsm_L2)

## ----Tabular data import and validation of file paths, eval=TRUE, message=FALSE, results='hide'----
samples <- read_tsv("samples.txt")
tx2gene <- read_csv("A_thal.tx2gene.csv")
dir <- #location of the data
files <- file.path(dir, paste0(samples$sample_id, ".isoforms.results"))
names(files) <- samples$sample_id
all(file.exists(files))
txi.rsem <- tximport(files, type = "rsem", tx2gene = tx2gene, txIn = TRUE, txOut = FALSE)

## ----Differential expression analysis with DESeq2, message=FALSE---------
dds <- DESeqDataSetFromTximport(txi.rsem, colData = samples, design = ~ Genotype + Treatment + Genotype:Treatment)
keep <- rowSums(counts(dds)) >= 10 #pre-filtering of large object to speed up DE analysis
dds <- dds[keep,]
dds$Genotype <- relevel(dds$Genotype, "WT")
dds <- DESeq(dds)
resultsNames(dds)

## ------------------------------------------------------------------------
rld <- rlog(dds, blind = FALSE)
ntd <- normTransform(dds)
library(vsn)
meanSdPlot(assay(rld))
meanSdPlot(assay(ntd))

## ------------------------------------------------------------------------
library(pheatmap)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Genotype", "Treatment")])
pheatmap(assay(rld)[select,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation = df)

## ------------------------------------------------------------------------
library(RColorBrewer)
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$Genotype, ntd$Treatment, sep = "_")
colnames(sampleDistMatrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, col=colours)

## ------------------------------------------------------------------------
plotPCA(ntd, intgroup=c("Genotype", "Treatment"))

## ----Printing of results tables------------------------------------------
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


## ----MA plots------------------------------------------------------------
p1 <- plotMA(res_1, ylim=c(-2,2), main="hsm13-leaf vs. WT-leaf")
p2 <- plotMA(res_2, ylim=c(-2,2), main="hsm13 roots (WCS365 vs buffer)")
p3 <- plotMA(res_3, ylim=c(-2,2), main="hsm13 vs WT (roots in buffer)")
p4 <- plotMA(res_4, ylim=c(-2,2), main="WT roots (WCS365 vs buffer)")
p5 <- plotMA(res_5, ylim=c(-2,2), main="hsm13 vs WT (roots with WCS365")

## ----Data export---------------------------------------------------------
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


