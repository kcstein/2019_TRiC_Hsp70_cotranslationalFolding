source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(data.table)
library(DESeq2)
library(dplyr)
library(plyr)


### TRiC
# Create a matrix with the read counts
rawcountData <- as.matrix(expression1[, .(T1X_total, T2X_total, R1X_total, R2X_total)])
row.names(rawcountData) <- as.character(expression1$orf)
rawcountData = rawcountData[rowSums(rawcountData) > 1, ]

# Define your variables
condition <- factor(c('TRiC','TRiC','Ribo','Ribo'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('T1X_total', 'T2X_total', 'R1X_total', 'R2X_total')
  
# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = rawcountData,
                              colData = colData,
                              design = ~ condition)
  
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "TRiC", "Ribo"), addMLE = T)


# This block will return a data.table with DESeq output
tric_deseq_gene <- as.data.frame(res)
setDT(tric_deseq_gene, keep.rownames = T)
colnames(tric_deseq_gene)[1] <- "orf"
setkeyv(tric_deseq_gene, c("orf"))
View(tric_deseq_gene[log2FoldChange > 1 & log2FoldChange < Inf & padj < 0.05])


### SSB
# Create a matrix with the read counts
rawcountData <- as.matrix(expression1[, .(ssb_Schx1_total, ssb_Schx2_total, ssb_Rchx1_total, ssb_Rchx2_total)])
row.names(rawcountData) <- as.character(expression1$orf)
rawcountData = rawcountData[rowSums(rawcountData) > 1, ]

# Define your variables
condition <- factor(c('SSB','SSB','Ribo','Ribo'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('ssb_Schx1_total', 'ssb_Schx2_total', 'ssb_Rchx1_total', 'ssb_Rchx2_total')

# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = rawcountData,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "SSB", "Ribo"), addMLE = T)


# This block will return a data.table with DESeq output to each codon position
ssb_deseq_gene <- as.data.frame(res)
setDT(ssb_deseq_gene, keep.rownames = T)
colnames(ssb_deseq_gene)[1] <- "orf"
setkeyv(ssb_deseq_gene, c("orf"))
View(ssb_deseq_gene[log2FoldChange > 1 & log2FoldChange < Inf & padj < 0.05])


### Atp2
# Create a matrix with the read counts
rawcountData <- as.matrix(expression_atp2[, .(T1_sum, T2_sum, R1_sum, R2_sum)])
row.names(rawcountData) <- as.character(expression_atp2$orf)
rawcountData = rawcountData[rowSums(rawcountData) > 1, ]

# Define your variables
condition <- factor(c('TRiC','TRiC','Ribo','Ribo'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('T1_sum', 'T2_sum', 'R1_sum', 'R2_sum')

# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = rawcountData,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "TRiC", "Ribo"), addMLE = T)


# This block will return a data.table with DESeq output
atp2_deseq_gene <- as.data.frame(res)
setDT(atp2_deseq_gene, keep.rownames = T)
colnames(atp2_deseq_gene)[1] <- "orf"
setkeyv(atp2_deseq_gene, c("orf"))
View(atp2_deseq_gene[log2FoldChange > 1 & log2FoldChange < Inf & padj < 0.05])



### Test for codon-level analysis for TRiC
# Create a matrix with the read counts
tric_deseq <- tric_dt[ribo_cor >= 0.5 & tric_cor >= 0.5 & 
                        R1X_rpc >= 0.5 & R2X_rpc >= 0.5 & T1X_rpc >= 0.5 & T2X_rpc >= 0.5 & 
                        ribo_total >= 128 & tric_total >= 128, c(1:7,167)]
tric_deseq <- tric_deseq[orf != "YDL133C-A" & orf != "YPL220W" & orf != "YDL184C"]

rawcountData <- as.matrix(tric_deseq[, .(T1X, T2X, R1X, R2X)])
row.names(rawcountData) <- as.character(tric_deseq$position1)
rawcountData = rawcountData[rowSums(rawcountData) > 1, ]

# Define your variables
condition <- factor(c('TRiC','TRiC','Ribo','Ribo'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('T1X', 'T2X', 'R1X', 'R2X')

# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = rawcountData,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "TRiC", "Ribo"), addMLE = T)


# This block will return a data.table with DESeq output
tric_deseq1 <- as.data.frame(res)
setDT(tric_deseq1, keep.rownames = T)
colnames(tric_deseq1)[1] <- "position1"
setkeyv(tric_deseq1, c("position1"))
setkeyv(tric_fishers1, c("position1"))
tric_deseq2 <- tric_deseq1[log2FoldChange > 0 & log2FoldChange < Inf & padj < 0.05]
View(tric_deseq1[tric_fishers1])

