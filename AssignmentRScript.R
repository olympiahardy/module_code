library(DESeq2)

# Creation of my experimental design table
colTable <- data.frame("ID" = c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12"), "Group" = c("A","A","A","A","B","B","B","B","C","C","C","C"), "Replicate" = c(1,2,3,4,1,2,3,4,1,2,3,4))
write.csv(colTable, "experimentaldesigntable.csv", row.names = F)
# Creation of gene count table and dds object for DeSeq2 removing all zero rows before differential analysis
geneCountTable <- read.csv("gene_count_matrix.csv", row.names = 1)
dds1 <- DESeqDataSetFromMatrix(countData = geneCountTable, colData = colTable, design = ~ Group)
notAllZero <- (rowSums(counts(dds1))>0)
dds1 <- dds1[notAllZero,]
dds1 <- DESeq(dds1)

# Generation of my transcript count table
transcriptCountTable <- read.csv("transcript_count_matrix.csv", row.names = 1)

# Creation of dds object for DeSeq2 removing all zero rows before differential analysis
dds2 <- DESeqDataSetFromMatrix(countData = transcriptCountTable, colData = colTable , design = ~ Group)
notAllZero <- (rowSums(counts(dds2))>0)
dds2 <- dds2[notAllZero,]
dds2 <- DESeq(dds2)

# Plotting of dispersion plots for gene counts and transcript counts
par(mfrow = c(1,1))
plotDispEsts(dds1, xlab = "Mean of Normalised Counts", ylab = "Dispersion", main = "Dispersion Plot of Gene Counts")
plotDispEsts(dds2, xlab = "Mean of Normalised Counts", ylab = "Dispersion", main = "Dispersion Plot of Transcript Counts")


# Creation of r-log based PCA plots for both gene count data and transcript level data
rld1 <- rlog(dds1, blind = FALSE)
PCAGene <- plotPCA(rld1, intgroup=c("Group"))
PCAGene + labs(title = "PCA for Gene Count Data")
??ggtitle
rld2 <- rlog(dds2, blind = FALSE)
plotPCA(rld2, intgroup=c("Group"))

# Differential Expression of gene count data

# Differential Expression of group A vs group B where LFC=0
res1AvB <- results(dds1, contrast = c("Group", "A","B"))
summary(res1AvB,alpha=0.05)
head(res1AvB)
res_sort1AvB <- res1AvB[order(res1AvB$padj),]
res_sig1AvB <- subset(res_sort1AvB, res_sort1AvB$padj < 0.001)
head(res_sig1AvB)
write.csv(res_sig1AvB, "Significant_Results_AvB.csv", row.names = TRUE)

# Differential Expression of group A vs group B where LFC=2
res.lfc2AvB <- results(dds1, contrast = c("Group", "A", "B"), lfcThreshold = 2)
summary(res.lfc2AvB,alpha=0.001)
res.lfc2AvB_sort <- res.lfc2AvB[order(res.lfc2AvB$padj),]
res.lfc2AvB_sig <- subset(res.lfc2AvB_sort, res.lfc2AvB_sort$padj < 0.001)
write.csv(res.lfc2AvB_sig, "Significant_Results_AvB_2.csv", row.names = TRUE)

# Differential Expression of group A vs group C where LFC=0
res1AvC <- results(dds1, contrast = c("Group", "A","C"))
summary(res1AvC,alpha=0.05)
head(res1AvC)
res_sort1AvC <- res1AvC[order(res1AvC$padj),]
res_sig1AvC <- subset(res_sort1AvC, res_sort1AvC$padj < 0.001)
head(res_sig1AvC)
write.csv(res_sig1AvC, "Significant_Results_AvC.csv", row.names = TRUE)

# Differential Expression of group A vs group C where LFC=2
res.lfc2AvC <- results(dds1, contrast = c("Group", "A", "C"), lfcThreshold = 2)
summary(res.lfc2AvC,alpha=0.001)
res.lfc2AvC_sort <- res.lfc2AvC[order(res.lfc2AvC$padj),]
res.lfc2AvC_sig <- subset(res.lfc2AvC_sort, res.lfc2AvC_sort$padj < 0.001)
write.csv(res.lfc2AvC_sig, "Significant_Results_AvC_2.csv", row.names = TRUE)

# Differential Expression of group C vs group B where LFC=0
res1CvB <- results(dds1, contrast = c("Group", "C","B"))
summary(res1CvB,alpha=0.05)
head(res1CvB)
res_sort1CvB <- res1CvB[order(res1CvB$padj),]
res_sig1CvB <- subset(res_sort1CvB, res_sort1CvB$padj < 0.001)
head(res_sig1CvB)
write.csv(res_sig1CvB, "Significant_Results_CvB.csv", row.names = TRUE)

# Differential Expression of group C vs group B where LFC=2
res.lfc2CvB <- results(dds1, contrast = c("Group", "C", "B"), lfcThreshold = 2)
summary(res.lfc2CvB,alpha=0.001)
res.lfc2CvB_sort <- res.lfc2CvB[order(res.lfc2CvB$padj),]
res.lfc2CvB_sig <- subset(res.lfc2CvB_sort, res.lfc2CvB_sort$padj < 0.001)
write.csv(res.lfc2CvB_sig, "Significant_Results_CvB_2.csv", row.names = TRUE)

par(mfrow = c(2,3))

# MA plots of Group A vs B where LFC=0 and LFC=2
DESeq2::plotMA(res1AvB, xlab = "Mean of Normalised Counts", ylab = "Log Fold Change", main = "MA Plot for Group A vs B where LFC=O")
DESeq2::plotMA(res.lfc2AvB, xlab = "Mean of Normalised Counts", ylab = "Log Fold Change", main = "MA Plot for Group A vs B where LFC=2")

# MA plots of Group A vs C where LFC=0 and LFC=2
DESeq2::plotMA(res1AvC, xlab = "Mean of Normalised Counts", ylab = "Log Fold Change", main = "MA Plot for Group A vs C where LFC=0")
DESeq2::plotMA(res.lfc2AvC, xlab = "Mean of Normalised Counts", ylab = "Log Fold Change", main = "MA Plot for Group A vs C where LFC=2")

# MA plots of Group C vs B where LFC=0 and LFC=2
DESeq2::plotMA(res1CvB, xlab = "Mean of Normalised Counts", ylab = "Log Fold Change", main = "MA Plot for Group C vs B where LFC=0")
DESeq2::plotMA(res.lfc2CvB, xlab = "Mean of Normalised Counts", ylab = "Log Fold Change", main = "MA Plot for Group C vs B where LFC=2")
