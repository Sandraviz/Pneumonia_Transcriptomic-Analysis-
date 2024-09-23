#### DE Analysis of RNA-Seq and Microarray data, correlation and Pathway analysis 

# Load necessary libraries
library(DESeq2)
library(RUVSeq)
library(EDASeq)
library(limma)
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(enrichplot)

# Differential expression analysis from RNA-Seq data with RUV-Seq and DESeq2
# Load your data 
# counts_merged <- read.xlsx("counts_file.xlsx")
# pheno_merged <- read.xlsx("pheno_file.xlsx")

# 1. Data Preprocessing
# Create DESeq2 dataset with RNA-Seq count data
ddsp2 <- DESeqDataSetFromMatrix(countData = counts_merged, colData = pheno_merged, design =~ sex + batch + tube + Condition)
ddsp2$Condition <- relevel(ddsp2$Condition, ref = "Control")

# 2. Initial DESeq2 analysis using RUVSeq
ddsp_for_ruv <- DESeq(ddsp2, test = "LRT", reduced = ~1)
res <- results(ddsp_for_ruv)

# 3. RUVSeq Normalization
set <- newSeqExpressionSet(counts(ddsp_for_ruv))
idx <- rowSums(counts(set) > 5) >= 2
set <- set[idx, ]
set <- betweenLaneNormalization(set, which = "upper")

# Extract empirical controls for RUV normalization
not_sig <- rownames(res)[which(res$pvalue > 0.1)]
empirical <- rownames(set)[rownames(set) %in% not_sig]
set <- RUVg(set, empirical, k = 2)
pdat <- pData(set)

# Update DESeq2 object with RUVSeq factors
ddsp2$W1 <- set$W_1
ddsp2$W2 <- set$W_2
design(ddsp2) <- ~ W1 + W2 + sex + batch + tube + Condition
ddsp3 <- DESeq(ddsp2)
ddsruvClean <- ddsp3[which(mcols(ddsp3)$betaConv), ]
resruv <- results(ddsruvClean, alpha = 0.05)

# 4. Selection of Differentially Expressed Genes (DEGs)
res_r <- as.data.frame(resruv)
res_r <- na.omit(res_r)
res05_r <- res_r[res_r$padj < 0.05 & (abs(res_r$log2FoldChange) > 1), ]

# Final table of DEGs
Pneumo_Toptable <- res05_r

# 5. Remove batch effect using variance-stabilizing transformation and limma
ddsp2 <- ddsp2[rowSums(counts(ddsp2)) > 5, ]
ddsp2 <- estimateSizeFactors(ddsp2)
ddsp2 <- estimateDispersionsGeneEst(ddsp2)
vsd <- varianceStabilizingTransformation(ddsp2, blind=FALSE)
# Create covariates for batch effect removal
batch <- vsd$batch
sex <- vsd$sex
tube <- vsd$tube
W1 <- vsd$W1
W2 <- vsd$W2

# Remove batch effects
covariates <- model.matrix(~batch + sex + tube + W1 + W2)
design.full <- model.matrix(~Condition, colData(vsd))
mat <- limma::removeBatchEffect(assay(vsd), covariates = covariates, design = design.full)
assay(vsd) <- mat

# 6. PCA plot visualization

pdf("PCA_plot.pdf", height = 8, width = 12)
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("batch", "Condition", "sex", "tube", "W1", "W2"), returnData = TRUE, ntop = 100)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()
print(pcaPlot)
dev.off()

# 7. Microarray DE Analysis

# Download and preprocess microarray data
gse <- getGEO("GSE103119", GSEMatrix = TRUE, getGPL = FALSE)
phenoGSE103119 <- pData(gse[[1]])

# Read and preprocess expression data
mydataRAW <- read.ilmn("GSE103119_non-normalized.txt", probeid="ID_REF", expr="836385")
mydata <- as.data.frame(mydataRAW$E)
mydata[mydata < 0] <- 0  # Set negative values to zero

# Normalize data
GSE103119exprs <- normalizeBetweenArrays(log2(mydata + 1), method="quantile")

# Probe annotation and averaging
probeID <- rownames(GSE103119exprs)
geneID <- data.frame(Gene = unlist(mget(x = probeID, envir = illuminaHumanv4SYMBOL)))
geneID <- na.omit(geneID)
GSE103119exprs <- GSE103119exprs[rownames(geneID), ]
rownames(GSE103119exprs) <- geneID$Gene
GSE103119exprs <- avereps(GSE103119exprs, ID = rownames(GSE103119exprs))

# 8. Limma DE Analysis
condition <- as.factor(phenoGSE103119$Control)
gender <- as.factor(phenoGSE103119$GENDER)
design.full <- model.matrix(~0 + condition + gender)
fit <- lmFit(GSE103119exprs, design.full)
contrast.matrix <- makeContrasts(PNEUMOvsCTL = conditionPNEUMO - conditionCTL, levels = colnames(design.full))
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

# DEGs from microarray
PNEUMO_top.table_Array <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
PNEUMO_top.table_Array <- PNEUMO_top.table_Array[abs(PNEUMO_top.table_Array$logFC) >= 1 & PNEUMO_top.table_Array$adj.P.Val < 0.05, ]


# 9. Correlation between RNA-Seq and Microarray DEGs

#Pneumo_Toptable
#PNEUMO_top.table_Array

PNEUMO_top.table_Array<-PNEUMO_top.table_Array[abs(PNEUMO_top.table_Array$logFC)>1.&PNEUMO_top.table_Array$adj.P.Val<0.05,]

table(rownames(Pneumo_Toptable)%in%rownames(PNEUMO_top.table_Array))

Genes_comunes<-rownames(Pneumo_Toptable[rownames(Pneumo_Toptable)%in%rownames(PNEUMO_top.table_Array),])


Pneumo_top_LFC_common<-Pneumo_Toptable[rownames(Pneumo_Toptable)%in%Genes_comunes,]
Pneumo_top_LFC_common<-Pneumo_top_LFC_common[,-c(1,3,4,5)]
colnames(Pneumo_top_LFC_common)<-c("LFC_RNASeq","pAdjust_RNA-Seq")
Pneumo_top_LFC_common$ID<-rownames(Pneumo_top_LFC_common)

dim(PNEUMO_top.table_Array)
dim(Pneumo_Toptable)

PNEUMO_top_LFC_common<-PNEUMO_top.table_Array[rownames(PNEUMO_top.table_Array)%in%Genes_comunes,]
PNEUMO_top_LFC_common<-PNEUMO_top_LFC_common[,-c(2,3,4,6)]
colnames(PNEUMO_top_LFC_common)<-c("LFC_Array","pAdjust_Array")
PNEUMO_top_LFC_common$ID<-rownames(PNEUMO_top_LFC_common)

Merged_common<-merge(Pneumo_top_LFC_common,PNEUMO_top_LFC_common, by= "ID")
Merged_common$Diff<-abs(Merged_common$LFC_RNASeq)-abs(Merged_common$LFC_Array)

#### Correlation  scatter

library(ggscatter)
library(ggpubr)

a<-ggscatter(Merged_common, x = "LFC_RNASeq", y = "LFC_Array",
             #add = "reg.line", 
             show.legend.text = "",
             #title="SEVERE VS Mild Aschenbrenner et al. (2021)",
             color="Diff",# Add regression line
             #conf.int = TRUE, 
             point=TRUE,
             cor.coef=F,
             label = Merged_common$ID,
             #label.select=labels_sel,
             font.label = c(12, "bold","gray44"),
             font.family = "",
             xlab = "LFC_RNASeq",
             ylab = "LFC_Array",
             repel = TRUE,
             label.rectangle = TRUE,# Add confidence interval
             alpha = 0.8,
             size=5)+
  #add.params = list(color = "blue",fill = "lightgray"))+
  gradient_color(c('dodgerblue','Black', 'yellow'))+
  #geom_point(aes(colour=padj_NonCovid_VS_Cont),size=5)+
  #scale_colour_gradient(low = "blue", high = "red")+
  ggpubr::stat_cor(method = "pearson", label.x = -5, label.y = 5,size=5.5)+
  scale_y_continuous(breaks=seq(-7, 7, 1))+
  scale_x_continuous(breaks=seq(-7, 7, 1))+
  theme_bw()+theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.text=element_text(size=15),axis.title=element_text(size=15))+geom_vline(aes(xintercept = 0))+geom_hline(aes(yintercept = 0))

#labs(title="Area Vs Population")        
#guides(colour=guide_legend("Padjusted RNA-Seq")) 
a<-ggpar(a,legend.title="LFC differences",ylim = c(-7, 7),xlim = c(-7, 7))
print(a)


# 10. Pathway Enrichment Analysis
# Gene annotation with biomaRt
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id'), mart = mart, uniqueRows = TRUE)
annotLookup_f <- annotLookup[annotLookup$hgnc_symbol %in% rownames(common_genes),]
annotLookup_f <- annotLookup_f[!duplicated(annotLookup_f$hgnc_symbol) & !is.na(annotLookup_f$entrezgene_id),]

# ORA using GO Biological Process (BP)
geneList_ORA <- sort(common_genes$log2FoldChange, decreasing = TRUE)
ego <- enrichGO(gene = names(geneList_ORA), keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH")
dat_ego <- as.data.frame(ego@result[ego@result$p.adjust < 0.05,])

# Reactome pathway enrichment
reactP <- enrichPathway(gene = names(geneList_ORA), organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)
dat_RE <- as.data.frame(reactP@result[reactP@result$p.adjust < 0.05,])






