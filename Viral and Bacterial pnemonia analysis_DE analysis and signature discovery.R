

##### DEGs, Pathway analysis and diagnostic Signature and validation 

# Load necessary packages
library(DESeq2)
library(RUVSeq)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(WGCNA)
library(flashClust)
library(dynamicTreeCut)
library(genefilter)
library(Hmisc)
library(biomaRt)

# Set working directory if necessary
# setwd("your/directory/path")

# Load your data (replace with your actual data)
# counts_merged <- read.xlsx("counts_file.xlsx")
# pheno_merged <- read.xlsx("pheno_file.xlsx")

# Define DESeqDataSet design
ddsp <- DESeqDataSetFromMatrix(countData = counts_merged,
                               colData = pheno_merged,
                               design = ~ sex + batch + pneumonia.type)

# Set reference level for pneumonia.type
ddsp$pneumonia.type <- relevel(ddsp$pneumonia.type, ref = "Viral")

# Run DESeq with likelihood ratio test (LRT)
ddsp_for_ruv <- DESeq(ddsp, test = "LRT", reduced = ~ 1)

# Get DESeq results
res <- results(ddsp_for_ruv)

# Create a SeqExpressionSet for RUVSeq
set <- newSeqExpressionSet(counts = counts(ddsp_for_ruv))

# Between-lane normalization (normalizing across samples)
idx <- rowSums(counts(set) > 5) >= 2
set <- set[idx,]
set <- betweenLaneNormalization(set, which = "upper")

# Apply RUVg method for removing unwanted variation
not_sig <- rownames(res)[which(res$pvalue > .1)]
empirical <- rownames(set)[rownames(set) %in% not_sig]
set <- RUVg(set, empirical, k = 2)

# Add RUVg W components to DESeqDataSet
ddsp$W1 <- set$W_1

# Redefine design including W components and other terms of interest
design(ddsp) <- ~ W1 + sex + batch + pneumonia.type

# Run DESeq again with the new design
ddsp2 <- DESeq(ddsp)

# Filter low-expressed genes
ddsp <- ddsp[rowSums(counts(ddsp)) > 5, ]

# Estimate size factors and dispersion estimates
ddsp2 <- estimateSizeFactors(ddsp2)
ddsp2 <- estimateDispersionsGeneEst(ddsp2)
# Variance-stabilizing transformation
vsd <- varianceStabilizingTransformation(ddsp2, blind = FALSE)
# Prepare data for batch effect removal and visualization 
vsd2 <- vsd 
covariates <- model.matrix(~ batch + sex + W1, data = colData(vsd2))
design.full <- model.matrix(~ pneumonia.type, data = colData(vsd2))
mat <- limma::removeBatchEffect(assay(vsd2), covariates = covariates, design = design.full)
# Update the matrix in vsd2
assay(vsd2) <- mat

# Save the transformed matrix
write.csv(mat, file = "mat_transformed.csv", row.names = FALSE)

# PCA analysis
pcaData <- plotPCA(vsd2, intgroup = c("batch", "pneumonia.type", "sex", "W1"), returnData = TRUE, ntop = 100)

# Save PCA plot
pdf("PCA_b_v.pdf", height = 8, width = 12)
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = pneumonia.type, label = name)) +
  geom_point(size = 12, alpha = 0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
print(pca_plot)
dev.off()

##### Get significant DEG results 

lfc.cutoff <- 1.5
resruv <- results(ddsp2, alpha = 0.05)
res_r<-as.data.frame(resruv)
res_r<-na.omit(res_r) 
res05_r<-res_r[res_r$padj < 0.05&(res_r$log2FoldChange>lfc.cutoff|res_r$log2FoldChange<(-lfc.cutoff)),]
Toptable<- res05_r
Toptable<-Toptable_B_V_05[Toptable_B_V_05$baseMean>100,]

# Save DEG results in CSV files
write.csv2(res05_r, "Pneumo_TopDEGs.csv") 

#############  MODULE ANALYSIS 
#### Load data
normalized_counts<-mat
dim(normalized_counts)
meta<-readRDS("pheno_merged.rds")
table(colnames(normalized_counts)==rownames(meta))
normalized_counts<-as.data.frame(normalized_counts)
dim(normalized_counts)

#### Filter non-protein coding genes

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    "ensembl_gene_id",
    "entrezgene_id"
  ),filters='biotype', values=c('protein_coding'),
  uniqueRows = TRUE)

annotLookup<-annotLookup[!annotLookup$hgnc_symbol=="",]
dim(annotLookup)
normalized_counts<-normalized_counts[rownames(normalized_counts)%in%annotLookup$hgnc_symbol,]
dim(normalized_counts)

#### Filter less variable genes
rv <- rowVars(data.matrix(normalized_counts))
q25 <- quantile(rowVars(data.matrix(normalized_counts)), .75)
normalized_counts_f <- data.matrix(normalized_counts)[rv > q25, ]
dim(normalized_counts_f)

#### Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
datExpr0 = as.data.frame(t(normalized_counts_f));
dim(datExpr0)

# Remove missing values and outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster samples to see outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# datExpr now contains the expression data ready for network analysis.
datExpr<-datExpr0

### Upload pheno info
datTraits<-meta
datTraits$ID<-rownames(datTraits)
table(rownames(datExpr)==rownames(datTraits))

#### Power selection
allowWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = "signed")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",#ylim=c(0,1),
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",#ylim=c(0,10000),
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#### Network construction and module detection
TOM = TOMsimilarityFromExpr(datExpr, TOMType = "signed",networkType="signed",power = 20)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9) 
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = T, pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods) # Give the module information
dynamicColors = labels2colors(dynamicMods) 
table(dynamicColors)
sizeGrWindow(8,6) 
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
MEList = moduleEigengenes(datExpr, colors = dynamicColors,softPower=20) 
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6) 
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors

#### Relating modules with the trait
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs = orderMEs(mergedMEs)
moduleTraitCor = WGCNA::cor(MEs, as.numeric(datTraits$Condition), use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
table_mod<-data.frame(moduleTraitCor,moduleTraitPvalue)
table_mod[table_mod$moduleTraitPvalue<0.05,]
sizeGrWindow(20,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation 
labeledHeatmap(Matrix = moduleTraitCor,
               yLabels = rownames(moduleTraitCor),
               xLabels="Viral-Bacterial",
               setStdMargins=F,
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 1.3,
               zlim = c(-1,1),
               cex.lab.x = 1.3,
               cex.lab.y = 1.7,
               main = paste("Module-trait relationships"))


#### Create a dataframe with GS, MMs and pMMs values 

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Condition)
names(weight) = "Phenotype"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, as.numeric(datTraits$Condition), use = "p"));

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

geneInfo0 = data.frame(
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, as.numeric(datTraits$Condition), use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Phenotype));
geneInfo = geneInfo0[geneOrder, ]


### One-step hub gene
hub<-chooseTopHubInEachModule(datExpr,moduleColors,type="signed", power=20)

###### Lasso models with most DE genes 
 
matrix<-t(mat)
matrix<-as.data.frame(matrix)  
matriz<-matriz[,colnames(matrix)%in%rownames(Toptable)] 
pheno_merged_filtered<-pheno_merged[,c(2,57)] #Select the condition colum of interest from the pheno dataset 
matrix<-merge(matrix,pheno_merged_filtered, by.x= 0,by.y=0 #Merge expresion and pheno data with the condition column 
rownames(matrix)<-matrix$Row.names 
matrix<-matrix[,-1] #remove Row.names column and the other column with pheno information 


set.seed(124)
cv.lassoModel <- cv.glmnet(
  x=data.matrix(matrix[,-37]),
  y=matriz$pneumonia.type,
  standardize=TRUE,
  alpha=1.0,
  nfolds=10,
  family="multinomial",
  parallel=TRUE) 
dev.off()

# plot variable deviances vs. shrinkage parameter, λ (lambda)

plot(cv.lassoModel) 

#Now we can identify best predictors via 2 metrics: min λ (lamba); 1 standard error of λ


idealLambda <- cv.lassoModel$lambda.min
idealLambda1se <- cv.lassoModel$lambda.1se

# derive coefficients for each gene
co <- coef(cv.lassoModel, s=idealLambda, exact=TRUE)
co

co.se <- coef(cv.lassoModel, s=idealLambda1se, exact=TRUE)
co.se

# identify predictors for each sub-type
rownames(co.se$viral)[which(co.se$viral != 0)]
rownames(co.se$ybacterial)[which(co.se$ybacterial != 0)] 

#######build a new model with our best predictors
matrix$pneumonia.type<-gsub("bacterial","1",matriz$pneumonia.type)
matrix$pneumonia.type<-gsub("viral","0",matriz$pneumonia.type)
matrix$pneumonia.type<-as.numeric(matriz$pneumonia.type)

finalLasso <- glm(matriz$pneumonia.type ~ FAM20A + BAG3 + MXRA7+ TDRD9+ KLF14 ,
                  data=matriz,
                  family=binomial(link="logit")) 

saveRDS(finalLasso,"finalLasso.rds") 


####### OPTIMAL CUTOFF for Binary Classification #########

library(cutpointr)
cp <- cutpointr(predictors_data, predictors, pneumonia.type, pos_class = "Bacterial", neg_class = "Viral", direction = ">=",
                method = maximize_metric, metric = sum_sens_spec) 
plot(cp)
summary(cp) 

roc <- roc(matriz$pneumonia.type, fitted(finalLasso), smooth=FALSE) 

##### Test the model in the dataset including probable samples (test group) or other validation datasets 

genes<- c("FAM20A","BAG3","MXRA7","TDRD9","KLF14","pneumonia.type")  

matriz_probable<-readRDS("normalized_probable.rds") #Load the matrix containing the neccesarry samples 
pheno_merged_filtered_probable<-readRDS("pheno_merged_filtered_probable.rds")
matriz_probable<-t(matprobable)
View(matriz_probable)
matriz_lasso_p<-matriz_probable[,colnames(matriz_probable)%in%genes]
matriz_lasso_p<-merge(matriz_lasso_p,pheno_merged_filtered_probable, by.x= 0,by.y=0)
View(matriz_lasso_p)
rownames(matriz_lasso_p)<-matriz_lasso_p$Row.names
matriz_lasso_p<-matriz_lasso_p[,-c(1,7)]
View(matriz_lasso_p) 

score<-(0.05964*matriz_lasso_p[,1])-(1.25018*matriz_lasso_p[,2])+(2.18318*matriz_lasso_p[,4])+(1.57463*matriz_lasso_p[,3])-(0.53523*matriz_lasso_p[,5])-16.31345

roc(matriz_lasso_p$pneumonia.type~score) 











