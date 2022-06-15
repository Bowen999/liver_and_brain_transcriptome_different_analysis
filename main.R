#0.data and package
library(pheatmap)
library(RColorBrewer)
library(BiocManager)
library(org.Mm.eg.db)
library(DESeq2)
library(genefilter)
library(edgeR)
library(EnhancedVolcano)

#1.Characterization of sample relationships
#1.1heatmaps
exp_counts <- read.table(file="~/Desktop/NGS_Proj/Output/Feature Count with Column name.csv", sep=",", header = T)
row.names(exp_counts) <- exp_counts$gene
exp_counts$gene <- NULL
#head(exp_counts)

#1.1.1Creating a DESeqDataSet from the Input Data
tissue <- c("brain", "brain", "brain", "liver", "liver", "liver")
rname <- c("brain1", "brain2", "brain3", "liver1", "liver2", "liver3")
sample_data <- data.frame(tissue, row.names=rname)
sample_data

data_deseq <- DESeqDataSetFromMatrix(countData = exp_counts, colData = sample_data, design = ~ 1)
#head(counts(data_deseq))

#1.1.2Filtering the DataSet
data_deseq <- data_deseq[ rowSums(counts(data_deseq)) > 1, ]
rld <- rlog(data_deseq, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
#sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$tissue, sep="-" )
colnames(sampleDistMatrix) <- paste( rld$tissue, sep="-"  )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#1.2.PCA
plotPCA(rld, intgroup = c("tissue"))

#1.3.MDS
y <- DGEList(counts=exp_counts[,1:6], gene=row.names(exp_counts)) #counts is the read counts col
#Organizing sample into groups
#Calculating Library Scaling Factors
y <- calcNormFactors(y)
y$samples$group = c("brain", "brain", "brain", "liver", "liver", "liver")
plotMDS(y)

#1.4.gene heatmap
#calculate the variance
geneVars <- rowVars(assay(rld))
#re-Order in decreasing
geneVarsOrdered <- order(geneVars, decreasing = TRUE)
topVarGenes <- head(geneVarsOrdered, 2000)
# Make a matrix
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)

#df <- as.data.frame(colData(rld)[,"tissue"])
#clear_col_names <- paste( rld$tissue, sep=".")
#topGenesHeatmap <- pheatmap(mat, annotation_col=df, labels_col = clear_col_names)
topGenesHeatmap <- pheatmap(mat, show_rownames = FALSE)

#2.DEG
tissue <- factor(tissue)
design <- model.matrix(~tissue)
rownames(design) <- colnames(y)
#design

Volcano_data <- estimateDisp(y, design, robust=TRUE)
Volcano_data$common.dispersion

#2.1 Num of up&down regualted
fit <- glmFit(Volcano_data, design)
lrt <- glmLRT(fit)
summary(de2 <- decideTestsDGE(lrt))

#2.2 Volcano plot
EnhancedVolcano(lrt$table,
                lab = rownames(lrt$table),
                x = 'logFC',y = 'PValue',
                pCutoff = 10e-30, FCcutoff = 2,
                pointSize = 1, labSize = 3,shape =16,
                col = c("azure3", "dodgerblue2", "slateblue", "royalblue4"))

diffExpGenes2 <- topTags(lrt, n=1000, p.value = 0.05)
head(diffExpGenes2$table)



#Differential Expression Analysis in Exact Test
y <- estimateDisp(y)
y$common.dispersion
#common dispersion is 0.2144381
#y$tagwise.dispersion
plotBCV(y)

et <- exactTest(y, pair=c("brain","liver"))
summary(de <- decideTestsDGE(et))

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

#Differentially Expressed Genes
library(org.Hs.eg.db)
diffExpGenes <- topTags(et, n=1000, p.value = 0.05)

head(diffExpGenes$table)
DEG <- diffExpGenes$table






# gene_names = rownames(diffExpGenes$table)
# mapping <- AnnotationDbi::select(org.Hs.eg.db, 
#                                    keys = gene_names,
#                                    columns = c("ENTREZID", "SYMBOL"),
#                                    keytype = "SYMBOL")
# 
# missing <- is.na(mapping$ENTREZID)
# sum(missing)
# mapping <- mapping[!missing,]
# 
# d2 <- duplicated(mapping$ENTREZID)
# sum(d2)
# mapping <- mapping[!d2,]
# 
# row.names(mapping) <- mapping$ENTREZID
# head(mapping) 
# gene_list = rownames(mapping) 
# head(gene_list)




# gene_list <- mapping["ENTREZID"]
# 
# DEG$ENTREZID <- gene_list
# 
# missing <- is.na(DEG$ENTREZID)
# sum(missing)
# DEG <- DEG[!missing,]
# head(DEG)
# d2 <- duplicated(DEG$ENTREZID)
# sum(d2)
# DEG <- DEG[!d2,]
# 
# #Exporting Differential Expression Results
write.table(DEG,
             file="DEG_Brain_vs_Liver.txt",
             sep = "\t", row.names=TRUE, col.names=NA)
# 
# head(gene_list)
# write.table(gene_list,
#             file="gene_list.txt",
#             sep = "\t", row.names=TRUE, col.names=NA)








#Enrichment Analysis

library(clusterProfiler)
#GOBP Enrichment Analysis

enrich_GO <- enrichGO(gene = gene_list,
                           OrgDb = org.Hs.eg.db,
                           ont = "BP", readable = T)

write.table(enrich_GO, file="GOBP_enrichment_results_R.txt", sep = "\t", row.names=TRUE, col.names=NA)


#barplot
barplot(enrich_GO, x = "GeneRatio", color = "p.adjust", showCategory = 20)
#dotplot
dotplot(enrich_GO, x = "Count", color = "pvalue", size = "GeneRatio",showCategory=20)



#KEGGG Enrichment Analysis
library(DOSE)
enrich_KEGG <- enrichKEGG(gene_list, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", 
                 qvalueCutoff=0.1)
head(enrich_KEGG)
write.table(enrich_KEGG, file="KEGG_enrichment_results_R.txt", sep = "\t", row.names=TRUE, col.names=NA)
#barplot
barplot(enrich_KEGG, x = "GeneRatio", color = "p.adjust", showCategory = 20)
#dotplot
dotplot(enrich_KEGG, x = "Count", color = "pvalue", size = "GeneRatio",showCategory=20)




#GO Network Visualization
fold_change <- diffExpGenes$table$logFC
names(fold_change) <- diffExpGenes$table$SYMBOL

cnetplot(enrich_KEGG, foldChange = fold_change,
         colorEdge = TRUE, cex_label_gene = 0.5, showCategory = 5)

#GO Heatmap
heatplot(enrich_KEGG, foldChange=fold_change, showCategory = 10)
library(GOexpress)
heatmap_GO(enrich_GO)


gseGO(gene_List, ont = "BP", OrgDb, keyType = "ENTREZID", exponent = 1,
      nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
      pAdjustMethod = "BH", verbose = TRUE, seed = FALSE, by = "fgsea")


head(gene_list)

write.table(gene_list, file="Gene_list.txt", sep = ",", row.names=TRUE, col.names=NA)
