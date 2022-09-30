library("DESeq2")
library("ggplot2")
library("fdrtool")
library("pheatmap")
library("ReportingTools") #per HTMLReport

# RNA-seq workflow: gene-level exploratory analysis and differential expression
## Starting from SummarizedExperiment to a DESeqDataSet object

coldataa = read.csv("coldata.csv")
class(coldataa)
coldataa$Condition

count_matrix = read.csv("count_matrix.csv", row.names = 1)
class(count_matrix)


Condition<- as.factor(coldataa[["Condition"]]);Condition
Gender<- as.factor(coldataa[["Gender"]]);Gender
Group<- as.factor(coldataa[["Group"]]);Group

coldatafinal<-data.frame(coldataa[1],Condition,Gender,coldataa["Age"],coldataa["ApoE"],Group,coldataa["Smell_test"])
coldatafinal

class(coldataa[[3]])
class(coldatafinal[[3]])

coldatafinal$Condition <- relevel(coldatafinal$Condition, "control")

result<-DESeqDataSetFromMatrix(countData = count_matrix, colData = coldatafinal, design = ~ Group + Condition)


#Exploratory analysis and visualization

#The first step is to pre-filter data. 
#In this case we removing rows of the DESeqDataSet that have no counts, 
#or only a single count across all samples
nrow(result)

#Basic filtering, (if equal to 0 is not expressed in any samples)
result <- result[rowSums(counts(result)) > 1, ]
nrow(result)

#Strict filtering, for each row I sum how many quanti values are > 10, I sum and I check if > 4
#at least 5 samples in which the gene are expressed more than 10
keep <- apply(counts(result), 1, function(x){
  sum(x>10) > 4
})
table(keep)

result <- result[keep,]

#stabilize the variance across the mean

#first way
rld <- rlog(result, blind=FALSE)
head(assay(rld), 3)

#second way, filtering more
vsd <- vst(result, blind=FALSE)
head(assay(vsd), 3)

#plot the result
par(mfrow=c(1, 3))
dds <- estimateSizeFactors(result)
lims <- c(-2, 20)

plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3, main="log2(x + 1)", xlim=lims, ylim=lims)

plot(assay(rld)[,1:2],
     pch=16, cex=0.3, main="rlog", xlim=lims, ylim=lims)

plot(assay(vsd)[,1:2],
     pch=16, cex=0.3, main="VST", xlim=lims, ylim=lims)


## Sample distance
#visualize sample-to-sample distances is a principal components analysis
# I can't do clustering on the result of plotPCA, only visualization
plotPCA(rld, intgroup = c("Group", "Condition"))
plotPCA(rld, intgroup = c("Group"))
plotPCA(rld, intgroup = c("Condition"))

## Other types of plot that you can use
## MDS plot
## t-SNE plot

## Differential expression analysis with DESeq2
#using pre-existing size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testingv
dds <- DESeq(dds)
res <- results(dds)
res

#show the significant genes with the strongest down or up regulation
#padj, must be small so indicate significance of the result
resSig <- subset(res, padj < 0.1)
nrow(resSig)

#padj contestualize the significance of the result (ex. if it's alone or with others genes)
head(resSig[ order(resSig$log2FoldChange), ])

#log2foldchange, is - or +, it's important to has an absolute big value
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])
summary(res)

#more strict threshold, more significant result
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
summary(res.05)

#raise the log2 fold change threshold to test for genes that show more substantial changes due to treatment
#expression of gene treat/expression of genes in control/untre
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

##compair 2 cell lines
#In general, the results for a comparison of any two levels of a variable can be extracted using the contrast argument to results.
results(dds, contrast=c("Group", "B", "C"))

## Log fold change shrinkage
#ranking and visualization, without the need for arbitrary filters on low count genes
#Genes with low counts or high dispersion values could lead to extremely high and not realistic LogFoldChange values. The shrunken log fold changes are useful for ranking and visualization, without the need for arbitrary filters on low count genes. Three different shrinkage methods are available
resultsNames(dds)
resAsh <- lfcShrink(dds, coef="Intercept", type="ashr")
table(resAsh$padj < .05)

## Plotting results
#A quick way to visualize the counts for a particular gene is to use the plotCounts function that takes as arguments the DESeqDataSet, a gene name, and the group over which to plot the counts
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("Group"))

#ggplot2
#We can also make custom plots using the ggplot function For example we can show the best gene based on treatment and cell line
geneCounts <- plotCounts(dds, gene=topGene, intgroup=c("Group","Condition"), returnData=TRUE)
ggplot(geneCounts, aes(x=Condition, y=count, color=Condition)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3)

ggplot(geneCounts, aes(x=Condition, y=count, fill=Condition)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")

ggplot(geneCounts, aes(x=Condition, y=count, color=Condition, group=Condition)) +
  scale_y_log10() + geom_point(size=3) + geom_line()

#MA-plot
#provides a useful overview for an experiment with a two-group comparison
plotMA(res, ylim=c(-5,5))

plotMA(resAsh, ylim=c(-3,3))

#hist
#nother useful diagnostic plot is the histogram of the p values - genes with mean normalized count larger than 1 - genes with small count are excluded
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")

## P-value distribution
hist(res$pvalue, col = "royalblue4")

#remove filtered out genes by independent filtering,
#they have NA adj. pvals
res <- res[ !is.na(res$padj),] 

#remove genes with NA pvals (outliers)
res <- res[ !is.na(res$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on
res <- res[, -which(names(res) == "padj")]

#we can re-estimate the variance inside the model and correct the p-value distribution.
# use z-scores as input to FDRtool to re-estimate the p-value
res_fdr <- fdrtool(res$stat, statistic= "normal")

#add values to the results data frame, also ad new BH- adjusted p-values
res[,"padj"] <- p.adjust(res_fdr$pval, method = "BH")

#plot correct p-value distribution 
hist(res_fdr$pval, col = "royalblue4", xlab = "CORRECTED p-values")

##Gene clustering
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20) #rowVars takes the variance of each gene  among the samples, in each row.

mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("Group","Condition")])

#data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
pheatmap(mat, annotation_col=df)



## Gene annotations

library("AnnotationDbi")
library("org.Hs.eg.db")

#this is to remove informations after the point, we don't need them

j<-row.names(res)

remove_upgrades <- function(j) {
  l<-unlist(strsplit(j, "[.]"))[1]
}

o<-sapply(j,remove_upgrades)
names(o)<-NULL

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=o,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=o,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

#show cols of org.Hs.eg.db
columns(org.Hs.eg.db)

## Exporting results
resOrderedDF <- as.data.frame(res)[1:100,]
write.csv(resOrderedDF, file="results.csv")

#export results in a dynamic HTML document
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")

publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)
