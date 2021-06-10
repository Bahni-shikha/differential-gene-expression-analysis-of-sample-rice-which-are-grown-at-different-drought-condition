setwd("C:/projects/rice/pipeline_result/heat_drought/data/workflow_SE/results/featureCounts")
counts=read.table("counts.txt", sep="", head=T, skip=1, row.names="Geneid")
countsnew <- counts [(c(6:11))]
countsnew
colnames(countsnew)
library(stringr)
colnames(countsnew)= str_split_fixed(colnames(countsnew),"\\.",6)[,1]
colnames(countsnew)


#define conditions for the samples
mycols = data.frame(row.name = (colnames(countsnew)))
coldata <- data.frame(mycols, condition = factor(c(rep("combined",3), rep("control",3))))
coldata

#check if row names and column names from sample and count matrix match
all(rownames(coldata) %in% colnames(countsnew))
all(rownames(coldata) %in% colnames(countsnew))


library(DESeq2)
library(gplots)
library(ggplot2)


#build the data frame to be used by DESeq2 package
dds=DESeqDataSetFromMatrix(countData = countsnew,colData = coldata,design = ~ condition)
dds


#remove rows with zero.
dds <- dds[rowSums(counts(dds)) > 10,]
dds <-DESeq(dds)


#define vst table
vst <- vst(dds, blind=FALSE)


#plot PCA plot
plotPCA(vst, intgroup="condition", ntop=nrow(counts(dds)))


#plot box plot
a <- DESeq2::plotPCA(vst, intgroup="condition")
a + geom_label(aes(label = coldata$condition),)
nudge <- position_nudge(y = 1)
a + geom_label(aes(label = coldata$condition), position = nudge)
a + geom_text(aes(label = coldata$condition), position = nudge, size=3)
boxplot(assay(vst), col=c("red", "red", "red", "green", "green", "green"), pch=".", vertical=TRUE, cex.axis=0.5, main = "Boxplot of heat and drought affected rice using vst method",
                         las=2, ylab="assay(vst)", xlab="samples", ylim=c(-10,30),
                         font.main= 5, font.axis=0.5, font.lab=2)




#plot correlation heatmap
cU <-cor(as.matrix(assay(vst)))
cols <- c("dodgerblue3", "firebrick3")[coldata$condition]
heatmap.2(cU, symm=TRUE, col= colorRampPalette(c("darkblue","white"))(100),
          labCol=colnames(cU),labRow=colnames(cU),
          distfun=function(c) as.dist(1 - c),
          trace="none",
          colv=TRUE, cexRow=0.9, cexCol=0.9, key=F,
          font=2,
          RowSideColors=cols, ColSideColors=cols)


#plot dispersion plot
plotDispEsts(dds)



#define filename and condition representation for output file generation
res <- results(dds, contrast=c("condition", "combined", "control"))
summary(res)
grp.mean <- sapply(levels(dds$condition),
                   function(lvl)
                     rowMeans(counts(dds,normalized=TRUE)[,dds$condition== lvl]))

norm.counts <- counts(dds,normalized=TRUE)
all <- data.frame(res,assay(vst))
nrow(all)
write.table(all, file="rice_combined_stree_main.csv",sep=",")



#sorting
q_value_sorting<- read.csv(file.choose())      #choose file "rice_combined_stree_main.csv"
library(dplyr)
#using filter to sort data.
q_value_sorting %>% filter(padj<=0.05)->sorted_q_values
nrow(sorted_q_values)

#arranging the table according to ascending order of padj values
sorted_q_values %>% arrange(padj)->ordered_values
top_n(ordered_values, 1000, padj)-> top_thousand_values

#creating a csv file of top differentially expressed gene
write.table(top_thousand_values, file="top_differentially_expressed_gene.csv",sep=",")


#finding top twenty downregulated gene
log2fold_sorting<- read.csv(file.choose())      # choose file "top_differentially_expressed_gene.csv"
log2fold_sorting %>% filter(log2FoldChange<=0)->downregulated_gene
downregulated_gene %>% arrange(log2FoldChange)->top_downregulated_gene
top_n(top_downregulated_gene, 20, log2FoldChange)->twenty_gene

#creating csv file
write.table(twenty_gene, file="top_twenty_downregulated_gene.csv", sep=",")

#finding to twenty downregulated gene
log2fold_sorting_1 <- read.csv(file.choose())      # choose file "top_differentially_expressed_gene.csv
log2fold_sorting_1 %>% filter(log2FoldChange>=0)->upregulated_gene
upregulated_gene %>% arrange(desc(log2FoldChange))->top_upregulated_gene
top_n(top_upregulated_gene, 20, log2FoldChange)->twenty_topgene

#creating csv file
write.table(twenty_topgene, file="top_twenty_upregulated_gene.csv", sep=",")





