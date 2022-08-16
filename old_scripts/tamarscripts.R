library(htmltools)
library("DESeq2")
library(ggplot2)
library(autoplotly)
library(factoextra)
library("pheatmap")
library("RColorBrewer")
library(gplots)
library(limma)
library("manhattanly")
library("ggrepel")
library(plotly)
library("heatmaply")

# loading the files
cts <- read.csv("chemo_count_matrix.csv", header = TRUE, sep = ",")
coldata <- read.csv("chemo_metadata.csv", header = TRUE, sep = ",")
# subsetting the problematic rows, if exist
cts <- subset(cts, Gene != "no_feature" & Gene != "ambiguous" & Gene != "too_low_aQual" &
                Gene != "not_aligned" & Gene != "alignment_not_unique")
# subsetting the non relevant data
coldata<-subset(coldata, exp_hours=="two")
ind<-colnames(cts)%in%coldata$ID
ind[1]<-TRUE
cts<-cts[,ind]
#create DESeq Dataset
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ treat + sex + exp + mouse, tidy=TRUE)
dds <- DESeq(dds)
#checking the names in order to analyze
resultsNames(dds)
#analyzing according to the two groups
res <- results(dds, name="treat_T_vs_NT")

#subsetting every gene that has less than 20 counts total
keep <- rowSums(counts(dds)) >= 20
cts<-cts[keep,]
dds <- dds[keep,]

#plots

vsdata <- vst(dds, blind=FALSE)

#PCA

vsdata$mouse<-as.character(vsdata$mouse)
vsdata$exp<-as.character(vsdata$exp)
#interactive PCA plots. need to upgrade: labels in boxes, the ability to reach a sample
#in the file from the plot in order to examine/change things
autoplotly::autoplotly((plotPCA(vsdata, intgroup="treat")+geom_text(aes(label=name),vjust=2)+labs(title = "PCA: treated vs. not treated"))
                       , data = vsdata, colour = 'Species', frame = TRUE)
autoplotly::autoplotly((plotPCA(vsdata, intgroup="sex")+geom_text(aes(label=name),vjust=2)+labs(title = "PCA: sex")
), data = vsdata, colour = 'Species', frame = TRUE)
autoplotly::autoplotly((plotPCA(vsdata, intgroup="exp")+geom_text(aes(label=name),vjust=2)+labs(title = "PCA: number of experiment")
), data = vsdata, colour = 'Species', frame = TRUE)
autoplotly::autoplotly((plotPCA(vsdata, intgroup="mouse")+geom_text(aes(label=name),vjust=2)+labs(title = "PCA: mouse")
), data = vsdata, colour = 'Species', frame = TRUE)

#create loadings arrows
num<-ncol(cts)
pca<-prcomp(cts[2:num], scale. = TRUE)  
fviz_pca_var(pca, col.var="contrib", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=TRUE)

#results check
#how many genes are there:
length(res$pvalue)
#how many of them are smaller than 0.05?
sum(res$pvalue <= 0.05, na.rm=TRUE)

#making volcano plots. blue if pval<0.05, red if log2FC>=1 and pval<0.05)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", sub="the results of T VS. NT in 2 hours", xlim=c(-3,3), col.sub="orchid4"))
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#creating a csv file with all the values with pval<.05 & log2FoldChange>=1 from T VS NT for up regulated genes
x<-subset(res, pvalue<.05 & log2FoldChange>=1)
x<-row.names(as.data.frame(x))
write.csv(x, "C:/Users/labyissachar/Okun Lab Dropbox/Yissachar Lab/Personal folders/Tamar Saad/R files/exercise files/RNASeq chemo/T VS NT upregulated in 2 hours.csv",
          row.names = FALSE)

#creating a csv file with all the values with pval<.05 & log2FoldChange>=1 from T VS NT for down regulated genes
x<-subset(res, pvalue<.05 & log2FoldChange<=-1)
x<-row.names(as.data.frame(x))
write.csv(x, "C:/Users/labyissachar/Okun Lab Dropbox/Yissachar Lab/Personal folders/Tamar Saad/R files/exercise files/RNASeq chemo/T VS NT downregulated in 2 hours.csv",
          row.names = FALSE)

#create dataframes & heatmaps of significant results 
#finding the index of the significant genes in cts
y<-subset(res, pvalue<.05 & abs(log2FoldChange)>=1)
idx<-which(rownames(vsdata[,1]) %in% row.names(y))
#finding the genes
significant<-cts[idx,]
#organize the data frame
row.names(significant)<-significant[,1] 
significant<-significant[,2:ncol(significant)]
colnames(significant)<-paste(coldata$mouse, coldata$treat, sep = ", ")
#create a matrix of the genes expression and their names
mat<-as.matrix(significant)

#create heatmap of most changed genes

df <- as.data.frame(colData(dds)[,"treat"])
colnames(df)<-"treatment"
row.names(df)<-colnames(mat)

pheatmap(mat, scale = "row", annotation_col =df)
heatmaply(mat, scale = "row", plot_method ="plotly")


