
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


library(dplyr)
cts <- (read.csv("D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/Salmon_TPM.csv",sep=",",row.names="GeneSymbol"))
coldata <- read.csv("D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/meta_test.csv", row.names=1)
coldata_not_factorial <- read.csv("D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/meta_test.csv")


data<-data.frame(row.names=rownames(cts))
num_of_exp =6
for (i in 1:num_of_exp){
  temp.df<-data.frame()
  temp.df<-cts %>%
    select(starts_with(coldata_not_factorial[i,1]))
  
  
  data[coldata_not_factorial[i,1]]<- rowMeans(temp.df)
}


coldata$condition <- factor(coldata$condition)


#didn't need to subset just change in the meta data
#data<-subset(data,select=-c(N831,N832))
#coldata<-subset(coldata,select=-c(N831,N832))


#rownames(coldata) <- sub("L00", "", rownames(coldata))
all(rownames(coldata) %in% colnames(data))
all(rownames(coldata) == colnames(data))
data <- data[, rownames(coldata)]
all(rownames(coldata) == colnames(data))



library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = coldata,
                              design = ~ condition)
dds

featureData <- data.frame(gene=rownames(data))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
#removing the outlier
library(dplyr)

#removing low count genes


#check this is put the condition right
#dds$condition <- factor(dds$condition, levels = c("males","females"))

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)

#keep <- rowSums(counts(dds)) > 0
#dds <- dds[keep,]

vsdata <-vst(dds,blind=FALSE)
#autoplotly::autoplotly((plotPCA(vsdata, intgroup="condition")+geom_text(aes(label=name),vjust=2)+labs(title = "PCA: male vs. female"))
#                      , data = vsdata, colour = 'Species', frame = TRUE)

autoplotly::autoplotly((plotPCA(vsdata, intgroup="condition")+geom_text(aes(label=name),vjust=2)+labs(title = "PCA:  male vs. female"))
                       , data = vsdata, colour = 'Species', frame = TRUE)
#create loadings arrows
num<-ncol(data)
pca<-prcomp(data[1:num], scale. = TRUE)  

fviz_pca_var(pca, col.var="contrib", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=TRUE)

#results check
#how many genes are there:
length(res$pvalue)
#how many of them are smaller than 0.05?
sum(res$pvalue <= 0.05, na.rm=TRUE)

#making volcano plots. blue if pval<0.05, red if log2FC>=1 and pval<0.05)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", sub="obp69", xlim=c(-3,3), col.sub="orchid4"))
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#creating a csv file with all the values with pval<.05 & log2FoldChange>=1 from T VS NT for up regulated genes
x<-subset(res, pvalue<.05 & log2FoldChange>=1)
x<-as.data.frame(x)
write.csv(x, "D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/upregulated.csv",
          row.names = TRUE)

#creating a csv file with all the values with pval<.05 & log2FoldChange>=1 from T VS NT for down regulated genes
x<-subset(res, pvalue<.05 & log2FoldChange<=-1)
x<-as.data.frame(x)
write.csv(x, "D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/downregulated.csv",
          row.names = TRUE)



#create dataframes & heatmaps of significant results 
#finding the index of the significant genes in cts
y<-subset(res, pvalue<.05 & abs(log2FoldChange)>=1)
idx<-which(rownames(vsdata[,1]) %in% row.names(y))
#finding the genes
significant<-data[idx,]
#organize the data frame
#row.names(significant)<-significant[,1] 
#significant<-significant[,2:ncol(significant)]
#colnames(significant)<-paste(coldata$mouse, coldata$treat, sep = ", ")
#create a matrix of the genes expression and their names

mat<-as.matrix(significant)

#create heatmap of most changed genes

df <- as.data.frame(colData(dds)[,"condition"])
colnames(df)<-"condition"
row.names(df)<-colnames(mat)

new_order<-data.frame()
new_order<-df[order(df$condition),,drop=FALSE]


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
library("dplyr")

#all that their sum of row is bigger than 50
#df_num_scale<-filter(significant,rowSums(significant)>50)
#df_num_scale = scale(df_num_scale)


#removed all the gene that have below low gene
data_without_zeros_sum<-filter(data,rowSums(data)>0)
df_num_scale<-scale(data)
df_num_scale_dafaframe<-as.data.frame(scale(data_without_zeros_sum))

sub_samp_ordered <- df_num_scale_dafaframe[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="all genes without zero")

#did scaling and removed the rows that have low expression of gene

library("dplyr")
df_num_scale_sum_above1<-filter(df_num_scale_dafaframe,rowSums(df_num_scale_dafaframe)>1)
df_num_scale_sum_below1<-filter(df_num_scale_dafaframe,rowSums(df_num_scale_dafaframe)<1)
sub_samp_ordered <- df_num_scale_sum_below1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="all genes until 1 without zero")
sub_samp_ordered <- df_num_scale_sum_above1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="all genes above 1 without zero")


####after removeiing all the gene that have low from 1 sum after scaling
dataonlybelow1rows <- data[intersect(rownames(data), rownames(df_num_scale_sum_below1)),]

dataonlybelow1rows_scaled <-as.data.frame(scale(dataonlybelow1rows))

sub_samp_ordered <- dataonlybelow1rows_scaled[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main = "without zero and without 1 ")




withoutzero_scale_sum_above1<-filter(dataonlybelow1rows_scaled,rowSums(dataonlybelow1rows_scaled)>1)
withoutzero_scale_sum_below1<-filter(dataonlybelow1rows_scaled,rowSums(dataonlybelow1rows_scaled)<1)
sub_samp_ordered <- withoutzero_scale_sum_below1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="all genes until 1 without zero after removed the values sumrow of 1")
sub_samp_ordered <- withoutzero_scale_sum_above1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="all genes above 1 without zero after removed the values sumrow of 1")




#for significant genes

#scale and remove zeros befor doing scaleing
data_without_zeros_sum<-filter(significant,rowSums(significant)>0)
df_num_scale_dafaframe<-as.data.frame(scale(data_without_zeros_sum))
#show after remove 0
sub_samp_ordered <- df_num_scale_dafaframe[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="significant without zero")

#show the gene that expressed the lowest
library("dplyr")
df_num_scale_sum_above1<-filter(df_num_scale_dafaframe,rowSums(df_num_scale_dafaframe)>1)
df_num_scale_sum_below1<-filter(df_num_scale_dafaframe,rowSums(df_num_scale_dafaframe)<1)
sub_samp_ordered <- df_num_scale_sum_below1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="significant until 1 without zero")
sub_samp_ordered <- df_num_scale_sum_above1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="significant above 1 without zero")

#keep only the lowes genes and return the to the origiinal data befor scale
####after removeiing all the gene that have low from 1 sum after scaling
dataonlybelow1rows <- data[intersect(rownames(significant), rownames(df_num_scale_sum_below1)),]
#doiing scale
dataonlybelow1rows_scaled <-as.data.frame(scale(dataonlybelow1rows))

sub_samp_ordered <- dataonlybelow1rows_scaled[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main = "significant without zero and without 1 ")




withoutzero_scale_sum_above1<-filter(dataonlybelow1rows_scaled,rowSums(dataonlybelow1rows_scaled)>1)
withoutzero_scale_sum_below1<-filter(dataonlybelow1rows_scaled,rowSums(dataonlybelow1rows_scaled)<1)
sub_samp_ordered <- withoutzero_scale_sum_below1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="significant until 1 without zero after removed the values sumrow of 1")
sub_samp_ordered <- withoutzero_scale_sum_above1[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="significant above 1 without zero after removed the values sumrow of 1")



##this part sperated between female and males 
sub_samp_ordered <- df_num_scale[,rownames(new_order) ]

pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order, cluster_cols = F,show_rownames=FALSE,main ="most significant")