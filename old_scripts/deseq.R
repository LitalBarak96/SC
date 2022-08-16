
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
up<-subset(res, pvalue<.05 & log2FoldChange>=1)
up<-as.data.frame(up)
write.csv(up, "D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/upregulated.csv",
          row.names = TRUE)

#creating a csv file with all the values with pval<.05 & log2FoldChange>=1 from T VS NT for down regulated genes
down<-subset(res, pvalue<.05 & log2FoldChange<=-1)
down<-as.data.frame(down)
write.csv(down, "D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq/downregulated.csv",
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

#here to change on what data we are doing
#if i want only significant write here instead data significant
df_num_scale<-scale(data)

sub_samp_ordered <- df_num_scale[,rownames(new_order) ]
pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order,  cutree_rows = 4,cluster_cols = F,show_rownames=FALSE,main ="all genes")


df_num_scale<-scale(significant)

##this part sperated between female and males 
sub_samp_ordered <- df_num_scale[,rownames(new_order) ]

pheatmap::pheatmap(sub_samp_ordered, annotation_col  = new_order,  cutree_rows = 4,cluster_cols = F,show_rownames=FALSE,main ="most significant")


#making avg for female and male
middle = ncol(sub_samp_ordered)/2
avg_feamle_and_male <-data.frame(female = rowSums(sub_samp_ordered[,1:middle]),male=rowSums(sub_samp_ordered[,(middle+1): ncol(sub_samp_ordered)]))

pheatmap::pheatmap(avg_feamle_and_male, cutree_rows = 4,cluster_cols = F,show_rownames=FALSE,main ="avg male vs female significant")

#pcc
library("ggpubr")
library("dplyr")
avg_feamle_and_male$id<-(row.names(avg_feamle_and_male) %in% row.names(up))
index<-which(row.names(avg_feamle_and_male) %in% row.names(up))
avg_feamle_and_male[index,]$id ="up"
index<-which(row.names(avg_feamle_and_male) %in% row.names(down))
avg_feamle_and_male[index,]$id ="down"
avg_feamle_and_male$id<-as.factor(avg_feamle_and_male$id)


avg_feamle_and_male$name <- rownames(avg_feamle_and_male)

ggscatter(avg_feamle_and_male, x = "female", y = "male", 
          add = "reg.line",
          conf.int = TRUE, color="id", palette = c("#00AFBB", "#E7B800"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "female", ylab = "male",label = "name")


#library(dplyr)
#test<-as.data.frame(rownames(data))
#after avg male and female we want only those who have diffrenece of 10
avg_feamle_and_male<-avg_feamle_and_male %>% filter(abs(female-male) >5)
pheatmap::pheatmap(avg_feamle_and_male, cutree_rows = 4,cluster_cols = F,main ="diffrence of 5 avg male vs female In significant")

#df
#to how many clusters cutree_rows
#fontsize_row=4
#cutree_cols  =2

#heatmaply(df_num_scale, plot_method ="plotly")


#not good tnse
#library(Rtsne)
#tsne_out <- Rtsne(as.matrix(df_num_scale), dims = 2, perplexity = 5)
#plot(tsne_out$Y, asp=1)
#tsne(df_num_scale,colvec=c('gold'))

