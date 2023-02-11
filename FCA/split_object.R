SERUAT_HEAD <-readRDS("D:/seq_data/FCA/seuruat_head5_12.rds")

library(edgeR) # BiocManager::install("edgeR")
library(limma)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(Matrix)
library(tidyr)
library(stringr)
library("gridExtra")
library(loomR)

SERUAT_HEAD$sample <- rownames(SERUAT_HEAD@meta.data)

SERUAT_HEAD$sample<-str_replace(SERUAT_HEAD$sample, "__", "_")

SERUAT_HEAD$sample<-str_replace(SERUAT_HEAD$sample, "-", "_")

View(SERUAT_HEAD@meta.data)
SERUAT_HEAD@meta.data <- separate(SERUAT_HEAD@meta.data, col = 'sample', into = c("test",'Barcode','identity', 'Sex'), 
                                    sep = '_')

new<-SplitObject(SERUAT_HEAD, split.by = "identity")

gene_name <-"NPF"
FCA1<-new$FCA1




all.genes_head<-as.data.frame(rownames(SERUAT_HEAD))

head_FCA<-function(gene_name,y_value,FCA){
  
  precnt_name<-paste("percent.",gene_name,sep="")
  FCA[[precnt_name]] <- PercentageFeatureSet(FCA,features =gene_name)
  
  gene_head<-as.data.frame(get(precnt_name,FCA@meta.data))
  colnames(gene_head)<-"head"

  
  gene_head_density <- hist(gene_head$head, breaks=200,plot = FALSE)         # Store output of hist function
  gene_head_density$density <- gene_head_density$counts /    # Compute density values
    sum(gene_head_density$counts) * 100
  
  

  
  
  main_title <-paste("head (red) of FCA1 of ",gene_name)
  
  plot(gene_head_density, freq = TRUE,col=alpha('red',0.5), main=main_title, xlab='%',ylim =c(0,y_value)) 
}

head_FCA("NPF",2,FCA1)

