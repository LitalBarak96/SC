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

FCA1<-new$FCA1
FCA4<-new$FCA4
FCA5<-new$FCA5
FCA6<-new$FCA6
FCA7<-new$FCA7
FCA8<-new$FCA8
FCA9<-new$FCA9
FCA10<-new$FCA10
FCA11<-new$FCA11
FCA12<-new$FCA12
FCA13<-new$FCA13
FCA14<-new$FCA14
FCA15<-new$FCA15



all.genes_head<-as.data.frame(rownames(SERUAT_HEAD))

head_FCA<-function(gene_name,y_value,x_value,FCA,NUMBER_OF_FCA){
  
  
  
  plotname<-paste("C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/12.02.23/new/",gene_name,"/",NUMBER_OF_FCA,".tiff",sep = "")
  precnt_name<-paste("percent.",gene_name,sep="")
  FCA[[precnt_name]] <- PercentageFeatureSet(FCA,features =gene_name)
  
  gene_head<-as.data.frame(get(precnt_name,FCA@meta.data))
  colnames(gene_head)<-"head"

  
  gene_head_density <- hist(gene_head$head, breaks=seq(0,10,l=1000),plot = FALSE)         # Store output of hist function
  #gene_head_density$density <- gene_head_density$counts /    # Compute density values
   # sum(gene_head_density$counts) * 100
  
  

  
  
  main_title <-paste("head (red) of " ,NUMBER_OF_FCA ,gene_name)
  tiff(file=plotname,units="in",height=11,width=14,res=600)
  
  plot(gene_head_density, freq = TRUE,col=alpha('red',0.5), main=main_title, xlab='%',ylim =c(0,y_value),xlim = c(0,x_value))
  dev.off()
  
}



number_of_fca<-c("1","4","5","6","7","8","9","10","11","12","13","14","15")
names_of_fca<-paste("FCA",number_of_fca,sep="")

for(i in 1:length(new)){
  head_FCA("trh",20,0.73,new[[i]],names_of_fca[i])
  
}


for(i in 1:length(new)){
  head_FCA("fru",400,0.25,new[[i]],names_of_fca[i])
  
}

for(i in 1:length(new)){
  head_FCA("NPFR",200,0.25,new[[i]],names_of_fca[i])
  
}

