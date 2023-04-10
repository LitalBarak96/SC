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

#gene_name <-"NPF"
#FCA1<-new$FCA1

# FCA1<-new$FCA1
# FCA4<-new$FCA4
# FCA5<-new$FCA5
# FCA6<-new$FCA6
# FCA7<-new$FCA7
# FCA8<-new$FCA8
# FCA9<-new$FCA9
# FCA10<-new$FCA10
# FCA11<-new$FCA11
# FCA12<-new$FCA12
# FCA13<-new$FCA13
# FCA14<-new$FCA14
# FCA15<-new$FCA15



all.genes_head<-as.data.frame(rownames(SERUAT_HEAD))

head_FCA<-function(gene_name,y_value,bottom_x_value,top_x_value,FCA,NUMBER_OF_FCA){
  
  #remmber to rename the path to the username of he computer
  
  plotname<-paste("C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/weekly presentation/16.2.23/",gene_name,"/",NUMBER_OF_FCA,".tiff",sep = "")
  precnt_name<-paste("percent.",gene_name,sep="")
  FCA[[precnt_name]] <- PercentageFeatureSet(FCA,features =gene_name)
  
  gene_had<-as.data.frame(get(precnt_name,FCA@meta.data))
  colnames(gene_head)<-"head"

  
  gene_head_density <- hist(gene_head$head, breaks=seq(0,10,l=1000),plot = FALSE)         # Store output of hist function
  #gene_head_density$density <- gene_head_density$counts /    # Compute density values
   # sum(gene_head_density$counts) * 100
  
  

  
  
  main_title <-paste("head (red) of " ,NUMBER_OF_FCA ,gene_name)
  tiff(file=plotname,units="in",height=13,width=18,res=600)
  
  t<-plot(gene_head_density, freq = TRUE,col=alpha('red',0.5), main=main_title, xlab='%',ylab='f',ylim =c(0,y_value),xlim =rev(c(bottom_x_value, top_x_value)),cex.axis=2.2,cex.lab=1.5, cex.main=2)
  dev.off()
  
}



number_of_fca<-c("1","4","5","6","7","8","9","10","11","12","13","14","15")
names_of_fca<-paste("FCA",number_of_fca,sep="")

for(i in 1:length(new)){
  head_FCA("trh",20,0.73,new[[i]],names_of_fca[i])
  
}



for(i in 1:length(new)){
  
  y_Vale_max<-round((2000*ncol(new[[i]]))/100000) # to find how many sample suppose to be
  y_Vale_max
  head_FCA("fru",y_Vale_max,0.58,5,new[[i]],names_of_fca[i])
  
}

for(i in 1:length(new)){
  y_Vale_max<-round((100*ncol(new[[i]]))/100000) # to find how many sample suppose to be
  y_Vale_max
  head_FCA("NPFR",y_Vale_max,0,2,new[[i]],names_of_fca[i])
  
}


for(i in 1:length(new)){
  y_Vale_max<-round((265*ncol(new[[i]]))/100000) # to find how many sample suppose to be
  y_Vale_max
  head_FCA("Tdc2",y_Vale_max,0.25,10,new[[i]],names_of_fca[i])
  
}

