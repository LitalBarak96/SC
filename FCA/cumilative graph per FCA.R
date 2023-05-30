library(spatstat.utils)
library(loomR)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)



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



SERUAT_HEAD <-readRDS("D:/seq_data/FCA/seuruat_head5_12.rds")
 
SERUAT_HEAD$sample <- rownames(SERUAT_HEAD@meta.data)

SERUAT_HEAD$sample<-str_replace(SERUAT_HEAD$sample, "__", "_")

SERUAT_HEAD$sample<-str_replace(SERUAT_HEAD$sample, "-", "_")

View(SERUAT_HEAD@meta.data)
SERUAT_HEAD@meta.data <- separate(SERUAT_HEAD@meta.data, col = 'sample', into = c("test",'Barcode','identity', 'Sex'), 
                                  sep = '_')

new<-SplitObject(SERUAT_HEAD, split.by = "identity")


all.genes_head<-as.data.frame(rownames(SERUAT_HEAD))

head_FCA_sum<-function(gene_name,y_value,bottom_x_value,top_x_value,FCA,NUMBER_OF_FCA,path_of_photos){
  
  #remmber to rename the path to the username of he computer
  plotname<- paste(path_of_photos,"/",NUMBER_OF_FCA,".tiff",sep = "")
  #plotname<-paste("C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/27.2.23/",gene_name,"/",NUMBER_OF_FCA,".tiff",sep = "")
  precnt_name<-paste("percent.",gene_name,sep="")
  FCA[[precnt_name]] <- PercentageFeatureSet(FCA,features =gene_name)
  
  gene_head<-as.data.frame(get(precnt_name,FCA@meta.data))
  colnames(gene_head)<-"head"
  
  
  #gene_head_density <- hist(gene_head$head, breaks=seq(0,10,l=1000),plot = FALSE)         # Store output of hist function
  #gene_head_density$density <- gene_head_density$counts /    # Compute density values
  # sum(gene_head_density$counts) * 100
  
  
  gene_sum_head <- hist(gene_head$head, breaks=seq(0,10,l=10000),plot = FALSE)                       # Store histogram info
  gene_sum_head$counts <- revcumsum((gene_sum_head$counts)) 
  
  
  main_title <-paste("cumilative garph in ",NUMBER_OF_FCA,"need to be " ,y_value," cells of ",gene_name)
  tiff(file=plotname,units="in",height=13,width=18,res=600)
  new_yval<-y_value + 0.5*(y_value)
  
  
  exp_per_counts<-gene_sum_head[["breaks"]][max(which(gene_sum_head[["counts"]] > y_value))]
  
  t<-plot(gene_sum_head,col=alpha('red',0.5), main=main_title, xlab='% expression',ylab='number of cells',ylim =c(0,new_yval),xlim =rev(c(bottom_x_value, top_x_value)),cex.axis=2.2,cex.lab=1.5, cex.main=2,xaxt = "n")
  abline(h = y_value, col = "darkgreen")
  abline(v =exp_per_counts, col = "blue") 
  
  axis(1, at = rev(seq(bottom_x_value,top_x_value, by = 0.01)),
       labels = rev(seq(bottom_x_value,top_x_value, by = 0.01)))
  
  
  dev.off()
  
}



#max(which(gene_sum_head[["counts"]] > 4))

number_of_fca<-c("1","4","5","6","7","8","9","10","11","12","13","14","15")
names_of_fca<-paste("FCA",number_of_fca,sep="")

path_of_photos<-"C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/27.2.23/"

gene<-"NPFR"
expceted_cell_num<-100


full_path_of_gene<-paste(path_of_photos,gene,"/",sep = "")

  
if(dir.exists(full_path_of_gene)){
  for(i in 1:length(new)){
    y_Vale_max_origin<-round((expceted_cell_num*ncol(new[[i]]))/100000) # to find how many sample suppose to be
    #y_Vale_max<-y_Vale_max_origin + 0.5 *(y_Vale_max_origin)
    head_FCA_sum(gene,y_Vale_max_origin,0,0.5,new[[i]],names_of_fca[i],full_path_of_gene)
  }
}else{
  dir.create(full_path_of_gene)
  for(i in 1:length(new)){
    y_Vale_max_origin<-round((expceted_cell_num*ncol(new[[i]]))/100000) # to find how many sample suppose to be
    #y_Vale_max<-y_Vale_max_origin + 0.5 *(y_Vale_max_origin)
    head_FCA_sum(gene,y_Vale_max_origin,0,0.3,new[[i]],names_of_fca[i],full_path_of_gene)
    
  }
}




path_of_photos<-"C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/27.2.23/"

gene<-"Dh44"

full_path_of_gene<-paste(path_of_photos,gene,"/",sep = "")

#y_Vale_max_origin<-round((expceted_cell_num*ncol(SERUAT_HEAD))/100000) # to find how many sample suppose to be

y_Vale_max_origin<-6

head_FCA_sum(gene,y_Vale_max_origin,0,0.5,SERUAT_HEAD,"ALL",full_path_of_gene)



#fru 200
#Tdc2 35
#NPF 5
#NPFR 13