
sub_Seted <-readRDS("D:/seq_data/FCA/seuruat_head18_5.rds")

library(Seurat)

genes_splited<-SplitObject(sub_Seted,split.by = "names" )

fru<-genes_splited[[1]]
Tbh<-genes_splited[[2]]
ple<-genes_splited[[3]]
Crz<-genes_splited[[4]]
NPFR<-genes_splited[[5]]
Tdc2<-genes_splited[[6]]
NPF<-genes_splited[[7]]

Idents(object = fru) <- fru@meta.data$identity

Idents(object = Tbh) <- Tbh@meta.data$identity
Idents(object = ple) <- ple@meta.data$identity
Idents(object = Crz) <- Crz@meta.data$identity
Idents(object = NPFR) <- NPFR@meta.data$identity
Idents(object = Tdc2) <- Tdc2@meta.data$identity
Idents(object = NPF) <- NPF@meta.data$identity


avg_fru<-as.data.frame(AverageExpression(object = fru,slot = 'counts' ))

avg_Tbh<-as.data.frame(AverageExpression(object = Tbh,slot = 'counts' ))
avg_ple<-as.data.frame(AverageExpression(object = ple,slot = 'counts' ))
avg_Crz<-as.data.frame(AverageExpression(object = Crz,slot = 'counts' ))
avg_NPFR<-as.data.frame(AverageExpression(object = NPFR,slot = 'counts' ))
avg_Tdc2<-as.data.frame(AverageExpression(object = Tdc2,slot = 'counts' ))
avg_NPF<-as.data.frame(AverageExpression(object = NPF,slot = 'counts' ))


colnames(avg_fru) = gsub("RNA.FCA", "fru_rep", colnames(avg_fru))
colnames(avg_Tbh) = gsub("RNA.FCA", "Tbh_rep", colnames(avg_Tbh))
colnames(avg_ple) = gsub("RNA.FCA", "ple_rep", colnames(avg_ple))
colnames(avg_Crz) = gsub("RNA.FCA", "Crz_rep", colnames(avg_Crz))
colnames(avg_NPFR) = gsub("RNA.FCA", "NPFR_rep", colnames(avg_NPFR))
colnames(avg_Tdc2) = gsub("RNA.FCA", "Tdc2_rep", colnames(avg_Tdc2))
colnames(avg_NPF) = gsub("RNA.FCA", "NPF._rep", colnames(avg_NPF))


meta_data<-NULL

meta_data<-colnames(avg_fru)
meta_data<-as.data.frame(meta_data)
colnames(meta_data)<-"sample"

meta_data$type<-"fru"


temp_<-NULL

temp_<-colnames(avg_NPF)
temp_<-as.data.frame(temp_)
colnames(temp_)<-"sample"
temp_$type<-"NPF"



meta_data<-rbind(meta_data,temp_)

# write.csv(meta_data, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/meta_data.csv", row.names=FALSE)

# write.csv(avg_fru, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/fru_avg_sc.csv", row.names=TRUE)
# write.csv(avg_Tbh, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/Tbh_avg_sc.csv", row.names=TRUE)
# write.csv(avg_ple, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/ple_avg_sc.csv", row.names=TRUE)
# write.csv(avg_Crz, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/Crz_avg_sc.csv", row.names=TRUE)
# write.csv(avg_NPFR, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/NPFR_avg_sc.csv", row.names=TRUE)
# write.csv(avg_Tdc2, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/Tdc2_avg_sc.csv", row.names=TRUE)
# write.csv(avg_NPF, "C:/Users/lital/OneDrive - Bar Ilan University/Lital/weekly presentation/28.5.23/NPF_avg_sc.csv", row.names=TRUE)

library("readxl")

all_genes<-read_excel("D:/seq_data/FCA/all_gene_avg_sc.xlsx",)

rowname_allgenes<-all_genes$...1


all_genes<-all_genes[,-1]

all_genes<-as.data.frame(all_genes)
row.names(all_genes)<-rowname_allgenes
