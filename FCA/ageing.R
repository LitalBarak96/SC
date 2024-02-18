# orgs <- as.data.frame(read.delim("/home/alu/aluguest/Lital_galitLab/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv", header=TRUE, allowEscapes=FALSE, sep="\t",  quote="", na.strings="", comment.char=""))

#tmp<-as.data.frame(Matrix::readMM('/home/alu/aluguest/Lital_galitLab/matrix.mtx'))

# expression_matrix <- ReadMtx(
#   mtx = "/home/alu/aluguest/Lital_galitLab/matrix.mtx", features = "/home/alu/aluguest/Lital_galitLab/genes.tsv",
#   cells = "/home/alu/aluguest/Lital_galitLab/barcodes.tsv",skip.feature = 
# )

GENES <- as.data.frame(read.delim("C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/pHd/18.2.24/genes.tsv", header=FALSE, allowEscapes=FALSE, sep="\t",  quote="", na.strings="", comment.char=""))
BARCODE <- as.data.frame(read.delim("C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/pHd/18.2.24/barcodes.tsv", header=FALSE, allowEscapes=FALSE, sep="\t",  quote="", na.strings="", comment.char=""))

MATRIX<-as.data.frame(Matrix::readMM('C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/pHd/18.2.24/matrix.mtx'))


rownames(MATRIX) <- GENES$V2
colnames(MATRIX) <- BARCODE$V1


meta <- as.data.frame(read.delim("C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/pHd/18.2.24/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv", header=TRUE, allowEscapes=FALSE, sep="\t",  quote="", na.strings="", comment.char="", row.names = 1))



test<-CreateSeuratObject( MATRIX, project = "SeuratProject", assay = "RNA",  meta.data = meta)



saveRDS(test,"C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/pHd/18.2.24/18_2_againg_not_normelzied.rds")




age_seruat <-readRDS("C:/Users/barakli8/OneDrive - Bar Ilan University/Lital/pHd/18.2.24/18_2_againg_not_normelzied.rds")

age_seruat <- NormalizeData(object = age_seruat)
age_seruat <- FindVariableFeatures(object = age_seruat)
age_seruat <- ScaleData(object = age_seruat)

age_seruat <- RunPCA(age_seruat,npcs=100)

age_seruat <- FindNeighbors(object = age_seruat, dims = 1:82)
age_seruat <- FindClusters(object = age_seruat, resolution = c(2,8))
#age_seruat <- RunUMAP(object = age_seruat, dims = 1:20)

age_seruat <- RunTSNE(object = age_seruat)


View(age_seruat@meta.data)

Idents(test) <- "sex"


#num_of_genes<-length(unique(age_seruat@meta.data$annotation))
#my_cols<-DiscretePalette(num_of_genes, palette = NULL)

Idents(age_seruat)<-'Age'

# listof_number<-0:80
# nums<-as.character(listof_number)
# test_nums<-subset(age_seruat, annotation != "0")
# View(test_nums@meta.data)

DimPlot(age_seruat,reduction = "tsne",split.by = 'sex',raster=FALSE)+guides(color = guide_legend(override.aes = list(size=1), ncol=2) )+theme(legend.title = element_text(size = 5), 
                                                                                                                                              legend.text = element_text(size = 7))
