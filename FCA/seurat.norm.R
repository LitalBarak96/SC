library(loomR)
library("ggplot2")
library("dplyr")
library("plotly")
setwd("~/Documents_general/Cloud/Bar-Ilan/Lab_files/Projects/Single_cell/FCA_Folder")

df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")
df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)

# let's read the anthena 
# list_files <- list(c("antenna","fca_biohub_antenna_ss2.loom"),c("body_wall","fca_biohub_body_wall_ss2.loom"),
#                    c("fat_body","fca_biohub_fat_body_ss2.loom"),c("gut","fca_biohub_gut_ss2.loom"),
#                    c("haltere","fca_biohub_haltere_ss2.loom"),c("heart","fca_biohub_heart_ss2.loom"),
#                    c("leg","fca_biohub_leg_ss2.loom"),c("male_reproductive_glands","fca_biohub_male_reproductive_glands_ss2.loom"),
#                    c("malpighian_tubule","fca_biohub_malpighian_tubule_ss2.loom"),c("oenocyte","fca_biohub_oenocyte_ss2.loom"),
#                    c("ovary","fca_biohub_ovary_ss2.loom"),c("proboscis_and_maxillary_palps","fca_biohub_proboscis_and_maxillary_palps_ss2.loom"),
#                    c("testis","fca_biohub_testis_ss2.loom"),c("trachea","fca_biohub_trachea_ss2.loom"),
#                    c("wing","fca_biohub_wing_ss2.loom")      

# loom_anthena <- connect(filename = "tissues/fca_biohub_antenna_ss2.loom", mode = "r",skip.validate = TRUE)
loom_anthena <- connect(filename = "tissues/fca_biohub_haltere_ss2.loom", mode = "r",skip.validate = TRUE)


df_tSNE <- as.data.frame(loom_anthena[["col_attrs/Embedding"]][])
df_tSNE$sex <- loom_anthena[["col_attrs/sex"]][]
df_tSNE$ID <- loom_anthena[["col_attrs/CellID"]][]
df_tSNE$attributes <- loom_anthena[["col_attrs/transf_annotation"]][]
names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
sum(is.na(df_tSNE.comb$run_accession))

matrix_data <- as.data.frame(loom_anthena[["matrix"]][,])
rownames(matrix_data) <- loom_anthena[["col_attrs/CellID"]][]
colnames(matrix_data) <- loom_anthena[["row_attrs/Gene"]][]

# each row is one cell. so... let's calculate how much/ if any have less than 200 genes (their cutof)
list_results <- c()
for (row in seq(1,nrow(matrix_data))) {
  list_results[row] <- sum(matrix_data[row,] > 0)
}


# sum(matrix_data[290,] > 0)
# colnames(matrix_data)
# 
hist(matrix_data[matrix_data[,"Adar"] > 0,"Adar"])


#strange.... there are low amount of genes, but they have value in tSNE
#let's check what happens in basal seuran normalization
library(dplyr)
library(Seurat)
library(patchwork)

fca_anthena_seurat <- CreateSeuratObject(counts = t(matrix_data), project = "FCA_anthena", min.cells = 0, min.features = 0)

fca_anthena_seurat[["percent.mt"]] <- PercentageFeatureSet(fca_anthena_seurat, pattern = "^MT-")
head(fca_anthena_seurat@meta.data, 5)
VlnPlot(fca_anthena_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# plot1 <- FeatureScatter(fca_anthena_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fca_anthena_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
plot2

#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

fca_anthena_seurat.norm <- NormalizeData(fca_anthena_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
head(fca_anthena_seurat.norm@meta.data, 5)


fca_anthena_seurat.norm <- FindVariableFeatures(fca_anthena_seurat.norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(fca_anthena_seurat.norm), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(fca_anthena_seurat.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

GetAssayData(object = fca_anthena_seurat.norm, slot = "counts")
# fca_anthena_seurat.norm[[]]@scale.data

fca_anthena_seurat.norm[["RNA"]]@data

df_anthena.norm <- as.data.frame(fca_anthena_seurat.norm[["RNA"]]@data)
df_anthena.norm.selected <- df_anthena.norm[c(grep("Adar",row.names(df_anthena.norm)),
                                              grep("elav",row.names(df_anthena.norm)),
                                              grep("Syt1",row.names(df_anthena.norm)),
                                              grep("para",row.names(df_anthena.norm)),
                                              grep("NPFR",row.names(df_anthena.norm)),
                                              grep("Shab",row.names(df_anthena.norm))
                                              ),] 
df_anthena.norm.selected <- as.data.frame(t(df_anthena.norm.selected))
df_anthena.norm.selected$Sample <- row.names(df_anthena.norm.selected)

df_editing_Sites_known <- read.csv("sample_summary_edit_586.csv")
df_index <- read.csv("devided4_editingIndex.csv")
df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")
df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)

df_anthena.norm.selected.comb <- left_join(df_anthena.norm.selected,df_info_table_selected, by=c("Sample" = "sample_title"))
df_anthena.norm.selected.comb.tSNE <- left_join(df_anthena.norm.selected.comb,df_tSNE, by=c("Sample" = "ID"))
df_anthena.norm.selected.comb.tSNE.index <- left_join(df_anthena.norm.selected.comb.tSNE,df_index, by=c("run_accession" = "Sample"))

#genes
df_anthena.norm.selected.comb.tSNE.index %>% select (Adar, elav, Syt1, para, NPFR, Shab, sex, attributes, X, Y, A2GEditingIndex  ) %>%
  tidyr::gather(Gene,norm.expression, -c(sex, attributes, X, Y, A2GEditingIndex )) %>%
  ggplot(aes(x= Gene, y=norm.expression, fill=Gene)) + geom_boxplot()+theme_bw()+facet_wrap(~attributes)

df_anthena.norm.selected.comb.tSNE.index %>% 
  select (Adar, elav, Syt1, para, NPFR, Shab, sex, attributes, X, Y, A2GEditingIndex  ) %>%
  filter(Adar >=0) %>%
  ggplot(aes(x= Adar, y=log(A2GEditingIndex), fill=Adar)) + geom_point()+theme_bw()+ 
  facet_wrap(~attributes) +
  ggpubr::stat_cor()


df_anthena.norm.selected.comb.tSNE.index %>% 
  select (Adar, elav, Syt1, para, NPFR, Shab, sex, attributes, X, Y, A2GEditingIndex  ) %>%
  filter(Adar >0 & elav > 0) %>%
  ggplot(aes(x= Adar, y=A2GEditingIndex, fill=Adar)) + geom_point()+theme_bw()+ 
  facet_wrap(~attributes) +
  ggpubr::stat_cor()

df_anthena.norm.selected.comb.tSNE.index %>% 
  select (Adar, elav, Syt1, para, NPFR, Shab, sex, attributes, X, Y, A2GEditingIndex  ) %>%
  filter(Adar >0 & elav > 0) %>%
  ggplot(aes(x= Adar, y=A2GEditingIndex, fill=Adar)) + geom_point()+theme_bw()+ 
  facet_wrap(~attributes) +
  ggpubr::stat_cor()


df_anthena.norm.selected.comb.tSNE.index %>% 
  select (Adar, elav, Syt1, para, NPFR, Shab, sex, attributes, X, Y, A2GEditingIndex  ) %>%
  filter(Adar >0 & para > 0 & Syt1>0 & elav >0) %>%
  ggplot(aes(x= Adar, y=A2GEditingIndex, fill=Adar)) + geom_point()+theme_bw()+ 
  facet_wrap(~attributes) +
  ggpubr::stat_cor()


# for a second i'll stop seurat workflow and recheck if h5ad contain something
library(SeuratDisk)
# Convert("tissues_h5ad/antenna.h5ad", "antenna.h5seurat",overwrite = T)
# seuratObj <- LoadH5Seurat("antenna.h5seurat",)

library(Seurat)
library(SeuratDisk)
path_to_fca<-"tissues_h5ad/antenna.h5ad"
Convert(path_to_fca,dest = "h5seurat",overwrite = TRUE)
seurat_anndata<-LoadH5Seurat("tissues_h5ad/antenna.h5seurat")
pbmc <- seurat_anndata

anndata::read_h5ad(path_to_fca)
Y