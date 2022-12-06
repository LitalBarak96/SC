# Install devtools from CRAN
install.packages("devtools")
# Use devtools to install hdf5r and loomR from GitHub
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

# Load loomR
library(loomR)

setwd("~/Documents_general/Cloud/Bar-Ilan/Lab_files/Projects/Single_cell/FCA_Folder")

lfile <- connect(filename = "fca_biohub_antenna_ss2.loom", mode = "r+", skip.validate = TRUE)
lfile


lfile[["matrix"]]
lfile[["col_attrs/cell_names"]]
lfile$row.attrs
lfile[["matrix"]][1:5, 1:5]

full.matrix <- lfile$matrix[, ]
dim(x = full.matrix)

dim(lfile[["matrix"]])


lfile[["matrix"]]$dims
gene.names <- lfile[["row_attrs/gene_names"]][]

lfile[["matrix"]]


l6.immune <- connect(filename = "fca_biohub_antenna_ss2.loom", mode = "r",skip.validate = TRUE)
l6.immune
gene_names <-  l6.immune[["row_attrs/Gene"]][]
gene_of_interest <- "NPFR"
gene_of_interest %in% gene_names
l6.seurat <- as.Seurat(l6.immune)



matrix <- l6.immune[["matrix"]][,]



l6.immune[["col_attrs/CellID"]][]
l6.immune[["col_attrs/"]]



l6.immune[["attrs/MetaData"]][] #python dict of all values
l6.immune[["col_attrs/CellID"]][] # cell ID
l6.immune[["col_attrs/sex"]][] # fly sex
l6.immune[["col_attrs/transf_annotation"]][] #cell annotations
l6.immune[["matrix"]][,] # matrix of expression of all genes
l6.immune[["row_attrs/Gene"]][] # All genes, 16310



lfile[["col_attrs/cell_names"]]
matrix_data <- as.data.frame(lfile[["matrix"]][,])
rownames(matrix_data) <- lfile[["col_attrs/CellID"]][]
colnames(matrix_data) <- lfile[["row_attrs/Gene"]][]


matrix_data[1:10,1:10]
rownames(matrix_data)
lfile[["col_attrs/CellID"]][]
lfile[["col_attrs/transf_annotation"]][]
lfile[["col_attrs/sex"]][]
lfile[["col_attrs/Embedding"]][]
matrix_data$annotatin<-lfile[["col_attrs/transf_annotation"]][]

matrix_data_annotaion_table <- as.data.frame(lfile[["col_attrs/Embedding"]][])
matrix_data_annotaion_table$CellID <- lfile[["col_attrs/CellID"]][]
matrix_data_annotaion_table$transf_annotation <- lfile[["col_attrs/transf_annotation"]][]
matrix_data_annotaion_table$sex <- lfile[["col_attrs/sex"]][]



"E"
library("dplyr")
matrix_data <- as.data.frame(l6.immune[["matrix"]][,])
rownames(matrix_data) <- l6.immune[["col_attrs/CellID"]][]
colnames(matrix_data) <- l6.immune[["row_attrs/Gene"]][]
matrix_data %>% 
  filter(elav >0 | NPFR >0)  %>%
  select("elav","NPFR")




library("ggplot2")
df_test <- as.data.frame(l6.immune[["col_attrs/Embedding"]][])
df_test$sex <- l6.immune[["col_attrs/sex"]][]
df_test$ID <- l6.immune[["col_attrs/CellID"]][]
df_test$attributes <- l6.immune[["col_attrs/transf_annotation"]][]

unique(df_test$attributes_shortened)
df_test$attributes_shortened <- df_test$attributes
df_test$attributes_shortened[grep("neuron",df_test$attributes_shortened)] <- "Neuron"
df_test$attributes_shortened[grep("epithelial cell",df_test$attributes_shortened)] <- "epithelia"
df_test$attributes_shortened[grep("adult antenna glial cell",df_test$attributes_shortened)] <- "glia"
df_test$attributes_shortened[grep("muscle cell",df_test$attributes_shortened)] <- "muscle"
unique(df_test$attributes_shortened)




df_editing_Sites_known <- read.csv("sample_summary_edit_586.csv")
df_index <- read.csv("editingIndex_group1_4_withErrors.tsv", sep = "\t")

df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")


df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)
df_test.comb <- left_join(df_test,df_info_table_selected, by=c("ID"="sample_title"))
# df_test.comb.known <- left_join(df_test.comb,df_editing_Sites_known, by=c("run_accession"="X"))
# df_test.comb.known[,9:ncol(df_test.comb.known)] %>%
#   summarise_all(funs(sum(is.na(.))))
df_test.comb.index <- left_join(df_test.comb,df_index, by=c("run_accession"="Sample"))



ggplot(df_test, aes(x = `_X`, y= `_Y`)) + geom_point(aes(color=attributes_shortened),size=3, alpha=0.3) + theme_bw()
df_test %>% filter(attributes_shortened == "Neuron") %>%
ggplot(aes(x = `_X`, y= `_Y`)) + geom_point(aes(color=attributes_shortened),size=3, alpha=0.3) + theme_bw()

df_test.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  filter(as.numeric(A2GEditingIndex) > 5) %>% 
  ggplot(aes(x = `_X`, y= `_Y` )) + geom_point(aes(color=attributes_shortened, 
                                                   alpha=as.numeric(A2GEditingIndex),
                                                   shape = sex),size=3) + theme_bw() 


df_test.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  filter(attributes_shortened == "Neuron") %>%
  filter(as.numeric(A2GEditingIndex) > 0) %>% 
  filter(attributes %in% c("Johnston organ neuron", "olfactory receptor neuron, coeloconics", "adult olfactory receptor neuron unknown type, Orco-") ) %>%
  ggplot(aes(x = `_X`, y= `_Y` )) + geom_point(aes(color=attributes, 
                                                   alpha=as.numeric(A2GEditingIndex),
                                                   shape = sex),size=3) + theme_bw() 





## allright. let's create the matrixes ----

library(loomR)
library(dplyr)
setwd("C:/Users/zbida/Downloads/FCA_Folder/")
df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")
df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)



## anthena ----
loom_anthena <- connect(filename = "tissues/fca_biohub_antenna_ss2.loom", mode = "r",skip.validate = TRUE)

df_tSNE <- as.data.frame(loom_anthena[["col_attrs/Embedding"]][])
df_tSNE$sex <- loom_anthena[["col_attrs/sex"]][]
df_tSNE$ID <- loom_anthena[["col_attrs/CellID"]][]
df_tSNE$attributes <- loom_anthena[["col_attrs/transf_annotation"]][]
names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
sum(is.na(df_tSNE.comb$run_accession))
write.csv(x = df_tSNE.comb, file="final_tables/tSNE_anthena_ss2.csv", quote = F, row.names = F)

# df_info_table <- bind_rows(df_info_table, select(df_tSNE.comb, run_accession, ID, attributes, sex))

matrix_data <- as.data.frame(loom_anthena[["matrix"]][,])
rownames(matrix_data) <- loom_anthena[["col_attrs/CellID"]][]
colnames(matrix_data) <- loom_anthena[["row_attrs/Gene"]][]
write.csv(x = matrix_data, file="final_tables/matrixGenes_anthena_ss2.csv", quote = F, 
          row.names = T)




## body wall ----


loom_body_wall <- connect(filename = "tissues/fca_biohub_body_wall_ss2.loom", mode = "r",skip.validate = TRUE)

df_tSNE <- as.data.frame(loom_body_wall[["col_attrs/Embedding"]][])
df_tSNE$sex <- loom_body_wall[["col_attrs/sex"]][]
df_tSNE$ID <- loom_body_wall[["col_attrs/CellID"]][]
df_tSNE$attributes <- loom_body_wall[["col_attrs/transf_annotation"]][]
names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
sum(is.na(df_tSNE.comb$run_accession))
write.csv(x = df_tSNE.comb, file="final_tables/tSNE_body_wall_ss2.csv", quote = F, row.names = F)

# df_info_table <- bind_rows(df_info_table, select(df_tSNE.comb, run_accession, ID, attributes, sex))


matrix_data <- as.data.frame(loom_body_wall[["matrix"]][,])
rownames(matrix_data) <- loom_body_wall[["col_attrs/CellID"]][]
colnames(matrix_data) <- loom_body_wall[["row_attrs/Gene"]][]
write.csv(x = matrix_data, file="final_tables/matrixGenes_body_wall_ss2.csv", quote = F, 
          row.names = T)



## fat_body ----

loom_fatBody <- connect(filename = "tissues/fca_biohub_fat_body_ss2.loom", mode = "r",skip.validate = TRUE)

df_tSNE <- as.data.frame(loom_fatBody[["col_attrs/Embedding"]][])
df_tSNE$sex <- loom_fatBody[["col_attrs/sex"]][]
df_tSNE$ID <- loom_fatBody[["col_attrs/CellID"]][]
df_tSNE$attributes <- loom_fatBody[["col_attrs/transf_annotation"]][]
names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
sum(is.na(df_tSNE.comb$run_accession))
write.csv(x = df_tSNE.comb, file="final_tables/tSNE_fatBody_ss2.csv", quote = F, row.names = F)

# df_info_table <- bind_rows(df_info_table, select(df_tSNE.comb, run_accession, ID, attributes, sex))


matrix_data <- as.data.frame(loom_fatBody[["matrix"]][,])
rownames(matrix_data) <- loom_fatBody[["col_attrs/CellID"]][]
colnames(matrix_data) <- loom_fatBody[["row_attrs/Gene"]][]
write.csv(x = matrix_data, file="final_tables/matrixGenes_fatBody_ss2.csv", quote = F, 
          row.names = T)


## gut ----

loom_gut <- connect(filename = "tissues/fca_biohub_gut_ss2.loom", mode = "r",skip.validate = TRUE)

df_tSNE <- as.data.frame(loom_gut[["col_attrs/Embedding"]][])
df_tSNE$sex <- loom_gut[["col_attrs/sex"]][]
df_tSNE$ID <- loom_gut[["col_attrs/CellID"]][]
df_tSNE$attributes <- loom_gut[["col_attrs/transf_annotation"]][]
names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
sum(is.na(df_tSNE.comb$run_accession))
write.csv(x = df_tSNE.comb, file="final_tables/tSNE_gut_ss2.csv", quote = F, row.names = F)

# df_info_table <- bind_rows(df_info_table, select(df_tSNE.comb, run_accession, ID, attributes, sex))


matrix_data <- as.data.frame(loom_gut[["matrix"]][,])
rownames(matrix_data) <- loom_gut[["col_attrs/CellID"]][]
colnames(matrix_data) <- loom_gut[["row_attrs/Gene"]][]
write.csv(x = matrix_data, file="final_tables/matrixGenes_gut_ss2.csv", quote = F, 
          row.names = T)



#---- general function ----
library(loomR)
library("ggplot2")
library("dplyr")
library("plotly")
library(dplyr)
library(Seurat)
library(patchwork)
setwd("~/Documents_general/Cloud/Bar-Ilan/Lab_files/Projects/Single_cell/FCA_Folder")
df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")


df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)


list_files <- list(c("antenna","fca_biohub_antenna_ss2.loom"),c("body_wall","fca_biohub_body_wall_ss2.loom"),
                   c("fat_body","fca_biohub_fat_body_ss2.loom"),c("gut","fca_biohub_gut_ss2.loom"),
                   c("haltere","fca_biohub_haltere_ss2.loom"),c("heart","fca_biohub_heart_ss2.loom"),
                   c("leg","fca_biohub_leg_ss2.loom"),c("male_reproductive_glands","fca_biohub_male_reproductive_glands_ss2.loom"),
                   c("malpighian_tubule","fca_biohub_malpighian_tubule_ss2.loom"),c("oenocyte","fca_biohub_oenocyte_ss2.loom"),
                   c("ovary","fca_biohub_ovary_ss2.loom"),c("proboscis_and_maxillary_palps","fca_biohub_proboscis_and_maxillary_palps_ss2.loom"),
                   c("testis","fca_biohub_testis_ss2.loom"),c("trachea","fca_biohub_trachea_ss2.loom"),
                   c("wing","fca_biohub_wing_ss2.loom")                   )

df_info_table <- data.frame(run_accession=character(), ID = character(), 
                            attributes=character(), sex = character(),
                            rowCounts.ADAR = numeric(),rowCounts.elav = numeric(), rowCounts.syt1 = numeric(), 
                            rowCounts.para = numeric(),rowCounts.NPFR = numeric(), rowCounts.shab = numeric(),
                            Adar.norm = numeric(),elav.norm = numeric(), Syt1.norm = numeric(), 
                            para.norm = numeric(),NPFR.norm = numeric(), Shab.norm = numeric() 
                          
)

for (item in list_files) {

  loom_file <- connect(filename = paste0("tissues/",item[[2]]), mode = "r",skip.validate = TRUE)
  
  df_tSNE <- as.data.frame(loom_file[["col_attrs/Embedding"]][])
  df_tSNE$sex <- loom_file[["col_attrs/sex"]][]
  df_tSNE$ID <- loom_file[["col_attrs/CellID"]][]
  df_tSNE$attributes <- loom_file[["col_attrs/transf_annotation"]][]
  df_tSNE$attributes <- gsub(",", ";",df_tSNE$attributes)
  names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
  names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
  
  df_tSNE$rowCounts.ADAR <- loom_file[["matrix"]][,grep("Adar", loom_file[["row_attrs/Gene"]][])]
  df_tSNE$rowCounts.elav <- loom_file[["matrix"]][,grep("elav", loom_file[["row_attrs/Gene"]][])]
  df_tSNE$rowCounts.syt1 <- loom_file[["matrix"]][,grep("Syt1$", loom_file[["row_attrs/Gene"]][])]
  df_tSNE$rowCounts.para <- loom_file[["matrix"]][,grep("para$", loom_file[["row_attrs/Gene"]][])]
  df_tSNE$rowCounts.NPFR <- loom_file[["matrix"]][,grep("NPFR", loom_file[["row_attrs/Gene"]][])]
  df_tSNE$rowCounts.shab <- loom_file[["matrix"]][,grep("Shab", loom_file[["row_attrs/Gene"]][])]
  
  #run seurat
  
  #let's check what happens in basal seuran normalization
  library(dplyr)
  library(Seurat)
  library(patchwork)
  
  
  matrix_data <- as.data.frame(loom_file[["matrix"]][,])
  rownames(matrix_data) <- loom_file[["col_attrs/CellID"]][]
  colnames(matrix_data) <- loom_file[["row_attrs/Gene"]][]
  fca_seurat <- CreateSeuratObject(counts = t(matrix_data), project = item[[1]], min.cells = 0, min.features = 0)
  
  fca_seurat[["percent.mt"]] <- PercentageFeatureSet(fca_seurat, pattern = "^MT-")
  
  fca_seurat.norm <- NormalizeData(fca_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  # fca_seurat.norm <- FindVariableFeatures(fca_seurat.norm, selection.method = "vst", nfeatures = 2000)
  # top10 <- head(VariableFeatures(fca_anthena_seurat.norm), 10)
  
  # plot variable features with and without labels

  # GetAssayData(object = fca_seurat.norm, slot = "counts")

  # fca_seurat.norm[["RNA"]]@data
  
  df_seurat.norm <- as.data.frame(fca_seurat.norm[["RNA"]]@data)
  df_seurat.norm.selected <- df_seurat.norm[c(grep("Adar",row.names(df_seurat.norm)),
                                                grep("elav",row.names(df_seurat.norm)),
                                                grep("Syt1",row.names(df_seurat.norm)),
                                                grep("para",row.names(df_anthena.norm)),
                                                grep("NPFR",row.names(df_seurat.norm)),
                                                grep("Shab",row.names(df_seurat.norm))
  ),] 
  df_seurat.norm.selected <- as.data.frame(t(df_seurat.norm.selected))
  df_seurat.norm.selected$Sample <- row.names(df_seurat.norm.selected)
  
  colnames(df_seurat.norm.selected)[colnames(df_seurat.norm.selected) == "Adar"] <- "Adar.norm"
  colnames(df_seurat.norm.selected)[colnames(df_seurat.norm.selected) == "elav"] <- "elav.norm"
  colnames(df_seurat.norm.selected)[colnames(df_seurat.norm.selected) == "Syt1"] <- "Syt1.norm"
  colnames(df_seurat.norm.selected)[colnames(df_seurat.norm.selected) == "para"] <- "para.norm"
  colnames(df_seurat.norm.selected)[colnames(df_seurat.norm.selected) == "NPFR"] <- "NPFR.norm"
  colnames(df_seurat.norm.selected)[colnames(df_seurat.norm.selected) == "Shab"] <- "Shab.norm"

  ## end run seurat
  
  df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
  df_tSNE.comb <- left_join(df_tSNE.comb,df_seurat.norm.selected, by=c("ID"="Sample"))
  
  sum(is.na(df_tSNE.comb$run_accession))
  write.csv(x = df_tSNE.comb, file=paste0("final_tables/tSNE_",item[[1]],"_ss2.csv"), quote = F, row.names = F)
  
  df_info_table <- bind_rows(df_info_table, select(df_tSNE.comb, run_accession, ID, attributes, sex,
                                                   rowCounts.ADAR, rowCounts.elav, rowCounts.syt1, 
                                                   rowCounts.para, rowCounts.NPFR, rowCounts.shab,
                                                   Adar.norm, elav.norm, Syt1.norm, para.norm, NPFR.norm,Shab.norm
                                                   
                                                   ))
  
  matrix_data <- as.data.frame(loom_file[["matrix"]][,])
  rownames(matrix_data) <- loom_file[["col_attrs/CellID"]][]
  colnames(matrix_data) <- loom_file[["row_attrs/Gene"]][]
  write.csv(x = matrix_data, file=paste0("final_tables/matrixGenes_",item[[2]],"_ss2.csv"), quote = F, 
            row.names = T)
  
  print(item[[1]])
  
}
df_info_table$attributes <- gsub(",", ";",df_info_table$attributes)
write.csv(x = df_info_table, file=paste0("final_tables/info_table_ss2_ver1_allLoom.csv"), quote = F, row.names = F)




df_tSNE$ADAR.rowCoumts <- loom_file[["matrix"]][,grep("Adar", loom_file[["row_attrs/Gene"]][])]
df_tSNE$elav.rowCoumts <- loom_file[["matrix"]][,grep("elav", loom_file[["row_attrs/Gene"]][])]
df_tSNE$syt1.rowCoumts <- loom_file[["matrix"]][,grep("Syt1$", loom_file[["row_attrs/Gene"]][])]
df_tSNE$para.rowCoumts <- loom_file[["matrix"]][,grep("para$", loom_file[["row_attrs/Gene"]][])]


#write all gene matrix
library(dplyr)
library(Seurat)
library(patchwork)
for (item in list_files) {

  loom_file <- connect(filename = paste0("tissues/",item[[2]]), mode = "r",skip.validate = TRUE)
  
  # df_tSNE <- as.data.frame(loom_file[["col_attrs/Embedding"]][])
  # df_tSNE$sex <- loom_file[["col_attrs/sex"]][]
  # df_tSNE$ID <- loom_file[["col_attrs/CellID"]][]
  # df_tSNE$attributes <- loom_file[["col_attrs/transf_annotation"]][]
  # df_tSNE$attributes <- gsub(",", ";",df_tSNE$attributes)
  # names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
  # names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"
  # 
  # df_tSNE$rowCounts.ADAR <- loom_file[["matrix"]][,grep("Adar", loom_file[["row_attrs/Gene"]][])]
  # df_tSNE$rowCounts.elav <- loom_file[["matrix"]][,grep("elav", loom_file[["row_attrs/Gene"]][])]
  # df_tSNE$rowCounts.syt1 <- loom_file[["matrix"]][,grep("Syt1$", loom_file[["row_attrs/Gene"]][])]
  # df_tSNE$rowCounts.para <- loom_file[["matrix"]][,grep("para$", loom_file[["row_attrs/Gene"]][])]
  # df_tSNE$rowCounts.NPFR <- loom_file[["matrix"]][,grep("NPFR", loom_file[["row_attrs/Gene"]][])]
  # df_tSNE$rowCounts.shab <- loom_file[["matrix"]][,grep("Shab", loom_file[["row_attrs/Gene"]][])]
  # 
  #run seurat
  
  #let's check what happens in basal seuran normalization
  
  
  
  matrix_data <- as.data.frame(loom_file[["matrix"]][,])
  rownames(matrix_data) <- loom_file[["col_attrs/CellID"]][]
  colnames(matrix_data) <- loom_file[["row_attrs/Gene"]][]
  fca_seurat <- CreateSeuratObject(counts = t(matrix_data), project = item[[1]], min.cells = 0, min.features = 0)
  
  fca_seurat[["percent.mt"]] <- PercentageFeatureSet(fca_seurat, pattern = "^MT-")
  
  fca_seurat.norm <- NormalizeData(fca_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  
  df_seurat.norm <- as.data.frame(fca_seurat.norm[["RNA"]]@data)

  
  write.csv(x = df_seurat.norm, file=paste0("final_tables/allGenes_norm_",item[[1]],"_ss2.csv"), quote = F, row.names = F)
  
  

  
  print(item[[1]])
  
}

# read matrix ----
library(loomR)
library("ggplot2")
library("dplyr")
library("plotly")
setwd("C:/Users/zbida/Downloads/FCA_Folder/")

list_files <- list(c("antenna","fca_biohub_antenna_ss2.loom"),c("body_wall","fca_biohub_body_wall_ss2.loom"),
                   c("fat_body","fca_biohub_fat_body_ss2.loom"),c("gut","fca_biohub_gut_ss2.loom"),
                   c("haltere","fca_biohub_haltere_ss2.loom"),c("heart","fca_biohub_heart_ss2.loom"),
                   c("leg","fca_biohub_leg_ss2.loom"),c("male_reproductive_glands","fca_biohub_male_reproductive_glands_ss2.loom"),
                   c("malpighian_tubule","fca_biohub_malpighian_tubule_ss2.loom"),c("oenocyte","fca_biohub_oenocyte_ss2.loom"),
                   c("ovary","fca_biohub_ovary_ss2.loom"),c("proboscis_and_maxillary_palps","fca_biohub_proboscis_and_maxillary_palps_ss2.loom"),
                   c("testis","fca_biohub_testis_ss2.loom"),c("trachea","fca_biohub_trachea_ss2.loom"),
                   c("wing","fca_biohub_wing_ss2.loom")                   )

df_editing_Sites_known <- read.csv("sample_summary_edit_586.csv")
df_index <- read.csv("devided4_editingIndex.csv")
df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")
df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)


loom_path <- list_files[[1]][2]
loom_name <-  list_files[[1]][1]
print(loom_name)
loom_file <- connect(filename = paste0("tissues/",loom_path), mode = "r",skip.validate = TRUE)

df_tSNE <- as.data.frame(loom_file[["col_attrs/Embedding"]][])
df_tSNE$sex <- loom_file[["col_attrs/sex"]][]
df_tSNE$ID <- loom_file[["col_attrs/CellID"]][]
df_tSNE$attributes <- loom_file[["col_attrs/transf_annotation"]][]
names(df_tSNE)[names(df_tSNE) == "_X"] <- "X"
names(df_tSNE)[names(df_tSNE) == "_Y"] <- "Y"


df_tSNE$ADAR.rowCoumts <- loom_file[["matrix"]][,grep("Adar", loom_file[["row_attrs/Gene"]][])]
df_tSNE$elav.rowCoumts <- loom_file[["matrix"]][,grep("elav", loom_file[["row_attrs/Gene"]][])]
df_tSNE$syt1.rowCoumts <- loom_file[["matrix"]][,grep("Syt1$", loom_file[["row_attrs/Gene"]][])]
df_tSNE$para.rowCoumts <- loom_file[["matrix"]][,grep("para$", loom_file[["row_attrs/Gene"]][])]

df_tSNE.comb <- left_join(df_tSNE,df_info_table_selected, by=c("ID"="sample_title"))
df_tSNE.comb.index <- left_join(df_tSNE.comb,df_index, by=c("run_accession"="Sample"))

df_tSNE.comb.index$attributes_shortened <- df_tSNE.comb.index$attributes
unique(df_tSNE.comb.index$attributes_shortened)
df_tSNE.comb.index$attributes_shortened[grep("neuron",df_tSNE.comb.index$attributes_shortened)] <- "Neuron"
df_tSNE.comb.index$attributes_shortened[grep("epithelial cell",df_tSNE.comb.index$attributes_shortened)] <- "epithelia"
df_tSNE.comb.index$attributes_shortened[grep("epidermal cell that specialized in antimicrobial response",df_tSNE.comb.index$attributes_shortened)] <- "epithelia"
df_tSNE.comb.index$attributes_shortened[grep("adult tracheal cell",df_tSNE.comb.index$attributes_shortened)] <- "epithelia"
df_tSNE.comb.index$attributes_shortened[grep("tendon cell",df_tSNE.comb.index$attributes_shortened)] <- "epithelia"


df_tSNE.comb.index$attributes_shortened[grep("adult antenna glial cell",df_tSNE.comb.index$attributes_shortened)] <- "glia"
df_tSNE.comb.index$attributes_shortened[grep("peripheral glial cell",df_tSNE.comb.index$attributes_shortened)] <- "glia"

df_tSNE.comb.index$attributes_shortened[grep("muscle cell",df_tSNE.comb.index$attributes_shortened)] <- "muscle"
df_tSNE.comb.index$attributes_shortened[grep("adult alary muscle",df_tSNE.comb.index$attributes_shortened)] <- "muscle"
df_tSNE.comb.index$attributes_shortened[grep("indirect flight muscle",df_tSNE.comb.index$attributes_shortened)] <- "muscle"

df_tSNE.comb.index$attributes_shortened[grep("adult ventral nervous system",df_tSNE.comb.index$attributes_shortened)] <- "nervous"
df_tSNE.comb.index$attributes_shortened[grep("adult peripheral nervous system",df_tSNE.comb.index$attributes_shortened)] <- "nervous"

unique(df_tSNE.comb.index$attributes_shortened)


#all

ggplot(df_tSNE.comb.index, aes(x = X, y= Y)) + 
  geom_point(aes(color=attributes),size=3, alpha=0.5) + theme_bw() +
  ggtitle(loom_name)# +  theme(legend.position="bottom")

table(df_tSNE.comb.index$attributes)[table(df_tSNE.comb.index$attributes)>10]


ggplot(df_tSNE.comb.index, aes(x = X, y= Y)) + 
  geom_point(aes(color=attributes_shortened),size=3, alpha=0.5) + 
  theme_bw() + ggtitle(loom_name) +  theme(legend.position="bottom") # +scale_color_brewer(palette="Dark2")

table(df_tSNE.comb.index$attributes_shortened)[table(df_tSNE.comb.index$attributes_shortened)>10]

#neurons
# df_tSNE.comb.index %>% filter(attributes_shortened == "Neuron") %>%
#   ggplot(aes(x = X, y= Y)) + 
#   geom_point(aes(color=attributes),size=3, alpha=0.3) + theme_bw()

#index
df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  mutate(A2G = as.numeric(A2GEditingIndex)) %>% 
  # filter(A2G > 0) %>% 
  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=A2G, shape = sex),size=3) + theme_bw() + scale_colour_gradient2()


# df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
#   mutate(A2G = as.numeric(A2GEditingIndex)) %>% 
#   filter(ADAR.rowCoumts > 0) %>%
#   ggplot(aes(x = X, y= Y )) + geom_point(aes(color=A2G, shape = sex),size=3) + theme_bw() + scale_colour_gradient2()


df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  mutate(A2G = as.numeric(A2GEditingIndex)) %>% 
  filter(ADAR.rowCoumts > 0) %>%
  # filter(A2G > 0) %>% 
  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=attributes_shortened, 
                                             alpha=A2G,
                                             shape = sex),size=5) + theme_bw() 



df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>%
  mutate(A2G = as.numeric(A2GEditingIndex)) %>%
  filter(ADAR.rowCoumts > 0) %>%
  # filter(A2G < 30) %>%
  ggplot(aes(x = attributes, y= A2G, fill=attributes_shortened )) +
  geom_boxplot() + geom_jitter(size=3, alpha=0.4)+ theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>%
  filter(elav.rowCoumts > 0) %>%
ggplot(aes(x = X, y= Y)) + 
  geom_point(aes(color=attributes_shortened),size=3, alpha=0.5) + 
  theme_bw() + ggtitle(loom_name) +  theme(legend.position="bottom") # +scale_color_brewer(palette="Dark2")

df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>%
  filter(elav.rowCoumts > 0) %>% 
  filter(ADAR.rowCoumts > 0) %>% 
  select(attributes, elav.rowCoumts,para.rowCoumts,A2GEditingIndex) %>% 
  group_by(attributes) %>% 
  summarise(min=min(elav.rowCoumts), max=max(elav.rowCoumts),
            mean = mean(elav.rowCoumts), median=median(elav.rowCoumts), count=length(elav.rowCoumts),
            A2G=median(A2GEditingIndex))


df_test.comb.known <- left_join(df_tSNE.comb.index,df_editing_Sites_known, by=c("run_accession"="X"))
minimal_cell_amount = 10
for (attribute in unique(df_test.comb.known$attributes_shortened)) {
  total_cells <- sum(df_test.comb.known$attributes_shortened == attribute & df_test.comb.known$ADAR.rowCoumts > 0)
  if(total_cells < minimal_cell_amount){next}
  print(attribute)
  
  sites_with_edit <- df_test.comb.known[df_test.comb.known$attributes_shortened == attribute & df_test.comb.known$ADAR.rowCoumts > 0,
                                        grep("chr",colnames(df_test.comb.known))] %>%
    summarise_all(funs(total_cells - sum(is.na(.)))) %>% 
    as.numeric() 
  sites_over_1 <- sum(sites_with_edit/total_cells >0.01)
  sites_over_5 <- sum(sites_with_edit/total_cells >0.05)
  sites_over_10 <- sum(sites_with_edit/total_cells >0.1)
  sites_over_20 <- sum(sites_with_edit/total_cells >0.2)
  sites_over_30 <- sum(sites_with_edit/total_cells >0.3)
  print(paste("Total", total_cells, "and amount cells with editing site - ",
              "1:",sites_over_1,"5:",sites_over_5,"10:",sites_over_10,"20:",sites_over_20,"30:",sites_over_30))
  
  
  df_test2 <- df_test.comb.known[df_test.comb.known$attributes == attribute,
                                 18:ncol(df_test.comb.known)]
  
  # ncol(df_test2)== length(sites_with_edit)
  # df_test2[,sites_with_edit>0]
}


  










 df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  mutate(A2G = as.numeric(A2GEditingIndex)) %>% 
  # filter(A2G > 0) %>% 
  ggplot(aes(x = attributes_shortened, y= A2G, fill=attributes_shortened )) + 
  geom_violin()+geom_boxplot() + theme_bw() + scale_fill_brewer(palette="Dark2") 



df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  mutate(A2G = as.numeric(A2GEditingIndex)) %>% 
  # filter(A2G > 0) %>% 
  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=attributes_shortened, 
                                             alpha=A2G,
                                             shape = sex),size=3) + theme_bw() 

df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  filter(as.numeric(A2GEditingIndex) > 0) %>% 
  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=attributes_shortened, 
                                                   alpha=as.numeric(A2GEditingIndex),
                                                   shape = sex),size=3) + theme_bw() 



df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>% 
  mutate(A2G = as.numeric(A2GEditingIndex)) %>%
  # filter(attributes_shortened == "Neuron") %>%
  mutate("A-G group" = if_else(A2G > 0, if_else(A2G > 10,"High: >10", "Low: <- 10"), "0")) %>%
  mutate(`A-G group`=factor(`A-G group`, levels = c("0","Low: <- 10","High: >10"))) %>%
  # filter(as.numeric(A2GEditingIndex) > 0) %>% 

  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=attributes_shortened, 
                                             alpha=`A-G group`,
                                             shape = sex),size=3) + theme_bw() 


df_tSNE.comb.index %>% filter(! is.na(A2GEditingIndex)) %>%
  mutate(A2G = as.numeric(A2GEditingIndex)) -> df_test
df_test$A2G[df_test$A2G >10] <- 20

df_test %>%
  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=A2G,
                                             shape = sex),size=3) + theme_bw() +
  scale_colour_gradient2() 

# df_test %>% filter(A2G >5) %>%
#   ggplot(aes(x = X, y= Y )) + geom_point(aes(color=A2G,
#                                              shape = sex),size=3) + theme_bw() +
#   scale_colour_gradient2() 

df_test %>%
  ggplot(aes(x = X, y= Y )) + geom_point(aes(color=A2G,
                                             shape = sex),size=3) + theme_bw() +
  scale_colour_gradient2() +facet_wrap(~attributes)




#All index files -----

df_info <- read.csv("final_tables/info_table_ss2_ver1_allLoom.csv")
df_info$tissue <- tidyr::separate(df_info,ID,c("FCA","P","Sex","tissue","B")) %>% .$tissue
df_index <- read.csv("devided4_editingIndex.csv")
df_comb <- left_join(df_info,df_index,by=c("run_accession"="Sample"))

sort(unique(df_comb$attributes))
table(df_comb$attributes_short,df_comb$tissue)
df_comb$attributes_short <- NA
df_comb$attributes_short[grep("neuron",df_comb$attributes)] <- "Neuron"
df_comb$attributes_short[grep("glia",df_comb$attributes)] <- "Glia"
df_comb$attributes_short[grep("muscle",df_comb$attributes)] <- "Muscle"
df_comb$attributes_short[grep("epithel",df_comb$attributes)] <- "Epithel"


df_comb %>% filter(! is.na(A2GEditingIndex)) %>%
  mutate(A2G = as.numeric(A2GEditingIndex)) %>%
  filter(!is.na(attributes_short)) %>%
  ggplot(aes(x = attributes_short, y= A2G, fill=attributes_short )) +
  geom_boxplot() + geom_jitter(size=3, alpha=0.4)+ theme_bw() +facet_wrap(~tissue) +ylim(0,40)
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df_comb %>% filter(! is.na(A2GEditingIndex)) %>%
  mutate(A2G = as.numeric(A2GEditingIndex)) %>%
  filter(!is.na(attributes_short)) %>%
  filter(attributes_short == "Neuron" ) %>%
  ggplot(aes(x = tissue, y= A2G, fill=tissue )) +
  geom_boxplot() + geom_jitter(size=3, alpha=0.4)+ theme_bw() 

df_comb %>% filter(! is.na(A2GEditingIndex)) %>%
  mutate(A2G = as.numeric(A2GEditingIndex)) %>%
  filter(rowCounts.elav > 0  & attributes_short != "Neuron") %>%
  ggplot(aes(x = tissue, y= A2G, fill=tissue )) +
  geom_boxplot() + geom_jitter(size=3, alpha=0.4)+ theme_bw() 


#read h5ad ----

# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_github(repo = "mojaveazure/seurat-disk")
# remotes::install_github('mojaveazure/seurat-disk', ref = 'feat/loom')
# remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
# remotes::install_github('mojaveazure/seurat-disk')

# Convert("akhCC.h5ad", dest = "h5seurat", overwrite = TRUE)
# pbmc3k <- LoadH5Seurat("akhCC.h5seurat",check.names=FALSE)
# 
# 
# library(SeuratDisk)
# Convert("tissues_h5ad/akhCC.h5ad", ".h5seurat")
# LoadH5Seurat("tissues_h5ad/akhCC.h5seurat")
# 
# # 
# library(SeuratDisk)
# Convert("cd10pos.h5ad", ".h5seurat")
# seuratObj <- LoadH5Seurat("cd10pos.h5Seurat")
# 
library(rhdf5)
h5ls("tissues_h5ad/akhCC.h5ad")
file <- h5file("tissues_h5ad/akhCC.h5ad")

file[["var/_index"]] #gene
file[["obs/_index"]][] #cols
file[["obs/sex"]] #sex



if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("zellkonverter")
library("zellkonverter")

single_cell_file <- zellkonverter::readH5AD("akhCC.h5ad")
single_cell_file@colData
names(single_cell_file@rowRanges[])
