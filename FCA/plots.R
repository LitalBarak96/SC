library(loomR)
library("ggplot2")
library("dplyr")
library("plotly")
setwd("~/Documents_general/Cloud/Bar-Ilan/Lab_files/Projects/Single_cell/FCA_Folder")

df_info_table <- read.csv("filereport_read_run_PRJEB45993_2_tsv.txt", sep = "\t")
df_info_table_selected <- select(df_info_table, run_accession, library_name, sample_title)


df_info_table_expression <- read.csv("final_tables/info_table_ss2_ver1_allLoom.csv", sep = ",")


df_editing_Sites_known <- read.csv("sample_summary_edit_586.csv")
df_index <- read.csv("devided4_editingIndex.csv")

df_comb  <- left_join(df_info_table_expression,df_index, by=c("run_accession"="Sample"))
df_comb$attributes_summerized <- df_comb$attributes

df_comb$attributes_summerized[grepl("neuron", df_comb$attributes_summerized)] <- "Neuron"
df_comb$attributes_summerized[grepl("nervous", df_comb$attributes_summerized)] <- "Neuron"
df_comb$attributes_summerized[grepl("epithel", df_comb$attributes_summerized)] <- "Epithel"
df_comb$attributes_summerized[grepl("glia", df_comb$attributes_summerized)] <- "Glia"
df_comb$attributes_summerized[grepl("ovary", df_comb$attributes_summerized)] <- "Ovary"
df_comb$attributes_summerized[grepl("sperm", df_comb$attributes_summerized)] <- "Sperm"
df_comb$attributes_summerized[grepl("cyst", df_comb$attributes_summerized)] <- "Cyst"
df_comb$attributes_summerized[grepl("muscle", df_comb$attributes_summerized)] <- "Muscle"


df_comb %>% 
  tidyr::separate(ID,into = ) 
  # filter(attributes_summerized %in% c("Neuron", "Epithel", "Glia", "Ovary", "Sperm", "Cyst", "Muscle")) %>%
  # filter(attributes_summerized %in% c("Neuron")) %>%
  # filter( grepl("olfactory",attributes)) %>%
  
  # filter(Adar.norm >1 & elav.norm > 0 ) %>% 
  # ggplot(aes(x=attributes_summerized , y=A2GEditingIndex )) +geom_boxplot() +geom_jitter()
  ggplot(aes(x= A2GEditingIndex, y=Adar.norm)) + geom_point()+theme_bw()+ggpubr::stat_cor() +facet_wrap(~attributes)
  





  # select (Adar, elav, Syt1, para, NPFR, Shab, sex, attributes, X, Y, A2GEditingIndex  ) %>%
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