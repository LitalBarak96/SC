#this script take the data after salmon and make dataframe with the group the gene and comoon name with MTP
library(dplyr)
setwd('D:/RNA_seq/20210608GalitOphir-270838574/SALMON_1.5.2/N831_L002-ds.9a3268502f2f4745827d7a83e6ccfac4')
q_gene<-read.csv("quant.sf",sep="\t")

setwd("D:/RNA_seq")
all_gene<-read.csv("dm6_refseq_common_names.tsv",sep="\t")

#

List_of_salmon <- list.files("D:/RNA_seq/20210608GalitOphir-270838574/SALMON_1.5.2", full.names=TRUE)

df_Salmon <- data.frame(matrix(ncol = length(List_of_salmon)+2, nrow = nrow(q_gene)))

#print(List_of_salmon)
group_name <<- tools::file_path_sans_ext(basename((List_of_salmon)))

df_Salmon[1]<-pull(q_gene,Name)

smaple_name<-tools::file_path_sans_ext(basename((List_of_salmon)))
colnames(df_Salmon)[1] <- "name"
colnames(df_Salmon)[2] <- "common_name"
colnames(df_Salmon)[3:34]<-smaple_name


for (sub_Salmon in List_of_salmon){
  setwd(sub_Salmon)
  smaple_current_name<-tools::file_path_sans_ext(basename((sub_Salmon)))
  q_gene<-read.csv("quant.sf",sep="\t")
  df_Salmon[smaple_current_name]= q_gene['TPM']
  
} 

all_gene_from_salamon<-filter(all_gene, X.name %in% df_Salmon$name)
df_Salmon[2]<-pull(all_gene_from_salamon,name2)




df_Salmon$sum <- apply(df_Salmon[3:34],1,sum)


most_exspressed<-df_Salmon$common_name[which.max(df_Salmon$sum)]

#present the 6 top expresses gene

df_Salmon %>%
  arrange(desc(df_Salmon$sum)) %>%
  head()