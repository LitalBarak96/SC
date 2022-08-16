library("dplyr")  
library("ggpubr")
library(psych)
rectan=2
mth="pearson"
res_dir="D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/Deseq"
element="genes"
df<-avg_feamle_and_male;
set_name<-a
#df<-spec_mat
#filter out some not interesting genes
cor_mat_sig <- corr.test(df, method=mth)
if (rectan==0){
  corrplot::corrplot(cor_mat_sig$r, p.mat = coset_namer_mat_sig$p, method = "square",#type = "lower",
                     insig = "label_sig", sig.level = c(.001, .01, .05),
                     title=paste("\nsignificant",element,"correlation", ),
                     pch.cex = .9, pch.col = "white", order="hclust")
  
}else {
  corrplot::corrplot(cor_mat_sig$r, p.mat = cor_mat_sig$p, method = "square",#type = "lower",
                     insig = "label_sig", sig.level = c(.001, .01, .05),
                     title=paste("\nsignificant",element,"correlation", set_name),
                     pch.cex = .9, pch.col = "white", order="hclust", addrect = rectan)
}
dev.copy(jpeg,filename=paste0(res_dir,"/correlation_sig_",element,"",mth,"",set_name, ".jpg"), 
         width=ncol(df)*20+200, height=ncol(df)*20+250)
dev.off()




