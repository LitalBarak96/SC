###DESCRIPTION: TRANSFER SALMON OUTPUT FROM TRANSCRIPT LEVEL TO GENE LEVEL ABUNDANCE + CALCULATE STATISTICS
## VER1: 18.7.2018
## VER2: 23.7.2018    removed sample filtering before statistics
## VER3: 18.8.2018    added filter by pattern and file prefix options, tailored for Salmon final script
## VER4: 21.8.2018    added argsparser, removed prefix + pattern

#!/usr/bin/env Rscript

# FUNCTIONS ---------------------------------------------------------------


#CBIND WHILE ALTERNATING
alternateCols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}


# ### SETUP ###
# #user_args are as written in the descriptions
# user_args <- cbind("/home/alu/dorony/gtf_files/hg38/transcriptIDtoGene_name", "/private7/projects/beta_cells/Oslo_data/Expression/Salmon/SE/samples.txt", "/private7/projects/beta_cells/Oslo_data/Expression/Salmon/SE/results", "/home/alu/fulther/BetaCells/SalmonTesting", "$", "DiViD2$")
# path_transcriptIDtoGene_name =user_args[1]
# path_barcode = user_args[2]
# path_tpm_results =user_args[3]
# out =user_args[4]

#CBIND WHILE ALTERNATING
alternateCols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}


#TPM FUNCTION
# path_transcriptIDtoGene_name: path of tanscript to gene ID file
# path_barcode: path for barcode file with sample names
# path_tpm_results: path for Salmon output to be processed
# out: path for directory for script output
tpm<- function(path_transcriptIDtoGene_name, path_barcode, path_tpm_results, out) {
  #notify of progress
  print(paste0("transcript ID to Gene path: ", path_transcriptIDtoGene_name), quote = F)
  print(paste0("barcode file path: ", path_barcode), quote = F)
  print(paste0("Salmon results path: ",path_tpm_results), quote = F)
  print(paste0("Output path: ",out), quote = F)
  
  #make sure all files are written to given path
  setwd(out)
  print(paste0("Workind directory (output dir): ",getwd()))
  
  library(data.table)
  #upload transcript to gene tabe
  tx2gene <-fread(path_transcriptIDtoGene_name,stringsAsFactors = F)
  tx2gene <-as.data.frame(tx2gene)
  
  #upload sample list = barcode file
  #prepare with:  (from directory of Salmon results)
  #               ls results > samples.txt
  #               sed -i -e '1i\Run'< samples.txt  samples.txt
  samples <- read.table(path_barcode, header = TRUE)
  files <- file.path(path_tpm_results,samples$Run)
  names(files) <- samples$Run
  
  
  #use tximport func to summarize to gene length
  library(tximport)
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  
  
  ###GET DATA
  #get TPM
  tpm<- txi$abundance
  colnames(tpm) <- gsub(".sf", "", colnames(tpm), fixed = T) #remove .sf from sample names
  #get counts
  counts <- txi$counts
  colnames(counts) <- gsub(".sf", "", colnames(counts), fixed = T) #remove .sf from sample names
  
  ### OUTPUT DATA
  out <- cbind(rownames(tpm), tpm) #make rownames a column
  colnames(out)[1] <- "GeneSymbol"
  print("writing abundance table...")
  write.csv(out,"Salmon_TPM.csv",row.names=F, quote = F)
  write.csv(out[grep("Obp69a", rownames(out), ignore.case = T), ],"Obp69a_Salmon_TPM.csv",row.names=F, quote = F)
  
  #add suffix to sample names to indeicate value
  colnames(tpm) <- lapply(colnames(tpm), paste0, "_TPM")
  colnames(counts) <- lapply(colnames(counts), paste0, "_counts")
  
  #create second table (with counts)
  out1 <- alternateCols(tpm, counts) #altenate columns
  out1 <- cbind(rownames(out1), out1) #make rownames a column
  colnames(out1)[1] <- "GeneSymbol"
  
  print("writing abundance and counts table...")
  write.csv(out1,"Salmon_TPM_wCounts.csv",row.names=F, quote = F)
  write.csv(out1[grep("Obp69a", rownames(out1), ignore.case = T), ],"Obp69a_Salmon_TPM_wCounts.csv",row.names=F, quote = F)
}


# SETUP -------------------------------------------------------------------

library("argparse")

# create parser object
parser <- ArgumentParser(description='Use tximport library to transfer transcript-level Salmon TPM output to gene-level TPM.')

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-i", action="store", dest = "indir", help="Path to Salmon output directory to be processed")
parser$add_argument("-o", action="store", dest = "outdir",  help="Output directory path")
parser$add_argument("-b", action="store", dest = "barcode", help="Path to barcode file with sample names")
parser$add_argument("-t", action="store", dest = "txToGeneID", default = NULL,  help="Path of tanscript to gene ID file")
parser$add_argument("-l", action="store", dest = "logPath", default = "log.txt",  help="Path of log file (appends)")





# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
user_args <- parser$parse_args()

#DEBUG

#_______________________________change folders from here,dm6_refseq_common_names for drozofila
user_args$txToGeneID <- "D:/RNA_seq/dm6_refseq_common_names.tsv"
user_args$barcode <- "D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/samples.txt"
user_args$indir <- "D:/RNA_seq/new_data/new/X201SC21111697-Z01-F001/SALMON_1.5.2/summery/"
user_args$outdir <- user_args$indir 
user_args$logPath <- paste0(user_args$indir,"/log.txt") 


#notify
print("", quote = F)
print("TXIMPORT SCRIPT PARAMETERS:", quote = F)
print(paste0("Log file: ",user_args$logPath), quote = F)
# direct stdoutput to a file 
sink(file = user_args$logPath, append=TRUE, split=FALSE)


#call function
library("tximport")
tpm(path_transcriptIDtoGene_name =user_args$txToGeneID, path_barcode = user_args$barcode,path_tpm_results =user_args$indir,out =user_args$outdir) 

%>% 