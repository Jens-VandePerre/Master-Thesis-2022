library("rpx")
library("mzR")
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
library("limma")
library("isobar")
library("devtools")
library("MSnbase")
library("Biobase")
library("dplyr")
library("tidyverse")
library("fs")
library("proxyC")
library("sjlabelled")
library("expss")
library("labelled")
library("patchwork")
library("msdata")
library("remotes")
library("janitor")
library("stringr")

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd) #all mzTabs as of now 

#Create Function readMzTab
  #Read an mzTab tab separated file 
readMzTab <- function(filename) {
  # read maximum number of columns in file
  ncol <- max(stats::na.omit(utils::count.fields(file=filename, sep = "\t")))
  print(ncol)
  mztab.table = utils::read.table(file=filename, header=FALSE,
                                  row.names=NULL, dec = ".", fill = TRUE,
                                  col.names = paste0("V", seq_len(ncol)),
                                  sep="\t", na.strings="null", quote = "")
  mztab.table
}
#Loop Reading in Files
  #File paths to direct the loop
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"))
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_names_short <- substring(file_name_long, 39, 46)) #Characters 86 untill 93 are uniqueue
  #Loop reading mzTabs
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
view(mzTab_files[[1]])

#extractMetadata
  #Extracting the MTD: only columns first 3 have inforamtion
extractMetadata <- function(mztab.table) {
  mztab.table[startsWith(as.character(mztab.table$V1), "MTD"), 1:3]
}
  #Loop
MTD <- list() #empty list
for (i in seq_along(mzTab_files)) {
  MTD[[i]] <- extractMetadata(mzTab_files[[i]])
}
mzTab_files_Metadata <- set_names(MTD, file_names_short) #names each file by file_names_short
view(mzTab_files_Metadata[[1]])

#extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
  #Loop for all files
PSM <- list() #empty list
for (i in seq_along(mzTab_files)) {
  PSM[[i]] <- extractFeaturesPSM(mzTab_files[[i]])
}
mzTab_files_PSM <- set_names(PSM, file_names_short) #names each file by file_names_short
view(mzTab_files_PSM[[3]])

#Create column with unmodified peptide sequences
#Make Tibble with PSMs and use PSH as column names
#Select sequence column + Removing brackets + Removing numbers
#Add new column sequence_no_mod behind the sequence column
        #Loop for all files
psm <- list() #empty list
for (i in seq_along(mzTab_files)) {
  psm[[i]] <- extractFeaturesPSM(mzTab_files[[i]]) %>%
  as_tibble() %>% row_to_names(row_number = 1) %>%
  mutate(sequence = trimws(str_remove_all(sequence, "n"))) %>% #remove n
  mutate(sequence = trimws(str_remove_all(sequence, "[0123456789]"))) %>% # remove numbers
  mutate(sequence = trimws(str_remove_all(sequence, "\\[|\\]"))) %>% #remove []
}
PSM <- set_names(psm, file_names_short) #names each file by file_names_short
view(PSM[[1]])
