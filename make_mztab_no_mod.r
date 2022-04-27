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
extractMetadata_long <- function(mztab.table) {
  mtd <- mztab.table[startsWith(as.character(mztab.table$V1), "MTD"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
    rbind(mtd,psh)

}
  #Loop
MTD_long <- list() #empty list
for (i in seq_along(mzTab_files)) {
  MTD_long[[i]] <- extractMetadata_long(mzTab_files[[i]])
}
mzTab_files_Metadata_long <- set_names(MTD_long, file_names_short) #names each file by file_names_short
view(mzTab_files_Metadata_long[[1]])

#extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  rbind(psm)
}
  #Loop for all files
PSM <- list() #empty list
for (i in seq_along(mzTab_files)) {
  PSM[[i]] <- extractFeaturesPSM(mzTab_files[[i]])
}
mzTab_files_PSM <- set_names(PSM, file_names_short) #names each file by file_names_short
view(mzTab_files_PSM[[3]])

#Remove PTMs from sequence collumn
        #Loop for all files
psm <- list() #empty list
for (i in seq_along(mzTab_files)) {
  psm[[i]] <- mzTab_files_PSM[[i]] %>%
  mutate(V2 = trimws(str_remove_all(V2, "n"))) %>% #remove n
  mutate(V2 = trimws(str_remove_all(V2, "[0123456789]"))) %>% # remove numbers
  mutate(V2 = trimws(str_remove_all(V2, "\\[|\\]"))) #remove []
}
PSM <- set_names(psm, file_names_short) #names each file by file_names_short
view(PSM[[1]])

#Joining MTD and PTM again
mzTab <- list()
for (i in 1:length(PSM)) {
    mzTab[[i]] <- rbind(MTD_long[[i]], psm[[i]])
}
view(mzTab[[1]])
MZTABBRO <- set_names(mzTab, file_name_long) #names each file by file_names_short


#Write these mzTabs, to be used as peptideindexer input
write.csv(MZTABBRO[[1]], file="/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq/test.mztab")
?writeMzTabData
as_data_frame()

writeMzTabData(MZTABBRO[[1]], what = "PEP", file = "/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq/test.mztab"
               )
