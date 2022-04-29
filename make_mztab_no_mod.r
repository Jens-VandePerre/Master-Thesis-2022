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
library("faahKO")


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

#Joining MTD and PTM again
mzTab <- list()
for (i in 1:length(PSM)) {
    mzTab[[i]] <- rbind(MTD_long[[i]], psm[[i]])
}
view(mzTab[[1]])
mzTab_no_mod <- set_names(mzTab, file_name_long) #names each file by file_names_short
length(mzTab_no_mod)
length(file_name_long)


#Write these mzTabs, to be used as peptideindexer input
for(i in seq_along(mzTab_no_mod)) {                             
  write.table(mzTab_no_mod[i],                              
             file = paste0("Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq/",
                    file_name_long[i],
                    ".mztab"),
             row.names = FALSE, quote=FALSE, sep='\t')
}

