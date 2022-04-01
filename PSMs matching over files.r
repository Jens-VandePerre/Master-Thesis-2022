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

#Loading in PSMs
    #01/04/22: Loading in the 23 mzTabs collected to date
wd <- setwd("~/Desktop/mzTab/Imported mzTab")
getwd() 
list.files(wd) #23 mzTabs in this wd
    #file paths to 23 files
(file_paths <- fs::dir_ls("~/Desktop/mzTab/Imported mzTab"))
    #Automate filename extraction
(file_names_wd <- list.files(wd)) #names of files in wd
(file_names_short <- substring(file_names_wd, 39, 46)) #Character 39 untill 46 are unique
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
  #Loop reading mzTabs
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
view(mzTab_files[[1]])
    #Selecting PSMs + create column names
        #extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
        #Loop for all files
psm <- list() #empty list
for (i in seq_along(mzTab_files)) {
  psm[[i]] <- extractFeaturesPSM(mzTab_files[[i]]) %>%
  as_tibble() %>% row_to_names(row_number = 1)
}
PSM <- set_names(psm, file_names_short) #names each file by file_names_short
PSM
view(mzTab_files_PSM[[1]])



