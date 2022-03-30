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

#mzML
mzML_WD <- setwd("~/Desktop/mzTab/mzML corresponding to mzTab")
getwd()
file_names_wd <- list.files(mzML_WD) #6 mzML files
    #The wanted files
file_names_short_typed <- c("B1S1_f10","B2S4_f10","B3S2_f09","B3S4_f04","B3S4_f06","B5S1_f08","B5S2_f04","B5S2_f07","B5S5_f04","B5S5_f08")
file_names_short <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
file_names_short #6 file names, both mzML and mzTab loaded




#Loading in matching mzMLs to mzTabs in wd
  #Names of wanted files

  #NOT all files here are matching
list.files(path = "~/Desktop/ALL mzML") #16 mzML files
list.files(path = "~/Desktop/mzTab/Imported mzTab") #23 mzTab files

  #Listing basted on string
list.files(path = "~/Desktop/ALL mzML", pattern = "B1S2_f10") #both present in mzTab and mzML
list.files(path = "~/Desktop/mzTab/Imported mzTab", pattern = "B1S2_f10") #both present in mzTab and mzML

  #Listing based on file_names_short
list.files(path = "~/Desktop/ALL mzML", pattern = file_names_short) #Gives only the first match
    #mzML
listed_mzMLs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzMLs[[i]] <- list.files(path = "~/Desktop/ALL mzML", pattern = file_names_short[[i]])
}
listed_mzMLs #from 16 mzML, the 6 from file_names_short are selected
    #mzTab
listed_mzTabs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzTabs[[i]] <- list.files(path = "~/Desktop/mzTab/Imported mzTab", pattern = file_names_short[[i]])
}
listed_mzTabs #from 23 mzTab, the 6 from file_names_short are selected

listed_mzTabs2 <- list()
for (i in seq_along(file_names_short)) {
  listed_mzTabs2[[i]] <- dir(pattern = file_names_short[[i]])
}
listed_mzTabs2 #from 23 mzTab, the 6 from file_names_short are selected



setwd("~/Desktop/mzTab/Imported mzTab")
mzTab <- list() #empty list
for (i in seq_along(listed_mzTabs2)) {
  mzTab[[i]] <- readMzTab(listed_mzTabs2[[i]])
}


mzML_files <- list() #empty list
for (i in seq_along(file_paths)) {
  mzML_files[[i]] <- readMSData(file_paths[[i]],
                                   msLevel = 2, verbose = FALSE, mode = "onDisk")
}
file_names_wd <- list.files(wd)
file_names_wd
mzML <- set_names(mzML_files, file_names_wd) #names each file by file_names_wd
mzML


