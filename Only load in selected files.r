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
file_names_short <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
file_names_short #6 file names, both mzML and mzTab loaded
    #Listing based on file_names_short
listed_mzMLs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzMLs[[i]] <- dir(pattern = file_names_short[[i]])
}
listed_mzMLs #from 16 mzML, the 6 from file_names_short are selected
    #Reading in listed mzMLs
mzML_files <- list() #empty list
for (i in seq_along(listed_mzMLs)) {
  mzML_files[[i]] <- readMSData(listed_mzMLs[[i]],
                                   msLevel = 2, verbose = FALSE, mode = "onDisk")
}
mzML_files



file_names_wd <- list.files(wd)
file_names_wd
mzML <- set_names(mzML_files, file_names_wd) #names each file by file_names_wd
mzML






    #mzTab

listed_mzTabs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzTabs[[i]] <- dir(path = "~/Desktop/mzTab/Imported mzTab", pattern = file_names_short[[i]])
}
listed_mzTabs #from 23 mzTab, the 6 from file_names_short are selected

?dir

setwd("~/Desktop/mzTab/Imported mzTab")
mzTab <- list() #empty list
for (i in seq_along(listed_mzTabs)) {
  mzTab[[i]] <- readMzTab(listed_mzTabs[[i]])
}




