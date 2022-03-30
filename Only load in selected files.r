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

#Working directory with all the wanted files
    #This wd has to contain all the files that have to be analyzed
mzML_WD <- setwd("~/Desktop/mzTab/mzML corresponding to mzTab")
getwd()
file_names_wd <- list.files(mzML_WD) #6 mzML files
    #The wanted files
file_names_short <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
file_names_short #6 file names, both mzML and mzTab loaded

#Wanted mzMLs
    #Listing based on file_names_short
listed_mzMLs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzMLs[[i]] <- dir(path = "~/Desktop/mzTab/mzML corresponding to mzTab", pattern = file_names_short[[i]])
}
listed_mzMLs #from 16 mzML, the 6 from file_names_short are selected
    #Reading in listed mzMLs
mzML_files <- list() #empty list
for (i in seq_along(listed_mzMLs)) {
  mzML_files[[i]] <- readMSData(listed_mzMLs[[i]], msLevel = 2, verbose = FALSE, mode = "onDisk")
}
Selected_mzML <- set_names(mzML_files, file_names_short)
Selected_mzML

#Wanted mzTabs
listed_mzTabs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzTabs[[i]] <- dir(path = "~/Desktop/mzTab/Imported mzTab", pattern = file_names_short[[i]])
}
listed_mzTabs #from 23 mzTab, the 6 from file_names_short are selected
    #Reading in listed mzTabs
mzTab_files <- list() #empty list
for (i in seq_along(listed_mzTabs)) {
  mzTab[[i]] <- readMzTab(listed_mzTabs[[i]])
}
Selected_mzTab <- set_names(mzML_files, file_names_short)
view(Selected_mzTab[[1]])