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
library("purrr")


#Working directory with all the wanted files
    #This wd has to contain all the files that have to be analyzed
mzTab_WD <- setwd("~/Desktop/Inputs/ALL_mzML")
getwd()
(file_names_wd <- list.files(mzTab_WD)) #6 mzML files
    #The wanted files
(file_names_short <- substring(fs::dir_ls("~/Desktop/Inputs/ALL_mzTab"), 86, 93)) #Character 86 untill 93 are unique

#Wanted mzMLs
    #Listing based on file_names_short
listed_mzMLs <- list()
for (i in seq_along(file_names_short)) {
  listed_mzMLs[[i]] <- dir(path = "~/Desktop/Inputs/ALL_mzML", pattern = file_names_short[[i]])
}
listed_mzMLs #from 16 mzML, the 6 from file_names_short are selected
    #Reading in listed mzMLs
mzML_files <- list() #empty list
for (i in seq_along(listed_mzMLs)) {
  mzML_files[[i]] <- readMSData(listed_mzMLs[[i]], msLevel = 2, verbose = FALSE, mode = "onDisk")
}
(Selected_mzML <- set_names(mzML_files, file_names_short))
 

#Wanted mzTabs
sel_mzTabs <- list()
for (i in seq_along(file_names_short)) {
  sel_mzTabs[[i]] <- dir(path = "~/Desktop/Inputs/ALL_mzTab", pattern = file_names_short[[i]])
}
sel_mzTabs

seleceted_mzTabs <- list[!sapply(list, identical, character(0))]
 
 
    #Reading in selected_mzTabs
mzTab_files <- list() #empty list
for (i in seq_along(seleceted_mzTabs)) {
  mzTab[[i]] <- readMzTab(seleceted_mzTabs[[i]])
}
Selected_mzTab <- set_names(mzML_files, file_names_short)
view(Selected_mzTab[[1]])