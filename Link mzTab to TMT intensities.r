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

#Load 6 mzTabs
mzTab_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/Test 6 mzTabs")

#Load 6 matching mzMLs
wd <- setwd("~/Desktop/mzTab/mzML corresponding to mzTab")
getwd() 
file_names_wd <- list.files(wd) #The first 6 mzML files
file_paths <- fs::dir_ls("~/Desktop/mzTab/mzML corresponding to mzTab")
file_paths 
mzML <- list() #empty list
for (i in seq_along(file_paths)) {
  mzML[[i]] <- readMSData(file_paths[[i]],
                         msLevel = 2, verbose = FALSE, mode = "onDisk")
}
mzML <- set_names(mzML, file_names_wd) #names each file by file_names_wd
mzML
#Extract TMT intensities from these 6 mzMLs
mzML_qnt <- list() #empty list
TMT <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt[[i]] <- 
    quantify(mzML[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    impute(method="MLE") %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE))
  for (j in seq_along(mzML_qnt)) {
    TMT[[j]] <- exprs(mzML_qnt[[j]]) #output all spectra, unclear in terminal
  }
}
TMT_Matched_mzML_6 <- set_names(TMT, file_names_wd) #names each file by file_names_wd
TMT_Matched_mzML_6

#Store this output
saveRDS(TMT_Matched_mzML_6, file = "~/Desktop/mzTab/Stored files/6 matched mzMLS")
Matched_mzML_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/6 matched mzMLS")