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
    #Contains column without modification = column 2
PSM <- readRDS(file = "~/Desktop/mzTab/Stored files/PSM column no modifications v2")
view(PSM[[1]][2]) #column with peptide sequences, no modifications

#For all files, only keep peptides that present in all files
filter(PSM[[1]][2] != PSM[[2]][2])

