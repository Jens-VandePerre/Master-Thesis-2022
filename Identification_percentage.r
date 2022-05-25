library("limma")
library("qvalue")
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
library("rpx")
library("mzR")
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
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


#TMT spectra, count rows
TMT <- readRDS("/Users/jensvandeperre/Desktop/Outputs/TMTs/ALL_TMTs_16.05.22")
spec_count <- list()
for (i in 1:264) {
    spec_count[[i]] <- TMT[[i]] %>% 
    n_distinct()
}


#PSMs from ANN-SoLo
AS <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSMs/ALL_PSMs_4.5.22)
"AS_count <- list()
for (i in 1:264) {
    AS_count[[i]] <- ANN_SoLo[[i]] %>% 
    nrow()
}


#PSMs original study
(file_paths_Original <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/"))
ORIG <- list()
for (i in 1:264) {
  ORIG[[i]] <- read.csv(file_paths_Original[[i]], header = FALSE, sep = ",")
}
view(ORIG[[1]])
nrow(ORIG[[1]])

OS_count <- list()
for (i in 1:264) {
    OS_count[[i]] <- ORIG[[i]] %>% 
    nrow()
}


#ANN-Solo identification %
AS_IP <- list()
for(i in 1:264) {
    AS_IP[[i]] <- AS_count[[i]]/TMT[[i]]*100%
}

#ANN-Solo identification %
OS_IP <- list()
for(i in 1:264) {
    OS_IP[[i]] <- OS_count[[i]]/TMT[[i]]*100%
}
