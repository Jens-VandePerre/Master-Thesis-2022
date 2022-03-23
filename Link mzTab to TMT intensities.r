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
mzTab <- readRDS(file = "~/Desktop/mzTab/Stored files/Test 6 mzTabs")


#Load 6 matching mzMLs

#Extract TMT intensities from these 6 mzMLs

#Store this output
saveRDS(file = "~/Desktop/mzTab/Stored files/6 matched mzMLS")