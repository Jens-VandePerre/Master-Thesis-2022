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

#Load 6 PSMs
mzTab_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/6 PSM")
mzTab_6 #This is a Tibble

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
view(mzML)
#Extract TMT intensities from these 6 mzMLs
    #impute: method="MLE"
    #no normalisation
mzML_qnt <- list() #empty list
TMT <- list() #empty list
tbl <- list()
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
Matched_mzML_6 <- set_names(TMT, file_names_wd) #names each file by file_names_wd
Matched_mzML_6
    #Store this output
saveRDS(Matched_mzML_6, file = "~/Desktop/mzTab/Stored files/6 matched mzMLS")
TMT_Matched_mzML_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/6 matched mzMLS")


#Look for matching scan numbers
view(mzTab_6[[4]])
view(TMT_Matched_mzML_6[[4]][,0])

#Make TMT tibble + add index column for matching
    #Extracting spectral number column
ind <- list()
for (i in seq_along(TMT_Matched_mzML_6)) { 
  ind[[i]] <- TMT_Matched_mzML_6[,1][[i]] %>%
  as_tibble(ind[[i]])
}
TMT_col_add <- set_names(ind, file_names_short)
TMT_col_add
view(ind[[1]])

    #Loop making  Tibble + adding column
tmt_tbl <- list()
for (i in seq_along(TMT_Matched_mzML_6)) { 
  tmt_tbl[[i]] <- as_tibble(TMT_Matched_mzML_6[[i]]) %>%
    add_column(index = TMT_col_add[[i]], .before = "126")
}
TMT_indexed <- set_names(tmt_tbl, file_names_short)
TMT_indexed
view(TMT_indexed[[1]])



#Make PSM index column for matching
    #Select column spectra_ref
select(mzTab_6[[4]], spectra_ref)
spec_ref <- list()
for (i in seq_along(mzTab_6)) {
    spec_ref[[i]] <- select(mzTab_6[[i]], spectra_ref)
}
spectra_ref <- set_names(spec_ref, file_names_short)
spectra_ref

    #Extract numbers from spectra_ref column
indx <- sub(".*= ", "", spectra_ref[[1]]) 
view(indx)
