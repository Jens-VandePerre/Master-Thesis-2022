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


class(TMT_Matched_mzML_6)

view(TMT_Matched_mzML_6[[6]][,0])


test <- matrix(unlist(TMT_Matched_mzML_6[[1]]), nrow=length(TMT_Matched_mzML_6[[1]]), byrow=TRUE)
view(test)

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


tmt_tbl <- list()
for (i in seq_along(TMT_Matched_mzML_6)) { 
  tmt_tbl[[i]] <- as_tibble(TMT_Matched_mzML_6[[i]]) %>%
    select
}
TMT_indexed <- set_names(tmt_tbl, file_names_short)
TMT_indexed
view(TMT_indexed[[1]])

view(as_tibble(TMT_Matched_mzML_6[[1]]))

#Make PSM index column for matching
    #Select column spectra_ref
ind_mzTab1 <- list()
for (i in seq_along(mzTab_6)) {
    ind_mzTab1[[i]] <- select(mzTab_6[[i]], spectra_ref)
}
spectra_ind <- set_names(ind_mzTab1, file_names_short)
spectra_ind
    #Extract numbers from spectra_ref column
ind_mzTab2 <- list()
for (i in seq_along(spectra_ind)) {
ind_mzTab2[[i]] <- spectra_ind[[i]] %>% 
                mutate(index = trimws(str_remove_all(spectra_ref, "index="))) %>%
                mutate(index = trimws(str_remove_all(index, "ms_run"))) %>%
                mutate(index = trimws(str_remove_all(index, "\\[|\\]"))) %>%
                mutate(index = trimws(str_remove_all(index, "1:"))) %>%
                select(index) 
}
mzTab_index_col <- set_names(ind_mzTab2, file_names_short)
mzTab_index_col
    #Adding this index column to mzTab_6
ind_mzTab3 <- list()
for (i in seq_along(mzTab_6)) {
    ind_mzTab3[[i]] <- mzTab_6[[i]] %>%
        add_column(index = mzTab_index_col[[i]], .before = "sequence")
}
mzTab_6_ready_for_matching <- set_names(ind_mzTab3, file_names_short)
mzTab_6_ready_for_matching
view(mzTab_6_ready_for_matching[[1]])



#Merging the 2
mzTab_6_ready_for_matching

