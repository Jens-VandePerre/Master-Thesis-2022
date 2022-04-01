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
PSM_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/PSM column no modifications v2")
PSM_6 #This is a Tibble
view(PSM_6[[1]])

#Load 6 matching mzMLs
wd <- setwd("~/Desktop/mzTab/mzML corresponding to mzTab")
getwd() 
file_names_wd <- list.files(wd) #The first 6 mzML files
file_names_short <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
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
view(PSM_6[[1]]["PSM_ID"]) #Column PSM_ID
view(TMT_Matched_mzML_6[[1]][, 0]) #These are row names

#Make TMT tibble + add index column for matching
  #selecting index column
ind_TMT1 <- list()
for (i in seq_along(TMT_Matched_mzML_6)) {
  ind_TMT1[[i]] <- tibble(index=rownames(TMT_Matched_mzML_6[[i]][, 0])) %>%
  mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
  mutate(index = trimws(str_remove_all(index, "^0"))) %>%
  mutate(index = trimws(str_remove_all(index, "^0"))) %>%
  mutate(index = trimws(str_remove_all(index, "^0"))) %>%
  mutate(index = trimws(str_remove_all(index, "^0"))) %>%
  select(index) %>%
  cbind(TMT_Matched_mzML_6[[i]]) %>% as_tibble
}
TMT_ready_for_machting <- set_names(ind_TMT1, file_names_short)
TMT_ready_for_machting[[1]]
view(TMT_ready_for_machting[[1]])

#Make PSM index column for matching
  #Extract numbers from PSM_ID column + Adding this index column to PSM_6
ind_mzTab5 <- list()
for (i in seq_along(PSM_6)) {
ind_mzTab5[[i]] <- select(PSM_6[[i]], PSM_ID) %>% 
                mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
                mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
                mutate(index = trimws(str_remove_all(index, "scan="))) %>%
                select(index) %>%
                cbind(PSM_6[[i]]) %>% as_tibble
}
PSM_ready_for_matching <- set_names(ind_mzTab5, file_names_short)
PSM_ready_for_matching
view(PSM_ready_for_matching[[1]])

#Merging the 2
view(PSM_ready_for_matching[[1]])
view(TMT_ready_for_machting[[1]])
merge1 <- merge(PSM_ready_for_matching[[1]], TMT_ready_for_machting[[1]], by="index")
view(merge1)

merging <- list()
for (i in seq_along(TMT_ready_for_machting)) {
  merging[[i]] <- merge(PSM_ready_for_matching[[i]], TMT_ready_for_machting[[i]], by="index") %>% 
  as_tibble %>% arrange(desc(index[[i]]))
}
Merged_PSM_TMT <- set_names(merging, file_names_short)
Merged_PSM_TMT
view(Merged_PSM_TMT[[1]])
  #Save outputs
saveRDS(Merged_PSM_TMT, file = "~/Desktop/mzTab/Stored files/PSMs linked to TMT intensities")
PSM_TMT <- readRDS("~/Desktop/mzTab/Stored files/PSMs linked to TMT intensities")


#Checking if length stay the same after matching PSMs and TMTs
  #Length PSM
l_PSM <- list()
for (i in seq_along(PSM_6)) {
  l_PSM[[i]] <- nrow(PSM_6[[i]])
}
(PSM_length <- set_names(l_PSM, file_names_short))
  #Length matched
l_matched <- list()
for (i in seq_along(Merged_PSM_TMT)) {
  l_matched[[i]] <- nrow(Merged_PSM_TMT[[i]])
}
(Merged_length <- set_names(l_matched, file_names_short)) 
  #Is there a difference?
l_diff <- list()
for (i in seq_along(PSM_6)) {
  l_diff[[i]] <- (PSM_length[[i]]-Merged_length[[i]])
}
(Length_difference <- set_names(l_diff, file_names_short)) #All lengts are the same. Merging SUCCESS!!!

#Selecting the collumn for relative quantification
selected <- list()
for (i in seq_along(Merged_PSM_TMT)) {
  selected[[i]] <- Merged_PSM_TMT[[i]] %>% 
  select("sequence_no_mod", "126":"131") %>%
  rename(Peptide_sequence=sequence_no_mod, 
        `Repoter intensity corrected 126` = `126`,
        `Repoter intensity corrected 127N` = `127N`,
        `Repoter intensity corrected 127C` = `127C`,
        `Repoter intensity corrected 128N` = `128N`,
        `Repoter intensity corrected 128C` = `128C`,
        `Repoter intensity corrected 129N` = `129N`,
        `Repoter intensity corrected 129C` = `129C`,
        `Repoter intensity corrected 130N` = `130N`,
        `Repoter intensity corrected 130C` = `130C`,
        `Repoter intensity corrected 113` = `131`
        )
}
view(selected[[1]])
