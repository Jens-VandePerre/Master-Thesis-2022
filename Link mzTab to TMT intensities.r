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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/mzTab_19_04_22")
getwd() 
list.files(wd)

#Load PSMs 
PSM_07_04_22 <- readRDS(file = "~/Desktop/Outputs/PSMs/20_04_22_PSMs") #PSMs stored 07/04/22
view(PSM_07_04_22[[1]]) 
#Load matching mzMLs
TMT_Intensities_06_04_22 <- readRDS(file = "~/Desktop/Outputs/TMTs/20.04.22_TMT") #TMT intensities stored 06/04/22
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/mzTab_19_04_22"))
(file_names_short <- substring(file_paths, 86, 93)) #Characters 86 untill 93 are uniqueue

#Look for matching scan numbers
view(PSM_07_04_22[[1]]["PSM_ID"]) #Column PSM_ID
view(TMT_Intensities_06_04_22[[1]][, 0]) #These are row names

#Make TMT tibble + add index column for matching
  #selecting index column
ind_TMT1 <- list()
for (i in seq_along(TMT_Intensities_06_04_22)) {
  ind_TMT1[[i]] <- tibble(index=rownames(TMT_Intensities_06_04_22[[i]][, 0])) %>%
    mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    select(index) %>%
    cbind(TMT_Intensities_06_04_22[[i]]) %>% as_tibble
}
TMT_ready_for_machting <- set_names(ind_TMT1, file_names_short)
TMT_ready_for_machting[[1]]
view(TMT_ready_for_machting[[1]])

#Make PSM index column for matching
  #Extract numbers from PSM_ID column + Adding this index column to PSM_6
ind_mzTab5 <- list()
for (i in seq_along(PSM_07_04_22)) {
  ind_mzTab5[[i]] <- select(PSM_07_04_22[[i]], PSM_ID) %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
    mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
    mutate(index = trimws(str_remove_all(index, "scan="))) %>%
    select(index) %>%
    cbind(PSM_07_04_22[[i]]) %>% as_tibble
}
PSM_ready_for_matching <- set_names(ind_mzTab5, file_names_short)
PSM_ready_for_matching
view(PSM_ready_for_matching[[1]]) #PSM column could be removed

#Merging by index
merging <- list()
for (i in seq_along(TMT_ready_for_machting)) {
  merging[[i]] <- merge(PSM_ready_for_matching[[i]], TMT_ready_for_machting[[i]], by="index") %>% 
  as_tibble
}
Merged_PSM_TMT <- set_names(merging, file_names_short)
Merged_PSM_TMT
view(Merged_PSM_TMT[[1]]) 
view(merging[[1]] %>% select(index, sequence_no_mod, 25:34) %>% arrange(sequence_no_mod)) 
  #Save outputs
saveRDS(Merged_PSM_TMT, file = "~/Desktop/Outputs/PSM_TMT_linked/07_04_22_PSM_TMT_Linked")
PSM_TMT_07_04_22 <- readRDS("~/Desktop/Outputs/PSM_TMT_linked/07_04_22_PSM_TMT_Linked")

#Checking if length stay the same after matching PSMs and TMTs
  #Length PSM
l_PSM <- list()
for (i in seq_along(PSM_TMT_07_04_22)) {
  l_PSM[[i]] <- nrow(PSM_TMT_07_04_22[[i]])
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
for (i in seq_along(PSM_TMT_07_04_22)) {
  l_diff[[i]] <- (PSM_length[[i]]-Merged_length[[i]])
}
(Length_difference <- set_names(l_diff, file_names_short))#All 0 Merging SUCCESS

#Selecting the collumn for relative quantification
selected <- list()
for (i in seq_along(Merged_PSM_TMT)) {
  selected[[i]] <- Merged_PSM_TMT[[i]] %>% 
  select("sequence_no_mod", "126":"131") %>%
  rename(Peptide_sequence = sequence_no_mod, 
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
