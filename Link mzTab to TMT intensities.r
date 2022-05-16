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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"))
(file_names_short <- substring(file_name_long, 39, 46)) 
length(file_names_short)

#Load PSMs 
PSM_all <- readRDS(file = "~/Desktop/Outputs/PSMs/ALL_PSMs_4.5.22") 
view(PSM_all[[1]]) 
length(PSM_all)
  #Create column no_mod
psm <- list() #empty list
for (i in seq_along(PSM_all)) {
  psm[[i]] <- PSM_all[[i]] %>%
  as_tibble() %>% 
  mutate(sequence_no_mod = trimws(str_remove_all(sequence, "n"))) %>% #remove n
  mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "[0123456789]"))) %>% # remove numbers
  mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "\\[|\\]"))) %>% #remove []
  relocate(sequence_no_mod, .after = sequence)
}
PSM <- set_names(psm, file_names_short) #names each file by file_names_short
view(PSM[[1]])

#Load matching TMTs
TMT_all <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/ALL_TMTs_16.05.22") 

#Look for matching scan numbers
view(PSM_all[[1]]["PSM_ID"]) #Column PSM_ID
view(TMT_all[[1]]$index) #These are row names

#Make PSM index column for matching
  #Extract numbers from PSM_ID column + Adding this index column to PSM_6
ind_PSM <- list()
for (i in seq_along(PSM)) {
  ind_PSM[[i]] <- select(PSM[[i]], PSM_ID) %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
    mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
    mutate(index = trimws(str_remove_all(index, "scan="))) %>%
    select(index) %>%
    cbind(PSM[[i]]) %>% as_tibble
}
PSM_ready_for_matching <- set_names(ind_PSM, file_names_short)
PSM_ready_for_matching
view(PSM_ready_for_matching[[1]])

#Merging by index
merging <- list()
for (i in seq_along(TMT_all)) {
  merging[[i]] <- merge(PSM_ready_for_matching[[i]], TMT_all[[i]], by="index") %>% 
  as_tibble
}
Merged_PSM_TMT <- set_names(merging, file_names_short)
view(Merged_PSM_TMT[[1]]) 
view(merging[[1]] %>% select(index, sequence_no_mod, 27:36) %>% arrange(sequence_no_mod)) 
  #Save outputs
saveRDS(Merged_PSM_TMT, file = "~/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
PSM_TMT_all <- readRDS("~/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")

#Checking if length stay the same after matching PSMs and TMTs
  #Length PSM
l_PSM <- list()
for (i in seq_along(PSM_TMT_all)) {
  l_PSM[[i]] <- nrow(PSM_TMT_all[[i]])
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
for (i in seq_along(PSM_TMT_all)) {
  l_diff[[i]] <- (PSM_length[[i]]-Merged_length[[i]])
}
(Length_difference <- set_names(l_diff, file_names_short))#All 0 Merging SUCCESS
