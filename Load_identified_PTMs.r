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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
    #Automate filename extraction
(file_name_long <- substring(list.files(wd), 1, 46))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"))
(file_names_short <- substring(file_name_long, 39, 46))
length(file_names_short)

#Load identified PTMs
(ptm_filepaths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PTM_identification_tol_10"))
PTM <- list()
for (i in 1:264) {
  PTM[[i]] <- read.csv(ptm_filepaths[[i]], sep = ",", header = TRUE)
}
view(PTM[[1]])
nrow(PTM[[1]])
length(PTM)
#Create index for matching to PSM_TMT
ptm_index <- list()
for (i in 1:264) {
  ptm_index[[i]] <- PTM[[i]] %>%
  as_tibble() %>%
  mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
  mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
  mutate(index = trimws(str_remove_all(index, "scan=")))
}
view(ptm_index[[1]])

#merge to PSM_TMT
PSM_TMT <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
nrow(PSM_TMT[[1]])
view(PSM_TMT[[1]])
PTM_PSM_TMT <- list()
for (i in 1:264) {
  PTM_PSM_TMT[[i]] <- merge(PSM_TMT[[i]], ptm_index[[i]], by = "index")
}
view(PTM_PSM_TMT[[1]])
nrow(PTM_PSM_TMT[[1]])

#Add a counter column
new_col <- list()
for (i in seq_along(PTM_PSM_TMT)) {
    new_col[[i]] <- PTM_PSM_TMT[[i]] %>%
    mutate(index_filename = paste0(index, " ",file_names_short[[i]]), .before = index)
}
view(new_col[[1]])
    #Selecting the sequence column + index_filename
seq <- list()
for (i in seq_along(new_col)) {
    seq[[i]] <- select(new_col[[i]], index_filename, sequence)
}



(all_seq <- bind_rows(seq) %>% arrange(sequence, index_filename))
view(all_seq) #all rows arrenged by sequence
Sum_all_seq <- all_seq %>% add_column(count = rep(1, nrow(all_seq))) %>% 
  group_by(sequence) %>%
  summarise(Count = sum(count))
view(Sum_all_seq) #sequence summarized, with a counter
view(Sum_all_seq %>% arrange(-Count)) #most frequent sequence at top


#Look for Cancer PTMs
#CRC
  #Phosphorylation
  #Citrunilation
  #Hydroxylation
#Cancer in general
  #Ubiquitination 
  #Methylation
  #Glycosilation
  #Acetylation




#Most prevalant PTMs




#How many types of PTMs
  #Count unique PTMs



#Most prevalent modified protein
AS_prot <- fread(file="/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/data_input.txt", sep="\t")
nrow(AS_prot)


(pia_filepaths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(pia_filepaths[[i]], sep = ",", header = TRUE)
}
view(PIA[[1]])
nrow(PIA[[1]])
length(PIA)




#Make PTM_PSM_TMT ready for merge to PIA output
#Select the needed columns
                    PTM_PSM_TMT_input <- list()
                    for (i in 1:264) {
                      PTM_PSM_TMT_input[[i]] <- PSM_TMT_ALL[[i]] %>% 
                      as_tibble() %>%
                      select(index, sequence, sequence_no_mod, "126":"131", charge) %>%
                      mutate(file_name = file_names_short[[i]]) 
                    }
                    view(PTM_PSM_TMT_input[[1]])
                    view(PTM_PSM_TMT_input[[2]])
                    view(PTM_PSM_TMT_input[[69]])
                    length(PTM_PSM_TMT_input)
