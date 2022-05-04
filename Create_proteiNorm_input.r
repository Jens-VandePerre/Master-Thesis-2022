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


wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq")
getwd() 
list.files(wd)
#Load PSMs 
#Load matching mzMLs
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab_pure_seq"))
(file_names_short <- substring(file_paths, 91, 98)) 

ALL_PSMs_4.5.22 <- readRDS(file = "~/Desktop/Outputs/PSMs/ALL_PSMs_4.5.22")
view(ALL_PSMs_4.5.22[1])

nrow(ALL_PSMs_4.5.22[1] %>% as_tibble)

  #Loop for all files
    #Make an index collumn for mathcing to TMTs
ind_mzTab <- list()
for (i in seq_along(PSM_no_mod)) {
  ind_mzTab[[i]] <- select(PSM_no_mod[[i]], PSM_ID) %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
    mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
    mutate(index = trimws(str_remove_all(index, "scan="))) %>%
    select(index) %>%
    cbind(PSM_no_mod[[i]]) %>% as_tibble
}
PSM_ready_for_matching <- set_names(ind_mzTab, file_names_short)
PSM_ready_for_matching
view(PSM_ready_for_matching[[1]]) #PSM column could be removed




names(ALL_PSMs_4.5.22[1])
PSM_B1S1_f01 <- ALL_PSMs_4.5.22[[1]] %>% 
    as_tibble %>%
    select(PSM_ID) %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
    mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
    mutate(index = trimws(str_remove_all(index, "scan="))) %>%
    select(index) %>%
    cbind(ALL_PSMs_4.5.22[[1]]) %>% as_tibble

#Match to TMTs
    #Load TMT saved
#TMT_06.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/06.04.22_TMT")
#TMT_20.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/20.04.22_TMT")
#TMT_22.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/22.04.22_TMT")
#TMT_28.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/28.04.22_TMT")

TMT_part1_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part1")

#TMT_part2_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part2")
#TMT_part3_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part3")


    #Create index collumn
ind_TMT1 <- list()
for (i in seq_along(TMT_part1_03.05.22)) {
  ind_TMT1[[i]] <- tibble(index=rownames(TMT_part1_03.05.22[[i]][, 0])) %>%
    mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    select(index) %>%
    cbind(TMT_part1_03.05.22[[i]]) %>% as_tibble
}
TMT_ready_for_machting <- set_names(ind_TMT1, file_names_short)
TMT_ready_for_machting[[1]]
view(TMT_ready_for_machting[[1]])


TMT_B1S1_f01 <- tibble(index=rownames(TMT_part1_03.05.22[[1]][, 0])) %>%
    mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    select(index) %>%
    cbind(TMT_part1_03.05.22[[1]]) %>% as_tibble



#Load protein info from PIA output
test1 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/test_1.csv', header = FALSE, sep = ",", stringsAsFactors = TRUE)
test2 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/test_2.csv', header = FALSE, sep = ",", stringsAsFactors = TRUE)
test3 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/test_3.csv', header = FALSE, sep = ",", stringsAsFactors = TRUE)
test4_PSM <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/PSM_exp_test_4.csv', header = FALSE, sep = ",", stringsAsFactors = TRUE)
test4_PRO <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/PROT_exp_test_4.csv', header = FALSE, sep = ",", stringsAsFactors = TRUE)
test5 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/PROT_exp_test_5.csv', header = FALSE, sep = ",", stringsAsFactors = TRUE)
 

view(test1)
view(test2)
view(test3)
view(test4_PSM)
view(test4_PRO)
view(test5)


protein <- test4_PRO %>% as_tibble %>% filter(str_detect(V1, "PROTEIN")) %>% row_to_names(row_number = 1)
peptide <- test4_PRO %>% as_tibble %>% filter(str_detect(V1, "PEPTIDE")) %>% row_to_names(row_number = 1)
psmset <- test4_PRO %>% as_tibble %>% filter(str_detect(V1, "PSMSET")) %>% row_to_names(row_number = 1)
psm <- test4_PRO %>% as_tibble  %>%  filter(str_detect(V1, "PSM")) %>% row_to_names(row_number = 1)


view(protein) #columns OK
view(peptide) #modifications column emprty
view(psmset) #decoy column not aligned
view(psm)




#Make TMT tibble + add index column for matching
  #selecting index column
ind_TMT1 <- list()
for (i in seq_along(TMT_Intensities_22_04_22)) {
  ind_TMT1[[i]] <- tibble(index=rownames(TMT_Intensities_22_04_22[[i]][, 0])) %>%
    mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    select(index) %>%
    cbind(TMT_Intensities_22_04_22[[i]]) %>% as_tibble
}
TMT_ready_for_machting <- set_names(ind_TMT1, file_names_short)
TMT_ready_for_machting[[1]]
view(TMT_ready_for_machting[[1]])



#Create Peptide.txt
PepSeq_ProAcc <- test5 %>% 
          as_tibble %>% 
          filter(str_detect(V1, "PSMSET")) %>% 
          row_to_names(row_number = 1) %>%
          select(Sequence, Accessions)
view(PepSeq_ProAcc)
nrow(PepSeq_ProAcc)
n_distinct(PepSeq_ProAcc)

ClusID_Des <- test5 %>% 
          as_tibble %>% 
          filter(str_detect(V1, "PROTEIN")) %>% 
          row_to_names(row_number = 1) %>%
          select(Proteins, ClusterID, Description) %>%
          rename(Accessions = Proteins)
view(ClusID_Des)
nrow(ClusID_Des)
n_distinct(ClusID_Des)

mergebro <- merge(PepSeq_ProAcc, ClusID_Des, by="Accessions")
nrow(mergebro)
view(mergebro)
n_distinct(mergebro)

mztab_TMT <- merge(TMT_B1S1_f01, PSM_B1S1_f01, by="index") %>%
        select(index, "126":"131", sequence, sequence_no_mod)
view(mztab_TMT)

mztab_TMT %>%
  rename(`Repoter intensity corrected 126` = `126`,
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


  #Leading rezor peptide
  #Gene names = Description
  #Reporter intensity corrected
  #Reverse = empty
  #Potential contaminant = empty
  #id = just number
  #Protein groups IDs = ClusterID????


#Create ProteinGroup.txt
  #ID = leading rezor peptide 
  #Reporter intensity corrected