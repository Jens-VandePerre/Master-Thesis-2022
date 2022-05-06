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




names(ALL_PSMs_4.5.22)
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
names(TMT_part1_03.05.22)
view(TMT_part1_03.05.22[[1]])

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
view(TMT_B1S1_f01)


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
mztab_TMT <- merge(TMT_B1S1_f01, PSM_B1S1_f01, by="index") %>%
        select(index, "126":"131", sequence, sequence_no_mod)
view(mztab_TMT)

mztab_TMT <- mztab_TMT %>%
  rename(`Reporter intensity corrected 1 TMT126` = `126`,
        `Reporter intensity corrected 1 TMT127N` = `127N`,
        `Reporter intensity corrected 2 TMT127C` = `127C`,
        `Reporter intensity corrected 2 TMT128N` = `128N`,
        `Reporter intensity corrected 2 TMT128C` = `128C`,
        `Reporter intensity corrected 2 TMT129N` = `129N`,
        `Reporter intensity corrected 2 TMT129C` = `129C`,
        `Reporter intensity corrected 2 TMT130N` = `130N`,
        `Reporter intensity corrected 2 TMT130C` = `130C`,
        `Reporter intensity corrected 3 TMT131` = `131`
        ) %>% as_tibble 
view(mztab_TMT)



#Create Peptide.txt
PepSeq_ProAcc <- test5 %>% 
          as_tibble %>% 
          filter(str_detect(V1, "PSMSET")) %>% 
          row_to_names(row_number = 1) %>%
          select(Sequence, Accessions) %>%
          rename(sequence_no_mod = Sequence)
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

ProteiNorm <- merge(PepSeq_ProAcc, ClusID_Des, by= "Accessions") %>%
            as_tibble %>%
            rename("Leading razor peptide" = Accessions, "Protein group IDs" = ClusterID) %>%
            merge(mztab_TMT, by = "sequence_no_mod") %>% distinct()


nrow(ProteiNorm)
view(ProteiNorm)
n_distinct(ProteiNorm)
nrow(mztab_TMT)
n_distinct(mztab_TMT)

rep(" ", nrow(ProteiNorm))

          add_column("Reverse" = rep(" ", nrow(ProteiNorm)), .after = "Reporter intensity corrected 3 TMT131") 
          add_column("Potential contaminant" = rep(" ", nrow(ProteiNorm)), .after = "Reverse") 
          add_column("id" = 1:nrow(ProteiNorm), .after="Potential contaminant") 


peptide_txt <- ProteiNorm %>%
          as_tibble %>%
          rename("Gene names" = Description) %>%
          select("Leading razor peptide", "Gene names", "Reporter intensity corrected 1 TMT126":"Reporter intensity corrected 3 TMT131", "Protein group IDs") %>%
          add_column("Reverse" = rep(" ", nrow(ProteiNorm)), .after = "Reporter intensity corrected 3 TMT131") %>%
          add_column("Potential contaminant" = rep(" ", nrow(ProteiNorm)), .after = "Reverse") %>%
          add_column("id" = 1:nrow(ProteiNorm), .after="Potential contaminant") 
write.table(protein_txt, file="/Users/jensvandeperre/Desktop/Inputs/ProteiNorm/peptide/peptide.txt", append = FALSE, sep = "\t", dec = ".",
             col.names = TRUE)

view(peptide_txt)
#add id, empty Reverse and empty Potential contaminant


proteinGroup_txt <- ProteiNorm %>%
           as_tibble %>%
           rename(id = "Leading razor peptide") %>%
           select(id, `Reporter intensity corrected 1 TMT126`:`Reporter intensity corrected 3 TMT131`)
view(proteinGroup_txt)
write.table(proteinGroup_txt, file="/Users/jensvandeperre/Desktop/Inputs/ProteiNorm/protein/proteinGroup.txt", append = FALSE, sep = "\t", dec = ".",
             col.names = TRUE)


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