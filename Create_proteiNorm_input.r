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



#Read in mzTab without modifications
readMzTab <- function(filename) {
  # read maximum number of columns in file
  ncol <- max(stats::na.omit(utils::count.fields(file=filename, sep = "\t")))
  print(ncol)
  mztab.table = utils::read.table(file=filename, header=FALSE,
                                  row.names=NULL, dec = ".", fill = TRUE,
                                  col.names = paste0("V", seq_len(ncol)),
                                  sep="\t", na.strings="null", quote = "")
  mztab.table
}
    #Loop
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
view(mzTab_files[[1]])
 #Save mzTabs to different location
saveRDS(mzTab_files, file = "~/Desktop/Outputs/mzTabs_imported/mzTab_nomod_02.05.22")
mzTab_no_mod <- readRDS(file = "~/Desktop/Outputs/mzTabs_imported/mzTab_nomod_02.05.22")

#Select PSMs
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
  #Loop for all files
PSM <- list() #empty list
for (i in seq_along(mzTab_no_mod)) {
  PSM[[i]] <- extractFeaturesPSM(mzTab_no_mod[[i]]) %>% row_to_names(row_number = 1)
}
mzTab_no_mod_PSM <- set_names(PSM, file_names_short) #names each file by file_names_short
view(mzTab_no_mod_PSM[[3]])
    #Save PSM output
saveRDS(mzTab_no_mod_PSM, file = "~/Desktop/Outputs/PSMs/All_PSMs_nomod_02.05.22")
PSM_no_mod <- readRDS(file = "~/Desktop/Outputs/PSMs/All_PSMs_nomod_02.05.22")
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


PSM_B1S1_f01 <- select(PSM_no_mod[[1]], PSM_ID) %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
    mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
    mutate(index = trimws(str_remove_all(index, "scan="))) %>%
    select(index) %>%
    cbind(PSM_no_mod[[1]]) %>% as_tibble

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

TMT_B1S1_f01 <- tibble(index=rownames(TMT_part1_03.05.22[[1]][, 0])) %>%
    mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    select(index) %>%
    cbind(TMT_part1_03.05.22[[1]]) %>% as_tibble


TMT_ready_for_machting <- set_names(ind_TMT1, file_names_short)
TMT_ready_for_machting[[1]]
view(TMT_ready_for_machting[[1]])

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


test4_PRO$V1

scan(text=PROTEIN, sep=",", what="", quiet=TRUE)[3]



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
  #Leading rezor peptide
  test4_PRO %>% as_tibble %>% filter(V1, str_detect("PROTEIN"))
  ?select

test4_PRO %>% rename(new_name = starts_with(V1))
as_tibble
 filter(str_detect("PROTEIN"))






  #Gene names
  #Reporter intensity corrected
  #Reverse
  #Potential contaminant
  #id = just number
  #Protein groups IDs


#Create ProteinGroup.txt
  #ID = leading rezor peptide 
  #Reporter intensity corrected