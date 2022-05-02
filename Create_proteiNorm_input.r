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

#Match to TMTs
    #Load TMT saved

    #Create index collumn
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

#Load protein info from PIA output
test1 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/test_1.csv', header = TRUE, sep = "\t")
test2 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/test_2.csv', header = TRUE, sep = "\t")
test3 <- read.csv(file = '/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/test_2.csv', header = TRUE, sep = "\t")
view(test1)
view(test2)
view(test3)



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