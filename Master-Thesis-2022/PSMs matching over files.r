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


#07/04/22 PSMs linked to TMT intensities for 23 files
PSM_TMT_07_04_22 <- readRDS("~/Desktop/Outputs/PSM_TMT_linked/07_04_22_PSM_TMT_Linked")
    #Make an extra index column that includes the file_names_short
new_col <- list()
for (i in seq_along(PSM_TMT_07_04_22)) {
    new_col[[i]] <- PSM_TMT_07_04_22[[i]] %>%
    mutate(index_filename = paste0(index, " ",file_names_short[[i]]), .before = index)
}
view(new_col[[4]])
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
    
    #Selecting the sequence_no_mod column + index_filename
      #Different result, because the modification are not always the same
seq_no_mod <- list()
for (i in seq_along(new_col)) {
    seq_no_mod[[i]] <- select(new_col[[i]], index_filename, sequence_no_mod)
}
(all_seq_no_mod <- bind_rows(seq_no_mod) %>% arrange(sequence_no_mod, index_filename))
view(all_seq_no_mod) #all rows arrenged by sequence_no_mod
Sum_all_seq_no_mod <- all_seq_no_mod %>%
  add_column(count = rep(1, nrow(all_seq_no_mod))) %>%
  group_by(sequence_no_mod) %>%
  summarise(Count = sum(count))
view(Sum_all_seq_no_mod) #sequence_no_mod summarized, with a counter
view(Sum_all_seq_no_mod %>% arrange(-Count))#most frequent sequence_no_mod at top

#Save the Outputs
  #Sequence_no_mod
saveRDS(all_seq_no_mod, file = "~/Desktop/Outputs/PSMs/07_04_22_Count_identical_seq_no_mod")
Count_identical_seq_no_mod_07_04_22 <- readRDS(file = "~/Desktop/Outputs/PSMs/07_04_22_Count_identical_seq_no_mod")
saveRDS(Sum_all_seq_no_mod, file = "~/Desktop/Outputs/PSMs/07_04_22_All_seq_no_mod")
All_seq_no_mod_07_04_22 <- readRDS(file = "~/Desktop/Outputs/PSMs/07_04_22_All_seq_no_mod")
  #Sequence
saveRDS(Sum_all_seq, file = "~/Desktop/Outputs/PSMs/07_04_22_Count_identical_seq_with_mod")
Count_identical_seq_with_mod_07_04_22 <- readRDS(file = "~/Desktop/Outputs/PSMs/07_04_22_Count_identical_seq_with_mod")
saveRDS(all_seq, file = "~/Desktop/Outputs/PSMs/07_04_22_All_seq_with_mod")
All_seq_with_mod_07_04_22 <- readRDS(file = "~/Desktop/Outputs/PSMs/07_04_22_All_seq_with_mod")


######################
#Extra: Starting from unchanged mzTab
######################
#Loading in PSMs
    #07/04/22: Loading in the 23 mzTabs collected to date
wd <- setwd("~/Desktop/Outputs/mzTab/Imported mzTab")
getwd() 
list.files(wd) #23 mzTabs in this wd
    #file paths to 23 files
(file_paths <- fs::dir_ls("~/Desktop/mzTab/Imported mzTab"))
    #Automate filename extraction
(file_names_wd <- list.files(wd)) #names of files in wd
(file_names_short <- substring(file_names_wd, 39, 46)) #Character 39 untill 46 are unique
  #Read an mzTab tab separated file 
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
  #Loop reading mzTabs
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
view(mzTab_files[[1]])
    #Selecting PSMs + create column names
        #extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
    #Loop for all files
    #Select sequence column + Removing brackets + Removing numbers
    #Add new column sequence_no_mod behind the sequence column
psm <- list() #empty list
for (i in seq_along(mzTab_files)) {
  psm[[i]] <- extractFeaturesPSM(mzTab_files[[i]]) %>%
  as_tibble() %>% row_to_names(row_number = 1) %>%
    mutate(sequence_no_mod = trimws(str_remove_all(sequence, "n"))) %>% #remove n
    mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "[0123456789]"))) %>% # remove numbers
    mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "\\[|\\]"))) %>% #remove []
    relocate(sequence_no_mod, .after = sequence)
}
PSM <- set_names(psm, file_names_short) #names each file by file_names_short
view(PSM[[1]])