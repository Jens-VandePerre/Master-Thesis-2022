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

#Loading in PSMs
    #01/04/22: Loading in the 23 mzTabs collected to date
wd <- setwd("~/Desktop/mzTab/Imported mzTab")
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

#Finding matching peptide sequences between files
#Work further with PSMs linked to TMTs, the index is needed
PSM_TMT <- readRDS("~/Desktop/mzTab/Stored files/PSMs linked to TMT intensities")
view(PSM_TMT[[1]])
    #Make an extra index column that includes the file_names_short
new_col <- list()
for (i in seq_along(PSM_TMT)) {
    new_col[[i]] <- PSM_TMT[[i]] %>%
    mutate(index_filename = paste0(index, " ",file_names_short[[i]]), .before = index)
}
view(new_col[[4]])
    #Selecting the sequence column + index_filename
seq <- list()
for (i in seq_along(new_col)) {
    seq[[i]] <- select(new_col[[i]], index_filename, sequence)
}
(all_seq <- bind_rows(seq))
view(all_seq)

?arrange

    #Selecting the sequence_no_mod column + index_filename
seq_no_mod <- list()
for (i in seq_along(new_col)) {
    seq_no_mod[[i]] <- select(new_col[[i]], index_filename, sequence_no_mod)
}
(all_seq_no_mod <- bind_rows(seq_no_mod) %>% arrange(sequence_no_mod, index_filename))
view(all_seq_no_mod)

