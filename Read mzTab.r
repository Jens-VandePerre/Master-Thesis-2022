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

wd <- setwd("~/Desktop/mzTab/Imported mzTab")
getwd() 
list.files(wd)

#Create Function readMzTab
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
(Test <- readMzTab("02CPTAC_COprospective_W_PNNL_20170123_B1S2_f10.mztab")) #The first file
#Loop Reading in Files
  #File paths to direct the loop
file_paths <- fs::dir_ls("~/Desktop/mzTab/Imported mzTab")
file_paths 
    #Automate filename extraction
file_names_wd <- list.files(wd) #names of files in wd
file_names_short <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
  #Loop reading mzTabs
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
mzTab_files[[1]]

 #Save mzMLs to different location: TMT outputs
saveRDS(mzTab_files, file = "~/Desktop/mzTab/Stored files/Test 6 mzTabs")
Test_6_mzTabs <- readRDS(file = "~/Desktop/mzTab/Stored files/Test 6 mzTabs")

##################
#Extract functions
##################

#extractMetadata_long
extractMetadata_long <- function(mztab.table) {
  mztab.table[startsWith(as.character(mztab.table$V1), "MTD"),]
}
  #Loop
MTD_long <- list() #empty list
for (i in seq_along(mzTab_files)) {
  MTD_long[[i]] <- extractMetadata_long(mzTab_files[[i]])
}
mzTab_files_Metadata_long <- set_names(MTD_long, file_names_short) #names each file by file_names_short
mzTab_files_Metadata_long

#extractMetadata
  #Extracting the MTD: only columns first 3 have inforamtion
extractMetadata <- function(mztab.table) {
  mztab.table[startsWith(as.character(mztab.table$V1), "MTD"), 1:3]
}
  #Loop
MTD <- list() #empty list
for (i in seq_along(mzTab_files)) {
  MTD[[i]] <- extractMetadata(mzTab_files[[i]])
}
mzTab_files_Metadata <- set_names(MTD, file_names_short) #names each file by file_names_short
mzTab_files_Metadata

#extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
  #Loop for all files
PSM <- list() #empty list
for (i in seq_along(mzTab_files)) {
  PSM[[i]] <- extractFeaturesPSM(mzTab_files[[i]])
}
mzTab_files_PSM <- set_names(PSM, file_names_short) #names each file by file_names_short
mzTab_files_PSM
view(mzTab_files_PSM[[1]])

nrow(mzTab_files_PSM[[4]])
n_distinct(mzTab_files_PSM[[4]])
nrow(unique(mzTab_files_PSM[[4]]))

#Determine identification percentage
  #Row count for mzTab_PSM = amount of identified peptides
nrow_PSM <- list() #empty list
for (i in seq_along(mzTab_files_PSM)) {
  nrow_PSM[[i]] <- as_tibble(mzTab_files_PSM[[i]]) %>%
    row_to_names(row_number = 1) %>%
    n_distinct() #count unique rows
}
Identified_peptides <- set_names(nrow_PSM, file_names_short) #names each file by file_names_short
Identified_peptides
  #Row count original
TMT_Matched_mzML_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/6 matched mzMLS")
n_distinct(TMT_Matched_mzML_6[[1]]) #count unique rows

    #Loop counting spectra
nrow_TMT <- list() #empty list
for (i in seq_along(TMT_Matched_mzML_6)) { #Only for the first 6 spectra, 
  nrow_TMT[[i]] <- n_distinct(TMT_Matched_mzML_6[[i]])
}
Spectral_count <- set_names(nrow_TMT, file_names_short) #names each file by file_names_short
Spectral_count

  #Loop percentage calculation
perc <- list() #empty list
for (i in seq_along(Identified_peptides)) { #Only for the first 6 spectra, 
  perc[[i]] <- (Identified_peptides[[i]]/Spectral_count[[i]])*100
}
Identification_percentage <- set_names(perc, file_names_short) #names each file by file_names_short
Identification_percentage

  #Making tibble of Identification_percentage, for plot
df_perc <- matrix(unlist(Identification_percentage), nrow=length(Identification_percentage), byrow=TRUE)
tbl_Identification_percentage <- tibble(File_Name=file_names_short , Identification_Percentage=df_perc)
tbl_Identification_percentage
  #Plot
p_perc <- ggplot(tbl_Identification_percentage, aes(x= File_Name , y= Identification_Percentage)) +
  geom_col() +
  labs(x="File Name", y="Percentage Identified Spectra (%)", title="Identification Percentage" , 
        subtitle="Percentage of identified peptides for each file with spectra") +
  geom_text(aes(label=round(df_perc, digits = 4)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 18))
  #Print p_perc
p_perc
pdf(file = "~/Desktop/mzTab/Plots/Identification Percentage.pdf")
   p_perc
dev.off()

  #Extra: create a tibble
cols_old <- colnames(df_perc)
Perc <- c("Identification Percentage")
file_names_6 <- file_names_short[c(1:6)]
df_newnames_perc <- setnames(df_perc, old = cols_old, new = Perc)
tbl_name <- tibble(File_Name = file_names_6)
tbl_perc <- as_tibble(df_newnames_perc)
tbl_id_perc <- cbind(tbl_name, tbl_perc)

#Create column with unmodified peptide sequences
#Make Tibble with PSMs and use PSH as column names
#Select sequence column + Removing brackets + Removing numbers
psms1 <- list()
for (i in seq_along(mzTab_files_PSM)) { #Only for the first 6 spectra
  psms1[[i]] <- as_tibble(mzTab_files_PSM[[i]]) %>%
    row_to_names(row_number = 1) %>%
    select(sequence) %>%
    mutate(sequence_no_mod = trimws(str_remove_all(sequence, "n"))) %>% #remove n
    mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "[0123456789]"))) %>% # remove numbers
    mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "\\[|\\]"))) %>% #remove []
    select(sequence_no_mod)
}
tbl_seq_no_mod <- set_names(psms1, file_names_short)
tbl_seq_no_mod
#Adding new column with unmodified peptide sequences
psms2 <- list()
for (i in 1:6) { #Only for the first 6 spectra
  psms2[[i]] <- as_tibble(mzTab_files_PSM[[i]]) %>%
    row_to_names(row_number = 1) %>%
    add_column(sequence_no_mod = tbl_seq_no_mod[[i]], .after = "sequence")
}
tbl_mzTab_PSM <- set_names(psms2, file_names_short)
tbl_mzTab_PSM
view(tbl_mzTab_PSM[[1]])

#Adding new column with unmodified peptide sequences with CBIND
psms3 <- list()
for (i in seq_along(mzTab_files_PSM)) {
  psms3[[i]] <- as_tibble(mzTab_files_PSM[[i]]) %>%
  row_to_names(row_number = 1)
}
tibbeltje <- set_names(psms3, file_names_short)
tibbeltje

psms4 <- list()
for (i in seq_along(tibbeltje)) { #Only for the first 6 spectra
  psms4[[i]] <- cbind(tibbeltje[[i]][,2] ,tbl_seq_no_mod[[i]], tibbeltje[[i]][,3:ncol(tibbeltje[[i]])]) #Column is all the way in the back
}
 <- set_names(psms4, file_names_short)
view([[1]])

  #Store PSM files
saveRDS(tbl_mzTab_PSM, file = "~/Desktop/mzTab/Stored files/6 PSM")
mzTab_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/6 PSM")
saveRDS(, file = "~/Desktop/mzTab/Stored files/6 PSM")
mzTab_6 <- readRDS(file = "~/Desktop/mzTab/Stored files/6 PSM")

