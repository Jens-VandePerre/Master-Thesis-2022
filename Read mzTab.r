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
list.files(wd) #all mzTabs as of now 

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
#Loop Reading in Files
  #File paths to direct the loop
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/mzTab_19_04_22"))
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_names_short <- substring(file_paths, 91, 98)) #Characters 86 untill 93 are uniqueue
  #Loop reading mzTabs
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
view(mzTab_files[[1]])
 #Save mzMLs to different location: TMT outputs
saveRDS(mzTab_files, file = "~/Desktop/Outputs/mzTabs_imported/20_04_22_mzTabs")
mzTabs_20_04_22 <- readRDS(file = "~/Desktop/Outputs/mzTabs_imported/20_04_22_mzTabs")

##################
#Extract functions
##################

#extractMetadata_long
extractMetadata_long <- function(mztab.table) {
  mztab.table[startsWith(as.character(mztab.table$V1), "MTD"),]
}
  #Loop
MTD_long <- list() #empty list
for (i in seq_along(mzTabs_20_04_22)) {
  MTD_long[[i]] <- extractMetadata_long(mzTabs_20_04_22[[i]])
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
for (i in seq_along(mzTabs_20_04_22)) {
  MTD[[i]] <- extractMetadata(mzTabs_20_04_22[[i]])
}
mzTab_files_Metadata <- set_names(MTD, file_names_short) #names each file by file_names_short
view(mzTab_files_Metadata[[1]])

#extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
  #Loop for all files
PSM <- list() #empty list
for (i in seq_along(mzTabs_20_04_22)) {
  PSM[[i]] <- extractFeaturesPSM(mzTabs_20_04_22[[i]])
}
mzTab_files_PSM <- set_names(PSM, file_names_short) #names each file by file_names_short
view(mzTab_files_PSM[[1]])

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
TMT_Intensities_20_04_22 <- readRDS(file = "~/Desktop/Outputs/TMTs/20.04.22_TMT")
    #Loop counting spectra
nrow_TMT <- list() #empty list
for (i in seq_along(TMT_Intensities_20_04_22)) { #Only for the first 6 spectra, 
  nrow_TMT[[i]] <- n_distinct(TMT_Intensities_20_04_22[[i]])
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
  geom_text(aes(label=round(df_perc, digits = 1)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 18))
  #Print p_perc
p_perc
pdf(file = "~/Desktop/Outputs/Plots/20_04_22_Identification_Percentage.pdf")
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
#Add new column sequence_no_mod behind the sequence column
        #Loop for all files
psm <- list() #empty list
for (i in seq_along(mzTabs_20_04_22)) {
  psm[[i]] <- extractFeaturesPSM(mzTabs_20_04_22[[i]]) %>%
  as_tibble() %>% row_to_names(row_number = 1) %>%
  mutate(sequence_no_mod = trimws(str_remove_all(sequence, "n"))) %>% #remove n
  mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "[0123456789]"))) %>% # remove numbers
  mutate(sequence_no_mod = trimws(str_remove_all(sequence_no_mod, "\\[|\\]"))) %>% #remove []
  relocate(sequence_no_mod, .after = sequence)
}
PSM <- set_names(psm, file_names_short) #names each file by file_names_short
view(PSM[[1]])
  #Store PSM files
saveRDS(PSM, file = "~/Desktop/Outputs/PSMs/20_04_22_PSMs") 
PSM_20_04_22 <- readRDS(file = "~/Desktop/Outputs/PSMs/20_04_22_PSMs")
