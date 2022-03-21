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
  #Names for the stored files
file_names_wd <- list.files(wd)
file_names_short_typed <- c("B1S2_f10","B1S3_f05","B1S3_f10","B1S4_f02","B1S4_f06","B2S4_f09")
    #Automate filename extraction
file_names_wd #Long file names
file_names_short <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
  #Loop reading mzTabs
mzTab <- list() #empty list
for (i in seq_along(file_paths)) {
  mzTab[[i]] <- readMzTab(file_paths[[i]])
}
mzTab_files <- set_names(mzTab, file_names_short) #names each file by file_names_short
mzTab_files
  #View first file 
view(mzTab_files[[3]])
view(mzTab_files[[1]] %>% as_tibble())

 #Save mzMLs to different location: TMT outputs
saveRDS(mzTab_files, file = "~/Desktop/mzTab/Stored files/Test 6 mzTabs")
Test_6_mzTabs <- readRDS(file = "~/Desktop/mzTab/Stored files/Test 6 mzTabs")

##################
#Extract functions
##################

#extractMetadata
extractMetadata <- function(mztab.table) {
  mztab.table[startsWith(as.character(mztab.table$V1), "MTD"),]
}
  #For the first file
(mzTab_files_Metadata_First_File <- extractMetadata(mzTab_files[[1]]))
  #Loop for all files
MTD <- list() #empty list
for (i in seq_along(mzTab_files)) {
  MTD[[i]] <- extractMetadata(mzTab_files[[i]])
}
mzTab_files_Metadata <- set_names(MTD, file_names_short) #names each file by file_names_short
mzTab_files

#extractMetadata.V1V2V3
  #Extracting the MTD: only columns first 3 have inforamtion
extractMetadata.V1V2V3 <- function(mztab.table) {
  mztab.table[startsWith(as.character(mztab.table$V1), "MTD"), 1:3]
}
  #For the first file
(mzTab_files_Metadata_V1V2V3_File <- extractMetadata.V1V2V3(mzTab_files[[1]]))
  #Loop for all files
MTD_V1V2V3 <- list() #empty list
for (i in seq_along(mzTab_files)) {
  MTD_V1V2V3[[i]] <- extractMetadata.V1V2V3(mzTab_files[[i]])
}
mzTab_files_Metadata_V1V2V3 <- set_names(MTD_V1V2V3, file_names_short) #names each file by file_names_short
mzTab_files_Metadata_V1V2V3

#extractFeaturesPSM 
extractFeaturesPSM <- function(mztab.table) {
  psm <- mztab.table[startsWith(as.character(mztab.table$V1), "PSM"),]
  psh <- mztab.table[startsWith(as.character(mztab.table$V1), "PSH"),]
  rbind(psh,psm)
}
  #For the first file
(mzTab_files_PSM_First_File <- extractFeaturesPSM(mzTab_files[[1]]))
view(mzTab_files_PSM_First_File)
  #Loop for all files
PSM <- list() #empty list
for (i in seq_along(mzTab_files)) {
  PSM[[i]] <- extractFeaturesPSM(mzTab_files[[i]])
}
mzTab_files_PSM <- set_names(PSM, file_names_short) #names each file by file_names_short
mzTab_files_PSM
view(mzTab_files_PSM[[4]])

#Determine identification percentage
  #Row count for mzTab_PSM = amount of identified peptides
nrow_PSM <- list() #empty list
for (i in seq_along(mzTab_files_PSM)) {
  nrow_PSM[[i]] <- nrow(mzTab_files_PSM[[i]])
}
Identified_peptides <- set_names(nrow_PSM, file_names_short) #names each file by file_names_short
Identified_peptides

  #Row count original
mzML1_10 <- readRDS(file = "~/Desktop/Read raw file/TMT outputs/Combined Files/mzML1-10") #Read in original mzML files
mzML1_10 #Output gives number of spectra
TMT_Intensities1_10 <- readRDS(file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10")
TMT_Intensities1_10 
nrow(TMT_Intensities1_10[[1]]) # Is the same as the number of spectra
    #Loop counting spectra
nrow_TMT <- list() #empty list
for (i in 1:6) { #Only for the first 6 spectra, 
  nrow_TMT[[i]] <- nrow(TMT_Intensities1_10[[i]])
}
Spectral_count <- set_names(nrow_TMT, file_names_short) #names each file by file_names_short
Spectral_count

  #Loop percentage calculation
nrow_TMT <- list() #empty list
for (i in 1:6) { #Only for the first 6 spectra, 
  nrow_TMT[[i]] <- nrow(TMT_Intensities1_10[[i]])
}
Spectral_count <- set_names(nrow_TMT, file_names_short) #names each file by file_names_short
Spectral_count