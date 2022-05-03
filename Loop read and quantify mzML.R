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
library("purrr")

#Run for 39 mzMLs dowloaded 19/04/22
#Working directory with all the wanted files
    #This wd has to contain all the files that have to be analyzed
mzML_WD <- setwd("/Users/jensvandeperre/Desktop/Inputs/mzML_3.5.22")
getwd()
(file_names_wd <- list.files(mzML_WD)) #39 mzML files
(mzML_file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/mzML_3.5.22"))
    #The wanted files
(file_names_short <- substring(fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/mzML_3.5.22"), 88, 95)) #Character 90 untill 97 are unique

#Not needed for this run
  #Wanted mzMLs
      #Listing based on file_names_short
  #listed_mzMLs <- list()
    #for (i in seq_along(file_names_short)) {
      #listed_mzMLs[[i]] <- dir(path = "/Users/jensvandeperre/Desktop/Inputs/mzML_19_04_22", pattern = file_names_short[[i]])
      #}

    #Reading in listed mzMLs
mzML_files <- list() #empty list
for (i in seq_along(mzML_file_paths)) {
  mzML_files[[i]] <- readMSData(mzML_file_paths[[i]], msLevel = 2, verbose = FALSE, mode = "onDisk")
}
(Selected_mzML <- set_names(mzML_files, file_names_short))
 #Save read in mzMLs to output location
saveRDS(Selected_mzML, file = "~/Desktop/Outputs/mzML_imported/03.05.22_mzML")
mzML <- readRDS(file = "~/Desktop/Outputs/mzML_imported/03.05.22_mzML")
view(mzML[[1]])

#Loop extracting TMT intensities + Printing TMT intensities
  #1. Loop that outputs intensities for ALL SPECTRA
      #not removing NAs
      #no imputation
      #no normalisation
mzML_qnt1 <- list() #empty list
TMT1 <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt1[[i]] <- 
    quantify(mzML[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE))
  for (j in seq_along(mzML_qnt1)) {
    TMT1[[j]] <- exprs(mzML_qnt1[[j]]) #outputs all spectra, unclear in terminal
  }
}
TMT_intensities1 <- set_names(TMT1, file_names_short) #names each file by file_names_wd
view(TMT_intensities1[1])
view(TMT1[[1]])
    #Save output to different location: TMT outputs/Combined Files
saveRDS(TMT_intensities1, file = "~/Desktop/Outputs/TMTs/03.05.22_TMT")
TMT_Intensities_28_04_22 <- readRDS(file = "~/Desktop/Outputs/TMTs/03.05.22_TMT")
view(TMT_Intensities_28_04_22[18])
         #1. Check missing data before imputation
missing1 <- list () #empty list
for (i in seq_along(TMT_Intensities_22_04_22)) {
   missing1[[i]] <- sum(is.na(TMT_Intensities_22_04_22[[i]]))}
missing_tot1 <- set_names(missing1, file_names_wd) #names each file by file_names_wd
missing_tot1 # Total missing for each file

######################################################
#Loops Exploring Normalization and Imputation
######################################################
  #2. Loop that outputs intensities for ALL spectra
    #impute: method="MLE"
    #no normalisation
mzML_qnt2 <- list() #empty list
TMT2 <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt2[[i]] <- 
    quantify(mzML[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    impute(method="MLE") %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE))
  for (j in seq_along(mzML_qnt2)) {
    TMT2[[j]] <- exprs(mzML_qnt2[[j]]) #output all spectra, unclear in terminal
  }
}
TMT_intensities2 <- set_names(TMT2, file_names_wd) #names each file by file_names_wd
TMT_intensities2
   #Save output to different location: TMT outputs/Combined Files
saveRDS(TMT_intensities2, file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation")
TMT_Intensities1_10_Imputation <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation")
TMT_Intensities1_10_Imputation
          #2. Check missing data after imputation
missing2 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation)) {
   missing2[[i]] <- sum(is.na(TMT_Intensities1_10_Imputation[[i]]))}
missing_tot2 <- set_names(missing2, file_names_wd) #names each file by file_names_wd
missing_tot2 # Total missing for each file


   #3. Loop that outputs intensities for ALL spectra
          #no imputation
          #normalise: method="center.median"
mzML_qnt3 <- list() #empty list
TMT3 <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt3[[i]] <- 
    quantify(mzML[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE)) %>%
    normalise(method="center.median")
  for (j in seq_along(mzML_qnt3)) {
    TMT3[[j]] <- exprs(mzML_qnt3[[j]]) #output all spectra, unclear in terminal
  }
}
TMT_intensities3 <- set_names(TMT3, file_names_wd) #names each file by file_names_wd
TMT_intensities3 #The Terminal output is unclear
   #Save output to different location: TMT outputs/Combined Files
saveRDS(TMT_intensities3, file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Normalization")
TMT_Intensities1_10_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Normalization")
TMT_Intensities1_10_Normalization  
          #3. Check missing data after normalization
missing3 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Normalization)) {
   missing3[[i]] <- sum(is.na(TMT_intTMT_Intensities1_10_Normalizationensities3[[i]]))}
missing_tot3 <- set_names(missing3, file_names_wd) #names each file by file_names_wd
missing_tot3 # Total missing for each file


    #4. Loop that outputs intensities for ALL spectra
        #Output for later analysis
          #impute: method="MLE"
          #normalise: method="center.median"
mzML_qnt4 <- list() #empty list
TMT4 <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt4[[i]] <- 
    quantify(mzML[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    impute(method="MLE") %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE)) %>%
    normalise(method="center.median")
  for (j in seq_along(mzML_qnt4)) {
    TMT4[[j]] <- exprs(mzML_qnt4[[j]]) #output all spectra, unclear in terminal
  }
}
TMT_intensities4 <- set_names(TMT4, file_names_wd) #names each file by file_names_wd
TMT_intensities4 #The Terminal output is unclear
   #Save output to different location: TMT outputs/Combined Files
saveRDS(TMT_intensities4, file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation+Normalization")
TMT_Intensities1_10_Imputation_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation+Normalization")
TMT_Intensities1_10_Imputation_Normalization
         #4. Check missing data after imputation and normalization
missing4 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation_Normalization)) {
   missing4[[i]] <- sum(is.na(TMT_Intensities1_10_Imputation_Normalization[[i]]))}
missing_tot4 <- set_names(missing4, file_names_wd) #names each file by file_names_wd
missing_tot4 # Total missing for each file



#Difference in missing values after imputation
  #total
diff <- mapply('-', missing_tot1, missing_tot4, SIMPLIFY = FALSE)
missing_diff <- set_names(diff, file_names_wd)
missing_diff
  #percentage lowered missing values
diff_perc <- mapply('/', missing_diff, missing_tot1, SIMPLIFY = FALSE)
diff_perc2 <- mapply('*', diff_perc, 100, SIMPLIFY = FALSE)
missing_diff_perc <- set_names(diff_perc2, file_names_wd)
missing_diff_perc

#####################################################################################
###########################
#Loops Saving output inside
###########################

 #5. Loop that outputs intensities for ALL SPECTRA
      #not removing NAs
      #no imputation
      #no normalisation
mzML_qnt5 <- list() #empty list
TMT5 <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt5[[i]] <- 
    quantify(mzML[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE))
  for (j in seq_along(mzML_qnt5)) {
    TMT5[[j]] <- exprs(mzML_qnt5[[j]]) #outputs all spectra, unclear in terminal
  }
}
TMT_intensities5 <- set_names(TMT5, file_names_wd) #names each file by file_names_wd
    #Save output to different location: TMT outputs
saveRDS(TMT_intensities5, file = "~/Desktop/Read raw file/TMT outputs/TMT_Intensities")

saved_files5 <- list() #empty list
for (i in 1:10) {
  assign(names[i],TMT_intensities5[[i]])
  saveRDS(list=names[i],file = paste(names[i],'~/Desktop/Read raw file/TMT outputs/TMT_Intensities' , sep=""), compress=TRUE)
}

#####################################################################################
########
#Extra's
########
# Loading in multiple mzML files
   #Method 2: purrr map
      #Automatically names the files
      #LONGER LOADING TIME
mzML_purrr <- file_paths %>%
  map(function(path) {
    readMSData(path)
  })
#Two seperate loops for extracting and printing intensities
    #Loop extracting TMT intensities
mzML_qnt_extra <- list() #empty list
for (i in seq_along(mzML)) {
  mzML_qnt_extra[[i]] <- 
    quantify(mzML[[i]] ,method = "max", #max is the only working method
                                   reporters = TMT10,
                                   strict = FALSE,
                                   verbose = FALSE) %>%
    impute(method="MLE") %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE)) %>%
    normalise(method="center.median")
}   
mzML_qnt_extre2 <- set_names(mzML_qnt_extra, file_names_wd) #names each file by file_names_wd
mzML_qnt_extra2
    #Printing TMT intensities
intensities_extra <- list() #empty list
for (i in seq_along(mzML_qnt_extra2)) {
  intensities[[i]] <-
  head(exprs(mzML_qnt_extra2[[i]]))
} 
intensities_extra2 <- set_names(intensities_extra, file_names_wd) #names each file
intensities_extra2
########################################################################################
#TESTING FILES
############################################
#1. Output intensities for ALL SPECTRA
      #not removing NAs
      #no imputation
      #no normalisation
TMT_Intensities1_10 <- readRDS(file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10")
TMT_Intensities1_10 

#2. Output intensities for ALL spectra
    #impute: method="MLE"
    #no normalisation
TMT_Intensities1_10_Imputation <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation")
TMT_Intensities1_10_Imputation

#3. Loop that outputs intensities for ALL spectra
    #no imputation
    #normalise: method="center.median"
TMT_Intensities1_10_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Normalization")
TMT_Intensities1_10_Normalization  

#4. Output intensities for ALL spectra
    #impute: method="MLE"
    #normalise: method="center.median"
TMT_Intensities1_10_Imputation_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation+Normalization")
TMT_Intensities1_10_Imputation_Normalization
########################################################################################
