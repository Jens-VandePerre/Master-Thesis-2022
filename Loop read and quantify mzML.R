library(rpx)
library(mzR)
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
library("limma")
library(isobar)
library(devtools)
library(MSnbase)
library(Biobase)
library(dplyr)
library(tidyverse)
library(fs)
library(proxyC)

wd <- setwd("~/Desktop/Read raw file/Data mzML")
getwd() 
list.files(wd) #The first 10 mzML files of CPTAC
# Loading in multiple mzML files
file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML")
file_paths #The first 10 mzML files paths of CPTAC

mzML_files <- list() #empty list
for (i in seq_along(file_paths)) {
  mzML_files[[i]] <- readMSData(file_paths[[i]],
                                   msLevel = 2, verbose = FALSE, mode = "onDisk")
}
file_names_wd <- list.files(wd)
file_names_wd
mzML <- set_names(mzML_files, file_names_wd) #names each file by file_names_wd
mzML

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
TMT_intensities1 <- set_names(TMT1, file_names_wd) #names each file by file_names_wd
TMT_intensities1
         #1. Check missing data before imputation
missing1 <- list () #empty list
for (i in seq_along(TMT_intensities1)) {
   missing1[[i]] <- sum(is.na(TMT_intensities1[[i]]))}
missing_tot1 <- set_names(missing1, file_names_wd) #names each file by file_names_wd
missing_tot1 # Total missing for each file


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
          #2. Check missing data after imputation
missing2 <- list () #empty list
for (i in seq_along(TMT_intensities2)) {
   missing2[[i]] <- sum(is.na(TMT_intensities2[[i]]))}
missing_tot2 <- set_names(missing2, file_names_wd) #names each file by file_names_wd
missing_tot2 # Total missing for each file

   #3. Loop that outputs intensities for ALL spectra
        #Output for later analysis
          #no imputation
          #normalise: method="center.median"
        #The Terminal output is unclear
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
          #3. Check missing data after imputation
missing3 <- list () #empty list
for (i in seq_along(TMT_intensities3)) {
   missing3[[i]] <- sum(is.na(TMT_intensities3[[i]]))}
missing_tot3 <- set_names(missing3, file_names_wd) #names each file by file_names_wd
missing_tot3 # Total missing for each file


    #4. Loop that outputs intensities for ALL spectra
        #Output for later analysis
          #impute: method="MLE"
          #normalise: method="center.median"
        #The Terminal output is unclear
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
         #4. Check missing data after imputation
missing4 <- list () #empty list
for (i in seq_along(TMT_intensities4)) {
   missing4[[i]] <- sum(is.na(TMT_intensities4[[i]]))}
missing_tot4 <- set_names(missing4, file_names_wd) #names each file by file_names_wd
missing_tot4 # Total missing for each file










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
