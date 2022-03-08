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
(file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML"))
file_paths #The first 10 mzML files paths of CPTAC
mzML_4files <- list() #empty list
for (i in seq_along(file_paths)) {
  mzML_4files[[i]] <- readMSData(file_paths [[i]],
                                   msLevel = 2, verbose = FALSE, mode = "onDisk")
}
file_names_wd <- list.files(wd)
file_names_wd
mzML_4files2 <- set_names(mzML_4files, file_names_wd) #names each file by file_names_wd
mzML_4files2

#Loop extracting TMT intensities + Printing TMT intensities
    #Loop that outputs intensities for the 10 first spectra
mzML_4files_qnt <- list() #empty list
TMT_intensities <- list() #empty list
for (i in seq_along(mzML_4files2)) {
  mzML_4files_qnt[[i]] <- 
    quantify(mzML_4files2[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    filterNA(pNA = 0) %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE)) %>%
    normalise("max")
  
  for (j in seq_along(mzML_4files_qnt)) {
    TMT_intensities[[j]] <- head(exprs(mzML_4files_qnt[[j]]),n=10) #only 10 first spectra
}
TMT_intensities2 <- set_names(TMT_intensities, file_names_wd) #names each file by file_names_wd
TMT_intensities2

    #Loop that outputs intensities for all spectra
        #The Terminal output is unclear
mzML_4files_qnt_2 <- list() #empty list
TMT_intensities_2 <- list() #empty list
for (i in seq_along(mzML_4files2)) {
  mzML_4files_qnt_2[[i]] <- 
    quantify(mzML_4files2[[i]], method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    filterNA(pNA = 0) %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE)) %>%
    normalise("max")

  for (j in seq_along(mzML_4files_qnt_2)) {
    TMT_intensities_2[[j]] <- exprs(mzML_4files_qnt_2[[j]]) 
    #all spectra, intensities for each TMT
  }
}
TMT_intensities3 <- set_names(TMT_intensities_2, file_names_wd) #names each file by file_names_wd
TMT_intensities3 #The Terminal output is unclear

########
#Extra's
########
# Loading in multiple mzML files
   #Method 2: purrr map
      #Automatically names the files
      #LONGER LOADING TIME
mzML_4files_purrr <- file_paths %>%
  map(function(path) {
    readMSData(path)
  })
#Two seperate loops for extracting and printing intensities
    #Loop extracting TMT intensities
mzML_4files_qnt <- list() #empty list
for (i in seq_along(mzML_4files)) {
  mzML_4files.qnt[[i]] <- 
    quantify(mzML_4files[[i]] ,method = "max", #max is the only working method
                                   reporters = TMT10,
                                   strict = FALSE,
                                   verbose = FALSE) %>%
    filterNA(pNA = 0) %>%
    purityCorrect(makeImpuritiesMatrix(10, edit = FALSE)) %>%
    normalise("max")
}
data_names_qnt <- c("mzML.qnt1", "mzML.qnt2", "mzML.qnt3", "mzML.qnt4")    
mzML_4files_qnt2 <- set_names(mzML_4files_qnt, file_names_wd) #names each file
mzML_4files_qnt2
    #Printing TMT intensities
intensities <- list() #empty list
for (i in seq_along(mzML_4files_qnt)) {
  intensities[[i]] <-
  head(exprs(mzML_4files_qnt[[i]]))
} 
intensities2 <- set_names(intensities, file_names_wd) #names each file
intensities2
