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


wd <-setwd("~/Desktop/Read raw file/Data mzML")
getwd() #the first 4 mzML files of CPTAC


# Loading in multiple mzML files

    #Method 1: for loop 1
file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML")
file_paths #The first 4 mzML files of CPTAC

mzML.4files <- list() #empty list
for (i in seq_along(file_paths)) {
  mzML.4files[[i]] <- readMSData(file = file_paths [[i]], 
                                   msLevel = 2, verbose = FALSE, mode="onDisk")
}

data_names <- c("mzML1", "mzML2", "mzML3", "mzML4")    
mzML.4files <- set_names(mzML.4files, data_names) #names each file
mzML.4files


    #Method 2: purrr map
      #Automatically names the files
      #LONGER LOADING TIME
mzML.4files2 <- file_paths %>%
  map(function(path) {
    readMSData(path)
  })

#Loop extracting TMT intensities
mzML.4files.qnt <- list() #empty list
mzML.4files.qnt

for (i in seq_along(mzML.4files)) {
  mzML.4files.qnt[[i]] <- 
    quantify(mzML.4files[[i]] ,method = "max", #max is the only working method
                                   reporters = TMT10,
                                   strict = FALSE,
                                   verbose = FALSE) %>%
    filterNA(pNA = 0) %>%
    purityCorrect(impurities) %>%
    normalise("max")
}
data_names.qnt <- c("mzML.qnt1", "mzML.qnt2", "mzML.qnt3", "mzML.qnt4")    
mzML.4files.qnt2 <- set_names(mzML.4files.qnt, data_names.qnt) #names each file
mzML.4files.qnt2

#Printing TMT intensities
intensities <- list() #empty list
intensities

for (i in seq_along(mzML.4files.qnt)) {
  intensities[[i]] <-
  head(exprs(mzML.4files.qnt[[i]]))
}
data_names.qnt <- c("mzML.qnt1", "mzML.qnt2", "mzML.qnt3", "mzML.qnt4")    
intensities2 <- set_names(intensities, data_names.qnt) #names each file
intensities2


#Loop extracting TMT intensities + Printing TMT intensities
mzML.4files.qnt2 <- list() #empty list
mzML.4files.qnt2
TMT.intensities <- list() #empty list
TMT.intensities

for (i in seq_along(mzML.4files)) {
  mzML.4files.qnt[[i]] <- 
    quantify(mzML.4files[[i]] ,method = "max", #max is the only working method
             reporters = TMT10,
             strict = FALSE,
             verbose = FALSE) %>%
    filterNA(pNA = 0) %>%
    purityCorrect(impurities) %>%
    normalise("max")
  
  for (j in seq_along(mzML.4files.qnt)) {
    TMT.intensities[[j]] <-
      head(exprs(mzML.4files.qnt[[j]]))
  }
}
data_names.qnt <- c("mzML.qnt1", "mzML.qnt2", "mzML.qnt3", "mzML.qnt4")    
TMT.intensities2 <- set_names(TMT.intensities, data_names.qnt) #names each file
TMT.intensities2


