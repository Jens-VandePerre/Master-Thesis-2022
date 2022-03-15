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

########################################################################################
#Load in outputs directly, not running loops
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

wd <- setwd("~/Desktop/Read raw file/Data mzML")
getwd() 
list.files(wd) #The first 10 mzML files of CPTAC
# Loading in multiple mzML files
file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML")
file_paths #The first 10 mzML files paths of CPTAC