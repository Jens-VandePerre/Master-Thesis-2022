library(tidyverse)
library(data.table)
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
library("rlist")

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
    #Automate filename extraction
(file_name_long <- substring(list.files(wd), 1, 46))
(file_names_short <- substring(file_name_long, 39, 46)) 

#Select the WHOLE string!!!
text_before <- #from here
'"considerModifications": false,
"createPSMsets": true,
"errorOnNoDecoys": false,
"psmLevelFileID": 0,
"calculateAllFDR": true,
"calculateCombinedFDRScore": true,
"topIdentifications": 0,
"infereProteins": true,
"inferenceMethod": "inference_spectrum_extractor",
"scoringMethod": "scoring_multiplicative",
"scoringBaseScore": "psm_combined_fdr_score",
"scoringPSMs": "best",
"proteinFilters": [
    "nr_unique_peptides_per_protein_filter >= 2"
],
"proteinExportFile": "/Users/jensvandeperre/Desktop/Outputs/PIA_analysis/PIA_output_' #until here

#Filename

text_after <-
'.csv",
"proteinExportWithPSMs": true,
"proteinExportWithPeptides": true,
"proteinExportWithProteinSequences": false`'

#Loop pasting filename in text
parameter_tex <- list()
for (i in 1:264) {
    parameter_tex[[i]] <- 
    paste(text_before, file_name_long[[i]], text_after, sep = "")
}
view(parameter_tex[[1]])
parameter_tex[[20]]
parameter_tex[[264]]

#Loop writing JSON files
for (i in 1:264) {
   write(parameter_tex[[i]], file = paste("/Users/jensvandeperre/Desktop/Inputs/PIA_parameter_files/", file_name_long[[i]], ".json", sep=""))
}
