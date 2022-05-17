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

#   command line (not needed for running in R)
#args <- commandArgs(trailingOnly = TRUE)
#csv <- args[1]
#bait <- args[2]
#mass_tolerance <- as.numeric(args[3])
#output_path <- paste(bait, "_tol_py.csv", sep = "")

bait <- name
mass_tolerance <- 20
output_path <- "/Users/jensvandeperre/Desktop/Inputs/PTM_mass_differences/PTM_mass_diff_"

#   PSM file (TMTs are included, but not needed)
PSM_TMT_all <- readRDS("~/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
view(PSM_TMT_all[[1]])
str(PSM_TMT_all[[1]])

#File naming 
wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
    #Automate filename extraction
(file_name_long <- substring(list.files(wd), 1, 46))
(file_names_short <- substring(file_name_long, 39, 46))

#   calculate the mass tolerance in Da for each psm
mztab_mtol <- PSM_TMT_all[[1]] %>% as_tibble %>%
    select(sequence, PSM_ID, exp_mass_to_charge, calc_mass_to_charge, charge) %>%
    mutate(exp_mass_to_charge = as.numeric(exp_mass_to_charge)) %>%
    mutate(calc_mass_to_charge = as.numeric(calc_mass_to_charge)) %>%
    mutate(charge = as.numeric(charge)) %>%
    mutate(mass_diff = (exp_mass_to_charge - calc_mass_to_charge)*charge) %>%
    mutate(exp_mass = exp_mass_to_charge * charge) %>%
    mutate(mass_tol = mass_tolerance * exp_mass * 10**-6) %>%   #something goes wrong here
    mutate(mass_tol_pos = mass_diff + mass_tol) %>%
    mutate(mass_tol_neg = mass_diff - mass_tol) %>%
    filter(!(mass_tol_neg < 0 & mass_tol_pos > 0)) %>%      # filter out unmodified psms
    mutate(file_name_long = file_name_long) %>%
    arrange(mass_tol_pos)

view(mztab_mtol)
str(mztab_mtol)

#   save to file for python
for (i in 1:264) {
fwrite(mztab_mtol[[i]], file = paste(output_path, file_name_long[[i]], ".csv", sep = ""), na = "NA", append = FALSE, col.names = TRUE)

}