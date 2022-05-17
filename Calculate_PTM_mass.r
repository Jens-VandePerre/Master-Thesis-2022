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


#   command line
#args <- commandArgs(trailingOnly = TRUE)
#csv <- args[1]
#bait <- args[2]
#mass_tolerance <- as.numeric(args[3])
#output_path <- paste(bait, "_tol_py.csv", sep = "")

args <- TRUE
csv <- csv
bait <- name
mass_tolerance <- 20
output_path <- outputpath

#   mztab file
mztab = fread(file=csv, sep = '\t')
PSM_TMT_all <- readRDS("~/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
view(PSM_TMT_all[[1]])
str(PSM_TMT_all[[1]])

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
    mutate(bait = bait) %>%
    arrange(mass_tol_pos)

view(mztab_mtol)
str(mztab_mtol)

#   save to file for python
fwrite(mztab_mtol, output_path, na = "NA", append = FALSE, col.names = TRUE)