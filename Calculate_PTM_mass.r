library("data.table")
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

#   PSM file (TMTs are included, but not needed)
PSM_TMT_all <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
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
    #Loop
mass_tolerance_20 <- 20
mztab_mtol_20 <- list()
for (i in 1:264) {
    mztab_mtol_20[[i]] <- PSM_TMT_all[[i]] %>% 
    as_tibble %>%
    dplyr::select(sequence, PSM_ID, exp_mass_to_charge, calc_mass_to_charge, charge) %>%
    mutate(exp_mass_to_charge = as.numeric(exp_mass_to_charge)) %>%
    mutate(calc_mass_to_charge = as.numeric(calc_mass_to_charge)) %>%
    mutate(charge = as.numeric(charge)) %>%
    mutate(mass_diff = (exp_mass_to_charge - calc_mass_to_charge)*charge) %>%
    mutate(exp_mass = exp_mass_to_charge * charge) %>%
    mutate(mass_tol = mass_tolerance_20 * exp_mass * 10**-6) %>%
    mutate(mass_tol_pos = mass_diff + mass_tol) %>%
    mutate(mass_tol_neg = mass_diff - mass_tol) %>%
    filter(!(mass_tol_neg < 0 & mass_tol_pos > 0)) %>% #filter out unmodified psms
    arrange(mass_tol_pos)
}
view(mztab_mtol_20[[1]])

mass_tolerance_10 <- 10
mztab_mtol_10 <- list()
for (i in 1:264) {
    mztab_mtol_10[[i]] <- PSM_TMT_all[[i]] %>% 
    as_tibble %>%
    dplyr::select(sequence, PSM_ID, exp_mass_to_charge, calc_mass_to_charge, charge) %>%
    mutate(exp_mass_to_charge = as.numeric(exp_mass_to_charge)) %>%
    mutate(calc_mass_to_charge = as.numeric(calc_mass_to_charge)) %>%
    mutate(charge = as.numeric(charge)) %>%
    mutate(mass_diff = (exp_mass_to_charge - calc_mass_to_charge)*charge) %>%
    mutate(exp_mass = exp_mass_to_charge * charge) %>%
    mutate(mass_tol = mass_tolerance_10 * exp_mass * 10**-6) %>%
    mutate(mass_tol_pos = mass_diff + mass_tol) %>%
    mutate(mass_tol_neg = mass_diff - mass_tol) %>%
    filter(!(mass_tol_neg < 0 & mass_tol_pos > 0)) %>% #filter out unmodified psms
    arrange(mass_tol_pos)
}
view(mztab_mtol_10[[1]])

#   save to file for python
    # mass_tolerance = 20
for (i in 1:264) {
    fwrite(mztab_mtol_20[[i]], file = paste("/Users/jensvandeperre/Desktop/Inputs/PTM_mass_differences/Mass_tolerance_20/"
    , file_name_long[[i]], ".csv", sep = ""), na = "NA", append = FALSE, col.names = TRUE)
}
    # mass_tolerance = 10
for (i in 1:264) {
    fwrite(mztab_mtol_10[[i]], file = paste("/Users/jensvandeperre/Desktop/Inputs/PTM_mass_differences/Mass_tolerance_10/"
    , file_name_long[[i]], ".csv", sep = ""), na = "NA", append = FALSE, col.names = TRUE)
}

#Leave all columns in
    #Loop
mass_tolerance_20 <- 20
PSM_20 <- list()
for (i in 1:264) {
    PSM_20[[i]] <- PSM_TMT_all[[i]] %>% 
    as_tibble %>%
    mutate(exp_mass_to_charge = as.numeric(exp_mass_to_charge)) %>%
    mutate(calc_mass_to_charge = as.numeric(calc_mass_to_charge)) %>%
    mutate(charge = as.numeric(charge)) %>%
    mutate(mass_diff = (exp_mass_to_charge - calc_mass_to_charge)*charge) %>%
    mutate(exp_mass = exp_mass_to_charge * charge) %>%
    mutate(mass_tol = mass_tolerance_20 * exp_mass * 10**-6) %>%
    mutate(mass_tol_pos = mass_diff + mass_tol) %>%
    mutate(mass_tol_neg = mass_diff - mass_tol) %>%
    filter(!(mass_tol_neg < 0 & mass_tol_pos > 0)) %>% #filter out unmodified psms
    arrange(mass_tol_pos)
}
view(PSM_20[[1]])
length(PSM_20)

mass_tolerance_10 <- 10
PSM_10 <- list()
for (i in 1:264) {
    PSM_10[[i]] <- PSM_TMT_all[[i]] %>% 
    as_tibble %>%
    mutate(exp_mass_to_charge = as.numeric(exp_mass_to_charge)) %>%
    mutate(calc_mass_to_charge = as.numeric(calc_mass_to_charge)) %>%
    mutate(charge = as.numeric(charge)) %>%
    mutate(mass_diff = (exp_mass_to_charge - calc_mass_to_charge)*charge) %>%
    mutate(exp_mass = exp_mass_to_charge * charge) %>%
    mutate(mass_tol = mass_tolerance_10 * exp_mass * 10**-6) %>%
    mutate(mass_tol_pos = mass_diff + mass_tol) %>%
    mutate(mass_tol_neg = mass_diff - mass_tol) %>%
    filter(!(mass_tol_neg < 0 & mass_tol_pos > 0)) %>% #filter out unmodified psms
    arrange(mass_tol_pos)
}
view(PSMl_10[[1]])
length(PSM_10)


# Save files for Relative Quantification
    # mass_tolerance = 20
for (i in 1:264) {
    fwrite(PSM_20[[i]], file = paste("/Users/jensvandeperre/Desktop/Inputs/PSM_TMT_mass_diff/Mass_tolerance_20/"
    , file_name_long[[i]], ".csv", sep = ""), na = "NA", append = FALSE, col.names = TRUE)
}
    # mass_tolerance = 10
for (i in 1:264) {
    fwrite(PSMl_10[[i]], file = paste("/Users/jensvandeperre/Desktop/Inputs/PSM_TMT_mass_diff/Mass_tolerance_10/"
    , file_name_long[[i]], ".csv", sep = ""), na = "NA", append = FALSE, col.names = TRUE)
}
