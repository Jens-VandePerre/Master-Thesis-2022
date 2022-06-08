library("limma")
library("qvalue")
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
library("rpx")
library("mzR")
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
    #Automate filename extraction
(file_name_long <- substring(list.files(wd), 1, 46))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"))
(file_names_short <- substring(file_name_long, 39, 46))
length(file_names_short)

#Load identified PTMs
(ptm_filepaths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PTM_identification_tol_10"))
PTM <- list()
for (i in 1:264) {
  PTM[[i]] <- read.csv(ptm_filepaths[[i]], sep = ",", header = TRUE)
}
view(PTM[[1]])
nrow(PTM[[1]])
length(PTM)
#Create index for matching to PSM_TMT
ptm_index <- list()
for (i in 1:264) {
  ptm_index[[i]] <- PTM[[i]] %>%
  as_tibble() %>%
  mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
  mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
  mutate(index = trimws(str_remove_all(index, "scan=")))
}
view(ptm_index[[1]])

#merge to PSM_TMT
PSM_TMT <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
nrow(PSM_TMT[[1]])
view(PSM_TMT[[1]])
PTM_PSM_TMT <- list()
for (i in 1:264) {
  PTM_PSM_TMT[[i]] <- merge(PSM_TMT[[i]], ptm_index[[i]], by = "index", all.x = TRUE)
}
view(PTM_PSM_TMT[[1]])
nrow(PTM_PSM_TMT[[1]])

#Add a counter column
new_col <- list()
for (i in seq_along(PTM_PSM_TMT)) {
    new_col[[i]] <- PTM_PSM_TMT[[i]] %>%
    #Selecting the right columns
    select(index, sequence.x, sequence_no_mod, "126":"131", mod, mod_mass) %>%
    #Making new index colmun
    mutate(index_filename = paste0(index, "_",file_names_short[[i]]), .before = index) %>%
    #Adding count column
    add_column(count = rep(1, nrow(PTM_PSM_TMT[[i]])))
}
view(new_col[[1]])
nrow(new_col[[1]])
  #Combine all file into one dataframe
ALL_seq <- bind_rows(new_col)
nrow(ALL_seq)
view(ALL_seq) 
fwrite(ALL_seq, "/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/PSM_TMT_PTM.csv")

#Most prevalant PTMs
Most_PTM <- ALL_seq %>%
    select(index, sequence.x, sequence_no_mod, mod, mod_mass, count) %>%
    group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
nrow(Most_PTM) #1413, same as unique PTMs
view(Most_PTM)
str(Most_PTM)
view(Most_PTM %>% arrange(-Count)) #most frequent mod at top

#Look for Cancer PTMs
CRC_PTMs <- ALL_seq %>%
    select(index, sequence.x, sequence_no_mod, mod, mod_mass, count)
  #CRC
    #Phosphorylation
    #79.966331 monoisotopic
    #79.9799 avg
Phos <- CRC_PTMs %>%
  filter(stringr::str_detect(mod, 'Phospho|phospho')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Phos)

    #Citrullination = Deamidation
    #0.984016 monoisotopic
    #0.9848 avg
Cit <- CRC_PTMs %>%
  filter(stringr::str_detect(mod_mass, '0.98|0.9848')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Cit)

    #Hydroxylation
    #15.994915 monoiso
    #15.9994 avg
Hyd <- CRC_PTMs %>%
  filter(stringr::str_detect(mod_mass, '15.99|15.9994')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Hyd)

  #Cancer in general
    #Ubiquitination
    #383.228103	monoiso
    #383.4460 avg
    #Overblijfsel 2xGly
    #massa 114
    #Vallen deze samen op zelfde peptdie? als volledige ub
Ubq <- CRC_PTMs %>%
  filter(stringr::str_detect(mod_mass, '383.|383.4460')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Ubq)

    #Methylation
    #14.015650 monoiso
    #14.0266 avg
Meth <- CRC_PTMs %>%
  filter(stringr::str_detect(mod_mass, '14.01|14.0266')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Meth)

    #Acetylation
    #42.010565 monoiso
    #42.0367 avg
Acet <- CRC_PTMs %>%
  filter(stringr::str_detect(mod_mass, '42.01|42.0367')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Acet)

    #Glycosilation
    # Overblijfsel
Glyc <- CRC_PTMs %>%
  filter(stringr::str_detect(mod_mass, 'Phospho|phospho')) %>%
  group_by(mod, mod_mass) %>%
    summarise(Count = sum(count))




#How many types of PTMs
  #Count unique PTMs
Unique_PTMs <- ALL_seq %>% 
    select(mod) %>%
    unique()
nrow(Unique_PTMs) #1413
view(Unique_PTMs)

#Most prevalent modified protein
(file_paths_PIA <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(file_paths_PIA[[i]], header = FALSE, sep = ",")
}
  #Merging: Proteins, Peptides and TMTs
    #Select peptide seq and proteins from PIA
PepSeq_ProAcc <- list()
for (i in 1:264) {
PepSeq_ProAcc[[i]] <- PIA[[i]] %>%
          as_tibble() %>%
          filter(str_detect(V1, "PSMSET")) %>%
          mutate(sequence_no_mod = V2) %>%
          mutate(Accessions = V3) %>%
          dplyr::select(sequence_no_mod, Accessions) %>%
          slice(-1) %>%
          as_tibble()
}
      #Select proteins and genes from PIA
Des <- list()
for (i in 1:264) {
Des[[i]] <- PIA[[i]] %>%
          as_tibble() %>%
          filter(str_detect(V1, "PROTEIN")) %>%
          mutate(Accessions = V2) %>%
          mutate(ClusterID = V8) %>%
          mutate(Description = V9) %>%
          dplyr::select(Accessions, ClusterID, Description) %>%
          slice(-1) %>%
          as_tibble()
}
    #Merge: PTMS, proteins, genes, peptide seq/PSM and TMTs
PRO_PTM_PSM_TMT <- list()
for (i in 1:264) {
PRO_PTM_PSM_TMT[[i]] <- merge(PepSeq_ProAcc[[i]], Des[[i]], by = "Accessions") %>%
            as_tibble() %>%
            merge(new_col[[i]], by = "sequence_no_mod", all.y = TRUE) %>%
            distinct() %>%
            as_tibble %>%
            select(Accessions, Description, sequence.x, sequence_no_mod,
            "126":"131", mod, mod_mass, index_filename, count)
}
view(PRO_PTM_PSM_TMT[[1]])
nrow(PRO_PTM_PSM_TMT[[1]])
PRO_PTM_PSM_TMT <- bind_rows(PRO_PTM_PSM_TMT)
fwrite(PRO_PTM_PSM_TMT, "/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/PRO_PTM_PSM_TMT.csv")
nrow(PRO_PTM_PSM_TMT)
      #Most prevalent protein
Most_prev_pro <- PRO_PTM_PSM_TMT %>%
    group_by(Accessions) %>%
    summarise(Count = sum(count))
view(Most_prev_pro)
nrow(Most_prev_pro)

      #Most prevalent modified protein
Most_prev_MOD_pro <- PRO_PTM_PSM_TMT %>%
    group_by(Accessions, mod, mod_mass) %>%
    summarise(Count = sum(count))
view(Most_prev_MOD_pro)
nrow(Most_prev_MOD_pro)


 



