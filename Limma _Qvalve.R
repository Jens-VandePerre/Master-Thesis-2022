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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/PSM_TMT_mass_diff/Mass_tolerance_10")
getwd() 
list.files(wd)
#Load PSMs 
#Load matching mzMLs
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/PSM_TMT_mass_diff/Mass_tolerance_10"))
(file_names_short <- substring(file_name_long, 55, 62))

#Load in all PSMs
  #PSM files linked to TMTs
PSM_TMT_ALL <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
nrow(PSM_TMT_ALL[[1]])
view(PSM_TMT_ALL[[1]])

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
                      PTM_PSM_TMT <- list()
                      for (i in 1:264) {
                        PTM_PSM_TMT[[i]] <- merge(PSM_TMT[[i]], ptm_index[[i]], by = "index")
                      }
                      view(PTM_PSM_TMT[[1]])
                      length(PTM_PSM_TMT)

                    #Make PTM_PSM_TMT ready for merge to PIA output
                    #Select the needed columns
                    PTM_PSM_TMT_input <- list()
                    for (i in 1:264) {
                      PTM_PSM_TMT_input[[i]] <- PSM_TMT_ALL[[i]] %>% 
                      as_tibble() %>%
                      select(index, sequence, sequence_no_mod, "126":"131", charge) %>%
                      mutate(file_name = file_names_short[[i]]) 
                    }
                    view(PTM_PSM_TMT_input[[1]])
                    view(PTM_PSM_TMT_input[[2]])
                    view(PTM_PSM_TMT_input[[69]])
                    length(PTM_PSM_TMT_input)

#Load protein info from PIA output
(file_paths_PIA <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(file_paths_PIA[[i]], header = FALSE, sep = ",")
}
view(PIA[[1]])
nrow(PIA[[1]])

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
view(PepSeq_ProAcc[[1]])
nrow(PepSeq_ProAcc[[1]])
length(PepSeq_ProAcc)
str(PepSeq_ProAcc[[1]])

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
str(Des[[1]])
view(Des[[1]])
nrow(Des[[1]])
length(Des)

  #Merge: proteins, genes, peptide seq/PSM and TMTs
data <- list()
for (i in 1:264) {
data[[i]] <- merge(PepSeq_ProAcc[[i]], Des[[i]], by = "Accessions") %>%
            as_tibble() %>%
            rename("Protein.Group.Accessions" = Accessions, "Protein.Descriptions" = Description) %>%
            merge(PSM_TMT_ALL[[i]], by = "sequence_no_mod") %>%
            distinct()
}
view(data[[1]])
length(data)

#Start creating input file
mydat <- list()
for (i in 1:264) {
  mydat[[i]] <- data[[i]] %>%
            as_tibble %>%
            select("Protein.Group.Accessions", "Protein.Descriptions", "sequence_no_mod",
            "126":"131") %>%
            rename("Sequence" = sequence_no_mod) %>%
            add_column("Quan.Usage" = "Used", .after = "Sequence") %>%
            add_column("Quan.Info" = "Unique", .after = "Quan.Usage") %>%
            add_column("Isolation.Interference" = 30, .after = "Quan.Info") %>%
            add_column("file_name" = file_names_short[[i]])
}
view(mydat[[2]])
length(mydat)

#Renaming TMTs based on sample type
B1S1_f01_f12 <- mydat[1:12]
B1S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S1_f01_f12_renamed[[i]] <- B1S1_f01_f12[[i]] %>%
  as_tibble() %>%
  rename(NAT_126_B1S1_f01_f12 = "126") %>%
  rename(NAT_127N_B1S1_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B1S1_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B1S1_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B1S1_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B1S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B1S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B1S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S1_f01_f12 = "130C") %>%
  rename(REF_131_B1S1_f01_f12 = "131")
}
B1S2_f01_f12 <- mydat[13:24]
B1S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S2_f01_f12_renamed[[i]] <- B1S2_f01_f12[[i]] %>%
  rename(NAT_126_B1S2_f01_f12 = "126") %>%
  rename(TUMOR_127N_B1S2_f01_f12 = "127N") %>%
  rename(NAT_127C_B1S2_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B1S2_f01_f12 = "128N") %>%
  rename(NAT_128C_B1S2_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B1S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B1S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B1S2_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S2_f01_f12 = "130C") %>%
  rename(REF_131_B1S2_f01_f12 = "131")
}
B1S3_f01_f12 <- mydat[25:36]
B1S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S3_f01_f12_renamed[[i]] <- B1S3_f01_f12[[i]] %>%
  rename(TUMOR_126_B1S3_f01_f12 = "126") %>%
  rename(NAT_127N_B1S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B1S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B1S3_f01_f12 = "128N") %>%
  rename(NAT_128C_B1S3_f01_f12 = "128C") %>%
  rename(NAT_129N_B1S3_f01_f12 = "129N") %>%
  rename(NAT_129C_B1S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B1S3_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S3_f01_f12 = "130C") %>%
  rename(REF_131_B1S3_f01_f12 = "131")
}
B1S4_f01_f12 <- mydat[37:48]
B1S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S4_f01_f12_renamed[[i]] <- B1S4_f01_f12[[i]] %>%
  rename(NAT_161_B1S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B1S4_f01_f12 = "127N") %>%
  rename(NAT_127C_B1S4_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B1S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B1S4_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B1S4_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B1S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B1S4_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S4_f01_f12 = "130C") %>%
  rename(REF_131__B1S4_f01_f12 = "131")
}
B2S1_f01_f12 <- mydat[49:60]
B2S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S1_f01_f12_renamed[[i]] <- B2S1_f01_f12[[i]] %>%
  rename(TUMOR_126_B2S1_f01_f12 = "126") %>%
  rename(TUMOR_127N_B2S1_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B2S1_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S1_f01_f12 = "128N") %>%
  rename(NAT_128C_B2S1_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B2S1_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B2S1_f01_f12 = "130C") %>%
  rename(REF_131_B2S1_f01_f12 = "131")
}
B2S2_f01_f12 <- mydat[61:72]
B2S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S2_f01_f12_renamed[[i]] <- B2S2_f01_f12[[i]] %>%
  rename(TUMOR_126_B2S2_f01_f12 = "126") %>%
  rename(TUMOR_127N_B2S2_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B2S2_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S2_f01_f12 = "128N") %>%
  rename(NAT_128C_B2S2_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B2S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S2_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B2S2_f01_f12 = "130C") %>%
  rename(REF_131_B2S2_f01_f12 = "131")
}
B2S3_f01_f12 <- mydat[73:84]
B2S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S3_f01_f12_renamed[[i]] <- B2S3_f01_f12[[i]] %>%
  rename(TUMOR_126_B2S3_f01_f12 = "126") %>%
  rename(NAT_127N_B2S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B2S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S3_f01_f12 = "128N") %>%
  rename(NAT_128C_B2S3_f01_f12 = "128C") %>%
  rename(NAT_129N_B2S3_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S3_f01_f12 = "130N") %>%
  rename(NAT_130C_B2S3_f01_f12 = "130C") %>%
  rename(REF_131_B2S3_f01_f12 = "131")
}
B2S4_f01_f12 <- mydat[85:96]
B2S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S4_f01_f12_renamed[[i]] <- B2S4_f01_f12[[i]] %>%
  rename(TUMOR_126_B2S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B2S4_f01_f12 = "127N") %>%
  rename(NAT_127C_B2S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B2S4_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B2S4_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S4_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B2S4_f01_f12 = "130C") %>%
  rename(REF_131_B2S4_f01_f12 = "131")
}
B3S1_f01_f12 <- mydat[97:108]
B3S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S1_f01_f12_renamed[[i]] <- B3S1_f01_f12[[i]] %>%
  rename(TUMOR_126_B3S1_f01_f12 = "126") %>%
  rename(NAT_127N_B3S1_f01_f12 = "127N") %>%
  rename(NAT_127C_B3S1_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B3S1_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B3S1_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B3S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B3S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B3S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B3S1_f01_f12 = "130C") %>%
  rename(REF_131_B3S1_f01_f12 = "131")
}
B3S2_f01_f12 <- mydat[109:120]
B3S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S2_f01_f12_renamed[[i]] <- B3S2_f01_f12[[i]] %>%
  rename(NAT_126_B3S2_f01_f12 = "126") %>%
  rename(NAT_127N_B3S2_f01_f12 = "127N") %>%
  rename(NAT_127C_B3S2_f01_f12 = "127C") %>%
  rename(NAT_128N_B3S2_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B3S2_f01_f12 = "128C") %>%
  rename(NAT_129N_B3S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B3S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B3S2_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B3S2_f01_f12 = "130C") %>%
  rename(REF_131_B3S2_f01_f12 = "131")
}
B3S3_f01_f12 <- mydat[121:132]
B3S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S3_f01_f12_renamed[[i]] <- B3S3_f01_f12[[i]] %>%
  rename(TUMOR_126_B3S3_f01_f12 = "126") %>%
  rename(NAT_127N_B3S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B3S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B3S3_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B3S3_f01_f12 = "128C") %>%
  rename(NAT_129N_B3S3_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B3S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B3S3_f01_f12 = "130N") %>%
  rename(NAT_130C_B3S3_f01_f12 = "130C") %>%
  rename(REF_131_B3S3_f01_f12 = "131")
}
B3S4_f01_f12 <- mydat[133:144]
B3S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S4_f01_f12_renamed[[i]] <- B3S4_f01_f12[[i]] %>%
  rename(NAT_126_B3S4_f01_f12 = "126") %>%
  rename(NAT_127N_B3S4_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B3S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B3S4_f01_f12 = "128N") %>%
  rename(NAT_128C_B3S4_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B3S4_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B3S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B3S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B3S4_f01_f12 = "130C") %>%
  rename(REF_131_B3S4_f01_f12 = "131")
}
B4S1_f01_f12 <- mydat[145:156]
B4S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S1_f01_f12_renamed[[i]] <- B4S1_f01_f12[[i]] %>%
  rename(NAT_126_B4S1_f01_f12 = "126") %>%
  rename(TUMOR_127N_B4S1_f01_f12 = "127N") %>%
  rename(NAT_127C_B4S1_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S1_f01_f12 = "128N") %>%
  rename(NAT_128C_B4S1_f01_f12 = "128C") %>%
  rename(NAT_129N_B4S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B4S1_f01_f12 = "129C") %>%
  rename(NAT_130N_B4S1_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S1_f01_f12 = "130C") %>%
  rename(REF_131_B4S1_f01_f12 = "131")
}
B4S2_f01_f12 <- mydat[157:168]
B4S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S2_f01_f12_renamed[[i]] <- B4S2_f01_f12[[i]] %>%
  rename(NAT_126_B4S2_f01_f12 = "126") %>%
  rename(TUMOR_127N_B4S2_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B4S2_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S2_f01_f12 = "128N") %>%
  rename(NAT_128C_B4S2_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B4S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B4S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B4S2_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S2_f01_f12 = "130C") %>%
  rename(REF_131_B4S2_f01_f12 = "131")
}
B4S3_f01_f12 <- mydat[169:180]
B4S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S3_f01_f12_renamed[[i]] <- B4S3_f01_f12[[i]] %>%
  rename(TUMOR_126_B4S3_f01_f12 = "126") %>%
  rename(NAT_127N_B4S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B4S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S3_f01_f12 = "128N") %>%
  rename(NAT_128C_B4S3_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B4S3_f01_f12 = "129N") %>%
  rename(NAT_129C_B4S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B4S3_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S3_f01_f12 = "130C") %>%
  rename(REF_131_B4S3_f01_f12 = "131")
}
B4S4_f01_f12 <- mydat[181:192]
B4S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S4_f01_f12_renamed[[i]] <- B4S4_f01_f12[[i]] %>%
  rename(TUMOR_126_B4S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B4S4_f01_f12 = "127N") %>%
  rename(NAT_127C_B4S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B4S4_f01_f12 = "128C") %>%
  rename(NAT_129N_B4S4_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B4S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B4S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S4_f01_f12 = "130C") %>%
  rename(REF_131_B4S4_f01_f12 = "131")
}
B5S1_f01_f12 <- mydat[193:204]
B5S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S1_f01_f12_renamed[[i]] <- B5S1_f01_f12[[i]] %>%
  rename(NAT_126_B5S1_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S1_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S1_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S1_f01_f12 = "128N") %>%
  rename(NAT_128C_B5S1_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B5S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B5S1_f01_f12 = "130C") %>%
  rename(REF_131_B5S1_f01_f12 = "131")
}
B5S2_f01_f12 <- mydat[205:216]
B5S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S2_f01_f12_renamed[[i]] <- B5S2_f01_f12[[i]] %>%
  rename(NAT_126_B5S2_f01_f12 = "126") %>%
  rename(NAT_127N_B5S2_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S2_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B5S2_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B5S2_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S2_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B5S2_f01_f12 = "130N") %>%
  rename(NAT_130C_B5S2_f01_f12 = "130C") %>%
  rename(REF_131_B5S2_f01_f12 = "131")
}
B5S3_f01_f12 <- mydat[217:228]
B5S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S3_f01_f12_renamed[[i]] <- B5S3_f01_f12[[i]] %>%
  rename(TUMOR_126_B5S3_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S3_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S3_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B5S3_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B5S3_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S3_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S3_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B5S3_f01_f12 = "130C") %>%
  rename(REF_131_B5S3_f01_f12 = "131")
}
B5S4_f01_f12 <- mydat[229:240]
B5S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S4_f01_f12_renamed[[i]] <- B5S4_f01_f12[[i]] %>%
  rename(TUMOR_126_B5S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S4_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B5S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B5S4_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S4_f01_f12 = "129N") %>%
  rename(NAT_129C_B5S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B5S4_f01_f12 = "130C") %>%
  rename(REF_131_B5S4_f01_f12 = "131")
}
B5S5_f01_f12 <- mydat[241:252]
B5S5_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S5_f01_f12_renamed[[i]] <- B5S5_f01_f12[[i]] %>%
  rename(TUMOR_126_B5S5_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S5_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B5S5_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S5_f01_f12 = "128N") %>%
  rename(NAT_128C_B5S5_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B5S5_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S5_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S5_f01_f12 = "130N") %>%
  rename(NAT_130C_B5S5_f01_f12 = "130C") %>%
  rename(REF_131_B5S5_f01_f12 = "131")
}
B5S6_f01_f12 <- mydat[253:264]
B5S6_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S6_f01_f12_renamed[[i]] <- B5S6_f01_f12[[i]] %>%
  rename(TUMOR_126_B5S6_f01_f12 = "126") %>%
  rename(NAT_127N_B5S6_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S6_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S6_f01_f12 = "128N") %>%
  rename(NAT_128C_B5S6_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S6_f01_f12 = "129N") %>%
  rename(NAT_129C_B5S6_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S6_f01_f12 = "130N") %>%
  rename(REF_130C_B5S6_f01_f12 = "130C") %>%
  rename(REF_131_B5S6_f01_f12 = "131")
}


dat <- bind_rows(
  B1S1_f01_f12_renamed,
  B1S2_f01_f12_renamed,
  B1S3_f01_f12_renamed,
  B1S4_f01_f12_renamed,
  B2S1_f01_f12_renamed,
  B2S2_f01_f12_renamed,
  B2S3_f01_f12_renamed,
  B2S4_f01_f12_renamed,
  B3S1_f01_f12_renamed,
  B3S2_f01_f12_renamed,
  B3S3_f01_f12_renamed,
  B3S4_f01_f12_renamed,
  B4S1_f01_f12_renamed,
  B4S2_f01_f12_renamed,
  B4S3_f01_f12_renamed,
  B4S4_f01_f12_renamed,
  B5S1_f01_f12_renamed,
  B5S2_f01_f12_renamed,
  B5S3_f01_f12_renamed,
  B5S4_f01_f12_renamed,
  B5S5_f01_f12_renamed,
  B5S6_f01_f12_renamed
)
view(dat)

  #Here something still need to happen
    #introduce NAs
    dat <- dat %>% mutate_all(na_if,"")
    dat %>% mutate_each(funs(empty_as_na)) 
    #OR only keep proteins present in all rows




cha <- c(
    "NAT_126_B1S1_f01_f12", 
"NAT_127N_B1S1_f01_f12", 
"TUMOR_127C_B1S1_f01_f12", 
"TUMOR_128N_B1S1_f01_f12", 
"TUMOR_128C_B1S1_f01_f12", 
"TUMOR_129N_B1S1_f01_f12", 
"TUMOR_129C_B1S1_f01_f12", 
"TUMOR_130N_B1S1_f01_f12", 
"TUMOR_130C_B1S1_f01_f12", 
"REF_131_B1S1_f01_f12", 
"NAT_126_B1S2_f01_f12", 
"TUMOR_127N_B1S2_f01_f12", 
"NAT_127C_B1S2_f01_f12", 
"TUMOR_128N_B1S2_f01_f12", 
"NAT_128C_B1S2_f01_f12", 
"TUMOR_129N_B1S2_f01_f12", 
"NAT_129C_B1S2_f01_f12", 
"TUMOR_130N_B1S2_f01_f12", 
"TUMOR_130C_B1S2_f01_f12", 
"REF_131_B1S2_f01_f12", 
"TUMOR_126_B1S3_f01_f12", 
"NAT_127N_B1S3_f01_f12", 
"TUMOR_127C_B1S3_f01_f12", 
"NAT_128N_B1S3_f01_f12", 
"NAT_128C_B1S3_f01_f12", 
"NAT_129N_B1S3_f01_f12", 
"NAT_129C_B1S3_f01_f12", 
"TUMOR_130N_B1S3_f01_f12", 
"TUMOR_130C_B1S3_f01_f12", 
"REF_131_B1S3_f01_f12", 
"NAT_161_B1S4_f01_f12", 
"TUMOR_127N_B1S4_f01_f12", 
"NAT_127C_B1S4_f01_f12", 
"TUMOR_128N_B1S4_f01_f12", 
"TUMOR_128C_B1S4_f01_f12", 
"TUMOR_129N_B1S4_f01_f12", 
"TUMOR_129C_B1S4_f01_f12", 
"NAT_130N_B1S4_f01_f12", 
"TUMOR_130C_B1S4_f01_f12", 
"REF_131__B1S4_f01_f12", 
"TUMOR_126_B2S1_f01_f12", 
"TUMOR_127N_B2S1_f01_f12", 
"TUMOR_127C_B2S1_f01_f12", 
"NAT_128N_B2S1_f01_f12", 
"NAT_128C_B2S1_f01_f12", 
"TUMOR_129N_B2S1_f01_f12", 
"NAT_129C_B2S1_f01_f12", 
"TUMOR_130N_B2S1_f01_f12", 
"TUMOR_130C_B2S1_f01_f12", 
"REF_131_B2S1_f01_f12", 
"TUMOR_126_B2S2_f01_f12", 
"TUMOR_127N_B2S2_f01_f12", 
"TUMOR_127C_B2S2_f01_f12", 
"NAT_128N_B2S2_f01_f12", 
"NAT_128C_B2S2_f01_f12", 
"TUMOR_129N_B2S2_f01_f12", 
"NAT_129C_B2S2_f01_f12", 
"TUMOR_130N_B2S2_f01_f12", 
"TUMOR_130C_B2S2_f01_f12", 
"REF_131_B2S2_f01_f12", 
"TUMOR_126_B2S3_f01_f12", 
"NAT_127N_B2S3_f01_f12", 
"TUMOR_127C_B2S3_f01_f12", 
"NAT_128N_B2S3_f01_f12", 
"NAT_128C_B2S3_f01_f12", 
"NAT_129N_B2S3_f01_f12", 
"NAT_129C_B2S3_f01_f12", 
"TUMOR_130N_B2S3_f01_f12", 
"NAT_130C_B2S3_f01_f12", 
"REF_131_B2S3_f01_f12", 
"TUMOR_126_B2S4_f01_f12", 
"TUMOR_127N_B2S4_f01_f12", 
"NAT_127C_B2S4_f01_f12", 
"NAT_128N_B2S4_f01_f12", 
"TUMOR_128C_B2S4_f01_f12", 
"TUMOR_129N_B2S4_f01_f12", 
"NAT_129C_B2S4_f01_f12", 
"TUMOR_130N_B2S4_f01_f12", 
"NAT_130C_B2S4_f01_f12", 
"REF_131_B2S4_f01_f12", 
"TUMOR_126_B3S1_f01_f12", 
"NAT_127N_B3S1_f01_f12", 
"NAT_127C_B3S1_f01_f12", 
"TUMOR_128N_B3S1_f01_f12", 
"TUMOR_128C_B3S1_f01_f12", 
"TUMOR_129N_B3S1_f01_f12", 
"TUMOR_129C_B3S1_f01_f12", 
"TUMOR_130N_B3S1_f01_f12", 
"TUMOR_130C_B3S1_f01_f12", 
"REF_131_B3S1_f01_f12", 
"NAT_126_B3S2_f01_f12", 
"NAT_127N_B3S2_f01_f12", 
"NAT_127C_B3S2_f01_f12", 
"NAT_128N_B3S2_f01_f12", 
"TUMOR_128C_B3S2_f01_f12", 
"NAT_129N_B3S2_f01_f12", 
"NAT_129C_B3S2_f01_f12", 
"TUMOR_130N_B3S2_f01_f12", 
"TUMOR_130C_B3S2_f01_f12", 
"REF_131_B3S2_f01_f12", 
"TUMOR_126_B3S3_f01_f12", 
"NAT_127N_B3S3_f01_f12", 
"TUMOR_127C_B3S3_f01_f12", 
"NAT_128N_B3S3_f01_f12", 
"TUMOR_128C_B3S3_f01_f12", 
"NAT_129N_B3S3_f01_f12", 
"TUMOR_129C_B3S3_f01_f12", 
"TUMOR_130N_B3S3_f01_f12", 
"NAT_130C_B3S3_f01_f12", 
"REF_131_B3S3_f01_f12", 
"NAT_126_B3S4_f01_f12", 
"NAT_127N_B3S4_f01_f12", 
"TUMOR_127C_B3S4_f01_f12", 
"NAT_128N_B3S4_f01_f12", 
"NAT_128C_B3S4_f01_f12", 
"TUMOR_129N_B3S4_f01_f12", 
"TUMOR_129C_B3S4_f01_f12", 
"NAT_130N_B3S4_f01_f12", 
"NAT_130C_B3S4_f01_f12", 
"REF_131_B3S4_f01_f12", 
"NAT_126_B4S1_f01_f12", 
"TUMOR_127N_B4S1_f01_f12", 
"NAT_127C_B4S1_f01_f12", 
"NAT_128N_B4S1_f01_f12", 
"NAT_128C_B4S1_f01_f12", 
"NAT_129N_B4S1_f01_f12", 
"TUMOR_129C_B4S1_f01_f12", 
"NAT_130N_B4S1_f01_f12", 
"NAT_130C_B4S1_f01_f12", 
"REF_131_B4S1_f01_f12", 
"NAT_126_B4S2_f01_f12", 
"TUMOR_127N_B4S2_f01_f12", 
"TUMOR_127C_B4S2_f01_f12", 
"NAT_128N_B4S2_f01_f12", 
"NAT_128C_B4S2_f01_f12", 
"TUMOR_129N_B4S2_f01_f12", 
"NAT_129C_B4S2_f01_f12", 
"TUMOR_130N_B4S2_f01_f12", 
"NAT_130C_B4S2_f01_f12", 
"REF_131_B4S2_f01_f12", 
"TUMOR_126_B4S3_f01_f12", 
"NAT_127N_B4S3_f01_f12", 
"TUMOR_127C_B4S3_f01_f12", 
"NAT_128N_B4S3_f01_f12", 
"NAT_128C_B4S3_f01_f12", 
"TUMOR_129N_B4S3_f01_f12", 
"NAT_129C_B4S3_f01_f12", 
"TUMOR_130N_B4S3_f01_f12", 
"NAT_130C_B4S3_f01_f12", 
"REF_131_B4S3_f01_f12", 
"TUMOR_126_B4S4_f01_f12", 
"TUMOR_127N_B4S4_f01_f12", 
"NAT_127C_B4S4_f01_f12", 
"NAT_128N_B4S4_f01_f12", 
"TUMOR_128C_B4S4_f01_f12", 
"NAT_129N_B4S4_f01_f12", 
"TUMOR_129C_B4S4_f01_f12", 
"NAT_130N_B4S4_f01_f12", 
"NAT_130C_B4S4_f01_f12", 
"REF_131_B4S4_f01_f12", 
"NAT_126_B5S1_f01_f12", 
"TUMOR_127N_B5S1_f01_f12", 
"NAT_127C_B5S1_f01_f12", 
"NAT_128N_B5S1_f01_f12", 
"NAT_128C_B5S1_f01_f12", 
"NAT_129N_B5S1_f01_f12", 
"TUMOR_129C_B5S1_f01_f12", 
"TUMOR_130N_B5S1_f01_f12", 
"TUMOR_130C_B5S1_f01_f12", 
"REF_131_B5S1_f01_f12", 
"NAT_126_B5S2_f01_f12", 
"NAT_127N_B5S2_f01_f12", 
"NAT_127C_B5S2_f01_f12", 
"TUMOR_128N_B5S2_f01_f12", 
"TUMOR_128C_B5S2_f01_f12", 
"NAT_129N_B5S2_f01_f12", 
"TUMOR_129C_B5S2_f01_f12", 
"TUMOR_130N_B5S2_f01_f12", 
"NAT_130C_B5S2_f01_f12", 
"REF_131_B5S2_f01_f12", 
"TUMOR_126_B5S3_f01_f12", 
"TUMOR_127N_B5S3_f01_f12", 
"NAT_127C_B5S3_f01_f12", 
"NAT_128N_B5S3_f01_f12", 
"TUMOR_128C_B5S3_f01_f12", 
"TUMOR_129N_B5S3_f01_f12", 
"TUMOR_129C_B5S3_f01_f12", 
"NAT_130N_B5S3_f01_f12", 
"TUMOR_130C_B5S3_f01_f12", 
"REF_131_B5S3_f01_f12", 
"TUMOR_126_B5S4_f01_f12", 
"TUMOR_126_B5S4_f01_f12", 
"TUMOR_127C_B5S4_f01_f12", 
"NAT_128N_B5S4_f01_f12", 
"TUMOR_128C_B5S4_f01_f12", 
"NAT_129N_B5S4_f01_f12", 
"NAT_129C_B5S4_f01_f12", 
"NAT_130N_B5S4_f01_f12", 
"NAT_130C_B5S4_f01_f12", 
"REF_131_B5S4_f01_f12", 
"TUMOR_126_B5S5_f01_f12", 
"TUMOR_127N_B5S5_f01_f12", 
"TUMOR_127C_B5S5_f01_f12", 
"NAT_128N_B5S5_f01_f12", 
"NAT_128C_B5S5_f01_f12", 
"TUMOR_129N_B5S5_f01_f12", 
"TUMOR_129C_B5S5_f01_f12", 
"NAT_130N_B5S5_f01_f12", 
"NAT_130C_B5S5_f01_f12", 
"REF_131_B5S5_f01_f12", 
"TUMOR_126_B5S6_f01_f12", 
"NAT_127N_B5S6_f01_f12", 
"NAT_127C_B5S6_f01_f12", 
"NAT_128N_B5S6_f01_f12", 
"NAT_128C_B5S6_f01_f12", 
"NAT_129N_B5S6_f01_f12", 
"NAT_129C_B5S6_f01_f12", 
"NAT_130N_B5S6_f01_f12", 
"REF_130C_B5S6_f01_f12", 
"REF_131_B5S6_f01_f12"
)

#Functions to load in
read.peptides <- function(dat, cha){
  output <- NULL
  
  dat$Sequence <- as.character(dat$Sequence)
  dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
  dat$Quan.Usage <- as.character(dat$Quan.Usage)
  dat$Quan.Info <- as.character(dat$Quan.Info)
  dat$Isolation.Interference <- as.numeric(as.character(dat$Isolation.Interference))
  
  dat <- subset(dat, Isolation.Interference<=30)  
  dat <- subset(dat, Quan.Usage=="Used")
  dat <- subset(dat, Protein.Group.Accessions!="")
  dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))  
}

quantify.proteins <- function(dat, cha){
  e.function <- function(x, seq) tapply(x, seq, median)
  output <- NULL
  
  dat$Sequence <- toupper(dat$Sequence) # Capital letters
  accessions <- as.character(unique(dat$Protein.Group.Accessions))
  n.proteins <- length(accessions)
  n.cha <- length(cha)
  
  for(k in 1:n.proteins){
    id <- accessions[k]
    sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
    sdat[cha] <- log2(sdat[cha])
    sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
    pdat <- sdat[, -1]
    n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
    temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])          
    n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)    
    if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
    pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
    output <- rbind(output, pdat)
  }
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- round(output[,1:n.cha],3)
  row.names(output) <- accessions
  output <- as.data.frame(output)
  return(output)
}

eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  log2FC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

eb.fit.mult <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  log2FC <- fit.eb$coef[, "tr2"]
  df.0 <- rep(fit.eb$df.prior, n)
  df.r <- fit.eb$df.residual
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
  t.mod <- fit.eb$t[, "tr2"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "tr2"]
  q.ord <- q.mod <- rep(NA,n)
  ids <- which(!is.na(p.ord))
  k <- 0
  q1 <- qvalue(p.ord[!is.na(p.ord)])$q
  q2 <- qvalue(p.mod[!is.na(p.mod)])$q
  for(i in ids){ 
    k <- k+1
    q.ord[i] <- q1[k] 
    q.mod[i] <- q2[k]
  }
  results.eb.mult <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
  return(results.eb.mult)
}

#Start analysis
(dat <- read.peptides(dat, cha))
dim(dat)
####HIER DELEN??? door ref
dat <- quantify.proteins(dat, cha)
view(dat)
dat.onehit <- subset(dat, dat$n.peptides == 1) 
view(dat)
######OF HIER DELEN??? door ref
dim(dat.onehit)
dat <- subset(dat, dat$n.peptides != 1)
dim(dat)
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat[, 1:length(cha)],  ylim = c(-3, 3), main="Boxplot normalized Intensities")






nor <- c("126", "127N")
tum <- c("127C","128N","128C","129N","129C","130N","130C")

design <- model.matrix(~factor(c(2,2,2,2,2,1,1,1,1)))
design
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(dat[, c(nor,tum)], design)
str(res.eb)
view(res.eb)

rx <- c(-1, 1)*max(abs(res.eb$log2FC))*1.1
ry <- c(0, ceiling(max(-log10(res.eb$p.ord), -log10(res.eb$p.mod))))

par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
par(las=1, xaxs="i", yaxs="i")

plot(res.eb$log2FC, -log10(res.eb$p.ord), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of ordinary p-values")

plot(res.eb$log2FC, -log10(res.eb$p.mod), pch=21, bg="lightgrey", cex=0.9,
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of moderated p-values")



