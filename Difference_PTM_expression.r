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

#Load proteins
#Load in all PSMs
  #PSM files linked to TMTs
PSM_TMT_ALL <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
dim(PSM_TMT_ALL[[1]])
view(head(PSM_TMT_ALL[[1]]))
  #Load protein info from PIA output
(file_paths_PIA <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(file_paths_PIA[[i]], header = FALSE, sep = ",")
}
view(head(PIA[[1]]))
dim(PIA[[1]])
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
dim(PepSeq_ProAcc[[1]])
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
dim(Des[[1]])
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
view(head(data[[1]]))
length(data)

#Start creating input file
mydat <- list()
for (i in 1:264) {
  mydat[[i]] <- data[[i]] %>%
            as_tibble %>%
            select("Protein.Group.Accessions", "Protein.Descriptions", "sequence",
            "sequence_no_mod","126":"131") %>%
            add_column("file_name" = file_names_short[[i]])
}
view(mydat[[2]])
length(mydat)

saveRDS(mydat, "/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/Input_diff_PTM_analysis")
mydat <- readRDS("/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/Input_diff_PTM_analysis")
view(head(mydat[[1]]))
str(mydat[[1]])
B1S1_f01_f12 <- mydat[1:12]
B1S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S1_f01_f12_renamed[[i]] <- B1S1_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
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
  rename(REF_131_B1S1_f01_f12 = "131") %>%
  drop_na()
}
nrow(B1S1_f01_f12_renamed[[1]])
view(head(B1S1_f01_f12_renamed[[1]]))
B1S2_f01_f12 <- mydat[13:24]
B1S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S2_f01_f12_renamed[[i]] <- B1S2_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B1S2_f01_f12 = "126") %>%
  rename(TUMOR_127N_B1S2_f01_f12 = "127N") %>%
  rename(NAT_127C_B1S2_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B1S2_f01_f12 = "128N") %>%
  rename(NAT_128C_B1S2_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B1S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B1S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B1S2_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S2_f01_f12 = "130C") %>%
  rename(REF_131_B1S2_f01_f12 = "131") %>%
  drop_na()
}
B1S3_f01_f12 <- mydat[25:36]
B1S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S3_f01_f12_renamed[[i]] <- B1S3_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B1S3_f01_f12 = "126") %>%
  rename(NAT_127N_B1S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B1S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B1S3_f01_f12 = "128N") %>%
  rename(NAT_128C_B1S3_f01_f12 = "128C") %>%
  rename(NAT_129N_B1S3_f01_f12 = "129N") %>%
  rename(NAT_129C_B1S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B1S3_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S3_f01_f12 = "130C") %>%
  rename(REF_131_B1S3_f01_f12 = "131") %>%
  drop_na()
}
B1S4_f01_f12 <- mydat[37:48]
B1S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S4_f01_f12_renamed[[i]] <- B1S4_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B1S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B1S4_f01_f12 = "127N") %>%
  rename(NAT_127C_B1S4_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B1S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B1S4_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B1S4_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B1S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B1S4_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B1S4_f01_f12 = "130C") %>%
  rename(REF_131__B1S4_f01_f12 = "131") %>%
  drop_na()
}
B2S1_f01_f12 <- mydat[49:60]
B2S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S1_f01_f12_renamed[[i]] <- B2S1_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B2S1_f01_f12 = "126") %>%
  rename(TUMOR_127N_B2S1_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B2S1_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S1_f01_f12 = "128N") %>%
  rename(NAT_128C_B2S1_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B2S1_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B2S1_f01_f12 = "130C") %>%
  rename(REF_131_B2S1_f01_f12 = "131") %>%
  drop_na()
}
B2S2_f01_f12 <- mydat[61:72]
B2S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S2_f01_f12_renamed[[i]] <- B2S2_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B2S2_f01_f12 = "126") %>%
  rename(TUMOR_127N_B2S2_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B2S2_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S2_f01_f12 = "128N") %>%
  rename(NAT_128C_B2S2_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B2S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S2_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B2S2_f01_f12 = "130C") %>%
  rename(REF_131_B2S2_f01_f12 = "131") %>%
  drop_na()
}
B2S3_f01_f12 <- mydat[73:84]
B2S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S3_f01_f12_renamed[[i]] <- B2S3_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B2S3_f01_f12 = "126") %>%
  rename(NAT_127N_B2S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B2S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S3_f01_f12 = "128N") %>%
  rename(NAT_128C_B2S3_f01_f12 = "128C") %>%
  rename(NAT_129N_B2S3_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S3_f01_f12 = "130N") %>%
  rename(NAT_130C_B2S3_f01_f12 = "130C") %>%
  rename(REF_131_B2S3_f01_f12 = "131") %>%
  drop_na()
}
B2S4_f01_f12 <- mydat[85:96]
B2S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B2S4_f01_f12_renamed[[i]] <- B2S4_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B2S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B2S4_f01_f12 = "127N") %>%
  rename(NAT_127C_B2S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B2S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B2S4_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B2S4_f01_f12 = "129N") %>%
  rename(NAT_129C_B2S4_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B2S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B2S4_f01_f12 = "130C") %>%
  rename(REF_131_B2S4_f01_f12 = "131") %>%
  drop_na()
}
B3S1_f01_f12 <- mydat[97:108]
B3S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S1_f01_f12_renamed[[i]] <- B3S1_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B3S1_f01_f12 = "126") %>%
  rename(NAT_127N_B3S1_f01_f12 = "127N") %>%
  rename(NAT_127C_B3S1_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B3S1_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B3S1_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B3S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B3S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B3S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B3S1_f01_f12 = "130C") %>%
  rename(REF_131_B3S1_f01_f12 = "131") %>%
  drop_na()
}
B3S2_f01_f12 <- mydat[109:120]
B3S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S2_f01_f12_renamed[[i]] <- B3S2_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B3S2_f01_f12 = "126") %>%
  rename(NAT_127N_B3S2_f01_f12 = "127N") %>%
  rename(NAT_127C_B3S2_f01_f12 = "127C") %>%
  rename(NAT_128N_B3S2_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B3S2_f01_f12 = "128C") %>%
  rename(NAT_129N_B3S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B3S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B3S2_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B3S2_f01_f12 = "130C") %>%
  rename(REF_131_B3S2_f01_f12 = "131") %>%
  drop_na()
}
B3S3_f01_f12 <- mydat[121:132]
B3S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S3_f01_f12_renamed[[i]] <- B3S3_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B3S3_f01_f12 = "126") %>%
  rename(NAT_127N_B3S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B3S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B3S3_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B3S3_f01_f12 = "128C") %>%
  rename(NAT_129N_B3S3_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B3S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B3S3_f01_f12 = "130N") %>%
  rename(NAT_130C_B3S3_f01_f12 = "130C") %>%
  rename(REF_131_B3S3_f01_f12 = "131") %>%
  drop_na()
}
B3S4_f01_f12 <- mydat[133:144]
B3S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B3S4_f01_f12_renamed[[i]] <- B3S4_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B3S4_f01_f12 = "126") %>%
  rename(NAT_127N_B3S4_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B3S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B3S4_f01_f12 = "128N") %>%
  rename(NAT_128C_B3S4_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B3S4_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B3S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B3S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B3S4_f01_f12 = "130C") %>%
  rename(REF_131_B3S4_f01_f12 = "131") %>%
  drop_na()
}
B4S1_f01_f12 <- mydat[145:156]
B4S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S1_f01_f12_renamed[[i]] <- B4S1_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B4S1_f01_f12 = "126") %>%
  rename(TUMOR_127N_B4S1_f01_f12 = "127N") %>%
  rename(NAT_127C_B4S1_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S1_f01_f12 = "128N") %>%
  rename(NAT_128C_B4S1_f01_f12 = "128C") %>%
  rename(NAT_129N_B4S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B4S1_f01_f12 = "129C") %>%
  rename(NAT_130N_B4S1_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S1_f01_f12 = "130C") %>%
  rename(REF_131_B4S1_f01_f12 = "131") %>%
  drop_na()
}
B4S2_f01_f12 <- mydat[157:168]
B4S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S2_f01_f12_renamed[[i]] <- B4S2_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B4S2_f01_f12 = "126") %>%
  rename(TUMOR_127N_B4S2_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B4S2_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S2_f01_f12 = "128N") %>%
  rename(NAT_128C_B4S2_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B4S2_f01_f12 = "129N") %>%
  rename(NAT_129C_B4S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B4S2_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S2_f01_f12 = "130C") %>%
  rename(REF_131_B4S2_f01_f12 = "131") %>%
  drop_na()
}
B4S3_f01_f12 <- mydat[169:180]
B4S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S3_f01_f12_renamed[[i]] <- B4S3_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B4S3_f01_f12 = "126") %>%
  rename(NAT_127N_B4S3_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B4S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S3_f01_f12 = "128N") %>%
  rename(NAT_128C_B4S3_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B4S3_f01_f12 = "129N") %>%
  rename(NAT_129C_B4S3_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B4S3_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S3_f01_f12 = "130C") %>%
  rename(REF_131_B4S3_f01_f12 = "131") %>%
  drop_na()
}
B4S4_f01_f12 <- mydat[181:192]
B4S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B4S4_f01_f12_renamed[[i]] <- B4S4_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B4S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B4S4_f01_f12 = "127N") %>%
  rename(NAT_127C_B4S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B4S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B4S4_f01_f12 = "128C") %>%
  rename(NAT_129N_B4S4_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B4S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B4S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B4S4_f01_f12 = "130C") %>%
  rename(REF_131_B4S4_f01_f12 = "131") %>%
  drop_na()
}
B5S1_f01_f12 <- mydat[193:204]
B5S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S1_f01_f12_renamed[[i]] <- B5S1_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B5S1_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S1_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S1_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S1_f01_f12 = "128N") %>%
  rename(NAT_128C_B5S1_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S1_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S1_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B5S1_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B5S1_f01_f12 = "130C") %>%
  rename(REF_131_B5S1_f01_f12 = "131") %>%
  drop_na()
}
B5S2_f01_f12 <- mydat[205:216]
B5S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S2_f01_f12_renamed[[i]] <- B5S2_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_126_B5S2_f01_f12 = "126") %>%
  rename(NAT_127N_B5S2_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S2_f01_f12 = "127C") %>%
  rename(TUMOR_128N_B5S2_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B5S2_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S2_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S2_f01_f12 = "129C") %>%
  rename(TUMOR_130N_B5S2_f01_f12 = "130N") %>%
  rename(NAT_130C_B5S2_f01_f12 = "130C") %>%
  rename(REF_131_B5S2_f01_f12 = "131") %>%
  drop_na()
}
B5S3_f01_f12 <- mydat[217:228]
B5S3_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S3_f01_f12_renamed[[i]] <- B5S3_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B5S3_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S3_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S3_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S3_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B5S3_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B5S3_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S3_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S3_f01_f12 = "130N") %>%
  rename(TUMOR_130C_B5S3_f01_f12 = "130C") %>%
  rename(REF_131_B5S3_f01_f12 = "131") %>%
  drop_na()
}
B5S4_f01_f12 <- mydat[229:240]
B5S4_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S4_f01_f12_renamed[[i]] <- B5S4_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B5S4_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S4_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B5S4_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S4_f01_f12 = "128N") %>%
  rename(TUMOR_128C_B5S4_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S4_f01_f12 = "129N") %>%
  rename(NAT_129C_B5S4_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S4_f01_f12 = "130N") %>%
  rename(NAT_130C_B5S4_f01_f12 = "130C") %>%
  rename(REF_131_B5S4_f01_f12 = "131") %>%
  drop_na()
}
B5S5_f01_f12 <- mydat[241:252]
B5S5_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S5_f01_f12_renamed[[i]] <- B5S5_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B5S5_f01_f12 = "126") %>%
  rename(TUMOR_127N_B5S5_f01_f12 = "127N") %>%
  rename(TUMOR_127C_B5S5_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S5_f01_f12 = "128N") %>%
  rename(NAT_128C_B5S5_f01_f12 = "128C") %>%
  rename(TUMOR_129N_B5S5_f01_f12 = "129N") %>%
  rename(TUMOR_129C_B5S5_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S5_f01_f12 = "130N") %>%
  rename(NAT_130C_B5S5_f01_f12 = "130C") %>%
  rename(REF_131_B5S5_f01_f12 = "131") %>%
  drop_na()
}
B5S6_f01_f12 <- mydat[253:264]
B5S6_f01_f12_renamed <- list()
for (i in 1:12) {
  B5S6_f01_f12_renamed[[i]] <- B5S6_f01_f12[[i]] %>%
  select(file_name, Protein.Group.Accessions, Protein.Descriptions, sequence, sequence_no_mod, "126":"131") %>%
  as_tibble() %>%
  rename(TUMOR_126_B5S6_f01_f12 = "126") %>%
  rename(NAT_127N_B5S6_f01_f12 = "127N") %>%
  rename(NAT_127C_B5S6_f01_f12 = "127C") %>%
  rename(NAT_128N_B5S6_f01_f12 = "128N") %>%
  rename(NAT_128C_B5S6_f01_f12 = "128C") %>%
  rename(NAT_129N_B5S6_f01_f12 = "129N") %>%
  rename(NAT_129C_B5S6_f01_f12 = "129C") %>%
  rename(NAT_130N_B5S6_f01_f12 = "130N") %>%
  rename(REF_130C_B5S6_f01_f12 = "130C") %>%
  rename(REF_131_B5S6_f01_f12 = "131") %>%
  drop_na()
}

#Calculate relative TMT intensities
dat_B1S1_f01_f12 <- bind_rows(B1S1_f01_f12_renamed) %>%
  group_by(Protein.Group.Accessions, sequence) %>% 
  summarise(NAT_126_B1S1_f01_f12 = sum(NAT_126_B1S1_f01_f12), 
  NAT_127N_B1S1_f01_f12 = sum(NAT_127N_B1S1_f01_f12), 
  TUMOR_127C_B1S1_f01_f12 = sum(TUMOR_127C_B1S1_f01_f12),
  TUMOR_128N_B1S1_f01_f12 =sum(TUMOR_128N_B1S1_f01_f12),
  TUMOR_128C_B1S1_f01_f12 =sum(TUMOR_128C_B1S1_f01_f12),
  TUMOR_129N_B1S1_f01_f12 =sum(TUMOR_129N_B1S1_f01_f12),
  TUMOR_129C_B1S1_f01_f12 =sum(TUMOR_129C_B1S1_f01_f12),
  TUMOR_130N_B1S1_f01_f12 =sum(TUMOR_130N_B1S1_f01_f12),
  TUMOR_130C_B1S1_f01_f12 =sum(TUMOR_130C_B1S1_f01_f12),
  REF_131_B1S1_f01_f12 =sum(REF_131_B1S1_f01_f12)
  ) %>%
  mutate (NAT_126_B1S1_f01_f12 = NAT_126_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (NAT_127N_B1S1_f01_f12 = NAT_127N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_127C_B1S1_f01_f12= TUMOR_127C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_128N_B1S1_f01_f12 = TUMOR_128N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_128C_B1S1_f01_f12 = TUMOR_128C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_129N_B1S1_f01_f12 = TUMOR_129N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_129C_B1S1_f01_f12 = TUMOR_129C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_130N_B1S1_f01_f12 = TUMOR_130N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_130C_B1S1_f01_f12 = TUMOR_130C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  select(-REF_131_B1S1_f01_f12)
dim(dat_B1S1_f01_f12) 
view(head(dat_B1S1_f01_f12))

dat_B1S2_f01_f12 <- bind_rows(B1S2_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B1S2_f01_f12 = sum(NAT_126_B1S2_f01_f12),
        TUMOR_127N_B1S2_f01_f12 = sum(TUMOR_127N_B1S2_f01_f12),
        NAT_127C_B1S2_f01_f12 = sum(NAT_127C_B1S2_f01_f12),
        TUMOR_128N_B1S2_f01_f12 = sum(TUMOR_128N_B1S2_f01_f12),
        NAT_128C_B1S2_f01_f12 = sum(NAT_128C_B1S2_f01_f12),
        TUMOR_129N_B1S2_f01_f12 = sum(TUMOR_129N_B1S2_f01_f12),
        NAT_129C_B1S2_f01_f12 = sum(NAT_129C_B1S2_f01_f12),
        TUMOR_130N_B1S2_f01_f12 = sum(TUMOR_130N_B1S2_f01_f12),
        TUMOR_130C_B1S2_f01_f12 = sum(TUMOR_130C_B1S2_f01_f12),
        REF_131_B1S2_f01_f12 = sum(REF_131_B1S2_f01_f12)
    ) %>%
  mutate(NAT_126_B1S2_f01_f12 = NAT_126_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_127N_B1S2_f01_f12 = TUMOR_127N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(NAT_127C_B1S2_f01_f12 = NAT_127C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_128N_B1S2_f01_f12 = TUMOR_128N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(NAT_128C_B1S2_f01_f12 = NAT_128C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_129N_B1S2_f01_f12 = TUMOR_129N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(NAT_129C_B1S2_f01_f12 = NAT_129C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_130N_B1S2_f01_f12 = TUMOR_130N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_130C_B1S2_f01_f12 = TUMOR_130C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  select(-REF_131_B1S2_f01_f12) 

dat_B1S3_f01_f12 <- bind_rows(B1S3_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B1S3_f01_f12 = sum(TUMOR_126_B1S3_f01_f12),
        NAT_127N_B1S3_f01_f12 = sum(NAT_127N_B1S3_f01_f12),
        TUMOR_127C_B1S3_f01_f12 = sum(TUMOR_127C_B1S3_f01_f12),
        NAT_128N_B1S3_f01_f12 = sum(NAT_128N_B1S3_f01_f12),
        NAT_128C_B1S3_f01_f12 = sum(NAT_128C_B1S3_f01_f12),
        NAT_129N_B1S3_f01_f12 = sum(NAT_129N_B1S3_f01_f12),
        NAT_129C_B1S3_f01_f12 = sum(NAT_129C_B1S3_f01_f12),
        TUMOR_130N_B1S3_f01_f12 = sum(TUMOR_130N_B1S3_f01_f12),
        TUMOR_130C_B1S3_f01_f12 = sum(TUMOR_130C_B1S3_f01_f12),
        REF_131_B1S3_f01_f12 = sum(REF_131_B1S3_f01_f12)
    ) %>%
  mutate(TUMOR_126_B1S3_f01_f12 = TUMOR_126_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_127N_B1S3_f01_f12 = NAT_127N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(TUMOR_127C_B1S3_f01_f12 = TUMOR_127C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_128N_B1S3_f01_f12 = NAT_128N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_128C_B1S3_f01_f12 = NAT_128C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_129N_B1S3_f01_f12 = NAT_129N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_129C_B1S3_f01_f12 = NAT_129C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(TUMOR_130N_B1S3_f01_f12 = TUMOR_130N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(TUMOR_130C_B1S3_f01_f12 = TUMOR_130C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  select(-REF_131_B1S3_f01_f12) 

dat_B1S4_f01_f12 <- bind_rows(B1S4_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B1S4_f01_f12 = sum(NAT_126_B1S4_f01_f12),
        TUMOR_127N_B1S4_f01_f12 = sum(TUMOR_127N_B1S4_f01_f12),
        NAT_127C_B1S4_f01_f12 = sum(NAT_127C_B1S4_f01_f12),
        TUMOR_128N_B1S4_f01_f12 = sum(TUMOR_128N_B1S4_f01_f12),
        TUMOR_128C_B1S4_f01_f12 = sum(TUMOR_128C_B1S4_f01_f12),
        TUMOR_129N_B1S4_f01_f12 = sum(TUMOR_129N_B1S4_f01_f12),
        TUMOR_129C_B1S4_f01_f12 = sum(TUMOR_129C_B1S4_f01_f12),
        NAT_130N_B1S4_f01_f12 = sum(NAT_130N_B1S4_f01_f12),
        TUMOR_130C_B1S4_f01_f12 = sum(TUMOR_130C_B1S4_f01_f12),
        REF_131__B1S4_f01_f12 =sum(REF_131__B1S4_f01_f12)
    ) %>%
  mutate(NAT_126_B1S4_f01_f12 = NAT_126_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_127N_B1S4_f01_f12 = TUMOR_127N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(NAT_127C_B1S4_f01_f12 = NAT_127C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_128N_B1S4_f01_f12 = TUMOR_128N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_128C_B1S4_f01_f12 = TUMOR_128C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_129N_B1S4_f01_f12 = TUMOR_129N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_129C_B1S4_f01_f12 = TUMOR_129C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(NAT_130N_B1S4_f01_f12 = NAT_130N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_130C_B1S4_f01_f12 = TUMOR_130C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  select(-REF_131__B1S4_f01_f12) 

dat_B2S1_f01_f12 <- bind_rows(B2S1_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B2S1_f01_f12 = sum(TUMOR_126_B2S1_f01_f12),
        TUMOR_127N_B2S1_f01_f12 = sum(TUMOR_127N_B2S1_f01_f12),
        TUMOR_127C_B2S1_f01_f12 = sum(TUMOR_127C_B2S1_f01_f12),
        NAT_128N_B2S1_f01_f12 = sum(NAT_128N_B2S1_f01_f12),
        NAT_128C_B2S1_f01_f12 = sum(NAT_128C_B2S1_f01_f12),
        TUMOR_129N_B2S1_f01_f12 = sum(TUMOR_129N_B2S1_f01_f12),
        NAT_129C_B2S1_f01_f12 = sum(NAT_129C_B2S1_f01_f12),
        TUMOR_130N_B2S1_f01_f12 = sum(TUMOR_130N_B2S1_f01_f12),
        TUMOR_130C_B2S1_f01_f12 = sum(TUMOR_130C_B2S1_f01_f12),
        REF_131_B2S1_f01_f12 = sum(REF_131_B2S1_f01_f12)        
    ) %>%
  mutate(TUMOR_126_B2S1_f01_f12 = TUMOR_126_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_127N_B2S1_f01_f12 = TUMOR_127N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_127C_B2S1_f01_f12 = TUMOR_127C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(NAT_128N_B2S1_f01_f12 = NAT_128N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(NAT_128C_B2S1_f01_f12 = NAT_128C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_129N_B2S1_f01_f12 = TUMOR_129N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(NAT_129C_B2S1_f01_f12 = NAT_129C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_130N_B2S1_f01_f12 = TUMOR_130N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_130C_B2S1_f01_f12 = TUMOR_130C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  select(-REF_131_B2S1_f01_f12) 

dat_B2S2_f01_f12 <- bind_rows(B2S2_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B2S2_f01_f12 = sum(TUMOR_126_B2S2_f01_f12),
        TUMOR_127N_B2S2_f01_f12 = sum(TUMOR_127N_B2S2_f01_f12),
        TUMOR_127C_B2S2_f01_f12 = sum(TUMOR_127C_B2S2_f01_f12),
        NAT_128N_B2S2_f01_f12 = sum(NAT_128N_B2S2_f01_f12),
        NAT_128C_B2S2_f01_f12 = sum(NAT_128C_B2S2_f01_f12),
        TUMOR_129N_B2S2_f01_f12 = sum(TUMOR_129N_B2S2_f01_f12),
        NAT_129C_B2S2_f01_f12 = sum(NAT_129C_B2S2_f01_f12),
        TUMOR_130N_B2S2_f01_f12 = sum(TUMOR_130N_B2S2_f01_f12),
        TUMOR_130C_B2S2_f01_f12 = sum(TUMOR_130C_B2S2_f01_f12),
        REF_131_B2S2_f01_f12 =sum(REF_131_B2S2_f01_f12)
    ) %>%
  mutate(TUMOR_126_B2S2_f01_f12 = TUMOR_126_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_127N_B2S2_f01_f12 = TUMOR_127N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_127C_B2S2_f01_f12 = TUMOR_127C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(NAT_128N_B2S2_f01_f12 = NAT_128N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(NAT_128C_B2S2_f01_f12 = NAT_128C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_129N_B2S2_f01_f12 = TUMOR_129N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(NAT_129C_B2S2_f01_f12 = NAT_129C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_130N_B2S2_f01_f12 = TUMOR_130N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_130C_B2S2_f01_f12 = TUMOR_130C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  select(-REF_131_B2S2_f01_f12) 

dat_B2S3_f01_f12 <- bind_rows(B2S3_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B2S3_f01_f12 = sum(TUMOR_126_B2S3_f01_f12),
        NAT_127N_B2S3_f01_f12 = sum(NAT_127N_B2S3_f01_f12),
        TUMOR_127C_B2S3_f01_f12 = sum(TUMOR_127C_B2S3_f01_f12),
        NAT_128N_B2S3_f01_f12 = sum(NAT_128N_B2S3_f01_f12),
        NAT_128C_B2S3_f01_f12 = sum(NAT_128C_B2S3_f01_f12),
        NAT_129N_B2S3_f01_f12 = sum(NAT_129N_B2S3_f01_f12),
        NAT_129C_B2S3_f01_f12 = sum(NAT_129C_B2S3_f01_f12),
        TUMOR_130N_B2S3_f01_f12 = sum(TUMOR_130N_B2S3_f01_f12),
        NAT_130C_B2S3_f01_f12 = sum(NAT_130C_B2S3_f01_f12),
        REF_131_B2S3_f01_f12 = sum(REF_131_B2S3_f01_f12)
    ) %>%
  mutate(TUMOR_126_B2S3_f01_f12 = TUMOR_126_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_127N_B2S3_f01_f12 = NAT_127N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(TUMOR_127C_B2S3_f01_f12 = TUMOR_127C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_128N_B2S3_f01_f12 = NAT_128N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_128C_B2S3_f01_f12 = NAT_128C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_129N_B2S3_f01_f12 = NAT_129N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_129C_B2S3_f01_f12 = NAT_129C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(TUMOR_130N_B2S3_f01_f12 = TUMOR_130N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_130C_B2S3_f01_f12 = NAT_130C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  select(-REF_131_B2S3_f01_f12) 

dat_B2S4_f01_f12 <- bind_rows(B2S4_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B2S4_f01_f12 = sum(TUMOR_126_B2S4_f01_f12),
        TUMOR_127N_B2S4_f01_f12 = sum(TUMOR_127N_B2S4_f01_f12),
        NAT_127C_B2S4_f01_f12 = sum(NAT_127C_B2S4_f01_f12),
        NAT_128N_B2S4_f01_f12 = sum(NAT_128N_B2S4_f01_f12),
        TUMOR_128C_B2S4_f01_f12 = sum(TUMOR_128C_B2S4_f01_f12),
        TUMOR_129N_B2S4_f01_f12 = sum(TUMOR_129N_B2S4_f01_f12),
        NAT_129C_B2S4_f01_f12 = sum(NAT_129C_B2S4_f01_f12),
        TUMOR_130N_B2S4_f01_f12 = sum(TUMOR_130N_B2S4_f01_f12),
        NAT_130C_B2S4_f01_f12 = sum(NAT_130C_B2S4_f01_f12),
        REF_131_B2S4_f01_f12 = sum(REF_131_B2S4_f01_f12)
    ) %>%
  mutate(TUMOR_126_B2S4_f01_f12 = TUMOR_126_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_127N_B2S4_f01_f12 = TUMOR_127N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_127C_B2S4_f01_f12 = NAT_127C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_128N_B2S4_f01_f12 = NAT_128N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_128C_B2S4_f01_f12 = TUMOR_128C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_129N_B2S4_f01_f12 = TUMOR_129N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_129C_B2S4_f01_f12 = NAT_129C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_130N_B2S4_f01_f12 = TUMOR_130N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_130C_B2S4_f01_f12 = NAT_130C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  select(-REF_131_B2S4_f01_f12) 

dat_B3S1_f01_f12 <- bind_rows(B3S1_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B3S1_f01_f12 = sum(TUMOR_126_B3S1_f01_f12),
        NAT_127N_B3S1_f01_f12 = sum(NAT_127N_B3S1_f01_f12),
        NAT_127C_B3S1_f01_f12 = sum(NAT_127C_B3S1_f01_f12),
        TUMOR_128N_B3S1_f01_f12 = sum(TUMOR_128N_B3S1_f01_f12),
        TUMOR_128C_B3S1_f01_f12 = sum(TUMOR_128C_B3S1_f01_f12),
        TUMOR_129N_B3S1_f01_f12 = sum(TUMOR_129N_B3S1_f01_f12),
        TUMOR_129C_B3S1_f01_f12 = sum(TUMOR_129C_B3S1_f01_f12),
        TUMOR_130N_B3S1_f01_f12 = sum(TUMOR_130N_B3S1_f01_f12),
        TUMOR_130C_B3S1_f01_f12 = sum(TUMOR_130C_B3S1_f01_f12),
        REF_131_B3S1_f01_f12 = sum(REF_131_B3S1_f01_f12)
    ) %>%
  mutate(TUMOR_126_B3S1_f01_f12 = TUMOR_126_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(NAT_127N_B3S1_f01_f12 = NAT_127N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(NAT_127C_B3S1_f01_f12 = NAT_127C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_128N_B3S1_f01_f12 = TUMOR_128N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_128C_B3S1_f01_f12 = TUMOR_128C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_129N_B3S1_f01_f12 = TUMOR_129N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_129C_B3S1_f01_f12 = TUMOR_129C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_130N_B3S1_f01_f12 = TUMOR_130N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_130C_B3S1_f01_f12 = TUMOR_130C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  select(-REF_131_B3S1_f01_f12) 

dat_B3S2_f01_f12 <- bind_rows(B3S2_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B3S2_f01_f12 = sum(NAT_126_B3S2_f01_f12),
        NAT_127N_B3S2_f01_f12 = sum(NAT_127N_B3S2_f01_f12),
        NAT_127C_B3S2_f01_f12 = sum(NAT_127C_B3S2_f01_f12),
        NAT_128N_B3S2_f01_f12 = sum(NAT_128N_B3S2_f01_f12),
        TUMOR_128C_B3S2_f01_f12 = sum(TUMOR_128C_B3S2_f01_f12),
        NAT_129N_B3S2_f01_f12 = sum(NAT_129N_B3S2_f01_f12),
        NAT_129C_B3S2_f01_f12 = sum(NAT_129C_B3S2_f01_f12),
        TUMOR_130N_B3S2_f01_f12 = sum(TUMOR_130N_B3S2_f01_f12),
        TUMOR_130C_B3S2_f01_f12 = sum(TUMOR_130C_B3S2_f01_f12),
        REF_131_B3S2_f01_f12 = sum(REF_131_B3S2_f01_f12)
    ) %>%
  mutate(NAT_126_B3S2_f01_f12 = NAT_126_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_127N_B3S2_f01_f12 = NAT_127N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_127C_B3S2_f01_f12 = NAT_127C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_128N_B3S2_f01_f12 = NAT_128N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(TUMOR_128C_B3S2_f01_f12 = TUMOR_128C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_129N_B3S2_f01_f12 = NAT_129N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_129C_B3S2_f01_f12 = NAT_129C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(TUMOR_130N_B3S2_f01_f12 = TUMOR_130N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(TUMOR_130C_B3S2_f01_f12 = TUMOR_130C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  select(-REF_131_B3S2_f01_f12) 

dat_B3S3_f01_f12 <- bind_rows(B3S3_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B3S3_f01_f12 = sum(TUMOR_126_B3S3_f01_f12),
        NAT_127N_B3S3_f01_f12 = sum(NAT_127N_B3S3_f01_f12),
        TUMOR_127C_B3S3_f01_f12 = sum(TUMOR_127C_B3S3_f01_f12),
        NAT_128N_B3S3_f01_f12 = sum(NAT_128N_B3S3_f01_f12),
        TUMOR_128C_B3S3_f01_f12 = sum(TUMOR_128C_B3S3_f01_f12),
        NAT_129N_B3S3_f01_f12 = sum(NAT_129N_B3S3_f01_f12),
        TUMOR_129C_B3S3_f01_f12 = sum(TUMOR_129C_B3S3_f01_f12),
        TUMOR_130N_B3S3_f01_f12 = sum(TUMOR_130N_B3S3_f01_f12),
        NAT_130C_B3S3_f01_f12 = sum(NAT_130C_B3S3_f01_f12),
        REF_131_B3S3_f01_f12 = sum(REF_131_B3S3_f01_f12)
    ) %>%
  mutate(TUMOR_126_B3S3_f01_f12 = TUMOR_126_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_127N_B3S3_f01_f12 = NAT_127N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_127C_B3S3_f01_f12 = TUMOR_127C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_128N_B3S3_f01_f12 = NAT_128N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_128C_B3S3_f01_f12 = TUMOR_128C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_129N_B3S3_f01_f12 = NAT_129N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_129C_B3S3_f01_f12 = TUMOR_129C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_130N_B3S3_f01_f12 = TUMOR_130N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_130C_B3S3_f01_f12 = NAT_130C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  select(-REF_131_B3S3_f01_f12) 

dat_B3S4_f01_f12 <- bind_rows(B3S4_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B3S4_f01_f12 = sum(NAT_126_B3S4_f01_f12),
        NAT_127N_B3S4_f01_f12 = sum(NAT_127N_B3S4_f01_f12),
        TUMOR_127C_B3S4_f01_f12 = sum(TUMOR_127C_B3S4_f01_f12),
        NAT_128N_B3S4_f01_f12 = sum(NAT_128N_B3S4_f01_f12),
        NAT_128C_B3S4_f01_f12 = sum(NAT_128C_B3S4_f01_f12),
        TUMOR_129N_B3S4_f01_f12 = sum(TUMOR_129N_B3S4_f01_f12),
        TUMOR_129C_B3S4_f01_f12 = sum(TUMOR_129C_B3S4_f01_f12),
        NAT_130N_B3S4_f01_f12 = sum(NAT_130N_B3S4_f01_f12),
        NAT_130C_B3S4_f01_f12 = sum(NAT_130C_B3S4_f01_f12),
        REF_131_B3S4_f01_f12 = sum(REF_131_B3S4_f01_f12)
    ) %>%
  mutate(NAT_126_B3S4_f01_f12 = NAT_126_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_127N_B3S4_f01_f12 = NAT_127N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(TUMOR_127C_B3S4_f01_f12 = TUMOR_127C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_128N_B3S4_f01_f12 = NAT_128N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_128C_B3S4_f01_f12 = NAT_128C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(TUMOR_129N_B3S4_f01_f12 = TUMOR_129N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(TUMOR_129C_B3S4_f01_f12 = TUMOR_129C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_130N_B3S4_f01_f12 = NAT_130N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_130C_B3S4_f01_f12 = NAT_130C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  select(-REF_131_B3S4_f01_f12) 

dat_B4S1_f01_f12 <- bind_rows(B4S1_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B4S1_f01_f12 = sum(NAT_126_B4S1_f01_f12),
        TUMOR_127N_B4S1_f01_f12 = sum(TUMOR_127N_B4S1_f01_f12),
        NAT_127C_B4S1_f01_f12 = sum(NAT_127C_B4S1_f01_f12),
        NAT_128N_B4S1_f01_f12 = sum(NAT_128N_B4S1_f01_f12),
        NAT_128C_B4S1_f01_f12 = sum(NAT_128C_B4S1_f01_f12),
        NAT_129N_B4S1_f01_f12 = sum(NAT_129N_B4S1_f01_f12),
        TUMOR_129C_B4S1_f01_f12 = sum(TUMOR_129C_B4S1_f01_f12),
        NAT_130N_B4S1_f01_f12 = sum(NAT_130N_B4S1_f01_f12),
        NAT_130C_B4S1_f01_f12 = sum(NAT_130C_B4S1_f01_f12),
        REF_131_B4S1_f01_f12 = sum(REF_131_B4S1_f01_f12)
    ) %>%
  mutate(NAT_126_B4S1_f01_f12 = NAT_126_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(TUMOR_127N_B4S1_f01_f12 = TUMOR_127N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_127C_B4S1_f01_f12 = NAT_127C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_128N_B4S1_f01_f12 = NAT_128N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_128C_B4S1_f01_f12 = NAT_128C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_129N_B4S1_f01_f12 = NAT_129N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(TUMOR_129C_B4S1_f01_f12 = TUMOR_129C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_130N_B4S1_f01_f12 = NAT_130N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_130C_B4S1_f01_f12 = NAT_130C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  select(-REF_131_B4S1_f01_f12) 

dat_B4S2_f01_f12 <- bind_rows(B4S2_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B4S2_f01_f12 = sum(NAT_126_B4S2_f01_f12),
        TUMOR_127N_B4S2_f01_f12 = sum(TUMOR_127N_B4S2_f01_f12),
        TUMOR_127C_B4S2_f01_f12 = sum(TUMOR_127C_B4S2_f01_f12),
        NAT_128N_B4S2_f01_f12 = sum(NAT_128N_B4S2_f01_f12),
        NAT_128C_B4S2_f01_f12 = sum(NAT_128C_B4S2_f01_f12),
        TUMOR_129N_B4S2_f01_f12 = sum(TUMOR_129N_B4S2_f01_f12),
        NAT_129C_B4S2_f01_f12 = sum(NAT_129C_B4S2_f01_f12),
        TUMOR_130N_B4S2_f01_f12 = sum(TUMOR_130N_B4S2_f01_f12),
        NAT_130C_B4S2_f01_f12 = sum(NAT_130C_B4S2_f01_f12),
        REF_131_B4S2_f01_f12 = sum(REF_131_B4S2_f01_f12)
    ) %>%
  mutate(NAT_126_B4S2_f01_f12 = NAT_126_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_127N_B4S2_f01_f12 = TUMOR_127N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_127C_B4S2_f01_f12 = TUMOR_127C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_128N_B4S2_f01_f12 = NAT_128N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_128C_B4S2_f01_f12 = NAT_128C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_129N_B4S2_f01_f12 = TUMOR_129N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_129C_B4S2_f01_f12 = NAT_129C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_130N_B4S2_f01_f12 = TUMOR_130N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_130C_B4S2_f01_f12 = NAT_130C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  select(-REF_131_B4S2_f01_f12) 

dat_B4S3_f01_f12 <- bind_rows(B4S3_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B4S3_f01_f12 = sum(TUMOR_126_B4S3_f01_f12),
        NAT_127N_B4S3_f01_f12 = sum(NAT_127N_B4S3_f01_f12),
        TUMOR_127C_B4S3_f01_f12 = sum(TUMOR_127C_B4S3_f01_f12),
        NAT_128N_B4S3_f01_f12 = sum(NAT_128N_B4S3_f01_f12),
        NAT_128C_B4S3_f01_f12 = sum(NAT_128C_B4S3_f01_f12),
        TUMOR_129N_B4S3_f01_f12 = sum(TUMOR_129N_B4S3_f01_f12),
        NAT_129C_B4S3_f01_f12 = sum(NAT_129C_B4S3_f01_f12),
        TUMOR_130N_B4S3_f01_f12 = sum(TUMOR_130N_B4S3_f01_f12),
        NAT_130C_B4S3_f01_f12 = sum(NAT_130C_B4S3_f01_f12),
        REF_131_B4S3_f01_f12 = sum(REF_131_B4S3_f01_f12)
    ) %>%
  mutate(TUMOR_126_B4S3_f01_f12 = TUMOR_126_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_127N_B4S3_f01_f12 = NAT_127N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(TUMOR_127C_B4S3_f01_f12 = TUMOR_127C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_128N_B4S3_f01_f12 = NAT_128N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_128C_B4S3_f01_f12 = NAT_128C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(TUMOR_129N_B4S3_f01_f12 = TUMOR_129N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_129C_B4S3_f01_f12 = NAT_129C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(TUMOR_130N_B4S3_f01_f12 = TUMOR_130N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_130C_B4S3_f01_f12 = NAT_130C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  select(-REF_131_B4S3_f01_f12) 

dat_B4S4_f01_f12 <- bind_rows(B4S4_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
TUMOR_126_B4S4_f01_f12 = sum(TUMOR_126_B4S4_f01_f12),
TUMOR_127N_B4S4_f01_f12 = sum(TUMOR_127N_B4S4_f01_f12),
NAT_127C_B4S4_f01_f12 = sum(NAT_127C_B4S4_f01_f12),
NAT_128N_B4S4_f01_f12 = sum(NAT_128N_B4S4_f01_f12),
TUMOR_128C_B4S4_f01_f12 = sum(TUMOR_128C_B4S4_f01_f12),
NAT_129N_B4S4_f01_f12 = sum(NAT_129N_B4S4_f01_f12),
TUMOR_129C_B4S4_f01_f12 = sum(TUMOR_129C_B4S4_f01_f12),
NAT_130N_B4S4_f01_f12 = sum(NAT_130N_B4S4_f01_f12),
NAT_130C_B4S4_f01_f12 = sum(NAT_130C_B4S4_f01_f12),
REF_131_B4S4_f01_f12 = sum(REF_131_B4S4_f01_f12)
    ) %>%
  mutate(TUMOR_126_B4S4_f01_f12 = TUMOR_126_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(TUMOR_127N_B4S4_f01_f12 = TUMOR_127N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_127C_B4S4_f01_f12 = NAT_127C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_128N_B4S4_f01_f12 = NAT_128N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(TUMOR_128C_B4S4_f01_f12 = TUMOR_128C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_129N_B4S4_f01_f12 = NAT_129N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(TUMOR_129C_B4S4_f01_f12 = TUMOR_129C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_130N_B4S4_f01_f12 = NAT_130N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_130C_B4S4_f01_f12 = NAT_130C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  select(-REF_131_B4S4_f01_f12) 

dat_B5S1_f01_f12 <- bind_rows(B5S1_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B5S1_f01_f12 = sum(NAT_126_B5S1_f01_f12),
        TUMOR_127N_B5S1_f01_f12 = sum(TUMOR_127N_B5S1_f01_f12),
        NAT_127C_B5S1_f01_f12 = sum(NAT_127C_B5S1_f01_f12),
        NAT_128N_B5S1_f01_f12 = sum(NAT_128N_B5S1_f01_f12),
        NAT_128C_B5S1_f01_f12 = sum(NAT_128C_B5S1_f01_f12),
        NAT_129N_B5S1_f01_f12 = sum(NAT_129N_B5S1_f01_f12),
        TUMOR_129C_B5S1_f01_f12 = sum(TUMOR_129C_B5S1_f01_f12),
        TUMOR_130N_B5S1_f01_f12 = sum(TUMOR_130N_B5S1_f01_f12),
        TUMOR_130C_B5S1_f01_f12 = sum(TUMOR_130C_B5S1_f01_f12),
        REF_131_B5S1_f01_f12 = sum(REF_131_B5S1_f01_f12)
    ) %>%
  mutate(NAT_126_B5S1_f01_f12 = NAT_126_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_127N_B5S1_f01_f12 = TUMOR_127N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_127C_B5S1_f01_f12 = NAT_127C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_128N_B5S1_f01_f12 = NAT_128N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_128C_B5S1_f01_f12 = NAT_128C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_129N_B5S1_f01_f12 = NAT_129N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_129C_B5S1_f01_f12 = TUMOR_129C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_130N_B5S1_f01_f12 = TUMOR_130N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_130C_B5S1_f01_f12 = TUMOR_130C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  select(-REF_131_B5S1_f01_f12) 

dat_B5S2_f01_f12 <- bind_rows(B5S2_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        NAT_126_B5S2_f01_f12 = sum(NAT_126_B5S2_f01_f12),
        NAT_127N_B5S2_f01_f12 = sum(NAT_127N_B5S2_f01_f12),
        NAT_127C_B5S2_f01_f12 = sum(NAT_127C_B5S2_f01_f12),
        TUMOR_128N_B5S2_f01_f12 = sum(TUMOR_128N_B5S2_f01_f12),
        TUMOR_128C_B5S2_f01_f12 = sum(TUMOR_128C_B5S2_f01_f12),
        NAT_129N_B5S2_f01_f12 = sum(NAT_129N_B5S2_f01_f12),
        TUMOR_129C_B5S2_f01_f12 = sum(TUMOR_129C_B5S2_f01_f12),
        TUMOR_130N_B5S2_f01_f12 = sum(TUMOR_130N_B5S2_f01_f12),
        NAT_130C_B5S2_f01_f12 = sum(NAT_130C_B5S2_f01_f12),
        REF_131_B5S2_f01_f12 = sum(REF_131_B5S2_f01_f12)
    ) %>%
  mutate(NAT_126_B5S2_f01_f12 = NAT_126_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_127N_B5S2_f01_f12 = NAT_127N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_127C_B5S2_f01_f12 = NAT_127C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_128N_B5S2_f01_f12 = TUMOR_128N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_128C_B5S2_f01_f12 = TUMOR_128C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_129N_B5S2_f01_f12 = NAT_129N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_129C_B5S2_f01_f12 = TUMOR_129C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_130N_B5S2_f01_f12 = TUMOR_130N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_130C_B5S2_f01_f12 = NAT_130C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  select(-REF_131_B5S2_f01_f12) 

dat_B5S3_f01_f12 <- bind_rows(B5S3_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B5S3_f01_f12 = sum(TUMOR_126_B5S3_f01_f12),
        TUMOR_127N_B5S3_f01_f12 = sum(TUMOR_127N_B5S3_f01_f12),
        NAT_127C_B5S3_f01_f12 = sum(NAT_127C_B5S3_f01_f12),
        NAT_128N_B5S3_f01_f12 = sum(NAT_128N_B5S3_f01_f12),
        TUMOR_128C_B5S3_f01_f12 = sum(TUMOR_128C_B5S3_f01_f12),
        TUMOR_129N_B5S3_f01_f12 = sum(TUMOR_129N_B5S3_f01_f12),
        TUMOR_129C_B5S3_f01_f12 = sum(TUMOR_129C_B5S3_f01_f12),
        NAT_130N_B5S3_f01_f12 = sum(NAT_130N_B5S3_f01_f12),
        TUMOR_130C_B5S3_f01_f12 = sum(TUMOR_130C_B5S3_f01_f12),
        REF_131_B5S3_f01_f12 = sum(REF_131_B5S3_f01_f12)
    ) %>%
  mutate(TUMOR_126_B5S3_f01_f12 = TUMOR_126_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_127N_B5S3_f01_f12 = TUMOR_127N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(NAT_127C_B5S3_f01_f12 = NAT_127C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(NAT_128N_B5S3_f01_f12 = NAT_128N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_128C_B5S3_f01_f12 = TUMOR_128C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_129N_B5S3_f01_f12 = TUMOR_129N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_129C_B5S3_f01_f12 = TUMOR_129C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(NAT_130N_B5S3_f01_f12 = NAT_130N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_130C_B5S3_f01_f12 = TUMOR_130C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  select(-REF_131_B5S3_f01_f12) 

dat_B5S4_f01_f12 <- bind_rows(B5S4_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
        TUMOR_126_B5S4_f01_f12 = sum(TUMOR_126_B5S4_f01_f12),
        TUMOR_127N_B5S4_f01_f12 = sum(TUMOR_127N_B5S4_f01_f12),
        TUMOR_127C_B5S4_f01_f12 = sum(TUMOR_127C_B5S4_f01_f12),
        NAT_128N_B5S4_f01_f12 = sum(NAT_128N_B5S4_f01_f12),
        TUMOR_128C_B5S4_f01_f12 = sum(TUMOR_128C_B5S4_f01_f12),
        NAT_129N_B5S4_f01_f12 = sum(NAT_129N_B5S4_f01_f12),
        NAT_129C_B5S4_f01_f12 = sum(NAT_129C_B5S4_f01_f12),
        NAT_130N_B5S4_f01_f12 = sum(NAT_130N_B5S4_f01_f12),
        NAT_130C_B5S4_f01_f12 = sum(NAT_130C_B5S4_f01_f12),
        REF_131_B5S4_f01_f12 = sum(REF_131_B5S4_f01_f12)
    ) %>%
  mutate(TUMOR_126_B5S4_f01_f12 = TUMOR_126_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(TUMOR_127N_B5S4_f01_f12 = TUMOR_127N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(TUMOR_127C_B5S4_f01_f12 = TUMOR_127C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_128N_B5S4_f01_f12 = NAT_128N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(TUMOR_128C_B5S4_f01_f12 = TUMOR_128C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_129N_B5S4_f01_f12 = NAT_129N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_129C_B5S4_f01_f12 = NAT_129C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_130N_B5S4_f01_f12 = NAT_130N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_130C_B5S4_f01_f12 = NAT_130C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  select(-REF_131_B5S4_f01_f12) 

dat_B5S5_f01_f12 <- bind_rows(B5S5_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
TUMOR_126_B5S5_f01_f12 = sum(TUMOR_126_B5S5_f01_f12),
TUMOR_127N_B5S5_f01_f12 = sum(TUMOR_127N_B5S5_f01_f12),
TUMOR_127C_B5S5_f01_f12 = sum(TUMOR_127C_B5S5_f01_f12),
NAT_128N_B5S5_f01_f12 = sum(NAT_128N_B5S5_f01_f12),
NAT_128C_B5S5_f01_f12 = sum(NAT_128C_B5S5_f01_f12),
TUMOR_129N_B5S5_f01_f12 = sum(TUMOR_129N_B5S5_f01_f12),
TUMOR_129C_B5S5_f01_f12 = sum(TUMOR_129C_B5S5_f01_f12),
NAT_130N_B5S5_f01_f12 = sum(NAT_130N_B5S5_f01_f12),
NAT_130C_B5S5_f01_f12 = sum(NAT_130C_B5S5_f01_f12),
REF_131_B5S5_f01_f12 = sum(REF_131_B5S5_f01_f12)
    ) %>%
  mutate(TUMOR_126_B5S5_f01_f12 = TUMOR_126_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_127N_B5S5_f01_f12 = TUMOR_127N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_127C_B5S5_f01_f12 = TUMOR_127C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_128N_B5S5_f01_f12 = NAT_128N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_128C_B5S5_f01_f12 = NAT_128C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_129N_B5S5_f01_f12 = TUMOR_129N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_129C_B5S5_f01_f12 = TUMOR_129C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_130N_B5S5_f01_f12 = NAT_130N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_130C_B5S5_f01_f12 = NAT_130C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  select(-REF_131_B5S5_f01_f12) 

dat_B5S6_f01_f12 <- bind_rows(B5S6_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, sequence) %>% 
    summarise(
TUMOR_126_B5S6_f01_f12 = sum(TUMOR_126_B5S6_f01_f12),
NAT_127N_B5S6_f01_f12 = sum(NAT_127N_B5S6_f01_f12),
NAT_127C_B5S6_f01_f12 = sum(NAT_127C_B5S6_f01_f12),
NAT_128N_B5S6_f01_f12 = sum(NAT_128N_B5S6_f01_f12),
NAT_128C_B5S6_f01_f12 = sum(NAT_128C_B5S6_f01_f12),
NAT_129N_B5S6_f01_f12 = sum(NAT_129N_B5S6_f01_f12),
NAT_129C_B5S6_f01_f12 = sum(NAT_129C_B5S6_f01_f12),
NAT_130N_B5S6_f01_f12 = sum(NAT_130N_B5S6_f01_f12),
REF_130C_B5S6_f01_f12 = sum(REF_130C_B5S6_f01_f12),
REF_131_B5S6_f01_f12 = sum(REF_131_B5S6_f01_f12)
    ) %>%
  mutate(TUMOR_126_B5S6_f01_f12 = TUMOR_126_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_127N_B5S6_f01_f12 = NAT_127N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_127C_B5S6_f01_f12 = NAT_127C_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_128N_B5S6_f01_f12 = NAT_128N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_128C_B5S6_f01_f12 = NAT_128C_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_129N_B5S6_f01_f12 = NAT_129N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_129C_B5S6_f01_f12 = NAT_129C_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_130N_B5S6_f01_f12 = NAT_130N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  select(-REF_130C_B5S6_f01_f12 ) %>%
  select(-REF_131_B5S6_f01_f12 ) 

all_batches <- list(
    dat_B1S1_f01_f12,
    dat_B1S2_f01_f12, 
    dat_B1S3_f01_f12, 
    dat_B1S4_f01_f12, 
    dat_B2S1_f01_f12, 
    dat_B2S2_f01_f12, 
    dat_B2S3_f01_f12, 
    dat_B2S4_f01_f12,
    dat_B3S1_f01_f12, 
    dat_B3S2_f01_f12, 
    dat_B3S3_f01_f12,
    dat_B3S4_f01_f12, 
    dat_B4S1_f01_f12,
    dat_B4S2_f01_f12, 
    dat_B4S3_f01_f12, 
    dat_B4S4_f01_f12, 
    dat_B5S1_f01_f12, 
    dat_B5S2_f01_f12, 
    dat_B5S3_f01_f12, 
    dat_B5S4_f01_f12, 
    dat_B5S5_f01_f12, 
    dat_B5S6_f01_f12
)
#Loop selecting wanted protein from each list
sel_pro <- list()
for (i in 1:length(all_batches)) {
  sel_pro[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P25815") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
view(sel_pro[[4]])
dim(sel_pro[[4]])

#Combine selected protein of all batches
Selected_protein <- sel_pro %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()
dim(Selected_protein)
view(Selected_protein)

#Load PTMs + peptide seq
ALL_seq <- fread("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/PSM_TMT_PTM.csv")
view(head(ALL_seq))

class(Selected_protein)
class(ALL_seq)

#Merge, based on modified peptide sequence
protein_PTM <- merge(
  Selected_protein, ALL_seq, by.x=sequence, by.y=sequence.x, all.x=TRUE
)
dim(protein_PTM)
view(protein_PTM)
  #Order columns
colorder <- c("Protein.Group.Accessions",
  "sequence", 
  "mod",
  "mod_mass",
  "TUMOR_127C_B1S1_f01_f12", 
  "TUMOR_128N_B1S1_f01_f12", 
  "TUMOR_128C_B1S1_f01_f12", 
  "TUMOR_129N_B1S1_f01_f12", 
  "TUMOR_129C_B1S1_f01_f12", 
  "TUMOR_130N_B1S1_f01_f12", 
  "TUMOR_130C_B1S1_f01_f12", 
  "TUMOR_127N_B1S2_f01_f12", 
  "TUMOR_128N_B1S2_f01_f12", 
  "TUMOR_129N_B1S2_f01_f12", 
  "TUMOR_130N_B1S2_f01_f12", 
  "TUMOR_130C_B1S2_f01_f12", 
  "TUMOR_126_B1S3_f01_f12", 
  "TUMOR_127C_B1S3_f01_f12", 
  "TUMOR_130N_B1S3_f01_f12", 
  "TUMOR_130C_B1S3_f01_f12", 
  "TUMOR_127N_B1S4_f01_f12", 
  "TUMOR_128N_B1S4_f01_f12", 
  "TUMOR_128C_B1S4_f01_f12", 
  "TUMOR_129N_B1S4_f01_f12", 
  "TUMOR_129C_B1S4_f01_f12", 
  "TUMOR_130C_B1S4_f01_f12", 
  "TUMOR_126_B2S1_f01_f12", 
  "TUMOR_127N_B2S1_f01_f12", 
  "TUMOR_127C_B2S1_f01_f12", 
  "TUMOR_129N_B2S1_f01_f12", 
  "TUMOR_130N_B2S1_f01_f12", 
  "TUMOR_130C_B2S1_f01_f12", 
  "TUMOR_126_B2S2_f01_f12", 
  "TUMOR_127N_B2S2_f01_f12", 
  "TUMOR_127C_B2S2_f01_f12", 
  "TUMOR_129N_B2S2_f01_f12", 
  "TUMOR_130N_B2S2_f01_f12", 
  "TUMOR_130C_B2S2_f01_f12", 
  "TUMOR_126_B2S3_f01_f12", 
  "TUMOR_127C_B2S3_f01_f12", 
  "TUMOR_130N_B2S3_f01_f12", 
  "TUMOR_126_B2S4_f01_f12", 
  "TUMOR_127N_B2S4_f01_f12", 
  "TUMOR_128C_B2S4_f01_f12", 
  "TUMOR_129N_B2S4_f01_f12", 
  "TUMOR_130N_B2S4_f01_f12", 
  "TUMOR_126_B3S1_f01_f12", 
  "TUMOR_128N_B3S1_f01_f12", 
  "TUMOR_128C_B3S1_f01_f12", 
  "TUMOR_129N_B3S1_f01_f12", 
  "TUMOR_129C_B3S1_f01_f12", 
  "TUMOR_130N_B3S1_f01_f12", 
  "TUMOR_130C_B3S1_f01_f12", 
  "TUMOR_128C_B3S2_f01_f12", 
  "TUMOR_130N_B3S2_f01_f12", 
  "TUMOR_130C_B3S2_f01_f12", 
  "TUMOR_126_B3S3_f01_f12", 
  "TUMOR_127C_B3S3_f01_f12", 
  "TUMOR_128C_B3S3_f01_f12", 
  "TUMOR_129C_B3S3_f01_f12", 
  "TUMOR_130N_B3S3_f01_f12", 
  "TUMOR_127C_B3S4_f01_f12", 
  "TUMOR_129N_B3S4_f01_f12", 
  "TUMOR_129C_B3S4_f01_f12", 
  "TUMOR_127N_B4S1_f01_f12", 
  "TUMOR_129C_B4S1_f01_f12", 
  "TUMOR_127N_B4S2_f01_f12", 
  "TUMOR_127C_B4S2_f01_f12", 
  "TUMOR_129N_B4S2_f01_f12", 
  "TUMOR_130N_B4S2_f01_f12", 
  "TUMOR_126_B4S3_f01_f12", 
  "TUMOR_127C_B4S3_f01_f12", 
  "TUMOR_129N_B4S3_f01_f12", 
  "TUMOR_130N_B4S3_f01_f12", 
  "TUMOR_126_B4S4_f01_f12", 
  "TUMOR_127N_B4S4_f01_f12", 
  "TUMOR_128C_B4S4_f01_f12", 
  "TUMOR_129C_B4S4_f01_f12", 
  "TUMOR_127N_B5S1_f01_f12", 
  "TUMOR_129C_B5S1_f01_f12", 
  "TUMOR_130N_B5S1_f01_f12", 
  "TUMOR_130C_B5S1_f01_f12", 
  "TUMOR_128N_B5S2_f01_f12", 
  "TUMOR_128C_B5S2_f01_f12", 
  "TUMOR_129C_B5S2_f01_f12", 
  "TUMOR_130N_B5S2_f01_f12", 
  "TUMOR_126_B5S3_f01_f12", 
  "TUMOR_127N_B5S3_f01_f12", 
  "TUMOR_128C_B5S3_f01_f12", 
  "TUMOR_129N_B5S3_f01_f12", 
  "TUMOR_129C_B5S3_f01_f12", 
  "TUMOR_130C_B5S3_f01_f12", 
  "TUMOR_126_B5S4_f01_f12", 
  "TUMOR_127N_B5S4_f01_f12", 
  "TUMOR_127C_B5S4_f01_f12", 
  "TUMOR_128C_B5S4_f01_f12", 
  "TUMOR_126_B5S5_f01_f12", 
  "TUMOR_127N_B5S5_f01_f12", 
  "TUMOR_127C_B5S5_f01_f12", 
  "TUMOR_129N_B5S5_f01_f12", 
  "TUMOR_129C_B5S5_f01_f12", 
  "TUMOR_126_B5S6_f01_f12",
  "NAT_126_B1S1_f01_f12", 
  "NAT_127N_B1S1_f01_f12", 
  "NAT_126_B1S2_f01_f12", 
  "NAT_127C_B1S2_f01_f12", 
  "NAT_128C_B1S2_f01_f12", 
  "NAT_129C_B1S2_f01_f12", 
  "NAT_127N_B1S3_f01_f12", 
  "NAT_128N_B1S3_f01_f12", 
  "NAT_128C_B1S3_f01_f12", 
  "NAT_129N_B1S3_f01_f12", 
  "NAT_129C_B1S3_f01_f12", 
  "NAT_126_B1S4_f01_f12", 
  "NAT_127C_B1S4_f01_f12", 
  "NAT_130N_B1S4_f01_f12", 
  "NAT_128N_B2S1_f01_f12", 
  "NAT_128C_B2S1_f01_f12", 
  "NAT_129C_B2S1_f01_f12", 
  "NAT_128N_B2S2_f01_f12", 
  "NAT_128C_B2S2_f01_f12", 
  "NAT_129C_B2S2_f01_f12", 
  "NAT_127N_B2S3_f01_f12", 
  "NAT_128N_B2S3_f01_f12", 
  "NAT_128C_B2S3_f01_f12", 
  "NAT_129N_B2S3_f01_f12", 
  "NAT_129C_B2S3_f01_f12", 
  "NAT_130C_B2S3_f01_f12", 
  "NAT_127C_B2S4_f01_f12", 
  "NAT_128N_B2S4_f01_f12", 
  "NAT_129C_B2S4_f01_f12", 
  "NAT_130C_B2S4_f01_f12", 
  "NAT_127N_B3S1_f01_f12", 
  "NAT_127C_B3S1_f01_f12", 
  "NAT_126_B3S2_f01_f12", 
  "NAT_127N_B3S2_f01_f12", 
  "NAT_127C_B3S2_f01_f12", 
  "NAT_128N_B3S2_f01_f12", 
  "NAT_129N_B3S2_f01_f12", 
  "NAT_129C_B3S2_f01_f12", 
  "NAT_127N_B3S3_f01_f12", 
  "NAT_128N_B3S3_f01_f12", 
  "NAT_129N_B3S3_f01_f12", 
  "NAT_130C_B3S3_f01_f12", 
  "NAT_126_B3S4_f01_f12", 
  "NAT_127N_B3S4_f01_f12", 
  "NAT_128N_B3S4_f01_f12", 
  "NAT_128C_B3S4_f01_f12", 
  "NAT_130N_B3S4_f01_f12", 
  "NAT_130C_B3S4_f01_f12", 
  "NAT_126_B4S1_f01_f12", 
  "NAT_127C_B4S1_f01_f12", 
  "NAT_128N_B4S1_f01_f12", 
  "NAT_128C_B4S1_f01_f12", 
  "NAT_129N_B4S1_f01_f12", 
  "NAT_130N_B4S1_f01_f12", 
  "NAT_130C_B4S1_f01_f12", 
  "NAT_126_B4S2_f01_f12", 
  "NAT_128N_B4S2_f01_f12", 
  "NAT_128C_B4S2_f01_f12", 
  "NAT_129C_B4S2_f01_f12", 
  "NAT_130C_B4S2_f01_f12", 
  "NAT_127N_B4S3_f01_f12", 
  "NAT_128N_B4S3_f01_f12", 
  "NAT_128C_B4S3_f01_f12", 
  "NAT_129C_B4S3_f01_f12", 
  "NAT_130C_B4S3_f01_f12", 
  "NAT_127C_B4S4_f01_f12", 
  "NAT_128N_B4S4_f01_f12", 
  "NAT_129N_B4S4_f01_f12", 
  "NAT_130N_B4S4_f01_f12", 
  "NAT_130C_B4S4_f01_f12", 
  "NAT_126_B5S1_f01_f12", 
  "NAT_127C_B5S1_f01_f12", 
  "NAT_128N_B5S1_f01_f12", 
  "NAT_128C_B5S1_f01_f12", 
  "NAT_129N_B5S1_f01_f12", 
  "NAT_126_B5S2_f01_f12", 
  "NAT_127N_B5S2_f01_f12", 
  "NAT_127C_B5S2_f01_f12", 
  "NAT_129N_B5S2_f01_f12", 
  "NAT_130C_B5S2_f01_f12", 
  "NAT_127C_B5S3_f01_f12", 
  "NAT_128N_B5S3_f01_f12", 
  "NAT_130N_B5S3_f01_f12", 
  "NAT_128N_B5S4_f01_f12", 
  "NAT_129N_B5S4_f01_f12", 
  "NAT_129C_B5S4_f01_f12", 
  "NAT_130N_B5S4_f01_f12", 
  "NAT_130C_B5S4_f01_f12", 
  "NAT_128N_B5S5_f01_f12", 
  "NAT_128C_B5S5_f01_f12", 
  "NAT_130N_B5S5_f01_f12", 
  "NAT_130C_B5S5_f01_f12", 
  "NAT_127N_B5S6_f01_f12", 
  "NAT_127C_B5S6_f01_f12", 
  "NAT_128N_B5S6_f01_f12", 
  "NAT_128C_B5S6_f01_f12", 
  "NAT_129N_B5S6_f01_f12", 
  "NAT_129C_B5S6_f01_f12", 
  "NAT_130N_B5S6_f01_f12"
)
protein_PTM <- protein_PTM[, colorder] %>%
    as_data_frame() %>%
    drop_na()

#Differential expression in analyis
  #Define sig difference in TMT intensities for all PTMs of a protein
ttestFunc <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
rawpvalue = apply(dat_col_ordered, 1, ttestFunc, grp1 = c(2:99), grp2 = c(100:198))
view(rawpvalue)



#Filter on accession number
#then only keep peptide sequences
#then reduce by sequence
#keeps 1 row per sequence