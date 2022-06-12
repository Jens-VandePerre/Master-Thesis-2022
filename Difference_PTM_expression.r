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


view(dat_B1S1_f01_f12)
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
saveRDS(all_batches, "/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/Differential_PTM_analysis")

CRC_proteins <- c("P25815",
"Q12884",
"P36952",
"P06702",
"P83881",
"Q96CG8",
"P06731",
"P05109",
"Q8NFJ5",
"P50454",
"P80511",
"Q9NR99",
"Q01650",
"P31949",
"P40261",
"P02788",
"O00469",
"P35442",
"O76021",
"P08246",
"P05164",
"Q9NR30",
"P30273",
"Q99715",
"P49913",
"P52926",
"Q8NCL4",
"O00425",
"P80188",
"P62760",
"P11388")


#Load PTMs + peptide seq
ALL_seq <- fread("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/PSM_TMT_PTM.csv") %>%
  rename(sequence = "sequence.x") %>%
  select(-c("126":"131")) %>%
  select(-index_filename, -index)
view(head(ALL_seq, 200))
dim(ALL_seq)


#Loop selecting wanted protein from each list
P25815 <- list()
for (i in 1:length(all_batches)) {
  P25815[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P25815") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P25815 <- P25815 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q12884 <- list()
for (i in 1:length(all_batches)) {
  Q12884[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q12884") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q12884 <- Q12884 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P36952 <- list()
for (i in 1:length(all_batches)) {
  P36952[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P36952") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P36952 <- P36952 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P06702 <- list()
for (i in 1:length(all_batches)) {
  P06702[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P06702") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P06702 <- P06702 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P83881 <- list()
for (i in 1:length(all_batches)) {
  P83881[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P83881") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P83881 <- P83881 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q96CG8 <- list() #NOT FOUND
for (i in 1:length(all_batches)) {
  Q96CG8[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q96CG8") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q96CG8 <- Q96CG8 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P06731 <- list() #NOT FOUND
for (i in 1:length(all_batches)) {
  P06731[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P06731") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P06731 <- P06731 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P05109 <- list()
for (i in 1:length(all_batches)) {
  P05109[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P05109") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P05109 <- P05109 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q8NFJ5 <- list() #NOT FOUND
for (i in 1:length(all_batches)) {
  Q8NFJ5[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q8NFJ5") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q8NFJ5 <- Q8NFJ5 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P50454 <- list()
for (i in 1:length(all_batches)) {
  P50454[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P50454") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P50454 <- P50454 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P80511 <- list()
for (i in 1:length(all_batches)) {
  P80511[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P80511") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P80511 <- P80511 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q9NR99 <- list()
for (i in 1:length(all_batches)) {
  Q9NR99[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q9NR99") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q9NR99 <- Q9NR99 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q01650 <- list()
for (i in 1:length(all_batches)) {
  Q01650[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q01650") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q01650 <- Q01650 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P31949 <- list()
for (i in 1:length(all_batches)) {
  P31949[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P31949") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P31949 <- P31949 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P40261 <- list()
for (i in 1:length(all_batches)) {
  P40261[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P40261") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P40261 <- P40261 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P02788 <- list()
for (i in 1:length(all_batches)) {
  P02788[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P02788") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P02788 <- P02788 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


O00469 <- list()
for (i in 1:length(all_batches)) {
  O00469[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "O00469") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
O00469 <- O00469 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P35442 <- list()
for (i in 1:length(all_batches)) {
  P35442[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P35442") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P35442 <- P35442 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


O76021 <- list()
for (i in 1:length(all_batches)) {
  O76021[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "O76021") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
O76021 <- O76021 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P08246 <- list()
for (i in 1:length(all_batches)) {
  P08246[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P08246") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P08246 <- P08246 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P05164 <- list()
for (i in 1:length(all_batches)) {
  P05164[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P05164") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P05164 <- P05164 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q9NR30 <- list()
for (i in 1:length(all_batches)) {
  Q9NR30[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q9NR30") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q9NR30 <- Q9NR30 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P30273 <- list()
for (i in 1:length(all_batches)) {
  P30273[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P30273") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P30273 <- P30273 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q99715 <- list()
for (i in 1:length(all_batches)) {
  Q99715[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q99715") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q99715 <- Q99715 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P49913 <- list() #NOT FOUND
for (i in 1:length(all_batches)) {
  P49913[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P49913") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P49913 <- P49913 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()

P52926 <- list()
for (i in 1:length(all_batches)) {
  P52926[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P52926") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P52926 <- P52926 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


Q8NCL4 <- list()
for (i in 1:length(all_batches)) {
  Q8NCL4[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q8NCL4") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
Q8NCL4 <- Q8NCL4 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


O00425 <- list()
for (i in 1:length(all_batches)) {
  O00425[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "O00425") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
O00425 <- O00425 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P80188 <- list()
for (i in 1:length(all_batches)) {
  P80188[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P80188") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P80188 <- P80188 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P62760 <- list()
for (i in 1:length(all_batches)) {
  P62760[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P62760") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P62760 <- P62760 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()


P11388 <- list()
for (i in 1:length(all_batches)) {
  P11388[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P11388") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P11388 <- P11388 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame()



#Merge, based on modified peptide sequence
#Differential expression in analyis
  #Define sig difference in TMT intensities for all PTMs of a protein
    #Bind t-test results to PTMs + BH correction
ttestFunc <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y, na.action=na.omit)
  results$p.value
}

P25815_PTM <- merge(P25815, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P25815_PTMrawpvalue = apply(P25815_PTM, 1, ttestFunc, grp1 = colnames(P25815_PTM)[grep("TUMOR",colnames(P25815_PTM))]
, grp2 = colnames(P25815_PTM)[grep("NAT",colnames(P25815_PTM))])
P25815_PTM_t_test <- cbind(P25815_PTM, P25815_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q12884_PTM <- merge(Q12884, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q12884_PTMrawpvalue = apply(Q12884_PTM, 1, ttestFunc, grp1 = colnames(Q12884_PTM)[grep("TUMOR",colnames(Q12884_PTM))]
, grp2 = colnames(Q12884_PTM)[grep("NAT",colnames(Q12884_PTM))])
Q12884_PTM_t_test <- cbind(Q12884_PTM, Q12884_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P36952_PTM <- merge(P36952, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P36952_PTMrawpvalue = apply(P36952_PTM, 1, ttestFunc, grp1 = colnames(P36952_PTM)[grep("TUMOR",colnames(P36952_PTM))]
, grp2 = colnames(P36952_PTM)[grep("NAT",colnames(P36952_PTM))])
P36952_PTM_t_test <- cbind(P36952_PTM, P36952_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P06702_PTM <- merge(P06702, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P06702_PTMrawpvalue = apply(P06702_PTM, 1, ttestFunc, grp1 = colnames(P06702_PTM)[grep("TUMOR",colnames(P06702_PTM))]
, grp2 = colnames(P06702_PTM)[grep("NAT",colnames(P06702_PTM))])
P06702_PTM_t_test <- cbind(P06702_PTM, P06702_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P83881_PTM <- merge(P83881, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P83881_PTMrawpvalue = apply(P83881_PTM, 1, ttestFunc, grp1 = colnames(P83881_PTM)[grep("TUMOR",colnames(P83881_PTM))]
, grp2 = colnames(P83881_PTM)[grep("NAT",colnames(P83881_PTM))])
P83881_PTM_t_test <- cbind(P83881_PTM, P83881_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q96CG8_PTM <- merge(Q96CG8, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q96CG8_PTMrawpvalue = apply(Q96CG8_PTM, 1, ttestFunc, grp1 = colnames(Q96CG8_PTM)[grep("TUMOR",colnames(Q96CG8_PTM))]
, grp2 = colnames(Q96CG8_PTM)[grep("NAT",colnames(Q96CG8_PTM))])
Q96CG8_PTM_t_test <- cbind(Q96CG8_PTM, Q96CG8_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P06731_PTM <- merge(P06731, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P06731_PTMrawpvalue = apply(P06731_PTM, 1, ttestFunc, grp1 = colnames(P06731_PTM)[grep("TUMOR",colnames(P06731_PTM))]
, grp2 = colnames(P06731_PTM)[grep("NAT",colnames(P06731_PTM))])
P06731_PTM_t_test <- cbind(P06731_PTM, P06731_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P05109_PTM <- merge(P05109, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P05109_PTMrawpvalue = apply(P05109_PTM, 1, ttestFunc, grp1 = colnames(P05109_PTM)[grep("TUMOR",colnames(P05109_PTM))]
, grp2 = colnames(P05109_PTM)[grep("NAT",colnames(P05109_PTM))])
P05109_PTM_t_test <- cbind(P05109_PTM, P05109_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q8NFJ5_PTM <- merge(Q8NFJ5, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q8NFJ5_PTMrawpvalue = apply(Q8NFJ5_PTM, 1, ttestFunc, grp1 = colnames(Q8NFJ5_PTM)[grep("TUMOR",colnames(Q8NFJ5_PTM))]
, grp2 = colnames(Q8NFJ5_PTM)[grep("NAT",colnames(Q8NFJ5_PTM))])
Q8NFJ5_PTM_t_test <- cbind(Q8NFJ5_PTM, Q8NFJ5_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P50454_PTM <- merge(P50454, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P50454_PTMrawpvalue = apply(P50454_PTM, 1, ttestFunc, grp1 = colnames(P50454_PTM)[grep("TUMOR",colnames(P50454_PTM))]
, grp2 = colnames(P50454_PTM)[grep("NAT",colnames(P50454_PTM))])
P50454_PTM_t_test <- cbind(P50454_PTM, P50454_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P80511_PTM <- merge(P80511, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P80511_PTMrawpvalue = apply(P80511_PTM, 1, ttestFunc, grp1 = colnames(P80511_PTM)[grep("TUMOR",colnames(P80511_PTM))]
, grp2 = colnames(P80511_PTM)[grep("NAT",colnames(P80511_PTM))])
P80511_PTM_t_test <- cbind(P80511_PTM, P80511_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q9NR99_PTM <- merge(Q9NR99, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q9NR99_PTMrawpvalue = apply(Q9NR99_PTM, 1, ttestFunc, grp1 = colnames(Q9NR99_PTM)[grep("TUMOR",colnames(Q9NR99_PTM))]
, grp2 = colnames(Q9NR99_PTM)[grep("NAT",colnames(Q9NR99_PTM))])
Q9NR99_PTM_t_test <- cbind(Q9NR99_PTM, Q9NR99_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q01650_PTM <- merge(Q01650, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q01650_PTMrawpvalue = apply(Q01650_PTM, 1, ttestFunc, grp1 = colnames(Q01650_PTM)[grep("TUMOR",colnames(Q01650_PTM))]
, grp2 = colnames(Q01650_PTM)[grep("NAT",colnames(Q01650_PTM))])
Q01650_PTM_t_test <- cbind(Q01650_PTM, Q01650_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P31949_PTM <- merge(P31949, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P31949_PTMrawpvalue = apply(P31949_PTM, 1, ttestFunc, grp1 = colnames(P31949_PTM)[grep("TUMOR",colnames(P31949_PTM))]
, grp2 = colnames(P31949_PTM)[grep("NAT",colnames(P31949_PTM))])
P31949_PTM_t_test <- cbind(P31949_PTM, P31949_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P40261_PTM <- merge(P40261, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P40261_PTMrawpvalue = apply(P40261_PTM, 1, ttestFunc, grp1 = colnames(P40261_PTM)[grep("TUMOR",colnames(P40261_PTM))]
, grp2 = colnames(P40261_PTM)[grep("NAT",colnames(P40261_PTM))])
P40261_PTM_t_test <- cbind(P40261_PTM, P40261_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P02788_PTM <- merge(P02788, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P02788_PTMrawpvalue = apply(P02788_PTM, 1, ttestFunc, grp1 = colnames(P02788_PTM)[grep("TUMOR",colnames(P02788_PTM))]
, grp2 = colnames(P02788_PTM)[grep("NAT",colnames(P02788_PTM))])
P02788_PTM_t_test <- cbind(P02788_PTM, P02788_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

O00469_PTM <- merge(O00469, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
O00469_PTMrawpvalue = apply(O00469_PTM, 1, ttestFunc, grp1 = colnames(O00469_PTM)[grep("TUMOR",colnames(O00469_PTM))]
, grp2 = colnames(O00469_PTM)[grep("NAT",colnames(O00469_PTM))])
O00469_PTM_t_test <- cbind(O00469_PTM, O00469_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P35442_PTM <- merge(P35442, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P35442_PTMrawpvalue = apply(P35442_PTM, 1, ttestFunc, grp1 = colnames(P35442_PTM)[grep("TUMOR",colnames(P35442_PTM))]
, grp2 = colnames(P35442_PTM)[grep("NAT",colnames(P35442_PTM))])
P35442_PTM_t_test <- cbind(P35442_PTM, P35442_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

O76021_PTM <- merge(O76021, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
O76021_PTMrawpvalue = apply(O76021_PTM, 1, ttestFunc, grp1 = colnames(O76021_PTM)[grep("TUMOR",colnames(O76021_PTM))]
, grp2 = colnames(O76021_PTM)[grep("NAT",colnames(O76021_PTM))])
O76021_PTM_t_test <- cbind(O76021_PTM, O76021_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P08246_PTM <- merge(P08246, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P08246_PTMrawpvalue = apply(P08246_PTM, 1, ttestFunc, grp1 = colnames(P08246_PTM)[grep("TUMOR",colnames(P08246_PTM))]
, grp2 = colnames(P08246_PTM)[grep("NAT",colnames(P08246_PTM))])
P08246_PTM_t_test <- cbind(P08246_PTM, P08246_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P05164_PTM <- merge(P05164, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P05164_PTMrawpvalue = apply(P05164_PTM, 1, ttestFunc, grp1 = colnames(P05164_PTM)[grep("TUMOR",colnames(P05164_PTM))]
, grp2 = colnames(P05164_PTM)[grep("NAT",colnames(P05164_PTM))])
P05164_PTM_t_test <- cbind(P05164_PTM, P05164_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q9NR30_PTM <- merge(Q9NR30, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q9NR30_PTMrawpvalue = apply(Q9NR30_PTM, 1, ttestFunc, grp1 = colnames(Q9NR30_PTM)[grep("TUMOR",colnames(Q9NR30_PTM))]
, grp2 = colnames(Q9NR30_PTM)[grep("NAT",colnames(Q9NR30_PTM))])
Q9NR30_PTM_t_test <- cbind(Q9NR30_PTM, Q9NR30_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P30273_PTM <- merge(P30273, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P30273_PTMrawpvalue = apply(P30273_PTM, 1, ttestFunc, grp1 = colnames(P30273_PTM)[grep("TUMOR",colnames(P30273_PTM))]
, grp2 = colnames(P30273_PTM)[grep("NAT",colnames(P30273_PTM))])
P30273_PTM_t_test <- cbind(P30273_PTM, P30273_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q99715_PTM <- merge(Q99715, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q99715_PTMrawpvalue = apply(Q99715_PTM, 1, ttestFunc, grp1 = colnames(Q99715_PTM)[grep("TUMOR",colnames(Q99715_PTM))]
, grp2 = colnames(Q99715_PTM)[grep("NAT",colnames(Q99715_PTM))])
Q99715_PTM_t_test <- cbind(Q99715_PTM, Q99715_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P49913_PTM <- merge(P49913, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P49913_PTMrawpvalue = apply(P49913_PTM, 1, ttestFunc, grp1 = colnames(P49913_PTM)[grep("TUMOR",colnames(P49913_PTM))]
, grp2 = colnames(P49913_PTM)[grep("NAT",colnames(P49913_PTM))])
P49913_PTM_t_test <- cbind(P49913_PTM, P49913_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P52926_PTM <- merge(P52926, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P52926_PTMrawpvalue = apply(P52926_PTM, 1, ttestFunc, grp1 = colnames(P52926_PTM)[grep("TUMOR",colnames(P52926_PTM))]
, grp2 = colnames(P52926_PTM)[grep("NAT",colnames(P52926_PTM))])
P52926_PTM_t_test <- cbind(P52926_PTM, P52926_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

Q8NCL4_PTM <- merge(Q8NCL4, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
Q8NCL4_PTMrawpvalue = apply(Q8NCL4_PTM, 1, ttestFunc, grp1 = colnames(Q8NCL4_PTM)[grep("TUMOR",colnames(Q8NCL4_PTM))]
, grp2 = colnames(Q8NCL4_PTM)[grep("NAT",colnames(Q8NCL4_PTM))])
Q8NCL4_PTM_t_test <- cbind(Q8NCL4_PTM, Q8NCL4_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

O00425_PTM <- merge(O00425, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
O00425_PTMrawpvalue = apply(O00425_PTM, 1, ttestFunc, grp1 = colnames(O00425_PTM)[grep("TUMOR",colnames(O00425_PTM))]
, grp2 = colnames(O00425_PTM)[grep("NAT",colnames(O00425_PTM))])
O00425_PTM_t_test <- cbind(O00425_PTM, O00425_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P80188_PTM <- merge(P80188, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P80188_PTMrawpvalue = apply(P80188_PTM, 1, ttestFunc, grp1 = colnames(P80188_PTM)[grep("TUMOR",colnames(P80188_PTM))]
, grp2 = colnames(P80188_PTM)[grep("NAT",colnames(P80188_PTM))])
P80188_PTM_t_test <- cbind(P80188_PTM, P80188_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P62760_PTM <- merge(P62760, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P62760_PTMrawpvalue = apply(P62760_PTM, 1, ttestFunc, grp1 = colnames(P62760_PTM)[grep("TUMOR",colnames(P62760_PTM))]
, grp2 = colnames(P62760_PTM)[grep("NAT",colnames(P62760_PTM))])
P62760_PTM_t_test <- cbind(P62760_PTM, P62760_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")

P11388_PTM <- merge(P11388, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count)
P11388_PTMrawpvalue = apply(P11388_PTM, 1, ttestFunc, grp1 = colnames(P11388_PTM)[grep("TUMOR",colnames(P11388_PTM))]
, grp2 = colnames(P11388_PTM)[grep("NAT",colnames(P11388_PTM))])
P11388_PTM_t_test <- cbind(P11388_PTM, P11388_PTMrawpvalue) %>%
  select(index_filename, sequence, sequence_no_mod, 
    mod, mod_mass, count, x) %>%
  rename(p.value = x) %>%
  BH = p.adjust(p.value, 
               method = "BH")











