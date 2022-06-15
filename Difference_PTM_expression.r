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
library("datawizard")


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
all_batches <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/Differential_PTM_analysis")


#Load PTMs + peptide seq
ALL_seq <- fread("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/PSM_TMT_PTM.csv") %>%
  rename(sequence = "sequence.x") %>%
  select(-c("126":"131")) %>%
  select(-index_filename, -index)
view(head(ALL_seq, 200))
dim(ALL_seq)


#Loop selecting wanted protein from each list
  #Create dataframe with all peptide seq for each selected protein
P05783 <- list()
for (i in 1:length(all_batches)) {
  P05783[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P05783") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P05783 <- P05783 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 
view(P05783_med)

  P61626 <- list()
for (i in 1:length(all_batches)) {
  P61626[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P61626") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P61626 <- P61626 %>% 
  reduce(full_join, by = "sequence") %>% 
  select_if(~!all(is.na(.))) %>%
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


Q96CG8 <- list() #Not Identified
for (i in 1:length(all_batches)) {
  Q96CG8[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q96CG8") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}
P06731 <- list() #Not Identified
for (i in 1:length(all_batches)) {
  P06731[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P06731") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


Q8NFJ5 <- list() #Not Identified
for (i in 1:length(all_batches)) {
  Q8NFJ5[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "Q8NFJ5") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 
view(P08246)

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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


P49913 <- list() #Not Identified
for (i in 1:length(all_batches)) {
  P49913[[i]] <- all_batches[[i]] %>%
    filter(Protein.Group.Accessions == "P49913") %>%
    ungroup() %>%
    select(-Protein.Group.Accessions)
}

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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 


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
  as_data_frame() %>%
  center(robust=TRUE, exclude = "sequence") #Median centring 



#Merge, based on modified peptide sequence
#Differential expression analyis
  #Define sig difference in TMT intensities for select PTMs of a protein
    #Bind wilcox.test results to PTMs + BH correction
WilcoxFunc <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = wilcox.test(x, y)
  results$p.value
}




P05783_PTM <- merge(P05783, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P05783_PTMrawpvalue = apply(P05783_PTM, 1, WilcoxFunc, grp1 = colnames(P05783_PTM)[grep("TUMOR",colnames(P05783_PTM))]
, grp2 = colnames(P05783_PTM)[grep("NAT",colnames(P05783_PTM))])
P05783_PTM_T_test <- cbind(P05783_PTM, P05783_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P05783_PTMrawpvalue) 
P05783_PTM_T_test$BH_adjusted_pval = p.adjust(P05783_PTMrawpvalue, 
               method = "BH")
P05783_PTM_T_test$expression_diff <- "No sig. expression difference"
P05783_PTM_T_test$expression_diff[P05783_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P05783_PTM_T_test <- P05783_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') | #phosphorylation
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na()
view(P05783_PTM_T_test)

P61626_PTM <- merge(P61626, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P61626_PTMrawpvalue = apply(P61626_PTM, 1, WilcoxFunc, grp1 = colnames(P61626_PTM)[grep("TUMOR",colnames(P61626_PTM))]
, grp2 = colnames(P61626_PTM)[grep("NAT",colnames(P61626_PTM))])
P61626_PTM_T_test <- cbind(P61626_PTM, P61626_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P61626_PTMrawpvalue) 
P61626_PTM_T_test$BH_adjusted_pval = p.adjust(P61626_PTMrawpvalue, 
               method = "BH")
P61626_PTM_T_test$expression_diff <- "No sig. expression difference"
P61626_PTM_T_test$expression_diff[P61626_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P61626_PTM_T_test <- P61626_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(P61626_PTM_T_test)




P25815_PTM <- merge(P25815, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P25815_PTMrawpvalue = apply(P25815_PTM, 1, WilcoxFunc, grp1 = colnames(P25815_PTM)[grep("TUMOR",colnames(P25815_PTM))]
, grp2 = colnames(P25815_PTM)[grep("NAT",colnames(P25815_PTM))])
P25815_PTM_T_test <- cbind(P25815_PTM, P25815_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P25815_PTMrawpvalue) 
 P25815_PTM_T_test$BH_adjusted_pval = p.adjust(P25815_PTMrawpvalue, 
               method = "BH")
P25815_PTM_T_test$expression_diff <- "No sig. expression difference"
P25815_PTM_T_test$expression_diff[P25815_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

Q12884_PTM <- merge(Q12884, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
Q12884_PTMrawpvalue = apply(Q12884_PTM, 1, WilcoxFunc, grp1 = colnames(Q12884_PTM)[grep("TUMOR",colnames(Q12884_PTM))]
, grp2 = colnames(Q12884_PTM)[grep("NAT",colnames(Q12884_PTM))])
Q12884_PTM_T_test <- cbind(Q12884_PTM, Q12884_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, Q12884_PTMrawpvalue) 
Q12884_PTM_T_test$BH_adjusted_pval = p.adjust(Q12884_PTMrawpvalue, 
               method = "BH")
Q12884_PTM_T_test$expression_diff <- "No sig. expression difference"
Q12884_PTM_T_test$expression_diff[Q12884_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P36952_PTM <- merge(P36952, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    drop_na %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P36952_PTMrawpvalue = apply(P36952_PTM, 1, WilcoxFunc, grp1 = colnames(P36952_PTM)[grep("TUMOR",colnames(P36952_PTM))]
, grp2 = colnames(P36952_PTM)[grep("NAT",colnames(P36952_PTM))])
P36952_PTM_T_test <- cbind(P36952_PTM, P36952_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P36952_PTMrawpvalue) 
 P36952_PTM_T_test$BH_adjusted_pval = p.adjust(P36952_PTMrawpvalue, 
               method = "BH")
P36952_PTM_T_test$expression_diff <- "No sig. expression difference"
P36952_PTM_T_test$expression_diff[P36952_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P36952_PTM_T_test <- P36952_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) 
view(P36952_PTM_T_test)
view(drop_na(P36952))


P06702_PTM <- merge(P06702, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P06702_PTMrawpvalue = apply(P06702_PTM, 1, WilcoxFunc, grp1 = colnames(P06702_PTM)[grep("TUMOR",colnames(P06702_PTM))]
, grp2 = colnames(P06702_PTM)[grep("NAT",colnames(P06702_PTM))])
P06702_PTM_T_test <- cbind(P06702_PTM, P06702_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P06702_PTMrawpvalue) 
 P06702_PTM_T_test$BH_adjusted_pval = p.adjust(P06702_PTMrawpvalue, 
               method = "BH")
P06702_PTM_T_test$expression_diff <- "No sig. expression difference"
P06702_PTM_T_test$expression_diff[P06702_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P06702_PTM_T_test <- P06702_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(P06702_PTM_T_test)


P83881_PTM <- merge(P83881, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P83881_PTMrawpvalue = apply(P83881_PTM, 1, WilcoxFunc, grp1 = colnames(P83881_PTM)[grep("TUMOR",colnames(P83881_PTM))]
, grp2 = colnames(P83881_PTM)[grep("NAT",colnames(P83881_PTM))])
P83881_PTM_T_test <- cbind(P83881_PTM, P83881_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P83881_PTMrawpvalue) 
 P83881_PTM_T_test$BH_adjusted_pval = p.adjust(P83881_PTMrawpvalue, 
               method = "BH")
P83881_PTM_T_test$expression_diff <- "No sig. expression difference"
P83881_PTM_T_test$expression_diff[P83881_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P05109_PTM <- merge(P05109, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P05109_PTMrawpvalue = apply(P05109_PTM, 1, WilcoxFunc, grp1 = colnames(P05109_PTM)[grep("TUMOR",colnames(P05109_PTM))]
, grp2 = colnames(P05109_PTM)[grep("NAT",colnames(P05109_PTM))])
P05109_PTM_T_test <- cbind(P05109_PTM, P05109_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P05109_PTMrawpvalue) 
P05109_PTM_T_test$BH_adjusted_pval = p.adjust(P05109_PTMrawpvalue, 
               method = "BH")
P05109_PTM_T_test$expression_diff <- "No sig. expression difference"
P05109_PTM_T_test$expression_diff[P05109_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P05109_PTM_T_test <- P05109_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  )  
view(P05109_PTM_T_test)



P50454_PTM <- merge(P50454, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P50454_PTMrawpvalue = apply(P50454_PTM, 1, WilcoxFunc, grp1 = colnames(P50454_PTM)[grep("TUMOR",colnames(P50454_PTM))]
, grp2 = colnames(P50454_PTM)[grep("NAT",colnames(P50454_PTM))])
P50454_PTM_T_test <- cbind(P50454_PTM, P50454_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P50454_PTMrawpvalue) 
 P50454_PTM_T_test$BH_adjusted_pval = p.adjust(P50454_PTMrawpvalue, 
               method = "BH")
P50454_PTM_T_test$expression_diff <- "No sig. expression difference"
P50454_PTM_T_test$expression_diff[P50454_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P50454_PTM_T_test <- P50454_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(P50454_PTM_T_test)


P80511_PTM <- merge(P80511, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P80511_PTMrawpvalue = apply(P80511_PTM, 1, WilcoxFunc, grp1 = colnames(P80511_PTM)[grep("TUMOR",colnames(P80511_PTM))]
, grp2 = colnames(P80511_PTM)[grep("NAT",colnames(P80511_PTM))])
P80511_PTM_T_test <- cbind(P80511_PTM, P80511_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P80511_PTMrawpvalue) 
 P80511_PTM_T_test$BH_adjusted_pval = p.adjust(P80511_PTMrawpvalue, 
               method = "BH")
P80511_PTM_T_test$expression_diff <- "No sig. expression difference"
P80511_PTM_T_test$expression_diff[P80511_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
Q9NR99_PTM <- merge(Q9NR99, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
Q9NR99_PTMrawpvalue = apply(Q9NR99_PTM, 1, WilcoxFunc, grp1 = colnames(Q9NR99_PTM)[grep("TUMOR",colnames(Q9NR99_PTM))]
, grp2 = colnames(Q9NR99_PTM)[grep("NAT",colnames(Q9NR99_PTM))])
Q9NR99_PTM_T_test <- cbind(Q9NR99_PTM, Q9NR99_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, Q9NR99_PTMrawpvalue) 
Q9NR99_PTM_T_test$BH_adjusted_pval = p.adjust(Q9NR99_PTMrawpvalue, 
               method = "BH")
Q9NR99_PTM_T_test$expression_diff <- "No sig. expression difference"
Q9NR99_PTM_T_test$expression_diff[Q9NR99_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

Q01650_PTM <- merge(Q01650, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
Q01650_PTMrawpvalue = apply(Q01650_PTM, 1, WilcoxFunc, grp1 = colnames(Q01650_PTM)[grep("TUMOR",colnames(Q01650_PTM))]
, grp2 = colnames(Q01650_PTM)[grep("NAT",colnames(Q01650_PTM))])
Q01650_PTM_T_test <- cbind(Q01650_PTM, Q01650_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, Q01650_PTMrawpvalue) 
Q01650_PTM_T_test$BH_adjusted_pval = p.adjust(Q01650_PTMrawpvalue, 
               method = "BH")
Q01650_PTM_T_test$expression_diff <- "No sig. expression difference"
Q01650_PTM_T_test$expression_diff[Q01650_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

P31949_PTM <- merge(P31949, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P31949_PTMrawpvalue = apply(P31949_PTM, 1, WilcoxFunc, grp1 = colnames(P31949_PTM)[grep("TUMOR",colnames(P31949_PTM))]
, grp2 = colnames(P31949_PTM)[grep("NAT",colnames(P31949_PTM))])
P31949_PTM_T_test <- cbind(P31949_PTM, P31949_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P31949_PTMrawpvalue) 
 P31949_PTM_T_test$BH_adjusted_pval = p.adjust(P31949_PTMrawpvalue, 
               method = "BH")
P31949_PTM_T_test$expression_diff <- "No sig. expression difference"
P31949_PTM_T_test$expression_diff[P31949_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P31949_PTM_T_test <- P31949_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(P31949_PTM_T_test)


#PROBLEM
P40261_PTM <- merge(P40261, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P40261_PTMrawpvalue = apply(P40261_PTM, 1, WilcoxFunc, grp1 = colnames(P40261_PTM)[grep("TUMOR",colnames(P40261_PTM))]
, grp2 = colnames(P40261_PTM)[grep("NAT",colnames(P40261_PTM))])
P40261_PTM_T_test <- cbind(P40261_PTM, P40261_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P40261_PTMrawpvalue) 
 P40261_PTM_T_test$BH_adjusted_pval = p.adjust(P40261_PTMrawpvalue, 
               method = "BH")
P40261_PTM_T_test$expression_diff <- "No sig. expression difference"
P40261_PTM_T_test$expression_diff[P40261_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

P02788_PTM <- merge(P02788, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P02788_PTMrawpvalue = apply(P02788_PTM, 1, WilcoxFunc, grp1 = colnames(P02788_PTM)[grep("TUMOR",colnames(P02788_PTM))]
, grp2 = colnames(P02788_PTM)[grep("NAT",colnames(P02788_PTM))])
P02788_PTM_T_test <- cbind(P02788_PTM, P02788_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P02788_PTMrawpvalue) 
 P02788_PTM_T_test$BH_adjusted_pval = p.adjust(P02788_PTMrawpvalue, 
               method = "BH")
P02788_PTM_T_test$expression_diff <- "No sig. expression difference"
P02788_PTM_T_test$expression_diff[P02788_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P02788_PTM_T_test <- P02788_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(P02788_PTM_T_test)


#TOO MUCH NA
O00469_PTM <- merge(O00469, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat) 
O00469_PTMrawpvalue = apply(O00469_PTM, 1, WilcoxFunc, grp1 = colnames(O00469_PTM)[grep("TUMOR",colnames(O00469_PTM))]
, grp2 = colnames(O00469_PTM)[grep("NAT",colnames(O00469_PTM))])
O00469_PTM_T_test <- cbind(O00469_PTM, O00469_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, O00469_PTMrawpvalue) 
 O00469_PTM_T_test$BH_adjusted_pval = p.adjust(O00469_PTMrawpvalue, 
               method = "BH")
O00469_PTM_T_test$expression_diff <- "No sig. expression difference"
O00469_PTM_T_test$expression_diff[O00469_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
O00469_PTM_T_test$expression_diff[O00469_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P35442_PTM <- merge(P35442, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P35442_PTMrawpvalue = apply(P35442_PTM, 1, WilcoxFunc, grp1 = colnames(P35442_PTM)[grep("TUMOR",colnames(P35442_PTM))]
, grp2 = colnames(P35442_PTM)[grep("NAT",colnames(P35442_PTM))])
P35442_PTM_T_test <- cbind(P35442_PTM, P35442_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P35442_PTMrawpvalue) 
 P35442_PTM_T_test$BH_adjusted_pval = p.adjust(P35442_PTMrawpvalue, 
               method = "BH")
P35442_PTM_T_test$expression_diff <- "No sig. expression difference"
P35442_PTM_T_test$expression_diff[P35442_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
O76021_PTM <- merge(O76021, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
O76021_PTMrawpvalue = apply(O76021_PTM, 1, WilcoxFunc, grp1 = colnames(O76021_PTM)[grep("TUMOR",colnames(O76021_PTM))]
, grp2 = colnames(O76021_PTM)[grep("NAT",colnames(O76021_PTM))])
O76021_PTM_T_test <- cbind(O76021_PTM, O76021_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, O76021_PTMrawpvalue) 
 O76021_PTM_T_test$BH_adjusted_pval = p.adjust(O76021_PTMrawpvalue, 
               method = "BH")
O76021_PTM_T_test$expression_diff <- "No sig. expression difference"
O76021_PTM_T_test$expression_diff[O76021_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P08246_PTM <- merge(P08246, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P08246_PTMrawpvalue = apply(P08246_PTM, 1, WilcoxFunc, grp1 = colnames(P08246_PTM)[grep("TUMOR",colnames(P08246_PTM))]
, grp2 = colnames(P08246_PTM)[grep("NAT",colnames(P08246_PTM))])
P08246_PTM_T_test <- cbind(P08246_PTM, P08246_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P08246_PTMrawpvalue) 
 P08246_PTM_T_test$BH_adjusted_pval = p.adjust(P08246_PTMrawpvalue, 
               method = "BH")
P08246_PTM_T_test$expression_diff <- "No sig. expression difference"
P08246_PTM_T_test$expression_diff[P08246_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P05164_PTM <- merge(P05164, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P05164_PTMrawpvalue = apply(P05164_PTM, 1, WilcoxFunc, grp1 = colnames(P05164_PTM)[grep("TUMOR",colnames(P05164_PTM))]
, grp2 = colnames(P05164_PTM)[grep("NAT",colnames(P05164_PTM))])
P05164_PTM_T_test <- cbind(P05164_PTM, P05164_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P05164_PTMrawpvalue) 
 P05164_PTM_T_test$BH_adjusted_pval = p.adjust(P05164_PTMrawpvalue, 
               method = "BH")
P05164_PTM_T_test$expression_diff <- "No sig. expression difference"
P05164_PTM_T_test$expression_diff[P05164_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P05164_PTM_T_test <- P05164_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(P05164_PTM_T_test)


Q9NR30_PTM <- merge(Q9NR30, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
Q9NR30_PTMrawpvalue = apply(Q9NR30_PTM, 1, WilcoxFunc, grp1 = colnames(Q9NR30_PTM)[grep("TUMOR",colnames(Q9NR30_PTM))]
, grp2 = colnames(Q9NR30_PTM)[grep("NAT",colnames(Q9NR30_PTM))])
Q9NR30_PTM_T_test <- cbind(Q9NR30_PTM, Q9NR30_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, Q9NR30_PTMrawpvalue) 
 Q9NR30_PTM_T_test$BH_adjusted_pval = p.adjust(Q9NR30_PTMrawpvalue, 
               method = "BH")
Q9NR30_PTM_T_test$expression_diff <- "No sig. expression difference"
Q9NR30_PTM_T_test$expression_diff[Q9NR30_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

P30273_PTM <- merge(P30273, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P30273_PTMrawpvalue = apply(P30273_PTM, 1, WilcoxFunc, grp1 = colnames(P30273_PTM)[grep("TUMOR",colnames(P30273_PTM))]
, grp2 = colnames(P30273_PTM)[grep("NAT",colnames(P30273_PTM))])
P30273_PTM_T_test <- cbind(P30273_PTM, P30273_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P30273_PTMrawpvalue) 
 P30273_PTM_T_test$BH_adjusted_pval = p.adjust(P30273_PTMrawpvalue, 
               method = "BH")
P30273_PTM_T_test$expression_diff <- "No sig. expression difference"
P30273_PTM_T_test$expression_diff[P30273_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

Q99715_PTM <- merge(Q99715, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
Q99715_PTMrawpvalue = apply(Q99715_PTM, 1, WilcoxFunc, grp1 = colnames(Q99715_PTM)[grep("TUMOR",colnames(Q99715_PTM))]
, grp2 = colnames(Q99715_PTM)[grep("NAT",colnames(Q99715_PTM))])
Q99715_PTM_T_test <- cbind(Q99715_PTM, Q99715_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, Q99715_PTMrawpvalue) 
 Q99715_PTM_T_test$BH_adjusted_pval = p.adjust(Q99715_PTMrawpvalue, 
               method = "BH")
Q99715_PTM_T_test$expression_diff <- "No sig. expression difference"
Q99715_PTM_T_test$expression_diff[Q99715_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
Q99715_PTM_T_test <- Q99715_PTM_T_test %>%
  filter(expression_diff == "Sig. expression difference") %>%
  filter(
      stringr::str_detect(mod, 'Phospho|phospho') |
      stringr::str_detect(mod_mass, '0.98') | #Citru
      stringr::str_detect(mod, 'Hydroxylation') |#Hyd
      stringr::str_detect(mod_mass, '114.04|114.042927') |#Ub residue
      stringr::str_detect(mod_mass, '383.|383.4460') | #Ub full
      stringr::str_detect(mod_mass, '14.01|14.0266') |  #Methyl
      stringr::str_detect(mod_mass, '42.01|42.0367') #Acetyl
  ) %>% drop_na() 
view(Q99715_PTM_T_test)


#PROBLEM
P52926_PTM <- merge(P52926, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P52926_PTMrawpvalue = apply(P52926_PTM, 1, WilcoxFunc, grp1 = colnames(P52926_PTM)[grep("TUMOR",colnames(P52926_PTM))]
, grp2 = colnames(P52926_PTM)[grep("NAT",colnames(P52926_PTM))])
P52926_PTM_T_test <- cbind(P52926_PTM, P52926_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P52926_PTMrawpvalue) 
 P52926_PTM_T_test$BH_adjusted_pval = p.adjust(P52926_PTMrawpvalue, 
               method = "BH")
P52926_PTM_T_test$expression_diff <- "No sig. expression difference"
P52926_PTM_T_test$expression_diff[P52926_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P52926_PTM_T_test$expression_diff[P52926_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

Q8NCL4_PTM <- merge(Q8NCL4, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
Q8NCL4_PTMrawpvalue = apply(Q8NCL4_PTM, 1, WilcoxFunc, grp1 = colnames(Q8NCL4_PTM)[grep("TUMOR",colnames(Q8NCL4_PTM))]
, grp2 = colnames(Q8NCL4_PTM)[grep("NAT",colnames(Q8NCL4_PTM))])
Q8NCL4_PTM_T_test <- cbind(Q8NCL4_PTM, Q8NCL4_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, Q8NCL4_PTMrawpvalue) 
 Q8NCL4_PTM_T_test$BH_adjusted_pval = p.adjust(Q8NCL4_PTMrawpvalue, 
               method = "BH")
Q8NCL4_PTM_T_test$expression_diff <- "No sig. expression difference"
Q8NCL4_PTM_T_test$expression_diff[Q8NCL4_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
Q8NCL4_PTM_T_test$expression_diff[Q8NCL4_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
O00425_PTM <- merge(O00425, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
O00425_PTMrawpvalue = apply(O00425_PTM, 1, WilcoxFunc, grp1 = colnames(O00425_PTM)[grep("TUMOR",colnames(O00425_PTM))]
, grp2 = colnames(O00425_PTM)[grep("NAT",colnames(O00425_PTM))])
O00425_PTM_T_test <- cbind(O00425_PTM, O00425_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, O00425_PTMrawpvalue) 
 O00425_PTM_T_test$BH_adjusted_pval = p.adjust(O00425_PTMrawpvalue, 
               method = "BH")
O00425_PTM_T_test$expression_diff <- "No sig. expression difference"
O00425_PTM_T_test$expression_diff[O00425_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
O00425_PTM_T_test$expression_diff[O00425_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P80188_PTM <- merge(P80188, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P80188_PTMrawpvalue = apply(P80188_PTM, 1, WilcoxFunc, grp1 = colnames(P80188_PTM)[grep("TUMOR",colnames(P80188_PTM))]
, grp2 = colnames(P80188_PTM)[grep("NAT",colnames(P80188_PTM))])
P80188_PTM_T_test <- cbind(P80188_PTM, P80188_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P80188_PTMrawpvalue) 
 P80188_PTM_T_test$BH_adjusted_pval = p.adjust(P80188_PTMrawpvalue, 
               method = "BH")
P80188_PTM_T_test$expression_diff <- "No sig. expression difference"
P80188_PTM_T_test$expression_diff[P80188_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P80188_PTM_T_test$expression_diff[P80188_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#PROBLEM
P62760_PTM <- merge(P62760, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P62760_PTMrawpvalue = apply(P62760_PTM, 1, WilcoxFunc, grp1 = colnames(P62760_PTM)[grep("TUMOR",colnames(P62760_PTM))]
, grp2 = colnames(P62760_PTM)[grep("NAT",colnames(P62760_PTM))])
P62760_PTM_T_test <- cbind(P62760_PTM, P62760_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P62760_PTMrawpvalue) 
 P62760_PTM_T_test$BH_adjusted_pval = p.adjust(P62760_PTMrawpvalue, 
               method = "BH")
P62760_PTM_T_test$expression_diff <- "No sig. expression difference"
P62760_PTM_T_test$expression_diff[P62760_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P62760_PTM_T_test$expression_diff[P62760_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"

#TOO MUCH NA
P11388_PTM <- merge(P11388, ALL_seq, by="sequence", all.x=TRUE) %>%
  unique()%>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, starts_with("TUM"), starts_with("NAT"), count) %>%
    mutate(Mean_tumor = rowMeans(.[grep("^TUM", names(.))])) %>%
    mutate(Mean_nat = rowMeans(.[grep("^NAT", names(.))])) %>%
    mutate(Difference_Tum_vs_Nat = Mean_tumor - Mean_nat)
P11388_PTMrawpvalue = apply(P11388_PTM, 1, WilcoxFunc, grp1 = colnames(P11388_PTM)[grep("TUMOR",colnames(P11388_PTM))]
, grp2 = colnames(P11388_PTM)[grep("NAT",colnames(P11388_PTM))])
P11388_PTM_T_test <- cbind(P11388_PTM, P11388_PTMrawpvalue) %>%
  select(sequence, sequence_no_mod, 
    mod, mod_mass, count, Mean_tumor, Mean_nat, Difference_Tum_vs_Nat, P11388_PTMrawpvalue) 
 P11388_PTM_T_test$BH_adjusted_pval = p.adjust(P11388_PTMrawpvalue, 
               method = "BH")
P11388_PTM_T_test$expression_diff <- "No sig. expression difference"
P11388_PTM_T_test$expression_diff[P11388_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"
P11388_PTM_T_test$expression_diff[P11388_PTM_T_test$BH_adjusted_pval < 0.05] <- "Sig. expression difference"











