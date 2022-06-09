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

#Load proteins
mydat <- readRDS("/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/data_before_name_change")
view(head(mydat[[1]]))
str(mydat[[1]])
B1S1_f01_f12 <- mydat[1:12]
B1S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S1_f01_f12_renamed[[i]] <- B1S1_f01_f12[[i]] %>%
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
  as_tibble() %>%
  rename(NAT_161_B1S4_f01_f12 = "126") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Sequence, Protein.Descriptions, "126":"131") %>%
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

dat_B1S1_f01_f12 <- bind_rows(B1S1_f01_f12_renamed) %>%
  group_by(Protein.Group.Accessions, Sequence) %>% 
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
nrow(dat_B1S1_f01_f12) 
view(dat_B1S1_f01_f12) 

dat_B1S2_f01_f12 <- bind_rows(B1S2_f01_f12_renamed) %>%
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
    summarise(
        NAT_161_B1S4_f01_f12 = sum(NAT_161_B1S4_f01_f12),
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
  mutate(NAT_161_B1S4_f01_f12 = NAT_161_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
    group_by(Protein.Group.Accessions, Sequence) %>% 
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
dat <- all_batches %>% 
  reduce(full_join, by='Protein.Group.Accessions')
dim(dat)
view(dat)


#Load PTMs + peptide seq