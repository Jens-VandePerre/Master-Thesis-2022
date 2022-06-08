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

#Input 
    #Proteins with their TMT intensities, labelled as TUMOR_batch and NAT_batch
    #The intensiteis are relative intensities:
        #Sum of all intensities from 1 sample (= one TMT channel), devided by reference sample
    #All NAs are removed, to do t.test
    #2819 quantified proteins remain
dat_col_ordered <- fread("/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/log2FC_input.txt")
view(dat_col_ordered)
dim(dat_col_ordered)

#Function to do t.test, comparing Tumor against NAT
ttestFunc <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
rawpvalue = apply(dat_col_ordered, 1, ttestFunc, grp1 = c(2:99), grp2 = c(100:198))
hist(rawpvalue)

#Log2 of tumor TMT relative intensities
    #followed by zero-centering
tum <- log2(dat_col_ordered %>%
    select("TUMOR_127C_B1S1_f01_f12":"TUMOR_126_B5S6_f01_f12")) %>%
    scale(scale = FALSE) 
dim(tum)

#Log2 of NAT TMT relative intensities
    #followed by zero-centering
nat <- log2(dat_col_ordered %>%
    select("NAT_126_B1S1_f01_f12":"NAT_130N_B5S6_f01_f12")) %>%
    scale(scale = FALSE) 
dim(nat)

#Median of log2 transformed relative intesities
Tum = apply(tum, 1, median)
NAT = apply(nat, 1, median) 

#Subtraction for log2FC
foldchange <- Tum - NAT 
hist(foldchange, xlab = "log2 Fold Change (NAT vs Tumor)")
view(foldchange)

#Log2FC is plotted against p.value from t.test
results = cbind(foldchange, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)
volcano = ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
volcano + geom_point()

view((dat_col_ordered))
view((results))

#Log2FC is attached to its Protein.Group.Accessions
FC <- cbind(dat_col_ordered, results) %>%
    select(Protein.Group.Accessions, foldchange, rawpvalue)
FC$BH = p.adjust(FC$rawpvalue, 
               method = "BH")
view(FC)
dim(FC)

#BH filtering of p.values
FC_BH_filtered <- FC %>%
filter(rawpvalue < BH)
view(FC_BH_filtered)
dim(FC_BH_filtered) #Only 1 removed???

























mydat <- readRDS("/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/data_before_name_change")
view(head(mydat[[1]]))
str(mydat[[1]])

B1S1_f01_f12 <- mydat[1:12]
B1S1_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S1_f01_f12_renamed[[i]] <- B1S1_f01_f12[[i]] %>%
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
B1S2_f01_f12 <- mydat[13:24]
B1S2_f01_f12_renamed <- list()
for (i in 1:12) {
  B1S2_f01_f12_renamed[[i]] <- B1S2_f01_f12[[i]] %>%
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  select(Protein.Group.Accessions, Protein.Descriptions, "126":"131") %>%
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
  group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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
    group_by(Protein.Group.Accessions) %>% 
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

colorder <- c("Protein.Group.Accessions",
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
"NAT_161_B1S4_f01_f12", 
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
dat_col_ordered <- dat[, colorder] %>%
    as_data_frame() %>%
    drop_na()

fwrite(dat_col_ordered, "/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/log2FC_input.txt")
dat_col_ordered <- fread("/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/log2FC_input.txt")
view(dat_col_ordered)
dim(dat_col_ordered)

ttestFunc <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
rawpvalue = apply(dat_col_ordered, 1, ttestFunc, grp1 = c(2:99), grp2 = c(100:198))
hist(rawpvalue)

tum <- log2(dat_col_ordered %>%
    select("TUMOR_127C_B1S1_f01_f12":"TUMOR_126_B5S6_f01_f12")) %>%
    scale(scale = FALSE) 
dim(tum)

nat <- log2(dat_col_ordered %>%
    select("NAT_126_B1S1_f01_f12":"NAT_130N_B5S6_f01_f12")) %>%
    scale(scale = FALSE) 
dim(nat)


Tum = apply(tum, 1, median)
NAT = apply(nat, 1, median) 

foldchange <- Tum - NAT 
hist(foldchange, xlab = "log2 Fold Change (NAT vs Tumor)")
view(foldchange)

results = cbind(foldchange, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)

library(ggplot2)
volcano = ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
volcano + geom_point()

FC <- merge(dat_col_ordered, results, by=0) %>%
    select(Protein.Group.Accessions, foldchange, rawpvalue, BH)
FC$BH = p.adjust(FC$rawpvalue, 
               method = "BH")
view(FC)
dim(FC)

FC_BH_filtered <- FC %>%
filter(rawpvalue < BH)
view(FC_BH_filtered)
dim(FC_BH_filtered) #Only 1 removed???
