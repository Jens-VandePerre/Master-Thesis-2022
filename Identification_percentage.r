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

#TMT spectra, count rows
TMT <- readRDS("/Users/jensvandeperre/Desktop/Outputs/TMTs/ALL_TMTs_16.05.22")
view(TMT[[1]])
spec_list_count <- list()
for (i in 1:264) {
    spec_list_count[[i]] <- TMT[[i]] %>% 
    n_distinct()
}
TMT_count <- Reduce("+", spec_list_count) #9547320
TMT_count <- 9547320

#PSMs from ANN-SoLo
AS <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
view(AS[[1]])
AS_list_count <- list()
for (i in 1:264) {
    AS_list_count[[i]] <- AS[[i]] %>% 
    n_distinct()
}
AS_count <- Reduce("+", AS_list_count) #4652819
AS_count <- 4652819

#PTMs
#Load identified PTMs
(ptm_filepaths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PTM_identification_tol_10"))
PTM <- list()
for (i in 1:264) {
  PTM[[i]] <- read.csv(ptm_filepaths[[i]], sep = ",", header = TRUE)
                      }
view(PTM[[1]])
nrow(PTM[[1]])
length(PTM)
PTM_count <- bind_rows(PTM)
nrow(PTM_count)
PTM_count <- 2464841



#PSMs original study
  #Loading in files
ORIG_PSM <- readRDS(file = "/Users/jensvandeperre/Desktop/Inputs/Original_study_saved_list/Original_study_PSMs")
view(ORIG_PSM[[1]])
  #Counting rows
OS_list_count <- list()
for (i in 1:264) {
    OS_list_count[[i]] <- ORIG_PSM[[i]] %>% 
    nrow() 
}
OS_count <- Reduce("+", OS_list_count) #3441708
OS_count <- 3441708


DIFF <- AS_count - OS_count

#Plot Peptide Identification Percentage
  #input tibble
AS_count_no_PTM <- AS_count-PTM_count
Study_IP <- c("Original Study", "ANN-SoLo")
Unmodified_IP <- c(OS_count, AS_count_no_PTM)
Modified_IP <- c(NA, PTM_count)
Total_spectra <- TMT_count
df_IP <- data.frame(Study_IP, Unmodified_IP, Modified_IP, Total_spectra) %>%
  mutate(Modified = Modified_IP/Total_spectra*100) %>%
  mutate(Unmodified = Unmodified_IP/Total_spectra*100)
tbl_Identification_IP <- pivot_longer(df_IP, Modified:Unmodified, names_to = "Spectra_Type", 
    values_to = "Identification_percentage") %>%
      mutate(Label = c(NA, "36.05%", "25.82%", "22.92%")) 
view(tbl_Identification_IP)
str(tbl_Identification_IP)
  #Plot
pep_IP <- ggplot(tbl_Identification_IP, 
    aes(fill=Spectra_Type, y=Identification_percentage, x=Study_IP, label=Label)) +
  geom_bar(position="stack", stat="identity", width = 0.65) +
  geom_text(size = 3.5, position = position_stack(vjust = 0.5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_color_manual(labels=c("Non-modified","Modified")) +
  labs(x="Study", y="Percentage Identified Spectra") +
  labs(fill="Spectra Type") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 12))
  #Print pep_IP
pep_IP
pdf(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Peptide_Identification_Percentage.pdf")
   pep_IP
dev.off()

#Total amount of peptides identified
  #Making tibble
OS_count
AS_count
DIFF
PTM_count
AS_count_no_PTM <- AS_count-PTM_count

Study <- c("Original Study", "ANN-SoLo")
Unmodified <- c(OS_count, AS_count_no_PTM)
Modified <- c(NA,PTM_count)
df_IT <- data.frame(Study, Unmodified, Modified)
tbl_Identification_total <- pivot_longer(df_IT, Unmodified:Modified, names_to = "Spectra_Type", 
    values_to = "Identification_count")
tbl_Identification_total
  #Plot
pep_IT <- ggplot(tbl_Identification_total, 
  aes(fill=Spectra_Type, y=Identification_count, x=Study, label=Identification_count)) + 
    geom_bar(position="stack", stat="identity", width = 0.65) +
    geom_text(size = 3.5, position = position_stack(vjust = 0.5)) +
    scale_color_manual(labels=c("Non-modified","Modified")) +
    labs(x = "Study", y = "Identified Peptides") +
    labs(fill='Spectra Type') +
    theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 12))
  #Print pep_IP
pep_IT
pdf(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Total_Identified_Peptides.pdf")
   pep_IT
dev.off()

#Total amount of proteins identified
  #Counts unique proteins
          #ALL ANN-SoLo + PIA identifiec proteins
    AS_prot <- fread(file="/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/data_input.txt", sep="\t")
    view(head(AS_prot))
    AS_prot_distinct <- AS_prot %>% 
    select(Protein.Group.Accessions) %>%
    unique() 
    dim(AS_prot)
    view(AS_prot)
        #ALL original study identifiec proteins
    OS_prot_NAT <- fread(file="/Users/jensvandeperre/Desktop/Inputs/Original_Proteins/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Normal_TMT_UnsharedLogRatio.cct.txt") %>%
            select(attrib_name)
    OS_prot_TUM <- fread(file="/Users/jensvandeperre/Desktop/Inputs/Original_Proteins/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct.txt") %>%
            select(attrib_name)
view(head(OS_prot_NAT))
OS_prot_distinct <- rbind(OS_prot_NAT, OS_prot_TUM) %>%
    unique() %>%
    nrow()
        #Unique protein counts
OS_pro_count <- 8067
AS_pro_count <- 8688

tbl_pro_Identification_total <- tibble(Study=c("Original Study", "ANN-SoLo") , Identification_count=c(OS_pro_count, AS_pro_count)) %>%
  mutate(Difference = AS_pro_count - OS_pro_count)
tbl_pro_Identification_total
  #Plot
pro_IT <- ggplot(tbl_pro_Identification_total, aes(x= Study , y= Identification_count)) +
  geom_col(width = 0.45, fill = "#0071b2") +
  labs(x = "Study", y = "Identified Proteins") +
  geom_text(aes(label=c(OS_pro_count, AS_pro_count)),
        position = position_dodge(width = 0.9), vjust = -0.25, size = 3.5) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 12))
pro_IT
pdf(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Unique_Identified_Proteins.pdf")
   pro_IT
dev.off()

#PTM identification
mod <- 1665076
mod_notfound <- 2464841 - 1665076
mod + mod_notfound
mod_notfound/2464841
ptm_label <- c("Indentified", "Not identified")
tbl_PTM <- tibble(Modfied_Peptide_Sepctra="Modfied Peptide Sepctra" , Spectral_Count=c(mod, mod_notfound), ptm_label) 
tbl_PTM
  #Plot
PLOTPTM <- ggplot(tbl_PTM, aes(x = Modfied_Peptide_Sepctra, y = Spectral_Count, fill = ptm_label)) +
  geom_col(position = "stack", width = 0.65) +
  labs(y = "Spectra Count", x="") +
  labs(fill='Post Translational Modifications') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  geom_text(label=c("1665076 (67.55%)", "799765 (32.45%)"), size = 3.5, vjust = 5) 
PLOTPTM
pdf(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/PTM_idif.pdf")
   PLOTPTM
dev.off()
