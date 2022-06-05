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
PTM_count <- bind_rows(PTM) %>%
filter(mod != "No direct match found in Unimod")
nrow(PTM_count)
PTM_count <- 1665076



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
  #ANN-Solo identification %
  ANN_SoLo <- (AS_count/TMT_count)*100
  #Original Study identification %
  Original_study <- (OS_count/TMT_count)*100
  #Making tibble 
tbl_Identification_percentage <- tibble(Study=c("Original Study", "ANN-SoLo") , Identification_Percentage=c(Original_study, ANN_SoLo))
tbl_Identification_percentage
  #Plot
pep_IP <- ggplot(tbl_Identification_percentage, aes(x= Study , y= Identification_Percentage)) +
  geom_col(width = 0.45, fill = c("#56B4E9", "#E69F00")) +
  labs(x="Study", y="Percentage Identified Spectra (%)", title="Peptide Identification Percentage" , 
        subtitle="Comparing ANN-SoLo and the Original Study") +
  geom_text(aes(label=round(c(Original_study, ANN_SoLo), digits = 1)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 3.5) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18)) 
  #Print pep_IP
pep_IP
png(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Peptide_Identification_Percentage.png")
   pep_IP
dev.off()

#Total amount of peptides identified
  #Making tibble
OS_count
AS_count
DIFF
PTM_count
AS_count_no_PTM <- AS_count-PTM_count
Study=c("Original Study", "ANN-SoLo")

Study <- c("Original Study", "ANN-SoLo")
Non_modified <- c(OS_count, AS_count_no_PTM)
Modified <- c(NA,PTM_count)
df_IT <- data.frame(Study, Non_modified, Modified)
tbl_Identification_total <- pivot_longer(df_IT, Non_modified:Modified, names_to = "Spectra_Type", values_to = "Identification_count")
tbl_Identification_total


  #Plot
pep_IT <- ggplot(tbl_Identification_total, 
  aes(fill=Spectra_Type, y=Identification_count, x=Study, label=Identification_count)) + 
    geom_bar(position="stack", stat="identity", width = 0.65) +
    geom_text(size = 3.5, position = position_stack(vjust = 0.5)) +
    scale_color_manual(labels=c("Non-modified","Modified")) +
    labs(x = "Study", y = "Identified Peptides", 
        title = "Total Identified Peptides", 
        subtitle = "Comparing ANN-SoLo and the Original Study") +
    labs(fill='Spectra Type') +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18))
  #Print pep_IP
pep_IT
png(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Total_Identified_Peptides.png")
   pep_IT
dev.off()

#Total amount of proteins identified
  #Counts unique proteins
  #TMT global proteomic analysis of the 96 tumor and NAT pairs identified a total of 8,067 proteins
    AS_prot <- read.csv(file="/Users/jensvandeperre/Desktop/Inputs/Limm_Qvalve/data_input.txt", sep="\t")
    view(head(AS_prot))
    AS_prot_distinct <- AS_prot %>% 
    select(Protein.Group.Accessions) %>%
    distinct() %>%
    nrow()

OS_pro_count <- 8067
AS_pro_count <- 8688

tbl_pro_Identification_total <- tibble(Study=c("Original Study", "ANN-SoLo") , Identification_count=c(OS_pro_count, AS_pro_count)) %>%
  mutate(Difference = AS_pro_count - OS_pro_count)
tbl_pro_Identification_total
  #Plot
pro_IT <- ggplot(tbl_pro_Identification_total, aes(x= Study , y= Identification_count)) +
  geom_col(width = 0.45, fill = c("#56B4E9", "#E69F00")) +
  labs(x = "Study", y = "Identified Proteins",
        title = "Unique Identified Proteins",
        subtitle = "Comparing ANN-SoLo and Original Study ") +
  geom_text(aes(label=c(OS_pro_count, AS_pro_count)),
        position = position_dodge(width = 0.9), vjust = -0.25, size = 3.5) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18))
  #Print pep_IP
pro_IT
png(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Unique_Identified_Proteins.png")
   pro_IT
dev.off()

#Comparing amount of unique peptides found

#Most common PTM

#Most found protein

    #Identification percentage of proteins
      #ANN-Solo identification %
      ANN_SoLo <- (AS_count/TMT_count)*100
      #Original Study identification %
      Original_study <- (OS_count/TMT_count)*100
      #Making tibble 
    tbl_Identification_percentage <- tibble(Study=c("ANN-SoLo", "Original Study") , Identification_Percentage=c(ANN_SoLo, Original_study))
    tbl_Identification_percentage
      #Plot
    pro_IP <- ggplot(tbl_Identification_percentage, aes(x= Study , y= Identification_Percentage)) +
      geom_col(width = 0.45, fill = c("#E69F00", "#56B4E9")) +
      labs(x="Study", y="Percentage Identified Proteins (%)", title="Protein Identification Percentage" , 
            subtitle="Comparing the protein identification percentages of ANN-SoLo and the Original Study") +
      geom_text(aes(label=round(Identification_Percentage, digits = 1)), 
                  position=position_dodge(width=0.9), vjust=-0.25, size = 3.5) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 18)) 
      
      #Print pep_IP
    pro_IP
    pdf(file = )
      pro_IP
    dev.off()
