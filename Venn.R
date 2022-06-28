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
library("ggvenn")
library("ggVennDiagram")
library("VennDiagram")
library("venn")


install.packages("venn", repos="https://www.freestatistics.org/cran/")

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
#Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"))
(file_names_short <- substring(file_name_long, 39, 46)) 
length(file_names_short)

#Original study psms
#Load in original studies
(file_paths_ori_study <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/Orig_study_PSMs"))
ori_study <- list()
for (i in seq_along(file_paths_ori_study)) {
  ori_study[[i]] <- fread(file_paths_ori_study[[i]])
}
str(ori_study[[1]])
length(ori_study)
#select indexes
OS_index <- list()
for (i in 1:264) {
  OS_index[[i]] <- ori_study[[i]] %>%
    as_tibble() %>%
    select(FileName)
}


#Select sequences
orig <- list()
for (i in 1:264) {
  orig[[i]] <- ori_study[[i]] %>%
    as_tibble() %>%
    mutate(sequence = trimws(str_remove_all(RTAtPrecursorHalfElution, "[[:digit:]]+"))) %>%
    mutate(sequence = trimws(str_remove_all(sequence, "[+.]"))) %>%
    mutate(PSM = paste(sequence, FileName, sep="_")) %>%
    select(PSM)
}
view(orig[[1]])
#Turn into 1 dataframe
OS_PSM <- bind_rows(orig)
OS_PSM_unique <- OS_PSM %>% distinct()
OS_index <- bind_rows(OS_index)
nrow(OS_PSM) #over 3 million peptides identified
nrow(OS_PSM_unique) #Over 20000 unique peptides identified
#Save as input for online Venny tool
write.table(OS_PSM, file="/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/Original_study/Original_study.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(OS_PSM_unique, file="/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/Original_study/Original_study_unique.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)


#Load PSMs 
PSM_all <- readRDS(file = "~/Desktop/Outputs/PSMs/ALL_PSMs_4.5.22") 
view(PSM_all[[1]]) 
length(PSM_all)

AS_index <- list() #empty list
for (i in seq_along(PSM_all)) {
  AS_index[[i]] <- PSM_all[[i]] %>%
    as_tibble() %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 controllerNumber=1 scan="))) %>% 
    select(index)
}

psm <- list() #empty list
for (i in seq_along(PSM_all)) {
  psm[[i]] <- PSM_all[[i]] %>%
    as_tibble() %>% 
    mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 controllerNumber=1 scan="))) %>% 
    mutate(PSM = paste(sequence_no_mod, index, sep="_")) %>%
    select(PSM)
}
view(psm[[1]])

AS_PSM <- bind_rows(psm)
AS_PSM_unique <- AS_PSM %>% distinct()
AS_index <- bind_rows(AS_index)
nrow(AS_PSM) #Over 4.5 million peptides identified
nrow(AS_PSM_unique) #Over 20000 unique peptides identified
#Save as input for online Venny tool
write.table(AS_PSM, file="/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/ANN_SoLo/ANN_SoLo_no_mod.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(AS_PSM_unique, file="/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/ANN_SoLo/ANN_SoLo_no_mod_unique.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)

write.table(AS_index, file="/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/AS_indexe.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(OS_index, file="/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/OS_index.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)


AS <- fread("/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/ANN_SoLo/ANN_SoLo_no_mod.txt", header=FALSE) %>%
  as_tibble 
view(head(AS))


OS <- fread("/Users/jensvandeperre/Desktop/Inputs/Venn_diagram/Original_study/Original_study.txt", header=FALSE) %>%
  as_tibble 
view(head(AS))


AS <- list(AS)
OS <- list(OS)
view(AS)

AS <- as.data.frame(AS)
OS <- as.data.frame(OS)

venn.diagram(
  x = list(AS, OS),
  category.names = c("ANN-SoLo", "Original Study"),
  filename = '/Users/jensvandeperre/Desktop/Outputs/Plots/venn_diagramm.png',
  output=TRUE
)

ANN_vs_ORG_unique <- list(
  `Reanalysis` = AS_index,
  `Original Study` = OS_index
)
Venn <- ggvenn(ANN_vs_ORG_unique, c("Reanalysis", "Original Study"),
                   fill_color = c("blue", "yellow"),
                   stroke_color = FALSE) +
  labs(title="Comparison Identified Peptide Sequences") +
  theme_void() +
  theme(plot.title = element_text(size = 20, face = "bold"))
Venn
png(file = "/Users/jensvandeperre/Desktop/Outputs/Plots/Seq_Venn.png")
Venn
dev.off()
