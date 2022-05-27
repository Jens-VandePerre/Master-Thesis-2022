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


#TMT spectra, count rows
TMT <- readRDS("/Users/jensvandeperre/Desktop/Outputs/TMTs/ALL_TMTs_16.05.22")
spec_count <- list()
for (i in 1:264) {
    spec_count[[i]] <- TMT[[i]] %>% 
    n_distinct()
}

#PSMs from ANN-SoLo
AS <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSMs/ALL_PSMs_4.5.22")
AS_count <- list()
for (i in 1:264) {
    AS_count[[i]] <- ANN_SoLo[[i]] %>% 
    nrow()
}

#PSMs original study
(file_paths_Original <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/"))
ORIG <- list()
for (i in 1:264) {
  ORIG[[i]] <- read.csv(file_paths_Original[[i]], header = FALSE, sep = ",")
}
OS_count <- list()
for (i in 1:264) {
    OS_count[[i]] <- ORIG[[i]] %>% 
    nrow()
}

#ANN-Solo identification %
AS_IP <- list()
for(i in 1:264) {
    AS_IP[[i]] <- paste((AS_count[[i]]/TMT[[i]])*100, "%")
}
mean(AS_IP)
median(AS_IP)

#ANN-Solo identification %
OS_IP <- list()
for(i in 1:264) {
    OS_IP[[i]] <- (OS_count[[i]]/TMT[[i]])*100
}
mean(OS_IP)
median(OS_IP)


perc <- list() #empty list
for (i in seq_along(Identified_peptides)) { #Only for the first 6 spectra, 
  perc[[i]] <- (Identified_peptides[[i]]/Spectral_count[[i]])*100
}
Identification_percentage <- set_names(perc, file_names_short) #names each file by file_names_short
Identification_percentage
  #Making tibble of Identification_percentage, for plot
df_perc <- matrix(unlist(Identification_percentage), nrow=length(Identification_percentage), byrow=TRUE)
tbl_Identification_percentage <- tibble(File_Name=file_names_short , Identification_Percentage=df_perc)
tbl_Identification_percentage
  #Plot
p_perc <- ggplot(tbl_Identification_percentage, aes(x= File_Name , y= Identification_Percentage)) +
  geom_col() +
  labs(x="File Name", y="Percentage Identified Spectra (%)", title="Identification Percentage" , 
        subtitle="Percentage of identified peptides for each file with spectra") +
  geom_text(aes(label=round(df_perc, digits = 1)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 18))
  #Print p_perc
p_perc
pdf(file = "~/Desktop/Outputs/Plots/Identification_Percentage.pdf")
   p_perc
dev.off()
