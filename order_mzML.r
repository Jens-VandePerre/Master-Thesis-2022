library("rpx")
library("mzR")
library("OrgMassSpecR")
library("biomaRt")
library("Hmisc")
library("gplots")
library("limma")
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
library("purrr")
library("rlist")


mzTab_WD <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd()
(file_names_wd <- list.files(mzTab_WD)) #6 mzML files
    #The wanted files
(file_names_short <- substring(fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"), 86, 93)) #Character 86 untill 93 are unique


saveRDS(TMT_part7_03.05.22, file = "~/Desktop/Outputs/mzML_imported/03.05.22_mzML_part7")



TMT_06.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/06.04.22_TMT")
TMT_20.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/20.04.22_TMT")
TMT_22.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/22.04.22_TMT")
TMT_28.04.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/28.04.22_TMT")
TMT_part1_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part1")
TMT_part2_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part2")
TMT_part3_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part3")
TMT_part4_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part4")
TMT_part5_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part5")
TMT_part6_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part6")
TMT_part7_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_part7")
TMT_f03_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_")

mzML <- c(TMT_06.04.22, TMT_20.04.22, TMT_22.04.22, TMT_28.04.22, TMT_part1_03.05.22, TMT_part2_03.05.22, 
        TMT_part3_03.05.22, TMT_part4_03.05.22, TMT_part5_03.05.22, TMT_part6_03.05.22, TMT_part7_03.05.22)
names(mzML)
distinct(mzML)
length(mzML)
file_names_short

mzTab <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/mzTabs_imported/ALL_mzTab_4.5.22")
names(mzTab)
length(mzTab)


#Reorder based on mzTab file order
    #Add filename column
order_mzML <- list()
namemzml <- names(mzML)
for (i in seq_along(mzML)) {
    order_mzML[[i]] <- mzML[[i]] %>% as_tibble() %>%
    add_column(namemzml= namemzml[[i]], .before = "126") 
}
view(order_mzML[[1]]) 
length(order_mzML)
view(mzML[[1]])


view(sort(namemzml))




    #Arrange with new column
all <- bind_rows(order_mzML)
arrange(all, namemzml)
all #264 files put in one file

all %>% select(namemzml) %>%
distinct()



    #Split file back up, in right order
data_list <- split(all, f = all$namemzml)                    
data_list[[1]]
data_list[[263]]
data_list[[12]]
data_list[[88]]
length(data_list)
length(namemzml)



