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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab")
getwd() 
list.files(wd)
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/ALL_mzTab"))
(file_names_short <- substring(file_name_long, 39, 46)) 
length(file_names_short)

#All TMT intensities
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
TMT_f03_03.05.22 <- readRDS(file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/03.05.22_TMT_f03") #One lost file
    #Combine into one list
mzML <- c(TMT_06.04.22, TMT_20.04.22, TMT_22.04.22, TMT_28.04.22, TMT_part1_03.05.22, TMT_part2_03.05.22, 
        TMT_part3_03.05.22, TMT_part4_03.05.22, TMT_part5_03.05.22, TMT_part6_03.05.22, TMT_part7_03.05.22, TMT_f03_03.05.22)
names(mzML) #B1S3_f02 is a duplicate
length(mzML) #Should be 264
n_distinct(mzML)

#Reorder based on mzTab file order
    #Add filename column
    #Add index column
order_mzML <- list()
namemzml <- names(mzML)
for (i in seq_along(mzML)) {
    order_mzML[[i]] <- tibble(index=rownames(mzML[[i]][, 0])) %>%
    mutate(index = trimws(str_remove_all(index, "F1.S"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    mutate(index = trimws(str_remove_all(index, "^0"))) %>%
    select(index) %>%
    cbind(mzML[[i]]) %>% 
    as_tibble %>%
    add_column(namemzml = namemzml[[i]], .before = "126") 

}
view(order_mzML[[1]]) 

    #Arrange with new column
all <- bind_rows(order_mzML)  %>%
    distinct() #Removes the duplicate
arrange(all, namemzml) #arrange in alphabetical order
    #Split file back up, in right order
TMT_264 <- split(all, f = all$namemzml)  
TMT_264 <- set_names(TMT_264, file_names_short)
    #Checking if it worked
TMT_264[[1]]
TMT_264[[264]]
TMT_264$B1S3_f02 #namemzml matches
TMT_264$B1S3_f03 #namemzml matches
TMT_264[[265]] #Does not excist = GOOD!
length(TMT_264) #264

#Save
saveRDS(TMT_264, file = "/Users/jensvandeperre/Desktop/Outputs/TMTs/ALL_TMTs_16.05.22")

