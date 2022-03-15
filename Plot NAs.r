library("rpx")
library(mzR)
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


########################################################################################
#Load in outputs directly, not running loops
############################################
#1. Output intensities for ALL SPECTRA
      #not removing NAs
      #no imputation
      #no normalisation
TMT_Intensities1_10 <- readRDS(file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10")
TMT_Intensities1_10 

#2. Output intensities for ALL spectra
    #impute: method="MLE"
    #no normalisation
TMT_Intensities1_10_Imputation <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation")
TMT_Intensities1_10_Imputation

#3. Loop that outputs intensities for ALL spectra
    #no imputation
    #normalise: method="center.median"
TMT_Intensities1_10_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Normalization")
TMT_Intensities1_10_Normalization  

#4. Output intensities for ALL spectra
    #impute: method="MLE"
    #normalise: method="center.median"
TMT_Intensities1_10_Imputation_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation+Normalization")
TMT_Intensities1_10_Imputation_Normalization
########################################################################################

wd <- setwd("~/Desktop/Read raw file/Data mzML")
getwd() 
(file_names_wd <- list.files(wd)) #The first 10 mzML files of CPTAC
# Loading in multiple mzML files
file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML")
file_paths #The first 10 mzML files paths of CPTAC


         #1. Check missing data before imputation
missing1 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing1[[i]] <- sum(is.na(TMT_Intensities1_10[[i]]))}
missing_tot1 <- set_names(missing1, file_names_wd) #names each file by file_names_wd
missing_tot1 # Total missing for each file


          #2. Check missing data after imputation
missing2 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation)) {
   missing2[[i]] <- sum(is.na(TMT_Intensities1_10_Imputation[[i]]))}
missing_tot2 <- set_names(missing2, file_names_wd) #names each file by file_names_wd
missing_tot2 # Total missing for each file


          #3. Check missing data after normalization
missing3 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Normalization)) {
   missing3[[i]] <- sum(is.na(TMT_Intensities1_10_Normalization[[i]]))}
missing_tot3 <- set_names(missing3, file_names_wd) #names each file by file_names_wd
missing_tot3 # Total missing for each file


         #4. Check missing data after imputation and normalization
missing4 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation_Normalization)) {
   missing4[[i]] <- sum(is.na(TMT_Intensities1_10_Imputation_Normalization[[i]]))}
missing_tot4 <- set_names(missing4, file_names_wd) #names each file by file_names_wd
missing_tot4 # Total missing for each file



#Missing per spectrum
missing_spec <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing_spec[[i]] <- is.na(TMT_Intensities1_10[[i]])
}
missing_spec #for ALL spectra, prints TRUE/FALSE per chanel

#Missing total over 10 files
missing2 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing2[[i]] <- sum(is.na(TMT_Intensities1_10[[i]]))}
missing_file_tot <- set_names(missing2, file_names_wd) #names each file by file_names_wd
missing_file_tot # Total missing for each file

#Mean missing values over 10 files
missing3 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing3[[i]] <- mean(is.na(TMT_Intensities1_10[[i]]))
}
missing_file_mean <- set_names(missing3, file_names_wd) #names each file by file_names_wd 
missing_file_mean #mean missing values per file

#Missing per row 
missing4 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing4[[i]] <- rowSums(is.na(TMT_Intensities1_10[[i]]))}
missing_row <- set_names(missing4, file_names_wd) #names each file by file_names_wd
missing_row #missing for each row

#Mean missing per row 
missing6 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing6[[i]] <- mean(rowSums(is.na(TMT_Intensities1_10[[i]])))}
missing_row_mean <- set_names(missing6, file_names_wd) #names each file by file_names_wd
missing_row_mean #missing for each row

#Missing total per column/TMT channel
missing5 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing5[[i]] <- colSums(is.na(TMT_Intensities1_10[[i]]))}
missing_col <- set_names(missing5, file_names_wd) #names each file by file_names_wd
missing_col #total missing for each col

#Mean missing per column/TMT channel 
missing7 <- list () #empty list
means <- list()
for (i in seq_along(TMT_Intensities1_10)) {
   missing7[[i]] <- is.na(TMT_Intensities1_10[[i]])
   for (j in seq_along(missing7)) { 
       means[[j]] <- colMeans(missing7[[j]])         
}}
missing_col_mean <- set_names(means, file_names_wd) #names each file by file_names_wd
missing_col_mean #mean missing for each col

#Max TMT intensity
missing8 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing8[[i]] <- max(TMT_Intensities1_10[[i]], na.rm=TRUE)
}
max <- set_names(missing8, file_names_wd) #names each file by file_names_wd
max

#Min TMT intensity
missing9 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing9[[i]] <- min(TMT_Intensities1_10[[i]], na.rm=TRUE)
}
min <- set_names(missing9, file_names_wd) #names each file by file_names_wd
min

#How many zero?
missing10 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing10[[i]] <- colZeros(TMT_Intensities1_10[[i]])
}
zero_col <- set_names(missing10, file_names_wd) #names each file by file_names_wd
zero_col #There are no zero intensities for the TMTs, so NA = 0. 

    #Not really need to run
missing11 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing11[[i]] <- rowZeros(TMT_Intensities1_10[[i]])
}
zero_row <- set_names(missing11, file_names_wd) #names each file by file_names_wd
zero_row


#################
#Making Plots
#################

meansgg <-missing_file_mean%>% as_tibble
barplot(TMT_intensities1)

row.names(TMT_intensities1)
colnames(TMT_intensities1)

length(TMT_intensities1)

data.frame(matrix(unlist(missing_file_mean), nrow=length(missing_file_mean), byrow=TRUE))
x %>% 
as_tibble %>%
var_label(x) 
<- c("Mean")


x <- unlist(missing_file_mean)
x %>% 
as_tibble

var_label(x$value) <- "Mean"

%>%
add_column(File_name=file_names) %>%






png(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per File.png",
width=1000, height=1000)
barplot(unlist(missing_file_mean), main="Mean Missing Values per File", xlab="Mean Missing Values", ylab=file_names, las=1)
dev.off()


file_names <- c("B1S1_f10","B2S4_f10","B3S2_f09","B3S4_f04","B3S4_f06","B5S1_f08","B5S2_f04","B5S2_f07","B5S5_f04","B5S5_f08")
file_names
file_names_wd
means <- c(0.01650336, 0.02377111, 0.01958903, 0.01656064, 0.0144778, 0.0236738, 0.02782404, 0.02605497, 0.02480649, 0.02118762)
df <- tibble(File_name=file_names , Mean=means)

png(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per File.png",
width=600, height=350)
ggplot(df, mapping = aes(x=File_name, y=Mean)) +
    geom_col() +
    labs(x="Mean Missing Values", y="File Name", title="Mean Missing Values per File")
dev.off()



