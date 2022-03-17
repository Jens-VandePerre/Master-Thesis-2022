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


########################################################################################
#Load in outputs directly, not running loops
############################################
#1. Output intensities for ALL SPECTRA
      #not removing NAs
      #no imputation
      #no normalisation
TMT_Intensities1_10 <- readRDS(file = "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10")
#2. Output intensities for ALL spectra
    #impute: method="MLE"
    #no normalisation
TMT_Intensities1_10_Imputation <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation")
#3. Loop that outputs intensities for ALL spectra
    #no imputation
    #normalise: method="center.median"
TMT_Intensities1_10_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Normalization")
#4. Output intensities for ALL spectra
    #impute: method="MLE"
    #normalise: method="center.median"
TMT_Intensities1_10_Imputation_Normalization <- readRDS(file= "~/Desktop/Read raw file/TMT outputs/Combined Files/TMT_Intensities1-10_Imputation+Normalization")
########################################################################################

wd <- setwd("~/Desktop/Read raw file/Data mzML")
getwd() 
(file_names_wd <- list.files(wd)) #The first 10 mzML files of CPTAC
# Loading in multiple mzML files
file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML")
file_paths #The first 10 mzML files paths of CPTAC


file_names <- c("B1S1_f10","B2S4_f10","B3S2_f09","B3S4_f04","B3S4_f06","B5S1_f08","B5S2_f04","B5S2_f07","B5S5_f04","B5S5_f08")
file_names #Short file names
file_names_wd #Long file names


         #1. Check missing data before imputation
missing1 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing1[[i]] <- sum(is.na(TMT_Intensities1_10[[i]]))}
missing_tot1_10 <- set_names(missing1, file_names_wd) #names each file by file_names_wd
missing_tot1_10 # Total missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot1_10), nrow=length(missing_tot1_10), byrow=TRUE)) 
Total_Missing_Values1_10 <- matrix(unlist(missing_tot1_10), nrow=length(missing_tot1_10), byrow=TRUE, ncol=1)
Total_Missing_Values1_10
df_missing <- tibble(File_name=file_names , Total_Missing_Values=Total_Missing_Values1_10)
df_missing

          #2. Check missing data after imputation
missing2 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation)) {
   missing2[[i]] <- sum(is.na(TMT_Intensities1_10_Imputation[[i]]))}
missing_tot2 <- set_names(missing2, file_names_wd) #names each file by file_names_wd
missing_tot2 # Total missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot2), nrow=length(missing_tot2), byrow=TRUE)) 
Total_Missing_Values1_10_Imputation <- matrix(unlist(missing_tot2), nrow=length(missing_tot2), byrow=TRUE, ncol=1)
Total_Missing_Values1_10_Imputation
df_missing_imp <- tibble(File_name=file_names , Total_Missing_Values=Total_Missing_Values1_10_Imputation)
df_missing_imp

          #3. Check missing data after normalization
missing3 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Normalization)) {
   missing3[[i]] <- sum(is.na(TMT_Intensities1_10_Normalization[[i]]))}
missing_tot3 <- set_names(missing3, file_names_wd) #names each file by file_names_wd
missing_tot3 # Total missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot3), nrow=length(missing_tot3), byrow=TRUE)) 
Total_Missing_Values1_10_Normalization <- matrix(unlist(missing_tot3), nrow=length(missing_tot3), byrow=TRUE, ncol=1)
Total_Missing_Values1_10_Normalization
df_missing_norm <- tibble(File_name=file_names , Total_Missing_Values=Total_Missing_Values1_10_Normalization)
df_missing_norm

         #4. Check missing data after imputation and normalization
missing4 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation_Normalization)) {
   missing4[[i]] <- sum(is.na(TMT_Intensities1_10_Imputation_Normalization[[i]]))}
missing_tot4 <- set_names(missing4, file_names_wd) #names each file by file_names_wd
missing_tot4 # Total missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot4), nrow=length(missing_tot4), byrow=TRUE)) 
Total_Missing_Values1_10_Imputation_Normalization <- matrix(unlist(missing_tot4), nrow=length(missing_tot4), byrow=TRUE, ncol=1)
Total_Missing_Values1_10_Imputation_Normalization
df_missing_imp_norm <- tibble(File_name=file_names , Total_Missing_Values=Total_Missing_Values1_10_Imputation_Normalization)
df_missing_imp_norm


#1-4 Plots
p1 <- ggplot(df_missing, mapping = aes(x=File_name, y=Total_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Total Missing Values", title="Total Missing Values: No Imputation", 
      subtitle="Total missing values before imputation", tag="A") +
   geom_text(aes(label=Total_Missing_Values1_10), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p2 <- ggplot(df_missing_imp, mapping = aes(x=File_name, y=Total_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Total Missing Values", title="Total Missing Values: Imputation", 
      subtitle="Total missing values after imputation", tag="B") +
   geom_text(aes(label=Total_Missing_Values1_10_Imputation), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p3 <- ggplot(df_missing_norm, mapping = aes(x=File_name, y=Total_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Total Missing Values", title="Total Missing Values: Normalization", 
      subtitle="Total missing values after normalization", tag="C") +
   geom_text(aes(label=Total_Missing_Values1_10_Normalization), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p4 <- ggplot(df_missing_imp_norm, mapping = aes(x=File_name, y=Total_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Total Missing Values", title="Total Missing Values: Normalization", 
      subtitle="Total missing values after imputation and normalization", tag="D") +
   geom_text(aes(label=Total_Missing_Values1_10_Imputation_Normalization), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
   #Putting the 4 plots on 1 page
p1 + p2 + p3 + p4
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Total Missing Values Plots.pdf")
   p1 + p2 + p3 + p4
dev.off()

#Mean missing values over 10 files
missing3 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing3[[i]] <- mean(is.na(TMT_Intensities1_10[[i]]))
}
missing_file_mean <- set_names(missing3, file_names_wd) #names each file by file_names_wd 
missing_file_mean #mean missing values per file
data.frame(matrix(unlist(missing_file_mean), nrow=length(missing_file_mean), byrow=TRUE)) 
Mean_Missing <- matrix(unlist(missing_file_mean), nrow=length(missing_file_mean), byrow=TRUE)
Mean_Missing
df_missing_mean <- tibble(File_name=file_names , Missing_Mean=Mean_Missing)
df_missing_mean
  #Plot 3 
p3 <- ggplot(df_missing_mean, mapping = aes(x=File_name, y=Missing_Mean)) +
   geom_col(group=TMT_Labels) +
   labs(x="File Name", y="Mean Missing Intensiteis", title="Mean Missing Values per File") +
   geom_text(aes(label=round(Mean_Missing, digits = 3)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), axis.text.y = element_text(size = 7)
            plot.title = element_text(size = 15))
   #Print 8-9 plots together
p3
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per File.pdf")
   p3
dev.off()

#Mean missing per row 
missing6 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing6[[i]] <- mean(rowSums(is.na(TMT_Intensities1_10[[i]])))}
missing_row_mean <- set_names(missing6, file_names_wd) #names each file by file_names_wd
missing_row_mean #missing for each row

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
data.frame(matrix(unlist(missing_col_mean), nrow=length(missing_col_mean), byrow=TRUE)) 
Mean_Missing_Channel <- matrix(unlist(missing_col_mean), nrow=length(missing_col_mean), byrow=TRUE, ncol=10)
Mean_Missing_Channel
df_missing_col_mean <- tibble(File_name=file_names , Missing_Channel=Mean_Missing_Channel)
df_missing_col_mean
   #Plot 7 NOT WORKING
TMT_Labels <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131")
p7 <- ggplot(df_missing_col_mean, mapping = aes(x=File_name, y=Mean_Missing_Channel)) +
   geom_col(group=TMT_Labels) +
   labs(x="File Name", y="Minimum Intensiteis", title="Mean Missing Values per TMT Channel") +
   geom_text(aes(label=TMT_Labels), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
   #Print 8-9 plots together
p7
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per TMT Channel.pdf")
   p7
dev.off()



#Plots Max and Min Intensities 
   #Max TMT intensity
missing8 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing8[[i]] <- max(TMT_Intensities1_10[[i]], na.rm=TRUE)
}
max <- set_names(missing8, file_names_wd) #names each file by file_names_wd
max
data.frame(matrix(unlist(max), nrow=length(max), byrow=TRUE)) 
Max_Values1_10 <- matrix(unlist(max), nrow=length(max), byrow=TRUE, ncol=1)
Max_Values1_10
df_max <- tibble(File_name=file_names , Max_Values=Max_Values1_10)
df_max
   #Min TMT intensity
missing9 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing9[[i]] <- min(TMT_Intensities1_10[[i]], na.rm=TRUE)
}
min <- set_names(missing9, file_names_wd) #names each file by file_names_wd
min
data.frame(matrix(unlist(min), nrow=length(min), byrow=TRUE)) 
Min_Values1_10 <- matrix(unlist(min), nrow=length(min), byrow=TRUE, ncol=1)
Min_Values1_10
df_min <- tibble(File_name=file_names , Min_Values=Min_Values1_10)
df_min
   #8-9 Plots
p8 <- ggplot(df_max, mapping = aes(x=File_name, y=Max_Values)) +
   geom_col() +
   labs(x="File Name", y="Maximun Intensities", title="Maximun TMT Intensities", 
      subtitle="Maximun measuered TMT intensities", tag="A") +
   geom_text(aes(label=round(Max_Values1_10, digits = 0)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p9 <- ggplot(df_min, mapping = aes(x=File_name, y=Min_Values)) +
   geom_col() +
   labs(x="File Name", y="Minimum Intensiteis", title="Minumim TMT Intensities", 
      subtitle="Minimum measuered TMT intensities", tag="B") +
   geom_text(aes(label=round(Min_Values1_10, digits = 3)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
   #Print 8-9 plots together
p8 + p9
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Plots Max and Min TMT Intensities.pdf")
   p8 + p9
dev.off()









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











means <- c(0.01650336, 0.02377111, 0.01958903, 0.01656064, 0.0144778, 0.0236738, 0.02782404, 0.02605497, 0.02480649, 0.02118762)
df1 <- tibble(File_name=file_names , Mean=means)
df1
png(file="~/Desktop/Read raw file/TMT outputs/Plots/Non Efficient Plot.png",
width=1000, height=750)
ggplot(df1, mapping = aes(x=File_name, y=Mean)) +
    geom_col() +
    labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values", subtitle="Mean missing values for the 10 first CPTAC files")
dev.off()

#Efficient Plot Making
data.frame(matrix(unlist(missing_file_mean), nrow=length(missing_file_mean), byrow=TRUE)) 
mean <- matrix(unlist(missing_file_mean), nrow=length(missing_file_mean), byrow=TRUE, ncol=1)
mean
df <- tibble(File_name=file_names_wd , Mean=mean)
df
ggplot(df, mapping = aes(x=File_name, y=Mean)) +
    geom_col() +
    labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values", 
      subtitle="Mean missing values for the 10 first CPTAC files")
   #png
png(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per File.png",
width=1000, height=750)
ggplot(df, mapping = aes(x=File_name, y=Mean)) +
    geom_col() +
    labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values", 
      subtitle="Mean missing values for the 10 first CPTAC files") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
   #pdf
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per File.pdf",
)
ggplot(df, mapping = aes(x=File_name, y=Mean)) +
    geom_col() +
    labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values", 
      subtitle="Mean missing values for the 10 first CPTAC files")
dev.off()



