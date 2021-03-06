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
############################################
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
############################################
#File Names
wd <- setwd("~/Desktop/Read raw file/Data mzML")
getwd() 
(file_names_wd <- list.files(wd)) #The first 10 mzML files of CPTAC
# Loading in multiple mzML files
file_paths <- fs::dir_ls("~/Desktop/Read raw file/Data mzML")
file_paths #The first 10 mzML files paths of CPTAC
file_names_short_typed <- c("B1S1_f10","B2S4_f10","B3S2_f09","B3S4_f04","B3S4_f06","B5S1_f08","B5S2_f04","B5S2_f07","B5S5_f04","B5S5_f08")
file_names <- substring(file_names_wd, 39, 46) #Character 39 untill 46 are unique
file_names #Short file names
file_names_wd #Long file names

#1-4 Plots Exploring missing data in relation to imputation and normalization
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
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Total Missing Values 4 Plots.pdf")
   p1 + p2 + p3 + p4
dev.off()

#1.1-4.1 Mean Values: Plots exploring missing data in relation to imputation and normalization
         #1.1 Check missing data before imputation: mean
missing1_mean <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing1_mean[[i]] <- mean(is.na(TMT_Intensities1_10[[i]]))}
missing_tot1_10_mean <- set_names(missing1_mean, file_names_wd) #names each file by file_names_wd
missing_tot1_10_mean # Mean missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot1_10_mean), nrow=length(missing_tot1_10_mean), byrow=TRUE)) 
Mean_Missing_Values1_10 <- matrix(unlist(missing_tot1_10_mean), nrow=length(missing_tot1_10_mean), byrow=TRUE, ncol=1)
Mean_Missing_Values1_10
df_perc_spec_missing <- tibble(File_name=file_names , Mean_Missing_Values=Mean_Missing_Values1_10)
df_perc_spec_missing
          #2.1 Check missing data after imputation: mean
missing2_mean <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation)) {
   missing2_mean[[i]] <- mean(is.na(TMT_Intensities1_10_Imputation[[i]]))}
missing_tot2_mean <- set_names(missing2_mean, file_names_wd) #names each file by file_names_wd
missing_tot2_mean # Mean missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot2_mean), nrow=length(missing_tot2_mean), byrow=TRUE)) 
Mean_Missing_Values1_10_Imputation <- matrix(unlist(missing_tot2_mean), nrow=length(missing_tot2_mean), byrow=TRUE, ncol=1)
Mean_Missing_Values1_10_Imputation
df_missing_imp_mean <- tibble(File_name=file_names , Mean_Missing_Values=Mean_Missing_Values1_10_Imputation)
df_missing_imp_mean
          #3.1 Check missing data after normalization: mean
missing3_mean <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Normalization)) {
   missing3_mean[[i]] <- mean(is.na(TMT_Intensities1_10_Normalization[[i]]))}
missing_tot3_mean <- set_names(missing3_mean, file_names_wd) #names each file by file_names_wd
missing_tot3_mean # Mean missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot3_mean), nrow=length(missing_tot3_mean), byrow=TRUE)) 
Mean_Missing_Values1_10_Normalization <- matrix(unlist(missing_tot3_mean), nrow=length(missing_tot3_mean), byrow=TRUE, ncol=1)
Mean_Missing_Values1_10_Normalization
df_missing_norm_mean <- tibble(File_name=file_names , Mean_Missing_Values=Mean_Missing_Values1_10_Normalization)
df_missing_norm_mean
         #4.1 Check missing data after imputation and normalization: mean
missing4_mean <- list () #empty list
for (i in seq_along(TMT_Intensities1_10_Imputation_Normalization)) {
   missing4_mean[[i]] <- mean(is.na(TMT_Intensities1_10_Imputation_Normalization[[i]]))}
missing_tot4_mean <- set_names(missing4_mean, file_names_wd) #names each file by file_names_wd
missing_tot4_mean # Mean missing for each file
            #Transforming List Data
data.frame(matrix(unlist(missing_tot4_mean), nrow=length(missing_tot4_mean), byrow=TRUE)) 
Mean_Missing_Values1_10_Imputation_Normalization <- matrix(unlist(missing_tot4_mean), nrow=length(missing_tot4_mean), byrow=TRUE, ncol=1)
Mean_Missing_Values1_10_Imputation_Normalization
df_missing_imp_norm_mean <- tibble(File_name=file_names , Mean_Missing_Values=Mean_Missing_Values1_10_Imputation_Normalization)
df_missing_imp_norm_mean
   #1.1-4.1 Plots
p1_mean <- ggplot(df_perc_spec_missing, mapping = aes(x=File_name, y=Mean_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values: No Imputation", 
      subtitle="Mean missing values before imputation", tag="A") +
   geom_text(aes(label=round(Mean_Missing_Values1_10, digits=3)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p2_mean <- ggplot(df_missing_imp_mean, mapping = aes(x=File_name, y=Mean_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values: Imputation", 
      subtitle="Mean missing values after imputation", tag="B") +
   geom_text(aes(label=round(Mean_Missing_Values1_10_Imputation, digits=3)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p3_mean <- ggplot(df_missing_norm_mean, mapping = aes(x=File_name, y=Mean_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values: Normalization", 
      subtitle="Mean missing values after normalization", tag="C") +
   geom_text(aes(label=round(Mean_Missing_Values1_10_Normalization, digits = 4)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p4_mean <- ggplot(df_missing_imp_norm_mean, mapping = aes(x=File_name, y=Mean_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values: Normalization", 
      subtitle="Mean missing values after imputation and normalization", tag="D") +
   geom_text(aes(label=round(Mean_Missing_Values1_10_Imputation_Normalization, digits = 4)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
   #Putting the 4 plots on 1 page
p1_mean + p2_mean + p3_mean + p4_mean
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values 4 Plots.pdf")
   p1_mean + p2_mean + p3_mean + p4_mean
dev.off()

#5. Perc of Spectra with at least one missing TMT intensity
missing5 <- list() #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing5[[i]] <- sum(rowSums(is.na(TMT_Intensities1_10[[i]])) > 0)/ nrow(TMT_Intensities1_10[[i]])*100
}
missing_row_mean <- set_names(missing5, file_names_wd) #names each file by file_names_wd
missing_row_mean #missing for each  spectrum 
data.frame(matrix(unlist(missing_row_mean), nrow=length(missing_row_mean), byrow=TRUE)) 
Perc_spec_missing <- matrix(unlist(missing_row_mean), nrow=length(missing_row_mean), byrow=TRUE)
Perc_spec_missing
df_perc_spec_missing <- tibble(File_name = file_names , Perc_spec_missing = Perc_spec_missing)
df_perc_spec_missing
  #Plot 5
p5 <- ggplot(df_perc_spec_missing, mapping = aes(x = File_name, y = Perc_spec_missing)) +
   geom_col() +
   labs(x = "File Name", y = "Percentage of Spectra", 
   title = "Percentage of Spectra with Missing Intensities", 
         subtitle = "Per file, the percentage of spectra with at least one missing TMT intensity") +
   geom_text(aes(label = round(Perc_spec_missing, digits = 2)), 
               position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 18), plot.subtitle = element_text(size = 10))
   #Print plot 5
p5
pdf(file = "~/Desktop/Read raw file/TMT outputs/Plots/Percentage of Spectra with Missing Intensities.pdf")
   p5
dev.off()

#6. Mean missing per column/TMT channel 
missing6 <- list () #empty list
means <- list()
for (i in seq_along(TMT_Intensities1_10)) {
   missing6[[i]] <- is.na(TMT_Intensities1_10[[i]])
   for (j in seq_along(missing6)) { 
       means[[j]] <- colMeans(missing6[[j]])
}}
missing_col_mean <- set_names(means, file_names_wd) #names each file by file_names_wd
missing_col_mean #mean missing for each col
data.frame(matrix(unlist(missing_col_mean), nrow=length(missing_col_mean), byrow=TRUE)) 
Mean_Missing_Channel <- data.frame(matrix(unlist(missing_col_mean), nrow=length(missing_col_mean), byrow=TRUE)) 
TMT_Labels <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131")
   #Example of creating new column names. Not needed to create plot 6
cols <- colnames(Mean_Missing_Channel)
df_newnames_Mean_Missing_Channel <- setnames(Mean_Missing_Channel, old = cols, new = TMT_Labels)
tbl_files <- tibble(File_name=file_names)
tbl_tmt <- as_tibble(df_newnames_Mean_Missing_Channel)
tbl_missing_col_mean <- cbind(tbl_files, tbl_tmt)
df_missing_col_mean <- tibble(File_name=file_names , Missing_Channel=tbl_tmt)
df_missing_col_mean

   #Create Data Frame for Plot 6
Labels = rep(TMT_Labels, each=10)
File = rep(file_names, times=10)
Missing_Values = (Mean_Missing_Channel)
df = data.frame(Labels,File,Missing_Values)
df
      #Lenghts have to be equal
length(Labels) #100
length(File) #100
length(Missing_Values) #100
   #Plot 6
p6 <- ggplot(df, mapping = aes(x=Labels, y=Missing_Values, fill=File)) +
   geom_bar(stat="identity", position="dodge") +
   labs(x="TMT Channel", y="Mean Missing Values", title="Mean Missing Values per TMT Channel") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8)) +
   geom_text(aes(label=round(Missing_Values, digits = 4)), 
               position=position_dodge(width=0.9), vjust=0.5, size = 1.55, angle = 90, hjust=-0.1) 
   #Print plot 6
p6
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Mean Missing Values per TMT Channel.pdf")
   p6
dev.off()

#7-8 Plots Max and Min Intensities
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
missing8 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing8[[i]] <- min(TMT_Intensities1_10[[i]], na.rm=TRUE)
}
min <- set_names(missing8, file_names_wd) #names each file by file_names_wd
min
data.frame(matrix(unlist(min), nrow=length(min), byrow=TRUE)) 
Min_Values1_10 <- matrix(unlist(min), nrow=length(min), byrow=TRUE, ncol=1)
Min_Values1_10
df_min <- tibble(File_name=file_names , Min_Values=Min_Values1_10)
df_min
   #7-8 Plots
p7 <- ggplot(df_max, mapping = aes(x=File_name, y=Max_Values)) +
   geom_col() +
   labs(x="File Name", y="Maximun Intensities", title="Maximun TMT Intensities", 
      subtitle="Maximun measuered TMT intensities", tag="A") +
   geom_text(aes(label=round(Max_Values1_10, digits = 0)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p8 <- ggplot(df_min, mapping = aes(x=File_name, y=Min_Values)) +
   geom_col() +
   labs(x="File Name", y="Minimum Intensiteis", title="Minumim TMT Intensities", 
      subtitle="Minimum measuered TMT intensities", tag="B") +
   geom_text(aes(label=round(Min_Values1_10, digits = 4)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
   #Print 8-9 plots together
p7 + p8
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Plots Max and Min TMT Intensities.pdf")
  p7 + p8
dev.off()

#9-10 Plot difference in missing values after imputation
  #total
diff <- mapply('-', missing_tot1_10, missing_tot4, SIMPLIFY = FALSE)
missing_diff <- set_names(diff, file_names_wd)
missing_diff
data.frame(matrix(unlist(missing_diff), nrow=length(missing_diff), byrow=TRUE)) 
Difference <- matrix(unlist(missing_diff), nrow=length(missing_diff), byrow=TRUE, ncol=1)
Difference
df_diff <- tibble(File_name=file_names , Decrease_Missing_Values=Difference)
df_diff
  #perc lowered missing values
diff_perc <- mapply('/', missing_diff, missing_tot1_10, SIMPLIFY = FALSE)
diff_perc2 <- mapply('*', diff_perc, 100, SIMPLIFY = FALSE)
missing_diff_perc <- set_names(diff_perc2, file_names_wd)
missing_diff_perc
data.frame(matrix(unlist(missing_diff_perc), nrow=length(missing_diff_perc), byrow=TRUE)) 
Difference_Perc <- matrix(unlist(missing_diff_perc), nrow=length(missing_diff_perc), byrow=TRUE, ncol=1)
Difference_Perc
df_diff_perc <- tibble(File_name=file_names , Perc_Decrease_Missing_Values=Difference_Perc)
df_diff_perc
   #Plots 9-10
p9 <- ggplot(df_diff, mapping = aes(x=File_name, y=Decrease_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Decrease missing values", title="Decrease Missing Values", 
      subtitle="Decrease in missing values after imputation and normalization", tag="A") +
   geom_text(aes(label=Difference), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 6))
p10 <- ggplot(df_diff_perc, mapping = aes(x=File_name, y=Perc_Decrease_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Perc decrease missing values", title="Perc Decrease Missing Values", 
      subtitle="Perc decrease in missing values after imputation and normalization", tag="B") +
   geom_text(aes(label=round(Difference_Perc, digits = 0)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.5) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 6))
  #Print 9-10 plots together
p9 + p10
pdf(file="~/Desktop/Read raw file/TMT outputs/Plots/Plots Differences Missing Values.pdf")
   p9 + p10
dev.off()

######################################
#Check Missing Data from PSMs
######################################

#11. Perc of Spectra with at least one missing TMT intensity
   #TMTs matched to PSMs
PSM_TMT <- readRDS("~/Desktop/mzTab/Stored files/PSMs linked to TMT intensities")
   #Selecting the TMT columns
PSM_TMT[[1]]
TMT_psm <- list()
for (i in seq_along(PSM_TMT)) {
   TMT_psm[[i]] <- PSM_TMT[[i]] %>% select(`126`:`131`)
}
TMT_psm #Tibblesof only TMT intensities
file_names_short_typed <- c("B1S2_f10", "B1S3_f05", "B1S3_f10", "B1S4_f02", "B1S4_f06", "B2S4_f09")
missing11 <- list() #empty list
for (i in seq_along(TMT_psm)) {
   missing11[[i]] <- sum(rowSums(is.na(TMT_psm[[i]])) > 0)/ nrow(TMT_psm[[i]])*100
}
PSM_missing_row <- set_names(missing11, file_names_short_typed) #names each file by file_names_short_typed
PSM_missing_row #missing for each  spectrum 
   #unlist
data.frame(matrix(unlist(PSM_missing_row), nrow=length(PSM_missing_row), byrow=TRUE)) 
PSM_Perc_spec_missing <- matrix(unlist(PSM_missing_row), nrow=length(PSM_missing_row), byrow=TRUE)
PSM_Perc_spec_missing
PSM_df_perc_spec_missing <- tibble(File_name = file_names_short_typed , Perc_spec_missing = PSM_Perc_spec_missing)
PSM_df_perc_spec_missing
  #Plot 11
p11 <- ggplot(PSM_df_perc_spec_missing, mapping = aes(x = File_name, y = PSM_Perc_spec_missing)) +
   geom_col() +
   labs(x = "File Name", y = "Percentage of PSM Spectra", 
   title = "Percentage of Spectra with Missing Intensities", 
         subtitle = "Afeter matchin TMTs to PSMS. Per file, the percentage of spectra with at least one missing TMT intensity") +
   geom_text(aes(label = round(PSM_Perc_spec_missing, digits = 2)), 
               position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 18), plot.subtitle = element_text(size = 10))
   #Print plot 11
p11
pdf(file = "~/Desktop/Read raw file/TMT outputs/Plots/Percentage_of_Spectra_Missing_Intensities_after_PSM_Matching.pdf")
   p11
dev.off()

         #1.1 Check missing data before imputation: mean
missing12_mean <- list () #empty list
for (i in seq_along(TMT_psm)) {
   missing12_mean[[i]] <- mean(is.na(TMT_psm[[i]]))}
PSM_mm12 <- set_names(missing12_mean, file_names_short_typed) #names each file by file_names_wd
PSM_mm12 # Mean missing for each file
            #Transforming List Data
data.frame(matrix(unlist(PSM_mm12), nrow=length(PSM_mm12), byrow=TRUE)) 
PSM_MM12 <- matrix(unlist(PSM_mm12), nrow=length(PSM_mm12), byrow=TRUE, ncol=1)
PSM_MM12
df_perc_spec_missing <- tibble(File_name=file_names_short_typed , Mean_Missing_Values=PSM_MM12)
df_perc_spec_missing

p12 <- ggplot(df_perc_spec_missing, mapping = aes(x=File_name, y=Mean_Missing_Values)) +
   geom_col() +
   labs(x="File Name", y="Mean Missing Values", title="Mean Missing Values after PSM TMT Matching") +
   geom_text(aes(label=round(PSM_MM12, digits=3)), 
               position=position_dodge(width=0.9), vjust=-0.25, size = 1.75) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
            plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
p12
pdf(file = "~/Desktop/Read raw file/TMT outputs/Plots/Mean_Missing_Values_after_PSM_TMT_Matching.pdf")
   p12
dev.off()












###################################
#Loops for exploring missing values
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
   #Missing per spectrum 
missing4 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing4[[i]] <- rowSums(is.na(TMT_Intensities1_10[[i]]))}
missing_row <- set_names(missing4, file_names_wd) #names each file by file_names_wd
missing_row #missing for each spectrum
   #Mean missing per spectrum 
missing5 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing5[[i]] <- mean(rowSums(is.na(TMT_Intensities1_10[[i]])))}
missing_row_mean <- set_names(missing5, file_names_wd) #names each file by file_names_wd
missing_row_mean #missing for each spectrum 
   #Missing total per column/TMT channel
missing5 <- list () #empty list
for (i in seq_along(TMT_Intensities1_10)) {
   missing5[[i]] <- colSums(is.na(TMT_Intensities1_10[[i]]))}
missing_col <- set_names(missing5, file_names_wd) #names each file by file_names_wd
missing_col #total missing for each col
   #Mean missing per column/TMT channel 
missing6 <- list () #empty list
means <- list()
for (i in seq_along(TMT_Intensities1_10)) {
   missing6[[i]] <- is.na(TMT_Intensities1_10[[i]])
   for (j in seq_along(missing6)) { 
       means[[j]] <- colMeans(missing6[[j]])         
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
missing8 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing8[[i]] <- min(TMT_Intensities1_10[[i]], na.rm=TRUE)
}
min <- set_names(missing8, file_names_wd) #names each file by file_names_wd
min
   #How many zero?
missing10 <- list()
for (i in seq_along(TMT_Intensities1_10)) {
    missing10[[i]] <- colZeros(TMT_Intensities1_10[[i]])
}
zero_col <- set_names(missing10, file_names_wd) #names each file by file_names_wd
zero_col #There are no zero intensities for the TMTs, so NA = 0. 