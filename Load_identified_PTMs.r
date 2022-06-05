#Load identified PTMs
(ptm_filepaths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PTM_identification_tol_10"))
PTM <- list()
for (i in 1:264) {
  PTM[[i]] <- read.csv(ptm_filepaths[[i]], sep = ",", header = TRUE)
}
view(PTM[[1]])
nrow(PTM[[1]])
length(PTM)
#Create index for matching to PSM_TMT
ptm_index <- list()
for (i in 1:264) {
  ptm_index[[i]] <- PTM[[i]] %>%
  as_tibble() %>%
  mutate(index = trimws(str_remove_all(PSM_ID, "controllerType=0 "))) %>%
  mutate(index = trimws(str_remove_all(index, "controllerNumber=1 "))) %>%
  mutate(index = trimws(str_remove_all(index, "scan=")))
}
view(ptm_index[[1]])

#merge to PSM_TMT
PSM_TMT <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
PTM_PSM_TMT <- list()
for (i in 1:264) {
  PTM_PSM_TMT[[i]] <- merge(PSM_TMT[[i]], ptm_index[[i]], by = "index")
}
view(PTM_PSM_TMT[[1]])
length(PTM_PSM_TMT)

#Look for Cancer PTMs




#Most prevalant PTMs




#How many types of PTMs



#Most modified protein
(pia_filepaths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(pia_filepaths[[i]], sep = ",", header = TRUE)
}
view(PIA[[1]])
nrow(PIA[[1]])
length(PIA)




#Make PTM_PSM_TMT ready for merge to PIA output
#Select the needed columns
                    PTM_PSM_TMT_input <- list()
                    for (i in 1:264) {
                      PTM_PSM_TMT_input[[i]] <- PSM_TMT_ALL[[i]] %>% 
                      as_tibble() %>%
                      select(index, sequence, sequence_no_mod, "126":"131", charge) %>%
                      mutate(file_name = file_names_short[[i]]) 
                    }
                    view(PTM_PSM_TMT_input[[1]])
                    view(PTM_PSM_TMT_input[[2]])
                    view(PTM_PSM_TMT_input[[69]])
                    length(PTM_PSM_TMT_input)
