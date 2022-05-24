#PSM
PSM_TMT_ALL <- readRDS("/Users/jensvandeperre/Desktop/Outputs/PSM_TMT_linked/ALL_PSM_TMT_Linked")
    #Don't forget to devide by REF
normalized <- list()
for (i in 1:264) {
    normalized[[i]] <- PSM_TMT_ALL[[i]] %>%
        mutate("126" = "126"/"131")
}

#Proteins PIA
(file_paths_PIA <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(file_paths_PIA[[i]], header = FALSE, sep = ",")
}

PepSeq_ProAcc <- list()
for (i in 1:264) {
PepSeq_ProAcc[[i]] <- PIA[[i]] %>%
          as_tibble() %>%
          filter(str_detect(V1, "PSMSET")) %>%
          mutate(sequence_no_mod = V2) %>%
          mutate(Accessions = V3) %>%
          dplyr::select(sequence_no_mod, Accessions) %>%
          slice(-1) %>%
          as_tibble()
}

Des <- list()
for (i in 1:264) {
Des[[i]] <- PIA[[i]] %>%
          as_tibble() %>%
          filter(str_detect(V1, "PROTEIN")) %>%
          mutate(Accessions = V2) %>%
          mutate(ClusterID = V8) %>%
          mutate(Description = V9) %>%
          dplyr::select(Accessions, ClusterID, Description) %>%
          slice(-1) %>%
          as_tibble()
}

PIA_PSM_TMT <- list()
for (i in 1:264) {
PIA_PSM_TMT[[i]] <- merge(PepSeq_ProAcc[[i]], Des[[i]], by = "Accessions") %>%
            as_tibble() %>%
            rename("Protein.Group.Accessions" = Accessions, "Protein.Descriptions" = Description) %>%
            merge(PSM_TMT_ALL[[i]], by = "sequence_no_mod") %>%
            distinct()
}







#PTMs
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
PTM_PSM_TMT <- list()
for (i in 1:264) {
    PTM_PSM_TMT[[i]] <- merge(PSM_TMT[[i]], ptm_index[[i]], by = "index")
}
view(PTM_PSM_TMT[[1]])
length(PTM_PSM_TMT)
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








