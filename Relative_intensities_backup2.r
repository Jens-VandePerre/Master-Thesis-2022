for (i in 1:12) {
    batch_B1S1_f01_f12[[i]] <- proteins_B1S1_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B1S2_f01_f12[[i]] <- proteins_B1S2_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B1S3_f01_f12[[i]] <- proteins_B1S3_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B1S4_f01_f12[[i]] <- proteins_B1S4_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B2S1_f01_f12[[i]] <- proteins_B2S1_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B2S2_f01_f12[[i]] <- proteins_B2S2_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B2S3_f01_f12[[i]] <- proteins_B2S3_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B2S4_f01_f12[[i]] <- proteins_B2S4_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B3S1_f01_f12[[i]] <- proteins_B3S1_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B3S2_f01_f12[[i]] <- proteins_B3S2_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B3S3_f01_f12[[i]] <- proteins_B3S3_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B3S4_f01_f12[[i]] <- proteins_B3S4_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B4S1_f01_f12[[i]] <- proteins_B4S1_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B4S2_f01_f12[[i]] <- proteins_B4S2_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B4S3_f01_f12[[i]] <- proteins_B4S3_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B4S4_f01_f12[[i]] <- proteins_B4S4_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B5S1_f01_f12[[i]] <- proteins_B5S1_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B5S2_f01_f12[[i]] <- proteins_B5S2_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B5S3_f01_f12[[i]] <- proteins_B5S3_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B5S4_f01_f12[[i]] <- proteins_B5S4_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B5S5_f01_f12[[i]] <- proteins_B5S5_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}
for (i in 1:12) {
    batch_B5S6_f01_f12[[i]] <- proteins_B5S6_f01_f12[[i]] %>% 
    rownames_to_column(var = "Accessions")
}




dat_B1S1_f01_f12 <- bind_rows(peptides_B1S1_f01_f12) %>%
  mutate (NAT_126_B1S1_f01_f12 = NAT_126_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (NAT_127N_B1S1_f01_f12 = NAT_127N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_127C_B1S1_f01_f12= TUMOR_127C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_128N_B1S1_f01_f12 = TUMOR_128N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_128C_B1S1_f01_f12 = TUMOR_128C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_129N_B1S1_f01_f12 = TUMOR_129N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_129C_B1S1_f01_f12 = TUMOR_129C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_130N_B1S1_f01_f12 = TUMOR_130N_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  mutate (TUMOR_130C_B1S1_f01_f12 = TUMOR_130C_B1S1_f01_f12/REF_131_B1S1_f01_f12) %>%
  select(-REF_131_B1S1_f01_f12)
nrow(dat_B1S1_f01_f12) 
view(dat_B1S1_f01_f12) 

dat_B1S2_f01_f12 <- bind_rows(peptides_B1S2_f01_f12) %>%
  mutate(NAT_126_B1S2_f01_f12 = NAT_126_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_127N_B1S2_f01_f12 = TUMOR_127N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(NAT_127C_B1S2_f01_f12 = NAT_127C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_128N_B1S2_f01_f12 = TUMOR_128N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(NAT_128C_B1S2_f01_f12 = NAT_128C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_129N_B1S2_f01_f12 = TUMOR_129N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(NAT_129C_B1S2_f01_f12 = NAT_129C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_130N_B1S2_f01_f12 = TUMOR_130N_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  mutate(TUMOR_130C_B1S2_f01_f12 = TUMOR_130C_B1S2_f01_f12/REF_131_B1S2_f01_f12) %>%
  select(-REF_131_B1S2_f01_f12) 

dat_B1S3_f01_f12 <- bind_rows(peptides_B1S3_f01_f12) %>%
    mutate(TUMOR_126_B1S3_f01_f12 = TUMOR_126_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_127N_B1S3_f01_f12 = NAT_127N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(TUMOR_127C_B1S3_f01_f12 = TUMOR_127C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_128N_B1S3_f01_f12 = NAT_128N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_128C_B1S3_f01_f12 = NAT_128C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_129N_B1S3_f01_f12 = NAT_129N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(NAT_129C_B1S3_f01_f12 = NAT_129C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(TUMOR_130N_B1S3_f01_f12 = TUMOR_130N_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  mutate(TUMOR_130C_B1S3_f01_f12 = TUMOR_130C_B1S3_f01_f12/REF_131_B1S3_f01_f12) %>%
  select(-REF_131_B1S3_f01_f12) 

dat_B1S4_f01_f12 <- bind_rows(peptides_B1S4_f01_f12) %>%
    mutate(NAT_161_B1S4_f01_f12 = NAT_161_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_127N_B1S4_f01_f12 = TUMOR_127N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(NAT_127C_B1S4_f01_f12 = NAT_127C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_128N_B1S4_f01_f12 = TUMOR_128N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_128C_B1S4_f01_f12 = TUMOR_128C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_129N_B1S4_f01_f12 = TUMOR_129N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_129C_B1S4_f01_f12 = TUMOR_129C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(NAT_130N_B1S4_f01_f12 = NAT_130N_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  mutate(TUMOR_130C_B1S4_f01_f12 = TUMOR_130C_B1S4_f01_f12/REF_131__B1S4_f01_f12) %>%
  select(-REF_131__B1S4_f01_f12) 

dat_B2S1_f01_f12 <- bind_rows(peptides_B2S1_f01_f12) %>%
    mutate(TUMOR_126_B2S1_f01_f12 = TUMOR_126_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_127N_B2S1_f01_f12 = TUMOR_127N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_127C_B2S1_f01_f12 = TUMOR_127C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(NAT_128N_B2S1_f01_f12 = NAT_128N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(NAT_128C_B2S1_f01_f12 = NAT_128C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_129N_B2S1_f01_f12 = TUMOR_129N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(NAT_129C_B2S1_f01_f12 = NAT_129C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_130N_B2S1_f01_f12 = TUMOR_130N_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  mutate(TUMOR_130C_B2S1_f01_f12 = TUMOR_130C_B2S1_f01_f12/REF_131_B2S1_f01_f12) %>%
  select(-REF_131_B2S1_f01_f12) 

dat_B2S2_f01_f12 <- bind_rows(peptides_B2S2_f01_f12) %>%
    mutate(TUMOR_126_B2S2_f01_f12 = TUMOR_126_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_127N_B2S2_f01_f12 = TUMOR_127N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_127C_B2S2_f01_f12 = TUMOR_127C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(NAT_128N_B2S2_f01_f12 = NAT_128N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(NAT_128C_B2S2_f01_f12 = NAT_128C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_129N_B2S2_f01_f12 = TUMOR_129N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(NAT_129C_B2S2_f01_f12 = NAT_129C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_130N_B2S2_f01_f12 = TUMOR_130N_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  mutate(TUMOR_130C_B2S2_f01_f12 = TUMOR_130C_B2S2_f01_f12/REF_131_B2S2_f01_f12) %>%
  select(-REF_131_B2S2_f01_f12) 

dat_B2S3_f01_f12 <- bind_rows(peptides_B2S3_f01_f12) %>%
    mutate(TUMOR_126_B2S3_f01_f12 = TUMOR_126_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_127N_B2S3_f01_f12 = NAT_127N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(TUMOR_127C_B2S3_f01_f12 = TUMOR_127C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_128N_B2S3_f01_f12 = NAT_128N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_128C_B2S3_f01_f12 = NAT_128C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_129N_B2S3_f01_f12 = NAT_129N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_129C_B2S3_f01_f12 = NAT_129C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(TUMOR_130N_B2S3_f01_f12 = TUMOR_130N_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  mutate(NAT_130C_B2S3_f01_f12 = NAT_130C_B2S3_f01_f12/REF_131_B2S3_f01_f12) %>%
  select(-REF_131_B2S3_f01_f12) 

dat_B2S4_f01_f12 <- bind_rows(peptides_B2S4_f01_f12) %>%
    mutate(TUMOR_126_B2S4_f01_f12 = TUMOR_126_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_127N_B2S4_f01_f12 = TUMOR_127N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_127C_B2S4_f01_f12 = NAT_127C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_128N_B2S4_f01_f12 = NAT_128N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_128C_B2S4_f01_f12 = TUMOR_128C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_129N_B2S4_f01_f12 = TUMOR_129N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_129C_B2S4_f01_f12 = NAT_129C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(TUMOR_130N_B2S4_f01_f12 = TUMOR_130N_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  mutate(NAT_130C_B2S4_f01_f12 = NAT_130C_B2S4_f01_f12/REF_131_B2S4_f01_f12) %>%
  select(-REF_131_B2S4_f01_f12) 

dat_B3S1_f01_f12 <- bind_rows(peptides_B3S1_f01_f12) %>%
    mutate(TUMOR_126_B3S1_f01_f12 = TUMOR_126_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(NAT_127N_B3S1_f01_f12 = NAT_127N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(NAT_127C_B3S1_f01_f12 = NAT_127C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_128N_B3S1_f01_f12 = TUMOR_128N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_128C_B3S1_f01_f12 = TUMOR_128C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_129N_B3S1_f01_f12 = TUMOR_129N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_129C_B3S1_f01_f12 = TUMOR_129C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_130N_B3S1_f01_f12 = TUMOR_130N_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  mutate(TUMOR_130C_B3S1_f01_f12 = TUMOR_130C_B3S1_f01_f12/REF_131_B3S1_f01_f12) %>%
  select(-REF_131_B3S1_f01_f12) 

dat_B3S2_f01_f12 <- bind_rows(peptides_B3S2_f01_f12) %>%
    mutate(NAT_126_B3S2_f01_f12 = NAT_126_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_127N_B3S2_f01_f12 = NAT_127N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_127C_B3S2_f01_f12 = NAT_127C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_128N_B3S2_f01_f12 = NAT_128N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(TUMOR_128C_B3S2_f01_f12 = TUMOR_128C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_129N_B3S2_f01_f12 = NAT_129N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(NAT_129C_B3S2_f01_f12 = NAT_129C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(TUMOR_130N_B3S2_f01_f12 = TUMOR_130N_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  mutate(TUMOR_130C_B3S2_f01_f12 = TUMOR_130C_B3S2_f01_f12/REF_131_B3S2_f01_f12) %>%
  select(-REF_131_B3S2_f01_f12) 

dat_B3S3_f01_f12 <- bind_rows(peptides_B3S3_f01_f12) %>%
    mutate(TUMOR_126_B3S3_f01_f12 = TUMOR_126_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_127N_B3S3_f01_f12 = NAT_127N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_127C_B3S3_f01_f12 = TUMOR_127C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_128N_B3S3_f01_f12 = NAT_128N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_128C_B3S3_f01_f12 = TUMOR_128C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_129N_B3S3_f01_f12 = NAT_129N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_129C_B3S3_f01_f12 = TUMOR_129C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(TUMOR_130N_B3S3_f01_f12 = TUMOR_130N_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  mutate(NAT_130C_B3S3_f01_f12 = NAT_130C_B3S3_f01_f12/REF_131_B3S3_f01_f12) %>%
  select(-REF_131_B3S3_f01_f12) 

dat_B3S4_f01_f12 <- bind_rows(peptides_B3S4_f01_f12) %>%
    mutate(NAT_126_B3S4_f01_f12 = NAT_126_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_127N_B3S4_f01_f12 = NAT_127N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(TUMOR_127C_B3S4_f01_f12 = TUMOR_127C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_128N_B3S4_f01_f12 = NAT_128N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_128C_B3S4_f01_f12 = NAT_128C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(TUMOR_129N_B3S4_f01_f12 = TUMOR_129N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(TUMOR_129C_B3S4_f01_f12 = TUMOR_129C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_130N_B3S4_f01_f12 = NAT_130N_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  mutate(NAT_130C_B3S4_f01_f12 = NAT_130C_B3S4_f01_f12/REF_131_B3S4_f01_f12) %>%
  select(-REF_131_B3S4_f01_f12) 

dat_B4S1_f01_f12 <- bind_rows(peptides_B4S1_f01_f12) %>%
    mutate(NAT_126_B4S1_f01_f12 = NAT_126_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(TUMOR_127N_B4S1_f01_f12 = TUMOR_127N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_127C_B4S1_f01_f12 = NAT_127C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_128N_B4S1_f01_f12 = NAT_128N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_128C_B4S1_f01_f12 = NAT_128C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_129N_B4S1_f01_f12 = NAT_129N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(TUMOR_129C_B4S1_f01_f12 = TUMOR_129C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_130N_B4S1_f01_f12 = NAT_130N_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  mutate(NAT_130C_B4S1_f01_f12 = NAT_130C_B4S1_f01_f12/REF_131_B4S1_f01_f12) %>%
  select(-REF_131_B4S1_f01_f12) 

dat_B4S2_f01_f12 <- bind_rows(peptides_B4S2_f01_f12) %>%
    mutate(NAT_126_B4S2_f01_f12 = NAT_126_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_127N_B4S2_f01_f12 = TUMOR_127N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_127C_B4S2_f01_f12 = TUMOR_127C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_128N_B4S2_f01_f12 = NAT_128N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_128C_B4S2_f01_f12 = NAT_128C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_129N_B4S2_f01_f12 = TUMOR_129N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_129C_B4S2_f01_f12 = NAT_129C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(TUMOR_130N_B4S2_f01_f12 = TUMOR_130N_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  mutate(NAT_130C_B4S2_f01_f12 = NAT_130C_B4S2_f01_f12/REF_131_B4S2_f01_f12) %>%
  select(-REF_131_B4S2_f01_f12) 

dat_B4S3_f01_f12 <- bind_rows(peptides_B4S3_f01_f12) %>%
    mutate(TUMOR_126_B4S3_f01_f12 = TUMOR_126_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_127N_B4S3_f01_f12 = NAT_127N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(TUMOR_127C_B4S3_f01_f12 = TUMOR_127C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_128N_B4S3_f01_f12 = NAT_128N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_128C_B4S3_f01_f12 = NAT_128C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(TUMOR_129N_B4S3_f01_f12 = TUMOR_129N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_129C_B4S3_f01_f12 = NAT_129C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(TUMOR_130N_B4S3_f01_f12 = TUMOR_130N_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  mutate(NAT_130C_B4S3_f01_f12 = NAT_130C_B4S3_f01_f12/REF_131_B4S3_f01_f12) %>%
  select(-REF_131_B4S3_f01_f12) 

dat_B4S4_f01_f12 <- bind_rows(peptides_B4S4_f01_f12) %>%
    mutate(TUMOR_126_B4S4_f01_f12 = TUMOR_126_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(TUMOR_127N_B4S4_f01_f12 = TUMOR_127N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_127C_B4S4_f01_f12 = NAT_127C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_128N_B4S4_f01_f12 = NAT_128N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(TUMOR_128C_B4S4_f01_f12 = TUMOR_128C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_129N_B4S4_f01_f12 = NAT_129N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(TUMOR_129C_B4S4_f01_f12 = TUMOR_129C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_130N_B4S4_f01_f12 = NAT_130N_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  mutate(NAT_130C_B4S4_f01_f12 = NAT_130C_B4S4_f01_f12/REF_131_B4S4_f01_f12) %>%
  select(-REF_131_B4S4_f01_f12) 

dat_B5S1_f01_f12 <- bind_rows(peptides_B5S1_f01_f12) %>%
    mutate(NAT_126_B5S1_f01_f12 = NAT_126_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_127N_B5S1_f01_f12 = TUMOR_127N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_127C_B5S1_f01_f12 = NAT_127C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_128N_B5S1_f01_f12 = NAT_128N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_128C_B5S1_f01_f12 = NAT_128C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(NAT_129N_B5S1_f01_f12 = NAT_129N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_129C_B5S1_f01_f12 = TUMOR_129C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_130N_B5S1_f01_f12 = TUMOR_130N_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  mutate(TUMOR_130C_B5S1_f01_f12 = TUMOR_130C_B5S1_f01_f12/REF_131_B5S1_f01_f12) %>%
  select(-REF_131_B5S1_f01_f12) 

dat_B5S2_f01_f12 <- bind_rows(peptides_B5S2_f01_f12) %>%
    mutate(NAT_126_B5S2_f01_f12 = NAT_126_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_127N_B5S2_f01_f12 = NAT_127N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_127C_B5S2_f01_f12 = NAT_127C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_128N_B5S2_f01_f12 = TUMOR_128N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_128C_B5S2_f01_f12 = TUMOR_128C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_129N_B5S2_f01_f12 = NAT_129N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_129C_B5S2_f01_f12 = TUMOR_129C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(TUMOR_130N_B5S2_f01_f12 = TUMOR_130N_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  mutate(NAT_130C_B5S2_f01_f12 = NAT_130C_B5S2_f01_f12/REF_131_B5S2_f01_f12) %>%
  select(-REF_131_B5S2_f01_f12) 

dat_B5S3_f01_f12 <- bind_rows(peptides_B5S3_f01_f12) %>%
    mutate(TUMOR_126_B5S3_f01_f12 = TUMOR_126_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_127N_B5S3_f01_f12 = TUMOR_127N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(NAT_127C_B5S3_f01_f12 = NAT_127C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(NAT_128N_B5S3_f01_f12 = NAT_128N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_128C_B5S3_f01_f12 = TUMOR_128C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_129N_B5S3_f01_f12 = TUMOR_129N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_129C_B5S3_f01_f12 = TUMOR_129C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(NAT_130N_B5S3_f01_f12 = NAT_130N_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  mutate(TUMOR_130C_B5S3_f01_f12 = TUMOR_130C_B5S3_f01_f12/REF_131_B5S3_f01_f12) %>%
  select(-REF_131_B5S3_f01_f12) 

dat_B5S4_f01_f12 <- bind_rows(peptides_B5S4_f01_f12) %>%
    mutate(TUMOR_126_B5S4_f01_f12 = TUMOR_126_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(TUMOR_127N_B5S4_f01_f12 = TUMOR_127N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(TUMOR_127C_B5S4_f01_f12 = TUMOR_127C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_128N_B5S4_f01_f12 = NAT_128N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(TUMOR_128C_B5S4_f01_f12 = TUMOR_128C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_129N_B5S4_f01_f12 = NAT_129N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_129C_B5S4_f01_f12 = NAT_129C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_130N_B5S4_f01_f12 = NAT_130N_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  mutate(NAT_130C_B5S4_f01_f12 = NAT_130C_B5S4_f01_f12/REF_131_B5S4_f01_f12) %>%
  select(-REF_131_B5S4_f01_f12) 

dat_B5S5_f01_f12 <- bind_rows(peptides_B5S5_f01_f12) %>%
    mutate(TUMOR_126_B5S5_f01_f12 = TUMOR_126_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_127N_B5S5_f01_f12 = TUMOR_127N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_127C_B5S5_f01_f12 = TUMOR_127C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_128N_B5S5_f01_f12 = NAT_128N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_128C_B5S5_f01_f12 = NAT_128C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_129N_B5S5_f01_f12 = TUMOR_129N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(TUMOR_129C_B5S5_f01_f12 = TUMOR_129C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_130N_B5S5_f01_f12 = NAT_130N_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  mutate(NAT_130C_B5S5_f01_f12 = NAT_130C_B5S5_f01_f12/REF_131_B5S5_f01_f12) %>%
  select(-REF_131_B5S5_f01_f12) 

dat_B5S6_f01_f12 <- bind_rows(peptides_B5S6_f01_f12) %>%
    mutate(TUMOR_126_B5S6_f01_f12 = TUMOR_126_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_127N_B5S6_f01_f12 = NAT_127N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_127C_B5S6_f01_f12 = NAT_127C_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_128N_B5S6_f01_f12 = NAT_128N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_128C_B5S6_f01_f12 = NAT_128C_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_129N_B5S6_f01_f12 = NAT_129N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_129C_B5S6_f01_f12 = NAT_129C_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  mutate(NAT_130N_B5S6_f01_f12 = NAT_130N_B5S6_f01_f12/REF_131_B5S6_f01_f12) %>%
  select(-REF_130C_B5S6_f01_f12 ) %>%
  select(-REF_131_B5S6_f01_f12 ) 
