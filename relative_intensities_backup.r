      batch_B1S1_f01_f12 <- list() #
      batch_B1S2_f01_f12 <- list() #
      batch_B1S3_f01_f12 <- list() #
      batch_B1S4_f01_f12 <- list() #
      batch_B2S1_f01_f12 <- list() #
      batch_B2S2_f01_f12 <- list() #
      batch_B2S3_f01_f12 <- list() #
      batch_B2S4_f01_f12 <- list() #
      batch_B3S1_f01_f12 <- list() #
      batch_B3S2_f01_f12 <- list() #
      batch_B3S3_f01_f12 <- list() #
      batch_B3S4_f01_f12 <- list() #
      batch_B4S1_f01_f12 <- list() #
      batch_B4S2_f01_f12 <- list() #
      batch_B4S3_f01_f12 <- list() #
      batch_B4S4_f01_f12 <- list() #
      batch_B5S1_f01_f12 <- list() #
      batch_B5S2_f01_f12 <- list() #
      batch_B5S3_f01_f12 <- list() #
      batch_B5S4_f01_f12 <- list() #
      batch_B5S5_f01_f12 <- list() #
      batch_B5S6_f01_f12 <- list() #



dat_B1S1_f01_f12 <- bind_rows(batch_B1S1_f01_f12) %>%
  group_by(Accessions) %>% 
  summarise(NAT_126_B1S1_f01_f12 = sum(NAT_126_B1S1_f01_f12), 
  NAT_127N_B1S1_f01_f12 = sum(NAT_127N_B1S1_f01_f12), 
  TUMOR_127C_B1S1_f01_f12 = sum(TUMOR_127C_B1S1_f01_f12),
  TUMOR_128N_B1S1_f01_f12 =sum(TUMOR_128N_B1S1_f01_f12),
  TUMOR_128C_B1S1_f01_f12 =sum(TUMOR_128C_B1S1_f01_f12),
  TUMOR_129N_B1S1_f01_f12 =sum(TUMOR_129N_B1S1_f01_f12),
  TUMOR_129C_B1S1_f01_f12 =sum(TUMOR_129C_B1S1_f01_f12),
  TUMOR_130N_B1S1_f01_f12 =sum(TUMOR_130N_B1S1_f01_f12),
  TUMOR_130C_B1S1_f01_f12 =sum(TUMOR_130C_B1S1_f01_f12),
  REF_131_B1S1_f01_f12 =sum(REF_131_B1S1_f01_f12),
  n.peptides =sum(n.peptides),
  n.spectra =sum(n.spectra)
  ) %>%
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

dat_B1S2_f01_f12 <- bind_rows(batch_B1S2_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B1S2_f01_f12 = sum(NAT_126_B1S2_f01_f12),
        TUMOR_127N_B1S2_f01_f12 = sum(TUMOR_127N_B1S2_f01_f12),
        NAT_127C_B1S2_f01_f12 = sum(NAT_127C_B1S2_f01_f12),
        TUMOR_128N_B1S2_f01_f12 = sum(TUMOR_128N_B1S2_f01_f12),
        NAT_128C_B1S2_f01_f12 = sum(NAT_128C_B1S2_f01_f12),
        TUMOR_129N_B1S2_f01_f12 = sum(TUMOR_129N_B1S2_f01_f12),
        NAT_129C_B1S2_f01_f12 = sum(NAT_129C_B1S2_f01_f12),
        TUMOR_130N_B1S2_f01_f12 = sum(TUMOR_130N_B1S2_f01_f12),
        TUMOR_130C_B1S2_f01_f12 = sum(TUMOR_130C_B1S2_f01_f12),
        REF_131_B1S2_f01_f12 = sum(REF_131_B1S2_f01_f12),
        n.peptides =sum(n.peptides),
        n.spectra =sum(n.spectra)
    ) %>%
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

dat_B1S3_f01_f12 <- bind_rows(batch_B1S3_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B1S3_f01_f12 = sum(TUMOR_126_B1S3_f01_f12),
        NAT_127N_B1S3_f01_f12 = sum(NAT_127N_B1S3_f01_f12),
        NAT_128N_B1S3_f01_f12 = sum(NAT_128N_B1S3_f01_f12),
        NAT_128C_B1S3_f01_f12 = sum(NAT_128C_B1S3_f01_f12),
        NAT_129N_B1S3_f01_f12 = sum(NAT_129N_B1S3_f01_f12),
        NAT_129C_B1S3_f01_f12 = sum(NAT_129C_B1S3_f01_f12),
        TUMOR_130N_B1S3_f01_f12 = sum(TUMOR_130N_B1S3_f01_f12),
        TUMOR_130C_B1S3_f01_f12 = sum(TUMOR_130C_B1S3_f01_f12),
        REF_131_B1S3_f01_f12 = sum(REF_131_B1S3_f01_f12),
        n.peptides =sum(n.peptides),
        n.spectra =sum(n.spectra)
    ) %>%
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

dat_B1S4_f01_f12 <- bind_rows(batch_B1S4_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_161_B1S4_f01_f12 = sum(NAT_161_B1S4_f01_f12),
        TUMOR_127N_B1S4_f01_f12 = sum(TUMOR_127N_B1S4_f01_f12),
        NAT_127C_B1S4_f01_f12 = sum(NAT_127C_B1S4_f01_f12),
        TUMOR_128N_B1S4_f01_f12 = sum(TUMOR_128N_B1S4_f01_f12),
        TUMOR_128C_B1S4_f01_f12 = sum(TUMOR_128C_B1S4_f01_f12),
        TUMOR_129N_B1S4_f01_f12 = sum(TUMOR_129N_B1S4_f01_f12),
        TUMOR_129C_B1S4_f01_f12 = sum(TUMOR_129C_B1S4_f01_f12),
        NAT_130N_B1S4_f01_f12 = sum(NAT_130N_B1S4_f01_f12),
        TUMOR_130C_B1S4_f01_f12 = sum(TUMOR_130C_B1S4_f01_f12),
        REF_131__B1S4_f01_f12 =sum(REF_131__B1S4_f01_f12),
        n.peptides =sum(n.peptides),
        n.spectra =sum(n.spectra)
    ) %>%
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

dat_B2S1_f01_f12 <- bind_rows(batch_B2S1_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B2S1_f01_f12 = sum(TUMOR_126_B2S1_f01_f12)
        TUMOR_127N_B2S1_f01_f12 = sum(TUMOR_127N_B2S1_f01_f12)
        TUMOR_127C_B2S1_f01_f12 = sum(TUMOR_127C_B2S1_f01_f12)
        NAT_128N_B2S1_f01_f12 = sum(NAT_128N_B2S1_f01_f12)
        NAT_128C_B2S1_f01_f12 = sum(NAT_128C_B2S1_f01_f12)
        TUMOR_129N_B2S1_f01_f12 = sum(TUMOR_129N_B2S1_f01_f12)
        NAT_129C_B2S1_f01_f12 = sum(NAT_129C_B2S1_f01_f12)
        TUMOR_130N_B2S1_f01_f12 = sum(TUMOR_130N_B2S1_f01_f12)
        TUMOR_130C_B2S1_f01_f12 = sum(TUMOR_130C_B2S1_f01_f12)
        REF_131_B2S1_f01_f12 = sum(REF_131_B2S1_f01_f12)
        n.peptides =sum(n.peptides),
                n.spectra =sum(n.spectra)
    ) %>%
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

dat_B2S2_f01_f12 <- bind_rows(batch_B2S2_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B2S2_f01_f12 = sum(TUMOR_126_B2S2_f01_f12),
        TUMOR_127N_B2S2_f01_f12 = sum(TUMOR_127N_B2S2_f01_f12),
        TUMOR_127C_B2S2_f01_f12 = sum(TUMOR_127C_B2S2_f01_f12),
        NAT_128N_B2S2_f01_f12 = sum(NAT_128N_B2S2_f01_f12),
        NAT_128C_B2S2_f01_f12 = sum(NAT_128C_B2S2_f01_f12),
        TUMOR_129N_B2S2_f01_f12 = sum(TUMOR_129N_B2S2_f01_f12),
        NAT_129C_B2S2_f01_f12 = sum(NAT_129C_B2S2_f01_f12),
        TUMOR_130N_B2S2_f01_f12 = sum(TUMOR_130N_B2S2_f01_f12),
        TUMOR_130C_B2S2_f01_f12 = sum(TUMOR_130C_B2S2_f01_f12),
        REF_131_B2S2_f01_f12 =sum(REF_131_B2S2_f01_f12),
            n.peptides =sum(n.peptides),
            n.spectra =sum(n.spectra)
    ) %>%
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

dat_B2S3_f01_f12 <- bind_rows(batch_B2S3_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B2S3_f01_f12 = sum(TUMOR_126_B2S3_f01_f12),
        NAT_127N_B2S3_f01_f12 = sum(NAT_127N_B2S3_f01_f12),
        TUMOR_127C_B2S3_f01_f12 = sum(TUMOR_127C_B2S3_f01_f12),
        NAT_128N_B2S3_f01_f12 = sum(NAT_128N_B2S3_f01_f12),
        NAT_128C_B2S3_f01_f12 = sum(NAT_128C_B2S3_f01_f12),
        NAT_129N_B2S3_f01_f12 = sum(NAT_129N_B2S3_f01_f12),
        NAT_129C_B2S3_f01_f12 = sum(NAT_129C_B2S3_f01_f12),
        TUMOR_130N_B2S3_f01_f12 = sum(TUMOR_130N_B2S3_f01_f12),
        NAT_130C_B2S3_f01_f12 = sum(NAT_130C_B2S3_f01_f12),
        REF_131_B2S3_f01_f12 = sum(REF_131_B2S3_f01_f12),
                    n.peptides =sum(n.peptides),
                    n.spectra =sum(n.spectra)
    ) %>%
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

dat_B2S4_f01_f12 <- bind_rows(batch_B2S4_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B2S4_f01_f12 = (TUMOR_126_B2S4_f01_f12),
        TUMOR_127N_B2S4_f01_f12 = (TUMOR_127N_B2S4_f01_f12),
        NAT_127C_B2S4_f01_f12 = (NAT_127C_B2S4_f01_f12),
        NAT_128N_B2S4_f01_f12 = (NAT_128N_B2S4_f01_f12),
        TUMOR_128C_B2S4_f01_f12 = (TUMOR_128C_B2S4_f01_f12),
        TUMOR_129N_B2S4_f01_f12 = (TUMOR_129N_B2S4_f01_f12),
        NAT_129C_B2S4_f01_f12 = (NAT_129C_B2S4_f01_f12),
        TUMOR_130N_B2S4_f01_f12 = (TUMOR_130N_B2S4_f01_f12),
        NAT_130C_B2S4_f01_f12 = (NAT_130C_B2S4_f01_f12),
        REF_131_B2S4_f01_f12 = sum(REF_131_B2S4_f01_f12),
            n.peptides =sum(n.peptides),
            n.spectra =sum(n.spectra)
    ) %>%
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

dat_B3S1_f01_f12 <- bind_rows(batch_B3S1_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B3S1_f01_f12 = sum(TUMOR_126_B3S1_f01_f12),
        NAT_127N_B3S1_f01_f12 = sum(NAT_127N_B3S1_f01_f12),
        NAT_127C_B3S1_f01_f12 = sum(NAT_127C_B3S1_f01_f12),
        TUMOR_128N_B3S1_f01_f12 = sum(TUMOR_128N_B3S1_f01_f12),
        TUMOR_128C_B3S1_f01_f12 = sum(TUMOR_128C_B3S1_f01_f12),
        TUMOR_129N_B3S1_f01_f12 = sum(TUMOR_129N_B3S1_f01_f12),
        TUMOR_129C_B3S1_f01_f12 = sum(TUMOR_129C_B3S1_f01_f12),
        TUMOR_130N_B3S1_f01_f12 = sum(TUMOR_130N_B3S1_f01_f12),
        TUMOR_130C_B3S1_f01_f12 = sum(TUMOR_130C_B3S1_f01_f12),
        REF_131_B3S1_f01_f12 = sum(REF_131_B3S1_f01_f12)
           n.peptides =sum(n.peptides),
            n.spectra =sum(n.spectra)
    ) %>%
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

dat_B3S2_f01_f12 <- bind_rows(batch_B3S2_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B3S2_f01_f12 = sum(NAT_126_B3S2_f01_f12),
        NAT_127N_B3S2_f01_f12 = sum(NAT_127N_B3S2_f01_f12),
        NAT_127C_B3S2_f01_f12 = sum(NAT_127C_B3S2_f01_f12),
        NAT_128N_B3S2_f01_f12 = sum(NAT_128N_B3S2_f01_f12),
        TUMOR_128C_B3S2_f01_f12 = sum(TUMOR_128C_B3S2_f01_f12),
        NAT_129N_B3S2_f01_f12 = sum(NAT_129N_B3S2_f01_f12),
        NAT_129C_B3S2_f01_f12 = sum(NAT_129C_B3S2_f01_f12),
        TUMOR_130N_B3S2_f01_f12 = sum(TUMOR_130N_B3S2_f01_f12),
        TUMOR_130C_B3S2_f01_f12 = sum(TUMOR_130C_B3S2_f01_f12),
        REF_131_B3S2_f01_f12 = sum(REF_131_B3S2_f01_f12),
                n.peptides =sum(n.peptides),
                    n.spectra =sum(n.spectra)
    ) %>%
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

dat_B3S3_f01_f12 <- bind_rows(batch_B3S3_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B3S3_f01_f12 = sum(TUMOR_126_B3S3_f01_f12),
        NAT_127N_B3S3_f01_f12 = sum(NAT_127N_B3S3_f01_f12),
        TUMOR_127C_B3S3_f01_f12 = sum(TUMOR_127C_B3S3_f01_f12),
        NAT_128N_B3S3_f01_f12 = sum(NAT_128N_B3S3_f01_f12),
        TUMOR_128C_B3S3_f01_f12 = sum(TUMOR_128C_B3S3_f01_f12),
        NAT_129N_B3S3_f01_f12 = sum(NAT_129N_B3S3_f01_f12),
        TUMOR_129C_B3S3_f01_f12 = sum(TUMOR_129C_B3S3_f01_f12),
        TUMOR_130N_B3S3_f01_f12 = sum(TUMOR_130N_B3S3_f01_f12),
        NAT_130C_B3S3_f01_f12 = sum(NAT_130C_B3S3_f01_f12),
                n.peptides =sum(n.peptides),
                    n.spectra =sum(n.spectra)
    ) %>%
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

dat_B3S4_f01_f12 <- bind_rows(batch_B3S4_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B3S4_f01_f12 = sum(NAT_126_B3S4_f01_f12),
        NAT_127N_B3S4_f01_f12 = sum(NAT_127N_B3S4_f01_f12),
        TUMOR_127C_B3S4_f01_f12 = sum(TUMOR_127C_B3S4_f01_f12),
        NAT_128N_B3S4_f01_f12 = sum(NAT_128N_B3S4_f01_f12),
        NAT_128C_B3S4_f01_f12 = sum(NAT_128C_B3S4_f01_f12),
        TUMOR_129N_B3S4_f01_f12 = sum(TUMOR_129N_B3S4_f01_f12),
        TUMOR_129C_B3S4_f01_f12 = sum(TUMOR_129C_B3S4_f01_f12),
        NAT_130N_B3S4_f01_f12 = sum(NAT_130N_B3S4_f01_f12),
        NAT_130C_B3S4_f01_f12 = sum(NAT_130C_B3S4_f01_f12),
        REF_131_B3S4_f01_f12 = sum(REF_131_B3S4_f01_f12),
                n.peptides =sum(n.peptides),
                    n.spectra =sum(n.spectra)
    ) %>%
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

dat_B4S1_f01_f12 <- bind_rows(batch_B4S1_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B4S1_f01_f12 = sum(NAT_126_B4S1_f01_f12),
        TUMOR_127N_B4S1_f01_f12 = sum(TUMOR_127N_B4S1_f01_f12),
        NAT_127C_B4S1_f01_f12 = sum(NAT_127C_B4S1_f01_f12),
        NAT_128N_B4S1_f01_f12 = sum(NAT_128N_B4S1_f01_f12),
        NAT_128C_B4S1_f01_f12 = sum(NAT_128C_B4S1_f01_f12),
        NAT_129N_B4S1_f01_f12 = sum(NAT_129N_B4S1_f01_f12),
        TUMOR_129C_B4S1_f01_f12 = sum(TUMOR_129C_B4S1_f01_f12),
        NAT_130N_B4S1_f01_f12 = sum(NAT_130N_B4S1_f01_f12),
        NAT_130C_B4S1_f01_f12 = sum(NAT_130C_B4S1_f01_f12),
        REF_131_B4S1_f01_f12 = sum(REF_131_B4S1_f01_f12),
                n.peptides =sum(n.peptides),
                    n.spectra =sum(n.spectra)
    ) %>%
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

dat_B4S2_f01_f12 <- bind_rows(batch_B4S2_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B4S2_f01_f12 = sum(NAT_126_B4S2_f01_f12),
        TUMOR_127N_B4S2_f01_f12 = sum(TUMOR_127N_B4S2_f01_f12),
        TUMOR_127C_B4S2_f01_f12 = sum(TUMOR_127C_B4S2_f01_f12),
        NAT_128N_B4S2_f01_f12 = sum(NAT_128N_B4S2_f01_f12),
        NAT_128C_B4S2_f01_f12 = sum(NAT_128C_B4S2_f01_f1),
        TUMOR_129N_B4S2_f01_f12 = sum(TUMOR_129N_B4S2_f01_f12),
        NAT_129C_B4S2_f01_f12 = sum(NAT_129C_B4S2_f01_f12),
        TUMOR_130N_B4S2_f01_f12 = sum(TUMOR_130N_B4S2_f01_f12),
        NAT_130C_B4S2_f01_f12 = sum(NAT_130C_B4S2_f01_f12),
        REF_131_B4S2_f01_f12 = sum(REF_131_B4S2_f01_f12)
                        n.peptides =sum(n.peptides),
                            n.spectra =sum(n.spectra)
    ) %>%
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

dat_B4S3_f01_f12 <- bind_rows(batch_B4S3_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B4S3_f01_f12 = sum(TUMOR_126_B4S3_f01_f12),
        NAT_127N_B4S3_f01_f12 = sum(NAT_127N_B4S3_f01_f12),
        TUMOR_127C_B4S3_f01_f12 = sum(TUMOR_127C_B4S3_f01_f12),
        NAT_128N_B4S3_f01_f12 = sum(NAT_128N_B4S3_f01_f12),
        NAT_128C_B4S3_f01_f12 = sum(NAT_128C_B4S3_f01_f12),
        TUMOR_129N_B4S3_f01_f12 = sum(TUMOR_129N_B4S3_f01_f12),
        NAT_129C_B4S3_f01_f12 = sum(NAT_129C_B4S3_f01_f12),
        TUMOR_130N_B4S3_f01_f12 = sum(TUMOR_130N_B4S3_f01_f12),
        NAT_130C_B4S3_f01_f12 = sum(NAT_130C_B4S3_f01_f12),
        REF_131_B4S3_f01_f12 = sum(REF_131_B4S3_f01_f12),
                                n.peptides =sum(n.peptides),
                                    n.spectra =sum(n.spectra)
    ) %>%
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

dat_B4S4_f01_f12 <- bind_rows(batch_B4S4_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
TUMOR_126_B4S4_f01_f12 = sum(TUMOR_126_B4S4_f01_f12),
TUMOR_127N_B4S4_f01_f12 = sum(TUMOR_127N_B4S4_f01_f12),
NAT_127C_B4S4_f01_f12 = sum(NAT_127C_B4S4_f01_f12),
NAT_128N_B4S4_f01_f12 = sum(NAT_128N_B4S4_f01_f12),
TUMOR_128C_B4S4_f01_f12 = sum(TUMOR_128C_B4S4_f01_f12),
NAT_129N_B4S4_f01_f12 = sum(NAT_129N_B4S4_f01_f12),
TUMOR_129C_B4S4_f01_f12 = sum(TUMOR_129C_B4S4_f01_f12),
NAT_130N_B4S4_f01_f12 = sum(NAT_130N_B4S4_f01_f12),
NAT_130C_B4S4_f01_f12 = sum(NAT_130C_B4S4_f01_f12),
REF_131_B4S4_f01_f12 = sum(REF_131_B4S4_f01_f12),
            n.peptides =sum(n.peptides),
            n.spectra =sum(n.spectra)

    ) %>%

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

dat_B5S1_f01_f12 <- bind_rows(batch_B5S1_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B5S1_f01_f12 = sum(NAT_126_B5S1_f01_f12),
        TUMOR_127N_B5S1_f01_f12 = sum(TUMOR_127N_B5S1_f01_f12),
        NAT_127C_B5S1_f01_f12 = sum(NAT_127C_B5S1_f01_f12),
        NAT_128N_B5S1_f01_f12 = sum(NAT_128N_B5S1_f01_f12),
        NAT_128C_B5S1_f01_f12 = sum(NAT_128C_B5S1_f01_f12),
        NAT_129N_B5S1_f01_f12 = sum(NAT_129N_B5S1_f01_f12),
        TUMOR_129C_B5S1_f01_f12 = sum(TUMOR_129C_B5S1_f01_f12),
        TUMOR_130N_B5S1_f01_f12 = sum(TUMOR_130N_B5S1_f01_f12),
        TUMOR_130C_B5S1_f01_f12 = sum(TUMOR_130C_B5S1_f01_f12),
        REF_131_B5S1_f01_f12 = sum(REF_131_B5S1_f01_f12),
                n.peptides =sum(n.peptides),
                n.spectra =sum(n.spectra)
    ) %>%
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

dat_B5S2_f01_f12 <- bind_rows(batch_B5S2_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        NAT_126_B5S2_f01_f12 = sum(NAT_126_B5S2_f01_f12),
        NAT_127N_B5S2_f01_f12 = sum(NAT_127N_B5S2_f01_f12),
        NAT_127C_B5S2_f01_f12 = sum(NAT_127C_B5S2_f01_f12),
        TUMOR_128N_B5S2_f01_f12 = sum(TUMOR_128N_B5S2_f01_f12),
        TUMOR_128C_B5S2_f01_f12 = sum(TUMOR_128C_B5S2_f01_f12),
        NAT_129N_B5S2_f01_f12 = sum(NAT_129N_B5S2_f01_f12),
        TUMOR_129C_B5S2_f01_f12 = sum(TUMOR_129C_B5S2_f01_f12),
        TUMOR_130N_B5S2_f01_f12 = sum(TUMOR_130N_B5S2_f01_f12),
        NAT_130C_B5S2_f01_f12 = sum(NAT_130C_B5S2_f01_f12),
        REF_131_B5S2_f01_f12 = sum(REF_131_B5S2_f01_f12),
                n.peptides =sum(n.peptides),
                n.spectra =sum(n.spectra)
    ) %>%
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

dat_B5S3_f01_f12 <- bind_rows(batch_B5S3_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B5S3_f01_f12 = sum(TUMOR_126_B5S3_f01_f12),
        TUMOR_127N_B5S3_f01_f12 = sum(TUMOR_127N_B5S3_f01_f12),
        NAT_127C_B5S3_f01_f12 = sum(NAT_127C_B5S3_f01_f12),
        NAT_128N_B5S3_f01_f12 = sum(NAT_128N_B5S3_f01_f12),
        TUMOR_128C_B5S3_f01_f12 = sum(TUMOR_128C_B5S3_f01_f12),
        TUMOR_129N_B5S3_f01_f12 = sum(TUMOR_129N_B5S3_f01_f12),
        TUMOR_129C_B5S3_f01_f12 = sum(TUMOR_129C_B5S3_f01_f12),
        NAT_130N_B5S3_f01_f12 = sum(NAT_130N_B5S3_f01_f12),
        TUMOR_130C_B5S3_f01_f12 = sum(TUMOR_130C_B5S3_f01_f12),
        REF_131_B5S3_f01_f12 = sum(REF_131_B5S3_f01_f12),
            n.peptides =sum(n.peptides),
            n.spectra =sum(n.spectra)
    ) %>%
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

dat_B5S4_f01_f12 <- bind_rows(batch_B5S4_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
        TUMOR_126_B5S4_f01_f12 = sum(TUMOR_126_B5S4_f01_f12),
        TUMOR_127N_B5S4_f01_f12 = sum(TUMOR_127N_B5S4_f01_f12),
        TUMOR_127C_B5S4_f01_f12 = sum(TUMOR_127C_B5S4_f01_f12),
        NAT_128N_B5S4_f01_f12 = sum(NAT_128N_B5S4_f01_f12),
        TUMOR_128C_B5S4_f01_f12 = sum(TUMOR_128C_B5S4_f01_f12),
        NAT_129N_B5S4_f01_f12 = sum(NAT_129N_B5S4_f01_f12),
        NAT_129C_B5S4_f01_f12 = sum(NAT_129C_B5S4_f01_f12),
        NAT_130N_B5S4_f01_f12 = sum(NAT_130N_B5S4_f01_f12),
        NAT_130C_B5S4_f01_f12 = sum(NAT_130C_B5S4_f01_f12),
        REF_131_B5S4_f01_f12 = sum(REF_131_B5S4_f01_f12),
                n.peptides =sum(n.peptides),
                n.spectra =sum(n.spectra)
    ) %>%
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

dat_B5S5_f01_f12 <- bind_rows(batch_B5S5_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
TUMOR_126_B5S5_f01_f12 = sum(TUMOR_126_B5S5_f01_f12),
TUMOR_127N_B5S5_f01_f12 = sum(TUMOR_127N_B5S5_f01_f12),
TUMOR_127C_B5S5_f01_f12 = sum(TUMOR_127C_B5S5_f01_f12),
NAT_128N_B5S5_f01_f12 = sum(NAT_128N_B5S5_f01_f12),
NAT_128C_B5S5_f01_f12 = sum(NAT_128C_B5S5_f01_f12),
TUMOR_129N_B5S5_f01_f12 = sum(TUMOR_129N_B5S5_f01_f12),
TUMOR_129C_B5S5_f01_f12 = sum(TUMOR_129C_B5S5_f01_f12),
NAT_130N_B5S5_f01_f12 = sum(NAT_130N_B5S5_f01_f12),
NAT_130C_B5S5_f01_f12 = sum(NAT_130C_B5S5_f01_f12),
REF_131_B5S5_f01_f12 = sum(REF_131_B5S5_f01_f12),
        n.peptides =sum(n.peptides),
        n.spectra =sum(n.spectra)
    ) %>%
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

dat_B5S6_f01_f12 <- bind_rows(batch_B5S6_f01_f12) %>%
    group_by(Accessions) %>% 
    summarise(
TUMOR_126_B5S6_f01_f12 = sum(TUMOR_126_B5S6_f01_f12),
NAT_127N_B5S6_f01_f12 = sum(NAT_127N_B5S6_f01_f12),
NAT_127C_B5S6_f01_f12 = sum(NAT_127C_B5S6_f01_f12),
NAT_128N_B5S6_f01_f12 = sum(NAT_128N_B5S6_f01_f12),
NAT_128C_B5S6_f01_f12 = sum(NAT_128C_B5S6_f01_f12),
NAT_129N_B5S6_f01_f12 = sum(NAT_129N_B5S6_f01_f12),
NAT_129C_B5S6_f01_f12 = sum(NAT_129C_B5S6_f01_f12),
NAT_130N_B5S6_f01_f12 = sum(NAT_130N_B5S6_f01_f12),
REF_130C_B5S6_f01_f12 = sum(REF_130C_B5S6_f01_f12),
REF_131_B5S6_f01_f12 = sum(REF_131_B5S6_f01_f12),
        n.peptides =sum(n.peptides),
        n.spectra =sum(n.spectra)
    ) %>%
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

