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

wd <- setwd("/Users/jensvandeperre/Desktop/Inputs/PSM_TMT_mass_diff/Mass_tolerance_10")
getwd() 
list.files(wd)
#Load PSMs 
#Load matching mzMLs
    #Automate filename extraction
(file_name_long <- list.files(wd))
(file_paths <- fs::dir_ls("/Users/jensvandeperre/Desktop/Inputs/PSM_TMT_mass_diff/Mass_tolerance_10"))
(file_names_short <- substring(file_name_long, 54, 62))

#Load in all PSMs
  #PSM files, with calulated mass differences at mass_tolerance = 10
PSM_TMT <- list()
for (i in 1:264) {
  PSM_TMT[[i]] <- read.csv(file_paths[[i]])
}
view(PSM_TMT[[1]])
nrow(PSM_TMT[[1]])
  #Select the needed columns
PSM_TMT_input <- list()
for (i in 1:264){
  PSM_TMT_input[[i]] <- PSM_TMT[[i]] %>% 
  select(index, "X126":"X131", sequence, sequence_no_mod, charge) %>%
  rename("126" = "X126") %>%
  rename("127N" = "X127N") %>%
  rename("127C" = "X127C") %>%
  rename("128N" = "X128N") %>%
  rename("128C" = "X128C") %>%
  rename("129N" = "X129N") %>%
  rename("129C" = "X129C") %>%
  rename("130N" = "X130N") %>%
  rename("130C" = "X130C") %>%
  rename("131" = "X131")
}
view(PSM_TMT_input[[1]])

#Load protein info from PIA output
(file_paths_PIA <- fs::dir_ls("/Users/jensvandeperre/Desktop/Outputs/PIA_analysis"))
PIA <- list()
for (i in 1:264) {
  PIA[[i]] <- read.csv(file_paths_PIA[[i]], header = FALSE, sep = ",")
}
view(PIA[[1]])
nrow(PIA[[1]])

#Merging: Proteins, Peptides and TMTs
  #Select peptide seq and proteins from PIA
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
view(PepSeq_ProAcc[[1]])
nrow(PepSeq_ProAcc[[1]])
length(PepSeq_ProAcc)
str(PepSeq_ProAcc[[1]])

  #Select proteins and genes from PIA
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
str(Des[[1]])
view(Des[[1]])
nrow(Des[[1]])
length(Des)

  #Merge: proteins, genes, peptide seq/PSM and TMTs
data <- list()
for (i in 1:264) {
data[[i]] <- merge(PepSeq_ProAcc[[i]], Des[[i]], by = "Accessions") %>%
            as_tibble() %>%
            rename("Protein.Group.Accessions" = Accessions, "Protein.Descriptions" = Description) %>%
            merge(PSM_TMT_input[[i]], by = "sequence_no_mod") %>%
            distinct()
}
view(data[[1]])
nrow(data[[1]])

#Start creating input file
mydat <- data %>%
          as_tibble %>%
          select("Protein.Group.Accessions", "Protein.Descriptions", "sequence_no_mod",
          "126":"131") %>%
          rename("Sequence" = sequence_no_mod) %>%
          add_column("Quan.Usage" = "Used", .after = "Sequence") %>%
          add_column("Quan.Info" = "Unique", .after = "Quan.Usage") %>%
          add_column("Isolation.Interference" = 30 , .after = "Quan.Info") 
view(mydat)
nrow(mydat)

dat <- read.csv("http://www.biostat.jhsph.edu/~kkammers/software/eupa/example_iTRAQ.csv")
view(dat)

(cha <- c("126","127N","127C","128N","128C","129N","129C","130N","130C","131"))
dim(mydat)
str(mydat)

read.peptides <- function(dat, cha){
  output <- NULL
  
  dat$Sequence <- as.character(dat$Sequence)
  dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
  dat$Quan.Usage <- as.character(dat$Quan.Usage)
  dat$Quan.Info <- as.character(dat$Quan.Info)
  dat$Isolation.Interference <- as.numeric(as.character(dat$Isolation.Interference))
  
  dat <- subset(dat, Isolation.Interference<=30)  
  dat <- subset(dat, Quan.Usage=="Used")
  dat <- subset(dat, Protein.Group.Accessions!="")
  dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))  
}

quantify.proteins <- function(dat, cha){
  e.function <- function(x, seq) tapply(x, seq, median)
  output <- NULL
  
  dat$Sequence <- toupper(dat$Sequence) # Capital letters
  accessions <- as.character(unique(dat$Protein.Group.Accessions))
  n.proteins <- length(accessions)
  n.cha <- length(cha)
  
  for(k in 1:n.proteins){
    id <- accessions[k]
    sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
    sdat[cha] <- log2(sdat[cha])
    sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
    pdat <- sdat[, -1]
    n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
    temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])          
    n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)    
    if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
    pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
    output <- rbind(output, pdat)
  }
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- round(output[,1:n.cha],3)
  row.names(output) <- accessions
  output <- as.data.frame(output)
  return(output)
}

eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  log2FC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

eb.fit.mult <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  log2FC <- fit.eb$coef[, "tr2"]
  df.0 <- rep(fit.eb$df.prior, n)
  df.r <- fit.eb$df.residual
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
  t.mod <- fit.eb$t[, "tr2"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "tr2"]
  q.ord <- q.mod <- rep(NA,n)
  ids <- which(!is.na(p.ord))
  k <- 0
  q1 <- qvalue(p.ord[!is.na(p.ord)])$q
  q2 <- qvalue(p.mod[!is.na(p.mod)])$q
  for(i in ids){ 
    k <- k+1
    q.ord[i] <- q1[k] 
    q.mod[i] <- q2[k]
  }
  results.eb.mult <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
  return(results.eb.mult)
}


(dat <- read.peptides(mydat, cha))
dim(dat)
####HIER DELEN??? door ref
dat <- quantify.proteins(dat, cha)
view(dat)
dat.onehit <- subset(dat, dat$n.peptides == 1) 
view(dat)
######OF HIER DELEN??? door ref
dim(dat.onehit)
dat <- subset(dat, dat$n.peptides != 1)
dim(dat)
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat[, 1:length(cha)],  ylim = c(-3, 3), main="Boxplot normalized Intensities")






nor <- c("126", "127N")
tum <- c("127C","128N","128C","129N","129C","130N","130C")

design <- model.matrix(~factor(c(2,2,2,2,2,1,1,1,1)))
design
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(dat[, c(nor,tum)], design)
str(res.eb)
view(res.eb)

rx <- c(-1, 1)*max(abs(res.eb$log2FC))*1.1
ry <- c(0, ceiling(max(-log10(res.eb$p.ord), -log10(res.eb$p.mod))))

par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
par(las=1, xaxs="i", yaxs="i")

plot(res.eb$log2FC, -log10(res.eb$p.ord), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of ordinary p-values")

plot(res.eb$log2FC, -log10(res.eb$p.mod), pch=21, bg="lightgrey", cex=0.9,
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of moderated p-values")



