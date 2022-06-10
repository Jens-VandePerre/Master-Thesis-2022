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
library(broom)
library(cowplot)
library(data.table)

maxquant_path <- "/Users/jensvandeperre/Downloads/proteinGroups.txt" 


maxquant <- fread(maxquant_path, sep="\t") %>%
	filter(!Reverse=="+") %>%
	filter(!`Potential contaminant`=="+") %>%
	filter(!`Only identified by site`=="+") %>%
	mutate(`Protein IDs` = str_replace_all(`Protein IDs`, ";", "_"))
view(maxquant)
dim(maxquant)


maxquant_l <- maxquant %>% select(`Protein IDs`, starts_with("LFQ intensity"))
maxquant_t <- data.frame(t(maxquant_l[-1]))
colnames(maxquant_t) <- maxquant_t[1, ]
maxquant_t <- maxquant_t[-1, ]
view(maxquant_l)
dim(maxquant_l)


maxquant_t <- rownames_to_column(maxquant_t, "Condition")

maxquant_pca_input <- maxquant_t %>% as_tibble %>% 
	separate(Condition, into=c("Condition", "sample"), sep="_rep") %>%
	select(-c(sample)) %>%
	mutate_at(c(2:3096), as.numeric)
view(maxquant_pca_input)
dim(maxquant_pca_input)


maxquant_pca_input %>% 
  select(-Condition) %>%    # remove Condition column
  prcomp() ->             	# do PCA
  pca                     	# store result as `pca`

pca_data <- data.frame(pca$x, Condition = maxquant_pca_input$Condition)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) + 
  geom_point() +
  scale_color_manual(values = c("#FFA185", "#FF521B", "#80A7EF", "#123F91", "#84C28D", "#2F6036"))
#   scale_color_manual(values = c("#071736", "#881D07", "#7FA7F0", "#F88D77", "#1956C8", "#F33B16"))

# capture the rotation matrix in a data frame
rotation_data <- data.frame(
  pca$rotation,
  variable = row.names(pca$rotation)
)

# define a pleasing arrow style
arrow_style <- arrow(
  length = unit(0.05, "inches"),
  type = "closed"
)

# now plot, using geom_segment() for arrows and geom_text for labels
pca_biplot <- ggplot(rotation_data) + 
  geom_segment(aes(xend = PC1, yend = PC2), x = 0, y = 0, arrow = arrow_style) + 
  geom_text(aes(x = PC1, y = PC2, label = variable), hjust = 1, size = 3, color = "red") + 
  xlim(-1.25, .5) + 
  ylim(-.5, 1) +
  coord_fixed() # fix aspect ratio to 1:1

percent <- 100*pca$sdev^2 / sum(pca$sdev^2)

perc_data <- data.frame(percent = percent, PC = 1:length(percent))
percent_plot <- ggplot(perc_data, aes(x = PC, y = percent)) + 
  geom_col() + 
  geom_text(aes(label = round(percent, 1)), size = 4, vjust = -.5) + 
  ylim(0, 80) +
  scale_x_continuous(breaks = 1:9) # make sure each PC gets an axis tick

ggsave("", pca_plot, width = 8, height = 3)
ggsave("", pca_biplot, width = 8, height = 8)
ggsave("", percent_plot, width = 6, height = 3)