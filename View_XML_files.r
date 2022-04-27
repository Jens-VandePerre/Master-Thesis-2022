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
library("XML")
library("methods")

wd <- setwd("/Users/jensvandeperre/Desktop/Outputs/PIA")
getwd()
list.files(wd)

#Load an XML file, PIA output
PIA <- xmlParse("test.xml")

# Exract the root node form the xml file.
rootnode <- xmlRoot(PIA)

# Find number of nodes in the root.
rootsize <- xmlSize(rootnode)

# Print the result.
print(rootsize)
print(rootnode[1])


xmldataframe <- xmlToDataFrame("test.xml")
print(xmldataframe)
