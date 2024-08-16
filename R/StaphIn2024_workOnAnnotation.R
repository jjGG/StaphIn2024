#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# read in protein annotation RData
load("p35269_firstAnalysis_SA6850_ProteinAnnotation.RData")
proteinAnnotation$description

library(stringr)
library(openxlsx)
library(dplyr)


# globals
fgczproject <- "p35269"
approach <- "firstAnalysis_SA6850"


# add more meta info columns to output
# blast related
#proteinAnnotation$description[100:103]
proteinAnnotation$blastOrthoAdd <- gsub(x = proteinAnnotation$description, pattern = ".*BLASTORTHO (.*)", replacement = "\\1")
proteinAnnotation$SAorthologue <- gsub(x = proteinAnnotation$blastOrthoAdd, pattern = "(.*) evalue=.*", replacement = "\\1")
proteinAnnotation$SAevalue <- as.numeric(gsub(x = proteinAnnotation$blastOrthoAdd, pattern = ".* evalue=(.*)", replacement = "\\1"))

# ncbi description related
proteinAnnotation$locusTag <- gsub(x = proteinAnnotation$description, pattern = ".*\\[locus_tag=(.*)\\] \\[protein=.*", replacement = "\\1")
proteinAnnotation$proteinDesc <- gsub(x = proteinAnnotation$description, pattern = ".*\\[protein=(.*)\\] \\[protein_id=.*", replacement = "\\1")
proteinAnnotation$SAprotID <- gsub(x = proteinAnnotation$description, pattern = ".*\\[protein_id=(.*)\\] \\[location=.*", replacement = "\\1")
proteinAnnotation$SAlocation <- gsub(x = proteinAnnotation$description, pattern = ".*\\[location=(.*)\\] \\[gbkey=.*", replacement = "\\1")

proteinAnnotation$blastOrthoAdd <- NULL
proteinAnnotation$UniprotAccOrtho <- unlist(lapply(str_split(string = proteinAnnotation$SAorthologue, pattern = "_"),FUN = "[",2))
proteinAnnotation$UniprotGNOrtho <- unlist(lapply(str_split(string = proteinAnnotation$SAorthologue, pattern = "_"),FUN = "[",3))

# read.xlsx(xlsxFile = totXlsx, sheet = "diff_exp_analysis")
eggNogg <- read.xlsx(xlsxFile = "../data/SA6850/eggNogg_SA6850_out.emapper.annotations.xlsx", sheet = "Sheet1", startRow = 3)

# combo <- left_join(x = phosRes, y = totRes, join_by("protein_Id" == "protAcc", "contrast" == "contrast"))
proteinAnnotationWEggNogg <- left_join(x = proteinAnnotation, y = eggNogg, join_by("protein_Id" == "query"))

(fN <- paste(fgczproject, "_", approach, "_proteinAnnotationWEggNoggNDBparsing", ".RData", sep = ""))
save(proteinAnnotationWEggNogg, file = fN)


