#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# here we develop an Rscript to analyze LFQ-DIA data analyzed w/ DIANN and using
# prolfqua and proflquapp

# here path to install packages
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("fgcz/prolfqua")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rentrez")

library(rentrez) #to fetch data from NCBI
library(prolfqua)
library(prolfquapp)
library(stringr)
library(dplyr)
library(rvest)
library(stringr)

# globals
fgczproject <- "p35269"
approach <- "firstAnalysis_SA6850"


load("p35269_firstAnalysis_SA6850_lfqDataObject.RData")
transformedData <- lfqdata

# what is the name of the working intensity
transformedData$config$table$get_response()


# what is the factors table
annotable <- transformedData$factors()
(annotable)
annotable$sampleName

transformedData$data

# define model with interaction
model <- strategy_lm("transformedIntensity ~ Genotype * GrowingCondition")

# what subject goes into the model?
# find out subject_Id name! transformedData$config$table$hkeysDepth()
mod <- build_model(transformedData$data, model, subject_Id = transformedData$config$table$hkeysDepth())

# where do we get small pvalues
mod$modelDF
mod$get_anova() # accept some warnings where some proteins have only few quants
mod$anova_histogram() # where do we get small pvalues?

annotable

# what contrasts are relevant?
##Focus first of comparison within Clone but between Environments
##Second compare clones with ancestor 6850
Contrast <- c("Growing_inPAvsTSB_GVancestor" = "`Genotype6850:GrowingConditionSN` - `Genotype6850:GrowingConditionTSB`",
              "Growing_inPAvsTSB_SB0403" = "`GenotypeSB0403:GrowingConditionSN` - `GenotypeSB0403:GrowingConditionTSB`",
              "Growing_inPAvsTSB_SB0804" = "`GenotypeSB0804:GrowingConditionSN` - `GenotypeSB0804:GrowingConditionTSB`",
              "Growing_inPAvsTSB_SB1002" = "`GenotypeSB1002:GrowingConditionSN` - `GenotypeSB1002:GrowingConditionTSB`",

              "Evolution_B0403vsAncestor_allSamples" = "`GenotypeSB0403` - `Genotype6850`",
              "Evolution_B1002vsAncestor_allSamples" = "`GenotypeSB1002` - `Genotype6850`")




# with moderation
contr <- Contrasts$new(mod, Contrast) %>% ContrastsModerated$new()
contrDF <- contr$get_contrasts()


# do imputing here
imputedContr_prot <- ContrastsMissing$new(lfqdata = transformedData, contrasts = Contrast)
colnames(imputedContr_prot$get_contrasts())

#merge
mergedResults_prot <- prolfqua::merge_contrasts_results(prefer = contr, add = imputedContr_prot)$merged
table(mergedResults_prot$contrast_result$modelName)
mergedResults_prot$contrast_result

# see where there are things going
# mod$anova_histogram()

pl <- contr$get_Plotter()
pl$volcano_plotly()

pl_imp <- mergedResults_prot$get_Plotter()
pl_imp$volcano_plotly()


(fN <- paste(fgczproject, "_", approach, "_ContrastResults", ".pdf",sep = ""))
pdf(fN, 15,15)
pl_imp$volcano()
pl_imp$ma_plot()
pl_imp$histogram_diff()
dev.off()


# write out results
(fN <- paste(fgczproject, "_", approach, "_ContrastResults", ".txt",sep = ""))
write.table(x = mergedResults_prot$contrast_result, file = fN, sep="\t", row.names = FALSE, col.names = TRUE)

contrDF <- mergedResults_prot$contrast_result
str(mergedResults_prot)

# merge with meta information
dim(contrDF)
(fN <- paste(fgczproject, "_", approach, "_ContrastObject", ".RData",sep = ""))
save(x = mergedResults_prot, file = fN)


#####Add additional protein info to data frame
# Function to extract AGU number from Protein_ID

extract_agu <- function(protein_id) {
  agu <- gsub(".*AGU([0-9]+)\\..*", "\\1", protein_id)
  return(paste0("AGU", agu))
}

# Function to fetch protein information from NCBI for S. aureus 6850
fetch_protein_info <- function(agu) {
  tryCatch({
    # Construct the URL for the GenBank format
    url <- paste0("https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=", agu, "&db=protein&report=genbank&conwithfeat=on&retmode=html&withmarkup=on")
    page <- read_html(url)
    content <- page %>% html_text()

    # Debugging: Print the first 500 characters of the page content
    print(substr(content, 1, 500))

    # Extract note, locus_tag, and gene using regular expressions
    note <- stringr::str_extract(content, "DEFINITION  (.+?)\\n")
    note <- ifelse(!is.na(note), stringr::str_trim(gsub("DEFINITION  ", "", note)), "Not found")

    locus_tag <- stringr::str_extract(content, "/locus_tag=\"([^\"]+)\"")
    locus_tag <- ifelse(!is.na(locus_tag), gsub("\"", "", stringr::str_extract(locus_tag, "\"([^\"]+)\"")), "Not found")

    gene <- stringr::str_extract(content, "/gene=\"([^\"]+)\"")
    gene <- ifelse(!is.na(gene), gsub("\"", "", stringr::str_extract(gene, "\"([^\"]+)\"")), "Not found")

    return(list(note=note, locus_tag=locus_tag, gene=gene))
  }, error = function(e) {
    print(paste("Error fetching data for AGU:", agu))
    return(list(note="Error", locus_tag="Error", gene="Error"))
  })
}



# Apply the functions to the test dataframe
contrDF <- contrDF %>%
  mutate(AGU = sapply(Protein_ID, extract_agu)) %>%
  rowwise() %>%
  mutate(protein_info = list(fetch_protein_info(AGU))) %>%
  mutate(
    note = protein_info$note,
    locus_tag = protein_info$locus_tag,
    gene = protein_info$gene
  ) %>%
  select(-AGU, -protein_info)

# Save the updated test dataframe to an Excel file
library(writexl)
write_xlsx(contrDF, "updated_contrDF_Saureus6850.xlsx")

# Print the results for inspection
print(contrDF_label)

# Print a summary of the results
print(paste("Total proteins processed:", nrow(contrDF_test)))
print(paste("Proteins found in S. aureus 6850:", sum(contrDF_test$note != "Not found" & contrDF_test$note != "Error")))
print(paste("Proteins not found:", sum(contrDF_test$note == "Not found")))
print(paste("Errors encountered:", sum(contrDF_test$note == "Error")))
#p3404_SAJE2_ResultsWithMeta <- left_join(contrDF, metaInfo_SAJE2)
#write.table(x = p3404_SAJE2_ResultsWithMeta, file = "p3404_SAJE2_fullResults_ancestor_growingInPAvsTSB_WithMetaInfoRedundant_wImputes_2022-08-17.txt", sep = "\t", row.names = FALSE)
#save(p3404_SAJE2_ResultsWithMeta, file = "p3404_SAJE2_results_ancestor_growingInPAvsTSB_wMetaInfo_wImputes_2022-08-17.RData")

