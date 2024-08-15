#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# here we develop an Rscript to analyze LFQ-DIA data analyzed w/ DIANN and using
# prolfqua and proflquapp

# here path to install packages

library(prolfqua)
library(prolfquapp)
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
Contrast <- c("Growing_inPAvsTSB_GVancestor" = "`Genotype6850:GrowingConditionSN` - `Genotype6850:GrowingConditionTSB`",
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


#p3404_SAJE2_ResultsWithMeta <- left_join(contrDF, metaInfo_SAJE2)
#write.table(x = p3404_SAJE2_ResultsWithMeta, file = "p3404_SAJE2_fullResults_ancestor_growingInPAvsTSB_WithMetaInfoRedundant_wImputes_2022-08-17.txt", sep = "\t", row.names = FALSE)
#save(p3404_SAJE2_ResultsWithMeta, file = "p3404_SAJE2_results_ancestor_growingInPAvsTSB_wMetaInfo_wImputes_2022-08-17.RData")




