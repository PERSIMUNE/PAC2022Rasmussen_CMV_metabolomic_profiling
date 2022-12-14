#------------------------------------------------------------------------------
# Script performs clustering of metabolites/lipids by use of the WGCNA 
# package. Code also aids in the selection of parameter values for the
# analysis and produces various figures and a final file containing all
# metabolites/lipids in the clusters.
#
# Author: Kirstine K. Rasmussen
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Load libraries and set up script
#------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())
cat("\014")

# Load packages
library(WGCNA)
library(tidyverse)
library(cowplot)
library(zeallot) # for %<-%
library(plotly)
library(ggplot2)
library(heatmaply)
library(htmlwidgets)
options(stringsAsFactors = FALSE);

# Enables multi-threating within WGCNA to speed up calculations
enableWGCNAThreads()

# Load custom functions for WGCNA
source('Cross_sectional_WGCNA_functions.R')

#------------------------------------------------------------------------------
# Select dataset and import data
# The script is run for metabolomics and lipidomics separately by selecting
# the dataType parameter
#------------------------------------------------------------------------------

# Selecct dataset
dataType <- "CS_met"
#dataType <- "CS_lip"

if (dataType == "CS_met") {
  
  # Import data
  metData <- read_tsv("Data/01_Metabolites_clinical_cross_sectional.tsv")
  
  # Transform to dataframe and remove columns "SAMPLE" 
  # and all clinical data
  log2data = as.data.frame(metData[, -c(1,924:952)]) %>% log2(.)
  patientIDs <- metData$PATIENT
  
  # Plot distributions
  plot_before <- data.frame(val = as.vector(as.matrix(metData[, -c(1,924:952)]))) %>%
    ggplot(., aes(val)) +
    geom_histogram() +
    theme_light() +
    xlab("metabolite abundances")
  
  plot_after <- data.frame(val = as.vector(as.matrix(log2data))) %>%
    ggplot(., aes(val)) +
    geom_histogram() +
    theme_light() +
    xlab("log2 transformed metabolite abundances")
  
  plot_grid(plot_before, plot_after, labels = c('A','B')) %>%
    ggsave("Figures/04_Data_distribution_met.png",., 
           width = 2700, height = 1100, units = "px")
  
  # Check for non-varied analytes
  datExpr0 <- data_init(log2data, patientIDs)
}
if (dataType == "CS_lip") {

  # Import data (log2 fold changes and clinical data)
  lipData <- read_tsv("Data/01_Lipids_clinical_cross_sectional.tsv")

  # Transform to dataframe and remove columns "SAMPLE", 
  # total lipids columns, and all clinical columns
  log2data <- as.data.frame(lipData[, -c(1,935:978)]) %>% log2(.)
  patientIDs <- lipData$PATIENT
  
  # Plot distributions
  plot_before <- data.frame(val = as.vector(as.matrix(lipData[, -c(1,935:978)]))) %>%
    ggplot(., aes(val)) +
    geom_histogram() +
    theme_light() +
    xlab("lipid abundances")
  
  plot_after <- data.frame(val = as.vector(as.matrix(log2data))) %>%
    ggplot(., aes(val)) +
    geom_histogram() +
    theme_light() +
    xlab("log2 transformed lipid abundances")
  
  plot_grid(plot_before, plot_after, labels = c('C','D')) %>%
    ggsave("Figures/04_Data_distribution_lip.png",., 
           width = 2700, height = 1100, units = "px")
  
  datExpr0 <- data_init(log2data, patientIDs)
}

# Show range of abundances
min(datExpr0);max(datExpr0)

#------------------------------------------------------------------------------
# Outlier detection of samples
#------------------------------------------------------------------------------

# Visualise samples and remove potential outliers
# - Run first time to visualise samples
# - Run second time with appropriate cutoff that removes outliers
# Samples cut by the line will be removed
datExpr <- check_outliers(datExpr0, cutoff = 72)

#------------------------------------------------------------------------------
# Include trait data
#------------------------------------------------------------------------------

# Name of file with clinical variables (samples in rows, variables in columns)
clinicalData <- "Data/01_Clinical_data.tsv"

# Clean up clinical data and save relevant clinical traits
datTraits <- importClinical(clinicalData, datExpr)

# Make dendrogram of samples with heatmap of clinical traits
sampleTraitDend(datExpr, datTraits)

# Save expression matrix and clinical traits matrix
save(datExpr, datTraits, 
     file = paste("Rdata/04_WGCNA_",dataType,"_datExpr_datTraits.RData", sep='')); print(
       paste("The following file was created: Rdata/04_WGCNA_",dataType,"_datExpr_datTraits.RData", sep=''))

#------------------------------------------------------------------------------
# Determine power parameters
#
# If previous part of the script has been run before,
# you can start here after defining the dataType variable
#------------------------------------------------------------------------------

# Load data files
load(file = paste("Rdata/04_WGCNA_",dataType,"_datExpr_datTraits.RData", sep=''))

# Takes a little while to run
# Selects the lowest power resulting in a model fit of >0.9
# OR if not possible, a power of 30
power <- determine_power(datExpr); sprintf('A power of %d was selected', power)

#------------------------------------------------------------------------------
# minModuleSize determination
#------------------------------------------------------------------------------

# Load data if analysis has been run before
load(file = paste("Rdata/04_WGCNA_",dataType,"_datExpr_datTraits.RData", sep=''))
load(file = paste("Rdata/04_WGCNA_",dataType,"_ModuleMembership.RData", sep=''))

# Set parameters: Biweight mid-correlation and minimum module sizes 5-20
corType = "bicor"
modList = c(5:20)

# Loop for making networks with varying minimum module sizes
MMavgAll <- vector()
numMod <- vector()
for(minMod in modList){
  # Make network
  c(moduleColors, moduleContent, MEs, tree) %<-% 
    network_construction(datExpr, corType, power, minMod, dataType)

  # Make module-trait correlation heatmap
  traitHeatmap(datExpr, datTraits, MEs, moduleColors, minMod, dataType)
  
  # Make module-membership plot
  MMmax <- moduleMembershipPlot(datExpr, MEs, moduleColors, minMod, dataType, plot=T)
  
  # Save average module membership
  MMavgAll <- c(MMavgAll,mean(MMmax))
  
  # Save number of modules produced by branch cutting
  numMod <- c(numMod, length(unique(moduleColors)))
}

# Make plot of average max module memberships
MMmat <- plotAvgMM(numMod, MMavgAll, dataType)
writexl::write_xlsx(data.frame(MMmat), paste("Data/MMmat_",dataType,".xlsx"))

# Select appropriate minMod parameter by assessing results from the just run section
minMod <- 10

save(MMavgAll, numMod, power, minMod, corType, 
     file = paste("Rdata/04_WGCNA_",dataType,"_ModuleMembership.RData", sep='')); print(
       paste("The following file was created: Rdata/04_WGCNA_",dataType,"_ModuleMembership.RData", sep=''))

#------------------------------------------------------------------------------
# Final WGCNA with selected parameters
#------------------------------------------------------------------------------

# Load data if analysis has been run before
load(file = paste("Rdata/04_WGCNA_",dataType,"_datExpr_datTraits.RData", sep=''))
load(file = paste("Rdata/04_WGCNA_",dataType,"_ModuleMembership.RData", sep=''))

# Make network
c(moduleColors, moduleContent, MEs, geneTree) %<-% 
  network_construction(datExpr, corType, power, minMod, dataType)

save(datExpr, datTraits, MEs, power, minMod, dataType, moduleColors, moduleContent,
     file = paste("Rdata/04_WGCNA_",dataType,"_all.RData", sep='')); print(
       paste("The following file was created: Rdata/04_WGCNA_",dataType,"_all.RData", sep=''))

#------------------------------------------------------------------------------
# Module assessment
#
# If final WGCNA has been run before, you can start here
#------------------------------------------------------------------------------

load(paste("Rdata/04_WGCNA_",dataType,"_all.RData", sep=''))

# Find hub-nodes in each module
hubs <- chooseTopHubInEachModule(datExpr, moduleColors, power = power, 
                                 type = "signed", );hubs

# Make annotation dataframe
if (dataType == "CS_met") {
  annotFile = "Data/01_Metabolites_annotation.tsv"}

if (dataType == "CS_lip") {
  annotFile = "Data/01_Lipids_annotation.tsv"}

c(annot, IDs2annot) %<-% readAnnotationFile(annotFile, datExpr)

# Make long-format table of module results
all_mods <- data.frame("Module"        = moduleColors,
                       "BIOCHEMICAL"   = annot$BIOCHEMICAL[IDs2annot],
                       "SUPER_PATHWAY" = annot$SUPER_PATHWAY[IDs2annot],
                       "SUB_PATHWAY"   = annot$SUB_PATHWAY[IDs2annot],
                       "HMDB_ID"       = annot$HMDB_ID[IDs2annot],
                       "KEGG_ID"       = annot$KEGG_ID[IDs2annot])
openxlsx::write.xlsx(all_mods, file = paste('Data/04_WGCNA_',dataType,'_results.xlsx',sep=''))     

# Number of analytes in each module
countModules <- all_mods %>% count("Module")

# Make module-trait correlation heatmaps
traitHeatmap(datExpr, datTraits, MEs, moduleColors, minMod, dataType)
interactiveHeatmap <- interactiveHeatmap(datExpr, datTraits, MEs, minMod, dataType)

# Make dataframe of superpathways in modules
c(df_sup, SUP_PATH) %<-% makeSuperPathwayDataframe(annot, IDs2annot, MEs, dataType)
openxlsx::write.xlsx(df_sup, file = paste('Data/04_PathwaysInMods_',dataType,'.xlsx',sep='')) 

# Make barplot of how pathways are distributed in modules
modsInPathwaysPlot(df_sup, SUP_PATH, minMod, dataType)

# Make barplot of which pathways are in each module 
pathwaysInModsPlot(df_sup, moduleColors, minMod, dataType)

# Make module-membership plot
MMmax <- moduleMembershipPlot(datExpr, MEs, moduleColors, minMod, dataType, plot = TRUE)

# Load TOM
load(file = 'blockwiseTOM-block.1.RData')
# Export networks to Cytoscape files
exportNetworkToCytoscape(TOM, 
                         edgeFile = paste('Data/04_',dataType,'_cytoscape_edge_',minMod,'.txt', sep=''),
                         nodeFile = paste('Data/04_',dataType,'_cytoscape_node_',minMod,'.txt', sep=''),
                         nodeAttr = moduleColors, altNodeNames = SUP_PATH, weighted = T,
                         threshold = 0.25, 
                         nodeNames = names(datExpr))

# Make excel files with analytes and IDs per module
exportModuleIDs(annot, IDs2annot, minMod, dataType)
