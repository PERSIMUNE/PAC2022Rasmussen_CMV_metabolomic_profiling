#------------------------------------------------------------------------------
# Script sets up the workspace and performs pre-processing of data files
# from Metabolon (metabolite and lipid abundances).
#
# Author: Kirstine K. Rasmussen
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup script with necessary packages to run all scripts
#------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())
cat("\014")

# Used packages
list_of_packages <- c('tidyverse','plyr','readxl','broom',
                      'cowplot','ggsignif','gtools',
                      'htmlwidgets','magrittr','caret',
                      'plotly','WGCNA', 'zeallot', 'ggplot2')

# Install needed packages only if not installed
new_packages <- list_of_packages[!(list_of_packages %in% 
                                     installed.packages ()[, "Package"])]
if(length(new_packages)) {
  install.packages(new_packages) }; rm(list_of_packages)

# Make folders for placing output files
ifelse(!dir.exists("Data"), dir.create("Data"), FALSE) # For constructed data files
ifelse(!dir.exists("Figures"), dir.create("Figures"), FALSE) # For figures
ifelse(!dir.exists("Rdata"), dir.create("Rdata"), FALSE) # For Rdata files

#------------------------------------------------------------------------------
# Data pre-processing of metabolomics and lipidomics
#------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())
cat("\014")

# Load packages
library(plyr)
library(readxl)
library(tidyverse)
library(caret) # nearZeroVar
library(cowplot)
library(broom)

# Read in files from Metabolon 
# (median scaled and not detected biochemicals replaced with the lowest value)
metabolites <- read_xlsx("../CMVmetabo/CORI-02-18MLCLP+ CDT_ANONYMIZED_median scaled_metabolites.xlsx")
lipids <- read_xlsx("../CMVmetabo/CORI-02-18MLCLP+ CDT_ANONYMIZED_median scaled_lipids.xlsx")

#------------------------------------------------------------------------------
# Extract pathway data
#------------------------------------------------------------------------------

# Dataframe of metabolites and corresponding super- and subpathways
data.frame(COMP_ID = metabolites$COMP_ID,
           BIOCHEMICAL = metabolites$BIOCHEMICAL,
           HMDB_ID = metabolites$HMDB,
           PUBCHEM_ID = metabolites$PUBCHEM,
           KEGG_ID = metabolites$KEGG...13,
           SUPER_PATHWAY = metabolites$SUPER_PATHWAY,
           SUB_PATHWAY = metabolites$SUB_PATHWAY) %>%
  write_tsv("Data/01_Metabolites_annotation.tsv")

# Dataframe of lipids and corresponding super- and subpathways
data.frame(COMP_ID = lipids$COMP_ID,
           BIOCHEMICAL = lipids$BIOCHEMICAL,
           HMDB_ID = lipids$`Sample   HMDB_ID`,
           KEGG_ID = lipids$KEGG,
           SUPER_PATHWAY = lipids$SUPER_PATHWAY,
           SUB_PATHWAY = lipids$SUB_PATHWAY) %>%
  write_tsv("Data/01_Lipids_annotation.tsv")

#------------------------------------------------------------------------------
# Clean up expression data
#------------------------------------------------------------------------------

# Select biochemical name and sample data
# Transpose to get patients as rows and biochemical as columns
m <- metabolites %>% select(2, 17:176) %>%
  gather(SAMPLE, valname, 2:161) %>% 
  spread(BIOCHEMICAL, valname)

l <- lipids %>% select(2, 10:169) %>%
  gather(SAMPLE, valname, 2:161) %>% 
  spread(BIOCHEMICAL, valname)

# Identify and remove zero and near-zero variance columns
m_var <- names(m[nearZeroVar(m)])
sprintf('These %d columns are removed:', length(m_var));m_var
m_2 <- dplyr::select(m, -c(all_of(m_var)))

# Identify and remove zero and near-zero variance columns
l_var <- names(l[nearZeroVar(l)])
sprintf('These %d columns were removed:', length(l_var));l_var
l_2 <- dplyr::select(l, -c(all_of(l_var)))

#------------------------------------------------------------------------------
# Merge with clinical data
#------------------------------------------------------------------------------

# Load baseline and other clinical data
# Remove PATIENT1 column as this is identical to PATIENT column
# Change all CMVrisk = 0 to CMVrisk = 1
clinical <- read_tsv("../CMVmetabo/metaboCMV_casesandcontrols_characteristics.tsv") %>%
  dplyr::rename("Original_SAMPLE_ID" = "SAMPLE_ID", "SAMPLE" = "SAMPLE_ID_CHOSEN_BY_BIOBANK",
                "Sex" = "GENDER") %>%
  .[-c(18)] %>%
  mutate(CMV_risk_old = CMV_risk) %>%
  mutate(CMV_risk = ifelse(CMV_risk == 0, 1, CMV_risk));names(clinical)

# Check if unknown CMV risk samples can be grouped with low risk
table(clinical$CMV_risk_old);table(clinical$CMV_risk)

logres <- clinical %>%
  mutate(CMV_risk_old = relevel(as.factor(CMV_risk_old), '1')) %>%
  do(CMV_risk_old_logres = glm(Group ~ CMV_risk_old, data = .)) %>%
  ungroup %>%
  mutate(mapping = map(CMV_risk_old_logres, tidy, exponentiate = T, conf.int = T)) %>%
  unnest(mapping) %>%
  mutate_if(is.numeric, round, digits=2) %>%
  write_tsv(., "01_CMV_risk_LogisticRegression.tsv")

# Import clinical data about aHSCT
# Select columns, rename column for date of bone marrow transplant,
aHSCT <- read_tsv("../CMVmetabo/haem_clindata_incGvHD.tsv") %>%
  select(PATIENT, `Date of BMT (dd/mm/yyyy)`, 
         DiseaseNameIBMTR, `Total dose of TBI in cGy`,
         `Bone marrow`, `Ubilical cord`, `Peripheral blood`,
         Myeloabl, `Date of death`, `Primary cause of death`) %>% 
  dplyr::rename("Date_BMT" = "Date of BMT (dd/mm/yyyy)")

# Join clinical and aHSCT data
clinical_aHSCT<- left_join(clinical, aHSCT, by="PATIENT")

# Join expression and clinical data
clinical_aHSCT$SAMPLE <- as.character(clinical$"SAMPLE")
m_3 <- left_join(m_2, clinical_aHSCT, by = "SAMPLE")
l_3 <- left_join(l_2, clinical_aHSCT, by = "SAMPLE")

#------------------------------------------------------------------------------
# Remove patients and samples not needed
#------------------------------------------------------------------------------

# Remove patients not fitting within study description
removePatients <- c("588781","588914","588927","588932",
                    "588835","588848","588743","588844")

m_4 <- m_3 %>%
  mutate(PATIENT = as.character(m_3$PATIENT)) %>%
  anti_join(tibble(PATIENT = removePatients))

l_4 <- l_3 %>%
  mutate(PATIENT = as.character(l_3$PATIENT)) %>%
  anti_join(tibble(PATIENT = removePatients))

# Remove duplicate sample found by manual assessment of data
# And remove extra controls (8 samples)
m_5 <- m_4 %>% filter(Original_SAMPLE_ID != '26936841') %>%
  filter(setnumber != 'n/a') %>%
  filter(days_post_T < 40)

l_5 <- l_4 %>% filter(Original_SAMPLE_ID != '26936841') %>%
  filter(setnumber != 'n/a') %>%
  filter(days_post_T < 40)

#------------------------------------------------------------------------------
# Select latest collected sample per patient
#------------------------------------------------------------------------------

m_6 <- ddply(m_5, .(PATIENT), function(x) x[which.max(x$days_post_T),]) 

l_6 <- ddply(l_5, .(PATIENT), function(x) x[which.max(x$days_post_T),])

# Plot sample collection days relative to aHSCT
ggplot(m_6) + 
  geom_bar(mapping = aes(x = days_post_T), fill = "darkblue", alpha = 0.8) +
  xlab("Days post aHSCT") +
  ggtitle("Sample collection in days post-aHSCT") +
  scale_x_continuous(breaks = seq(-5, 60, by = 5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,22)) +
  theme_light()

# Save clinical data (whole cohort)
as.data.frame(m_6[924:952]) %>% write_tsv("Data/01_Clinical_data.tsv")

#------------------------------------------------------------------------------
# Write final data files
#------------------------------------------------------------------------------

write_tsv(m_6, "Data/01_Metabolites_clinical_cross_sectional.tsv")

write_tsv(l_6, "Data/01_Lipids_clinical_cross_sectional.tsv")
