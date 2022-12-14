#------------------------------------------------------------------------------
# Script calculates logistic regression models for metabolites and lipids
# previously found associated with CMV infection. Creates result files with
# estimates, p-values, and confidence intervals.
#
# Author: Kirstine K. Rasmussen
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Logistic regression models of previously found associations with CMV
#------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())
cat("\014")

# Load packages
library(tidyverse)
library(magrittr)
library(broom)
library(cowplot)
library(ggsignif)

#------------------------------------------------------------------------------
# Load data and prepare for analysis
#------------------------------------------------------------------------------

# Load in metabolites and select known metabolites and clinical data
m <- read_tsv("Data/01_Metabolites_clinical_cross_sectional.tsv") %>% 
  dplyr::select("SAMPLE", "PATIENT", "Sex", "ttype", "age_at_T", "CMV_risk", "Myeloabl",
                "Group", "Group description", "setnumber",
                "glutamine", "phenylalanine", "tryptophan", "kynurenine", 
                "quinolinate", "taurine", "trimethylamine N-oxide", "lactate", 
                "lysine", "choline", "alanine")

# Load in lipids and select known lipids
l <- read_tsv("Data/01_Lipids_clinical_cross_sectional.tsv") %>%
  dplyr::select("SAMPLE", "Total FFA"  )

# Join known metabolites and lipids
ml <- left_join(m, l, by = "SAMPLE")

# Known analytes
analytes <- c("glutamine", "phenylalanine", "tryptophan", "kynurenine", 
              "quinolinate", "taurine", "trimethylamine N-oxide", "lactate", 
              "lysine", "choline", "alanine", "Total FFA")

# Create long-format data frame
ml_long <- ml %>% gather(
  "analyte",
  "measurement",
  all_of(analytes)) %>%
  group_by(analyte) %>%
  mutate(analyte = factor(analyte))

#------------------------------------------------------------------------------
# Multivariable model
#------------------------------------------------------------------------------

# Multivariable model adjusted for conditioning regimen, sex, age, and CMV risk score
multiv <- lapply(levels(ml_long$analyte),
                  function(x){
                    glm(Group ~ measurement  + factor(ttype) + factor(Sex) + age_at_T + factor(CMV_risk), 
                        family=binomial, data=subset(ml_long,analyte==x)) %>%
                      tidy(., conf.int=T,exponentiate=T) %>%
                      mutate(analyte=x) %>%
                      select(analyte,everything()) %>%
                      .[.$term == "measurement",]
                  })

# Gather results in table
multiv_res <- bind_rows(multiv) %>%
  write_tsv('Data/03_MultivariableReg_results.tsv')

# Find significant results
multiv_res[multiv_res$p.value < 0.1,]
