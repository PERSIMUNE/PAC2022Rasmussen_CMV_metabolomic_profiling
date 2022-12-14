#------------------------------------------------------------------------------
# Script loads cleaned up data and calculates summary statistics for the
# first table in the paper.
#
# Author: Kirstine K. Rasmussen
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Load cleaned data
#------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())
cat("\014")

# Load packages
library(tidyverse)

# Load cohort (one sample for 68 patients)
metData <- read_tsv("Data/01_Metabolites_clinical_cross_sectional.tsv")

# Sample collection dates
summary(metData$CEP_D)

# Make separate datasets for the cases and controls
cases <- subset(metData, Group==1)
controls <- subset(metData, Group==0)

#------------------------------------------------------------------------------ 
# Sex
#------------------------------------------------------------------------------
# 1 = male, 2 = female
table(metData$Sex);table(metData$Sex)[[1]]/dim(metData)[1]*100 #male
table(cases$Sex);table(cases$Sex)[[1]]/dim(cases)[1]*100 #male
table(controls$Sex);table(controls$Sex)[[1]]/dim(controls)[1]*100 #male

Sex = ggplot(metData) + 
  aes(x = `Group description`, fill = as.character(Sex)) +
  geom_bar() + 
  ggtitle("Sex of patients") +
  labs(fill = "Sex"); Sex

# Chi square as Sex is categorical
chisq.test(metData$Sex, metData$`Group description`)

#------------------------------------------------------------------------------
# Age at transplantation
#------------------------------------------------------------------------------
summary(metData$age_at_T)
summary(cases$age_at_T)
summary(controls$age_at_T)

# Test if age is normally distributed
shapiro.test(metData$age_at_T)$p.value #p<0.05 -> NOT normally distributed

age = ggplot(metData) +
  aes(x = age_at_T, fill = `Group description`) +
  geom_bar() +
  ggtitle("Age at aHSCT");age

# Statistical test of non-normally distributed continuous variable
wilcox.test(age_at_T ~ `Group description`, data = metData)

#------------------------------------------------------------------------------
# Sample collected days post aHSCT
#------------------------------------------------------------------------------
summary(metData$days_post_T)
summary(cases$days_post_T)
summary(controls$days_post_T)

# Test if sample collection-times are normally distributed
shapiro.test(metData$days_post_T)$p.value #p<0.05 -> NOT normally distributed

samples = ggplot(metData) +
  aes(x = days_post_T, fill = `Group description`) +
  geom_bar() +
  theme_light() +
  xlab("days post-aHSCT") +
  theme(plot.title=element_text(size=16,face="bold")) +
  theme(axis.title=element_text(size=13,face="bold")) +
  theme(axis.text=element_text(size=12)) +
  theme(legend.text=element_text(size=13)) +
  theme(legend.title=element_text(size=13,face="bold")) +
  scale_fill_discrete(name="Group",labels=c("CMV infection","No CMV infection")); samples

# Wilcoxon test of non-normally distributed continuous variable
wilcox.test(days_post_T ~ `Group description`, data = metData)

#------------------------------------------------------------------------------
# Fist CMV in days post-aHSCT
#------------------------------------------------------------------------------
summary(cases$firstCMV)
summary(controls$firstCMV)

# Test if CMV occurrence is normally distributed
shapiro.test(metData$firstCMV)$p.value #p<0.05 -> NOT normally distributed

CMV = ggplot(metData) +
  aes(x = firstCMV, fill = `Group description`) +
  geom_bar() +
  ggtitle("First CMV occurrence post-aHSCT"); CMV

# Plot days between sample collection and CMV onset
pCMVpostaHSCT = ggplot(cases) + 
  aes(x = firstCMV-days_post_T) +
  geom_bar() +
  ggtitle("Days between sample collection and CMV onset") +
  scale_x_continuous(breaks = seq(5,60,5)); pCMVpostaHSCT

summary(cases$firstCMV-cases$days_post_T)

#------------------------------------------------------------------------------
# Graft origin
#------------------------------------------------------------------------------
print("Bone marrow");dim(metData[metData$`Bone marrow`=='yes',])[1];dim(metData[metData$`Bone marrow`=='yes',])[1]/length(metData$SAMPLE)*100
print("Peripheral blood");dim(metData[metData$`Peripheral blood`=='yes',])[1];dim(metData[metData$`Peripheral blood`=='yes',])[1]/length(metData$SAMPLE)*100

print("Bone marrow");dim(cases[cases$`Bone marrow`=='yes',])[1];dim(cases[cases$`Bone marrow`=='yes',])[1]/length(cases$SAMPLE)*100
print("Peripheral blood");dim(cases[cases$`Peripheral blood`=='yes',])[1];dim(cases[cases$`Peripheral blood`=='yes',])[1]/length(cases$SAMPLE)*100

print("Bone marrow");dim(controls[controls$`Bone marrow`=='yes',])[1];dim(controls[controls$`Bone marrow`=='yes',])[1]/length(controls$SAMPLE)*100
print("Peripheral blood");dim(controls[controls$`Peripheral blood`=='yes',])[1];dim(controls[controls$`Peripheral blood`=='yes',])[1]/length(controls$SAMPLE)*100

chisq.test(metData$`Peripheral blood`, metData$`Group description`)
chisq.test(metData$`Bone marrow`, metData$`Group description`)

# Fishers exact test because 
fisher.test(metData$`Bone marrow`, metData$`Group description`)

#------------------------------------------------------------------------------
# Conditioning
#------------------------------------------------------------------------------
# 1 = myeloablative, 2 = nonmyeloablative
table(metData$CEP_SPEC)
table(metData$CEP_SPEC)[[1]]/dim(metData)[1]*100 #myeloablative
table(metData$CEP_SPEC)[[2]]/dim(metData)[1]*100 #nonmyeloablative

table(cases$CEP_SPEC)
table(cases$CEP_SPEC)[[1]]/dim(cases)[1]*100 #myeloablative
table(cases$CEP_SPEC)[[2]]/dim(cases)[1]*100 #nonmyeloablative

table(controls$CEP_SPEC)
table(controls$CEP_SPEC)[[1]]/dim(controls)[1]*100 #myeloablative
table(controls$CEP_SPEC)[[2]]/dim(controls)[1]*100 #nonmyeloablative

ttype = ggplot(metData) +
  aes(x = CEP_SPEC, fill = `Group description`) +
  geom_bar() +
  ggtitle("Conditioning regimen pre-aHSCT"); ttype

# Chi square as ttype is categorical
chisq.test(metData$CEP_SPEC, metData$`Group description`)

#------------------------------------------------------------------------------
# Serostatus
#------------------------------------------------------------------------------
table(cases$serostatus)
table(controls$serostatus)

#------------------------------------------------------------------------------
# CMVrisk
#------------------------------------------------------------------------------
# 1 = low, 2 = intermediate, 3 = high
# CMV risks of 0 have been merged with CMV risk = 1 group
table(metData$CMV_risk)
table(metData$CMV_risk)[[1]]/length(metData$SAMPLE)*100 #low
table(metData$CMV_risk)[[2]]/length(metData$SAMPLE)*100 #intermediate
table(metData$CMV_risk)[[3]]/length(metData$SAMPLE)*100 #high

table(cases$CMV_risk)
table(cases$CMV_risk)[[1]]/length(cases$SAMPLE)*100 #low
table(cases$CMV_risk)[[2]]/length(cases$SAMPLE)*100 #intermediate
table(cases$CMV_risk)[[3]]/length(cases$SAMPLE)*100 #high

table(controls$CMV_risk)
table(controls$CMV_risk)[[1]]/length(controls$SAMPLE)*100 #low
table(controls$CMV_risk)[[2]]/length(controls$SAMPLE)*100 #intermediate
table(controls$CMV_risk)[[3]]/length(controls$SAMPLE)*100 #high

CMVrisk = ggplot(metData) +
  aes(x = as.character(CMV_risk), fill = `Group description`) +
  geom_bar() +
  ggtitle("CMVrisk - based on serostatus"); CMVrisk

# Fisher's exact test as cmvrisk is categorical w. small groups
fisher.test(metData$CMV_risk, metData$`Group description`)
