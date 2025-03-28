library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(ggthemes)
library(viridis)
library(MASS)
library(lme4)
library(tidyverse)
library(multcomp) 
library(multcompView)
library(FSA)
library(stats)
library(PMCMRplus)
library(PMCMR)
library(pscl)
library(emmeans)

library(DHARMa)
library(glmmTMB)

library(performance)

source("../graphing functions.r")
source("../Data Processing.r")

# Reformatted version of Study 1. 
# Maybe just use for some of the analysis? Live scale analysis?
# Make study 1 a folder with sub-files instead?

#### (1) Data Pre-processing ####
trt_labels <- c("May 2-3", "May 17", "May 28", "No Treatment")

# Read in file
scalecount <- read.csv(file="EHS count 2024 v7 (study 1).csv", strip.white=TRUE)
colnames(scalecount) <- c("Label", "County","Twigab","Date","Counter","Livescale1","Deadscale1",
                          "Livescale2","Deadscale2","Livescale3","Deadscale3","Prespara","Presfungus",
                          "Presscalenewgr","encarsia")

scalecount <- process_scalecount(scalecount)

# also change county to factor
scalecount$County <- as.factor(scalecount$County)

# Getting first count - July
scalecount_july <- subset(scalecount, grepl('July', Date)) # July data
scalecount_nov <- subset(scalecount, grepl('November', Date)) # November data

# Treatment survival means
tmeans_july <- get_treatment_survival_means(scalecount_july)
tmeans_nov <- get_treatment_survival_means(scalecount_nov)

##### Specialized parasitism and fungus preprocessing #####
# Drop twigs with 0 scale
scalecount_para_july <- scalecount_july[rowSums(scalecount_july[, 6:11] == 0) < 2,]
scalecount_para_july <- 
  scalecount_para_july %>% 
  drop_na(Label) #some samples only have one twig

scalecount_para_nov <- scalecount_nov[rowSums(scalecount_nov[, 6:11] == 0) < 2,]
scalecount_para_nov<- 
  scalecount_para_nov %>% 
  drop_na(Label) #some samples only have one twig

#### (2) Parasitism models ####
# Use mixed models to account for random block and tree effects

##### July #####
para_mod_july <- glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                         data=scalecount_para_july,
                         family=binomial)
simr_para_mod_july <- simulateResiduals(para_mod_july)

##### November #####
para_mod_nov <-  glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                        data=scalecount_para_nov,
                        family=binomial)
simr_para_mod_nov <- simulateResiduals(para_mod_nov)

#### (3) Parasitism Analysis ####

##### July #####
emm_para_july <- emmeans(para_mod_july, "Treatment")
pairs(emm_para_july, type="response", adjust="BH")
july_para_means_comp <- pairs(regrid(emm_para_july), adjust="BH")
# CI for pairwise comp
confint(july_para_means_comp)
# CI for means
confint(emm_para_july, type="response")

##### November #####
emm_para_nov <- emmeans(para_mod_nov, "Treatment")
pairs(emm_para_nov, type="response", adjust="BH")
nov_para_means_comp <- pairs(regrid(emm_para_nov), adjust="BH")
# Confint for pairwise comparison
confint(nov_para_means_comp)
# Confint for means
confint(emm_para_nov, type="response")

#### (4) Fungus models ####
##### July #####
fung_mod_july <- glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                         data=scalecount_para_july,
                         family=binomial)
simr_fung_mod_july <- simulateResiduals(para_mod_july)

##### November #####
fung_mod_nov <-  glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                         data=scalecount_para_nov,
                         family=binomial)
simr_fung_mod_nov <- simulateResiduals(fung_mod_nov)

#### (5) Fungus Analysis ####
##### July #####
emm_fung_july <- emmeans(fung_mod_july, "Treatment")
pairs(emm_fung_july, type="response", adjust="BH")
july_fung_means_comp <- pairs(regrid(emm_fung_july), adjust="BH")
# CI for pairwise comp
confint(july_fung_means_comp)
# CI for means
confint(emm_fung_july, type="response")

##### November #####
emm_fung_nov <- emmeans(fung_mod_nov, "Treatment")
pairs(emm_fung_nov, type="response", adjust="BH")
nov_fung_means_comp <- pairs(regrid(emm_fung_nov), adjust="BH")
# Confint for pairwise comparison
confint(nov_fung_means_comp)
# Confint for means
confint(emm_fung_nov, type="response")

