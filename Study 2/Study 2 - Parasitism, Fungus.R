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

source("../graphing functions.r")
source("../Data Processing.r")

# Reformatted version of Study 2.

#### (1) Data Pre-processing ####

# Reading in the data
# Not sure why but there is an extra column of missing values at the end of the data set
# This is removed below along with the 'notes' column
scalecount2 <- read.csv(file = "distance numbers v7 (study 2).csv", strip.white=TRUE)
colnames(scalecount2) <- c("Label","Type", "Twigab","Date","Livescale1","Deadscale1","Livescale2",
                           "Deadscale2","Livescale3","Deadscale3","Prespara",
                           "Presscalenewgr","Presfungus","encarsia","notes","extra")

# Remove notes column and random extra column
scalecount2 <- subset(scalecount2, select = -c(notes, extra))

# Process data set
scalecount2 <- process_scalecount(scalecount2)

# July and November data
scalecount2_july <- subset(scalecount2, grepl('July', Date))
scalecount2_nov <- subset(scalecount2, grepl('November', Date))

# Treatment survival means
tmeans2_july <- get_treatment_survival_means(scalecount2_july)
tmeans2_nov <- get_treatment_survival_means(scalecount2_nov)

##### Specialized parasitism and fungus preprocessing #####

# Should look at Meanlivescale and Meandeadscale
# Drop twigs with 0 scale

scalecount2_para_july <- 
  scalecount2_july %>%
  filter(Meanlivescale > 0 | Meandeadscale > 0) %>%
  drop_na(Label)

scalecount2_para_nov <-
  scalecount2_nov %>%
  filter(Meanlivescale > 0 | Meandeadscale > 0) %>%
  drop_na(Label)

#### (2) Parasitism Models ####

##### July #####

para_mod_july <- glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                         data=scalecount2_para_july,
                         family=binomial)

simr_para_mod_july <- simulateResiduals(para_mod_july)
plot(simr_para_mod_july)

##### November #####

para_mod_nov <-  glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                         data=scalecount2_para_nov,
                         family=binomial)

simr_para_mod_nov <- simulateResiduals(para_mod_nov)
plot(simr_para_mod_nov)

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
