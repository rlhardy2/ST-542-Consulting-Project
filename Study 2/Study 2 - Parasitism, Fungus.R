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
library(glmmTMB)

source("../graphing functions.r")
source("../Data Processing.r")

# Reformatted version of Study 2.

#### (1) Data Pre-processing ####

trt_labels2 <- c("Pyri/May", "Ace/May", "Pyri/May + Ace/June", "No Treatment")

# Reading in the data
# Not sure why but there is an extra column of missing values at the end of the data set
# This is removed below along with the 'notes' column
scalecount2 <- read.csv(file = "distance numbers v7 (study 2).csv", strip.white=TRUE)
colnames(scalecount2) <- c("Label","Type", "Twigab","Date","Livescale1","Deadscale1","Livescale2",
                           "Deadscale2","Livescale3","Deadscale3","Prespara",
                           "Presscalenewgr","Presfungus","encarsia","notes","extra")

# Remove notes column and random extra column and Presscalenewgr (unused)
scalecount2 <- subset(scalecount2, select = -c(notes, extra, Presscalenewgr))

# Drop uncollected labels
scalecount2 <- scalecount2 %>% drop_na(Livescale1)

# Extract treatment - needed for merging step
# Study 2-only processing - Add EHS and cryptomeria scale together
scalecount2$Treatment<- extract_treatment(scalecount2)

scalecount2$Prespara <- replace_strings_with_binary(scalecount2$Prespara)
scalecount2$Presfungus <- replace_strings_with_binary(scalecount2$Presfungus)
scalecount2$Prespara <- as.numeric(scalecount2$Prespara)
scalecount2$Presfungus <- as.numeric(scalecount2$Presfungus)

# Expand scalecount2 to "wide" format, creating vars for Live/Dead scale
# and presence to add together more easily
scalecount2_wide  <- 
  scalecount2 %>%
  pivot_wider(names_from=Type, values_from=c(Livescale1, Livescale2, Livescale3,
                                             Deadscale1, Deadscale2, Deadscale3,
                                             Prespara, Presfungus))

# Add variables together
scalecount2_wide <- create_live_deadscale_pres_vars_wide(scalecount2_wide)

# Process data set
scalecount2 <- process_scalecount(scalecount2_wide)

# July and November data
scalecount2_july <- subset(scalecount2, grepl('July', Date))
scalecount2_nov <- subset(scalecount2, grepl('November', Date))

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
confint_emm_para_july_orig_scale <- confint(emm_para_july, type="response")
confint_emm_para_july_orig_scale

##### November #####

emm_para_nov <- emmeans(para_mod_nov, "Treatment")
pairs(emm_para_nov, type="response", adjust="BH")
nov_para_means_comp <- pairs(regrid(emm_para_nov), adjust="BH")
# Confint for pairwise comparison
confint(nov_para_means_comp)
# Confint for means
confint_emm_para_nov_orig_scale <- confint(emm_para_nov, type="response")
confint_emm_para_nov_orig_scale

#### Parasitism Graphs ####

##### July #####

pairs_para_july_df <- as.data.frame(pairs(regrid(emm_para_july), adjust="BH"))
confint_para_july_df <- as.data.frame(confint_emm_para_july_orig_scale)
y_positions <- seq(.8, 1.3, by=.1)
get_cis_marginal_means_plot(ci_df=confint_para_july_df, pairs_df=pairs_para_july_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - July Parasitism, Study 2")

##### November #####

pairs_para_nov_df <- as.data.frame(pairs(regrid(emm_para_nov), adjust="BH"))
confint_para_nov_df <- as.data.frame(confint_emm_para_nov_orig_scale)
y_positions <- seq(1, 1.5, by=.1)
get_cis_marginal_means_plot(ci_df=confint_para_nov_df, pairs_df=pairs_para_nov_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - Nov Parasitism, Study 2")

#### (4) Fungus models ####

##### July #####

fung_mod_july <- glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                         data=scalecount2_para_july,
                         family=binomial)

simr_fung_mod_july <- simulateResiduals(fung_mod_july)
plot(simr_fung_mod_july)

##### November #####

fung_mod_nov <-  glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                         data=scalecount2_para_nov,
                         family=binomial)
# Warning about "model convergence problem"

simr_fung_mod_nov <- simulateResiduals(fung_mod_nov)
plot(simr_fung_mod_nov)

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

