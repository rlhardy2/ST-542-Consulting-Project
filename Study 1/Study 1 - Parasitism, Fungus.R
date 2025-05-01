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
library(ggpubr)

library(DHARMa)
library(glmmTMB)

library(performance)

source("../graphing functions.r")
source("../Data Processing.r")

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

##### Specialized parasitism and fungus preprocessing #####
# Drop twigs with 0 scale

# Should look at Meanlivescale and Meandeadscale
scalecount_para_july <- 
  scalecount_july %>%
  filter(Meanlivescale > 0 | Meandeadscale > 0) %>%
  drop_na(Label)

scalecount_para_nov <-
  scalecount_nov %>%
  filter(Meanlivescale > 0 | Meandeadscale > 0) %>%
  drop_na(Label)

tmeans_july_para <- scalecount_para_july %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes_para = mean(Prespara == 1),
    percent_yes_fung = mean(Presfungus == 1),
    percent_yes_para_100 = round(percent_yes_para*100),
    percent_yes_fung_100 = round(percent_yes_fung*100)
  )
tmeans_july_para$Collection<-"July"

tmeans_nov_para <- scalecount_para_nov %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes_para = mean(Prespara == 1),
    percent_yes_fung = mean(Presfungus == 1),
    percent_yes_para_100 = round(percent_yes_para*100),
    percent_yes_fung_100 = round(percent_yes_fung*100)
  )
tmeans_nov_para$Collection<-"Nov"

para_table <- dplyr::bind_rows(tmeans_july_para, tmeans_nov_para)

#### (2) Exploratory Analysis ####

##### Parasitism #####
plot_means_by_collection(data=para_table, 
                         title="Study 1 - Mean Presence Percentage of Parasitism",
                         x_str="Treatment",
                         y_str="percent_yes_para_100",
                         labels=trt_labels,
                         y_lab="Percent")

##### Fungus #####
plot_means_by_collection(data=para_table, 
                         title="Study 1 - Mean Presence Percentage of Fungus",
                         x_str="Treatment",
                         y_str="percent_yes_fung_100",
                         labels=trt_labels,
                         y_lab="Percent")


#### (3) Parasitism models ####

# Use mixed models to account for random block and tree effects

##### July #####
para_mod_july <- glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                         data=scalecount_para_july,
                         family=binomial)
simr_para_mod_july <- simulateResiduals(para_mod_july)
plot(simr_para_mod_july)

##### November #####
para_mod_nov <-  glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                        data=scalecount_para_nov,
                        family=binomial)
simr_para_mod_nov <- simulateResiduals(para_mod_nov)
plot(simr_para_mod_nov)

#### (4) Parasitism Analysis ####

##### July #####
emm_para_july <- emmeans(para_mod_july, "Treatment")
pairs(emm_para_july, type="response", adjust="BH")
july_para_means_comp <- pairs(regrid(emm_para_july), adjust="BH")
july_para_means_comp
# CI for pairwise comp
confint(july_para_means_comp)
# CI for means
confint_emm_para_july_orig_scale <- confint(emm_para_july, type="response")
confint_emm_para_july_orig_scale

##### November #####
emm_para_nov <- emmeans(para_mod_nov, "Treatment")
pairs(emm_para_nov, type="response", adjust="BH")
nov_para_means_comp <- pairs(regrid(emm_para_nov), adjust="BH")
nov_para_means_comp
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
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - July Parasitism, Study 1")

##### November #####
pairs_para_nov_df <- as.data.frame(pairs(regrid(emm_para_nov), adjust="BH"))
confint_para_nov_df <- as.data.frame(confint_emm_para_nov_orig_scale)
y_positions <- seq(1, 1.5, by=.1)
get_cis_marginal_means_plot(ci_df=confint_para_nov_df, pairs_df=pairs_para_nov_df, 
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - Nov Parasitism, Study 1")

#### (5) Fungus models ####
##### July #####
fung_mod_july <- glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                         data=scalecount_para_july,
                         family=binomial)
simr_fung_mod_july <- simulateResiduals(fung_mod_july)
plot(simr_fung_mod_july)

##### November #####
fung_mod_nov <-  glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                         data=scalecount_para_nov,
                         family=binomial)
simr_fung_mod_nov <- simulateResiduals(fung_mod_nov)
plot(simr_fung_mod_nov)

#### (6) Fungus Analysis ####

##### July #####
emm_fung_july <- emmeans(fung_mod_july, "Treatment")
pairs(emm_fung_july, type="response", adjust="BH")
july_fung_means_comp <- pairs(regrid(emm_fung_july), adjust="BH")
july_fung_means_comp
# CI for pairwise comp
confint(july_fung_means_comp)
# CI for means
confint_emm_fung_july_orig_scale <- confint(emm_fung_july, type="response")
confint_emm_fung_july_orig_scale

##### November #####
emm_fung_nov <- emmeans(fung_mod_nov, "Treatment")
pairs(emm_fung_nov, type="response", adjust="BH")
nov_fung_means_comp <- pairs(regrid(emm_fung_nov), adjust="BH")
nov_fung_means_comp
# Confint for pairwise comparison
confint(nov_fung_means_comp)
# Confint for means
confint_emm_fung_nov_orig_scale <- confint(emm_fung_nov, type="response")
confint_emm_fung_nov_orig_scale

#### (7) Fungus Graphs ####

##### July #####
pairs_fung_july_df <- as.data.frame(pairs(regrid(emm_fung_july), adjust="BH"))
confint_fung_july_df <- as.data.frame(confint_emm_fung_july_orig_scale)
y_positions <- seq(.25, .75, by=.1)

get_cis_marginal_means_plot(ci_df=confint_fung_july_df, pairs_df=pairs_fung_july_df, 
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - July Fungus, Study 1")

##### November #####
pairs_fung_nov_df <- as.data.frame(pairs(regrid(emm_fung_nov), adjust="BH"))
confint_fung_nov_df <- as.data.frame(confint_emm_fung_nov_orig_scale)
y_positions <- seq(.35, .85, by=.1)

get_cis_marginal_means_plot(ci_df=confint_fung_nov_df, pairs_df=pairs_fung_nov_df, 
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - Nov Fungus, Study 1")
