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
library(ggpubr)

library(performance)

source("../graphing functions.r")
source("../Data Processing.r")

#### (1) Data Pre-processing ####
trt_labels <- c("May 2-3", "May 17", "May 28", "No Treatment")

# Read in file
scalecount <- read.csv(file="EHS count 2024 v7 (study 1).csv", strip.white=TRUE)
colnames(scalecount) <- c("Label", "County","Twigab","Date","Counter","Livescale1","Deadscale1",
                          "Livescale2","Deadscale2","Livescale3","Deadscale3",
                          "Prespara","Presfungus",
                          "Presscalenewgr","encarsia")

scalecount <- process_scalecount(scalecount)

# also change county to factor
scalecount$County <- as.factor(scalecount$County)

# Getting first count - July
scalecount_july <- subset(scalecount, grepl('July', Date)) # July data
scalecount_nov <- subset(scalecount, grepl('November', Date)) # November data

##### Specialized Encarsia pre-processing #####
# Drop counts where encarsia is NA
scalecountencar_july <-scalecount_july %>% drop_na(encarsia)
scalecountencar_nov <-scalecount_nov %>% drop_na(encarsia)

# Study 1's data set has encarsia count for each twig
# Combine by tree, as encarsia is by tree
scalecountencar_july <- average_counts_by_tree(scalecountencar_july)
scalecountencar_nov <- average_counts_by_tree(scalecountencar_nov)

# Treatment means
tmeans_encar_july <- scalecountencar_july %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(encarsia, na.rm = T),3),
  )
tmeans_encar_july$Collection<-"July"

tmeans_encar_nov <- scalecountencar_nov %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(encarsia, na.rm = T),3),
  )
tmeans_encar_nov$Collection<-"Nov"

encar_table_trt <- dplyr::bind_rows(tmeans_encar_july, 
                                    tmeans_encar_nov)

#### (2) Exploratory Graphs ####
get_hist_all_trt(scalecountencar_july, x_str="encarsia", x_lab="Encarsia Count",
                 title="Study 1 - Encarsia, All Treatments, July")
get_hist_all_trt(scalecountencar_nov, x_str="encarsia", x_lab="Encarsia Count",
                 title="Study 1 - Encarsia, All Treatments, Nov")

get_hist_encarsia(data=scalecountencar_july, collection_date="July",
                  study=1, labels=trt_labels)

get_hist_encarsia(data=scalecountencar_nov, collection_date="November",
                  study=1, labels=trt_labels)

plot_means_by_collection(data=encar_table_trt, 
                         title="Study 1 - Mean Encarsia Count Per Tree",
                         x_str="Treatment",
                         y_str="Mean")

#### (3) Encarsia models ####

# Since these are by tree, only need random intercept for Block

##### July #####

# Poisson model
encar_pois_july <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar_july, ziformula = ~0,
                          family = poisson)
simr_encar_pois_july <- simulateResiduals(encar_pois_july)
plot(simr_encar_pois_july)

# Negative binomial, linear overdispersion
encar_nb1_july <- glmmTMB(encarsia ~ Treatment + (1| Block),
                        data=scalecountencar_july, ziformula = ~0,
                        family = nbinom1)
simr_encar_nb1_july <- simulateResiduals(encar_nb1_july)
plot(simr_encar_nb1_july)

##### November #####

# Poisson model
encar_pois_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar_nov, ziformula = ~0,
                          family = poisson)
simr_encar_pois_nov <- simulateResiduals(encar_pois_nov)
plot(simr_encar_pois_nov)

# Negative binomial - type 1 (linearly overdispersed)
encar_nb1_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                        data=scalecountencar_nov, ziformula = ~0,
                        family = nbinom1)
# Fits pretty well
simr_encar_nb1_nov <- simulateResiduals(encar_nb1_nov)
plot(simr_encar_nb1_nov)
testCategorical(simr_encar_nb1_nov, scalecountencar_nov$Treatment)

# Negative binomial - type 2 (quadratically overdispersed)
encar_nb2_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar_nov, ziformula = ~0,
                         family = nbinom2)
simr_encar_nb2_nov <- simulateResiduals(encar_nb2_nov)
plot(simr_encar_nb2_nov)
testCategorical(simr_encar_nb2_nov, scalecountencar_nov$Treatment)

###### Zero Inflation Test ######

# negative binomial 2
encar_zinb2_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar_nov, ziformula = ~1,
                         family = nbinom2)
simr_encar_zinb2_nov <- simulateResiduals(encar_zinb2_nov)
plot(simr_encar_zinb2_nov)

# Poisson
encar_zip_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                           data=scalecountencar_nov, ziformula = ~1,
                           family = poisson)
simr_encar_zip_nov <- simulateResiduals(encar_zip_nov)
plot(simr_encar_zip_nov)

#### (4) Analysis ####

##### July #####
# Treatment means - July
emm_encar_july <- emmeans(encar_nb1_july, "Treatment")
emm_encar_july_orig_scale <- emmeans(encar_nb1_july, "Treatment", type="response")
confint(emm_encar_july_orig_scale)

# Treatment means comparisons - July
pairs_encar_july <- pairs(regrid(emm_encar_july), adjust="BH")
pairs_encar_july
confint(pairs_encar_july)


##### November #####
# Treatment means - Nov
emm_encar_nov <- emmeans(encar_nb1_nov, "Treatment")
emm_encar_nov_orig_scale <- emmeans(encar_nb1_nov, "Treatment", type="response")
confint(emm_encar_nov_orig_scale)

# Treatment means comparisons - Nov
pairs_encar_nov <- pairs(regrid(emm_encar_nov), adjust="BH")
pairs_encar_nov
confint(pairs_encar_nov)

#### (5) Graphs ####
##### July #####
pairs_encar_july_df <- as.data.frame(pairs(regrid(emm_encar_july), adjust="BH"))
confint_encar_july_df <- as.data.frame(confint(emm_encar_july_orig_scale))
y_positions <- seq(19, 25, by=1.2)

get_cis_marginal_means_plot(ci_df=confint_encar_july_df, pairs_df=pairs_encar_july_df, 
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="response", 
                            y_lab="Encarsia Count (Per Tree)",
                            title="CIs of Estimated Treatment Means - July Encarsia, Study 1")

##### November #####
pairs_encar_nov_df <- as.data.frame(pairs(regrid(emm_encar_nov), adjust="BH"))
confint_encar_nov_df <- as.data.frame(confint(emm_encar_nov_orig_scale))
y_positions <- seq(19, 25, by=1.2)

get_cis_marginal_means_plot(ci_df=confint_encar_nov_df, pairs_df=pairs_encar_nov_df, 
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="response", 
                            y_lab="Encarsia Count (Per Tree)",
                            title="CIs of Estimated Treatment Means - Nov Encarsia, Study 1")

