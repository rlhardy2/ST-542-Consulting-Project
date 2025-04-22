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
library(rcompanion)
library(performance)
library(ggpubr)

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

# Treatment survival means
tmeans2_july <- get_treatment_survival_means(scalecount2_july)
tmeans2_july$Collection<-"July"
tmeans2_nov <- get_treatment_survival_means(scalecount2_nov)
tmeans2_nov$Collection<-"Nov"
tmeans2_table <- dplyr::bind_rows(tmeans2_july, tmeans2_nov)

# Getting the data by tree
scalecount2_avg_july <- average_counts_across_twigs(scalecount2_july)
scalecount2_avg_nov <- average_counts_across_twigs(scalecount2_nov)

#### (2) Graphs - Exploratory ####

##### Distribution #####
get_hist_livescale(scalecount2_july, collection_date="July", study=2, 
                   labels=trt_labels2)
get_hist_livescale(scalecount2_nov, collection_date="Nov", study=2, 
                   labels=trt_labels2)

get_hist_all_trt(scalecount2_july, x_str="Meanlivescale", x_lab="Mean Live Scale",
                 title="Study 2 - Mean Live Scale, All Treatments, July",
                 labels=trt_labels2)

get_hist_all_trt(scalecount2_nov, x_str="Meanlivescale", x_lab="Mean Live Scale",
                 title="Study 2 - Mean Live Scale, All Treatments, Nov",
                 labels=trt_labels2)

##### Mean Live Scale #####
plot_means_by_collection(data=tmeans2_table, 
                         title="Study 2 - Mean Live EHS & Crypto Scale",
                         x_str="Treatment", y_str="Mean", labels=trt_labels2)


#### (3) Friedman Test + Nemenyi Test ####

##### July #####

# Friedman test only supports unreplicated complete block designs
# This needs to be averaged across all treatments and will be less precise

scalecount2_avg_across_block_trt <-
  average_counts_across_block_trt(scalecount2_avg_july)

# For some reason, this won't work without converting to a matrix
scalecount2_avg_block_trt_matrix <-
  as.matrix(scalecount2_avg_across_block_trt)

#### Friedman test for July -- mean live scale
# For some reason this only works if you use as.matrix lol
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount2_avg_block_trt_matrix)
friedman

# P-value of 0.5779 from Friedman test above is not significant!
# Therefore we do not need to conduct the Nemenyi (or other) pairwise test!

##### November #####

# This needs to be averaged across all treatments and will be less precise
scalecount2_avg_across_block_trt <-
  average_counts_across_block_trt(scalecount2_avg_nov)

# For some reason, this won't work without converting to a matrix
scalecount2_avg_block_trt_matrix <-
  as.matrix(scalecount2_avg_across_block_trt)

#### Friedman test for November -- mean live scale
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount2_avg_block_trt_matrix)
friedman

# P-value of 0.1562 from Friedman test above is not significant!
# Therefore we do not need to conduct the Nemenyi (or other) pairwise test!

#### (4) Scheirer–Ray–Hare Test ####

##### July #####

# It appears that we don't need an unreplicated complete block design here
# In order to use this test the observations need to be independent, so I am
# using the by tree data here

scheirer_july <- scheirerRayHare(Meanlivescale ~ Treatment | Block,
                                 data = scalecount2_avg_july)
scheirer_july

# Note that the p-value for Treatment is not significant
# The p-value for Treatment:Block interaction is significant

##### November #####

scheirer_nov <- scheirerRayHare(Meanlivescale ~ Treatment | Block,
                                 data = scalecount2_avg_nov)
scheirer_nov

# Note that the p-value for Treatment is significant
# The p-value for Treatment:Block interaction is significant

#### (5) Poisson & Negative Binomial Models ####

##### July #####

# Mixed Poisson, no zero inflation
pois_july2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount2_july, ziformula = ~0,
                     family = poisson)
# Residuals
simr_pois_july2 <- simulateResiduals(pois_july2)
plot(simr_pois_july2)

# Mixed NB model, no zero inflation
# Negative Binomial 2 (typical)
nb2_july2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount2_july, ziformula = ~0,
                    family = nbinom2)

simr_nb2_july2 <- simulateResiduals(nb2_july2)
plot(simr_nb2_july2)

#testCategorical(simr_nb2_july, scalecount2_july$Treatment) # gives error...

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_july2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount2_july, ziformula = ~0,
                    family = nbinom1)

simr_nb1_july2 <- simulateResiduals(nb1_july2)
plot(simr_nb1_july2)

#testCategorical(simr_nb1_july, scalecount2_july$Treatment) # gives error...

##### November #####

# Mixed Poisson, no zero inflation
pois_nov2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount2_nov, ziformula = ~0,
                    family = poisson)

simr_pois_nov2 <- simulateResiduals(pois_nov2)
plot(simr_pois_nov2)
check_overdispersion(pois_nov2)

# Mixed NB model, no zero inflation
# Using type 2 nbinom as is typical
nb2_nov2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount2_nov, ziformula = ~0,
                   family = nbinom2)

simr_nb2_nov2 <- simulateResiduals(nb2_nov2)
plot(simr_nb2_nov2)

#testCategorical(simr_nb2_nov, scalecount2_nov$Treatment) # gives error...

# Mixed NB1 model, no zero inflation
nb1_nov2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount2_nov, ziformula = ~0,
                   family = nbinom1)

simr_nb1_nov2 <- simulateResiduals(nb1_nov2)
plot(simr_nb1_nov2)
#testCategorical(simr_nb1_nov, scalecount2_nov$Treatment) # gives error...

#### (6) Zero-Inflated & Hurdle Models ####

##### July #####

# Zero-inflated negative binomial 2
zinb2_july2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount2_july, ziformula = ~1,
                      family = nbinom2)
# Got a warning about "model convergence problem"

simr_zinb2_july2 <- simulateResiduals(zinb2_july2)
plot(simr_zinb2_july2)

# Zero-inflated negative binomial 1
zinb1_july2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount2_july, ziformula = ~1,
                      family = nbinom1)
# Got a warning about "model convergence problem"

simr_zinb1_july2 <- simulateResiduals(zinb1_july2)
plot(simr_zinb1_july2)

##### November #####
# Diagnostics easier with glmmTMB than PSCL due to DHARMa compatibility

# Using negative binomial 2 - quadratic overdispersion
zinb2_nov2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                     data=scalecount2_nov, ziformula = ~1,
                     family = nbinom2)

simr_zinb2_nov2 <- simulateResiduals(zinb2_nov2)
plot(simr_zinb2_nov2)

# Using negative binomial 1 - linear overdispersion
zinb1_nov2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                     data=scalecount2_nov, ziformula = ~1,
                     family = nbinom1)

simr_zinb1_nov2 <- simulateResiduals(zinb1_nov2)
plot(simr_zinb1_nov2)

# Zero-inflated Poisson
zip_nov2 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                   data=scalecount2_nov, ziformula = ~1,
                   family = poisson)

simr_zip_nov2 <- simulateResiduals(zip_nov2)
plot(simr_zip_nov2)

# Hurdle negative binomial
hnbinom_nov2 <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                        data=scalecount2_nov,
                        ziformula = ~1,
                        family=truncated_nbinom2)
# Got a warning about "model convergence problem"

simr_hnb2_nov2 <- simulateResiduals(hnbinom_nov2)
plot(simr_hnb2_nov2)

# Hurdle Poisson
hpois_nov2 <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount2_nov,
                      ziformula = ~1,
                      family=truncated_poisson) # Got a weird warning message?

simr_hpois_nov2 <- simulateResiduals(hpois_nov2)
plot(simr_hpois_nov2)

#### (7) Analysis ####

##### July #####

###### Means ######

# Estimated marginal means
emm_nb1_july2 <- emmeans(nb1_july2, "Treatment")
emm_nb1_july2_orig_scale <- emmeans(nb1_july2, 
                                     "Treatment", type="response")
confint(emm_nb1_july2_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_july2, 
         sigma=sigma(nb1_july2), edf=df.residual(nb1_july2))

###### Pairwise Comparisons ######

# EMMs differ if arrows don't overlap
plot(emm_nb1_july2_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_july2_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_july2, adjust="BH"))
# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_nb1_july2), adjust="BH"))

##### November #####

###### Means ######

# Estimated marginal means
emm_nb1_nov2 <- emmeans(nb1_nov2, "Treatment")
emm_nb1_nov2_orig_scale <- emmeans(nb1_nov2, 
                                    "Treatment", type="response")
confint(emm_nb1_nov2_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_nov2, 
         sigma=sigma(nb1_nov2), edf=df.residual(nb1_nov2))

###### Pairwise Comparisons ######

# EMMs differ if arrows don't overlap
plot(emm_nb1_nov2_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_nov2_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_nov2, adjust="BH"))
# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_nb1_nov2), adjust="BH"))

#### (8) Graphs ####

##### July #####

pairs_july2_df <- as.data.frame(pairs(regrid(emm_nb1_july2), adjust="BH"))
confint_july2_df <- as.data.frame(confint(emm_nb1_july2_orig_scale))
y_positions <- seq(2.5, 7.5, by=1)

get_cis_marginal_means_plot(ci_df=confint_july2_df, pairs_df=pairs_july2_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="response", 
                            y_lab="Sum of Live Scale Across Twig (Estimated By Mean of Shoots)",
                            title="CIs of Estimated Treatment Means - July Live Scale Count, Study 2")

##### November #####

pairs_nov2_df <- as.data.frame(pairs(regrid(emm_nb1_nov2), adjust="BH"))
confint_nov2_df <- as.data.frame(confint(emm_nb1_nov2_orig_scale))
y_positions <- seq(3.5, 8.5, by=1)

# This looks quite strange lol, because the May acetamiprid trt has a mean of 0...
get_cis_marginal_means_plot(ci_df=confint_nov2_df, pairs_df=pairs_nov2_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="response", 
                            y_lab="Sum of Live Scale Across Twig (Estimated By Mean of Shoots)",
                            title="CIs of Estimated Treatment Means - November Live Scale Count, Study 2")

