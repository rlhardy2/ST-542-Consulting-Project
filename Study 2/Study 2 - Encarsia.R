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

##### Specialized Encarsia pre-processing #####

# Drop counts where encarsia is NA
scalecountencar2_july <-scalecount2_july %>% drop_na(encarsia)
scalecountencar2_nov <-scalecount2_nov %>% drop_na(encarsia)

# Study 1's data set has encarsia count for each twig
# Combine by tree, as encarsia is by tree
scalecountencar2_july <- average_counts_across_twigs(scalecountencar2_july)
scalecountencar2_nov <- average_counts_across_twigs(scalecountencar2_nov)

#### (2) Encarsia models ####

##### July #####

trt_labels <- c("1 ?", "2 ?", "3 ?", "4 ?") # I don't think these are right lol

# Graph histogram of encarsia -- this needs to be double checked...
get_hist_all_trt(scalecountencar2_july, 
                 "encarsia", "Encarsia", 
                 "Study 2 Encarsia Count - July",
                 trt_labels)

# Poisson model
encar_pois_july <- glmmTMB(encarsia ~ Treatment + (1| Block),
                           data=scalecountencar2_july, ziformula = ~0,
                           family = poisson)

simr_encar_pois_july <- simulateResiduals(encar_pois_july)
plot(simr_encar_pois_july)

# Negative binomial, linear overdispersion
encar_nb1_july <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar2_july, ziformula = ~0,
                          family = nbinom1)

simr_encar_nb1_july <- simulateResiduals(encar_nb1_july)
plot(simr_encar_nb1_july)

##### November #####

# Graph histogram of encarsia -- this needs to be double checked...
get_hist_all_trt(scalecountencar2_nov, 
                 "encarsia", "Encarsia", 
                 "Study 2 Encarsia Count - Nov", trt_labels)

# Poisson model
# Since this is by tree, don't need label within block
encar_pois_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar2_nov, ziformula = ~0,
                          family = poisson)

simr_encar_pois_nov <- simulateResiduals(encar_pois_nov)
plot(simr_encar_pois_nov)

# Negative binomial - type 1 (linearly overdispersed)
encar_nb1_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar2_nov, ziformula = ~0,
                         family = nbinom1)

simr_encar_nb1_nov <- simulateResiduals(encar_nb1_nov)
plot(simr_encar_nb1_nov)
testCategorical(simr_encar_nb1_nov, scalecountencar2_nov$Treatment) # gives error...

# Negative binomial - type 2 (quadratically overdispersed)
encar_nb2_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar2_nov, ziformula = ~0,
                         family = nbinom2)

simr_encar_nb2_nov <- simulateResiduals(encar_nb2_nov)
plot(simr_encar_nb2_nov)
testCategorical(simr_encar_nb2_nov, scalecountencar2_nov$Treatment) # gives error...

# Check zero inflated just in case...
# No evidence of zero inflation
encar_zinb2_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                           data=scalecountencar2_nov, ziformula = ~1,
                           family = nbinom2)
# Got a warning about "model convergence problem"

simr_encar_zinb2_nov <- simulateResiduals(encar_zinb2_nov)
plot(simr_encar_zinb2_nov)

# Test zero inflated Poisson
encar_zip_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar2_nov, ziformula = ~1,
                         family = poisson)

simr_encar_zip_nov <- simulateResiduals(encar_zip_nov)
plot(simr_encar_zip_nov)

#### (3) Analysis ####

##### July #####

# Treatment means - July
emm_encar_july <- emmeans(encar_nb1_july, "Treatment")
emm_encar_july_orig_scale <- emmeans(encar_nb1_july, "Treatment", type="response")
confint(emm_encar_july_orig_scale)

# Treatment means comparisons - July
pairs_encar_july <- pairs(regrid(emm_encar_july), adjust="BH")
confint(pairs_encar_july)

##### November #####

# Treatment means - Nov
emm_encar_nov <- emmeans(encar_nb1_nov, "Treatment")
emm_encar_nov_orig_scale <- emmeans(encar_nb1_nov, "Treatment", type="response")
confint(emm_encar_nov_orig_scale)

# Treatment means comparisons - Nov
pairs_encar_nov <- pairs(regrid(emm_encar_nov), adjust="BH")
confint(pairs_encar_nov)
