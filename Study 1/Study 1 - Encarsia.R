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
scalecount_nov <- subset(scalecount, grepl('julyember', Date)) # November data

#### (2) Encarsia models ####

# drop counts where NA
scalecountencar_july <-scalecount_july %>% drop_na(encarsia)
scalecountencar_nov <-scalecount_nov %>% drop_na(encarsia)

# Study 1's data set has encarsia count for each twig
# Combine by tree, as encarsia is by tree
scalecountencar_july <- average_counts_across_twigs(scalecountencar_july)
scalecountencar_nov <- average_counts_across_twigs(scalecountencar_nov)

##### July #####
# Graph histogram of encarsia
get_hist_all_trt(scalecountencar_july, 
                 "encarsia", "Encarsia", 
                 "Study 1 Encarsia Count - July", trt_labels)

encar_pois_july <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar_july, ziformula = ~0,
                          family = poisson)
# Looks like there's more data at extremes
simr_encar_pois_july <- simulateResiduals(encar_pois_july)

# Negative binomial, linear overdispersion
encar_nb_july <- glmmTMB(encarsia ~ Treatment + (1| Block),
                        data=scalecountencar_july, ziformula = ~0,
                        family = nbinom1)
# Better looking residuals than Poisson
simr_encar_nb_july <- simulateResiduals(encar_nb_july)
plot(simr_encar_nb_july)

##### November #####
get_hist_all_trt(scalecountencar_nov, 
                 "encarsia", "Encarsia", 
                 "Study 1 Encarsia Count - Nov", trt_labels)

# Since this is by tree, don't need label within block
encar_pois_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar_nov, ziformula = ~0,
                          family = poisson)
# Some deviation issues...looks like there's more data at extremes
simr_encar_pois_nov <- simulateResiduals(encar_pois_nov)

# Negative binomial - type 1 (linearly overdispersed)
encar_nb_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                        data=scalecountencar_nov, ziformula = ~0,
                        family = nbinom1)
# Fits pretty well
simr_encar_nb_nov <- simulateResiduals(encar_nb_nov)

# Negative binomial - type 2 (quadratically overdispersed)
encar_nb2_nov <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar_nov, ziformula = ~0,
                         family = nbinom2)
# Not as good as nb type 1
simr_encar_nb2_nov <- simulateResiduals(encar_nb2_nov)

# Best one looks to be type 1 NB

#### (3) Analysis ####

##### July #####
# Treatment means - July
emm_encar_july <- emmeans(encar_nb_july, "Treatment")
emm_encar_july_orig_scale <- emmeans(encar_nb_july, "Treatment", type="response")
confint(emm_encar_july_orig_scale)

# Treatment means comparisons - july
pairs_encar_july <- pairs(regrid(emm_encar_july), adjust="BH")
confint(pairs_encar_july)


##### November #####
# Treatment means - nov
emm_encar_nov <- emmeans(encar_nb_nov, "Treatment")
emm_encar_nov_orig_scale <- emmeans(encar_nb_nov, "Treatment", type="response")

# Treatment means comparisons - nov
pairs_encar_nov <- pairs(regrid(emm_encar_nov), adjust="BH")
confint(pairs_encar_nov)
