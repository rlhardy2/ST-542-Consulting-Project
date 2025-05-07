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
library(emmeans)
library(DHARMa)
library(glmmTMB)
library(performance)

source("../graphing functions.r")
source("../Data Processing.r")

# Reformatted version of Study 1, for scale count analysis

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


##### New Processing for Per-Shoot Modeling ######

# Get long scalecount for July
scalecount_july_long <- get_per_shoot_scalecount(scalecount_july)

# Get only live scalecount from long, without NAs
scalecount_july_long_live <- get_per_shoot_scalecount_live(scalecount_july_long)


# Get long scalecount for November
scalecount_nov_long <- get_per_shoot_scalecount(scalecount_nov)

# Get only live scalecount from long, without NAs
scalecount_nov_long_live <- get_per_shoot_scalecount_live(scalecount_nov_long)


#### Models - Non-Zero Inflated ####

# Mixed Poisson, no zero inflation
pois_nov <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                     data=scalecount_nov_long_live, ziformula = ~0,
                     family = poisson)
simr_pois_nov <- simulateResiduals(pois_nov)
plot(simr_pois_nov, title="Poisson")

# Mixed NB model, no zero inflation
# Negative Binomial 2 (typical) - quadratic overdispersion
nb2_nov <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount_nov_long_live, ziformula = ~0,
                    family = nbinom2)
simr_nb2_nov <- simulateResiduals(nb2_nov)
plot(simr_nb2_nov, title="Negative Binomial 2")

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_nov <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount_nov_long_live, ziformula = ~0,
                    family = nbinom1)
simr_nb1_nov <- simulateResiduals(nb1_nov, title="Negative Binomial 1")
plot(simr_nb1_nov)

#### Models - Zero-Inflated ####

zinb2_nov <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                     data=scalecount_nov_long_live, ziformula = ~1,
                     family = nbinom2)
simr_zinb2_nov <- simulateResiduals(zinb2_nov)
plot(simr_zinb2_nov)

# Zero-inflated Poisson
zip_nov <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                   data=scalecount_nov_long_live, ziformula = ~1,
                   family = poisson)
simr_zip_nov <- simulateResiduals(zip_nov)
plot(simr_zip_nov)

# Hurdle negative binomial
hnbinom_nov <-  glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                        data=scalecount_nov_long_live,
                        ziformula = ~1,
                        family=truncated_nbinom2)
simr_hnb2_nov <- simulateResiduals(hnbinom_nov)
plot(simr_hnb2_nov)

# Hurdle Poisson
hpois_nov <-  glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                      data=scalecount_nov_long_live,
                      ziformula = ~1,
                      family=truncated_poisson)
simr_hpois_nov <- simulateResiduals(hpois_nov)
plot(simr_hpois_nov)


#### Analysis ####

##### Means #####
# Estimated marginal means
emm_nb1_nov <- emmeans(nb1_nov, "Treatment")
emm_nb1_nov_orig_scale <- emmeans(nb1_nov, 
                                  "Treatment", type="response")

# Per-shoot scale
confint(emm_nb1_nov_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_nov, 
         sigma=sigma(nb1_nov), edf=df.residual(nb1_nov))

##### Pairwise Comparisons #####
plot(emm_nb1_nov_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_nov_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_nov, adjust="BH"))


# CI for pairwise comparison on original scale (per-shoot)
confint(pairs(regrid(emm_nb1_nov), adjust="BH"))

# Pairwise comparisons per-shoot (p-values)
pairs(regrid(emm_nb1_nov), adjust="BH")


#### Graphs ####
pairs_nov_df <- as.data.frame(pairs(regrid(emm_nb1_nov), adjust="BH"))
confint_nov_df <- as.data.frame(confint(emm_nb1_nov_orig_scale))
y_positions <- seq(3, 9, by=1.2)

get_cis_marginal_means_plot(ci_df=confint_nov_df, pairs_df=pairs_nov_df, 
                            y_positions=y_positions, trt_labels=trt_labels,
                            y_str="response", 
                            y_lab="Scale Count (Per Shoot)",
                            title="CIs of Estimated Treatment Means - November Live Scale Count, Study 1")

