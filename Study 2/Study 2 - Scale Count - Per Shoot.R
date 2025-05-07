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
# Preferred so that we can tie scale count to parasitism/fungus properly
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


##### New Processing for Per-Shoot Modeling ######

# Get long scalecount for July
scalecount2_july_long <- get_per_shoot_scalecount(scalecount2_july)

# Get only live scalecount from long, without NAs
scalecount2_july_long_live <- get_per_shoot_scalecount_live(scalecount2_july_long)


# Get long scalecount for November
scalecount2_nov_long <- get_per_shoot_scalecount(scalecount2_nov)

# Get only live scalecount from long, without NAs
scalecount2_nov_long_live <- get_per_shoot_scalecount_live(scalecount2_nov_long)


#### Model ####

##### July #####
pois_july2 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                     data=scalecount2_july_long_live, ziformula = ~0,
                     family = poisson)
simr_pois_july2 <- simulateResiduals(pois_july2)
plot(simr_pois_july2, title="Poisson")

# Mixed NB model, no zero inflation
# Negative Binomial 2 (typical) - quadratic overdispersion
nb2_july2 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount2_july_long_live, ziformula = ~0,
                    family = nbinom2)
simr_nb2_july2 <- simulateResiduals(nb2_july2)
plot(simr_nb2_july2, title="Negative Binomial 2")

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_july2 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount2_july_long_live, ziformula = ~0,
                    family = nbinom1)
simr_nb1_july2 <- simulateResiduals(nb1_july2, title="Negative Binomial 1")
plot(simr_nb1_july2)


##### November #####
# Mixed Poisson, no zero inflation
pois_nov2 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount2_nov_long_live, ziformula = ~0,
                    family = poisson)
simr_pois_nov2 <- simulateResiduals(pois_nov2)
plot(simr_pois_nov2, title="Poisson")

# Mixed NB model, no zero inflation
# Negative Binomial 2 (typical) - quadratic overdispersion
nb2_nov2 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                   data=scalecount2_nov_long_live, ziformula = ~0,
                   family = nbinom2)
simr_nb2_nov2 <- simulateResiduals(nb2_nov2)
plot(simr_nb2_nov2, title="Negative Binomial 2")

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_nov2 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                   data=scalecount2_nov_long_live, ziformula = ~0,
                   family = nbinom1)
simr_nb1_nov2 <- simulateResiduals(nb1_nov2, title="Negative Binomial 1")
plot(simr_nb1_nov2)

#### Analysis - July ####

##### Means #####
# Estimated marginal means
emm_nb1_july2 <- emmeans(nb1_july2, "Treatment")
emm_nb1_july2_orig_scale <- emmeans(nb1_july2, 
                                   "Treatment", type="response")

# Per-shoot scale
confint(emm_nb1_july2_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_july2, 
         sigma=sigma(nb1_july2), edf=df.residual(nb1_july2))

##### Pairwise Comparisons #####
plot(emm_nb1_july2_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_july2_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_july2, adjust="BH"))

# CI for pairwise comparison on original scale (per-shoot)
confint(pairs(regrid(emm_nb1_july2), adjust="BH"))

# Pairwise comparisons per-shoot (p-values)
pairs(regrid(emm_nb1_july2), adjust="BH")


##### Graphs #####
pairs_july2_df <- as.data.frame(pairs(regrid(emm_nb1_july2), adjust="BH"))
confint_july2_df <- as.data.frame(confint(emm_nb1_july2_orig_scale))
y_positions <- seq(3, 9, by=1.2)

get_cis_marginal_means_plot(ci_df=confint_july2_df, pairs_df=pairs_july2_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="response", 
                            y_lab="Scale Count (Per Shoot)",
                            title="CIs of Estimated Treatment Means - July Live Scale Count, Study 2")


#### Analysis - November ####

##### Means #####
# Estimated marginal means
emm_nb1_nov2 <- emmeans(nb1_nov2, "Treatment")
emm_nb1_nov2_orig_scale <- emmeans(nb1_nov2, 
                                  "Treatment", type="response")

# Per-shoot scale
confint(emm_nb1_nov2_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_nov2, 
         sigma=sigma(nb1_nov2), edf=df.residual(nb1_nov2))

##### Pairwise Comparisons #####
plot(emm_nb1_nov2_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_nov2_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_nov2, adjust="BH"))

# CI for pairwise comparison on original scale (per-shoot)
confint(pairs(regrid(emm_nb1_nov2), adjust="BH"))

# Pairwise comparisons per-shoot (p-values)
pairs(regrid(emm_nb1_nov2), adjust="BH")


##### Graphs #####
pairs_nov2_df <- as.data.frame(pairs(regrid(emm_nb1_nov2), adjust="BH"))
confint_nov2_df <- as.data.frame(confint(emm_nb1_nov2_orig_scale))
y_positions <- seq(3, 9, by=1.2)

get_cis_marginal_means_plot(ci_df=confint_nov2_df, pairs_df=pairs_nov2_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="response", 
                            y_lab="Scale Count (Per Shoot)",
                            title="CIs of Estimated Treatment Means - November Live Scale Count, Study 2")

