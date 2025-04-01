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
tmeans2_nov <- get_treatment_survival_means(scalecount2_nov)

#### (2) Nonparametric Analysis ####

##### July #####

# Friedman test
# Friedman test only supports unreplicated complete block designs
# This needs to be averaged across all treatments and will be less precise
scalecount2_avg_july <- average_counts_across_twigs(scalecount2_july)
scalecount2_avg_across_block_trt <- average_counts_across_block_trt(scalecount2_avg_july)

# Note that after grouping data by treatment-block combination, there are only 16 observations!
# Is that going to be a problem...??

# For some reason, this won't work without converting to a matrix
scalecount2_avg_block_trt_matrix <- as.matrix(scalecount2_avg_across_block_trt)

#### Friedman test for July -- mean live scale
# For some reason this only works if you use as.matrix lol
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount2_avg_block_trt_matrix)
friedman

# P-value of 0.5456 from Friedman test above is not significant!
# Therefore we do not need to conduct the Nemenyi (or other) pairwise test!

##### November #####

# Friedman test
# This needs to be averaged across all treatments and will be less precise
scalecount2_avg_nov <- average_counts_across_twigs(scalecount2_nov)
scalecount2_avg_across_block_trt <- average_counts_across_block_trt(scalecount2_avg_nov)

# Again, note that the data grouped by treatment-block combination only has 16 observations

# For some reason, this won't work without converting to a matrix
scalecount2_avg_block_trt_matrix <- as.matrix(scalecount2_avg_across_block_trt)

#### Friedman test for November -- mean live scale
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount2_avg_block_trt_matrix)
friedman

# P-value of 0.1562 from Friedman test above is not significant!
# Therefore we do not need to conduct the Nemenyi (or other) pairwise test!

#### (3) Poisson & Negative Binomial Models ####

##### July #####

# Mixed Poisson, no zero inflation
pois_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount2_july, ziformula = ~0,
                     family = poisson)
# Residuals
simr_pois_july <- simulateResiduals(pois_july)
plot(simr_pois_july)

# Mixed NB model, no zero inflation
# Negative Binomial 2 (typical)
nb2_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount2_july, ziformula = ~0,
                    family = nbinom2)

simr_nb2_july <- simulateResiduals(nb2_july)
plot(simr_nb2_july)

#testCategorical(simr_nb2_july, scalecount2_july$Treatment) # gives error...

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount2_july, ziformula = ~0,
                    family = nbinom1)

simr_nb1_july <- simulateResiduals(nb1_july)
plot(simr_nb1_july)

#testCategorical(simr_nb1_july, scalecount2_july$Treatment) # gives error...

##### November #####

# Mixed Poisson, no zero inflation
pois_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount2_nov, ziformula = ~0,
                    family = poisson)

simr_pois_nov <- simulateResiduals(pois_nov)
plot(simr_pois_nov)
check_overdispersion(pois_nov)

# Mixed NB model, no zero inflation
# Using type 2 nbinom as is typical
nb2_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount2_nov, ziformula = ~0,
                   family = nbinom2)

simr_nb2_nov <- simulateResiduals(nb2_nov)
plot(simr_nb2_nov)

#testCategorical(simr_nb2_nov, scalecount2_nov$Treatment) # gives error...

# Mixed NB1 model, no zero inflation
nb1_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount2_nov, ziformula = ~0,
                   family = nbinom1)

simr_nb1_nov <- simulateResiduals(nb1_nov)
plot(simr_nb1_nov)
#testCategorical(simr_nb1_nov, scalecount2_nov$Treatment) # gives error...

#### (5) Zero-Inflated & Hurdle Models ####

##### July #####

# Zero-inflated negative binomial 2
zinb2_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount2_july, ziformula = ~1,
                      family = nbinom2)
# Got a warning about "model convergence problem"

simr_zinb2_july <- simulateResiduals(zinb2_july)
plot(simr_zinb2_july)

# Zero-inflated negative binomial 1
zinb1_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount2_july, ziformula = ~1,
                      family = nbinom1)
# Got a warning about "model convergence problem"

simr_zinb1_july <- simulateResiduals(zinb1_july)
plot(simr_zinb1_july)

##### November #####
# Diagnostics easier with glmmTMB than PSCL due to DHARMa compatibility

# Using negative binomial 2 - quadratic overdispersion
zinb2_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                     data=scalecount2_nov, ziformula = ~1,
                     family = nbinom2)

simr_zinb2_nov <- simulateResiduals(zinb2_nov)
plot(simr_zinb2_nov)

# Using negative binomial 1 - linear overdispersion
zinb1_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                     data=scalecount2_nov, ziformula = ~1,
                     family = nbinom1)

simr_zinb1_nov <- simulateResiduals(zinb1_nov)
plot(simr_zinb1_nov)

# Zero-inflated Poisson
zip_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                   data=scalecount2_nov, ziformula = ~1,
                   family = poisson)

simr_zip_nov <- simulateResiduals(zip_nov)
plot(simr_zip_nov)

# Hurdle negative binomial
hnbinom_nov <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                        data=scalecount2_nov,
                        ziformula = ~1,
                        family=truncated_nbinom2)
# Got a warning about "model convergence problem"

simr_hnb2_nov <- simulateResiduals(hnbinom_nov)
plot(simr_hnb2_nov)

# Hurdle Poisson
hpois_nov <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount2_nov,
                      ziformula = ~1,
                      family=truncated_poisson) # Got a weird warning message?

simr_hpois_nov <- simulateResiduals(hpois_nov)
plot(simr_hpois_nov)

#### (6) Analysis ####

##### July #####

###### Means ######

# Estimated marginal means
emm_zinb2_july <- emmeans(zinb2_july, "Treatment")
emm_zinb2_july_orig_scale <- emmeans(zinb2_july, 
                                     "Treatment", type="response")
confint(emm_zinb2_july_orig_scale)
# effect size - Cohen's d
eff_size(emm_zinb2_july, 
         sigma=sigma(zinb2_july), edf=df.residual(zinb2_july))

###### Pairwise Comparisons ######

# EMMs differ if arrows don't overlap
plot(emm_zinb2_july_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_zinb2_july_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_zinb2_july, adjust="BH"))
# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_zinb2_july), adjust="BH"))

##### November #####

###### Means ######

# Estimated marginal means
emm_zinb2_nov <- emmeans(zinb2_nov, "Treatment")
emm_zinb2_nov_orig_scale <- emmeans(zinb2_nov, 
                                    "Treatment", type="response")
confint(emm_zinb2_nov_orig_scale)
# effect size - Cohen's d
eff_size(emm_zinb2_nov, 
         sigma=sigma(zinb2_nov), edf=df.residual(zinb2_nov))

###### Pairwise Comparisons ######

# EMMs differ if arrows don't overlap
plot(emm_zinb2_nov_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_zinb2_nov_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_zinb2_nov, adjust="BH"))
# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_zinb2_nov), adjust="BH"))

#### (7) Graphs ####

##### July #####

# Making the confidence intervals a data frame
confint_july <- as.data.frame(confint(emm_zinb2_july_orig_scale))

pairs_july_df<- as.data.frame(pairs(regrid(emm_zinb2_july), adjust="BH"))

# Extracts treatments from pairs for graphing
# Must be called group1 and group2 for stat_pvalue_manual
pairs_july_df$group1 <- substring(pairs_july_df$contrast, first=10, last=10)
pairs_july_df$group2 <- substring(pairs_july_df$contrast, first=23, last=23)
pairs_july_df$p.value_round <- round(pairs_july_df$p.value, digits=3)
# y position of each pairwise comparison bar
pairs_july_df$y.position <- seq(1.5, 2.5, by=.2)

# Plot of the CIs for the estimated marginal means for each treatment
ggplot(confint_july, aes(x = Treatment, y = response)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(title = "CIs of the Estimated Marginal Means by Treatment",
       x = "Treatment",
       y = "Sumlivescale_from_mean") +
  stat_pvalue_manual(
    data = pairs_july_df, label = "p.value_round",
    xmin = "group1", xmax = "group2",
    y.position = "y.position"
  )


# Making the confidence intervals a data frame (for the pairs)
confint_pairs_july <- as.data.frame(confint(pairs(regrid(emm_zinb2_nov), adjust="BH")))


# Plot of the CIs for the estimated marginal means for each pair
ggplot(confint_pairs_july, aes(x = contrast, y = estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(title = "CIs of the Estimated Marginal Means by Treatment",
       x = "Treatment",
       y = "Sumlivescale_from_mean") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##### November #####
