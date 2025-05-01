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

##### Specialized Encarsia pre-processing #####

# Drop counts where encarsia is NA
scalecountencar2_july <-scalecount2_july %>% drop_na(encarsia)
scalecountencar2_nov <-scalecount2_nov %>% drop_na(encarsia)

# Study 1's data set has encarsia count for each twig
# Combine by tree, as encarsia is by tree
scalecountencar2_july <- average_counts_by_tree(scalecountencar2_july)
scalecountencar2_nov <- average_counts_by_tree(scalecountencar2_nov)

# July and November have same data?
scalecountencar2_july$encarsia
scalecountencar2_nov$encarsia

# Encarsia treatment means
tmeans_encar_july <- scalecountencar2_july %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(encarsia, na.rm = T),3)
  )


#### (2) Exploratory Graphs ####
get_hist_all_trt(scalecountencar2_july, x_str="encarsia", x_lab="Encarsia Count",
                 title="Study 2 - Encarsia, All Treatments, July")

get_hist_encarsia(data=scalecountencar2_july, collection_date="July",
                  study=2, labels=trt_labels2)

plot_means(data=tmeans_encar_july, 
                         title="Study 2 - Mean Encarsia Count Per Tree",
                         x_str="Treatment",
                         y_str="Mean", labels=trt_labels2)

#### (3) Encarsia models ####

##### July #####

# Poisson model
encar_pois_july2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                           data=scalecountencar2_july, ziformula = ~0,
                           family = poisson)

simr_encar_pois_july2 <- simulateResiduals(encar_pois_july2)
plot(simr_encar_pois_july2)

# Negative binomial, linear overdispersion
encar_nb1_july2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar2_july, ziformula = ~0,
                          family = nbinom1)

simr_encar_nb1_july2 <- simulateResiduals(encar_nb1_july2)
plot(simr_encar_nb1_july2)

##### November #####

# Poisson model
# Since this is by tree, don't need label within block
encar_pois_nov2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                          data=scalecountencar2_nov, ziformula = ~0,
                          family = poisson)

simr_encar_pois_nov2 <- simulateResiduals(encar_pois_nov2)
plot(simr_encar_pois_nov2)

# Negative binomial - type 1 (linearly overdispersed)
encar_nb1_nov2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar2_nov, ziformula = ~0,
                         family = nbinom1)

simr_encar_nb1_nov2 <- simulateResiduals(encar_nb1_nov2)
plot(simr_encar_nb1_nov2)
testCategorical(simr_encar_nb1_nov2, scalecountencar2_nov$Treatment)

# Negative binomial - type 2 (quadratically overdispersed)
encar_nb2_nov2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar2_nov, ziformula = ~0,
                         family = nbinom2)

simr_encar_nb2_nov2 <- simulateResiduals(encar_nb2_nov2)
plot(simr_encar_nb2_nov2)
testCategorical(simr_encar_nb2_nov2, scalecountencar2_nov$Treatment) # gives error...

###### Zero-Inflated ######

# Negative binomial 2
encar_zinb2_nov2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                           data=scalecountencar2_nov, ziformula = ~1,
                           family = nbinom2)
simr_encar_zinb2_nov2 <- simulateResiduals(encar_zinb2_nov2)
plot(simr_encar_zinb2_nov2)

# Poisson
encar_zip_nov2 <- glmmTMB(encarsia ~ Treatment + (1| Block),
                         data=scalecountencar2_nov, ziformula = ~1,
                         family = poisson)

simr_encar_zip_nov2 <- simulateResiduals(encar_zip_nov2)
plot(simr_encar_zip_nov2)

#### (4) Analysis ####

##### July #####

# Treatment means - July
emm_encar_july2 <- emmeans(encar_nb1_july2, "Treatment")
emm_encar_july2_orig_scale <- emmeans(encar_nb1_july2, "Treatment", type="response")
confint(emm_encar_july2_orig_scale)

# Treatment means comparisons - july2
pairs_encar_july2 <- pairs(regrid(emm_encar_july2), adjust="BH")
pairs_encar_july2
confint(pairs_encar_july2)

##### November #####

# Treatment means - nov2
emm_encar_nov2 <- emmeans(encar_nb1_nov2, "Treatment")
emm_encar_nov2_orig_scale <- emmeans(encar_nb1_nov2, "Treatment", type="response")
confint(emm_encar_nov2_orig_scale)

# Treatment means comparisons - nov2
pairs_encar_nov2 <- pairs(regrid(emm_encar_nov2), adjust="BH")
pairs_encar_nov2
confint(pairs_encar_nov2)


#### (5) Graphs ####
##### July #####
pairs_encar_july2_df <- as.data.frame(pairs(regrid(emm_encar_july2), adjust="BH"))
confint_encar_july2_df <- as.data.frame(confint(emm_encar_july2_orig_scale))
y_positions <- seq(4, 9, by=1)

get_cis_marginal_means_plot(ci_df=confint_encar_july2_df, pairs_df=pairs_encar_july2_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="response", 
                            y_lab="Encarsia Count (Per Tree)",
                            title="CIs of Estimated Treatment Means - July Encarsia Count, Study 2")

##### November #####
pairs_encar_nov2_df <- as.data.frame(pairs(regrid(emm_encar_nov2), adjust="BH"))
confint_encar_nov2_df <- as.data.frame(confint(emm_encar_nov2_orig_scale))
y_positions <- seq(4, 9, by=1)

get_cis_marginal_means_plot(ci_df=confint_encar_nov2_df, pairs_df=pairs_encar_nov2_df, 
                            y_positions=y_positions, trt_labels=trt_labels2,
                            y_str="response", 
                            y_lab="Encarsia Count (Per Tree)",
                            title="CIs of Estimated Treatment Means - November Encarsia Count, Study 2")
