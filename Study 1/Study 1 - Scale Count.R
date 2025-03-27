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

# Treatment survival means
tmeans_july <- get_treatment_survival_means(scalecount_july)
tmeans_nov <- get_treatment_survival_means(scalecount_nov)

#### (2) Friedman + Nemenyi ####

# Friedman test for block design
# Study was technically blocked, but this may not be totally necessary given
# that the blocking was performed randomly?

# Friedman test only supports unreplicated complete block designs
# This needs to be averaged across all treatments and will be less precise
scalecount_tree_mean_july <- average_counts_across_twigs(scalecount_july)
scalecount_tree_mean_across_block_trt <- 
  average_counts_across_block_trt(scalecount_tree_mean_july)

# For some reason, this won't work without converting to a matrix
scalecount_tree_mean_block_trt_matrix <- as.matrix(scalecount_tree_mean_across_block_trt)

#### Friedman test for July
# For some reason this only works if you use as.matrix lol
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
              data = scalecount_tree_mean_block_trt_matrix)
friedman

# Trying Nemenyi test because it doesn't require restructing the data
# Lower power than Conover
nemenyi <- frdAllPairsNemenyiTest(Meanlivescale ~ Treatment | Block, 
                       data = scalecount_tree_mean_block_trt_matrix)
PMCMRTable(nemenyi)

# Extract p-values from the Nemenyi test
p_values_july <- nemenyi$p.value

# Graphing the p-values using a heat map... weird but I found this online lol
heatmap(p_values_july, main = "Nemenyi Test P-Values",
        col = heat.colors(10),
        scale = "none")

# Doesn't work for November 
# Not unreplicated complete block design due to uncollected labels
scalecount_tree_mean_nov <- average_counts_across_twigs(scalecount_nov)


#### (3) Kruskal-Wallis + Dunn ####
# Since the blocks were not truly specified (randomly assigned)
# we may just be able to collapse this to one factor
# if we use the by-tree data set? Ask Dr. Harris...

# Insignificant results for July?
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_tree_mean_july)

# Significant results for November
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_tree_mean_nov)
dunn_nov <- dunnTest(Meanlivescale ~ Treatment, 
                     data = scalecount_tree_mean_nov, method = "bh")

print(dunn_nov, dunn.test.results=TRUE)

#### (4) Poisson/Negative binomial ####

# Mixed Poisson, no zero inflation
pois_nov_mm <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount_nov, ziformula = ~0,
                     family = poisson)
# Doesn't look very good...
simr_pois_nov <- simulateResiduals(pois_nov_mm)

# Mixed NB model, no zero inflation
nb_nov_mm <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount_nov, ziformula = ~0,
                     family = nbinom2)
# Lots of zero dispersion, deviation, etc - 
# probably better to just use zero inflation model
simr_nb_nov <- simulateResiduals(nb_nov_mm)

#### (5) Zero-Inflated & Hurdle Models ####

# Diagnostics easier with glmmTMB than PSCL due to DHARMa compatibility

# Random block, zero inflation, over each twig
# Random block seems to make the most sense in that we don't
# exactly know the effect of each block and it would vary across time, etc
# (soil quality changes?)
# Actually, since the tree is within the block, I guess this theoretically should be nested?
# Tree within block - doesn't matter in this case since each block has a different name
# Will be the same if you do (1 | Block) + (1 | Label) ("crossed factors")
zinb_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                                  data=scalecount_nov, ziformula = ~1,
                                  family = nbinom2)
# looks quite nice
simr_zinb_nov <- simulateResiduals(zinb_nov)


# Mixed NB model, mixed eff for tree with no block, zero inflation
# May be okay because block indicates location too?
# Maybe ask Dr. Harris
zinb_nov_noblock <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1|Label),
                               data=scalecount_nov, ziformula = ~1,
                               family = nbinom2)
# No issues
simr_nov_noblock <- simulateResiduals(zinb_nov_noblock)


# Zero-inflated Poisson
zip_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                   data=scalecount_nov, ziformula = ~1,
                   family = poisson)
# some residual variance issues
simr_zip_nov <- simulateResiduals(zip_nov)

# Hurdle negative binomial
hnbinom_nov <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                        data=scalecount_nov,
                        ziformula = ~1,
                        family=truncated_nbinom1)
# residual variance issues
simr_hnb_nov <- simulateResiduals(hnbinom_nov)

# Hurdle Poisson
hpois_nov <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                        data=scalecount_nov,
                        ziformula = ~1,
                        family=truncated_poisson)
# Homogeneity of variance issues
simr_hpois_nov <- simulateResiduals(hpois_nov)

# Check performance
check_model(zinb_nov)
performance(zinb_nov)
plot(simr_zinb_nov)


#### (6) Estimates, CIs, Treatment Means ####

# Use zero-inflated NB mixed model with nested tree within block for testing

# Estimated marginal means
emm_zinb_nov <- emmeans(zinb_nov, "Treatment")
emm_zinb_nov_orig_scale <- emmeans(zinb_nov, 
                                      "Treatment", type="response")

# effect size - Cohen's d
eff_size(emm_zinb_nov, 
         sigma=sigma(zinb_nov), edf=df.residual(zinb_nov))


# EMMs differ if arrows don't overlap
# https://cran.r-project.org/web/packages/emmeans/vignettes/xplanations.html#arrows
plot(emm_zinb_nov_orig_scale, comparison=TRUE)

# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_zinb_nov_orig_scale, adjust="BH"))

# CI for pairwise comparison on log scale
confint(pairs(emm_zinb_nov, adjust="BH"))

# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_zinb_nov), adjust="BH"))
