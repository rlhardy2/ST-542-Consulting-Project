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

##### July #####
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

##### November #####
# Doesn't work for November 
# Not unreplicated complete block design due to uncollected labels
scalecount_tree_mean_nov <- average_counts_across_twigs(scalecount_nov)


#### (3) Kruskal-Wallis + Dunn ####
# Since the blocks were not truly specified (randomly assigned)
# we may just be able to collapse this to one factor if 
# we use the by-tree data set?
# However, blocks seem to have quite an effect in Study 1 
# (may be bc of diff counties)

# Insignificant results for July?
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_tree_mean_july)

# Significant results for November
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_tree_mean_nov)
dunn_nov <- dunnTest(Meanlivescale ~ Treatment, 
                     data = scalecount_tree_mean_nov, method = "bh")

print(dunn_nov, dunn.test.results=TRUE)

#### (4) Poisson & Negative Binomial Models ####

# Just to test fit and necessity of zero inflation

##### July #####
# Mixed Poisson, no zero inflation
pois_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                       data=scalecount_july, ziformula = ~0,
                       family = poisson)
# Doesn't look very good...
simr_pois_july <- simulateResiduals(pois_july)

# Mixed NB model, no zero inflation
# Negative Binomial 2 (typical)
nb2_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount_july, ziformula = ~0,
                     family = nbinom2)
# Has some zero inflation issues
# Also deviation issues
simr_nb2_july <- simulateResiduals(nb2_july)
testCategorical(simr_nb2_july, scalecount_july$Treatment)

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount_july, ziformula = ~0,
                   family = nbinom1)
# Doesn't strictly have zero inflation
# Actually looks quite good?
# Should probably have similar assumptions as November though??
simr_nb1_july <- simulateResiduals(nb1_july)
# Homogeneity of variance looks okay...
testCategorical(simr_nb1_july, scalecount_july$Treatment)

##### November #####
# Mixed Poisson, no zero inflation
pois_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount_nov, ziformula = ~0,
                     family = poisson)
# Doesn't look very good...
simr_pois_nov <- simulateResiduals(pois_nov)
check_overdispersion(pois_nov)

# Mixed NB model, no zero inflation
# Using type 2 nbinom as is typical
nb2_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                     data=scalecount_nov, ziformula = ~0,
                     family = nbinom2)
# Doesn't technically fail zero inflation test
# But lots of issues
simr_nb2_nov <- simulateResiduals(nb2_nov)
testCategorical(simr_nb2_nov, scalecount_nov$Treatment)

# Mixed NB1 model, no zero inflation
nb1_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                  data=scalecount_nov, ziformula = ~0,
                  family = nbinom1)
# Doesn't technically fail zero inflation test
# Issue with homogeneity of variance
simr_nb1_nov <- simulateResiduals(nb1_nov)
plot(simr_nb1_nov)
testCategorical(simr_nb1_nov, scalecount_nov$Treatment)

#### (5) Zero-Inflated & Hurdle Models ####

##### July #####
zinb2_july <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                    data=scalecount_july, ziformula = ~1,
                    family = nbinom2)
# looks quite nice
simr_zinb2_july <- simulateResiduals(zinb2_july)
plot(simr_zinb2_july)

##### November #####
# Diagnostics easier with glmmTMB than PSCL due to DHARMa compatibility

# Random block, zero inflation, over each twig
# Random block seems to make the most sense in that we don't
# exactly know the effect of each block and it would vary across time, etc
# (soil quality changes?)
# Actually, since the tree is within the block, I guess this theoretically should be nested?
# Tree within block - doesn't matter in this case since each block has a different name
# Will be the same if you do (1 | Block) + (1 | Label) ("crossed factors")
# Using negative binomial 2 - quadratic overdispersion
zinb2_nov <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                                  data=scalecount_nov, ziformula = ~1,
                                  family = nbinom2)
# looks quite nice
simr_zinb2_nov <- simulateResiduals(zinb2_nov)

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
                        family=truncated_nbinom2)
# residual variance issues
simr_hnb2_nov <- simulateResiduals(hnbinom_nov)

# Hurdle Poisson
hpois_nov <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                        data=scalecount_nov,
                        ziformula = ~1,
                        family=truncated_poisson)
# Homogeneity of variance issues
simr_hpois_nov <- simulateResiduals(hpois_nov)

# Check performance
check_model(zinb2_nov)
performance(zinb2_nov)
plot(simr_zinb2_nov)

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
# https://cran.r-project.org/web/packages/emmeans/vignettes/xplanations.html#arrows
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
emm_nb1_nov <- emmeans(nb1_nov, "Treatment")
emm_nb1_nov_orig_scale <- emmeans(nb1_nov, 
                                      "Treatment", type="response")
confint(emm_nb1_nov_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_nov, 
         sigma=sigma(nb1_nov), edf=df.residual(nb1_nov))


###### Pairwise Comparisons ######
# EMMs differ if arrows don't overlap
# https://cran.r-project.org/web/packages/emmeans/vignettes/xplanations.html#arrows
plot(emm_nb1_nov_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_nov_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_nov, adjust="BH"))
# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_nb1_nov), adjust="BH"))

# Convert to data frame for graphing
regrid_pairs <- as.data.frame(pairs(regrid(emm_nb1_nov), adjust="BH"))

#### (7) Graphs ####

# Graph of sum of live scale per twig
ggboxplot(scalecount_nov, x="Treatment", y="Sumlivescale_from_mean", 
          ylab="Sum of Live Scale", title="Study 1 - Sum of Live Scale Per Twig", ggtheme=theme_gray(base_size=12))
