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

#### Data Pre-processing ####
trt_labels3 <- c("Acetamiprid", "Dinotefuron", "Sulfoxaflor", "Flupyradifurone", "No Treatment" )
scalecount3 <- read.csv(file = "acety counting v3 (study 3).csv", strip.white=TRUE)

colnames(scalecount3) <- c("Label", "Twigab","EHSLivescale1","EHSDeadscale1","EHSLivescale2","EHSDeadscale2",
                          "EHSLivescale3","EHSDeadscale3","Prespara","CryptoLivescale1","CryptoDeadscale1",
                          "CryptoLivescale2","CryptoDeadscale2","CryptoLivescale3","CryptoDeadscale3",
                          "Presfungus")

##### Special - study 3 only #####
# Create Livescale{x}, Deadscale{x} variables
scalecount3 <- create_live_deadscale_vars(scalecount3)
## End study 3 only processing

# Remove unused variables (crypto, ehs scale)
scalecount3 <-subset(scalecount3, select = -c(EHSLivescale1,EHSLivescale2,EHSLivescale3,EHSDeadscale1,
                                           EHSDeadscale2,EHSDeadscale3,CryptoLivescale1,
                                           CryptoLivescale2,CryptoLivescale3,CryptoDeadscale1,
                                           CryptoDeadscale2,
                                           CryptoDeadscale3))

scalecount3 <- process_scalecount(scalecount3)


#### Poisson/Negative binomial ####
# Mixed Poisson, no zero inflation
pois_nov3 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                    data=scalecount3, ziformula = ~0,
                    family = poisson)
# Doesn't look very good...variance issues
simr_pois_nov3 <- simulateResiduals(pois_nov3)
check_overdispersion(pois_nov3)

# Mixed NB model, no zero inflation
# Using type 2 nbinom as is typical
nb2_nov3 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount3, ziformula = ~0,
                   family = nbinom2)
# Variance looks a little weird, but actually pretty decent
# Underfitting zeroes a little
simr_nb2_nov3 <- simulateResiduals(nb2_nov3)
testCategorical(simr_nb2_nov3, scalecount3$Treatment)

# Mixed NB1 model, no zero inflation
nb1_nov3 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1| Block / Label),
                   data=scalecount3, ziformula = ~0,
                   family = nbinom1)

# Variance looks a little weird still, but not bad? (very narrow for Trt 4)
simr_nb1_nov3 <- simulateResiduals(nb1_nov3)
plot(simr_nb1_nov3)
# Zeroes look ok
testCategorical(simr_nb1_nov3, scalecount3$Treatment)


#### Zero-Inflated & Hurdle Models ####

# Zero inflation with nb2
zinb2_nov3 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount3, ziformula = ~1,
                      family = nbinom2)
simr_zinb2_nov3 <- simulateResiduals(zinb2_nov3)
plot(simr_zinb2_nov3)

# Zero inflation with nb 1
zinb1_nov3 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount3, ziformula = ~1,
                      family = nbinom1)
simr_zinb1_nov3 <- simulateResiduals(zinb1_nov3)
plot(simr_zinb1_nov3)

# Zero inflation with poisson
zip_nov3 <- glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                      data=scalecount3, ziformula = ~1,
                      family = poisson)
simr_zip_nov3 <- simulateResiduals(zip_nov3)
# Homogeneity of variance issues
plot(simr_zip_nov3)

# Hurdle negative binomial
hnbinom2_nov3 <-  glmmTMB(Sumlivescale_from_mean ~ Treatment + (1 | Block / Label),
                        data=scalecount3,
                        ziformula = ~1,
                        family=truncated_nbinom2)
# residual variance issues, don't use
simr_hnb2_nov3 <- simulateResiduals(hnbinom2_nov3)
plot(simr_hnb2_nov3)

# Might prefer the standard models for study 3? IDK

#### Analysis ####

# Maybe use nb1...?

# Estimated marginal means
emm_nb1_nov3 <- emmeans(nb1_nov3, "Treatment")
emm_nb1_nov3_orig_scale <- emmeans(nb1_nov3, 
                                     "Treatment", type="response")
confint(emm_nb1_nov3_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_nov3, 
         sigma=sigma(nb1_nov3), edf=df.residual(nb1_nov3))


###### Pairwise Comparisons ######
# EMMs differ if arrows don't overlap
# https://cran.r-project.org/web/packages/emmeans/vignettes/xplanations.html#arrows
plot(emm_nb1_nov3_orig_scale, comparison=TRUE)
# Pairwise comparisons for ratios 
# (happens if you take the pairs from emm on the original scale
# due to needing to perform tests on log)
confint(pairs(emm_nb1_nov3_orig_scale, adjust="BH"))
# CI for pairwise comparison on log scale
confint(pairs(emm_nb1_nov3, adjust="BH"))
# CI for pairwise comparison on original scale
confint(pairs(regrid(emm_nb1_nov3), adjust="BH"))

pairs(regrid(emm_nb1_nov3), adjust="BH")
