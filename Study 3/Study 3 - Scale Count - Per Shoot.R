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
library(ggpubr)

source("../graphing functions.r")
source("../Data Processing.r")

#### (1) Data Pre-processing ####

trt_labels3 <- c("Acetamiprid", "Dinotefuran", "Sulfoxaflor", "Flupyradifurone", "No Treatment")

scalecount3 <- read.csv(file = "acety counting v3 (study 3).csv", strip.white=TRUE)

colnames(scalecount3) <- c("Label", "Twigab","EHSLivescale1","EHSDeadscale1",
                           "EHSLivescale2","EHSDeadscale2",
                           "EHSLivescale3","EHSDeadscale3",
                           "Prespara","CryptoLivescale1","CryptoDeadscale1",
                           "CryptoLivescale2","CryptoDeadscale2",
                           "CryptoLivescale3","CryptoDeadscale3",
                           "Presfungus")

##### Special - Study 3 only #####

# Create Livescale{x}, Deadscale{x} variables
scalecount3 <- create_live_deadscale_vars(scalecount3)
## End study 3 only processing

# Remove unused variables (crypto, ehs scale)
scalecount3 <-subset(scalecount3, select = -c(EHSLivescale1,EHSLivescale2,
                                              EHSLivescale3,EHSDeadscale1,
                                              EHSDeadscale2,EHSDeadscale3,CryptoLivescale1,
                                              CryptoLivescale2,CryptoLivescale3,CryptoDeadscale1,
                                              CryptoDeadscale2,
                                              CryptoDeadscale3))

scalecount3 <- process_scalecount(scalecount3)

##### New Processing for Per-Shoot Modeling ######

# Get long scalecount
scalecount3_long <- get_per_shoot_scalecount(scalecount3)

# Get only live scalecount from long, without NAs
scalecount3_long_live <- get_per_shoot_scalecount_live(scalecount3_long)


#### Model ####

# Mixed Poisson, no zero inflation
pois_nov3 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                     data=scalecount3_long_live, ziformula = ~0,
                     family = poisson)
simr_pois_nov3 <- simulateResiduals(pois_nov3)
plot(simr_pois_nov3, title="Poisson")

# Mixed NB model, no zero inflation
# Negative Binomial 3 (typical) - quadratic overdispersion
nb2_nov3 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount3_long_live, ziformula = ~0,
                    family = nbinom2)
simr_nb2_nov3 <- simulateResiduals(nb2_nov3)
plot(simr_nb2_nov3, title="Negative Binomial 2")

# Mixed NB model, no zero inflation
# Negative Binomial 1 (linear dispersion)
nb1_nov3 <- glmmTMB(scalecount_shoot ~ Treatment + (1| Block / Label / Twigab),
                    data=scalecount3_long_live, ziformula = ~0,
                    family = nbinom1)
simr_nb1_nov3 <- simulateResiduals(nb1_nov3)
plot(simr_nb1_nov3, title="Negative Binomial 1")

#### (7) Analysis ####


##### Means #####
# Using NB1
# Estimated marginal means
emm_nb1_nov3 <- emmeans(nb1_nov3, "Treatment")
emm_nb1_nov3_orig_scale <- emmeans(nb1_nov3, 
                                     "Treatment", type="response")
confint(emm_nb1_nov3_orig_scale)
# effect size - Cohen's d
eff_size(emm_nb1_nov3, 
         sigma=sigma(nb1_nov3), edf=df.residual(nb1_nov3))


##### Pairwise Comparisons #####
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

##### Graphs #####
pairs_nov_df3 <- as.data.frame(pairs(regrid(emm_nb1_nov3), adjust="BH"))
confint_nov_df3 <- as.data.frame(confint(emm_nb1_nov3_orig_scale))
y_positions <- seq(8, 18.8, by=1.2)

get_cis_marginal_means_plot(ci_df=confint_nov_df3, pairs_df=pairs_nov_df3, 
                            y_positions=y_positions, trt_labels=trt_labels3,
                            y_str="response", 
                            y_lab="Scale Count (Per Shoot)",
                            title="CIs of Estimated Treatment Means - Live Scale Count, Study 3")
