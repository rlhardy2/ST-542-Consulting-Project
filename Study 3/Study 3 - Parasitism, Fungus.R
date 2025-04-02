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

# Only keep samples where scale (either live or dead) was found
scalecount3_para <- scalecount3 %>% filter(Meanlivescale > 0 | Meandeadscale > 0)

#### (2) Parasitism models ####

para_mod_nov3 <-  glmmTMB(Prespara ~ Treatment + (1 | Block / Label),
                         data=scalecount3_para,
                         family=binomial)
simr_para_mod_nov3 <- simulateResiduals(para_mod_nov3)
plot(simr_para_mod_nov3)

#### (3) Parasitism analysis ####
emm_para_nov3 <- emmeans(para_mod_nov3, "Treatment")
pairs(emm_para_nov3, type="response", adjust="BH")
nov_para_means_comp <- pairs(regrid(emm_para_nov3), adjust="BH")
# Confint for pairwise comparison
confint(nov_para_means_comp)
# Confint for means
confint(emm_para_nov3, type="response")


#### (4) Fungus models ####

fung_mod_nov3 <-  glmmTMB(Presfungus ~ Treatment + (1 | Block / Label),
                          data=scalecount3_para,
                          family=binomial)
simr_fung_mod_nov3 <- simulateResiduals(fung_mod_nov3)
plot(simr_fung_mod_nov3)

#### (5) Fungus analysis ####
emm_fung_nov3 <- emmeans(fung_mod_nov3, "Treatment")
pairs(emm_fung_nov3, type="response", adjust="BH")
nov_fung_means_comp <- pairs(regrid(emm_fung_nov3), adjust="BH")
# Confint for pairwise comparison
confint(nov_fung_means_comp)
# Confint for means
confint_emm_fung_nov3_orig_scale <- confint(emm_fung_nov3, type="response")


#### Fungus Graphs ####
pairs_fung_nov_df3 <- as.data.frame(pairs(regrid(emm_fung_nov3), adjust="BH"))
confint_fung_nov_df3 <- as.data.frame(confint_emm_fung_nov3_orig_scale)
y_positions <- seq(1, 1.9, by=.1)

get_cis_marginal_means_plot(ci_df=confint_fung_nov_df3, pairs_df=pairs_fung_nov_df3, 
                            y_positions=y_positions, trt_labels=trt_labels3,
                            y_str="prob", 
                            y_lab="Probability",
                            title="CIs of Estimated Treatment Means - Entomopathogenic Fungus, Study 3")


