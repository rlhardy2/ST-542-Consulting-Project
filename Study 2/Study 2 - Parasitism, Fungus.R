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

source("../graphing functions.r")
source("../Data Processing.r")

# Reformatted version of Study 2.

#### (1) Data Pre-processing ####

# Reading in the data
# Not sure why but there is an extra column of missing values at the end of the data set
# This is removed below along with the 'notes' column
scalecount2 <- read.csv(file = "distance numbers v7 (study 2).csv", strip.white=TRUE)
colnames(scalecount2) <- c("Label","Type", "Twigab","Date","Livescale1","Deadscale1","Livescale2",
                           "Deadscale2","Livescale3","Deadscale3","Prespara",
                           "Presscalenewgr","Presfungus","encarsia","notes","extra")

# Remove notes column and random extra column
scalecount2 <- subset(scalecount2, select = -c(notes, extra))

# Process data set
scalecount2 <- process_scalecount(scalecount2)

# July and November data
scalecount2_july <- subset(scalecount2, grepl('July', Date))
scalecount2_nov <- subset(scalecount2, grepl('November', Date))

# Treatment survival means
tmeans2_july <- get_treatment_survival_means(scalecount2_july)
tmeans2_nov <- get_treatment_survival_means(scalecount2_nov)

#### (2) Parasitism Analysis ???? ####

##### July #####

# Client dropped twigs with 0 scale
scalecount2_para_july <- scalecount2_july[rowSums(scalecount_july[, 5:10] == 0) < 2,]
scalecount2_para_july <- 
  scalecount2_para_july %>% 
  drop_na(Label) # Some samples only have one twig

# Compress to presence per tree
scalecount2_para_july_tree <- get_presence_across_twigs(scalecount2_para_july)

# Per tree model
parasitism_mod2_july <- glm(formula = Prespara ~ Treatment + Block,
                            data = scalecount2_para_july_tree,
                            family = binomial)

emm <- emmeans(parasitism_mod2_july, "Treatment")
pairs(emm)

##### November #####

# Umm, I'm not sure if we should be dropping twigs with 0 scale?
# Wouldn't 0 scale mean the pesticide worked (95% of the time like the client estimated)?
# We have A LOT of observations with 0 scale here, most likely because the pesticide worked
# If we remove all of these observations, we are left with almost nothing to use in the
# binomial model below...
# That's why the code below only shows one pairwise comparison...

# Client dropped twigs with 0 scale
scalecount2_para_nov <- scalecount2_nov[rowSums(scalecount2_nov[, 5:10] == 0) < 2,]
scalecount2_para_nov <- 
  scalecount2_para_nov %>% 
  drop_na(Label) # Some samples only have one twig

# Compress to presence per tree
scalecount2_para_nov_tree <- get_presence_across_twigs(scalecount2_para_nov)

# Per tree model
parasitism_mod2_nov <- glm(formula = Prespara ~ Treatment + Block,
                           data = scalecount2_para_nov_tree,
                           family = binomial)

emm <- emmeans(parasitism_mod2_nov, "Treatment")
pairs(emm)

#### (3) Fungus Analysis ???? ####

##### July #####




##### November #####





