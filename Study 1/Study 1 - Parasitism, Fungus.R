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

# Reformatted version of Study 1. 
# Maybe just use for some of the analysis? Live scale analysis?
# Make study 1 a folder with sub-files instead?

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


#### (2) Binomial Model - Parasitism ####

# Client dropped twigs with 0 scale
scalecount_para_july <- scalecount_july[rowSums(scalecount_july[, 6:11] == 0) < 2,]
scalecount_para_july <- 
  scalecount_para_july %>% 
  drop_na(Label) #some samples only have one twig,

# compress to presence per tree
scalecount_para_july_tree <- get_presence_across_twigs(scalecount_para_july)

# This doesn't converge if we include Label:Twigab (tree+twig) with the per-twig observations
# Hopefully we can say for the purposes of parasitism
# it's okay to just treat each twig separately

# Per-tree model
parasitism_mod <- glm(formula = Prespara ~ Treatment + Block,
                      data = scalecount_para_july,
                      family = binomial)

emm <- emmeans(parasitism_mod, "Treatment")
pairs(emm)

#### (3) Binomial Mixed Model - Parasitism ####

# Testing a generalized linear mixed model where we remove the block
# but instead include the tree as a mixed effect
test_glmer <- glmer(Prespara ~ Treatment + (1|Label),
                    data=scalecount_para_july, family = binomial(link = "logit"))
