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

source("graphing functions.r")
source("Data Processing.r")

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

#### (2) Nonparametric tests? ####

# Testing Friedman test
# Friedman test only supports unreplicated complete block designs
# This needs to be averaged across all treatments and will be less precise
scalecount_avg_july <- average_counts_across_twigs(scalecount_july)
scalecount_avg_across_block_trt <- 
  average_counts_across_block_trt(scalecount_avg_july)

# For some reason, this won't work without converting to a matrix
scalecount_avg_block_trt_matrix <- as.matrix(scalecount_avg_across_block_trt)

# Friedman test for July
# For some reason this only works if you use as.matrix lol
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
              data=scalecount_avg_block_trt_matrix)

# Trying Nemenyi test because it doesn't require restructing the data
# Lower power than Conover
nemenyi <- frdAllPairsNemenyiTest(Meanlivescale ~ Treatment | Block, 
                       data = scalecount_avg_block_trt_matrix)
PMCMRTable(nemenyi)

#### (3) Binomial model test ####

# This doesn't converge if we include Label:Twigab (tree+twig) lol...
# Hopefully we can say for the purposes of parasitism
# it's okay to just treat each twig separately
# Alternatively, we could process this to be per-tree - 1 if either twig has a 1
parasitism_mod <- glm(formula = Prespara ~ Treatment + Block,
                   data = scalecount_july,
                   family = binomial)

emm <- emmeans(parasitism_mod, "Treatment")
pairs(emm)
