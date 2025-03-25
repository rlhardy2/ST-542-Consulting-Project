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

source("graphing functions.r")
source("Data Processing.r")

# Reformatted version of Study 1. 
# Maybe just use for some of the analysis? Live scale analysis?
# Make study 1 a folder with sub-files instead?

#### (1) Data Pre-processing ####
trt_labels <- c("May 2-3", "May 17", "May 28", "No Treatment") # Are these correct?

# Read in file
scalecount <- read.csv(file="EHS count 2024 v7 (study 1).csv", strip.white=TRUE)
colnames(scalecount) <- c("Label", "County","Twigab","Date","Counter","Livescale1","Deadscale1",
                          "Livescale2","Deadscale2","Livescale3","Deadscale3","Prespara",
                          "Presfungus","Presscalenewgr","encarsia")

# Process data set
scalecount <- process_scalecount(scalecount)

# Also change county to factor
scalecount$County <- as.factor(scalecount$County)

# July and November data
scalecount_july <- subset(scalecount, grepl('July', Date))
scalecount_nov <- subset(scalecount, grepl('November', Date))

# Treatment survival means
tmeans_july <- get_treatment_survival_means(scalecount_july)
tmeans_nov <- get_treatment_survival_means(scalecount_nov)

#### (2) Nonparametric Analysis ####

#### July analysis -- mean live scale

# Friedman test
# Friedman test only supports unreplicated complete block designs
# This needs to be averaged across all treatments and will be less precise
scalecount_avg_july <- average_counts_across_twigs(scalecount_july)
scalecount_avg_across_block_trt <- average_counts_across_block_trt(scalecount_avg_july)

# For some reason, this won't work without converting to a matrix
scalecount_avg_block_trt_matrix <- as.matrix(scalecount_avg_across_block_trt)

#### Friedman test for July -- mean live scale
# For some reason this only works if you use as.matrix lol
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
              data = scalecount_avg_block_trt_matrix)
friedman

# Trying Nemenyi test because it doesn't require restructing the data
# Lower power than Conover
nemenyi <- frdAllPairsNemenyiTest(Meanlivescale ~ Treatment | Block, 
                       data = scalecount_avg_block_trt_matrix)
PMCMRTable(nemenyi)

# Extract p-values from the Nemenyi test
p_values_july <- nemenyi$p.value

# Graphing the p-values using a heat map... weird but I found this online lol
heatmap(p_values_matrix, main = "Nemenyi Test P-Values",
        col = heat.colors(10),
        scale = "none")

#### November analysis -- mean live scale (doesn't work)

# We do not have an unreplicated complete block design due to missing data!!
# Therefore we cannot use the Friedman test!!

# Again, but for November (doesn't work!)
scalecount_avg_nov <- average_counts_across_twigs(scalecount_nov)
scalecount_avg_across_block_trt <- 
  average_counts_across_block_trt(scalecount_avg_nov)

# For some reason, this won't work without converting to a matrix
scalecount_avg_block_trt_matrix <- as.matrix(scalecount_avg_across_block_trt)

#### Friedman test for November -- doesn't work!
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount_avg_block_trt_matrix)

#### (3) Binomial Model ####

#### July analysis -- parasitism

# Client dropped twigs with 0 scale
scalecount_para_july <- scalecount_july[rowSums(scalecount_july[, 6:11] == 0) < 2,]
scalecount_para_july <- 
  scalecount_para_july %>% 
  drop_na(Label) #some samples only have one twig,

# Compress to presence per tree
scalecount_para_july_tree <- get_presence_across_twigs(scalecount_para_july)

# This doesn't converge if we include Label:Twigab (tree+twig) with the per-twig observations
# Hopefully we can say for the purposes of parasitism
# it's okay to just treat each twig separately

# Per tree model
parasitism_mod_july <- glm(formula = Prespara ~ Treatment + Block,
                           data = scalecount_para_july_tree,
                           family = binomial)

emm <- emmeans(parasitism_mod_july, "Treatment")
pairs(emm)

#### November analysis -- parasitism






#### (4) Zero-Inflated Models ####

# We can use the "pscl" package for zero-inflated models using the zeroinfl() function!
# However, I'm not exactly sure of the correct way to input the formula for the function...

