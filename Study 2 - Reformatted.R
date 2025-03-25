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

#### (2) Nonparametric Analysis ####

#### July analysis -- mean live scale

# Friedman test
# Friedman test only supports unreplicated complete block designs
# This needs to be averaged across all treatments and will be less precise
scalecount2_avg_july <- average_counts_across_twigs(scalecount2_july)
scalecount2_avg_across_block_trt <- average_counts_across_block_trt(scalecount2_avg_july)

# Note that after grouping data by treatment-block combination, there are only 16 observations!
# Is that going to be a problem...??

# For some reason, this won't work without converting to a matrix
scalecount2_avg_block_trt_matrix <- as.matrix(scalecount2_avg_across_block_trt)

#### Friedman test for July -- mean live scale
# For some reason this only works if you use as.matrix lol
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount2_avg_block_trt_matrix)
friedman

# P-value of 0.5456 from Friedman test above is not significant!
# Therefore we do not need to conduct the Nemenyi (or other) pairwise test!

#### November analysis -- mean live scale

# Friedman test
# This needs to be averaged across all treatments and will be less precise
scalecount2_avg_nov <- average_counts_across_twigs(scalecount2_nov)
scalecount2_avg_across_block_trt <- average_counts_across_block_trt(scalecount2_avg_nov)

# Again, note that the data grouped by treatment-block combination only has 16 observations

# For some reason, this won't work without converting to a matrix
scalecount2_avg_block_trt_matrix <- as.matrix(scalecount2_avg_across_block_trt)

#### Friedman test for November -- mean live scale
friedman <- friedman.test(Meanlivescale ~ Treatment | Block,
                          data = scalecount2_avg_block_trt_matrix)
friedman

# P-value of 0.1562 from Friedman test above is not significant!
# Therefore we do not need to conduct the Nemenyi (or other) pairwise test!

#### (3) Binomial Model ####

#### July analysis -- parasitism





#### (4) Zero-Inflated Models ####

# We can use the "pscl" package for zero-inflated models using the zeroinfl() function!
# However, I'm not exactly sure of the correct way to input the formula for the function...





