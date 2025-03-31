# Loading necessary packages
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
library(biostat3)
library(multcompView)
library(FSA)

source("../graphing functions.r")

#### (1) Data Pre-processing ####

# Reading in the data
scalecount <- read.csv(file = "acety counting v3 (study 3).csv", strip.white=TRUE)

# Adding column names
colnames(scalecount) <- c("Label", "Twigab","EHSLivescale1","EHSDeadscale1","EHSLivescale2","EHSDeadscale2",
                          "EHSLivescale3","EHSDeadscale3","Prespara","CryptoLivescale1","CryptoDeadscale1",
                          "CryptoLivescale2","CryptoDeadscale2","CryptoLivescale3","CryptoDeadscale3",
                          "Presfungus")

# Creating live-scale and dead-scale variables
scalecount$Livescale1 <- scalecount$EHSLivescale1 + scalecount$CryptoLivescale1
scalecount$Livescale2 <- scalecount$EHSLivescale2 + scalecount$CryptoLivescale2
scalecount$Livescale3 <- scalecount$EHSLivescale3 + scalecount$CryptoLivescale3
scalecount$Deadscale1 <- scalecount$EHSDeadscale1 + scalecount$CryptoDeadscale1
scalecount$Deadscale2 <- scalecount$EHSDeadscale2 + scalecount$CryptoDeadscale2
scalecount$Deadscale3 <- scalecount$EHSDeadscale3 + scalecount$CryptoDeadscale2

scalecount<-subset(scalecount, select = -c(EHSLivescale1,EHSLivescale2,EHSLivescale3,EHSDeadscale1,
                                           EHSDeadscale2,EHSDeadscale3,CryptoLivescale1,
                                           CryptoLivescale2,CryptoLivescale3,CryptoDeadscale1,
                                           CryptoDeadscale2,
                                           CryptoDeadscale3))

# Getting 0s and 1s for presence and absence
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="yes", 1)
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="y", 1)
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="Y", 1)
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="", 0)
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="n", 0)
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="N", 0)
scalecount$Prespara <- replace(scalecount$Prespara, scalecount$Prespara=="no", 0)

scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="yes", 1)
scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="y", 1)
scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="Y", 1)
scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="", 0)
scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="n", 0)
scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="N", 0)
scalecount$Presfungus <- replace(scalecount$Presfungus, scalecount$Presfungus=="no", 0)

# Making certain variables numeric
scalecount <- scalecount %>% 
  mutate(across(c(Livescale1, Deadscale1, Livescale2, Deadscale2, Livescale3, Deadscale3), as.numeric))

# New variables for mean of live-scale and dead-scale (averages the three shoots)
scalecount<-scalecount %>% 
  mutate(Meanlivescale = rowMeans(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
scalecount<-scalecount %>% 
  mutate(Meandeadscale = rowMeans(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

# New variables for sum of live-scale and dead-scale (sums the three shoots)
scalecount<-scalecount %>% 
  mutate(Sumlivescale = rowSums(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
scalecount<-scalecount %>% 
  mutate(Sumdeadscale = rowSums(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

# Removing labels that were not collected
scalecount <- scalecount %>% drop_na(Livescale1)

# Making new data set
scalecountall <- scalecount[ -c(5:10) ]

# Extract and create new column with treatment variable, taking on values 1 through 5
scalecountall$Treatment <- substring(scalecountall$Label, first=2, last=2)

# Making the treatment variable a factor variable
scalecountall <- scalecountall %>% mutate(across(c(Treatment), as.factor))

#### End of data cleaning

#### (2) Graphical and Table Summaries ####

# Treatment survival means (mean of the three shoots per sample)
tmeans_mean <- scalecountall %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_mean)

# Treatment survival sums (sum of the live insects on the three shoots per sample)
tmeans_sum <- scalecountall %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(Sumlivescale, na.rm = T),3),
    sd = sd(Sumlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_sum)

#### Treatment means as bar plots (use this one)
# Technically this says EHS, but there are *some* cryptomeria scale? about 5 observations with
labels <- c("Acetamiprid", "Dinotefuron", "Sulfoxaflor", "Flupyradifurone", "No Treatment" )
tmeans_mean$Treatment <- factor(tmeans_mean$Treatment, levels=c("2","1","3","5","4"))
livescale_plot_mean <- ggplot(data=tmeans_mean, 
                         aes(x=Treatment, y=Mean
                         ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x = "Treatment",
       y = "Mean",
       title = "Study 3 - Mean Live EHS & Crypto Scale") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_discrete(label = labels) +
  scale_fill_manual(values = c("#CC6600","red","#9999FF","#0072B2","green4"))  
livescale_plot_mean


tmeans_sum$Treatment <- factor(tmeans_sum$Treatment, levels=c("2","1","3","5","4"))
livescale_plots_sum <- ggplot(data=tmeans_sum, 
                         aes(x=Treatment, y=Mean
                         ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x = "Treatment",
       y = "Mean",
       title = "Study 3 - Mean Sums of Live Scale on 10 needles \n of This Year's Growth") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_discrete(label = labels) +
  scale_fill_manual(values = c("#CC6600","red","#9999FF","#0072B2","green4")) 
livescale_plots_sum

# Get histograms
scalecountall$Treatment <- factor(scalecountall$Treatment,
                                  levels=c("2","1","3","5","4"))
get_hist(data=scalecountall, 
         x_str="Meanlivescale", 
         x_lab="Mean Live Scale",
         title="Study 3 - Mean Numbers of Live Scale Insects",
         labels=labels)

get_hist(data=scalecountall, 
         x_str="Sumlivescale", 
         x_lab="Sum of Live Scale",
         title="Study 3 - Sum of Live Scale Insects",
         labels=labels)

get_hist_all_trt(scalecountall, x_str="Meanlivescale", x_lab="Mean Live Scale",
                 title="Study 3 - Mean Live Scale, All Treatments",
                 labels=labels)

#### (3) Shapiro-Wilk tests (code by Rachel) ####

shapiro.test(scalecountall$Meanlivescale)
shapiro.test(scalecountall$Sumlivescale)

#### (4) Kruskal-Wallis tests (code by Rachel) ####

kruskal.test(Meanlivescale ~ Treatment, data = scalecountall)
kruskal.test(Sumlivescale ~ Treatment, data = scalecountall)

#### (5) Dunn's Test (code by Rachel) ####

dunnTest(Meanlivescale ~ Treatment, data = scalecountall, method = "bh")
dunnTest(Sumlivescale ~ Treatment, data = scalecountall, method = "bh")

#### (6) Parasitism and Fungus Analysis ####

# Making Prespara and Presfungus variables numeric
scalecountall <- scalecountall %>% 
  mutate(across(c(Prespara, Presfungus), as.numeric))

# Now with no zeros: taking out all samples in which no scale was found
scalecount_nozero <- scalecountall

# Note - this is from the original code
# This looks at Prespara, Presfungus, Meanlivescale, Meandeadscale?
# scalecount_nozero <- scalecount_nozero[rowSums(scalecount_nozero[, 3:6] == 0) < 2,]
scalecount_nozero <- 
  scalecount_nozero %>%
  filter(Meanlivescale>0 | Meandeadscale>0)

# Percentage of parasitism presence in a group
tmeans_para <- scalecount_nozero %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Prespara == 1),
    percent_yes100 = round(percent_yes*100),
  )

# Treatment means as bar plots for parasitism no zeros
labels <- c("Acetamiprid", "Dinotefuron", "Sulfoxaflor", "Flupyradifurone", "No Treatment" )
tmeans_para$Treatment <- factor(tmeans_para$Treatment, levels=c("2","1","3","5","4"))
livescale_plots_para <- ggplot(data=tmeans_para,
                              aes(x=Treatment, y=percent_yes100
                              ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +

  labs(x = "Treatment",
       y = "Percent",
       title = "Study 3 - Mean Presence Percentage of Parasitism") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_discrete(label = labels)
livescale_plots_para
ggsave(livescale_plots_para, file="livescale_plots_para.pdf", 
       width = 6, height=4)

# Percentage of fungus in a treatment plot no zeros
scalecount_nozero <- scalecount_nozero %>% 
  mutate(across(c(Presfungus), as.numeric))

# Percentage of fungus presence in a group (Rachel added this code)
tmeans_fungus <- scalecount_nozero %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )

# Treatment means as bar plots for fungus no zeros
# Simplified "entofungus" to just "fungus" for ST542 report
tmeans_fungus$Treatment <- factor(tmeans_fungus$Treatment, levels=c("2","1","3","5","4"))
livescale_plots_fungus <- ggplot(data=tmeans_fungus, 
                               aes(x=Treatment, y=percent_yes100
                               ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  labs(x="Treatment",
       y = "Percent",
       title = "Study 3 - Mean Presence Percentage of Fungus") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_discrete(label = labels)
livescale_plots_fungus
ggsave(livescale_plots_fungus, file="livescale_plots_fungus.pdf", 
       width = 6, height=4)
