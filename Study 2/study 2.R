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

source("graphing functions.r")

#### (1) Data Pre-processing ####

# Reading in the data
scalecount1 <- read.csv(file = "distance numbers v7 (study 2).csv", strip.white=TRUE)

# Adding column names
colnames(scalecount1) <- c("Label","Type", "Twigab","Date","Livescale1","Deadscale1","Livescale2",
                           "Deadscale2","Livescale3","Deadscale3","Prespara",
                          "Presscalenewgr","Presfungus","encarsia","notes")

# Sub-setting the data
scalecount1 <- subset(scalecount1, select = c("Label","Type", "Twigab","Date","Livescale1",
                                              "Deadscale1","Livescale2","Deadscale2",
                                              "Livescale3","Deadscale3","Prespara","Presfungus",
                                              "encarsia"))

# Extract and create new column with treatment
scalecount1$Treatment <- substring(scalecount1$Label, first=4, last=4)

scalecount1 <- scalecount1 %>% drop_na(Livescale1)

# Getting 0s and 1s for presence and absence
scalecount1$Prespara <- replace(scalecount1$Prespara, scalecount1$Prespara=="yes", 1)
scalecount1$Prespara <- replace(scalecount1$Prespara, scalecount1$Prespara=="no", 0)
scalecount1$Presfungus <- replace(scalecount1$Presfungus, scalecount1$Presfungus=="yes", 1)
scalecount1$Presfungus <- replace(scalecount1$Presfungus, scalecount1$Presfungus=="no", 0)

# Making the chosen variables numeric
scalecount1 <- scalecount1 %>% 
  mutate(across(c(Livescale1, Deadscale1, Livescale2, Deadscale2, Livescale3, Deadscale3), as.numeric))
str(scalecount1)

# Making a new table that adds the EHS and cryptomeria
scalecountboth <- scalecount1 %>% 
  group_by(Label,Twigab, Date, Treatment) %>% 
  dplyr::summarize(across(where(is.numeric), sum))
scalecountboth

# Getting new column 
# Get another data frame to send to no NA- this is the column we will use
# No NA means no missing values -- this is a data set with no missing values
scalecountboth_NoNA <- scalecountboth
scalecountboth$Treatment <- as.factor(scalecountboth$Treatment)
scalecountboth$Twigab <- as.factor(scalecountboth$Twigab)

# Must use across NOT SELECT - this is for getting mean scale

# Mean live-scale
scalecountboth_NoNA <- scalecountboth_NoNA %>% 
  mutate(Meanlivescale = rowMeans(across(c(Livescale1, Livescale2, Livescale3)), na.rm = TRUE))
scalecountboth_NoNA

# Mean dead-scale
scalecountboth_NoNA <- scalecountboth_NoNA %>% 
  mutate(Meandeadscale = rowMeans(across(c(Deadscale1, Deadscale2, Deadscale3)), na.rm = TRUE))
scalecountboth_NoNA

# New variables for sum of live-scale and dead-scale (sums the three shoots)
scalecountboth_NoNA <-scalecountboth_NoNA %>% 
  mutate(Sumlivescale = rowSums(across(c(Livescale1, Livescale2, Livescale3)), na.rm = TRUE))
scalecountboth_NoNA <-scalecountboth_NoNA  %>% 
  mutate(Sumdeadscale = rowSums(across(c(Deadscale1, Deadscale2, Deadscale3)), na.rm = TRUE))

# Making Treatment a factor variable
scalecountboth_NoNA$Treatment <- as.factor(scalecountboth_NoNA$Treatment)

#### End of data cleaning

#### Below starts analysis by month (July and November)

# Data sets for July and November
scalecount_July <- subset(scalecountboth_NoNA, grepl('July', Date))
scalecount_Nov <- subset(scalecountboth_NoNA, grepl('Nov', Date))

# Making Treatment a factor variable
scalecount_July$Treatment <- as.factor(scalecount_July$Treatment)
scalecount_Nov$Treatment <- as.factor(scalecount_Nov$Treatment)

#### (2) Shapiro-Wilk tests (code by Rachel) ####

## For July

shapiro.test(scalecount_July$Meanlivescale)
shapiro.test(scalecount_July$Sumlivescale)
shapiro.test(scalecount_July$encarsia)

## For November

shapiro.test(scalecount_Nov$Meanlivescale)
shapiro.test(scalecount_Nov$Sumlivescale)
shapiro.test(scalecount_Nov$encarsia)

#### (3) Kruskal-Wallis tests (code by Rachel) ####

## For whole data set (not sure if this is necessary)

kruskal.test(Meanlivescale ~ Treatment, data = scalecountboth_NoNA)
kruskal.test(Sumlivescale ~ Treatment, data = scalecountboth_NoNA)
kruskal.test(encarsia ~ Treatment, data = scalecountboth_NoNA)

## For July

kruskal.test(Meanlivescale ~ Treatment, data = scalecount_July)
kruskal.test(Sumlivescale ~ Treatment, data = scalecount_July)
kruskal.test(encarsia ~ Treatment, data = scalecount_July)

## For November

kruskal.test(Meanlivescale ~ Treatment, data = scalecount_Nov)
kruskal.test(Sumlivescale ~ Treatment, data = scalecount_Nov)
kruskal.test(encarsia ~ Treatment, data = scalecount_Nov)

#### (4) Dunn's Test (code by Rachel) ####

## For whole data set (not sure if this is necessary)

dunnTest(Meanlivescale ~ Treatment, data = scalecountboth_NoNA,
         method = "bh")
dunnTest(Sumlivescale ~ Treatment, data = scalecountboth_NoNA,
         method = "bh")
dunnTest(encarsia ~ Treatment, data = scalecountboth_NoNA,
         method = "bh") # Warning about missing data

## For July

dunnTest(Meanlivescale ~ Treatment, data = scalecount_July, method = "bh")
dunnTest(Sumlivescale ~ Treatment, data = scalecount_July, method = "bh")
dunnTest(encarsia ~ Treatment, data = scalecount_July, method = "bh") #Warning about missing data

## For November

dunnTest(Meanlivescale ~ Treatment, data = scalecount_Nov, method = "bh")
dunnTest(Sumlivescale ~ Treatment, data = scalecount_Nov, method = "bh")
dunnTest(encarsia ~ Treatment, data = scalecount_Nov, method = "bh") #Warning about missing data

#### (5) Table Summaries ####

# Treatment means for November
tmeans_nov <- scalecount_Nov %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(Meanlivescale, na.rm = T), 3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
tmeans_nov$Collection<-"Nov"
print(tmeans_nov)

tmeans_july <- scalecount_July %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(Meanlivescale, na.rm = T), 3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
tmeans_july$Collection<-"July"

#### (6) Original Plots ####

# Original bar plots from Dr. Bookwalter
# Treatment means as bar plots for November
labels <- c("Pyri/May", "Ace/May", "Pyri/May + Ace/June", "No Treatment")
Meanlivescale_plot_Nov <- ggplot(data = tmeans_nov, 
                                 aes(x = Treatment, y = Mean), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x = "Treatment",
       y = "Mean",
       title = "Mean Numbers of Live EHS & Crypto Scale \n on 10 needles") +
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
  scale_fill_brewer(palette = "Dark2")+ scale_x_discrete(label = labels)
Meanlivescale_plot_Nov
ggsave(Meanlivescale_plot_Nov, file="Meanlivescale_plot_Nov.pdf", 
       width = 6, height=4)

tmeans_table <- dplyr::bind_rows(tmeans_july, tmeans_nov)
plot_means_by_collection(data=tmeans_table, 
                         title="Study 2 - Mean Live EHS & Crypto Scale",
                         x_str="Treatment", y_str="Mean", labels=labels)

#### (7) Graphical Summaries ####

# Histograms of means
get_hist_livescale(scalecount_July, "July", 2, labels)
get_hist_livescale(scalecount_Nov, "Nov", 2, labels)

get_hist_all_trt(scalecount_July, x_str="Meanlivescale", x_lab="Mean Live Scale",
                 title="Study 2 - Mean Live Scale, All Treatments, July",
                 labels=labels)

get_hist_all_trt(scalecount_Nov, x_str="Meanlivescale", x_lab="Mean Live Scale",
                 title="Study 2 - Mean Live Scale, All Treatments, Nov",
                 labels=labels)

# Histograms of sums
get_hist(data=scalecount_July, 
         x_str="Sumlivescale", 
         x_lab="Sum",
         title="Study 2 - Sum of Live EHS & Crypto Scale, July",
         labels=labels)
get_hist(data=scalecount_Nov, 
         x_str="Sumlivescale", 
         x_lab="Sum",
         title="Study 2 - Sum of Live EHS & Crypto Scale, Nov",
         labels=labels)


#### (8) Encarsia Analysis ####

# Note - this should already be 1 per tree
scalecountencar_nov <- scalecount_Nov
scalecountencar_july <- scalecount_July

# Remove twigs with no scale
# Some samples only have one twig, dropping those
scalecountencar_nov <- scalecountencar_nov %>% drop_na(encarsia)
scalecountencar_july <- scalecountencar_july %>% drop_na(encarsia)

# Kruskal-Wallis test for encarsia
kruskal.test(encarsia ~ Treatment, data = scalecountencar_nov)
# Kruskal-Wallis chi-squared = 2.0473, df = 3, p-value = 0.5626
# No difference among treatments

# Treatment means for encarsia - Nov, July
tmeans_encar_nov <- scalecountencar_nov %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(encarsia, na.rm = T),3),
    sd = sd(encarsia),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
tmeans_encar_nov$Collection<-"Nov"
print(tmeans_encar_nov)

tmeans_encar_july <- scalecountencar_july %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean = round(mean(encarsia, na.rm = T),3),
    sd = sd(encarsia),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
tmeans_encar_july$Collection<-"July"

# Put both data sets in one table for graphing
encar_table_trt <- dplyr::bind_rows(tmeans_encar_july, 
                                    tmeans_encar_nov)

# Treatment means for encarsia as bar plots
labels <- c("Pyri/May", "Ace/May", "Pyri/May + Ace/June", "No Treatment")
encarsia_plot<- ggplot(data = tmeans_encar_nov, 
                       aes(x = Treatment, y = Mean), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x = "Treatment",
       y = "Mean",
       title = "Mean numbers of Encarsia") +
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
  )+
  scale_fill_brewer(palette = "Dark2")+ scale_x_discrete(label = labels)
encarsia_plot
ggsave(encarsia_plot, file="encarsia_plot.pdf", 
       width = 6, height=4)

# Plot of both July and November means
tmeans_encar_plot_trt <- 
  plot_means_by_collection(data=encar_table_trt, 
                                title="Study 2 - Mean Encarsia Count",
                                x_str="Treatment",
                                y_str="Mean",
                                labels = labels)
tmeans_encar_plot_trt

# Histograms
get_hist_encarsia(data=scalecountencar_july, collection_date="July", study=2, labels=labels)
get_hist_encarsia(data=scalecountencar_nov, collection_date="Nov", study=2, labels=labels)

get_hist_all_trt(scalecountencar_july, x_str="encarsia", 
                 x_lab="Count",
                 title="Study 2 - Encarsia Count, All Treatments, July",
                 labels=labels)

get_hist_all_trt(scalecountencar_nov, x_str="Meanlivescale", 
                 x_lab="Count",
                 title="Study 2 - Encarsia Count, All Treatments, Nov",
                 labels=labels)

#### (9) Parasitism Analysis ####

#### Fungus (November) (This is actually parasitism!!!)

# Note - this won't work. I realized the count with Prespara isn't quite correct.
# Commented this out so it doesn't affect running the file...
# Must use data set without summed scale counts here instead!

# scalecountfung_july <- scalecount_July
# scalecountfung_nov <- scalecount_Nov
# 
# # Remove twigs with no scale
# # Some samples only have one twig, dropping those
# scalecountfung_july <- scalecountfung_july[rowSums(scalecountfung_july[, 6:11] == 0) < 2,]
# scalecountfung_july <- scalecountfung_july %>% drop_na(Label)
# scalecountfung_nov <- scalecountfung_nov[rowSums(scalecountfung_nov[, 6:11] == 0) < 2,]
# scalecountfung_nov <- scalecountfung_nov %>% drop_na(Label)
# 
# scalecountfung_july$Prespara <- as.numeric(scalecountfung_july$Prespara)
# scalecountfung_nov$Prespara <- as.numeric(scalecountfung_nov$Prespara)
# 
# # Percentage of presence in a group
# # Fixed code after adding Prespara to scalecount processing
# # Wait this is actually parasitism, not fungus?
# tmeans_fung_nov <- scalecountfung_nov %>% 
#   group_by(Treatment) %>% 
#   summarise(
#     percent_yes = mean(Prespara == 1),
#     percent_yes100 = round(percent_yes*100),
#   )
# tmeans_fung_nov$Collection <- "Nov"
# 
# tmeans_fung_july <- scalecountfung_july %>% 
#   group_by(Treatment) %>% 
#   summarise(
#     percent_yes = mean(Prespara == 1),
#     percent_yes100 = round(percent_yes*100),
#   )
# tmeans_fung_july$Collection <- "July"
# 
# fungus_table_trt <- dplyr::bind_rows(tmeans_fung_july, 
#                                      tmeans_fung_nov)
# 
# plot_means_by_collection(data=fungus_table_trt, 
#                          title="Study 2 - Mean Presence Percentage of Parasitism",
#                          x_str="Treatment",
#                          y_str="percent_yes100",
#                          y_lab="Percent",
#                          labels=labels)
