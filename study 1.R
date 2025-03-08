# Load the package, read data, and manipulate
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

source("graphing functions.r")


### PREPROCESSING ###

# read in file
scalecount <- read.csv(file="EHS count 2024 v7 (study 1).csv", strip.white=TRUE)
colnames(scalecount) <- c("Label", "County","Twigab","Date","Counter","Livescale1","Deadscale1",
                          "Livescale2","Deadscale2","Livescale3","Deadscale3","Prespara","Presfungus",
                          "Presscalenewgr","encarsia")

#delete columns with no info
scalecount <- scalecount[ -c(16:19) ]

#extract and create new column with treatment
scalecount$Treatment<-
  substring(scalecount$Label, first=4, last=4)

#getting zero and 1s as presence and absence
# parasitism
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="yes", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="y", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="Y", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="n", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="N", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="no", 0)

# scale new growth?
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="yes", 1)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="y", 1)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="Y", 1)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="n", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="N", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="no", 0)

# presence of fungus
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="yes", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="y", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="Y", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="n", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="N", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="no", 0)

# change Treatment and County to factors
scalecount <- scalecount %>% 
  mutate(across(c(Treatment, County), as.factor))

# change counts to numeric
scalecount <- scalecount %>% 
  mutate(across(c(Livescale1, Deadscale1, Livescale2, Deadscale2, Livescale3, Deadscale3), as.numeric))

# Add mean live scale and Mean dead scale columns
scalecount<-scalecount %>% 
  mutate(Meanlivescale = rowMeans(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
scalecount<-scalecount %>% 
  mutate(Meandeadscale = rowMeans(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

# New variables for sum of live-scale and dead-scale (sums the three shoots)
scalecount<-scalecount %>% 
  mutate(Sumlivescale = rowSums(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
scalecount<-scalecount %>% 
  mutate(Sumdeadscale = rowSums(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

#removing labels that were not collected
scalecount<-scalecount %>% drop_na(Livescale1)

#scalecount has both July and November counts in it

# getting first count- July
# "scalecount_july" henceforth refers to July collection
scalecount_july <- subset(scalecount, grepl('July',Date))

# Treatment survival means
tmeans_all <- scalecount_july %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_all)

# Treatment labels
trt_labels <- c("May 2-3", "May 17", "May 28", "No Treatment")

get_hist_livescale(scalecount_july, collection_date="July", study=1, 
                   labels=trt_labels)
get_hist_livescale(scalecount_nov, collection_date="Nov", study=1, 
                   labels=trt_labels)

get_hist(data=scalecount_july, 
         x_str="Sumlivescale", 
         x_lab="Sum",
         title="Study 1 - Sum of Live EHS, July",
         labels=trt_labels)

get_hist(data=scalecount_nov, 
         x_str="Sumlivescale", 
         x_lab="Sum",
         title="Study 1 - Sum of Live EHS, Nov",
         labels=trt_labels)

# treatment means as bar plots
Meanlivescale_plot_all <- ggplot(data=tmeans_all, 
                                   aes(x=Treatment, y=Mean
                                   ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment Date", y="Mean Live EHS", title="Study 1 - Mean Live EHS on 10 needles of Last Year's Growth") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+  ylim(0, 3)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2") + scale_x_discrete(label = trt_labels)
Meanlivescale_plot_all
ggsave(Meanlivescale_plot_all, file="Meanlivescale_plot_all.pdf", 
       width = 6, height=4)


##


####by county#####
# Treatment survival means for July scale count
tmeans_county_july <- scalecount_july %>% 
  group_by(Treatment, County) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_county_july)

#by county
# Treatment survival means for second (November) scale count
scalecount_nov <- subset(scalecount, grepl('November', Date))
tmeans_county_nov <- scalecount_nov %>% 
  group_by(Treatment, County) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_county_nov)

# Mean live scale plot by county for July
Meanlivescale_plot_county_july <- ggplot(data=tmeans_county_july, 
                                    aes(x=Treatment, y=Mean, fill=County, colour=County
                                    ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Live Count", title="Study 1 - Mean Live Scale By County - July") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+  ylim(0, 5)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2") + scale_x_discrete(label = trt_labels)
Meanlivescale_plot_county_july

# save as a PDF image July collection 
ggsave(Meanlivescale_plot_county, file="Meanlivescale_plot_county.pdf", 
       width = 6, height=4)

# November collection plot, by county
Meanlivescale_plot_county_nov <- ggplot(data=tmeans_county_nov, 
                                     aes(x=Treatment, y=Mean, fill=County, colour=County
                                     ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Live Count", title="Study 1 - Mean Live Scale By County - Nov") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+
  ylim(0, 5)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2") + scale_x_discrete(label = trt_labels)
Meanlivescale_plot_county_nov

# save as a PDF image  - plot of mean live scale, november collection
ggsave(Meanlivescale_plot_county2, file="Meanlivescale_plot_county2.pdf", 
       width = 6, height=4)


##### ##Kruskal wallace, will use: ####
shapiro.test(scalecount_july$Livescale1)
#W = 0.58685, p-value < 2.2e-16 . less than .05, thus data is non normal
#first count:
kruskal.test(Livescale1 ~ Treatment, data = scalecount_july)

pairwise.wilcox.test(scalecount_july$Livescale1, scalecount_july$Treatment,
                     p.adjust.method = "BH")


#second count:
shapiro.test(scalecount_nov$Livescale1)
#W = 0.65372, p-value < 2.2e-16 less than .05, thus data is non normal
kruskal.test(Livescale1 ~ Treatment, data = scalecount_nov)

pairwise.wilcox.test(scalecount_nov$Livescale1, scalecount_nov$Treatment,
                     p.adjust.method = "BH")


###looking at percentages of parasitized ####

#first count:

#extract and create new column with treatment-
scalecountfung_july<-scalecount_july


#remove twigs with no scale
# Note: Ask client about meaning of 0 counts!! 
scalecountfung_july<-scalecountfung_july[rowSums(scalecountfung_july[, 6:11] == 0) < 2,]
scalecountfung_july<-scalecountfung_july %>% drop_na(Label) #some samples only have one twig, dropping those

# Treatment means - percentage of parasitism in July collection
tmeans_july_para <-scalecountfung_july %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Prespara == 1),
    percent_yes100 = round(percent_yes*100),
  )

tmeans_july_para$Collection<-"July"
#checking to make sure we are getting the presence absence mean percentage 
#write.csv(scalecountfung_july,"~/Downloads/scalecountfung_july.csv", row.names = FALSE)

# Convert presence of parasitism to numeric and perform normality check
# definitely not normal, binary data
scalecountfung_july$Prespara<-as.numeric(scalecountfung_july$Prespara)
shapiro.test(scalecountfung_july$Prespara)

#W = 0.58143, p-value < 2.2e-16 . less than .05, thus data is non normal

# perform Kruskal-Wallis test and Wilcox text
#first count:
kruskal.test(Prespara ~ Treatment, data = scalecountfung_july)

# Note - cannot compute exact value due to ties
pairwise.wilcox.test(scalecountfung_july$Prespara, scalecountfung_july$Treatment,
                     p.adjust.method = "BH")

#  1       2     3    4
#  2 0.313 -     -    
#  3 0.161 0.032 -    
#  4 0.766 0.388 0.124


#second count (November):

#extract and create new column with treatment-
scalecountfung2<-scalecount_nov

#remove twigs with no scale and only one shoot
# Remove twigs without scale - if there wasn't EHS this was already counted as 0
scalecountfung2<-scalecountfung2[rowSums(scalecountfung2[, 6:11] == 0) < 2,]
scalecountfung2<-scalecountfung2 %>% drop_na(Label) #some samples only have one twig, dropping those

#percentage of parasitism in a group
tmeans_nov_para<-scalecountfung2 %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Prespara == 1),
    percent_yes100 = round(percent_yes*100),
  )

tmeans_nov_para$Collection<-"Nov"
scalecountfung2$Prespara<-as.numeric(scalecountfung2$Prespara)

# check normality
# Should not be normal data? Binary
shapiro.test(scalecountfung2$Prespara)

#W = 0.58143, p-value < 2.2e-16 . less than .05, thus data is non normal
#first count:
kruskal.test(Prespara ~ Treatment, data = scalecountfung2)

# This one does run! Probably because the 0s were removed?
pairwise.wilcox.test(scalecountfung2$Prespara, scalecountfung2$Treatment,
                     p.adjust.method = "BH")

#1       2     3    
#2 0.256 -     -    
#3 0.021 0.205 -    
#4 0.238 0.885 0.205

# combine parasitic tables
para_table<-dplyr::bind_rows(tmeans_july_para, tmeans_nov_para)

#make plot for parasitism count
# label with treatment dates
trt_labels <- c("May 2-3", "May 17", "May 28", "No Treatment")

para_plot <- plot_means_by_collection(data=para_table, 
                                          title="Study 1 - Parasitism Means Percentages",
                                          x_str="Treatment",
                                          y_str="percent_yes100")
para_plot
ggsave(para_plot, file="Pres_mean_plot.pdf", 
       width = 6, height=10)


####presence of fungi####
#by county

# Treatment survival means - November, by county (not by treatment?)
tmeans_fungus_july <- scalecount_july %>% 
  group_by(County) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )
tmeans_fungus_july$Collection<-"July"

tmeans_fungus_nov <- scalecount_nov %>% 
  group_by(County) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )
tmeans_fungus_nov$Collection<-"Nov"
print(tmeans_fungus_nov)

fungus_table_county <- dplyr::bind_rows(tmeans_fungus_july, tmeans_fungus_nov)

# Plot for county with both July and November
county_labels <- c("Ashe", "Avery", "Mitchell", "Watauga")
plot_means_by_collection(data=fungus_table_county, 
                         title="Study 1 - Mean Presence Percentage of Fungus", 
                         x_str="County",
                         y_str="percent_yes100",
                         labels=county_labels)

# November-only plot for County
tmeans_fungus_plot <- ggplot(data=tmeans_fungus_nov, 
                             aes(x=County, y=percent_yes100
                             ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  labs(x="County", y="Mean Presence %", title="Study 1 - Mean Presence Percentage of Fungus, Nov") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2") 
tmeans_fungus_plot
ggsave(tmeans_fungus_plot, file="tmeans_fungus_plot.pdf", 
       width = 4, height=6)


# Get by treatment
tmeans_fungus_trt_july <- scalecount_july %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )
tmeans_fungus_trt_july$Collection<-"July"

tmeans_fungus_trt_nov <- scalecount_nov %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )
tmeans_fungus_trt_nov$Collection<-"Nov"

fungus_table_trt <- dplyr::bind_rows(tmeans_fungus_trt_july, 
                                        tmeans_fungus_trt_nov)
tmeans_fungus_plot_trt <- plot_means_by_collection(data=fungus_table_trt, 
                                                   title="Study 1 - Mean Presence Percentage of Fungus",
                                                   x_str="Treatment",
                                                   y_str="percent_yes100")
ggsave(tmeans_fungus_plot_trt, file="tmeans_fungus_plot_trt.pdf", 
       width = 4, height=6)


# Fungus count for Mitchell county only
scalecount_novmitchell <- subset(scalecount_nov, grepl('Mitchell',County))
scalecount_novmitchell<-scalecount_novmitchell[rowSums(scalecount_novmitchell[, 6:11] == 0) < 2,]
scalecount_novmitchell<-scalecount_novmitchell %>% drop_na(Label)

tmeans_nov_para_mitchell<-scalecount_novmitchell %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )

scalecount_novmitchell$Prespara<-as.numeric(scalecount_novmitchell$Prespara)

#not significant- no difference in fungus among treatments in mitchell county 
kruskal.test(Prespara ~ Treatment, data = scalecount_novmitchell)
#Kruskal-Wallis chi-squared = 3.9077, df = 3, p-value = 0.2716
pairwise.wilcox.test(scalecount_novmitchell$Prespara, scalecountfung2$Treatment,
                     p.adjust.method = "BH")

###encarsia########

# Encarsia originally only looked at for Nov?
scalecountencar_nov <- scalecount_nov
scalecountencar_july <- scalecount_july


#remove twigs with no scale
scalecountencar_nov <-scalecountencar_nov %>% drop_na(encarsia) #some samples only have one twig, dropping those
scalecountencar_july <- scalecountencar_july %>% drop_na(encarsia)

tmeans_encar_july <- scalecountencar_july %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(encarsia, na.rm = T),3),
    sd = sd(encarsia),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
tmeans_encar_july$Collection<-"July"

tmeans_encar_nov <- scalecountencar_nov %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(encarsia, na.rm = T),3),
    sd = sd(encarsia),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
tmeans_encar_nov$Collection<-"Nov"
print(tmeans_encar)

encar_table_trt <- dplyr::bind_rows(tmeans_encar_july, 
                                     tmeans_encar_nov)

# plot Encarsia means by collection date
tmeans_encar_plot_trt <- plot_means_by_collection(data=encar_table_trt, 
                                                   title="Study 1 - Mean Encarsia Count",
                                                   x_str="Treatment",
                                                   y_str="Mean")

# Encarsia by treatment, November only
scalecountencar_plot <- ggplot(data=tmeans_encar_nov, 
                               aes(x=Treatment, y=Mean
                               ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment Date", y="Mean Count", title="Mean Numbers of Encarsia - November") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")  + scale_x_discrete(label = trt_labels)
scalecountencar_plot
ggsave(scalecountencar_plot, file="Encarsia means.pdf", 
       width = 6, height=4)


# TODO: should convert into function
ggplot(data=scalecountencar_july, aes(x=encarsia, fill=Treatment)) + 
  geom_histogram(binwidth=2) + 
  labs(x="Encarsia Count") + 
  facet_wrap(~Treatment) + 
  scale_fill_discrete(labels=trt_labels) + 
  labs(title=paste("Study 1 - Encarsia Count,",
                   "July"))

ggplot(data=scalecountencar_nov, aes(x=encarsia, fill=Treatment)) + 
  geom_histogram(binwidth=2) + 
  labs(x="Encarsia Count") + 
  facet_wrap(~Treatment) + 
  scale_fill_discrete(labels=trt_labels) + 
  labs(title=paste("Study 1 - Encarsia Count,",
                   "Nov"))



#scalecountfung2$Prespara<-as.numeric(scalecountfung2$Prespara)
shapiro.test(scalecountencar$encarsia)

#W = 0.64279, p-value < 2.2e-16 not normal
kruskal.test(encarsia ~ Treatment, data = scalecountencar)
#Kruskal-Wallis chi-squared = 10.057, df = 3, p-value = 0.01809 significant

pairwise.wilcox.test(scalecountencar$encarsia, scalecountencar$Treatment,
                     p.adjust.method = "BH")
#1         2      3     
#2 0.4788 -      -     
#3 0.3075 0.7669 -     
#4 0.1603 0.1114 0.0065