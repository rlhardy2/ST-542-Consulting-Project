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

setwd("~/Documents/nc state/Statistics/2024 Acetamiprid")
scalecount <- read.csv(file="/Users/jbookwa/Documents/nc state/Statistics/2024 Acetamiprid/EHS count 2024 v7.csv",strip.white=TRUE)
colnames(scalecount) <- c("Label", "County","Twigab","Date","Counter","Livescale1","Deadscale1","Livescale2","Deadscale2",
                          "Livescale3","Deadscale3","Prespara","Presfungus","Presscalenewgr","encarsia")

#delete columns with no info
scalecount <- scalecount[ -c(16:19) ]

#extract and create new column with treatment
scalecount$Treatment<-
  substring(scalecount$Label, first=4, last=4)

#getting zero and 1s as pres and abs
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="yes", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="y", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="Y", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="n", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="N", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="no", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="yes", 1)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="y", 1)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="Y", 1)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="n", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="N", 0)
scalecount$Presscalenewgr<-replace(scalecount$Presscalenewgr, scalecount$Presscalenewgr=="no", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="yes", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="y", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="Y", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="n", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="N", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="no", 0)


scalecount <- scalecount %>% 
  mutate(across(c(Treatment, County), as.factor))

scalecount <- scalecount %>% 
 mutate(across(c(Livescale1, Deadscale1, Livescale2, Deadscale2, Livescale3, Deadscale3), as.numeric))


str(scalecount)



#getting new column 
scalecount<-scalecount %>% 
  mutate(Meanlivescale = rowMeans(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))

scalecount<-scalecount %>% 
  mutate(Meandeadscale = rowMeans(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

#removing labels that were not collected
scalecount<-scalecount %>% drop_na(Livescale1)

str(scalecount)


#scalecount has both July and November counts in it

#getting first count- July 

scalecount1 <- subset(scalecount, grepl('July',Date))

# Treatment survival means
tmeans <- scalecount1 %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans)


# treatment means as bar plots
labels <- c("May 2-3", "May 17", "May 28", "No Treatment")
Meanlivescale_plot_first <- ggplot(data=tmeans, 
                              aes(x=Treatment, y=Mean
                              ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment Date", y="Mean Numbers of Live EHS on 10 needles of Last Year's Growth") +
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
  scale_color_brewer(palette = "Dark2") + scale_x_discrete(label = labels)
Meanlivescale_plot_first
ggsave(Meanlivescale_plot_first, file="Meanlivescale_plot_first.pdf", 
       width = 6, height=4)


##






####by county#####
# Treatment survival means
tmeans <- scalecount1 %>% 
  group_by(Treatment, County) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans)

#by county
# Treatment survival means
tmeans_county2 <- scalecount2 %>% 
  group_by(Treatment, County) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_county2)


labels <- c("May 2-3", "May 17", "May 28", "No Treatment")
Meanlivescale_plot_county <- ggplot(data=tmeans, 
                     aes(x=Treatment, y=Mean, fill=County, colour=County
                     ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Live Count") +
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
  scale_color_brewer(palette = "Dark2") + scale_x_discrete(label = labels)
Meanlivescale_plot_county

# save as a PDF image July collection 
ggsave(Meanlivescale_plot_county, file="Meanlivescale_plot_county.pdf", 
       width = 6, height=4)

#november collection
labels <- c("May 2-3", "May 17", "May 28", "No Treatment")
Meanlivescale_plot_county2 <- ggplot(data=tmeans_county2, 
                                    aes(x=Treatment, y=Mean, fill=County, colour=County
                                    ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Live Count 2") +
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
  scale_color_brewer(palette = "Dark2") + scale_x_discrete(label = labels)
Meanlivescale_plot_county2

# save as a PDF image  november collection
ggsave(Meanlivescale_plot_county2, file="Meanlivescale_plot_county2.pdf", 
       width = 6, height=4)













## working on stats, do not use####

p <- ggplot(scalecount1, aes(x=Meanlivescale)) +
  geom_histogram(bins=10)
p

qqnorm(scalecount1$Meanlivescale, pch = 1, frame = FALSE)
#distribution squewed to the right
shapiro.test(scalecount1$Meanlivescale)
#the null hypothesis that the data follow a normal distribution. this result says
#we can reject the null.. so data is not normal


library(coxed)
library(DHARMa)
source("functions 2024 acetamiprid.r")

# To store the results.
res <- list()


# Regressions with a truncated-Poisson (1), truncated negative binomial (2) and (3) and
# generalized Poisson (4) distributional models.
res[[1]] <- twig_function(scalecount1, 1) 
res[[2]] <- twig_function(scalecount1, 2)
res[[3]] <- twig_function(scalecount1, 3)
res[[4]] <- twig_function(scalecount1, 4)

summary(res[[1]])
summary(res[[2]])
summary(res[[3]])
summary(res[[4]])

# Let's first have a look at how well those models perform. First, we calculate the Pearson, Kendall
# and Spearman correlation coefficients. The higher, the better.
print(c('Truncated Poisson'=cor(predict(res[[1]],type="response"),scalecount1$Meanlivescale),
        'Truncated negative binomial 1'=cor(predict(res[[2]],type="response"),scalecount1$Meanlivescale),
        'Truncated negative binomial 2'=cor(predict(res[[3]],type="response"),scalecount1$Meanlivescale),
        'Truncated generalized Poisson'=cor(predict(res[[4]],type="response"),scalecount1$Meanlivescale)))
#all about the same, maybe truncated generalized poisson is the best

print(c('Truncated Poisson'=cor(predict(res[[1]],type="response"),scalecount1$Meanlivescale,method="kendall"),
        'Truncated negative binomial 1'=cor(predict(res[[2]],type="response"),scalecount1$Meanlivescale,method="kendall"),
        'Truncated negative binomial 2'=cor(predict(res[[3]],type="response"),scalecount1$Meanlivescale,method="kendall"),
        'Truncated generalized Poisson'=cor(predict(res[[4]],type="response"),scalecount1$Meanlivescale,method="kendall")))
# Truncated Poisson is best
print(c('Truncated Poisson'=cor(predict(res[[1]],type="response"),scalecount1$Meanlivescale,method="spearman"),
        'Truncated negative binomial 1'=cor(predict(res[[2]],type="response"),scalecount1$Meanlivescale,method="spearman"),
        'Truncated negative binomial 2'=cor(predict(res[[3]],type="response"),scalecount1$Meanlivescale,method="spearman"),
        'Truncated generalized Poisson'=cor(predict(res[[4]],type="response"),scalecount1$Meanlivescale,method="spearman")))
# Truncated Poisson is best


## Next, we calculate the root-mean-square-error and the mean-absolute deviation. #JAMIE: Both MAE and RMSE 
#express average model prediction error in units of the variable of interest. 
#Both metrics can range from 0 to ∞ and are indifferent to the direction of errors. 
#They are negatively-oriented scores, which means lower values are better.


# Let's have a look at a plot of predicted vs. observed values. We use log-transformed axes.
par(mfrow=c(2,2))
plot(predict(res[[1]],type="conditional"),scalecount1$Meanlivescale,log="xy")
points(c(1,1000),c(1,1000),type="l")
plot(predict(res[[2]],type="conditional"),scalecount1$Meanlivescale,log="xy")
points(c(1,1000),c(1,1000),type="l")
plot(predict(res[[3]],type="conditional"),scalecount1$Meanlivescale,log="xy")
points(c(1,1000),c(1,1000),type="l")
plot(predict(res[[4]],type="conditional"),scalecount1$Meanlivescale,log="xy")
points(c(1,1000),c(1,1000),type="l")
par(mfrow=c(1,1))



# Residuals with the DHARMa package.
sim_res <- lapply(res,function(x) DHARMa::simulateResiduals(x, n = 1000))
plot(sim_res[[1]])
plot(sim_res[[2]])
plot(sim_res[[3]])
plot(sim_res[[4]])


############second count#####


scalecount2 <- subset(scalecount, grepl('November',Date))


# Treatment survival means
tmeans <- scalecount2 %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
    sd = sd(Meanlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans)


# treatment means as bar plots
labels <- c("May 2-3", "May 17", "May 28", "No Treatment")
Meanlivescale_plot_second <- ggplot(data=tmeans, 
                             aes(x=Treatment, y=Mean
                             ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment Date", y="Mean Numbers of Live EHS on 10 needles of This Year's Growth in November") +
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
  scale_color_brewer(palette = "Dark2")  + scale_x_discrete(label = labels)
Meanlivescale_plot_second
ggsave(Meanlivescale_plot_second, file="Meanlivescale_plot_second.pdf", 
       width = 6, height=4)




##### ##Kruskal wallace, will use: ####
shapiro.test(scalecount1$Livescale1)
#W = 0.58685, p-value < 2.2e-16 . less than .05, thus data is non normal
#first count:
kruskal.test(Livescale1 ~ Treatment, data = scalecount1)

pairwise.wilcox.test(scalecount1$Livescale1, scalecount1$Treatment,
                     p.adjust.method = "BH")


#second count:
shapiro.test(scalecount2$Livescale1)
#W = 0.65372, p-value < 2.2e-16 less than .05, thus data is non normal
kruskal.test(Livescale1 ~ Treatment, data = scalecount2)

pairwise.wilcox.test(scalecount2$Livescale1, scalecount2$Treatment,
                     p.adjust.method = "BH")









###looking at percentages of parasitized ####

#first count:

#extract and create new column with treatment-
scalecountfung1<-scalecount1


#remove twigs with no scale

scalecountfung1<-scalecountfung1[rowSums(scalecountfung1[, 6:11] == 0) < 2,]
scalecountfung1<-scalecountfung1 %>% drop_na(Label) #some samples only have one twig, dropping those

#percentage of presence in a group
tmeans1<-scalecountfung1 %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Prespara == 1),
    percent_yes100 = round(percent_yes*100),
  )

tmeans1$Collection<-"July"
#checking to make sure we are getting the presence absence mean percentage 
#write.csv(scalecountfung1,"~/Downloads/scalecountfung1.csv", row.names = FALSE)

scalecountfung1$Prespara<-as.numeric(scalecountfung1$Prespara)
shapiro.test(scalecountfung1$Prespara)

#W = 0.58143, p-value < 2.2e-16 . less than .05, thus data is non normal
#first count:
kruskal.test(Prespara ~ Treatment, data = scalecountfung1)

pairwise.wilcox.test(scalecountfung1$Prespara, scalecountfung1$Treatment,
                     p.adjust.method = "BH")

#  1       2     3    4
#  2 0.313 -     -    
#  3 0.161 0.032 -    
#  4 0.766 0.388 0.124

####parastism just in July plot not used####
labels <- c("May 2-3", "May 17", "May 28", "No Treatment")
scalecountfung1_para_plot <- ggplot(data=tmeans, 
                                    aes(x=Treatment, y=percent_yes100
                                    ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment Date", y="Mean Numbers of Live EHS on 10 needles of This Year's Growth in November") +
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
  scale_color_brewer(palette = "Dark2")  + scale_x_discrete(label = labels)
Meanlivescale_plot_second
ggsave(Meanlivescale_plot_second, file="Parasitism_plot_first.pdf", 
       width = 6, height=4)

#####


#second count:

#extract and create new column with treatment-
scalecountfung2<-scalecount2


#remove twigs with no scale

scalecountfung2<-scalecountfung2[rowSums(scalecountfung2[, 6:11] == 0) < 2,]
scalecountfung2<-scalecountfung2 %>% drop_na(Label) #some samples only have one twig, dropping those

#percentage of presence in a group
tmeans2<-scalecountfung2 %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Prespara == 1),
    percent_yes100 = round(percent_yes*100),
  )

tmeans2$Collection<-"Nov"
scalecountfung2$Prespara<-as.numeric(scalecountfung2$Prespara)
shapiro.test(scalecountfung2$Prespara)

#W = 0.58143, p-value < 2.2e-16 . less than .05, thus data is non normal
#first count:
kruskal.test(Prespara ~ Treatment, data = scalecountfung2)

pairwise.wilcox.test(scalecountfung2$Prespara, scalecountfung2$Treatment,
                     p.adjust.method = "BH")

#1       2     3    
#2 0.256 -     -    
#3 0.021 0.205 -    
#4 0.238 0.885 0.205


#combine parasitic tables


para_table<-dplyr::bind_rows(tmeans1, tmeans2)

#make plot for both parasitic count
labels <- c("May 2-3", "May 17", "May 28", "No Treatment")

para_plot <- ggplot(data=para_table, 
                                    aes(x=Treatment, y=percent_yes100, fill=Collection, colour=Collection
                                    ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  labs(x="Treatment", y="Presence Absence Mean Percentage") +
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
  scale_color_brewer(palette = "Dark2")  + scale_x_discrete(label = labels)
para_plot

ggsave(para_plot, file="Pres_mean_plot.pdf", 
       width = 6, height=10)















####presence of fungi####
#by county

# Treatment survival means
tmeans_fungus <- scalecount2 %>% 
  group_by(County) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )
print(tmeans_fungus)

tmeans_fungus_plot <- ggplot(data=tmeans_fungus, 
                             aes(x=County, y=percent_yes100
                             ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  labs(x="County", y="Mean Presence Percentage of Fungus Infestation") +
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




scalecount2mitchell <- subset(scalecount2, grepl('Mitchell',County))
scalecount2mitchell<-scalecount2mitchell[rowSums(scalecount2mitchell[, 6:11] == 0) < 2,]
scalecount2mitchell<-scalecount2mitchell %>% drop_na(Label)

tmeans2<-scalecount2mitchell %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Presfungus == 1),
    percent_yes100 = round(percent_yes*100),
  )

scalecount2mitchell$Prespara<-as.numeric(scalecount2mitchell$Prespara)

#not significant- no difference in fungus among treatments in mitchell county 
kruskal.test(Prespara ~ Treatment, data = scalecount2mitchell)
#Kruskal-Wallis chi-squared = 3.9077, df = 3, p-value = 0.2716
pairwise.wilcox.test(scalecount2mitchell$Prespara, scalecountfung2$Treatment,
                     p.adjust.method = "BH")

######encarsia########

scalecountencar<-scalecount2


#remove twigs with no scale

scalecountencar<-scalecountencar %>% drop_na(encarsia) #some samples only have one twig, dropping those

tmeans <- scalecountencar %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(encarsia, na.rm = T),3),
    sd = sd(encarsia),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans)


labels <- c("May 2-3", "May 17", "May 28", "No Treatment")
scalecountencar_plot <- ggplot(data=tmeans, 
                                    aes(x=Treatment, y=Mean
                                    ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment Date", y="Mean Numbers of Encarsia") +
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
  scale_color_brewer(palette = "Dark2")  + scale_x_discrete(label = labels)
scalecountencar_plot
ggsave(scalecountencar_plot, file="Encarsia means.pdf", 
       width = 6, height=4)


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



