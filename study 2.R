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

setwd("~/Documents/nc state/Statistics/2024 Distance")
scalecount1 <- read.csv(file="/Users/jbookwa/Documents/nc state/Statistics/2024 Distance/distance numbers v7.csv",strip.white=TRUE)
colnames(scalecount1) <- c("Label","Type", "Twigab","Date","Livescale1","Deadscale1","Livescale2","Deadscale2",
                          "Livescale3","Deadscale3","Prespara","Presscalenewgr","Presfungus","encarsia","notes")
#just getting columns I need
scalecount1 = subset(scalecount1, select = c("Label","Type", "Twigab","Date","Livescale1","Deadscale1","Livescale2","Deadscale2",
                                           "Livescale3","Deadscale3","Prespara","Presfungus","encarsia"))

#extract and create new column with treatment
scalecount1$Treatment<-
  substring(scalecount1$Label, first=4, last=4)

scalecount1<-scalecount1 %>% drop_na(Livescale1)


#getting zero and 1s as pres and abs
scalecount1$Prespara<-replace(scalecount1$Prespara, scalecount1$Prespara=="yes", 1)
scalecount1$Prespara<-replace(scalecount1$Prespara, scalecount1$Prespara=="no", 0)
scalecount1$Presfungus<-replace(scalecount1$Presfungus, scalecount1$Presfungus=="yes", 1)
scalecount1$Presfungus<-replace(scalecount1$Presfungus, scalecount1$Presfungus=="no", 0)


scalecount1 <- scalecount1 %>% 
  mutate(across(c(Livescale1, Deadscale1, Livescale2, Deadscale2, Livescale3, Deadscale3), as.numeric))
str(scalecount1)

#removing labels that were not collected
#scalecount1<-scalecount1 %>% drop_na(Livescale1)
#scalecount1<-scalecount1 %>% drop_na(Livescale2)
#scalecount1<-scalecount1 %>% drop_na(Livescale3)

#figuring out percentage of EHS to cryptomeria 
#scalecount1<-scalecount1 %>% 
#  mutate(SUMlivescale = rowSums(across(c(Livescale1,Livescale2,Livescale3)),na.rm = TRUE))

#tmeans <- scalecount1 %>% 
#  group_by(Type) %>% 
#  dplyr::summarize(
#    Mean=round(mean(Livescale1, na.rm = T),3),
#    sd = sd(Livescale1),
#    n = n(),
#    se = sd / sqrt(n),
#   cv = sd/Mean * 100
#  )
#print(tmeans)


#making a new table that adds the EHS and cryptomeria 
scalecountboth<- scalecount1 %>% 
  group_by(Label,Twigab, Date, Treatment) %>% 
  dplyr::summarize(across(where(is.numeric), sum))


#getting new column 
#get another dataframe to send to no NA- this is the column we will use
scalecountboth_NoNA<-scalecountboth
scalecountboth$Treatment<-as.factor(scalecountboth$Treatment)
scalecountboth$Twigab<-as.factor(scalecountboth$Twigab)

####must use across NOT SELECT - this is for getting mean scale#####
scalecountboth_NoNA<-scalecountboth_NoNA %>% 
  mutate(Meanlivescale = rowMeans(across(c(Livescale1,Livescale2,Livescale3)),na.rm = TRUE))


scalecountboth_NoNA<-scalecountboth_NoNA %>% 
  mutate(Meandeadscale = rowMeans(across(c(Deadscale1, Deadscale2, Deadscale3)), na.rm = TRUE))




#scalecountboth<-as.data.frame(scalecountboth)
#str(scalecountboth)
#gm <- glm(Livescale1 ~ Treatment , 
#            data = scalecountboth,
#             family = "poisson")                                                         
#summary(gm)  
########


#gm <- glm(Livescale1 ~ Treatment , 
#          data = scalecountboth,
#          family = "poisson")                                                         
#summary(gm)  


scalecount_July <- subset(scalecountboth_NoNA, grepl('July',Date))
scalecount_Nov <- subset(scalecountboth_NoNA, grepl('Nov',Date))

#####July kruskal#####
kruskal.test(Livescale1 ~ Treatment, data = scalecount_July)
  #Kruskal-Wallis chi-squared = 3.3811, df = 3, p-value = 0.3365

#modelg_scale<-glht(gm, mcp(Treatment="Tukey"))
#table_glht(modelg_scale)



###Nov kruskal
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_Nov)

#Kruskal-Wallis chi-squared = 17.059, df = 3, p-value = 0.0006873
pairwise.wilcox.test(scalecount_Nov$Meanlivescale, scalecount_Nov$Treatment,
                     p.adjust.method = "BH")

         #1      2        3     
#2     0.0087   -         -     
#3    0.0159   0.3992   -     
#4      0.7103 0.0086.  0.0104

tmeans <- scalecount_Nov %>% 
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
labels <- c("Pyri/May", "Ace/May", "Pyri/May + Ace/June", "No Treatment")
Meanlivescale_plot_Nov <- ggplot(data=tmeans, 
                             aes(x=Treatment, y=Mean
                             ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Numbers of Live EHS & Crypto Scale on 10 needles") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+
  scale_fill_brewer(palette = "Dark2")+ scale_x_discrete(label = labels)
Meanlivescale_plot_Nov
ggsave(Meanlivescale_plot_Nov, file="Meanlivescale_plot_Nov.pdf", 
       width = 6, height=4)



##encarsia####
scalecountencar<-scalecount_Nov


#remove twigs with no scale

scalecountencar<-scalecountencar %>% drop_na(encarsia) #some samples only have one twig, dropping those

kruskal.test(encarsia ~ Treatment, data = scalecountencar)
#Kruskal-Wallis chi-squared = 2.0473, df = 3, p-value = 0.5626
#no difference among treatments





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


# treatment means as bar plots
labels <- c("Pyri/May", "Ace/May", "Pyri/May + Ace/June", "No Treatment")
encarsia_plot<- ggplot(data=tmeans, 
                                 aes(x=Treatment, y=Mean
                                 ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean numbers of Encarsia ") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+
  scale_fill_brewer(palette = "Dark2")+ scale_x_discrete(label = labels)
encarsia_plot
ggsave(encarsia_plot, file="encarsia_plot.pdf", 
       width = 6, height=4)







##fungus#########

#extract and create new column with treatment-
scalecountfung2<-scalecount_Nov


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



