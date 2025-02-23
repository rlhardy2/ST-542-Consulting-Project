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

setwd("~/Documents/nc state/Statistics/2024 Acetyl")
scalecount <- read.csv(file="/Users/jbookwa/Documents/nc state/Statistics/2024 Acetyl/acety counting v2.csv",strip.white=TRUE)
colnames(scalecount) <- c("Label", "Twigab","EHSLivescale1","EHSDeadscale1","EHSLivescale2","EHSDeadscale2",
                          "EHSLivescale3","EHSDeadscale3","Prespara","CryptoLivescale1","CryptoDeadscale1",
                          "CryptoLivescale2","CryptoDeadscale2","CryptoLivescale3","CryptoDeadscale3","Presfungus")
scalecount$Livescale1<-scalecount$EHSLivescale1 + scalecount$CryptoLivescale1
scalecount$Livescale2<-scalecount$EHSLivescale2 + scalecount$CryptoLivescale2
scalecount$Livescale3<-scalecount$EHSLivescale3 + scalecount$CryptoLivescale3
scalecount$Deadscale1<-scalecount$EHSDeadscale1 + scalecount$CryptoDeadscale1
scalecount$Deadscale2<-scalecount$EHSDeadscale2 + scalecount$CryptoDeadscale2
scalecount$Deadscale3<-scalecount$EHSDeadscale3 + scalecount$CryptoDeadscale2

scalecount<-subset(scalecount, select=-c(EHSLivescale1,EHSLivescale2,EHSLivescale3,EHSDeadscale1,
                                         EHSDeadscale2,EHSDeadscale3,CryptoLivescale1,
                                         CryptoLivescale2,CryptoLivescale3,CryptoDeadscale1,CryptoDeadscale2,
                                         CryptoDeadscale3))

#getting zero and 1s as pres and abs
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="yes", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="y", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="Y", 1)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="n", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="N", 0)
scalecount$Prespara<-replace(scalecount$Prespara, scalecount$Prespara=="no", 0)

scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="yes", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="y", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="Y", 1)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="n", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="N", 0)
scalecount$Presfungus<-replace(scalecount$Presfungus, scalecount$Presfungus=="no", 0)


scalecount <- scalecount %>% 
  mutate(across(c(Livescale1, Deadscale1, Livescale2, Deadscale2, Livescale3, Deadscale3), as.numeric))

str(scalecount)


#getting new column 
scalecount<-scalecount %>% 
  mutate(Meanlivescale = rowMeans(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))

scalecount<-scalecount %>% 
  mutate(Meandeadscale = rowMeans(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

scalecount<-scalecount %>% 
  mutate(Sumlivescale = rowSums(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))

scalecount<-scalecount %>% 
  mutate(Sumdeadscale = rowSums(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

#removing labels that were not collected
scalecount<-scalecount %>% drop_na(Livescale1)

str(scalecount)# 
write.csv(scalecount,"/Users/jbookwa/Documents/nc state/Statistics/2024 Acetyl/scalecount.csv", row.names = FALSE)
#this is all data. to get replicates at shoots, delete columns

scalecountall <- scalecount[ -c(5:10) ]


#extract and create new column with treatment
scalecountall$Treatment<-
  substring(scalecountall$Label, first=2, last=2)

scalecountall <- scalecountall %>% 
  mutate(across(c(Treatment), as.factor))

write.csv(scalecountall,"/Users/jbookwa/Documents/nc state/Statistics/2024 Acetyl/scalecountall.csv", row.names = FALSE)


# Treatment survival means (mean of the three shoots per sample)
tmeans_mean <- scalecountall %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(Meanlivescale, na.rm = T),3),
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
    Mean=round(mean(Sumlivescale, na.rm = T),3),
    sd = sd(Sumlivescale),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans_sum)




###### treatment means as bar plots- use this one####
labels <- c("Acetamiprid","Dinotefuron","Sulfoxaflor", "Flupyradifurone","No Treatment" )
tmeans_mean$Treatment <- factor(tmeans_mean$Treatment, levels=c("2","1","3","5","4"))
livescale_plot_mean <- ggplot(data=tmeans_mean, 
                         aes(x=Treatment, y=Mean
                         ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Numbers of Live EHS on 10 needles of This Year's Growth") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+scale_x_discrete(label = labels)+
scale_fill_manual(values = c("#CC6600","red","#9999FF","#0072B2","green4"))  
livescale_plot_mean


#### treatment means as bar plots don't use ####
livescale_plots_sum <- ggplot(data=tmeans_sum, 
                         aes(x=Treatment, y=Mean
                         ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), position = position_dodge(0.9), width = 0,
                show.legend = FALSE, color="black") +
  labs(x="Treatment", y="Mean Numbers of Live EHS on 10 needles of This Year's Growth") +
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
  scale_fill_manual(values = c("#CC6600","red","#9999FF","#0072B2","green4")) 
livescale_plots_sum
######

#Shapiro-Wilk normality test
shapiro.test(scalecountall$Meanlivescale)
#p value lower than .05, therefore non normal distribution. use non parametric test


library(biostat3)
##this is what we are using
kruskal.test(Meanlivescale ~ Treatment, data = scalecountall)

letters<-pairwise.wilcox.test(scalecountall$Meanlivescale, scalecountall$Treatment,
                     p.adjust.method = "BH")
letters
library(multcompView)
comparison_letters <- multcompLetters(letters$p.value) 

comparison_letters


#1   dinote    
#2    acetampirid  
#3    sulflo 
#4    no app 
#5.    flupyradifurone


#1          2       3       4      
#2 0.47865 -       -       -      
#3 0.25351 0.05070 -       -      
#4 1.4e-06 1.9e-07 2.6e-05 -      
#5 0.00231 0.00015 0.03909 0.05070

multcompLetters2(letters)


#error due to ties, ties are thrown out. 

#looking at poisson as well
poisson.model <- glm(Sumlivescale ~ Treatment,
                     family = 'poisson', data = scalecountall)
summary(poisson.model)



###looking at percentages of fungus####
scalecountall <- scalecountall %>% 
  mutate(across(c(Prespara), as.numeric))

#percentage of parasitism in a treatment
tmeans_para <- scalecountall %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(Prespara, na.rm = T),2),
    sd = sd(Prespara),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans)

#now with no zeros- taking out all samples in which no EHS were found#####
scalecount_nozero<-scalecountall

scalecount_nozero<-scalecount_nozero[rowSums(scalecount_nozero[, 3:6] == 0) < 2,]

#percentage of presence in a group

tmeans1<-scalecount_nozero %>% 
  group_by(Treatment) %>% 
  summarise(
    percent_yes = mean(Prespara == 1),
    percent_yes100 = round(percent_yes*100),
  )


kruskal.test(Prespara ~ Treatment, data = scalecount_nozero)

pairwise.wilcox.test(scalecount_nozero$Prespara, scalecount_nozero$Treatment,
                     p.adjust.method = "BH")

#### treatment means as bar plots for parasitism no zeros####
labels <- c("Acetamiprid","Dinotefuron","Sulfoxaflor", "Flupyradifurone","No Treatment" )
tmeans1$Treatment <- factor(tmeans1$Treatment, levels=c("2","1","3","5","4"))
livescale_plots_para <- ggplot(data=tmeans1,
                              aes(x=Treatment, y=percent_yes
                              ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +

  labs(x="Treatment", y="Percentage of parasitism on 10 needles of This Year's Growth") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+scale_x_discrete(label = labels)
livescale_plots_para
ggsave(livescale_plots_para, file="livescale_plots_para.pdf", 
       width = 6, height=4)






#####percentage of fungus in a treatment plot no zeros####
scalecount_nozero <- scalecount_nozero %>% 
  mutate(across(c(Presfungus), as.numeric))


tmeans_fungus <- scalecount_nozero %>% 
  group_by(Treatment) %>% 
  dplyr::summarize(
    Mean=round(mean(Presfungus, na.rm = T),2),
    sd = sd(Presfungus),
    n = n(),
    se = sd / sqrt(n),
    cv = sd/Mean * 100
  )
print(tmeans)


kruskal.test(Presfungus ~ Treatment, data = scalecount_nozero)

pairwise.wilcox.test(scalecount_nozero$Presfungus, scalecount_nozero$Treatment,
                     p.adjust.method = "BH")

tmeans_fungus$Treatment <- factor(tmeans_fungus$Treatment, levels=c("2","1","3","5","4"))
livescale_plots_fungus <- ggplot(data=tmeans_fungus, 
                               aes(x=Treatment, y=Mean
                               ), na.rm = T) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"))  +
  labs(x="Treatment", y="Percentage of entofungus on 10 needles of this Year's Growth") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle = 35, vjust = 0.5, size=18), 
        axis.title.x = element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)
  )+scale_x_discrete(label = labels)
  
livescale_plots_fungus
ggsave(livescale_plots_fungus, file="livescale_plots_fungus.pdf", 
       width = 6, height=4)



