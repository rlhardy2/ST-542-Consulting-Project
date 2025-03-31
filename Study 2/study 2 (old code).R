## This document includes code from "study 2.R" that was deemed either old or unusable ##
## I didn't want to just outright delete the code! ##

## If you currently have the necessary data sets and such in your environment, this code can be run ##
## Please move unnecessary stuff here so the original file will be more organized ##


#### This is some code that was commented out of the original file

# Removing labels that were not collected
#scalecount1<-scalecount1 %>% drop_na(Livescale1)
#scalecount1<-scalecount1 %>% drop_na(Livescale2)
#scalecount1<-scalecount1 %>% drop_na(Livescale3)

# Figuring out percentage of EHS to cryptomeria 
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

# This is a Poisson model, but the response is live-scale1
gm <- glm(Livescale1 ~ Treatment,
          data = scalecountboth,
          family = "poisson")
summary(gm)

# Not sure what's going on here??
modelg_scale <- glht(gm, mcp(Treatment = "Tukey"))
#table_glht(modelg_scale)

# Pairwise Wilcox test -- gives same error about "ties"
pairwise.wilcox.test(scalecount_Nov$Meanlivescale,
                     scalecount_Nov$Treatment,
                     p.adjust.method = "BH")





