## This document includes code from "study 3.R" that was deemed either old or unusable ##
## I didn't want to just outright delete the code! ##

## If you currently have the necessary data sets and such in your environment, this code can be run ##
## Please move unnecessary stuff here so the original file will be more organized ##


#### This is some code that was commented out of the original file

# THIS DOESN'T WORK -- gives error message about ties
#letters < -pairwise.wilcox.test(scalecountall$Meanlivescale, scalecountall$Treatment,
#                               p.adjust.method = "BH")
#letters

#comparison_letters <- multcompLetters(letters$p.value)
#comparison_letters

#1 Dinotefuran   
#2 Acetampirid  
#3 Sulfoxaflor
#4 No Treatment
#5 Flupyradifurone

#1     1      2       3       4   
#2 0.47865 -       -       -      
#3 0.25351 0.05070 -       -      
#4 1.4e-06 1.9e-07 2.6e-05 -      
#5 0.00231 0.00015 0.03909 0.05070

#multcompLetters2(letters)

# Error due to ties, ties are thrown out ??

# Looking at poisson as well
#poisson.model <- glm(Sumlivescale ~ Treatment,
#                     family = 'poisson', data = scalecountall)
#summary(poisson.model)






