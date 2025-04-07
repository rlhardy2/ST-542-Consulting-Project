## This document includes code from "study 1.R" that was deemed either old or unusable ##
## I didn't want to just outright delete the code! ##

## If you currently have the necessary data sets and such in your environment, this code can be run ##
## Please move unnecessary stuff here so the original file will be more organized ##



#### Old code for Kruskal-Wallis and Wilcoxon tests below

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


# scalecount_para_july <- scalecount_july[rowSums(scalecount_july[, 12:13] == 0) < 2,]
# scalecount_para_july <- 
#   scalecount_para_july %>% 
#   drop_na(Label) #some samples only have one twig
# 
# scalecount_para_nov <- scalecount_nov[rowSums(scalecount_nov[, 12:13] == 0) < 2,]
# scalecount_para_nov<- 
#   scalecount_para_nov %>% 
#   drop_na(Label) #some samples only have one twig

#### Kruskal-Wallis + Dunn

# Since the blocks were not truly specified (randomly assigned)
# we may just be able to collapse this to one factor if 
# we use the by-tree data set?
# However, blocks seem to have quite an effect in Study 1 
# (may be bc of diff counties)

# Insignificant results for July?
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_tree_mean_july)

# Significant results for November
kruskal.test(Meanlivescale ~ Treatment, data = scalecount_tree_mean_nov)
dunn_nov <- dunnTest(Meanlivescale ~ Treatment, 
                     data = scalecount_tree_mean_nov, method = "bh")

print(dunn_nov, dunn.test.results=TRUE)



