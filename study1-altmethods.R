## working on stats, do not use####
# These are the alternative models that Jamie was pursuing
# Truncated Poisson, truncated negative binomial, truncated generalized Poisson

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


