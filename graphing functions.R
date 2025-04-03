library(dplyr)
library(tidyr)
library(ggplot2)

# plot means by collection date and an X value (treatment, county)
# "data" should have the tmeans format illustrated in Study 1
plot_means_by_collection <- function(data, title, x_str, y_str, labels=trt_labels, y_lab=y_str) {
  x_sym <- ensym(x_str)
  y_sym <- ensym(y_str)
  trt_coll_plot <- ggplot(data=data, 
                          aes(x={{x_sym}}, y={{y_sym}}, 
                              fill=Collection, colour=Collection
                          ), na.rm = T) +
    geom_bar(stat="identity", 
             position = position_dodge2(width = 0.9, preserve = "single"))  +
    labs(x=x_str, y=y_lab, title=title) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x=element_text(angle = 45, vjust = 0.5, size=15), 
          axis.title.x = element_text(size=18),
          axis.title.y=element_text(size=18),
          axis.text.y=element_text(size=18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18)
    )+
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2")  + scale_x_discrete(label = labels)
  
  return(trt_coll_plot)
}

# Histograms - one plot for each treatment
get_hist <- function(data, x_str, x_lab, title, labels, binwidth=1) {
  x_sym <- ensym(x_str)
  hist <- ggplot(data=data, aes(x={{x_sym}}, fill=Treatment)) + 
    geom_histogram(binwidth=binwidth) + 
    labs(x=x_lab) + 
    facet_wrap(~Treatment) + 
    scale_fill_discrete(labels=labels) + 
    labs(title=title)
  return(hist)
}

# Histogram of mean livescale
# outputs one plot for each treatment
get_hist_livescale <- function(data, collection_date, study, labels) {
  title <- paste("Study", study, "- Mean Numbers of Live Scale Insects,", collection_date)
  return (get_hist(data=data, x_str="Meanlivescale", x_lab="Mean Live Scale",
                  title=title, labels=labels))
}

get_hist_encarsia <- function(data, collection_date, study, labels) {
  title <- paste("Study", study, "- Encarsia Count,", collection_date)
  return (get_hist(data=data, x_str="encarsia", 
                   x_lab="Encarsia Count", title=title, labels=labels))
}

# Histogram, all treatments merged
# Change bin width as needed...
get_hist_all_trt <- function(data, x_str, x_lab, title, labels, binwidth=1) {
  x_sym <- ensym(x_str)
  hist <- ggplot(data=data, aes(x={{x_sym}})) + 
    geom_histogram(binwidth=binwidth) + 
    labs(x=x_lab) + 
    scale_fill_discrete(labels=labels) + 
    labs(title=title)
  return(hist)
}

get_cis_marginal_means_plot <- function(ci_df, pairs_df, y_positions, trt_labels,
                                        y_str, y_lab, title) {
  y_sym <- ensym(y_str)
  pairs_df_graphing <- pairs_df
  
  # Extracts treatments from pairs for graphing
  # Must be called group1 and group2 for stat_pvalue_manual
  pairs_df_graphing$group1 <- substring(pairs_df$contrast, first=10, last=10)
  pairs_df_graphing$group2 <- substring(pairs_df$contrast, first=23, last=23)
  pairs_df_graphing$p.value_round <- round(pairs_df$p.value, digits=5)
  # y position of each pairwise comparison bar
  pairs_df_graphing$y.position <- y_positions
  
  emmeans_ci_plot <- 
    ggplot(ci_df, aes(x = Treatment, y = {{y_sym}})) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
    labs(title = title,
         x = "Treatment",
         y = y_lab) +
    # Add treatment means comparison bars
    stat_pvalue_manual(
      data = pairs_df_graphing, label = "p.value_round",
      xmin = "group1", xmax = "group2",
      y.position = "y.position"
    ) +
    scale_x_discrete(label = trt_labels)
  return (emmeans_ci_plot)
}
  
