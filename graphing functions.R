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
get_hist <- function(data, x_str, x_lab, title, labels) {
  x_sym <- ensym(x_str)
  hist_livescale <- ggplot(data=data, aes(x={{x_sym}}, fill=Treatment)) + 
    geom_histogram(binwidth=1) + 
    labs(x=x_lab) + 
    facet_wrap(~Treatment) + 
    scale_fill_discrete(labels=labels) + 
    labs(title=title)
  return(hist_livescale)
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

