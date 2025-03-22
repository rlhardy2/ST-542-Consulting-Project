library(dplyr)
library(stringr)

# Check for na in live, dead scale counts
check_na_scale_counts <- function(scalecount) {
  has_na <- scalecount %>% 
    filter(if_any(c(Livescale1, Livescale2, Livescale3, 
                    Deadscale1, Deadscale2, Deadscale3), is.na))
  return (has_na)
}

# check counts by label and twig 
# Should be unique if all vals from dataset within same timeframe (July, Nov)
check_counts <- function(scalecount) {
  count <- scalecount %>% group_by(Label, Twigab) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count))
  return (count)
}

# Label denotes the same tree
# Average count data (live scale, dead scale, encarsia) 
# across both twigs to get by tree
average_counts_across_twigs <- function(scalecount) {
  scalecount_tree_mean <-
    scalecount %>%
    group_by(Label, Date, Treatment) %>% 
    summarize(across(where(is.numeric), mean)) %>%
    # drop rows that shouldn't use means
    dplyr::select(-one_of(c('Sumlivescale', 'Sumdeadscale')))
  
  return (scalecount_tree_mean)
}

sum_counts_across_twigs <- function(scalecount) {
  scalecount_tree_sum <-
    scalecount %>%
    group_by(Label, Date, Treatment) %>% 
    summarize(across(where(is.numeric), sum)) %>%
    # drop rows that shouldn't use sums
    dplyr::select(-one_of(c('encarsia', 'Meanlivescale', 'Meandeadscale')))
  
  return (scalecount_tree_sum)
}

# replace binary strings with 0s and 1s
replace_strings_with_binary <- function(vars) {
  vars_with_nums <- replace(vars, vars %in% c("yes", "y", "Y"), 1)
  
  # This is the criteria for both study 1 and study 3
  # study 2 doesn't include "", but there's no data with a blank anyway
  # that isn't also missing Livescale1
  vars_with_nums <- replace(vars_with_nums, 
                            vars_with_nums %in% c("", "n", "N", "no"), 0)
  return (vars_with_nums)
}

# Extracts block from treatment label string - first letter, hyphen 2nd letter
# First letter denotes county, 2nd block within the county
extract_block <- function(scalecount) {
  scalecount_block <- str_match(scalecount$Label, "[A-Z]-[A-Z]")
  scalecount_block <- as.factor(scalecount_block)
  return (scalecount_block)
}

# Common processing across all 3 experiments
process_scalecount <- function(scalecount_raw) {
  # Delete columns with no info
  scalecount <- scalecount_raw[ -c(16:19) ]
  
  # Extract and create new column with treatment
  scalecount$Treatment<-
    substring(scalecount$Label, first=4, last=4)
  
  # Encode presence and absence as 0s and 1s
  scalecount$Prespara <- replace_strings_with_binary(scalecount$Prespara)
  scalecount$Presfungus <- replace_strings_with_binary(scalecount$Presfungus)
  
  # Extract block
  scalecount$Block <- extract_block(scalecount)
  
  # Convert Treatment to factor - need to also add County for Study 1
  scalecount$Treatment <- as.factor(scalecount$Treatment)
  
  # Convert twig to factor
  scalecount$Twigab <- as.factor(scalecount$Twigab)
  
  # Change counts to numeric
  scalecount <- scalecount %>% 
    mutate(across(c(Livescale1, Deadscale1, Livescale2, 
                    Deadscale2, Livescale3, Deadscale3), 
                  as.numeric))
  
  # Add mean live scale and Mean dead scale columns
  scalecount<-scalecount %>% 
    mutate(Meanlivescale = rowMeans(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
  scalecount<-scalecount %>% 
    mutate(Meandeadscale = rowMeans(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))
  
  # Add sum live scale and sum dead scale
  scalecount<-scalecount %>% 
    mutate(Sumlivescale = rowSums(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
  scalecount<-scalecount %>% 
    mutate(Sumdeadscale = rowSums(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))
  
  # Removing labels that were not collected
  scalecount<-scalecount %>% drop_na(Livescale1)
  
  return (scalecount)
}

# Extracts treatment means by treatment
get_treatment_survival_means <- function(scalecount) {
  tmeans <- scalecount %>% 
    group_by(Treatment) %>% 
    dplyr::summarize(
      Mean=round(mean(Meanlivescale, na.rm = T),3),
      sd = sd(Meanlivescale),
      n = n(),
      se = sd / sqrt(n),
      cv = sd/Mean * 100
    )
  return (tmeans)
}

