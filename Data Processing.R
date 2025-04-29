library(dplyr)
library(stringr)

# Add EHS and Cryptomeria scale together into Livescale{x}, Deadscale{x} cols
create_live_deadscale_vars <- function(scalecount) {
  summed <- 
    scalecount %>%
    mutate(Livescale1=EHSLivescale1 + CryptoLivescale1,
           Livescale2=scalecount$EHSLivescale2 + scalecount$CryptoLivescale2,
           Livescale3=scalecount$EHSLivescale3 + scalecount$CryptoLivescale3,
           Deadscale1=scalecount$EHSDeadscale1 + scalecount$CryptoDeadscale1,
           Deadscale2=scalecount$EHSDeadscale2 + scalecount$CryptoDeadscale2,
           Deadscale3=scalecount$EHSDeadscale3 + scalecount$CryptoDeadscale3)
  return (summed)
}

# Create live and deadscale variables out of pivot wide table
create_live_deadscale_pres_vars_wide <- function(scalecount_wide) {
  summed <- 
    scalecount_wide %>%
    mutate(Livescale1=rowSums(across(c(Livescale1_both, Livescale1_EHS, Livescale1_Crypto)), na.rm=TRUE),
           Livescale2=rowSums(across(c(Livescale2_both, Livescale2_EHS, Livescale2_Crypto)), na.rm=TRUE),
           Livescale3=rowSums(across(c(Livescale3_both, Livescale3_EHS, Livescale3_Crypto)), na.rm=TRUE),
           Deadscale1=rowSums(across(c(Deadscale1_both, Deadscale1_EHS, Deadscale1_Crypto)), na.rm=TRUE),
           Deadscale2=rowSums(across(c(Deadscale2_both, Deadscale2_EHS, Deadscale2_Crypto)), na.rm=TRUE),
           Deadscale3=rowSums(across(c(Deadscale3_both, Deadscale3_EHS, Deadscale3_Crypto)), na.rm=TRUE),
           
           # add prespara and presfungus together - will replace all >= 1 with 1 after
           Prespara=rowSums(across(c(Prespara_both, Prespara_EHS, Prespara_Crypto)), na.rm=TRUE),
           Presfungus=rowSums(across(c(Presfungus_both, Presfungus_EHS, Presfungus_Crypto)), na.rm=TRUE),
           ) %>%
    # Replace Prespara and Presfungus >=1 with 1 as it's a binary variable
    # Sometimes there's data across crypto/ehs/both observations where several had presence marked as 1
    # So rowSums adds this to 1
    # There is probably a better way to do this
    mutate(across(.cols = c(Prespara, Presfungus), .fns = function(x) ifelse(x >= 1, 1, 0)))
  
  # If all entries are NA, then rowSums will calculate a 0...
  # Check if all are NA, and change to NA
  # There's probably a better way to do this
  summed[is.na(summed$Livescale1_both) & is.na(summed$Livescale1_Crypto) 
         & is.na(summed$Livescale1_EHS), "Livescale1"] <- NA
  summed[is.na(summed$Livescale2_both) & is.na(summed$Livescale2_Crypto) 
         & is.na(summed$Livescale2_EHS), "Livescale2"] <- NA
  summed[is.na(summed$Livescale3_both) & is.na(summed$Livescale3_Crypto) 
         & is.na(summed$Livescale3_EHS), "Livescale3"] <- NA
  summed[is.na(summed$Deadscale1_both) & is.na(summed$Deadscale1_Crypto) 
         & is.na(summed$Deadscale1_EHS), "Deadscale1"] <- NA
  summed[is.na(summed$Deadscale2_both) & is.na(summed$Deadscale2_Crypto) 
         & is.na(summed$Deadscale2_EHS), "Deadscale2"] <- NA
  summed[is.na(summed$Deadscale3_both) & is.na(summed$Deadscale3_Crypto) 
         & is.na(summed$Deadscale3_EHS), "Deadscale3"] <- NA
  summed[is.na(summed$Prespara_both) & is.na(summed$Prespara_Crypto) 
         & is.na(summed$Prespara_EHS), "Prespara"] <- NA
  summed[is.na(summed$Presfungus_both) & is.na(summed$Presfungus_Crypto) 
         & is.na(summed$Presfungus_EHS), "Presfungus"] <- NA
  
  # Remove temp variables
  summed <- 
    summed %>% 
    dplyr::select(-one_of(c('Livescale1_both', 'Livescale1_EHS', 'Livescale1_Crypto',
                          'Livescale2_both', 'Livescale2_EHS', 'Livescale2_Crypto',
                          'Livescale3_both', 'Livescale3_EHS', 'Livescale3_Crypto',
                          'Deadscale1_both', 'Deadscale1_EHS', 'Deadscale1_Crypto',
                          'Deadscale2_both', 'Deadscale2_EHS', 'Deadscale2_Crypto',
                          'Deadscale3_both', 'Deadscale3_EHS', 'Deadscale3_Crypto',
                          'Prespara_both', 'Prespara_EHS', 'Prespara_Crypto',
                          'Presfungus_both', 'Presfungus_EHS', 'Presfungus_Crypto'
  ))) 
  return (summed)
}


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
  count <- scalecount %>% group_by(Label, Twigab, Date) %>%
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
    group_by(Label, Date, Treatment, Block) %>% 
    summarize(across(where(is.numeric), mean)) %>%
    # drop rows that shouldn't use means
    dplyr::select(-one_of(c('Sumlivescale_from_mean',
                            'Sumdeadscale_from_mean')))
  
  return (scalecount_tree_mean)
}

#### For Study 3 ONLY!!!
# Label denotes the same tree
# Average count data (live scale, dead scale, encarsia) 
# across both twigs to get by tree
average_counts_across_twigs_study3 <- function(scalecount) {
  scalecount_tree_mean <-
    scalecount %>%
    group_by(Label, Treatment, Block) %>% 
    summarize(across(where(is.numeric), mean)) %>%
    # drop rows that shouldn't use means
    dplyr::select(-one_of(c('Sumlivescale_from_mean',
                            'Sumdeadscale_from_mean')))
  
  return (scalecount_tree_mean)
}

# Label denotes the same tree
# Sum count data (live scale, dead scale, encarsia) 
# across both twigs to get by tree
sum_counts_across_twigs <- function(scalecount) {
  scalecount_tree_sum <-
    scalecount %>%
    group_by(Label, Date, Treatment, Block) %>% 
    summarize(across(where(is.numeric), sum)) %>%
    # drop rows that shouldn't use sums
    dplyr::select(-one_of(c('encarsia', 'Meanlivescale', 'Meandeadscale')))
  
  return (scalecount_tree_sum)
}

# Average counts across block-treatment combo
# Must have already averaged across twigs first...
average_counts_across_block_trt <- function(scalecount) {
  scalecount_block_trt_mean <-
    scalecount %>%
    group_by(Date, Treatment, Block) %>% 
    summarize(across(where(is.numeric), mean)) %>%
    arrange(Block, Treatment)
  
  return (scalecount_block_trt_mean)
}

#### For Study 3 ONLY!!!
# Average counts across block-treatment combo
# Must have already averaged across twigs first...
average_counts_across_block_trt_study3 <- function(scalecount) {
  scalecount_block_trt_mean <-
    scalecount %>%
    group_by(Treatment, Block) %>% 
    summarize(across(where(is.numeric), mean)) %>%
    arrange(Block, Treatment)
  
  return (scalecount_block_trt_mean)
}

# replace binary strings with 0s and 1s
replace_strings_with_binary <- function(vars) {
  vars_with_nums <- replace(vars, vars %in% c("yes", "y", "Y"), 1)
  
  # This is the criteria for both study 1 and study 3
  # study 2 doesn't include "", but there's no data with a blank anyway
  # that isn't also missing Livescale1 (gets filtered out)
  vars_with_nums <- replace(vars_with_nums, 
                            vars_with_nums %in% c("", "n", "N", "no"), 0)
  return (vars_with_nums)
}

# Extracts block from treatment label string - first letter, hyphen 2nd letter
# First letter denotes county, 2nd block within the county
extract_block <- function(scalecount) {
  # Check first element - if doesn't include county, only take 1st char
  if (str_length(scalecount[1,]$Label) == 3) {
    scalecount_block <- substring(scalecount$Label, 1, 1)
  } else {
    scalecount_block <- str_match(scalecount$Label, "[A-Z]-[A-Z]")
  }
  
  scalecount_block <- as.factor(scalecount_block)
  return (scalecount_block)
}

# mean live and dead scale are not integers for count data
# Because we have missing shoots
# need to multiply all means by 3 instead of summing for accuracy...
get_sum_scale_from_mean <- function(mean_scale) {
  sum_scale_from_mean <- round(mean_scale *3)
  return (sum_scale_from_mean)
}

extract_treatment <- function(scalecount) {
  # if label doesn't include county
  if (str_length(scalecount[1,]$Label) == 3) {
    treatment <- substring(scalecount$Label, first=2, last=2)
  } else {
    treatment <- substring(scalecount$Label, first=4, last=4)
  }
  
  treatment <- as.factor(treatment)
  return (treatment)
}

# Common processing across all experiments
process_scalecount <- function(scalecount_raw, in_group=FALSE) {
  scalecount <- scalecount_raw
  
  # Extract and create new column with treatment
  # Use 2nd to last character
  scalecount$Treatment<- extract_treatment(scalecount)
  
  # Encode presence and absence as 0s and 1s
  scalecount$Prespara <- replace_strings_with_binary(scalecount$Prespara)
  scalecount$Presfungus <- replace_strings_with_binary(scalecount$Presfungus)
  scalecount$Prespara <- as.numeric(scalecount$Prespara)
  scalecount$Presfungus <- as.numeric(scalecount$Presfungus)
  
  # Extract block
  scalecount$Block <- extract_block(scalecount)
  
  # Convert twig to factor
  scalecount$Twigab <- as.factor(scalecount$Twigab)
  
  # Convert label to factor
  scalecount$Label <- as.factor(scalecount$Label)
  
  # Change counts to numeric
  scalecount <- scalecount %>% 
    mutate(across(c(Livescale1, Deadscale1, Livescale2, 
                    Deadscale2, Livescale3, Deadscale3), 
                  as.numeric))
  
  # Add mean live scale and Mean dead scale columns
  # This will technically retain labels with 2 shoots?
  scalecount<-scalecount %>%
    mutate(Meanlivescale = rowMeans(dplyr::select(., Livescale1, Livescale2, Livescale3), na.rm = TRUE))
  scalecount<-scalecount %>%
    mutate(Meandeadscale = rowMeans(dplyr::select(., Deadscale1, Deadscale2, Deadscale3), na.rm = TRUE))

  # Get the scaled "sum: of the mean live scale for use with count models
  scalecount$Sumlivescale_from_mean <- 
    get_sum_scale_from_mean(scalecount$Meanlivescale)
  scalecount$Sumdeadscale_from_mean <- 
    get_sum_scale_from_mean(scalecount$Meandeadscale)
  
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
      Sum_from_mean=mean(Meanlivescale * 3),
      sd = sd(Meanlivescale),
      n = n(),
      se = sd / sqrt(n),
      cv = sd/Mean * 100
    )
  return (tmeans)
}

# Extracts treatment means by county
get_treatment_survival_means_county <- function(scalecount) {
  tmeans <- scalecount %>% 
    group_by(Treatment, County) %>% 
    dplyr::summarize(
      Mean=round(mean(Meanlivescale, na.rm = T),3),
      Sum_from_mean=mean(Meanlivescale * 3),
      sd = sd(Meanlivescale),
      n = n(),
      se = sd / sqrt(n),
      cv = sd/Mean * 100
    )
  return (tmeans)
}

