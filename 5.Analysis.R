#================================
# Analysis
#================================

# Load packages
library(dplyr)
library(broom) # For models 
library(broom.mixed) # For mixed models 
library(parallelly)
library(ggplot2) # Plotting
library(data.table)
library(lme4)
library(stringr)
library(effects)
library(lmerTest)
library(scales)
library(ggpubr)
library(standardize)
library(egg)
library(sjPlot)
library(forcats)
library(vegan)
library(sjmisc)

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file again
All_regions <- read.csv('final_abundance_changes.csv', sep = ',') 

# Correct some naming formats
All_regions <- All_regions %>%
  mutate(Habitat = case_when(Habitat == "benthopelagic" ~ "Benthopelagic", Habitat == "demersal" ~ "Demersal", Habitat == "reef-associated" ~ "Reef-associated", Habitat == "pelagic" ~ "Pelagic", TRUE ~ Habitat)) %>%
  mutate(PlaceofDevelopment = case_when(PlaceofDevelopment == 'planktonic' ~ 'Pelagic', TRUE ~ PlaceofDevelopment))

# make habitat a factor
All_regions$Habitat <- as.factor(All_regions$Habitat)

# make place of larval development a factor
All_regions$PlaceofDevelopment <- as.factor(All_regions$PlaceofDevelopment)

#-------------------------------------
# 2. Define regions
#-------------------------------------

# Create region subsets
Oceania <- All_regions %>% filter(region == 'Australia' | region == 'New Zealand' | region == 'Indonesia')
Eastern_US <- All_regions %>% filter(region == 'Eastern_US' | region == 'Belize' | region == 'Panama' | region == 'Canada' | region == 'Caribbean')
Europe <- All_regions %>% filter(region == 'Europe' | region == 'Spain')
Western_US <- All_regions %>% filter(region == 'Western_US')
Other <- All_regions %>% filter(region == 'Antartica' | region == 'Pacific' | region == 'South_America') 

#-------------------------------------
# 3. Run basic model for various regions and durations
#-------------------------------------

# Create model function
run_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = data)
  print(summary(model))
  freq_plot <- ggplot(data, aes(position_in_range,fill = fct_rev(factor(abundance_temp_binary)), alpha=0.8)) +
    geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) +
    scale_fill_manual(values=c("orange", "steelblue2")) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    ggtitle(deparse(substitute(data))) + theme_classic() + theme(legend.position="none")
  return(freq_plot)
 }

#-------------------------------------
# 3.1 Time series of all durations
#-------------------------------------

# All regions
run_model(All_regions)

# Run model for each region
run_model(Eastern_US)
run_model(Europe)
run_model(Oceania)
run_model(Western_US)
run_model(Other)

#-------------------------------------
# 3.2 Time series of 20+ year duration
#-------------------------------------

# Filter for 20+ year duration
All_twenty <- All_regions %>% filter(duration > 19)

# Create region subsets
Oceania_twenty <- Oceania %>% filter(duration > 19)
Eastern_US_twenty <- Eastern_US %>% filter(duration > 19)
Europe_twenty <- Europe %>% filter(duration > 19)
Western_US_twenty <- Western_US %>% filter(duration > 19)

# All regions
run_model(All_twenty)

# Run model for each region
run_model(Oceania_twenty)
run_model(Eastern_US_twenty)
run_model(Europe_twenty)
run_model(Western_US_twenty)

#-------------------------------------
# 3.3 Time series of 30+ year duration
#-------------------------------------

# Filter for 30+ year duration
All_thirty <- All_regions %>% filter(duration > 29)

# Create region subsets
Oceania_thirty <- Oceania %>% filter(duration > 29)
Eastern_US_thirty <- Eastern_US %>% filter(duration > 29)
Europe_thirty <- Europe %>% filter(duration > 29)
Western_US_thirty <- Western_US %>% filter(duration > 29)

# All regions
run_model(All_thirty)

# Run model for each region
run_model(Oceania_thirty)
run_model(Eastern_US_thirty)
run_model(Europe_thirty)
run_model(Western_US_thirty)

#-------------------------------------
# 3.4 Time series of 40+ year duration
#-------------------------------------

# Filter for 40+ year duration
All_forty <- All_regions %>% filter(duration > 39)

# Create region subsets
Oceania_forty <- Oceania %>% filter(duration > 39)
Eastern_US_forty <- Eastern_US %>% filter(duration > 39)
Europe_forty <- Europe %>% filter(duration > 39)
Western_US_forty <- Western_US %>% filter(duration > 39)

# All regions
run_model(All_forty)

# Run model for each region
run_model(Oceania_forty)
run_model(Eastern_US_forty)
run_model(Europe_forty)
run_model(Western_US_forty)

#-------------------------------------
# 3.5 Survey as a random factor
#-------------------------------------

# Create survey model function
run_survey_model <- function(data) {
  model <- glmer(abundance_mean_temp_binary ~ position_in_range + (1|accepted_name) + (1|survey), family = binomial("logit"), data = data)
  print(summary(model))
  freq_plot <- ggplot(data, aes(position_in_range,fill = fct_rev(factor(abundance_mean_temp_binary)), alpha=0.8)) +
    geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) +
    scale_fill_manual(values=c("orange", "steelblue2")) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    ggtitle(deparse(substitute(data))) + theme_classic() + theme(legend.position="none")
  plot(freq_plot)
}

# All regions
run_survey_model(All_regions)

# Run model for each region
run_survey_model(Oceania)
run_survey_model(Eastern_US)
run_survey_model(Europe)
run_survey_model(Western_US)

#-------------------------------------
# 4. Ecological trait effects
#-------------------------------------

#-------------------------------------
# 4.1 Ocean habitat
#-------------------------------------

# Create habitat model function
run_habitat_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:Habitat + (1|accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  habitat_plot_avg <- plot_model(model, type = "int", terms = c("position_in_range [all]", "habitat"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate", "purple")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate", "purple")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))  
  print(habitat_plot_avg)
}

# All regions
run_habitat_model(All_regions)  

# Run model for each region
run_habitat_model(Eastern_US)
run_habitat_model(Europe)
run_habitat_model(Oceania)
run_habitat_model(Western_US)

#-------------------------------------
# 4.2 Larval phase
#-------------------------------------

# Create habitat model function
run_larval_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:PlaceofDevelopment + (1|accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  larval_plot_avg <- plot_model(model, type = "int", terms = c("position_in_range [all]", "larval"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate", "purple")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate", "purple")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))  
  print(larval_plot_avg)
}

# All regions
run_larval_model(All_regions)

# Run model for each region
run_larval_model(Eastern_US)
run_larval_model(Europe)
run_larval_model(Oceania)
run_larval_model(Western_US)

#-------------------------------------
# 4.3 Trophic level
#-------------------------------------

# Create trophic level model function
run_trophic_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:Troph + (1 | accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  plot_avg <- plot_model(model, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "trophic level"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate"), labels = c("Low", "Middle", "High")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) 
  print(plot_avg)
}

# All regions
run_trophic_model(All_regions)

# Run model for each region
run_trophic_model(Eastern_US)
run_trophic_model(Europe)
run_trophic_model(Oceania)
run_trophic_model(Western_US)

#-------------------------------------
# 4.4 Depth
#-------------------------------------

# Create depth model function
run_depth_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logDepth + (1 | accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  depth_plot_avg <- plot_model(model, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "logdepth"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate"), labels = c("Low depth", "Medium depth", "Large depth")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))  
  print(depth_plot_avg)
}

# All regions
run_depth_model(All_regions)

# Run model for each region
run_depth_model(Eastern_US)
run_depth_model(Europe)
run_depth_model(Oceania)
run_depth_model(Western_US)

#-------------------------------------
# 5. Life history effects
#-------------------------------------

#-------------------------------------
# 5.1 Growth rate
#-------------------------------------

# Create growth rate model function
run_k_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logK + (1 | accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  plot_avg <- plot_model(model, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "generation time"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate"), labels = c("Long", "Medium", "Short")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))  
  print(plot_avg)
}

# All regions
run_k_model(All_regions)

#Run model for each region
run_k_model(Eastern_US)
run_k_model(Europe)
run_k_model(Oceania)
run_k_model(Western_US)

#-------------------------------------
# 5.2 Fecundity
#-------------------------------------

# Create fecundity model function
run_fecundity_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logFecundity + (1 | accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  fecundity_plot_avg <- plot_model(model, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "logdepth"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate"), labels = c("Low fecundity", "Medium fecundity", "High fecundity")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01))   
  print(fecundity_plot_avg)
}

# All regions
run_fecundity_model(All_regions)

# Run model for each region
run_fecundity_model(Eastern_US)
run_fecundity_model(Europe)
run_fecundity_model(Oceania)
run_fecundity_model(Western_US)

#-------------------------------------
# 5.3 Body length
#-------------------------------------

# Create length model function
run_length_model <- function(data) {
  title <- deparse(substitute(data))
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logLength + (1 | accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  length_plot_avg <- plot_model(model, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "logLength"), title = title) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate")) +
    theme_classic() + scale_color_manual(values = c("#0072B2", "black", "chocolate"), labels = c("Small body length", "Medium body length", "Large body length")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) 
  print(length_plot_avg) 
} 

# All regions
run_length_model(All_regions)

# Run model for each region
run_length_model(Eastern_US)
run_length_model(Europe)
run_length_model(Oceania)
run_length_model(Western_US)

#-------------------------------------
# 6. Duration effect model
#-------------------------------------

# Create duration model function
run_duration_model <- function(data) {
  model <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:duration + (1 | accepted_name_numeric), family = binomial("logit"), data = data)
  print(summary(model))
  plot_avg <- plot_model(model, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "duration"), title = deparse(substitute(data))) + ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "black", "chocolate")) +
    theme_classic() + labs(color = "Duration of Timeseries", fill = "Duration of Timeseries") + scale_color_manual(values = c("#0072B2", "black", "chocolate"), labels = c("Short", "Medium", "Long")) +
    theme(legend.justification = c("right", "bottom"), legend.position = c(.95, .05), legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) 
  return(plot_avg)
}

# All regions
run_duration_model(All_regions)

# Run model for each region
run_duration_model(Oceania)
run_duration_model(Eastern_US)
run_duration_model(Europe)
run_duration_model(Western_US)

# end of script 