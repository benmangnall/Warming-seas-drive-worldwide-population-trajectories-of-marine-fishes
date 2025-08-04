#================================
# Create figures
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
library(viridis)
library(vegan)
library(ggpubr)
library(grid)         
library(maps)   

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file again
All_regions <- read.csv('final_abundance_changes.csv', sep = ',') 

# Correct some naming formats
All_regions <- All_regions %>%
  mutate(Habitat = case_when(Habitat == "benthopelagic" ~ "Benthopelagic", Habitat == "demersal" ~ "Demersal", Habitat == "reef-associated" ~ "Reef-associated", Habitat == "pelagic" ~ "Pelagic", TRUE ~ Habitat)) %>%
  mutate(PlaceofDevelopment = case_when(PlaceofDevelopment == 'planktonic' ~ 'Pelagic', TRUE ~ PlaceofDevelopment))

# Make habitat a factor
All_regions$Habitat <- as.factor(All_regions$Habitat)

# Make place of larval development a factor
All_regions$PlaceofDevelopment <- as.factor(All_regions$PlaceofDevelopment)

#-------------------------------------
# 2. Define regions
#-------------------------------------

# Create region subsets
Oceania <- All_regions %>% filter(region == 'Australia' | region == 'New Zealand' | region == 'Indonesia')
Eastern_US <- All_regions %>% filter(region == 'Eastern_US' | region == 'Belize' | region == 'Panama' | region == 'Canada' | region == 'Caribbean')
Europe <- All_regions %>% filter(region == 'Europe' | region == 'Spain')
Western_US <- All_regions %>% filter(region == "Western_US")
Other <- All_regions %>% filter(region == 'Antartica' | region == 'Pacific' | region == 'South_America') 

#-------------------------------------
# 3. Figure 1
#-------------------------------------

#-------------------------------------
# 3.1 Fig. 1A - Time series duration on world map
#-------------------------------------

# Get the map outline
world <- map_data("world")

# Create hex map
hex_map <- ggplot(All_regions, aes(x = longitude_gridcell, y = latitude_gridcell)) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill="grey", colour="grey70", linewidth=0.2, alpha=0.3) +
  geom_hex(bins=100) +
  scale_fill_gradientn(
    colors = c("#D9E1F8", "#8094CC", "#DA3A3A"),
    trans = "log10",
    name = "Number of time series",
    guide = guide_colorbar(
      barheight = unit(5, "mm"),   # Increase bar height for visibility
      barwidth = unit(100, "mm"),  # Keep the bar wide for horizontal appearance
      direction = "horizontal",    # Ensures the legend remains horizontal
      title.theme = element_text(margin = margin(r = 17, b = 19))
    )
  ) +
  theme_void() +
  theme(
    legend.position = c(0.39, 0.2),  # 0.39, 0.2 Adjust these values to move the legend
    legend.justification = c(0.0, 0.0),  # Aligns the legend within the plot
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  )

# Display hex map
hex_map

# Export hex map
ggsave(hex_map, file="hex_map.pdf",width = 40, height = 20, unit = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 3.2 Fig. 1B - Taxonomic order frequency plot
#-------------------------------------

# Filter to include only orders that occur in more than 2000 time series
taxa <- All_regions %>%
  group_by(Order) %>%
  summarise(n = n()) %>%
  filter(n > 2000)

# Create order plot
orderbarplot <- ggplot(data = taxa, aes(x = reorder(Order, -n), y=n))+ 
  geom_bar(stat = "identity", width=0.5,fill='red4') +
  coord_flip() + theme_classic() +
  xlab("Taxonomic Order")+ylab("Number of time series")+
  scale_y_continuous(limits = c(0,42000), expand = c(0,0)) +
  theme(axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 15, color = '#2C2323'),
        axis.title.x = element_text(margin = margin(t = 10)))

# Display order plot
orderbarplot

# Export order plot
ggsave(orderbarplot, file="orderbarplot.pdf",width = 40, height = 20, unit = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 3.3 Fig. 1C - Temperature change over time
#-------------------------------------

# Load in temperature and year values for grid cells 
temp <- fread('abundances.csv') %>%
  select(latitude_gridcell, longitude_gridcell, mean_temperature, year) %>%
  mutate(grid_cell = paste(latitude_gridcell, longitude_gridcell, sep = '_')) %>%
  select(-latitude_gridcell, -longitude_gridcell) 

# Get the post-filtered grid cells 
correct_grid_cells <- fread('final_abundance_changes.csv') %>%
  select(latitude_gridcell, longitude_gridcell) %>%
  mutate(grid_cell = paste(latitude_gridcell, longitude_gridcell, sep = '_')) %>%
  select(-latitude_gridcell, -longitude_gridcell) %>%
  unique()

# Join the temperature data to the grid cell data
temp_data <- correct_grid_cells %>%
  left_join(temp, by = 'grid_cell') %>%
  unique()

# Obtain the estimate of temperature change over time
temp_change_model <- lmer(mean_temperature ~ year + (1|grid_cell),data=temp_data)
summary(temp_change_model)

# Convert estimate to an effect for plotting 
effect_temp_change <- as.data.frame(effect("year", temp_change_model))

# Define theme for plot
temp_plot_theme <- theme(
  axis.title.x = element_text(size = 40, vjust=1.5, hjust=0.5, colour='black'),
  axis.title.y = element_text(size = 35, colour='black'),
  axis.text = element_text(size = 32, colour='black'),
  axis.text.x = element_text(size = 32, vjust = 0.5, colour='black'),
  plot.title = element_text(lineheight=.8, face="bold", size = 10, hjust=0.5, vjust=0.2),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 13),
  legend.position.inside = c(0.9,0.1),
  legend.justification = c("right","bottom"))

# Create plot
temp_plot <- ggplot(data=effect_temp_change, aes(x=year,y=fit)) + 
  geom_line(colour="#C7340C") +
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=0, alpha=0.3, colour="#EBC220", fill="#EBC220") +
  scale_y_continuous(name="Mean SST (°C)", limits = c(9, 12)) +
  scale_x_continuous(name="Year") +
  theme_bw() + temp_plot_theme

# Display plot
temp_plot

# Export plot
ggsave(temp_plot, file="temp_plot.pdf",width = 30, height = 30, unit = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 4. Fig. 2 
#-------------------------------------

#-------------------------------------
# 4.1 Defining plot themes
#-------------------------------------

# Define theme for plots
effect_theme_fig2 <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 21, vjust=0.4),
  axis.title.y = element_text(size = 21),
  axis.text = element_text(size = 20, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 19, hjust=0.5, vjust=0.000001),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c(1, 0), legend.position = c(0.9,0.1),legend.text = element_text(size = 15))

# Define theme for duration plot
effect_theme_duration <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 21, vjust=0.4),
  axis.title.y = element_text(size = 19),
  axis.text = element_text(size = 20, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 19, hjust=0.5, vjust=0.000001),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c(1, 0), legend.position = c(0.9,0.1),legend.text = element_text(size = 15))

# Define theme to strip away x-axis labels
strip_x <- theme(
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank())

#-------------------------------------
# 4.2 Run duration model for plot
#-------------------------------------

# Duration model
duration <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:duration + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(duration)

#-------------------------------------
# 4.3 Fig 2.A, B, C, D, E, F - Time series density plots
#-------------------------------------

# Create duration plot
duration_plot <- plot_model(duration, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "duration"), title = '') +
  ylab("Positive trend probability") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Short", "Medium", "Long")) +
  labs(color = "Timeseries Duration", fill  = "Timeseries Duration") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  effect_theme_duration +
  theme(legend.position = c(0.1, 0.96), legend.justification = c(0, 1), legend.title = element_text(size = 19), legend.text = element_text(size = 20))
duration_plot

# Create plot for all regions 
All_regions_plot <- ggplot(All_regions, aes(position_in_range,fill = fct_rev(factor(abundance_temp_binary)), alpha=0.8)) +
  geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) + ylab('Density of time series') +
  scale_fill_manual(values=c("orange", "steelblue2")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  ggtitle('') + effect_theme_fig2 + theme(legend.position="none") 
All_regions_plot

# Create plot for the Northwest Atlantic 
Eastern_US_plot <- ggplot(Eastern_US, aes(position_in_range,fill = fct_rev(factor(abundance_temp_binary)), alpha=0.8)) +
  geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) + ylab('Density of time series') +
  scale_fill_manual(values=c("orange", "steelblue2")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  ggtitle('') + effect_theme_fig2 + theme(legend.position="none") 
Eastern_US_plot

# Create plot for Europe
Europe_plot <- ggplot(Europe, aes(position_in_range,fill = fct_rev(factor(abundance_temp_binary)), alpha=0.8)) +
  geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) + ylab('Density of time series') +
  scale_fill_manual(values=c("orange", "steelblue2")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  ggtitle('') + effect_theme_fig2 + theme(legend.position="none") 
Europe_plot

# Create plot for Oceania
Oceania_plot <- ggplot(Oceania, aes(position_in_range,fill = fct_rev(factor(abundance_temp_binary)), alpha=0.8)) +
  geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) + ylab('Density of time series') +
  scale_fill_manual(values=c("orange", "steelblue2")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  ggtitle('') + effect_theme_fig2 + theme(legend.position="none") 
Oceania_plot

# Create plot for the Northeast Pacific 
Western_US_plot <- ggplot(Western_US, aes(position_in_range,fill = fct_rev(factor(abundance_temp_binary)), alpha=0.8)) +
  geom_density(adjust = 1/5,linewidth =0.15, alpha=.5) + ylab('Density of time series') +
  scale_fill_manual(values=c("orange", "steelblue2")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  ggtitle('') + effect_theme_fig2 + theme(legend.position="none") 
Western_US_plot

# Strip x-axis labels
All_regions_plot   <- All_regions_plot   + strip_x
duration_plot      <- duration_plot      + strip_x
Europe_plot        <- Europe_plot        + strip_x
Oceania_plot       <- Oceania_plot       + strip_x
Eastern_US_plot    <- Eastern_US_plot    + strip_x
Western_US_plot    <- Western_US_plot    + strip_x

# Arrange plots
fig2 <- ggarrange(All_regions_plot, duration_plot, Europe_plot, Oceania_plot, Eastern_US_plot, Western_US_plot, nrow = 3, ncol = 2)

# Apply one x-axis label
fig2 <- annotate_figure(fig2, bottom = text_grob("Position in Range (0 = equatorward limit, 1 = poleward limit)", size = 21, vjust = 0.3))

# Export figure
ggsave(fig2, file="fig2.pdf", width = 50, height = 25, units = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 5. Fig. 3 
#-------------------------------------

#-------------------------------------
# 5.1 Define theme for plots
#-------------------------------------

# Define theme for plots
effect_theme_fig3 <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 23, vjust=0.1),
  axis.title.y = element_text(size = 22),
  axis.text = element_text(size = 25, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 19, hjust=0.5, vjust=0.000001),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c(1, 0), legend.position = c(0.9,0.1),legend.text = element_text(size = 15))

#-------------------------------------
# 5.2 Run models needed for plots
#-------------------------------------

# Habitat model
habitat <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:Habitat + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(habitat)

# Larval habitat model
larval_habitat <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:PlaceofDevelopment + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(larval_habitat)

# Trophic level model
trophic <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:Troph + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(trophic)

# Depth model
depth <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logDepth + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(depth)

#-------------------------------------
# 5.3 Fig 3. A, B, C, D - Ecological trait effect size plots
#-------------------------------------

# Create plot for different habitats
habitat_plot <- plot_model(habitat, type = "int", terms = c("position_in_range [all]", "Habitat"), title = '') + ylab("") + xlab("") + scale_fill_manual(name = "Habitat type", values = c("#0072B2", "darkgrey", "chocolate", "purple")) +
  theme_classic() + scale_color_manual(name = "Habitat type", values = c("#0072B2", "darkgrey", "chocolate", "purple")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) + effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.59, .65), legend.title = element_text(size = 21), legend.text = element_text(size = 19))
print(habitat_plot)

# Create plot for different larval phase habitats
larval_plot <- plot_model(larval_habitat, type = "int", terms = c("position_in_range [all]", "PlaceofDevelopment"), title = '') + ylab("") + xlab("") + scale_fill_manual(name = "Habitat type", values = c("#0072B2", "darkgrey", "chocolate", "purple")) +
  theme_classic() + scale_color_manual(name = "Larval habitat type", values = c("#0072B2", "darkgrey", "chocolate", "purple")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) + effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.96, .08), legend.title = element_text(size = 23), legend.text = element_text(size = 25)) 
print(larval_plot)

# Create plot for different trophic levels
trophic_plot <- plot_model(trophic, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "trophic"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Low", "Medium", "High")) +
  labs(color = "Trophic level", fill  = "Trophic level") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.9, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
trophic_plot

# Create plot for different depths
depth_plot <- plot_model(depth, mdrt.values = "quart", type = "pred", terms = c("position_in_range [all]", "logDepth"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Low", "Medium", "High")) +
  labs(color = "Depth", fill  = "Depth") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.9, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
depth_plot

# Arrange plots 
fig3 <- ggarrange(habitat_plot, larval_plot, trophic_plot, depth_plot, nrow=2, ncol = 2)

# Add x and y axis labels
fig3 <- annotate_figure(fig3,
  bottom = grid::textGrob("Position in Range (0 = equatorward limit, 1 = poleward limit)", gp = gpar(fontsize = 25), vjust = 0.1),
  left = grid::textGrob("Abundance change with warming (Probability of positive trend)", rot = 90, gp = gpar(fontsize = 25), vjust = 0.9))

# Export figure
ggsave(fig3, file="fig3.pdf", width = 35, height = 30, units = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 6. Fig. 4
#-------------------------------------

#-------------------------------------
# 6.1 Define theme for plots
#-------------------------------------

# Define theme for plots
effect_theme_fig4 <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 15, vjust=0.1),
  axis.title.y = element_text(size = 15),
  axis.text = element_text(size = 22, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 19, hjust=0.5, vjust=0.000001),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c(1, 0), legend.position = c(0.9,0.1),legend.text = element_text(size = 15))

# Define theme to strip x and y axis labels
strip_ <- theme(
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank())

#-------------------------------------
# 6.2 Run models needed for plots
#-------------------------------------

# Growth_rate model
growth_rate <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logK + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(growth_rate)

# Fecundity model
fecundity <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logFecundity + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(fecundity)

# Body_length model
body_length <- glmer(abundance_temp_binary ~ position_in_range + position_in_range:logLength + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(body_length)

#-------------------------------------
# 6.3 Fig 3. A, B, C - Life history trait effect size plots
#-------------------------------------

# Create plot for different growth rates
growth_rate_plot <- plot_model(growth_rate, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "growth rate"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Slow", "Medium", "Fast")) +
  labs(color = "Growth Rates", fill  = "Growth Rates") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + effect_theme_fig4 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.9, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
growth_rate_plot

# Create plot for different fecundities
fecundity_plot <- plot_model(fecundity, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "Fecundity"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Low", "Medium", "High")) +
  labs(color = "Fecundity", fill  = "Fecundity") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + effect_theme_fig4 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.84, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
fecundity_plot

# Create plot for different body lengths
body_length_plot <- plot_model(body_length, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "Body length"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Short", "Medium", "Long")) +
  labs(color = "Body length", fill  = "Body length") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + effect_theme_fig4 + 
  theme(legend.justification = c("right", "bottom"), legend.position = c(.865, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
body_length_plot

# Strip away x and y axis labels
growth_rate_plot <- growth_rate_plot + strip_
fecundity_plot <- fecundity_plot + strip_
body_length_plot <- body_length_plot + strip_

# Arrange plots
fig4 <- ggarrange(growth_rate_plot, fecundity_plot, body_length_plot, nrow=3, ncol = 1)

# Annotate figure with x and y axis labels 
fig4 <- annotate_figure(fig4,
  bottom = grid::textGrob("Position in Range (0 = equatorward limit, 1 = poleward limit)", gp = gpar(fontsize = 24), vjust = 0.45),
  left = grid::textGrob("Abundance change with warming (Probability of positive trend)", rot = 90, gp = gpar(fontsize = 28), vjust = 0.5))

# Export figure
ggsave(fig4, file="fig4.pdf", width = 25, height = 37.5, units = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 7. Fig. 5 
#-------------------------------------

#-------------------------------------
# 7.1 Defining theme for plots
#-------------------------------------

# Define theme for plots
effect_theme_fig5 <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 41, vjust=0.3),
  axis.title.y = element_text(size = 43),
  axis.text = element_text(size = 38, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 19, hjust=0.5, vjust=0.000001),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c(1, 0), legend.position = c(0.9,0.1),legend.text = element_text(size = 15))

#-------------------------------------
# 7.2 Fig. 5F - Community response magnitude on world map
#-------------------------------------

# Get the community responses
community <- All_regions %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  reframe(mean_response = glm(abundance_temp_binary ~ position_in_range) %>% tidy() %>% filter(term == "position_in_range") %>% pull(estimate)) %>%
  filter(mean_response > -1 & mean_response < 1)

# Get the world outline 
world <- map_data("world")

# Plot the community responses onto the world map
hex_map_mean <- ggplot(community, aes(x = longitude_gridcell, y = latitude_gridcell)) +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill   = "grey",
               colour = "grey70",
               linewidth = 0.2,
               alpha  = 0.3) +
  stat_summary_hex(aes(z = mean_response), fun  = mean, bins = 100) +
  scale_fill_gradientn(
    colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"),
    limits = c(-1, 1),                  
    oob    = scales::squish,            
    breaks = c(-1, 0, 1),               
    labels = c("-1", "0", "1"),         
    name   = "Response\nMagnitude",
    guide  = guide_colorbar(
      barheight   = unit(7,  "mm"),
      barwidth    = unit(100, "mm"),
      direction   = "horizontal",
      title.theme = element_text(margin = margin(r = 13, b = 23))
    )) +
  theme_void() +
  theme(
    legend.position      = c(0.37, 0.2),
    legend.justification = c(0.0, 0.0),
    legend.text          = element_text(size = 15),
    legend.title         = element_text(size = 16),
    plot.title           = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle        = element_text(hjust = 0.5, size = 14))

# Display hex map
print(hex_map_mean)

# Export hex map
ggsave(hex_map_mean, file="hex_map_mean.pdf", width = 40, height = 20, units = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 7.3 Fig. 5A, B, C, D, E - Community response histograms by region
#-------------------------------------

#-------------------------------------
# 7.3.1 Obtaining the community response for each region
#-------------------------------------

# Europe 
Europe_community <- Europe %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  reframe(mean_response = glm(abundance_temp_binary ~ position_in_range) %>% tidy() %>% filter(term == "position_in_range") %>% pull(estimate)) %>%
  filter(mean_response > -1 & mean_response < 1)

# Oceania
Oceania_community <- Oceania %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  reframe(mean_response = glm(abundance_temp_binary ~ position_in_range) %>% tidy() %>% filter(term == "position_in_range") %>% pull(estimate)) %>%
  filter(mean_response > -1 & mean_response < 1)

# Northwest Atlantic
Eastern_US_community <- Eastern_US %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  reframe(mean_response = glm(abundance_temp_binary ~ position_in_range) %>% tidy() %>% filter(term == "position_in_range") %>% pull(estimate)) %>%
  filter(mean_response > -1 & mean_response < 1)

# Northwest Pacific
Western_US_community <- Western_US %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  reframe(mean_response = glm(abundance_temp_binary ~ position_in_range) %>% tidy() %>% filter(term == "position_in_range") %>% pull(estimate)) %>%
  filter(mean_response > -1 & mean_response < 1)

#-------------------------------------
# 7.3.2 Plotting the community response for each region
#-------------------------------------

# All regions
All_community <- ggplot(community, aes(x = mean_response)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
    limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  labs(x = "Response magnitude", y = "Number of communities") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60), axis.text = element_text(size = 45))

# Europe
Europe_community <- ggplot(Europe_community, aes(x = mean_response)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))

# Oceania
Oceania_community <- ggplot(Oceania_community, aes(x = mean_response)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))

# Northwest Atlantic
Eastern_US_community <- ggplot(Eastern_US_community, aes(x = mean_response)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))

# Northwest Pacific
Western_US_community <- ggplot(Western_US_community, aes(x = mean_response)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))

# Export plots
ggsave(All_community, file="All_community.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Europe_community, file="Europe_community.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Oceania_community, file="Oceania_community.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Eastern_US_community, file="Eastern_US_community.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Western_US_community, file="Western_US_community.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 7.4 Fig. 5G, H - Example calculation plots
#-------------------------------------

# Get the positive example grid cell
community_example_pos <- All_regions %>% filter(latitude_gridcell == '59' & longitude_gridcell == '-169')

# Plot the data from the positive example grid cell
community_example_pos <- ggplot(community_example_pos, aes(x = position_in_range, y = abundance_mean_temp_correlation)) +
  geom_point(size = 5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
  labs(
    x     = "Position in Latitudinal Range",
    y     = "Abundance response to temperature") +
  effect_theme_fig5

# Get the negative example grid cell
community_example_neg <- All_regions %>% filter(latitude_gridcell == '48' & longitude_gridcell == '-8')

# Plot the data from the negative example grid cell
community_example_neg <- ggplot(community_example_neg, aes(x = position_in_range, y = abundance_mean_temp_correlation)) +
  geom_point(size = 5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
  labs(
    x     = "Position in Latitudinal Range",
    y     = "Abundance response to temperature") +
  effect_theme_fig5

# Export plots 
ggsave(community_example_pos, file="community_example_pos.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(community_example_neg, file="community_example_neg.pdf", width = 30, height = 30, units = "cm", dpi = 600, device = cairo_pdf)

# end of script