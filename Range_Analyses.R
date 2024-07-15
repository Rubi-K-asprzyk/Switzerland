#!/usr/bin/env Rscript

# --------------------------- #
# Swiss Patterns / K.Thibault #
# --------------------------- #

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT SWISSPATTERNS.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 
cat(rule(left = "INITIALISATION", line_col = "green", line = "-", col = "br_green"))

#------------------------#
##### INITIALISATION #####
#------------------------#

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})

# Install/load tons of packages.
p_load(doParallel, # Allow parallel computation
       ape,        # Multiple tools of phylogenetic analyses
       dplyr,      # Allow the use of the plyr syntax
       tidyr,      # Allow the use of function "pivot_longer"
       ggplot2,
       purrr,
       ggpubr,
       ggnewscale)

# Set the parallel backend
registerDoParallel(cores=2)

# Set a global option to disable dplyr::summarise() warnings when grouping. 
options(dplyr.summarise.inform = TRUE)

# Call the functions needed to transform the PW matrices into vector
source(paste0(getwd(),"/Utilities/PW_to_Vector.R"))

# Message
cat(rule(left = "- Packages loaded - ", line_col = "white", line = " ", col = "green"))

# --- Theme settings --- #

# --- #
cat(rule(left = "PLOT THEME SETTINGS", line_col = "green", line = "-", col = "br_green"))
# --- #

# Create a color-scheme based on the taxa to use for all the plots. 
color_scheme <-  c("Liverworts" = "yellowgreen", "Mosses" = "darkgreen")

# Choose a palette (of the rColorBrewer package, max 2 breaks with the "Paired" palette. 
Pal_col <- "Paired"

# This theme extends the 'theme_light' that comes with ggplot2.
# The "Lato" font is used as the base font. This is similar
# to the original font in Cedric's work, Avenir Next Condensed.
# theme_set(theme_light(base_family = "Lato"))

# Set one known theme
theme_set(theme_light())
# Modify it
theme_update(
  # Remove title for both x and y axes
  # axis.title = element_blank(),
  # Axes labels are grey
  axis.text = element_text(color = "grey40"),
  # The size of the axes labels are different for x and y.
  axis.text.x = element_text(size = 10, margin = margin(t = 5)),
  axis.text.y = element_text(size = 10, margin = margin(r = 5)),
  # Also, the ticks have a very light grey color
  axis.ticks = element_line(color = "grey91", linewidth = .5),
  # The length of the axis ticks is increased.
  axis.ticks.length.x = unit(.5, "lines"),
  axis.ticks.length.y = unit(.5, "lines"),
  # Remove the grid lines that come with ggplot2 plots by default
  # panel.grid = element_blank(),
  # Customize margin values (top, right, bottom, left)
  plot.margin = margin(5, 5, 5, 5),
  # Use a light grey color for the background of both the plot and the panel
  plot.background = element_rect(fill = "grey98", color = "grey98"),
  panel.background = element_rect(fill = "grey98", color = "grey98"),
  # Customize title appearence
  plot.title = element_text(
    color = "grey10", 
    size = 10, 
    face = "bold",
    margin = margin(t = 10,b = 10)
  ),
  # Customize subtitle appearence
  plot.subtitle = element_text(
    color = "grey30", 
    size = 8,
    lineheight = 1.35,
    margin = margin(t = 10, b = 20)
  ),
  # Title and caption are going to be aligned
  plot.title.position = "panel",
  plot.caption.position = "panel",
  plot.caption = element_text(
    color = "grey30", 
    size = 8,
    lineheight = 1.2, 
    hjust = 0,
    margin = margin(t = 20) # Large margin on the top of the caption.
  ),
  # Remove legend
  # legend.position = "none"
  # Change the background of the legend
  legend.background = element_rect(fill = "grey98", color = "grey98"),
  # Add a square around the legend
  # legend.box.background = element_rect(colour = "black")
)

# Create a second theme for the ggarrange plots
arrange_theme <- function() {
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.margin = margin(5, 5, 5, 5),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ),
        # Customize subtitle appearence
        plot.subtitle = element_text(
          color = "grey30", 
          size = 8,
          lineheight = 1.35,
          margin = margin(l = 10 ,t = 10, b = 20),
          hjust = 0.5
        ),
        # Title and caption are going to be aligned
        plot.caption = element_text(
          color = "grey30", 
          size = 8,
          lineheight = 1.2,
          hjust = 0.1,
          margin = margin(t = 10, b = 10) # Large margin on the top of the caption.
        ))
  
} # End of arrange_theme

# Message
cat(rule(left = "- Theme set - ", line_col = "white", line = " ", col = "green"))

  #### . Data Preparation . ##### 

# -- Bryophytes -- # 

# Load the datasets
  # Occurence
bryoData <- read.csv("FinalMatrixLv95_SynChangedV3.csv",row.names=1) 
  # Traits
bryoTraits <- read.csv("Bryophytes_Traits_SpRichness.csv",row.names=1)
  # Environment Variables
bryoEnv <- read.csv("EnvironmentalValues_Bryophytes.csv",row.names=1) # mnt_mean was changed directly in the dataset to "z"

# Load the phylotrees
Liver_Tree <- read.tree("PhyloTree/timetree50mod-liverwortsV2.nwk")
Mosses_Tree <- read.tree("PhyloTree/timetree50mod-mossesV2.nwk")

# -- Split between Mosses and Liverworts -- #

# Load the names of the mosses present in the dataframe of Occ_Data_Moss for further splitting between Mosses and Liverworts.
Bryo_Names <- read.csv(file = "Utilities/Bryophyte_Names.csv", row.names = 1)
Moss_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Moss"),Species) %>% as.matrix()
Liver_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Liverwort"),Species) %>% as.matrix()

# Extract the meta data names
Bryo_Meta <- colnames(bryoData)[1:4]

# Split the Occ_Data_Moss between Mosses and Liverworts
liverData <- dplyr::select(bryoData,any_of(c(Bryo_Meta,Liver_Names)))
mossData <- dplyr::select(bryoData,any_of(c(Bryo_Meta,Moss_Names)))

# Get the names of the liverworts and mosses present
Moss_Names <- colnames(mossData)[-c(1:4)]
Liver_Names <- colnames(liverData)[-c(1:4)]

# Message
cat(rule(left = "- Data loaded - ", line_col = "white", line = " ", col = "green"))

#-------------------------------------------#
##### SR COMPUTATION AND PLOT FILTERING #####
#-------------------------------------------#

# -- Computation of the wanted Metrics -- #

# Input a threshold for the minimum number of species mandatory for the analyses
Thresh <- 5

# Total Bryophytes
Bryo.filtered <- bryoData %>%
  # Compute the Species Richness
  mutate(SR = apply(bryoData[,5:ncol(bryoData)],1,sum), .after = y) %>%
  # Select the plots based on the Species richness threshold.
  dplyr::filter(SR > Thresh)

# Number of plot lost : 224 / Number of plot left: 433
cat(paste0("BRYOPHYTES - Number of plot left: ",nrow(Bryo.filtered)))

# Total Mosses
Moss.filtered <- mossData %>%
  # Compute the Species Richness
  mutate(SR = apply(mossData[,5:ncol(mossData)],1,sum), .after = y) %>%
  # Select the plots based on the Species richness threshold.
  dplyr::filter(SR > Thresh)

# Number of plot lost : 234 / Number of plot left: 423
cat(paste0("MOSSES - Number of plot left: ",nrow(Moss.filtered)))

# Total Liverworts
Liver.filtered <- liverData %>%
  # Compute the Species Richness
  mutate(SR = apply(liverData[,5:ncol(liverData)],1,sum), .after = y) %>%
  # Select the plots based on the Species richness threshold.
  dplyr::filter(SR > Thresh)

# Number of plot lost : 631 / Number of plot left: 26
cat(paste0("LIVERWORTS - Number of plot left: ",nrow(Liver.filtered)))

# Message
cat(rule(left = "- Data filtered based on species richness - ", line_col = "white", line = " ", col = "green"))

#------------------------------------#
##### ENVIRONMENTAL DATA MERGING #####
#------------------------------------#

# Joining of the environmental data and move the z column
Bryo.filtered <- left_join(Bryo.filtered,bryoEnv) %>% relocate(z,.after = y)
Moss.filtered <- left_join(Moss.filtered,bryoEnv) %>% relocate(z,.after = y)
Liver.filtered <- left_join(Liver.filtered,bryoEnv) %>% relocate(z,.after = y)

#----------------------#
##### RANGE CHOICE #####
#----------------------#

# Set the number breaks
breaks = 10

# Choose between weighted or unweighted altitudinal stripes "UW_Break" or "W_Break" for the analyses. 
Break <- "W_Break"   

# Breaks have to be changed from the intervals to the simple number
#if(Break == "UW_Break") {Break_S <- "UW_Stripe"
#} else if (Break == "W_Break") {Break_S <- "W_Stripe"}

# Total Bryophytes
Bryo.filtered <- mutate(Bryo.filtered,
  # Add the Taxa
  "Taxa" = "Bryophytes",
  # Add the unweighted break (Range) 
  "UW_Break" = cut(Bryo.filtered$z, breaks, dig.lab = 4), 
  # Add the weighted break (Range) 
  "W_Break" = cut_number(Bryo.filtered$z, breaks, dig.lab = 4),
  # Add the unweighted break (Number)  
  "UW_Stripe" = cut(Bryo.filtered$z, breaks, dig.lab = 4, labels = F), # Set labels = F to have categorical values as simple numbers. 
  # Add the weighted break (Number) 
  "W_Stripe" = cut_number(Bryo.filtered$z, breaks, dig.lab = 4, labels = F), 
  .after = z)

# Total Mosses
Moss.filtered <- mutate(Moss.filtered,
  # Add the Taxa
  "Taxa" = "Mosses",
  # Add the unweighted break (Range) 
  "UW_Break" = cut(Moss.filtered$z, breaks, dig.lab = 4), 
  # Add the weighted break (Range) 
  "W_Break" = cut_number(Moss.filtered$z, breaks, dig.lab = 4),
  # Add the unweighted break (Number)  
  "UW_Stripe" = cut(Moss.filtered$z, breaks, dig.lab = 4, labels = F), # Set labels = F to have categorical values as simple numbers. 
  # Add the weighted break (Number) 
  "W_Stripe" = cut_number(Moss.filtered$z, breaks, dig.lab = 4, labels = F), 
  .after = z)

# Total Liverworts
Liver.filtered <- mutate(Liver.filtered,
  # Add the Taxa
  "Taxa" = "Liverworts",
  # Add the unweighted break (Range) 
  "UW_Break" = cut(Liver.filtered$z, breaks, dig.lab = 4), 
  # Add the weighted break (Range) 
  "W_Break" = cut_number(Liver.filtered$z, breaks, dig.lab = 4),
  # Add the unweighted break (Number)  
  "UW_Stripe" = cut(Liver.filtered$z, breaks, dig.lab = 4, labels = F), # Set labels = F to have categorical values as simple numbers. 
  # Add the weighted break (Number) 
  "W_Stripe" = cut_number(Liver.filtered$z, breaks, dig.lab = 4, labels = F), 
  .after = z)


# Add the stripe distances 
# Liverworts_Beta <- Liverworts_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)
# Mosses_Beta <- Mosses_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)
# Bryophytes_Beta <- Bryophytes_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)
# Tracheo_Beta <- Tracheo_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)

# Message
cat(rule(left = "- Altitudinal stripes added - ", line_col = "white", line = " ", col = "green"))

#-------------------------------#
##### PLOT SPECIES RICHNESS #####
#-------------------------------#

# Create a global dataframe containing the results for all taxa. 
Total.filtered <- 
  # Bind the all the wanted taxa
  purrr::reduce(list(Bryo.filtered,Moss.filtered,Liver.filtered),full_join) %>%
  # Select only the wanted data
  dplyr::select("Sites_notre_num_rotation_":"SR") %>%
  # Group the data
  group_by(Taxa, .data[[Break]]) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  mutate(Sum = sum(SR), Count=n(), Var_Value = var(SR), Mean_Value = mean(SR)) %>%
  # Sort the values by "z" across for both Taxa.
  group_by(Taxa) %>%
  arrange(z, .by_group = TRUE)


# Compute and create the summary table of the metric(s) splitted between stripes for each taxa
SR.summary <- foreach(Taxon = unique(Total.filtered$Taxa)) %dopar% {

  # Draw the plot 
  Plot <- Total.filtered %>%
    # Select the wanted taxa
    dplyr::filter(Taxa == Taxon) %>%
    # Group the data
    group_by(.data[[Break]]) %>%
    # Compute the summary stat
    get_summary_stats(SR,type = "common") %>% # Compute the summary statistics for each groups 
    # Draw the summary statistics
    ggsummarytable(
      x = Break,                           # Split by stripes
      y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
      digits = 2,                            # Number of digits 
      size = 12,                              # Size of the text
      color = Break,                       # Color by stripes
      ggtheme = arrange_theme() +            # Theme
       theme(legend.position = "none")
      )
      
  # Return the plot  
  return(Plot)

} %>% set_names(unique(Total.filtered$Taxa))
    
# Compute and create the summary table of the metric(s) unsplitted between stripes
SR.summary.unsplitted <- Total.filtered %>%
# Group the data
group_by(Taxa) %>%
# Compute the summary stat
get_summary_stats(SR,type = "common") %>% # Compute the summary statistics for each groups 
# Draw the summary statistics
ggsummarytable(x = "Taxa",                          # Split by stripes
               y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
               digits = 2,                            # Number of digits 
               size = 12,                              # Size of the text
               ggtheme = arrange_theme() +            # Theme
                 theme(legend.position = "none")
)

    
# - Distribution of SR ~ Altitudinal stripes - # 

SR_Density <- foreach(Taxon = unique(Total.filtered$Taxa)) %dopar% {    

  # Find the count of plots in each stripe
SR.mean <- Total.filtered %>%
  # Select the wanted taxa
  dplyr::filter(Taxa == Taxon[2]) %>%
  # Group the data
  group_by(.data[[Break]] ) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  summarise(z = mean(z), Mean_Value = mean(SR)) %>%
    arrange(z)


Plot <- Total.filtered %>% 
  # Select the wanted taxa
  dplyr::filter(Taxa == Taxon[2]) %>%
  # Aes
  ggplot(aes(x = z, y = SR, color = .data[[Break]])) +
  # Plot all the values
  geom_point() +
  # Plot the connected scatter plot of the mean values for each altitudinal stripe.
  new_scale_color() + # Define scales before initiating a new one
  # Use the new scale
  geom_point(data = SR.mean, aes(y = Mean_Value), size = 1) +
  geom_line(aes(y = Mean_Value),linewidth = 1, alpha = 0.7, group = 1) +
  # Guides
  guides(size = "none") +
  # Labels
  xlab("Altitude") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("ScatterPlot of each metric values ~ Altitude"),
    subtitle = paste0(
      "Black dotted line represent the mean for each metric of the altitudinal stripes.\nAltitudinal stripes are",{ifelse(Break == "W_Break",paste(" weighted by plot number i.e represent an equal number of plots."),paste(" unweighted by plot number i.e represent an equal geographical distance."))}))

  # Return the plot
  return(Plot)
   
} %>% set_names(unique(Total.filtered$Taxa))


    
  # - Density plots of the Metric(s) values - #
    
Plot2_Glob <- Alpha_Data %>% 
  # Group the data
  group_by(Metric, Taxa) %>%
  # Aes
  ggplot(aes(x = Value)) +
  # Plot
  geom_density(fill = "grey",alpha = 0.5) +
  # Split between Taxa and Metric
  facet_grid(Metric ~ Taxa, scales = "free") +
  # Labels
  labs(
    title = paste0("DensityPlot of each Metrics"))
    
    # ----- #
    
  # - BoxPlot of the Metric(s) values - #

Plot3_Glob <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Aes
  ggplot(aes(x = .data[[Break_S]], y = Value, color = .data[[Break]], group = .data[[Break]])) +
  # Color Palette    
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_boxplot() + 
  # Split between Taxa and Metric
  facet_grid(Metric ~ Taxa, scales = "free_y") +
  # Add the significativity labels
  stat_pvalue_manual(WT_S_Adj, 
                     label = "p.adj.signif",
                     hide.ns = T) +
  # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcox-tests were realized between adjacent stripes and significant results are displayed.",
                          "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

    # ----- #

  # - BoxPlot of the Metric(s) values ~ Taxa - #

Plot4_1_Glob <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Aes
  ggplot(aes(x = Taxa, y = Value, color = Taxa)) +
  # Color Palette    
  scale_color_manual(values = color_scheme) +
  # Plot all the values
  geom_boxplot() + 
  # Split between Taxa and Metric
  facet_grid(Metric ~ .data[[Break_S]], scales = "free_y") +
  # Add the significativity labels
  # stat_pvalue_manual(WT_Taxa, 
  #                    label = "p.adj.signif",
  #                    hide.ns = T,
  #                    step.increase = 0.1,
  #                    step.group.by = c("Metric",Break_S)) +
  # Rotate
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # Labels
  xlab("Taxa") +
  ylab("Metric values") +
  labs(
    color = "Taxa",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes ~ Comparison of taxa "),
    subtitle = paste0("Wilcox-tests were realized between adjacent stripes and significant results are displayed.",
                      "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

  # - Connected scatterplots of the Metric(s) values ~ Taxa - #

Plot4_2_Glob <- Summ_Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Aes
  ggplot(aes(x = .data[[Break]], y = Mean_Value, color = Taxa, group = Taxa)) +
  # Color Palette    
  scale_color_manual(values = color_scheme) +
  # Plot all the values
  geom_point() + 
  geom_line() +
  # Split between Taxa and Metric
  facet_grid(Metric ~ ., scales = "free_y") +
  # Rotate
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Taxa",
    title = paste0("Scatterplots of metric values by taxa ~ Altitudinal stripes")
  )

  # - Combination of the two previous plots of the Metric(s) values ~ Taxa - #    

Plot4_Glob <- 
    # Arrange the two plots
    ggarrange(Plot4_1_Glob,Plot4_2_Glob, ncol = 2, legend = "right", common.legend = T, widths = c(3,1)) %>% 
  as.ggplot() +
  arrange_theme()

##### ------------ 1.A.b: Statistics ----------------------- #####

##### ----- 1.A.b.1: Normality tests ----- #####

# Verification of the normality of the distribution of the values for each altitudinal stripe.
# We'll use the of test of Shapiro-Wilk because we have most of the time less than 50 values.
# H0: Data are normally distributed: if p-value < 0.05, we reject that hypothesis.

# Message
cat(rule(left = "- Normality tests: Shapiro-Wilk - ", line_col = "white", line = "_", col = "grey"))

  #### . Shapiro-Wilk (Metric + Taxa) . ####

SW_US <- Alpha_Data %>% 
  group_by(Metric,Taxa)  %>% 
  do(tidy(shapiro.test(.$Value))) %>% 
  ungroup() %>% 
  select(-method) %>%
  # do an assignment (:=) and pass variables as column names by unquoting (!!) to not evaluate it
  mutate(!!Break_S := NA)  # This column is added for future joining with the next dataset. 

# Result:
cat(paste0("- Normality tests grouped by Metric and Taxa: ",sum(SW_US$p.value > 0.05),"/",nrow(SW_US)," tests passed. - "))
  
    # ---- #
  
  #### . Shapiro-Wilk (Metric + Taxa + Stripes) . ####

SW_S <- Alpha_Data %>% 
  group_by(Metric, Taxa, .data[[Break_S]]) %>%
  # Find the groups where the Variance is 0
  mutate(Var_Value = var(Value)) %>%
  # Remove these groups because they make the shapiro.test bug
  filter(Var_Value != 0) %>%
  # Do the shapiro_test
  do(tidy(shapiro.test(.$Value))) %>%
  ungroup() %>% 
  select(-method)

# Result:
cat(paste0("- Normality tests grouped by Metric and Taxa and Stripes: ",sum(SW_S$p.value > 0.05),"/",nrow(SW_S)," tests passed. - "))
  
# Merge the two datasets
SW_Result <- full_join(SW_US,SW_S)

  # ------------------------------------------- # 
  # /!\ Hypothesis of normality are not met /!\ # 
  # ------------------------------------------- # 

##### ----- 1.A.b.2: Equality of mean tests. ----- #####

  # We will use Kruskall-Wallis tests to determine if means of the metrics are different from each other. 

# --- Create a Dataset of Mean Values --- #

# Find the count of plots in each stripe
Summ_Alpha_Data <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa, .data[[Break_S]],.data[[Break]] ) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  summarise(z = mean(z), Mean_Value = mean(Value))

# --- Compute the significance values with Kruskall-Wallis tests --- #

  # Is there a significant difference between the groups ? 
  # If p-value <0.05, there is a significant difference between the groups. 
  
# Message
cat(rule(left = "- Equality of means tests: Kruskall-Wallis - ", line_col = "white", line = "_", col = "grey"))
  
  #### . Kruskall-Wallis (Metric + Taxa) . ####

KW_US <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Compute the kruskall-test
  do(tidy(kruskal.test(x = .$Value, g = .data[[Break_S]]))) %>%
  # Remove the method
  select(-c(method,parameter))

    # ----- # 

  #### . Kruskall-Wallis effect size (Metric + Taxa) . ####

  # Effect size tells you how meaningful the relationship between variables or the difference between groups is. It indicates the practical significance of a research outcome.
  # A large effect size means that a research finding has practical significance, while a small effect size indicates limited practical applications.
  # The value multiplied by 100 represent the percentage of variance of the dependent variable explained by the independent one. 
  
# Message
cat(rule(left = "- Kruskall-Wallis Effect-Size - ", line_col = "white", line = "_", col = "grey"))
  
KW_ES <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Compute the kruskall-test
  kruskal_effsize(formula = as.formula(paste("Value","~",Break_S))) %>%
  # Remove the method
  select(-c(method))
  
    # ----- # 

  #### . Wilcoxon test (Metric + Taxa ~ Stripes). ####
  
# Message
cat(rule(left = "- Wilcox test (Metric + Taxa ~ Stripes) - ", line_col = "white", line = "_", col = "grey"))
  
WT_S <- Alpha_Data %>%
  # Group and filter the data with variance == 0
  group_by(Metric,Taxa,!!Break_S) %>%
  filter(Var_Value != 0) %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Compute the kruskall-test
  wilcox_test(formula = as.formula(paste("Value","~",Break_S)), p.adjust.method = "bonferroni") %>%
  # Transform "group1" and "group2" into numeric
  mutate_at(c("group1", "group2"), as.numeric ) %>%
  # Add a column that is the "stripe distance" between the distribution compared
  mutate(Stripe_Distance = abs(group1 - group2))

    # ----- # 

  #### . Wilcoxon test / Adjacent stripes (Metric + Taxa ~ Stripes ). ####

# Message
cat(rule(left = "- . Wilcoxon test for adjacent stripes (Metric + Taxa ~ Stripes ). - ", line_col = "white", line = "_", col = "grey"))

# We need a column "y.position" for the plotting of the significance brackets
y.position <- Alpha_Data %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Add the y.position with an increase of X%. 
  summarise(y.position = max(Value) * 1.2) 
  
# Filter the precedent data_frame to only keep the values of stripe distance == 1
WT_S_Adj <- WT_S %>%
  # Filter the data 
  filter(Stripe_Distance == 1) %>%
  # Add the values of y.position 
  left_join(y.position, by=c("Metric","Taxa"))

    # ----- # 

  #### . Wilcoxon test (Metric + Stripes ~ Taxa). ####

# Message
cat(rule(left = "- Kruskall-Wallis (Metric + Stripes ~ Taxa ) - ", line_col = "white", line = "_", col = "grey"))

# We need a column "y.position" for the plotting of the significance brackets
y.position.taxa <- Alpha_Data %>%
  # Group the data
  group_by(Metric,.data[[Break_S]]) %>%
  # Add the y.position with an increase of X%. 
  summarise(y.position = max(Value) * 1.2) 

# Make the test
WT_Taxa <- Alpha_Data %>%
  # Group the data
  group_by(Metric,.data[[Break_S]]) %>%
  # Compute the kruskall-test
  wilcox_test(formula = as.formula(paste("Value","~","Taxa")), p.adjust.method = "bonferroni") %>%
  # Add the values of y.position 
  left_join(y.position.taxa, by = c("Metric",Break_S))

  # ----- #

# ##### ----- 1.A.b.3: Correlations ----- #####
# 
#   #### . Correlation (Metric ~ Taxa). ####
# 
# # Compute the correlation between the metrics results and the variables. 
# Cor_1a <- Alpha_Data %>% 
#   # Group by the metric and taxa
#   group_by(Metric,.data[[Break]]) %>%
#   # only select the wanted columns
#   select(all_of(c("Metric","Taxa",Break,"Value"))) %>%
#   # Nest the data
#     # The goal of that is to keep the outer-grouping (Metric and eventually the stripes) to work on the inner data_frames
#   nest() %>%
#   # On each of the "data", apply functions with map()
#   mutate(
#     correlations = map(data, 
#                        ~ pivot_wider(data = .x,
#                                      names_from = Taxa, 
#                                      values_from = Value))
#   ) %>%
#   # Unnest to return to a dataframe shape
#   unnest(correlations) %>%
#   # Remove the "data" column
#   select(-data) %>%
#   # Keep the rows that correspond the the correlation between the actual "Value" and the Environnemental Variables
#   filter(term == "Value") %>% 
#   # Remove the "Values" and "Terms" columns
#   select(-c(Value,term)) %>%
#   # transform into a long format. 
#   pivot_longer(!c(Taxa,Metric), names_to = "Variable", values_to = "Correlation") %>%
#   # Sort the correlation column by decreasing order
#   arrange(desc(Correlation),.by_group = TRUE)
# 
# # --- Compute the correlations between the metrics values and the different wanted variables --- #
# 
# # Create a global dataframe containing the results for all taxas as well as the environnemental variables
# Alpha_Cor_Data <- Data_List %>%
#   # Join everything together
#   reduce(full_join) %>%
#   # Combine the multiple columns of metric results into two : Metric and Value
#   gather(key = "Metric", value = "Value", !!Metric) %>%
#   # Transform into factor
#   mutate(across(!!Break_S, as_factor)) %>%
#   # Remove the columns that correspond to the species names
#   select(-all_of(Sp_Names$Sp_name)) %>% 
#   # Relocate Metric, Taxa and Values
#   relocate(any_of(c("Metric","Value")),.after  = Taxa) %>%
#   # Remove columns that contains Any NA
#   select_if(~ !any(is.na(.)))
# 
# # Compute the correlation between the metrics results and the variables. 
# Alpha_Cor <- Alpha_Cor_Data %>% 
#   # Group the data
#   group_by(Taxa,Metric) %>%
#   # Create a list column of each grouped dataframes
#   nest() %>%
#   # Apply the correlations on each of the groups
#   mutate(
#     correlations = map(data, correlate)
#   ) %>%
#   # Unnest to return to a dataframe shape
#   unnest(correlations) %>%
#   # Remove the "data" column
#   select(-data) %>%
#   # Keep the rows that correspond the the correlation between the actual "Value" and the Environnemental Variables
#   filter(term == "Value") %>% 
#   # Remove the "Values" and "Terms" columns
#   select(-c(Value,term)) %>%
#   # transform into a long format. 
#   pivot_longer(!c(Taxa,Metric), names_to = "Variable", values_to = "Correlation") %>%
#   # Sort the correlation column by decreasing order
#   arrange(desc(Correlation),.by_group = TRUE)
# 
# # Plot the correlations values
# Alpha_Cor_Plot <- Alpha_Cor %>%
#   # Group the data
#   group_by(Metric, Taxa) %>%
#   # Aes
#   ggplot(aes(x = fct_inorder(Variable), y = Correlation, color = Taxa)) +
#   # Palette 
#   scale_color_brewer(palette = Pal_col) +
#   # Plot all the values
#   geom_point() +
#   # Split between Taxa and Metric
#   facet_grid(Taxa ~ Metric , scales = "free_y") +
#   # Reoder the x-axis
#   # scale_x_discrete(limits=Variable)
#   # Rotate
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   # Guides
#   guides(size = "none") +
#   # Labels
#   xlab("Altitude") +
#   ylab("Metric values") +
#   labs(
#     color = "Stripe altitudinal limits",
#     title = paste0("ScatterPlot of each metric values ~ Altitude"),
#     subtitle = paste0(
#       "Black dotted line represent the mean for each metric of the altitudinal stripes.\nAltitudinal stripes are",{ifelse(Break == "W_Break",paste(" weighted by plot number i.e represent an equal number of plots."),paste(" unweighted by plot number i.e represent an equal geographical distance."))}))


##### ----- 1.A.2.a : Global Plots ------------------- #####

    # ----- # 



    # ----- #
    
  # - ScatterPlot of the Kruskall_Wallis tests results - #
    
# Create a dataset of significance for the next plots.
Group_Signif <- WT_S %>%
  # Group.by the wanted columns
  group_by(Metric,Taxa,Stripe_Distance,p.adj.signif) %>%
  # Change the stripe distance into factors
  mutate_at(vars(Stripe_Distance),factor) %>%
  # Create a column: Significant VS Non-Significant
  mutate(Significance = ifelse(startsWith(p.adj.signif, "*"), "Significant", "Non_Significant")) %>%
  # Change the grouping
  ungroup() %>% group_by(Metric,Taxa,Stripe_Distance) %>%
  # Add the number of stripe distances replicates for each of the stripes distances
  mutate(Stripe_count = n()) %>%
  # Change the grouping
  ungroup() %>% group_by(Metric,Taxa,Stripe_Distance,Significance) %>%
  # Add the number of significant results for each of the stripes distances
  mutate(Signif_count = n()) %>%
  # Compute the ratio of significant and non-significant tests (SPLITTED)
  mutate(Signif_Percent = (Signif_count/Stripe_count)*100)
    
# Plot the comparisons of CONSECUTIVE altitudinal stripes.
Plot5_Glob  <- Group_Signif %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Aes
  ggplot(aes(x = Stripe_Distance, y = Signif_Percent, group = Significance, color = Significance)) +
  # Split between Taxa and Metric
  facet_grid(Metric ~ Taxa, scales = "free_y") +
  # Draw the plot
  geom_point() + 
  geom_line() + 
  geom_hline(yintercept = 50) +
  # Color Palette
  scale_color_brewer(palette=Pal_col) +
  # Labels
  xlab("Stripe distances") +
  ylab("Percentage of tests") +
  labs(
    title = paste0("Percentage of Significant Wilcoxon tests of metric values ~ Stripe distance"),
    subtitle = "Stripe distance is an articifial unit to describe how far the stripes are from each other.\nThe higher the stripe distance, the lower the number of possible pairwise comparisons.\nAll significance levels are combined into the two modalities."
  )


#--------------------------#
##### PIST COMPUTATION #####
#--------------------------#






# Creation of a function to compute the phybeta-diversity indices.
# Take as parameter : 
# - Data : A dataframe containing the Occurrence data.
# - Sp_names : A character vector of species names

Beta_Phylo <- function(Data,Sp_names) {
  
  # ----- #
  
  # --- Creation of the dataframe --- # 
  
  Beta_Occ <- Data %>% 
      dplyr::select(all_of(as.vector(Sp_names))) # Remove the meta_data columns.

  
  # ----- #
  
  # --- Computation of the Taxonomic Beta-Diversity Metrics --- # 
  
  # Sorensen
  Beta_sor <- 
    beta.pair(Beta_Occ, index.family = "sorensen") %>%
    lapply(as.matrix)
  
  # Transform it into a vector
  Beta_sor <-
    lapply(names(Beta_sor), function(x) {
      PW_to_Vector(Beta_sor[[x]], Colname=x)}) %>%
    # Merge the results altogether
    purrr::reduce(merge) 
  
  # Jaccard
  Beta_jac <- 
    beta.pair(Beta_Occ, index.family = "jaccard") %>%
    lapply(as.matrix)
  
  # Transform it into a vector
  Beta_jac <-lapply(names(Beta_jac), function(x) {
    PW_to_Vector(Beta_jac[[x]], Colname=x)}) %>%
    # Merge the results altogether
    purrr::reduce(merge) 
  
  # Bind the two metrics together 
  Taxo_Beta <- merge(Beta_sor,Beta_jac)
  
  # ----- #
  
  # --- Addition of Meta-Data --- # 
  
  # Add the z-distance and 3D-distances between the plots as a variable. 
  Taxo_Beta <- 
    if(Type == "B"){ 
      
      # -- Plot A: Add x, y and z -- # 
      
      left_join(Taxo_Beta,dplyr::select(Data,"Site_VDP","x","y","z"), c("PlotA" = "Site_VDP")) %>% # Use a left join
        relocate(any_of(c('x',"y","z")),.after = PlotB) %>% # Move the column
        data.table::setnames(c('x',"y","z"),c('PlotA_x','PlotA_y','PlotA_z')) %>% # rename the column
        
        # -- Plot B: Add x, y and z -- # 
        
        left_join(dplyr::select(Data,"Site_VDP","x","y","z"), c("PlotB" = "Site_VDP")) %>% # Use a left join
        relocate(any_of(c('x',"y","z")),.after = PlotA_z) %>% # Move the column
        data.table::setnames(c('x',"y","z"),c('PlotB_x','PlotB_y','PlotB_z')) %>% # rename the column
        
        # -- Add the absolute difference of altitude between the two plots (delta_z)
        mutate(delta_z = abs(PlotA_z - PlotB_z), .after = PlotB_z) %>%
        # -- Add the complete 3D-distance between the two plots (delta_xyz)
        mutate(delta_xyz = sqrt((PlotB_x - PlotA_x)^2 + (PlotB_y - PlotA_y)^2 + (PlotB_z - PlotA_z)^2), .after = delta_z)
      
    } else {
      
      # -- Plot A: Add x, y and z -- # 
      
      left_join(Taxo_Beta,dplyr::select(Data,"X","x","y","z"), c("PlotA" = "X")) %>% # Use a left join
        relocate(any_of(c('x',"y","z")),.after = PlotB) %>% # Move the column
        data.table::setnames(c('x',"y","z"),c('PlotA_x','PlotA_y','PlotA_z')) %>% # rename the column
        
        # -- Plot B: Add x, y and z -- # 
        
        left_join(dplyr::select(Data,"X","x","y","z"), c("PlotB" = "X")) %>% # Use a left join
        relocate(any_of(c('x',"y","z")),.after = PlotA_z) %>% # Move the column
        data.table::setnames(c('x',"y","z"),c('PlotB_x','PlotB_y','PlotB_z')) %>% # rename the column
        
        # -- Add the absolute difference of altitude between the two plots (delta_z)
        mutate(delta_z = abs(PlotA_z - PlotB_z), .after = PlotB_z) %>%
        # -- Add the complete 3D-distance between the two plots (delta_xyz)
        mutate(delta_xyz = sqrt((PlotB_x - PlotA_x)^2 + (PlotB_y - PlotA_y)^2 + (PlotB_z - PlotA_z)^2), .after = delta_z)
    }
  
  # Return the wanted data
  return(Taxo_Beta)
  
}