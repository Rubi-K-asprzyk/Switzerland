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

{

# Set Working Directory
setwd("/home/thibault/Documents/PhD-Thesis/Research/Switzerland")

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
       ggnewscale,
       gtools,
       forcats,
       rstatix,
       betapart,
       tibble,
       spacodiR,
       vegan,
       abdiv)

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
bryoData <- read.csv("FinalMatrixLv95_SynChangedV3.csv",row.names=1) %>%
  # Rename "Sites_notre_num_rotation_"
  rename("Site_VDP" = "Sites_notre_num_rotation_")
  # Traits
bryoTraits <- read.csv("Bryophytes_Traits_SpRichness.csv",row.names=1)
  # Environment Variables
bryoEnv <- read.csv("EnvironmentalValues_Bryophytes.csv",row.names=1) %>% # mnt_mean was changed directly in the dataset to "z"
  # Rename "Sites_notre_num_rotation_"
  rename("Site_VDP" = "Sites_notre_num_rotation_")
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

} # END OF INITIALISATION

#--------------------------------------------------#
##### RANGE-ANALYSES: CREATION OF THE FUNCTION #####
#--------------------------------------------------#

  # . Creation of a function "range-analyses" to realize all the wanted analyses for each taxa and threhold wanted . # 

  # - ARGUMENTS - #
    # Data: A community matrix of occurence with plot in rows and species in column (with of without metadata associated).
    # Sp_Names: A vector of species names. 
    # Phylo_Tree : Phylogenetic Tree associated with the dataset. 
    # Threshold: The minimum number of species wanted in each plot. (Default = 2).
    # Alt_limits: A vector of altitudinal limits to split the altitudinal gradient into stripes. (Default = NULL).
    # N_stripes: If {Alt_limits == NULL}, the number of stripes to split the altitudinal gradient. (Default = 6).
    # Weighted: If {Alt_limits == NULL}, Are the stripes weighted by the number of plots (i.e with varying size) or not (i.e varying in the number of plots). (Default = T). 

  # - OUTPUTS - #

  # - TEST - #
  Data <- mossData
  Sp_Names <- Moss_Names
  Phylo_Tree <- Mosses_Tree
  Threshold <- 4 # Mosses only accept a threshold of 4
  Alt_limits <- c(0,1000,1400,1800,2000,2200,Inf) # That makes 6 altitudinal # That makes 6 altitudinal ranges  /  Alt_limits <- c(1400,2000) # For the Liverworts
  N_stripes <- 6
  Weighted <- TRUE

  # - COMPUTED METRICS - # 
    # - Species Richness (SR)
    # - Gini-Simpson Index (GS)
    # - Phylogenetic Diversity (PD)
    # - Mean Phylogenetic Diversity (MPD)
    # - Mean Nearest Neighbour Distance (MNTD)

range_analyses <- function(Data, Sp_Names, Threshold, Alt_limits = NULL, N_stripes = 6, Weighted = TRUE) {

  # ----- METRIC COMPUTATION ----- #

  # Computation of the metrics and filtering of the data based on the specie richness Threshold
  Data.filtered <- Data %>%
    # Compute the Species Richness
    mutate(SR = apply(Data[,Sp_Names], 1, sum), .after = y) %>%
    # Compute the Gini-Simpson Index
    mutate(GS = apply(Data[,Sp_Names], 1, simpson), .after = y) %>%
    # Select the plots based on the Species richness threshold.
    dplyr::filter(SR > Threshold) %>%
    # Compute the Phylogenetic Diversity
    mutate(PD = picante::pd(samp = .[,Sp_Names], tree = Phylo_Tree)$PD, .after = y)

    # Compute and reorder the cophenetic distances
    Data.Diss <- cophenetic(Phylo_Tree)
    # Only keep the species that are present in the tree
    Sp_Names_Present <- Sp_Names[Sp_Names %in% colnames(Data.Diss)]
    # Reorder the distance matrix
    Data.Diss <- Data.Diss[Sp_Names_Present,Sp_Names_Present]
  
  # Add the MPD to the dataframe
  Data.filtered <- Data.filtered %>%
    # Compute the mean phylogenetic distance
    mutate(MPD = picante::mpd(samp = .[,colnames(Data.Diss)], dis = Data.Diss), .after = y) %>%
    # Compute the mean nearest neighbour distance
    mutate(MNTD = picante::mntd(samp = .[,colnames(Data.Diss)], dis = Data.Diss), .after = y) %>%
    # Pivot longer the metrics
    pivot_longer(cols = c(SR,GS,PD,MPD,MNTD), values_to = "Value", names_to = "Metric")

  # Message
  cat(rule(left = paste0("- Data filtered based on species richness / THRESHOLD = ",Threshold," - "), line_col = "white", line = " ", col = "green"))
  cat(paste0(" - Number of initial plots: ",length(unique(Data$Site_VDP))))
  cat(paste0(" - Number of plot left: ",length(unique(Data.filtered$Site_VDP))))
  cat(paste0(" - Number of plot eliminated: ",(length(unique(Data$Site_VDP)) - length(unique(Data.filtered$Site_VDP)))))

  # ----- ENVIRONMENTAL DATA MERGING ----- #####

  # Joining of the environmental data and move the z column
  Data.filtered <- left_join(Data.filtered,bryoEnv) %>% relocate(z,Metric,Value,.after = y)

  # --------------------------------------------------------------------------------------------------------------------------------- #

  # ----- RANGE CHOICE ----- #

    # Creation of the ranges based on the parameters entered in the function.

  # If Alt_limits are entered, create the dataframe based on them.
  if (is.null(Alt_limits) == FALSE) {
    # Create a dataframe of intervals as well as corresponding stripe number
    Range_Intervals <- paste0("[", Alt_limits[-length(Alt_limits)], "-", Alt_limits[-1], "[") %>%
      # Transform as dataframe
      data.frame(.) %>%
      dplyr::rename("Range_Intervals" = ".") %>%
      # Add the interval number corresponding
      cbind(Range_Number = 1:nrow(.))

    # Create a fonction to find in which interval the altitude "Z" falls. Z is a unique value (not a whole vector).
    Interval.finder <- function(Z) {
      # Scan all intervals.
      for (Int in 1:(length(Alt_limits) - 1)) {
        # Find in which interval the Z fall.
        if (Z >= Alt_limits[Int] & Z < Alt_limits[Int + 1]) {
          # Return "Int", i.e the number of the interval
          return(Int)
        }
      }
    }

    # Apply this function on all plots altitude.
    Range_Number <- foreach(Z = Data.filtered$z, .combine = c) %do% {
      # Apply the function
      Interval.finder(Z = Z)
    }

    # Add the range to the data
    Data.filtered <- Data.filtered %>%
      # Add it
      mutate(Range_Number = Range_Number, .after = z) %>%
      # Add the range intervals.
      right_join(Range_Intervals, by = "Range_Number") %>%
      # Move the column
      relocate(Range_Intervals, .before = Range_Number)

    # If Alt_limits are not entered, create N_stripes, weighted or not by the number of plots they contains
  } else {
    # Create the ranges
    Data.filtered <- mutate(Data.filtered,
      # Add the desired Range_Intervals
      Range_Intervals = case_when(Weighted == TRUE ~ cut_number(Data.filtered$z, n = N_stripes, dig.lab = 4, right = FALSE), # Right = FALSE is used to have all breaks that start with a "[" to help for the correct ordering of the levels of the factor.
        Weighted == FALSE ~ cut(Data.filtered$z, N_stripes, dig.lab = 4, right = FALSE),
        .default = "OSKOUR"
      ),
      # Add the desired Range_Number
      Range_Number = case_when(Weighted == TRUE ~ as.character(cut_number(Data.filtered$z, n = N_stripes, dig.lab = 4, labels = F, right = FALSE)),
        Weighted == FALSE ~ as.character(cut(Data.filtered$z, N_stripes, dig.lab = 4, labels = F, right = FALSE)),
        .default = "OSKOUR"
      ),
      .after = z
    )
  }

# Message
cat(rule(left = "- Altitudinal stripes added - ", line_col = "white", line = " ", col = "green"))

  # --------------------------------------------------------------------------------------------------------------------------------- #

  # ----- SUMMARISATION BY RANGE ----- #

# Create a global dataframe containing the results for all taxa. 
Data.filtered.summarize <- Data.filtered %>%
  # Select only the wanted data
  dplyr::select("Site_VDP":"Value") %>%
  # Group the data
  group_by(Range_Number, Metric) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  mutate(Sum_Value = sum(Value), Mean_Value = mean(Value),Var_Value = var(Value), Count=n()) %>%
  # -!- Reorder the breaks -!- #
  group_by(Metric) %>% 
  # Arrange by the altitude to order in the same time the breaks
  arrange(z, .by_group = TRUE) %>%
  mutate(Break = factor(Range_Number, levels = unique(.$Range_Number)))

  
# Compute and create the summary table of the metric(s) splitted between stripes for each taxa

Alpha.summary <- Total.filtered %>%
  # Group the data
  group_by(Taxa, Break, Metric) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(
    x = "Break",                           # Split by stripes
    y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
    digits = 2,                            # Number of digits 
    size = 4,                             # Size of the text
    color = "Break",                       # Color by stripes
    facet.by = c("Metric","Taxa"),
    scales = "free_x", 
    ggtheme = arrange_theme() +            # Theme
     theme(legend.position = "none")
    ) +
  # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0))
 
    
# Compute and create the summary table of the metric(s) unsplitted between stripes
Alpha.summary.unsplitted <- Total.filtered %>%
# Group the data
group_by(Taxa, Metric) %>%
# Compute the summary stat
get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
# Draw the summary statistics
ggsummarytable(x = "Taxa",                          # Split by stripes
               y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
               digits = 2,                            # Number of digits 
               size = 5,                             # Size of the text
               facet.by = c("Metric","Taxa"),
               scales = "free_x", 
               ggtheme = arrange_theme() +            # Theme
                 theme(legend.position = "none")
)

# - Distribution of SR points ~ Altitudinal stripes - # 

# Find the mean value of z for x_axis centering and Value for y-axis
Value.Mean <- Total.filtered %>%
  # Group the data
  group_by(Taxa,Metric,Break) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  summarise(z = mean(z), Mean_Value = mean(Value)) %>%
  # reorder z by group to have the breaks in the good order
  arrange(z, .by_group = TRUE) %>%
  # Transform the breaks into factors
  mutate(Break = factor(Break))

# Draw the plots.
Alpha.scatter <- Total.filtered %>% 
  # Aes
  ggplot(aes(x = z, y = Value, color = Break)) + # Reorder the factors correctly
  # Plot all the values
  geom_point() +
  # Facet the plot
  facet_grid(Metric ~ Taxa,
            scales = "free_y") +
  # Plot the connected scatter plot of the mean values for each altitudinal stripe.
  new_scale_color() + # Define scales before initiating a new one
  # Use the new scale
  geom_point(data = Value.Mean, aes(y = Mean_Value), size = 1) +
  geom_line(aes(y = Mean_Value),linewidth = 1, alpha = 0.7, group = 1) +
  # Guides
  guides(size = "none") +
  # Labels
  xlab("Altitude") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("ScatterPlot of metrics ~ Altitude"),
    subtitle = paste0(
      "Black dotted line represent the mean metric value for each of the altitudinal stripes.\nAltitudinal stripes are",{ifelse(breaks_type == "Weighted",paste(" weighted by plot number i.e represent an equal number of plots."),paste(" unweighted by plot number i.e represent an equal geographical distance."))}))

# - Boxplots of Metrics Values ~ Altitudinal stripes - # 

# FIRST STEP / Compute the wilcoxon tests between all the stripes to after display them on the boxplot.
  
Metric.Wilcoxon <- Total.filtered %>%
  # Group and filter the data with variance == 0
  group_by(Taxa,Break,Metric) %>%
  filter(Var_Value != 0) %>%
  # Group the data
  group_by(Taxa,Metric) %>%
  # Compute the kruskall-test
  wilcox_test(formula = Value ~ Stripe, p.adjust.method = "bonferroni") %>%
  # Transform "group1" and "group2" into numeric
  mutate_at(c("group1", "group2"), as.numeric ) %>%
  # Add a column that is the "stripe distance" between the distribution compared
  mutate(Stripe_Distance = abs(group1 - group2))

# We need a column "y.position" for the plotting of the significance brackets
y.position <- Total.filtered %>%
  # Group the data
  group_by(Taxa,Metric) %>%
  # Add the y.position with an increase of X%. 
  summarise(y.position = max(Value) * 1.05) 
  
# Filter the precedent data_frame to only keep the values of stripe distance == 1
Metric.Wilcoxon.Adj <- Metric.Wilcoxon %>%
  # Filter the data to only keep adjacent stripes eventually 
  # filter(Stripe_Distance == 1) %>%
  # Add the values of y.position 
  left_join(y.position)

# We will use ddboxplot instead of geom boxplot because it is easier to plot the significance brackets on each facet. 

Alpha.Boxplot <- Total.filtered %>%
  # Draw the plot
  ggboxplot(
    x = "Break", y = "Value", fill = "Break",
    facet = c("Metric", "Taxa"),
    scales = "free",
    ggtheme = arrange_theme()
    ) + 
  # Add the significance levels
  # stat_pvalue_manual(Metric.Wilcoxon.Adj, hide.ns = TRUE, step.increase = 0.1, step.group.by = "Metric") +
    # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcox-tests were realized between adjacent stripes and significant results are displayed.",
                          "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  ) + # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0))

#-------------------------------#
##### BETA METRICS COMPUTATION #####
#-------------------------------#

Taxo.beta <- foreach(Taxa = list(Moss.filtered,Liver.filtered)) %dopar% { 

  # Add the rowname as a column to keep this info
  Data <- Taxa %>%
    # Pivot back the Alpha metrics computed to remove them (Or all the lines will be multiplied)
    pivot_wider(names_from = Metric, values_from = Value) %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric)


  # Get the Metadata
  MetaData <- Data %>% 
    # Select only the Metadata (based on the species names )
    dplyr::select(!any_of(c(Moss_Names,Liver_Names)))

  # -- Sorensen -- # 

  Sorensen <- Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Compute the Sorensen index
    beta.pair(x = ., index.family = "sorensen") %>%
    # Transform into a matrix
    lapply(as.matrix) 

      # Transform it into a vector
  Sorensen <-
    lapply(names(Sorensen), function(x) {
    PW_to_Vector(Sorensen[[x]], Colname=x)}) %>%
    # Merge the results altogether
    purrr::reduce(merge) %>%
    # Add back the taxa we worked with
    mutate(Taxa = unique(Data$Taxa), .before = 1)

  # -- Jaccard -- # 

  Jaccard <- Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Compute the Sorensen index
    beta.pair(x = ., index.family = "jaccard") %>%
    # Transform into a matrix
    lapply(as.matrix) 

      # Transform it into a vector
  Jaccard <-
    lapply(names(Jaccard), function(x) {
    PW_to_Vector(Jaccard[[x]], Colname=x)}) %>%
    # Merge the results altogether
    purrr::reduce(merge) %>%
    # Add back the taxa we worked with
    mutate(Taxa = unique(Data$Taxa), .before = 1)

  # Change the metadata to NOT contains the alpha metrics results
  MetaData <- dplyr::select(MetaData, !c(SR,GS,PD))

  # -- Bind the metrics together -- #
  Taxo_Beta <- merge(Sorensen,Jaccard) %>%
    # Transform the wanted columns into numeric
    mutate_at(c('PlotA', 'PlotB'), as.numeric) %>%
    # Join the metadata for both plots used in the pairwise computation. 
    left_join(MetaData, c("PlotA" = "Plot", "Taxa")) %>%  # Plot A
    left_join(MetaData, c("PlotB" = "Plot", "Taxa"), suffix = c("_A", "_B")) %>% # Plot B + add suffixes "A" and "B" to distinguish the two plots.
    # -- Add the absolute difference of altitude between the two plots (delta_z)
    mutate(delta_z = abs(z_A - z_B), .after = z_B) %>%
    # - Add the complete 3D-distance between the two plots (delta_xyz) - #
    mutate(delta_xyz = sqrt((x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2), .after = delta_z)

  # Return the results
  return(Taxo_Beta)

} %>% set_names(c("Mosses","Liverworts"))

# Combine both mosses and liverworts results
Taxo.beta <- purrr::reduce(list(Taxo.beta$Mosses,Taxo.beta$Liverworts),full_join) %>%
  # Pivot the metrics
  pivot_longer(cols = c(beta.sim,beta.sne,beta.sor,beta.jtu,beta.jne,beta.jac), values_to = "Value", names_to = "Metric") %>%
  # Move these columns
  relocate(Metric,Value,.after = PlotB)

#--------------------------------#
##### PHYLO BETA COMPUTATION #####
#--------------------------------#

Phylo.beta <- foreach(Taxa = list(Moss.filtered,Liver.filtered)) %dopar% {

  # Add the rowname as a column to keep this info
  Data <- Taxa %>%
      # Pivot back the Alpha metrics computed to remove them (Or all the lines will be multiplied)
    pivot_wider(names_from = Metric, values_from = Value) %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric) %>%
    # Select only the occurence data
    dplyr::select(any_of(c(Moss_Names,Liver_Names))) %>%
    # Transpose the data for the analysis
    t() %>%
    # Add the colnames as the plot numbers
    `colnames<-`(1:ncol(.))

  # Get the Metadata
  MetaData <- Taxa %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric) %>% 
    # Select only the Metadata (based on the species names )
    dplyr::select(!any_of(c(Moss_Names,Liver_Names,"Metric","Value")))

  # Select the adequate phylogenetic tree between the moss and liverworts ones
  ifelse(Taxa$Taxa == "Mosses", Phylo <- Mosses_Tree, Phylo <- Liver_Tree)

# -- Compute the PIst metrics -- # 

  Hardy_Metrics <- spacodiR::spacodi.calc(
      sp.plot = Data,         # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
      phy = Phylo,             # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                   # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
      all.together = TRUE,     # whether to treat all traits together or separately
      prune = TRUE,
      pairwise = TRUE)

# -- Extract all Metrics -- #
  
  ## PST ##

Pst <- Hardy_Metrics$pairwise.Pst %>% 
  # Change colnames and rownames to be a simple number (and not "plt.X")
  `colnames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
  `rownames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
  # Transform the pairwise 
  PW_to_Vector(Colname = "Pst") %>%
  mutate_at(c("PlotA","PlotB"), as.numeric) %>%
  # Join the metadata for both plots used in the pairwise computation. 
  left_join(MetaData, c("PlotA" = "Plot")) %>% # Plot A
  left_join(MetaData, c("PlotB" = "Plot", "Taxa"), suffix = c("_A", "_B")) %>% # Plot B + add suffixes "A" and "B" to distinguish the two plots. 
  # - Add the absolute difference of altitude between the two plots (delta_z) - #
  mutate(delta_z = abs(z_A - z_B), .after = z_B) %>%
  # - Add the complete 3D-distance between the two plots (delta_xyz) - #
  mutate(delta_xyz = sqrt((x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2), .after = delta_z) 

  ## BST ##

Bst <- Hardy_Metrics$pairwise.Bst %>% 
  # Change colnames and rownames to be a simple number (and not "plt.X")
  `colnames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
  `rownames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
  # Transform the pairwise 
  PW_to_Vector(Colname = "Bst") %>%
  mutate_at(c("PlotA","PlotB"), as.numeric) %>%
  # Join the metadata for both plots used in the pairwise computation. 
  left_join(MetaData, c("PlotA" = "Plot")) %>% # Plot A
  left_join(MetaData, c("PlotB" = "Plot", "Taxa"), suffix = c("_A", "_B")) %>% # Plot B + add suffixes "A" and "B" to distinguish the two plots. 
  # - Add the absolute difference of altitude between the two plots (delta_z) - #
  mutate(delta_z = abs(z_A - z_B), .after = z_B) %>%
  # - Add the complete 3D-distance between the two plots (delta_xyz) - #
  mutate(delta_xyz = sqrt((x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2), .after = delta_z) 

  ## PIST ##

PIst <- Hardy_Metrics$pairwise.PIst %>% 
  # Change colnames and rownames to be a simple number (and not "plt.X")
  `colnames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
  `rownames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
  # Transform the pairwise 
  PW_to_Vector(Colname = "PIst") %>%
  mutate_at(c("PlotA","PlotB"), as.numeric) %>%
  # Join the metadata for both plots used in the pairwise computation. 
  left_join(MetaData, c("PlotA" = "Plot")) %>% # Plot A
  left_join(MetaData, c("PlotB" = "Plot", "Taxa"), suffix = c("_A", "_B")) %>% # Plot B + add suffixes "A" and "B" to distinguish the two plots. 
  # - Add the absolute difference of altitude between the two plots (delta_z) - #
  mutate(delta_z = abs(z_A - z_B), .after = z_B) %>%
  # - Add the complete 3D-distance between the two plots (delta_xyz) - #
  mutate(delta_xyz = sqrt((x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2), .after = delta_z)
    
# --- Combine the three altogether and pivot longer the metrics --- # 

  # There is no use of Bst here because we have only occurrence data.

Phylo.beta <- purrr::reduce(list(Pst,PIst),full_join) %>%
  # Pivot the metrics
  pivot_longer(cols = c(Pst,PIst), values_to = "Value", names_to = "Metric") %>%
  # Relocate columns
  dplyr::select(order(colnames(.))) %>%
  relocate(any_of(c("Taxa","Metric","Value")),.before = 1) %>%
  relocate(any_of(c("delta_z","delta_xyz")),.after = z_B)

# Return the good values
return(Phylo.beta)

} %>% set_names(c("Mosses","Liverworts"))

# Combine both mosses and liverworts results
Phylo.beta <- full_join(Phylo.beta$Mosses,Phylo.beta$Liverworts)

# ------------- #

# Combine the taxonomic and phylogenetic beta results

Beta.results <- full_join(Phylo.beta,Taxo.beta)  %>%
  # Create a column that combines the two stripes
  mutate(Stripe_AB = paste0(Stripe_A,"-",Stripe_B)) %>%
  # Create a column that is the stripe distance between the two plots.
  mutate(Stripe_distance = abs(as.numeric(Stripe_A)-as.numeric(Stripe_B)))

# ------------- #

#---------------------------#
##### PLOT BETA RESULTS #####
#---------------------------#

# -- Boxplot of Metrics ~ Stripe distance -- #

# FIRST STEP / Compute the wilcoxon tests between all the stripe distances to after display them on the boxplot.
  
# Beta.Wilcoxon <- Beta.results %>%
#   # Group the data
#   group_by(Taxa,Metric) %>%
#   # Filter the dataset to remove the computation intra plot
#   filter(PlotA != PlotB) %>%
#   # Compute the kruskall-test
#   wilcox_test(formula = Value ~ Stripe_distance, p.adjust.method = "bonferroni") %>%
#   # Transform "group1" and "group2" into numeric
#   mutate_at(c("group1", "group2"), as.numeric)

# # We need a column "y.position" for the plotting of the significance brackets
# y.position <- Beta.results%>%
#   # Group the data
#   group_by(Taxa,Metric) %>%
#   # Add the y.position with an increase of X%. 
#   summarise(y.position = max(Value) * 1.05) 
  
# # Filter the precedent data_frame to only keep the adjacent stripe distances (1-2-3 ... )
# Beta.Wilcoxon.Adj <- Beta.Wilcoxon %>%
#   # Filter the data to have eventually only the adjacent stripes
#   filter(group2 == group1 + 1) %>%
#   # Add the values of y.position 
#   left_join(y.position, by=c("Taxa","Metric")) %>%
#   # WARNING: Brackets are moved to the left, therefore, we will move them to the right
#   mutate(group1 = group1 + 1) %>%
#   mutate(group2 = group2 + 1)

# SECOND STEP / Draw the boxplots.

Beta_Boxplot_StripeDistance <- Beta.results %>%
  # Draw the plot
  ggboxplot(
    x = "Stripe_distance", y = "Value", fill = "Stripe_distance",
    facet = c("Metric", "Taxa"),
    scales = "free",
    ggtheme = arrange_theme()
    ) + 
  # Add the significance levels
  # stat_pvalue_manual(Beta.Wilcoxon.Adj, hide.ns = TRUE, step.increase = 0, step.group.by = "Metric") +
    # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcox-tests were realized between adjacent stripes and significant results are displayed.",
                          "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

# -- Boxplot of Metrics ~ Intra Stripe -- #

# FIRST STEP / Compute the wilcoxon tests between all the stripe distances to after display them on the boxplot.
  
# Beta.Wilcoxon <- Beta.results %>%
#   # Filter the dataset to only keep the pairwise plots from the same stripe
#   filter(Stripe_A == Stripe_B) %>%
#   # Filter the dataset to remove the computation intra plot
#   filter(PlotA != PlotB) %>%
#   # Group the data
#   group_by(Taxa,Metric) %>%
#   # Compute the kruskall-test
#   wilcox_test(formula = Value ~ Stripe_AB, p.adjust.method = "bonferroni")

# # We need a column "y.position" for the plotting of the significance brackets
# y.position <- Beta.results %>%
#   # Filter the dataset to only keep the pairwise plots from the same stripe
#   filter(Stripe_A == Stripe_B) %>%
#   # Group the data
#   group_by(Taxa,Metric) %>%
#   # Add the y.position with an increase of X%. 
#   summarise(y.position = max(Value) * 1.05)
  
# # Filter the precedent data_frame to only keep the adjacent stripe distances (1-2-3 ... )
# Beta.Wilcoxon.Adj <- Beta.Wilcoxon %>%
#   # Filter the data 
#   filter(group2 == group1 + 1) %>%
#   # Add the values of y.position 
#   left_join(y.position, by=c("Taxa","Metric")) 

# SECOND STEP / Draw the boxplots.

Beta_Boxplot_IntraStripe <- Beta.results %>%
  # Filter the dataset to only keep the pairwise plots from the same stripe
  dplyr::filter(Stripe_A == Stripe_B) %>%
  dplyr::filter(PlotA != PlotB) %>%
  mutate_at("Stripe_AB", as.character) %>%
  # Reorder the stripes to have 10-10 effectively last
  mutate(Stripe_AB = fct_relevel(Stripe_AB,mixedsort(unique(.$Stripe_AB)))) %>%
  # Aes
  ggboxplot(
    x = "Stripe_AB", y = "Value", fill = "Stripe_AB",
    facet = c("Metric", "Taxa"),
    scales = "free",
    ggtheme = arrange_theme()
    ) + 
  # Add the significance levels
  # stat_pvalue_manual(Beta.Wilcoxon.Adj, hide.ns = TRUE, step.increase = 0, step.group.by = "Metric") +
    # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcox-tests were realized between adjacent stripes and significant results are displayed.",
                          "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

# -- Scatterplot of Metrics ~ z-distance -- #

# Draw the scatterplot.
  # It's a big one, we better not do it everytime.

# Beta_Scatterplot <- Beta.results %>%
#   # Select the wanted taxa (and remove intra plot computation)
#   dplyr::filter(PlotA != PlotB) %>%
#   mutate_at("Stripe_distance", as.character) %>%
#   # Aes
#   ggplot(aes(x = delta_z, y = Value, color = Stripe_distance, group = Stripe_distance)) +
#   # Plot all the values
#   geom_point() +
#   # facet by Metric
#   facet_grid(Taxa ~ Metric) +
#   # Labels
#   xlab("Z-Distance between plots.") +
#   ylab("Metric values.") +
#   labs(
#     color = "Stripe distances between the two plots.",
#     title = paste0("Scatterplot of metric values ~ Altitudinal distance")
#   )

#---------------------#
##### MANTEL TEST #####
#---------------------#

# Mantel Test are only used for beta metrics (to compare matrices against matrices)

# Compute a Mantel test between the PIst and the z-distance between the two plots used used in the comparison.

## We need to recompute the metrics without vectorizing them to compute mantel tests

Mantel.Data <- foreach(Taxa = list(Moss.filtered,Liver.filtered)) %dopar% { 

  # Add the rowname as a column to keep this info
  Data <- Taxa %>%
    # Pivot back the Alpha metrics computed to remove them (Or all the lines will be multiplied)
    pivot_wider(names_from = Metric, values_from = Value) %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric)

  # Get the Metadata
  MetaData <- Data %>% 
    # Select only the Metadata (based on the species names )
    dplyr::select(!any_of(c(Moss_Names,Liver_Names)))

  # -- Sorensen -- # 

  Sorensen <- Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Compute the Sorensen index
    beta.pair(x = ., index.family = "sorensen") %>%
    # Transform into a matrix
    lapply(as.matrix) 

  # -- Jaccard -- # 

  Jaccard <- Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Compute the Sorensen index
    beta.pair(x = ., index.family = "jaccard") %>%
    # Transform into a matrix
    lapply(as.matrix) 

  # -- Hardy Metrics -- #

  # Prepare the data
    Data_Hardy <- Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Transpose the data for the analysis
    t() %>%
    # Add the colnames as the plot numbers
    `colnames<-`(rownames(Data))

  # Select the good phylogenetic tree.
  Phylo <- Mosses_Tree
  if(unique(Data$Taxa) == "Liverworts") {Phylo <- Liver_Tree}

  # Compute the PIst
  Hardy_Metrics <- spacodiR::spacodi.calc(
    sp.plot = Data_Hardy,    # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
    phy = Phylo,       # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                   # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
    all.together = TRUE,     # whether to treat all traits together or separately
    prune = TRUE,
    pairwise = TRUE)

  # -- Extract all the Metrics we want to keep -- #

    # Hardy Metrics
  PIst <- Hardy_Metrics$pairwise.PIst %>% as.dist()
  Pst <- Hardy_Metrics$pairwise.Pst %>% as.dist()
  Bst <- Hardy_Metrics$pairwise.Bst %>% as.dist()
    # Sorensen Metrics
  beta.sim <- Sorensen$beta.sim %>% as.dist()
  beta.sne <- Sorensen$beta.sne %>% as.dist()
  beta.sor <- Sorensen$beta.sor %>% as.dist()
    # Jaccard Metrics
  beta.jtu <- Jaccard$beta.jtu %>% as.dist()
  beta.jne <- Jaccard$beta.jne %>% as.dist()
  beta.jac <- Jaccard$beta.jac %>% as.dist()

  # Return the results
  return(list(PIst,Pst,Bst,beta.jac,beta.jne,beta.jtu,beta.sim,beta.sne,beta.sor,Data$z) %>%
            setNames(c("PIst","Pst","Bst","beta.jac","beta.jne","beta.jtu","beta.sim","beta.sne","beta.sor","Zdist")))

} %>% set_names(c("Mosses","Liverworts"))

#### --- Computation of the Mantel Test --- # 

# Execute the tests on both Mosses & Liverworts

Mantel.Results <- foreach(i = 1) %do% {

  # Because only mosses are working for now 
  Data <- Mantel.Data[[i]]
  
  # Extract the Zdist and remove it from the Data
  Zdist <- dist(Data$Zdist)
  Taxa.Mantel <- Data
  Taxa.Mantel$Zdist <- NULL

  # Compute all the Mantel tests
  Mantel <- foreach(Data2 = Taxa.Mantel, .combine = rbind) %do% {

    # Compute the test
    Test <- mantel(Data2,Zdist, method = "spearman", permutations = 9, na.rm = TRUE)

    # Extract the name of the metric, the significance and the p-value
    Result <- data.frame(Statistic = Test$statistic, Pvalue = Test$signif)
    # Return the results

  } %>% mutate(Metric = names(Taxa.Mantel), .before = Statistic)

} # %>% set_names(names(Mantel.Data))

#-----------------------------#
##### PEARSON CORRELATION #####
#-----------------------------#

# Pearson correlations are used with alpha metrics against a variable (z)

Pearson_Correlation <- Total.filtered %>%
  group_by(Taxa, Metric) %>%
  summarise(correlation = cor(Value, z))


#-----------------------------#
##### TEST #####
#-----------------------------#

  Data <- Mantel.Data$Liverworts

  # Extract the Zdist and remove it from the Data
  Zdist <- dist(Data$Zdist)
  Taxa.Mantel <- Data
  Taxa.Mantel$Zdist <- NULL

  # Compute all the Mantel tests
  Mantel <- foreach(Data2 = Taxa.Mantel, .combine = rbind) %do% {

    # Compute the test
    Test <- mantel(Data2,Zdist, method = "spearman", permutations = 9, na.rm = TRUE)

    # Extract the name of the metric, the significance and the p-value
    Result <- data.frame(Statistic = Test$statistic, Pvalue = Test$signif)
    # Return the results

  }