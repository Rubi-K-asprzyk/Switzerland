#!/usr/bin/env Rscript

# --------------------------- #
# Swiss Patterns / K.Thibault #
# --------------------------- #

# This script is the step 0 of the Swiss Patterns project. It will create the global data frames used in the main script. 

# Working directory
setwd("~/Documents/PhD-Thesis/Research/Switzerland")

#------------------------#
##### INITIALISATION #####
#------------------------#

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})

# Install/load tons of packages.
p_load("doParallel", # Allow parallel computation
       "fields",     # Allow the display of the grid of plots
       "ape",        # Multiple tools of phylogenetic analyses
       "untb",       # Allow the creation of the fisherian ecosystem
       "TreeSim",    # Allow the creation of the phylogenetic trees.
       "dplyr",      # Allow the use of the plyr syntax
       "tidyr",      # Allow the use of function "pivot_longer"
       "stringr",    # Allow the use of function "str_remove"
       "ggplot2",    # Graphical representation tools
       "magrittr",   # Allow the use of function "set_colnames"
       "argparser",   # Add a command_line parser 
       "tidyverse",
       "vegan",
       "raster",
       "ggplotify",
       "ggpubr",
       "patchwork",
       "RColorBrewer",
       "scales",
       "betapart",
       "purrr",
       "gridExtra",
       "ggpointdensity",
       "hexbin",
       "gdm",
       "scales",
       "FactoMineR",
       "factoextra",
       "broom",
       "spaa",
       "LaplacesDemon",
       "rstatix",
       "ufs",
       "ggnewscale",
       "multcompView",
       "varhandle",
       "caret",
       "corrr",
       "ggtree",
       "cli"
)

#### . Plot Theme . ####

# Print a message
cli_alert_info("Setting the theme ... ")

# This theme extends the 'theme_light' that comes with ggplot2.
# The "Lato" font is used as the base font. This is similar
# to the original font in Cedric's work, Avenir Next Condensed.
# theme_set(theme_light(base_family = "Lato"))
theme_set(theme_light())

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
  plot.margin = margin(10, 10, 10, 10),
  # Use a light grey color for the background of both the plot and the panel
  plot.background = element_rect(fill = "grey98", color = "grey98"),
  panel.background = element_rect(fill = "grey98", color = "grey98"),
  # Customize title appearence
  plot.title = element_text(
    color = "grey10",
    size = 10,
    face = "bold",
    margin = margin(t = 10)
  ),
  # Customize subtitle appearence
  plot.subtitle = element_text(
    color = "grey30",
    size = 8,
    lineheight = 1.35,
    margin = margin(t = 10, b = 20)
  ),
  # Title and caption are going to be aligned
  plot.title.position = "plot",
  plot.caption.position = "plot",
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

##### ----- Load the data ----- #####

# Load the environmental data (Bryo)
Env_Data_Bryo <- read.csv(file = "Env_Var_Bryo.csv",row.names = 1)

# Load the occurrence data (Bryo)
Occ_Data_Bryo <- read.csv(file = "FinalMatrixLv95_SynChanged.csv", row.names = 1)

# Load the environmental data (Tracheo)
Env_Data <- read.csv(file = "Env_Var_Tracheo.csv",row.names = 1)

# Load the occurrence data (Tracheo)
Occ_Data_Tracheo <- read.csv(file = "Occurence_Tracheo.csv", row.names = 1)

# Load the names of the mosses present in the dataframe of Occ_Data_Moss for further splitting between Mosses and Liverworts.
Bryo_Names <- read.csv(file = "Bryophytes/Bryophyte_Names.csv", row.names = 1)
Moss_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Moss"),Species) %>% as.matrix()
Liver_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Liverwort"),Species) %>% as.matrix()

# Extract the meta data names
Bryo_Meta <- colnames(Occ_Data_Bryo)[1:4]

# ----- RANGES ----- # 

# Add the range to the environmental dataframe
# [0,1000[ / [1000,1400[ / [1400,1800[ / [1800,2000[ / [2000,2200[ / [2200,+Inf[ 

Test <- Env_Data %>%
  mutate(Range = case_when( 
    z >= 0 & z < 1000 ~ 1,
    z >= 1000 & z < 1400 ~ 2,
    z >= 1400 & z < 1800 ~ 3,
    z >= 1800 & z < 2000 ~ 4,
    z >= 2000 & z < 2200 ~ 5,
    z >= 2200 ~ 6))

#-----------------------------------------#
#####  Creation of a global dataframe #####
#-----------------------------------------#

# ----- BRYOPHYTES ----- # 

# Split the Occ_Data_Moss between Mosses and Liverworts
Liver_Occ <- dplyr::select(Occ_Data_Bryo,all_of(c(Bryo_Meta,Liver_Names)))
Moss_Occ <- dplyr::select(Occ_Data_Bryo,all_of(c(Bryo_Meta,Moss_Names)))

# Create two global data frames for Mosses and Liverworts
Liver_Total <- inner_join(Liver_Occ,Test)
Moss_Total <- full_join(Env_Data_Bryo,Moss_Occ)

# Save the two data frames 
write_csv(Liver_Total,"Bryophytes/Liverworts_Total.csv")
write_csv(Moss_Total,"Bryophytes/Mosses_Total.csv")

#--------------------------#
#####  Check Tree tips #####
#--------------------------#

# Load the phylotrees
Liver_Tree <- read.tree("Bryophytes/timetree50mod-liverworts.nwk")
Mosses_Tree <- read.tree("Bryophytes/timetree50mod-mosses.nwk")

# Get the tips clean
Liver_Tips <- Liver_Tree$tip.label %>% gsub("'","",.)
Moss_Tips <- Mosses_Tree$tip.label %>% gsub("'","",.)  

# Find the species present or not in each tree
Int_Liver_Tips <-  intersect(Liver_Tips,Liver_Names)   # Found species present in the tree
Miss_Liver_Tips <- Liver_Names[!Liver_Names %in% intersect(Liver_Tips,Liver_Names)] # Found species not present in the tree

Int_Moss_Tips <-  intersect(Moss_Tips,Moss_Names) # Found species present in the tree
Miss_Moss_Tips <- Moss_Names[!Moss_Names %in% intersect(Moss_Tips,Moss_Names)] # Found species not present in the tree

# Saving a bunch of things to create a Excel for Alain.
# write_csv(data.frame(Liver_Tips),"Bryophytes/Liver_Tips.csv")
# write_csv(data.frame(Moss_Tips),"Bryophytes/Moss_Tips.csv")
# write_csv(data.frame(Miss_Liver_Tips),"Bryophytes/Miss_Liver_Tips.csv")
# write_csv(data.frame(Miss_Moss_Tips),"Bryophytes/Miss_Moss_Tips.csv")



