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
       ggplot2)

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