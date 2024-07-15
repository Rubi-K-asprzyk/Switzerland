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
nrow(bryoData) - nrow(Bryo.filtered)

# Total Mosses
Moss.filtered <- mossData %>%
  # Compute the Species Richness
  mutate(SR = apply(mossData[,5:ncol(mossData)],1,sum), .after = y) %>%
  # Select the plots based on the Species richness threshold.
  dplyr::filter(SR > Thresh)

# Number of plot lost : 234 / Number of plot left: 423
nrow(mossData) - nrow(Moss.filtered)

# Total Liverworts
Liver.filtered <- liverData %>%
  # Compute the Species Richness
  mutate(SR = apply(liverData[,5:ncol(liverData)],1,sum), .after = y) %>%
  # Select the plots based on the Species richness threshold.
  dplyr::filter(SR > Thresh)

# Number of plot lost : 631 / Number of plot left: 26
nrow(liverData) - nrow(Liver.filtered)

# Message
cat(rule(left = "- Data filtered based on species richness - ", line_col = "white", line = " ", col = "green"))

#----------------------#
##### RANGE CHOICE #####
#----------------------#

# Choose between weighted or unweighted altitudinal stripes "UW_Break" or "W_Break"
Break <- "UW_Break"   

# Breaks have to be changed from the intervals to the simple number
if(Break == "UW_Break") {Break_S <- "UW_Stripe"
} else if (Break == "W_Break") {Break_S <- "W_Stripe"}

# Set the breaks
breaks = 10

# The tag *_Nb refers to labels that are not the interval but rather a simple number. 

# We need to have the same scales for both the bryophytes and the tracheophytes but their Z range are different.
# Therefore, we will add the global min and max(z) to both dataframes and remove them just after.
Added_Range <- c(range(Tracheo_Data$z, na.rm = T), range(Bryophytes_Data$z, na.rm = T))

    # ----- # 

  # --- Taxonomic Alpha Data --- # 

    # -- Bryophytes, Mosses, Liverworts -- # 

# Un-weighted breaks 
UW_Break_Bryo <- cut(c(Bryophytes_Data$z,Added_Range), breaks, dig.lab = 4) %>% head(-4) # Remove the last 4 values using head(-X)
UW_Break_Bryo_Nb <- cut(c(Bryophytes_Data$z,Added_Range), breaks, dig.lab = 4, labels = F)  %>% head(-4) # Set labels = F to have categorical values as simple numbers. 

# Weighted breaks
W_Break_Bryo <- cut_number(c(Bryophytes_Data$z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
W_Break_Bryo_Nb <- cut_number(c(Bryophytes_Data$z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 

# Add it to the initial data
Liverworts_Data <- mutate(Liverworts_Data,"Taxa" = "Liverworts", "UW_Break" = UW_Break_Bryo, "W_Break" = W_Break_Bryo, "UW_Stripe" = UW_Break_Bryo_Nb, "W_Stripe" = W_Break_Bryo_Nb, .after = z)
Mosses_Data <- mutate(Mosses_Data,"Taxa" = "Mosses", "UW_Break" = UW_Break_Bryo, "W_Break" = W_Break_Bryo, "UW_Stripe" = UW_Break_Bryo_Nb, "W_Stripe" = W_Break_Bryo_Nb, .after = z)
Bryophytes_Data <- mutate(Bryophytes_Data,"Taxa" = "Bryophytes", "UW_Break" = UW_Break_Bryo, "W_Break" = W_Break_Bryo, "UW_Stripe" = UW_Break_Bryo_Nb, "W_Stripe" = W_Break_Bryo_Nb, .after = z)


    # ----- #  

  # --- Taxonomic Beta Data --- # 

    # -- Bryophytes, Mosses, Liverworts -- # 

# Un-weighted breaks 
UW_Break_Bryo_A <- cut(c(Bryophytes_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
UW_Break_Bryo_B <- cut(c(Bryophytes_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 

UW_Break_Bryo_A_Nb <- cut(c(Bryophytes_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 
UW_Break_Bryo_B_Nb <- cut(c(Bryophytes_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 

# Weighted breaks
W_Break_Bryo_A <- cut(c(Bryophytes_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
W_Break_Bryo_B <- cut(c(Bryophytes_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 

W_Break_Bryo_A_Nb <- cut_number(c(Bryophytes_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 
W_Break_Bryo_B_Nb <- cut_number(c(Bryophytes_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 

# Add them to the initial data 
Liverworts_Beta <- mutate(Liverworts_Beta,"Taxa" = "Liverworts", "UW_BreakA" = UW_Break_Bryo_A,"UW_BreakB" = UW_Break_Bryo_B, "W_BreakA" = W_Break_Bryo_A,"W_BreakB" = W_Break_Bryo_B, "UW_StripeA" = UW_Break_Bryo_A_Nb,"UW_StripeB" = UW_Break_Bryo_B_Nb, "W_StripeA" = W_Break_Bryo_A_Nb,"W_StripeB" = W_Break_Bryo_B_Nb, .after = delta_xyz)
Mosses_Beta <- mutate(Mosses_Beta,"Taxa" = "Mosses", "UW_BreakA" = UW_Break_Bryo_A,"UW_BreakB" = UW_Break_Bryo_B, "W_BreakA" = W_Break_Bryo_A,"W_BreakB" = W_Break_Bryo_B, "UW_StripeA" = UW_Break_Bryo_A_Nb,"UW_StripeB" = UW_Break_Bryo_B_Nb, "W_StripeA" = W_Break_Bryo_A_Nb,"W_StripeB" = W_Break_Bryo_B_Nb, .after = delta_xyz)
Bryophytes_Beta <- mutate(Bryophytes_Beta,"Taxa" = "Bryophytes", "UW_BreakA" = UW_Break_Bryo_A,"UW_BreakB" = UW_Break_Bryo_B, "W_BreakA" = W_Break_Bryo_A,"W_BreakB" = W_Break_Bryo_B, "UW_StripeA" = UW_Break_Bryo_A_Nb,"UW_StripeB" = UW_Break_Bryo_B_Nb, "W_StripeA" = W_Break_Bryo_A_Nb,"W_StripeB" = W_Break_Bryo_B_Nb, .after = delta_xyz)

    # -- Tracheophytes -- # 

# Un-weighted breaks 
UW_Break_Tracheo_A <- cut(c(Tracheo_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
UW_Break_Tracheo_B <- cut(c(Tracheo_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 

UW_Break_Tracheo_A_Nb <- cut(c(Tracheo_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 
UW_Break_Tracheo_B_Nb <- cut(c(Tracheo_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 

# Weighted breaks
W_Break_Tracheo_A <- cut(c(Tracheo_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
W_Break_Tracheo_B <- cut(c(Tracheo_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 

W_Break_Tracheo_A_Nb <- cut_number(c(Tracheo_Beta$PlotA_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 
W_Break_Tracheo_B_Nb <- cut_number(c(Tracheo_Beta$PlotB_z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 

# Add them to the initial data 
Tracheo_Beta <- mutate(Tracheo_Beta,"Taxa" = "Tracheophytes", "UW_BreakA" = UW_Break_Tracheo_A,"UW_BreakB" = UW_Break_Tracheo_B, "W_BreakA" = W_Break_Tracheo_A,"W_BreakB" = W_Break_Tracheo_B, "UW_StripeA" = UW_Break_Tracheo_A_Nb,"UW_StripeB" = UW_Break_Tracheo_B_Nb, "W_StripeA" = W_Break_Tracheo_A_Nb,"W_StripeB" = W_Break_Tracheo_B_Nb, .after = delta_xyz)

    # ----- #

      # -- Stripe distances -- # 

# Add the stripe distances 
Liverworts_Beta <- Liverworts_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)
Mosses_Beta <- Mosses_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)
Bryophytes_Beta <- Bryophytes_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)
Tracheo_Beta <- Tracheo_Beta %>% mutate(UW_Stripe_dist = abs(UW_StripeA - UW_StripeB), W_Stripe_dist = abs(W_StripeA - W_StripeB), .after = W_StripeB)

# Create a list of these 4 datasets
Beta_List <- list(Mosses_Beta,Liverworts_Beta,Bryophytes_Beta,Tracheo_Beta) %>%
  setNames(c("Mosses","Liverworts","Bryophytes","Tracheophytes"))

# Clean the environnment
rm(list=ls(pattern="_Break_"))

# Message
cat(rule(left = "- Altitudinal stripes added - ", line_col = "white", line = " ", col = "green"))


#------------------------------------#
##### ENVIRONMENTAL DATA MERGING #####
#------------------------------------#











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