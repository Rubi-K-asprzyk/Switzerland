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
       "ggtree"
)

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
color_scheme <-  c("Bryophytes" = "orange1", "Liverworts" = "yellowgreen", "Mosses" = "darkgreen", "Tracheophytes" = "orchid1")

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

# --- Load the data --- #

# Load the environmental mosses data
Mosses_Data <- read.csv(file = "Data/Bryophytes/Mosses_Total.csv")

# Load the environmental liverworts data
Liverworts_Data <- read.csv(file = "Data/Bryophytes/Liverworts_Total.csv")

# Create a global "Bryophyte Data" by merging the Moss and Liverworts data
Bryophytes_Data <- merge(Mosses_Data,Liverworts_Data)

# Load the environmental tracheophyte data
Tracheo_Data <- read.csv(file = "Data/Angiosperm/Tracheo_Total.csv")

# Load the species names
Sp_Names <- read.csv(file = "Utilities/sp_names.csv")
# Split them between the different taxons
  Liver_names <- Sp_Names[which(Sp_Names$Taxa == "Liverworts"),1]
  Moss_names <- Sp_Names[which(Sp_Names$Taxa == "Mosses"),1]
  Bryo_names <- Sp_Names[which(Sp_Names$Taxa != "Tracheophytes"),1]
  Tracheo_names <- Sp_Names[which(Sp_Names$Taxa == "Tracheophytes"),1]
  
# Load the phylotrees
Liver_Tree <- read.tree("Data/Bryophytes/timetree50mod-liverworts.nwk")
Mosses_Tree <- read.tree("Data/Bryophytes/timetree50mod-mosses.nwk")

# Tree tips modifications 
ggtree(Liver_Tree) + geom_tiplab() + ggplot2::xlim(0, 0.5)
ggtree(Mosses_Tree) + geom_tiplab()

# Message
cat(rule(left = "- Data loaded - ", line_col = "white", line = " ", col = "green"))

#------------------------------------------------#
##### STEP 0: Creation of a global dataframe #####
#------------------------------------------------#

##### ----- 0.A.1: Taxonomic alpha diversity -------------------- #####

# Initialize the metric we are working with.
Metric <- c("Sp_Richness","Sp_Simpson","Sp_Shannon")
Metric <- c("Sp_Richness")
Metric_B <- c("beta.sim","beta.sne","beta.sor","beta.jtu","beta.jne","beta.jac")

# Computation of the taxonomic alpha diversity across the altitudinal gradient.
  # Creation of a function to compute the species richness, the Simpson and Shannon diversity. 
  # Take as parameter : 
    # - Data : A dataframe containing the Occurrence data 
    # - Sp_names : A character vector of species names
    
Alpha_Taxo <- function(Data,Sp_names) {
  
  # Species Richness
  if("Sp_Richness" %in% Metric) {
    # Computation
    Value <- apply(dplyr::select(Data,all_of(Sp_names)),1,sum)
    # Add the results to the initial dataframe
    Data <- mutate(Data,"Sp_Richness" = Value,.after = z)
  }
  
  # Simpson Index
  if("Sp_Simpson" %in% Metric) {
    # Computation
    Value <- vegan::diversity(dplyr::select(Data,all_of(Sp_names)),index = "simpson")
    # Add the results to the initial dataframe
    Data <- mutate(Data,"Sp_Simpson" = Value,.after = z)
  }
  
  # Simpson Index
  if("Sp_Shannon" %in% Metric) {
    # Computation
    Value <- vegan::diversity(dplyr::select(Data,all_of(Sp_names)),index = "shannon")
    # Add the results to the initial dataframe
    Data <- mutate(Data,"Sp_Simpson" = Value,.after = z)
  }
  
  # Return the dataframe containing the results
  return(Data)
  
}

# Mosses
Mosses_Data <- Alpha_Taxo(Data = Mosses_Data,Moss_names)
# Liverworts
Liverworts_Data <- Alpha_Taxo(Data = Liverworts_Data,Liver_names)
# Bryophytes
Bryophytes_Data <- Alpha_Taxo(Data = Bryophytes_Data,Bryo_names)
# Tracheophytes
Tracheo_Data <- Alpha_Taxo(Data = Tracheo_Data,Tracheo_names)

# Message
cat(rule(left = "- Taxonomic alpha diversity added - ", line_col = "white", line = " ", col = "green"))

##### ----- 0.A.2: Taxonomic Beta diversity -------------------- #####

# Creation of a function to compute the beta-diversity indices.
# Take as parameter : 
# - Data : A dataframe containing the Occurrence data.
# - Type : The type of data we are working on : "B" for bryophytes or "T" for tracheophytes. 
# - Sp_names : A character vector of species names

Beta_Taxo <- function(Data,Sp_names,Type = "B") {
  
    # ----- #
  
  # --- Creation of the dataframe --- # 
  
  Beta_Occ <-
    if(Type == "B"){ 
      Data %>% 
        column_to_rownames(var = "Site_VDP") %>% # Change the row_names to be the plot number/name.
        dplyr::select(all_of(Sp_names)) # Remove the meta_data columns.
    } else {
      rownames(Data) <- NULL
      Data %>% 
        column_to_rownames(var = "X") %>% 
        dplyr::select(all_of(Sp_names))
    }
  
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

# Mosses
Mosses_Beta <- Beta_Taxo(Data = Mosses_Data,Moss_names)
# Liverworts
Liverworts_Beta <- Beta_Taxo(Data = Liverworts_Data,Liver_names)
# Bryophytes
Bryophytes_Beta <- Beta_Taxo(Data = Bryophytes_Data,Bryo_names)
# Tracheophytes
Tracheo_Beta <- Beta_Taxo(Data = Tracheo_Data,Tracheo_names,Type = "T")

# Message
cat(rule(left = "- Taxonomic beta diversity added - ", line_col = "white", line = " ", col = "green"))

##### ----- 0.B: Altitudinal stripes -------------------- #####

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

    # -- Tracheophytes -- # 

# Un-weighted breaks 
UW_Break_Tracheo <- cut(c(Tracheo_Data$z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
UW_Break_Tracheo_Nb <- cut(c(Tracheo_Data$z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4)

# Weighted breaks
W_Break_Tracheo <- cut_number(c(Tracheo_Data$z,Added_Range), breaks, dig.lab = 4) %>% head(-4) 
W_Break_Tracheo_Nb <- cut_number(c(Tracheo_Data$z,Added_Range), breaks, dig.lab = 4, labels = F) %>% head(-4) 

# Add it to the initial data
Tracheo_Data <- mutate(Tracheo_Data,"Taxa" = "Tracheophytes","UW_Break" = UW_Break_Tracheo, "W_Break" = W_Break_Tracheo,"UW_Stripe" = UW_Break_Tracheo_Nb, "W_Stripe" = W_Break_Tracheo_Nb,.after = z)

# /!\ WARNING /!\
# Three plots of the TracheoData are removed because they lack some data.
Tracheo_Data <- Tracheo_Data[-c(179,182,1057),]

# Create a list of these 4 datasets
Data_List <- list(Mosses_Data,Liverworts_Data,Bryophytes_Data,Tracheo_Data) %>%
  setNames(c("Mosses","Liverworts","Bryophytes","Tracheophytes"))

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

##### ----- 0.C: PhyloTree tips -------------------- #####

# PhyloTree tips are not congruent with species names in the community, make them congruent.
# Liver_names

#--------------------------------------#
##### STEP 1 : TAXONOMIC DIVERSITY #####
#--------------------------------------#

##### ______________________________________________________________________________________________________________________________________ #####

##### ------------ 1.A: TAXONOMIC ALPHA -------------------- #####

##### ______________________________________________________________________________________________________________________________________ #####

##### ------------ 1.A.a : Global dataframe(s)  ------------ #####

# Create a global dataframe containing the results for all taxa. 
Alpha_Data <- Data_List %>%
  # Bind the all the wanted taxa
  reduce(full_join) %>%
  # Keep the wanted columns
  # select(c("x", "y", "z", "Taxa", "UW_Break", "W_Break", "UW_Stripe", "W_Stripe", Metric)) %>%
  # Combine all the metrics into one column
  gather(key = "Metric", value = "Value", !!Metric) %>%
  # Relocate the Metric and Value column
  relocate(any_of(c("Metric","Value")),.after  = Taxa) %>%
  # Transform some columns into factors
  mutate(across(c("UW_Break", "W_Break", "UW_Stripe", "W_Stripe"), as_factor)) %>%
  # Group the data
  group_by(Metric, Taxa, .data[[Break_S]]) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  # This will be used in the future tests, to remove groups with null variance 
  mutate(Sum = sum(Value), Count=n(), Var_Value = var(Value), Mean_Value = mean(Value))

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

  # METRICS SUMMARY STATISTICS --- #

# Message
cat(rule(left = "- Summary statistics - ", line_col = "white", line = "_", col = "grey"))
    
# Compute and create the summary table of the metric(s) splitted between stripes
Summary_splitted <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa, .data[[Break_S]]) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = Break_S,                           # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 color = Break_S,                       # Color by stripes
                 facet.by = c("Metric", "Taxa"),        # Split by Metric and Taxa
                 palette  = Pal_col,                    # Color Palette
                 ggtheme = arrange_theme() +            # Theme
                   theme(legend.position = "none"))
    
    
# Compute and create the summary table of the metric(s) unsplitted between stripes
Summary_unsplitted <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "Metric",                           # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 facet.by = "Taxa",        # Split by Metric and Taxa
                 ggtheme = arrange_theme() +            # Theme
                   theme(legend.position = "none"))

    # ----- # 

  # GLOBAL PLOTTING OF TAXONOMIC ALPHA METRICS --- #

    # ----- # 
    
  # - Distribution of Metric(s) ~ Altitudinal stripes - # 
    
Plot1_Glob <- Alpha_Data %>% 
  # Group the data
  group_by(Metric, Taxa) %>%
  # Aes
  ggplot(aes(x = z, y = Value, color = .data[[Break]])) +
  # Palette 
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_point() +
  # Split between Taxa and Metric
  facet_grid(Metric ~ Taxa, scales = "free_y") +
  # Plot the connected scatter plot of the mean values for each altitudinal stripe.
  new_scale_color() + # Define scales before initiating a new one
  # Use the new scale
  geom_point(data = Summ_Alpha_Data, aes(y = Mean_Value), size = 1) +
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
   
    # ----- #
    
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
    
##### ----- 1.A.2.b : Individual Plots ------------------- #####

# Create a loop to draw and save a mix-up of the previous global plots created splitted by taxa and metric. 
# Take as parameter : 
  # - Taxa : A character vector of the name of the Taxa we are currently working on..
  # - Metric : A character vector of name of the Metric we are currently working on. 
   
# Create these two character vectors 
Metric <- unique(Alpha_Data$Metric)
Taxa <- unique(Alpha_Data$Taxa)

    # ------------------ #

# Use a double loop of foreach
Alpha_Taxo_Ind <- 
  foreach(j = 1:length(Taxa)) %:%
    foreach(i = 1:length(Metric)) %do% {
      
    # ----- # 
      
  # METRICS SUMMARY STATISTICS --- #
      
# Compute and create the summary table of the metric(s) splitted between stripes
Summary_splitted_ind <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i])) %>%
  # Group the data
  group_by(.data[[Break_S]]) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = Break,                             # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 color = Break_S,                       # Color by stripes
                 palette  = Pal_col,                    # Color Palette
                 ggtheme = arrange_theme() +            # Theme
                  theme(legend.position = "none"))
      
      
# Compute and create the summary table of the metric(s) unsplitted between stripes
Summary_unsplitted_ind <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i])) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "Metric",                          # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 ggtheme = arrange_theme() +            # Theme
                 theme(legend.position = "none"))
      
    # ----- # 
      
  # INDIVIDUAL PLOTTING OF TAXONOMIC ALPHA METRICS --- #
      
    # ----- # 
      
# - Distribution of Metric(s) ~ Altitudinal stripes - # 

# Filter Summ_Alpha_Data for the wanted metric. 
Ind_Summ_Alpha <- Summ_Alpha_Data %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i]))

# Create the density plot
Plot1_Ind <- Alpha_Data %>% 
  # Group the data
  group_by(Metric, Taxa) %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i])) %>%
  # Aes
  ggplot(aes(x = z, y = Value, color = .data[[Break]])) +
  # Palette 
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_point() +
  # Plot the connected scatter plot of the mean values for each altitudinal stripe.
  new_scale_color() + # Define scales before initiating a new one
  # Use the new scale
  geom_point(data = Ind_Summ_Alpha, aes(y = Mean_Value), size = 1) +
  geom_line(aes(y = Mean_Value),linewidth = 1, alpha = 0.7, group = 1) +
  # Labels
  xlab("Altitude") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("ScatterPlot of ",Metric[i]," values ~ Altitude"),
    subtitle = paste0(
      "Black dotted line represent the mean for each metric of the altitudinal stripes.\nAltitudinal stripes are",{ifelse(Break == "W_Break",paste(" weighted by plot number i.e represent an equal number of plots."),paste(" unweighted by plot number i.e represent an equal geographical distance."))}))

    # ----- #

  # - Density plots of the Metric(s) values - #

Plot2_Ind <- Alpha_Data %>% 
  # Group the data
  group_by(Metric, Taxa) %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i])) %>%
  # Aes
  ggplot(aes(x = Value)) +
  # Plot
  geom_density(fill = "grey",alpha = 0.5) +
  # Labels
  labs(
    title = paste0("DensityPlot of ",Metric[i]))

    # ----- #

  # - BoxPlot of the Metric(s) values - #

# Filter KW_S for the wanted metric. 
Ind_WT_S_Adj  <- WT_S %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i]) & Stripe_Distance == 1) %>%
  # Add the values of y.position 
  left_join(y.position, by=c("Metric","Taxa"))

# Create the boxplot
Plot3_Ind <- Alpha_Data %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i])) %>%
  # Aes
  ggplot(aes(x = .data[[Break_S]], y = Value, color = .data[[Break]], group = .data[[Break]])) +
  # Color Palette    
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_boxplot() +
  # Add the significativity labels (if some are significant)
  {if (sum(str_detect(Ind_WT_S_Adj$p.adj.signif,'ns')) != length(Ind_WT_S_Adj$p.adj.signif)) {
    # Add the significativity
    stat_pvalue_manual(Ind_WT_S_Adj, 
                     label = "p.adj.signif",
                     hide.ns = T)
    }} + 
  # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("BoxPlot of ",Metric[i]," ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcoxon tests were realized between adjacent stripes and significant results are displayed.",
                      "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

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
Plot4_Ind  <- Group_Signif %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Filter the data
  filter(str_detect(Taxa, !!Taxa[j]) & str_detect(Metric, !!Metric[i])) %>%
  # Aes
  ggplot(aes(x = Stripe_Distance, y = Signif_Percent, group = Significance, color = Significance)) +
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
    title = paste0("Percentage of Significant Wilcoxon tests of ",Metric[i]," ~ Stripe distance"),
    subtitle = "Stripe distance is an articifial unit to describe how far the stripes are from each other.\nThe higher the stripe distance, the lower the number of possible pairwise comparisons.\nAll significance levels are combined into the two modalities."
  )  

    # ----- #

  # - Arrange the plots all_together - #

    # PLOT1_ind = Distribution of Metric(s) ~ Altitudinal stripes 
    # PLOT2_ind = Density plots of the Metric(s) values 
    # PLOT3_ind = BoxPlot of the Metric(s) values
    # PLOT4_ind = ScatterPlot of the Kruskall_Wallis tests results
    # summary_splitted
    # summary_unsplitted


# Arrange the three graphs altogether with the statistics below.
Plot_Total <- 
  ggarrange(
    # Arrange the three main plots
    ggarrange(Plot1_Ind,Plot3_Ind,Plot4_Ind, ncol = 3, legend = "right", common.legend = T), 
    # Arrange the summary stats and plot2
    ggarrange(Summary_splitted_ind,Summary_unsplitted_ind,Plot2_Ind, ncol = 3, widths = c(3,1,1),legend = "none"), 
    # Set the global parameters
    nrow = 2, heights = c(3,1)) %>% 
  as.ggplot() +
  arrange_theme() +
  labs(
    title = paste0("Alpha Taxonomic diversity of ",Taxa[j])
  )

# Return this global plot
return(Plot_Total)


} # END OF THE INDIVIDUAL PLOTTING

    # ------------------ #
  
##### ----- 1.A.2.c : Final Plot ------------------- #####

  # Combination of the global plots and the individual plot inside one pdf output.

# Append and flatten all the plots inside one list.
pl <- append(
  # List of global plots
  list(Summary_splitted,Summary_unsplitted,Plot1_Glob,Plot2_Glob,Plot3_Glob,Plot4_Glob,Plot5_Glob),
  # Flattened list of plots mix-up by metrics / taxa
  flatten(Alpha_Taxo_Ind))

# Transform them into grobs
pl <- lapply(pl, as.grob)
# Create the pdf
ml <- marrangeGrob(pl,nrow=1,ncol = 1)
# Save a global plot
ggsave(filename = paste0("STEP1_",breaks,Break,"_Richness.pdf"),
       plot = ml,
       device = "pdf",
       path = paste0(getwd(),"/Plots"),
       width = 35,
       height = 25,
       units = "cm")


# -------------------------------------------------------- #
# -------------------------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

##### ----- 1.B: TAXONOMIC BETA ------------------- #####

##### ______________________________________________________________________________________________________________________________________ #####


# Breaks have to be changed from the intervals to the simple number
if(Break == "UW_Break") {Break_S <- "UW_Stripe_dist" ; BreakIntra <- "UW_StripeA"
} else if (Break == "W_Break") {Break_S <- "W_Stripe_dist" ; BreakIntra <- "W_StripeA"}

# Create a global dataframe containing the results for all taxas. 
Alpha_Beta <- Beta_List %>% 
  reduce(full_join, by = c("Sample", "PlotA", "PlotB", "PlotA_x", "PlotA_y", "PlotA_z", "PlotB_x", "PlotB_y", "PlotB_z", "delta_z", "delta_xyz", "Taxa", "UW_BreakA", "UW_BreakB","W_BreakA", "W_BreakB", "UW_StripeA", "UW_StripeB", "W_StripeA", "W_StripeB", "UW_Stripe_dist", "W_Stripe_dist", "beta.sim", "beta.sne", "beta.sor", "beta.jtu", "beta.jne","beta.jac")) %>%
  gather(key = "Metric", value = "Value", !!Metric_B) %>%
  mutate(across(c("UW_StripeA","UW_StripeB","W_StripeA","W_StripeB","UW_Stripe_dist","W_Stripe_dist"), as_factor)) %>%
  # Remove the NaN values
  filter(!is.nan(Value)) %>%
  # Relocate Metric, Taxa and Values
  relocate(any_of(c("Taxa","Metric","Value")),.before  = Sample)

##### ----- 1.B.1: Statistics ------------------- #####

# First step, verify the normality of the distribution of the values for each altitudinal stripe

# We'll use the of test of Kolmogorov-Smirnoff because we have a lot of pairwise values per group. 
# H0: Data are normally distributed: if p-value < 0.05, we reject that hypothesis.

  #### . Shapiro-Wilk (Metric + Taxa) . ####

# Are the metrics values UNSPLITTED between stripes distributed normally ? (Grouped by Metric and Taxa).
SW_US <- Alpha_Beta %>%
  # Grouped the data
  group_by(Metric,Taxa) %>%
  # Execute the test
  do(tidy(ks.test(.$Value, y = "pnorm", mean = mean(.$Value), sd=sd(.$Value)))) %>%
  # Remove unwanted columns
  select(-c(method,alternative))

  # RESULT:
cat(paste0("Normality tests grouped by Metric and Taxa: ",sum(SW_US$p.value > 0.05),"/",nrow(SW_US)," tests passed."))
    
  # ---- #

  #### . Shapiro-Wilk (Metric + Taxa) . ####

# Are the metrics values SPLITTED between stripes distributed normally ? (Grouped by Metric, Taxa and Stripe Distance).
SW_S <- Alpha_Beta %>%
  # Grouped the data
  group_by(Metric, Taxa, UW_Stripe_dist) %>%
  # Remove the NaN values
  filter(!is.nan(Value)) %>%
  # Execute the test
  do(tidy(ks.test(.$Value, y = "pnorm", mean = mean(.$Value), sd=sd(.$Value)))) %>%
  # Remove unwanted columns
  select(-c(method,alternative))

  # RESULT:
cat(paste0("Normality tests grouped by Metric, Taxa and stripe distances: ",sum(SW_S$p.value > 0.05),"/",nrow(SW_S)," tests passed."))

# Merge the two datasets
SW_Result <- full_join(SW_US,SW_S)

##### ----- 1.B.b.2: Equality of mean tests. ----- #####

# - Hypothesis of normality are not met, therefore, we'll use kruskall-wallis tests to determine if means of the metrics are different - # 

# A/ Between stripes for the same taxon, with computation intra-stripes.
# B/ Between stripes for the same taxon, with increasing distances between the plots.
# C/ Between taxon for the same stripes (Intra and Inter)

# We will use the simple factor of breaks to be able to compute a "stripe" distance. 

# --- Trim the global data --- #

  # -- Remove the groups where the variance is equal to 0.

# Find the count of plots in each stripe
Alpha_Beta <- Alpha_Beta %>%
  # Group the data
  group_by(Metric, Taxa, .data[[Break_S]]) %>%
  # Remove the NaN values
  filter(!is.nan(Value)) %>%
  # Compute the Sum of the Metric and the number of plots for each group.
  mutate(Sum = sum(Value), Count=n(), Var_Value = var(Value), Mean_Value = mean(Value)) %>%
  # Remove the lines with 0 variance.
  filter(Var_Value != 0)

  # -- Create a dataset of mean values INTRA.

# Find the count of plots in each stripe
Summ_Alpha_Beta_Intra <- Alpha_Beta %>%
  # Filter the data to keep only the intra-stripes comparisons 
  filter(.data[[Break_S]] == 0) %>%
  # Group the data
  group_by(Metric,Taxa,.data[[BreakIntra]]) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  summarise(Mean_Value = mean(Value), delta_z = mean(delta_z), Sd_Value = sd(Value))

# --- Compute the significance values with Kruskall-Wallis tests --- #

# Is there a significant difference between the groups ? 
# If p-value <0.05, there is a significant difference between the groups. 

  #### . Kruskall-Wallis (Metric + Taxa) . ####

# Message
cat(rule(left = "- Kruskall-Wallis Unsplitted - ", line_col = "white", line = "_", col = "grey"))

  # KRUSKALL-WALLIS TOTAL --- #

# Is there a significant difference between metrics values UNSPLITTED (Groupped by Metric and Taxa) ? 

KW_US <- Alpha_Beta %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Compute the kruskall-test
  do(tidy(kruskal.test(x = .$Value, g = .data[[Break_S]]))) %>%
  # Remove the method
  select(-c(method,parameter))

  # RESULT:
cat(paste0("KRUSKALL-WALLIS TOTAL: Significant differences grouped by Metric and Taxa: ",sum(KW_US$p.value < 0.05),"/",nrow(KW_US)," tests passed."))

    # ----- # 

  # KRUSKALL-WALLIS EFFECT-SIZE --- #

# Effect size tells you how meaningful the relationship between variables or the difference between groups is. It indicates the practical significance of a research outcome.
# A large effect size means that a research finding has practical significance, while a small effect size indicates limited practical applications.
# The value multiplied by 100 represent the percentage of variance of the dependent variable explained by the independent one. 

  #### . Kruskall-Wallis effect size (Metric + Taxa) . ####

# Message
cat(rule(left = "- Kruskall-Wallis Effect-Size - ", line_col = "white", line = "_", col = "grey"))

KW_ES <- Alpha_Beta %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Compute the kruskall-test
  kruskal_effsize(formula = as.formula(paste("Value","~",Break_S))) %>%
  # Remove the method
  select(-c(method))

    # ----- # 

  # KRUSKALL-WALLIS INTRA --- #

  #### . Kruskall-Wallis intra stripes (Metric + Taxa) . ####

# Message
cat(rule(left = "- Kruskall-Wallis Intra -  ", line_col = "white", line = "_", col = "grey"))

# The goal here is to determine if there is differences between the beta-comparisons intra-stripes. 
# So only stripe distances between plots compared pairwise equal to zero. And we have that for all the different stripes

KW_Intra <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Filter the data to keep only the intra-stripes comparisons 
  filter(.data[[Break_S]] == 0) %>%
  # Perform ANOVA tests
  do(tidy(kruskal.test(x = .$Value, g = .data[[BreakIntra]]))) %>%
  # Remove the method
  select(-c(method,parameter))

  # RESULT:
cat(paste0("KRUSKALL-WALLIS INTRA: Significant differences grouped by Metric and Taxa: ",sum(KW_Intra$p.value < 0.05),"/",nrow(KW_Intra)," tests passed."))

    # ----- # 

  # KRUSKALL-WALLIS INTER --- #

#### . Kruskall-Wallis inter stripes (Metric + Taxa + Stripe distances) . ####

# Message
cat(rule(left = "- Kruskall-Wallis Inter -  ", line_col = "white", line = "_", col = "grey"))

# The goal here is to determine if there is differences between the beta-comparisons inter-stripes ~ Stripe distance 

KW_Inter <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Do the test
  do(tidy(kruskal.test(x = .$Value, g = .data[[BreakIntra]]))) %>%
  # Remove the method
  select(-c(method,parameter))

  # RESULT:
cat(paste0("KRUSKALL-WALLIS INTER: Significant differences grouped by Metric and Taxa: ",sum(KW_Inter$p.value < 0.05),"/",nrow(KW_Inter)," tests passed."))

    # ----- # 

  # WILCOX TESTS INTRA --- #

#### . Wilcoxon tests intra stripes (Metric + Taxa) . ####

# Message
cat(rule(left = "- Wilcox tests Intra-Stripes -  ", line_col = "white", line = "_", col = "grey"))

# Comparison of the mean values of the metrics results from intra-pairwise plots computations. 
# i.e: Metrics are computed between plots from the same stripe. The results for each stripes are compared to determine if they are significantly different from each other. 

WT_Intra <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Filter the data to keep only the intra-stripes comparisons 
  filter(.data[[Break_S]] == 0) %>%
  # There is some trouble with the test, we need to remove some groups giving problems
  # Group the data with the stripe values 
  group_by(Metric,Taxa,.data[[BreakIntra]]) %>%
  # Recompute the mean and the variance of the groups
  mutate(Sum = sum(Value), Count=n(), Var_Value = var(Value), Mean_Value = mean(Value)) %>%
  # remove the groups with NullVariance
  filter(Var_Value != 0) %>%
  # Regroup again with the good groups 
  ungroup() %>%
  group_by(Metric,Taxa) %>%
  # pairwise comparisons
  wilcox_test(formula = as.formula(paste("Value","~",BreakIntra)), p.adjust.method = "bonferroni") %>%
  # Transform "group1" and "group2" into numeric
  mutate_at(c("group1", "group2"), as.numeric ) %>%
  # Add a column that is the "stripe distance" between the distribution compared
  mutate(Stripe_Distance = abs(group1 - group2))

    # ----- # 

  # WILCOX TESTS INTRA STRIPES / ONLY COMPARISONS BETWEEN ADJACENT STRIPES --- #

#### . Wilcoxon tests intra stripes / Adjacent Stripes (Metric + Taxa) . ####

# Message
cat(rule(left = "- Wilcox tests Intra-Stripes / comparisons between adjacent stripes - ", line_col = "white", line = "_", col = "grey"))

# We need a column "y.position" for the plotting of the significance brackets
y.position <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Add the y.position with an increase of X%. 
  summarise(y.position = max(Value) * 1.2) 

# Filter the precedent data_frame to only keep the values of stripe distance == 1
WT_Intra_Adj <- WT_Intra %>%
  # Filter the data 
  filter(Stripe_Distance == 1) %>%
  # Add the values of y.position 
  left_join(y.position, by=c("Metric","Taxa"))

    # ----- # 

  # WILCOX TESTS INTER --- #

#### . Wilcoxon tests inter stripes (Metric + Taxa) . ####

# Message
cat(rule(left = "- Wilcox tests Inter-Stripes -  ", line_col = "white", line = "_", col = "grey"))

# Comparison of the mean values of the metrics results from inter-pairwise plots computations. 
# i.e: Metrics are computed between plots from different stripes. The results for each stripe distances are compared to determine if they are significantly different from each other. 

WT_Inter <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa,.data[[Break_S]]) %>%
  # Recompute the mean and the variance of the groups
  mutate(Sum = sum(Value), Count=n(), Var_Value = var(Value), Mean_Value = mean(Value)) %>%
  # remove the groups with NullVariance
  filter(Var_Value != 0) %>%
  # Regroup again with the good groups 
  ungroup() %>%
  group_by(Metric,Taxa) %>%
  # We compare all the stripes distances to a ref_group, that is the Intra Group: Stripe_Dist = "O"
  wilcox_test(formula = as.formula(paste("Value","~",Break_S)), p.adjust.method = "bonferroni", ref.group = "0") %>%
  # Add the values of y.position 
  left_join(y.position, by=c("Metric","Taxa"))


##### ----- 1.B.2.a : Global Plots ------------------- #####

    # ----- # 

  # METRICS SUMMARY STATISTICS --- #

# Message
cat(rule(left = "- Summary statistics - ", line_col = "white", line = "_", col = "grey"))

# Compute and create the summary table of the metric(s) splitted between stripes distances between the two plots of the pairwise comparisons.
Summary_splitted <- Alpha_Beta %>%
  # Group the data
  group_by(Metric, Taxa, .data[[Break_S]]) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = Break_S,                           # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 color = Break_S,                       # Color by stripes
                 facet.by = c("Metric", "Taxa"),        # Split by Metric and Taxa
                 palette  = Pal_col,                    # Color Palette
                 ggtheme = arrange_theme() +            # Theme
                   theme(legend.position = "none")) +
  # Labels
  xlab("Stripe distance") +
  labs(
    title = ("Inter-Stripes metrics values ~ Altitudinal stripes number"),
    subtitle = paste0("Pairwise Beta-metrics are computed between plots from different stripe distances.\nAltitudinal stripes are",{ifelse(Break == "W_Break",paste(" weighted by plot number i.e represent an equal number of plots."),paste(" unweighted by plot number i.e represent an equal geographical distance."))}))


  # ----- #

# Compute and create the summary table of the metric(s) unsplitted between stripes distances of the pairwise comparisons.
Summary_unsplitted <- Alpha_Beta %>%
  # Group the data
  group_by(Metric, Taxa) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "Metric",                           # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 facet.by = "Taxa",        # Split by Metric and Taxa
                 ggtheme = arrange_theme() +            # Theme
                   theme(legend.position = "none"))

  # ----- #

# Compute and create the summary table of the metric(s) splitted between stripes. Plots from the pairwise comparisons being from the same stripe. (INTRA)
Summary_splitted_Intra <- Alpha_Beta %>%
  # Filter the data to keep only the intra-stripes comparisons 
  filter(.data[[Break_S]] == 0) %>%
  # Group the data
  group_by(Metric,Taxa,.data[[BreakIntra]]) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = BreakIntra,                           # Split by stripes
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 2,                            # Number of digits 
                 size = 3,                              # Size of the text
                 color = BreakIntra,                       # Color by stripes
                 facet.by = c("Metric", "Taxa"),        # Split by Metric and Taxa
                 palette  = Pal_col,                    # Color Palette
                 ggtheme = arrange_theme() +            # Theme
                   theme(legend.position = "none")) + 
  xlab("Stripe Number") +
  labs(
    title = ("Intra-Stripes metrics values ~ Altitudinal stripes number"),
    subtitle = paste0("Pairwise Beta-metrics are computed between plots from the same altitudinal stripe.\nAltitudinal stripes are",{ifelse(Break == "W_Break",paste(" weighted by plot number i.e represent an equal number of plots."),paste(" unweighted by plot number i.e represent an equal geographical distance."))}))

    # ----- # 

  # GLOBAL PLOTTING OF TAXONOMIC BETA METRICS --- #

    # ----- # 

# - Distribution of Metric(s) ~ Altitudinal stripes ~ Intra Pairwise Plots - # 

  # BoxPlots and significativity levels of the Intra-Stripes metrics values, splitted between metrics and taxas. 

Plot1_a_Glob <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Filter the data to keep only the intra-stripes comparisons 
  filter(.data[[Break_S]] == 0) %>%
  # Aes
  ggplot(aes(x = .data[[BreakIntra]], y = Value, color = .data[[BreakIntra]])) +
  # Palette 
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_boxplot() +
  # Split between Taxa and Metric
  facet_grid(Metric ~ Taxa, scales = "free_y") +
  # Add the significativity labels
  stat_pvalue_manual(WT_Intra_Adj, 
                     label = "p.adj.signif",
                     hide.ns = T) +
  # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  guides(color = "none") +
  labs(
    color = "Stripes",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcoxon tests were realized between adjacent stripes and significant results are displayed.",
                      "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

    # ----- #

  # Connected scatterplots and significativity levels of the Intra-Stripes metrics values, splitted between metrics

Plot1_b_Glob <- Summ_Alpha_Beta_Intra %>%
  # Aes
  ggplot(aes(x = .data[[BreakIntra]], y = Mean_Value, color = Taxa, group = Taxa)) +
  # Palette 
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_line() + 
  geom_point() +
  # Split between Taxa and Metric
  facet_grid(Metric ~ ., scales = "free_y") +
  # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Taxa",
    title = paste0("Connected scatterplot of metric(s) mean values")
  )

    # ----- #

# - Combination of the two previous plots of the distribution of Metric(s) ~ Altitudinal stripes ~ Intra Pairwise Plots - #   

Plot1_Glob <- 
  # Arrange the two plots
  ggarrange(Plot1_a_Glob,Plot1_b_Glob, ncol = 2, legend = "right", common.legend = T, widths = c(3,1)) %>% 
  as.ggplot() +
  arrange_theme()

    # ----- #
    # ----- # 

  # - Distribution of Metric(s) ~ Altitudinal stripes ~ Inter Pairwise Plots stripe distance - # 

# BoxPlots and significativity levels of the Inter Pairwise Plots stripe distance metrics values, splitted between metrics and taxas. 

Plot2_a_Glob <- Alpha_Beta %>%
  # Group the data
  group_by(Metric,Taxa) %>%
  # Aes
  ggplot(aes(x = .data[[Break_S]], y = Value, color = .data[[Break_S]])) +
  # Palette 
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_boxplot() +
  # Split between Taxa and Metric
  facet_grid(Metric ~ Taxa, scales = "free_y") +
  # # Add the significativity labels
  # stat_pvalue_manual(WT_Inter, 
  #                    label = "p.adj.signif",
  #                    hide.ns = T,
  #                    step.increase = 0.1,
  #                    step.group.by = c("Metric","Taxa")) +
  # Labels
  xlab("Altitudinal stripe distance") +
  ylab("Metric values") +
  guides(color = "none") +
  labs(
    color = "Stripes",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripe distance between plots"),
    subtitle = paste0("Wilcoxon tests were realized between intra-plot metrics (0) and all possible stripe distances. Significant results are displayed.",
                      "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  )

    # ----- #

  # Connected scatterplots and significativity levels of the Intra-Stripes metrics values, splitted between metrics

Plot2_b_Glob <- Alpha_Beta %>%
  # Aes
  ggplot(aes(x = .data[[Break_S]], y = Mean_Value, color = Taxa, group = Taxa)) +
  # Palette 
  scale_color_brewer(palette = Pal_col) +
  # Plot all the values
  geom_line() + 
  geom_point() +
  # Split between Taxa and Metric
  facet_grid(Metric ~ ., scales = "free_y") +
  # Labels
  xlab("Stripe distance") +
  ylab("Metric values") +
  labs(
    color = "Taxa",
    title = paste0("Connected scatterplot of metric(s) mean values")
  )

    # ----- #

  # - Combination of the two previous plots of the distribution of Metric(s) ~ Altitudinal stripes ~ Intra Pairwise Plots - #   

Plot2_Glob <- 
  # Arrange the two plots
  ggarrange(Plot2_a_Glob,Plot2_b_Glob, ncol = 2, legend = "right", common.legend = T, widths = c(3,1)) %>% 
  as.ggplot() +
  arrange_theme()

  # ----- #

# - ScatterPlot of the Kruskall_Wallis tests results - #

# Create a dataset of significance for the next plots.
Group_Signif <- KW_US %>%
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
Plot4_Glob  <- Group_Signif %>%
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
    title = paste0("Percentage of Significant Kruskall_Wallis tests of metric values ~ Stripe distance"),
    subtitle = "Stripe distance is an articifial unit to describe how far the stripes are from each other.\nThe higher the stripe distance, the lower the number of possible pairwise comparisons.\nAll significance levels are combined into the two modalities."
  )




# --- Plot the data --- #

# Create a function to create, draw and save a global ggplot object for the taxonomic alpha diversity of a specific taxonomic group
# Take as parameter : 
# - Data : The taxoBeta 
# - Taxa : The taxa we are working with as a character


Beta_Taxo_Plot <- function(Data,Taxa){

# Different possibilities:
  # The best but long and computational heavy.
   # geom_pointdensity() + 
   # scale_color_viridis_c() +
  
  # Much faster
   # stat_binhex() +
   # scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 8000)) +
  
p1_sor <- Data %>% 
  ggplot(aes(x = delta_z, y = beta.sor)) +
  stat_binhex(na.rm = T,
              bins = 50) +
  scale_fill_viridis_c() +
  xlab("Altitude difference") +
  ylab("Sorensen beta-diversity") +
  labs(
    title = paste0(Taxa,": Sorensen beta-diversity VS altitude difference between the plots.")
  )

# Summary
Summ_p1_sor <- Data %>%
  get_summary_stats(beta.sor,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "beta.sor",           # Sub scenario 
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 4,                        # Number of digits 
                 size = 4,                          # Size of the text
                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank(),
                                                   legend.position = "right")) %>% as.ggplot() 

# Arrange them together
sor <- ggarrange(p1_sor, Summ_p1_sor, nrow = 2, heights = c(3,1), legend = "right")

# --- #

p2_sor <- Data %>% 
  ggplot(aes(x = delta_z, y = beta.sim)) +
  stat_binhex(na.rm = T,
              bins = 50) +
  scale_fill_viridis_c() +
  xlab("Altitude difference") +
  ylab("Sorensen turnover component") +
  labs(
    title = paste0(Taxa,": Sorensen turnover VS altitude difference between the plots.")
  )

# Summary
Summ_p2_sor <- Data %>%
  get_summary_stats(beta.sim,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "beta.sim",           # Sub scenario 
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 4,                        # Number of digits 
                 size = 4,                          # Size of the text
                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank(),
                                                   legend.position = "right")) %>% as.ggplot() 

# Arrange them together
sim <- ggarrange(p2_sor, Summ_p2_sor, nrow = 2, heights = c(3,1), legend = "right")

# --- # 

p3_sor <- Data %>% 
  ggplot(aes(x = delta_z, y = beta.sne)) +
  stat_binhex(na.rm = T,
              bins = 50) +
  scale_fill_viridis_c() +
  xlab("Altitude difference") +
  ylab("Sorensen nestedness component") +
  labs(
    title = paste0(Taxa,": Sorensen nestedness VS altitude difference between the plots.")
  )

# Summary
Summ_p3_sor <- Data %>%
  get_summary_stats(beta.sne,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "beta.sne",           # Sub scenario 
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 4,                        # Number of digits 
                 size = 4,                          # Size of the text
                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank(),
                                                   legend.position = "right")) %>% as.ggplot() 

# Arrange them together
sne <- ggarrange(p3_sor, Summ_p3_sor, nrow = 2, heights = c(3,1), legend = "right")

# Arrange the three plots together
plot_sor <- ggarrange(sor,sim,sne,ncol = 3) 

# ----------------------------- #

p1_jac <- Data %>% 
  ggplot(aes(x = delta_z, y = beta.jac)) +
  stat_binhex(na.rm = T,
              bins = 50) +
  scale_fill_viridis_c() +
  xlab("Altitude difference") +
  ylab("Jaccard beta-diversity") +
  labs(
    title = paste0(Taxa,": Jaccard beta-diversity VS altitude difference between the plots.")
  )

# Summary
Summ_p1_jac <- Data %>%
  get_summary_stats(beta.jac,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "beta.jac",           # Sub scenario 
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 4,                        # Number of digits 
                 size = 4,                          # Size of the text
                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank(),
                                                   legend.position = "none")) %>% as.ggplot() 

# Arrange them together
jac <- ggarrange(p1_jac, Summ_p1_jac, nrow = 2, heights = c(3,1), legend = "right")

# --- #

p2_jac <- Data %>% 
  ggplot(aes(x = delta_z, y = beta.jtu)) +
  stat_binhex(na.rm = T,
              bins = 50) +
  scale_fill_viridis_c() +
  xlab("Altitude difference") +
  ylab("Jaccard turnover component") +
  labs(
    title = paste0(Taxa,": Jaccard turnover VS altitude difference between the plots.")
  )

# Summary
Summ_p2_jac <- Data %>%
  get_summary_stats(beta.jtu,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "beta.jtu",           # Sub scenario 
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 4,                        # Number of digits 
                 size = 4,                          # Size of the text
                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank(),
                                                   legend.position = "none")) %>% as.ggplot() 

# Arrange them together
jtu <- ggarrange(p2_jac, Summ_p2_jac, nrow = 2, heights = c(3,1), legend = "right")

# --- # 

p3_jac <- Data %>% 
  ggplot(aes(x = delta_z, y = beta.jne)) +
  stat_binhex(na.rm = T,
              bins = 50) +
  scale_fill_viridis_c() +
  xlab("Altitude difference") +
  ylab("Jaccard nestedness component") +
  labs(
    title = paste0(Taxa,": Jaccard nestedness VS altitude difference between the plots.;")
  )

# Summary
Summ_p3_jac <- Data %>%
  get_summary_stats(beta.jne,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(x = "beta.jne",           # Sub scenario 
                 y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                 digits = 4,                        # Number of digits 
                 size = 4,                          # Size of the text
                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank(),
                                                   legend.position = "right")) %>% as.ggplot() 

# Arrange them together
jne <- ggarrange(p3_jac, Summ_p3_jac, nrow = 2, heights = c(3,1), legend = "right")

# Arrange the three plots together
plot_jac <- ggarrange(jac,jtu,jne,ncol = 3)

# Return the two pages of plots
Total <- list(plot_sor,plot_jac) %>%
  return()

}

# Mosses
Mosses_Beta_Plot <- Beta_Taxo_Plot(Data = Mosses_Beta, Taxa = "Mosses")
# Liverworts
Liverworts_Beta_Plot <- Beta_Taxo_Plot(Data = Liverworts_Beta, Taxa = "Liverworts")
# Bryophytes
Bryophytes_Beta_Plot <- Beta_Taxo_Plot(Data = Bryophytes_Beta, Taxa = "Bryophytes")
# Tracheophytes
Tracheo_Beta_Plot <- Beta_Taxo_Plot(Data = Tracheo_Beta, Taxa = "Tracheophytes")

##### ----- 1.B.2.b : Final Plot ------------------- #####

# --- Save the plots --- #

# Append all the plots together in the same list
pl <- flatten(list(Mosses_Beta_Plot,Liverworts_Beta_Plot,Bryophytes_Beta_Plot,Tracheo_Beta_Plot))
# Transform them into grobs
pl <- lapply(pl, as.grob)
# Create the pdf
ml <- marrangeGrob(pl,nrow=1,ncol = 1)
# Save a global plot
ggsave(filename = "STEP1_Beta_Taxo_Richness.pdf",
       plot = ml,
       device = "pdf",
       path = paste0(getwd(),"/Plots"),
       width = 45,
       height = 25,
       units = "cm")

# _________________________________________ ####

#-----------------------------------------#
##### STEP 2 : PHYLOGENETIC DIVERSITY #####
#-----------------------------------------#

##### -- Creation of the global dataframe -- #####

##### ----- 1.A: Phylogenetic alpha ------------------- #####

# Creation of a function to compute the PD, the MPD and the MNTD
# Take as parameter : 
# - Data : A dataframe containing the Occurrence data 
# - Sp_names : A character vector of species names
# - Tree : The phylogenetic tree associated with the communities
# - Type : The type of data we are working on : "B" for bryophytes or "T" for tracheophytes. 

Alpha_Phylo <- function(Data,Sp_names,Tree,Type = "B") {
  
  # Subsetting of the global dataset to only keep the community values
  Community <- if(Type == "B"){ 
      Data %>% 
        column_to_rownames(var = "Site_VDP") %>% # Change the row_names to be the plot number/name.
        dplyr::select(all_of(Sp_names)) # Remove the meta_data columns.
    } else {
      rownames(Data) <- NULL
      Data %>% 
        column_to_rownames(var = "X") %>% 
        dplyr::select(all_of(Sp_names))
    }
  
  # Computation of the Faith PD
  Faith_PD <- picante::pd(Community, tree = Tree, include.root = T) %>%
    subset(select = -SR) 
  # Computation of the Shannon index
  Sp_Shannon <- vegan::diversity(dplyr::select(Data,all_of(Sp_names)),index = "shannon")
  
  # Add the results to the initial dataframe
  Data <- mutate(Data,Sp_Richness,Sp_Simpson,Sp_Shannon,.after = z)
  
  # Return the dataframe containing the results
  return(Data)
  
}

# Mosses
Mosses_Data <- Alpha_Phylo(Data = Mosses_Data,Moss_names,Tree = Mosses_Tree,Type = "B")

# Liverworts
Liverworts_Data <- Alpha_Taxo(Data = Liverworts_Data,Liver_names)
# Bryophytes
Bryophytes_Data <- Alpha_Taxo(Data = Bryophytes_Data,Bryo_names)
# Tracheophytes
Tracheo_Data <- Alpha_Taxo(Data = Tracheo_Data,Tracheo_names)

# Message
cat(rule(left = "- Taxonomic alpha diversity added - ", line_col = "white", line = " ", col = "green"))



# Load the environmental data
ACP_Data <- read.csv(file = "Data/Bryophytes/Env_Var_Bryo.csv", sep = ";", row.names = 1)

# Split the data between the different "types" of variables (the +3 allows to avoid the two first columns that are x, y and z coordinates)
# Types available on Env_Var_Names.csv BUT careful that mnt_mean was moved as "z". 
ACP_Edaphic <- ACP_Data[,c(1:3,6:8)+3]
ACP_Clim <- ACP_Data[,c(4,21:44,67:68)+3]
ACP_Topo <- ACP_Data[,c(5,45:52,66)+3]
ACP_LC <- ACP_Data[,c(9:20,53:65)+3]

# Make the ACPs
ACP_Edaphic <- PCA(ACP_Edaphic, scale.unit = TRUE, ncp = 5, graph = TRUE)  # Auto correlated 
ACP_Clim <- PCA(ACP_Clim, scale.unit = TRUE, ncp = 5, graph = TRUE) # Auto correlated
ACP_Topo <- PCA(ACP_Topo, scale.unit = TRUE, ncp = 5, graph = TRUE) # Not auto correlated
ACP_LC <- PCA(ACP_LC, scale.unit = TRUE, ncp = 5, graph = TRUE) # Not auto correlated

# Create the data set to make the analysis on. 
ACP_Data <- data.frame(
  ACP_Data[,1:3],
  "Edaphic" = ACP_Edaphic$ind$coord[,1],
  "Climatic" = ACP_Clim$ind$coord[,1],
  "Topo" = ACP_Topo$ind$coord[,1],
  "LC" = ACP_LC$ind$coord[,1])


##### ----- 3.A: GDM data ------------------- #####

# Create a dataframe of identifiers
IDs <- data.frame(Occ_Data[,1:2],"ID" = c(1:nrow(Occ_Data)))

# Create a fake distance matrix for testing
GDM_dist <- as.matrix(dist(Occ_Data[,5:ncol(Occ_Data)]))
GDM_dist <- cbind("ID" = IDs$ID,GDM_dist)

# Without ACPs 
GDM_Env_Data <- as.matrix(Env_Data)
GDM_Env_Data <- cbind("ID" = IDs$ID,GDM_Env_Data)

# With ACPs
GDM_ACP_Data <- as.matrix(ACP_Data)
GDM_ACP_Data <- cbind("ID" = IDs$ID,GDM_ACP_Data)

# --- 3D Distances between plots --- ###

# Load a needed function
source("/home/thibault/Documents/Scripts/Useful_function/PW_to_Vector.R")

# We have the x, y and z (Mnt_Mean) coordinates of each plots, we need to compute the distances between the plots in this 3D environnement. 

# Create the matrix
XYZ_Coord <- data.frame("x" = Env_Data$x, "y" = Env_Data$y, "z" = Env_Data$mnt_mean)
# Compute the distance
XYZ_Dist <- dist(XYZ_Coord, diag = FALSE, upper = FALSE)
# Create a matrix version to remove the upper.triangle and the diagonal
XYZ_Dist_matrix <- as.matrix(XYZ_Dist)                      # Create the matrix
# Transform it into a vector without the diagonal and the duplicated values
XYZ_Dist_Lower <- PW_to_Vector(XYZ_Dist_matrix,Colname = "3DDistance", Diag = F)
# Remove the "Index" column
XYZ_Dist_Lower <- as.matrix(select(XYZ_Dist_Lower,-Index))
# Sort the data frame to be compatible with format site pairs. 
XYZ_Dist_Lower <- XYZ_Dist_Lower[order(XYZ_Dist_Lower[,"PlotA"]),]

# Mathematical formula of the distance between 3points in space
ThreeDist <- sqrt((XYZ_Coord[2,1] - XYZ_Coord[1,1])^2 + (XYZ_Coord[2,2] - XYZ_Coord[1,2])^2 + (XYZ_Coord[2,3] - XYZ_Coord[1,3])^2) # IT WORKS.

# --- STEP 1: Format the data for the GDM --- #

# bioData: The response variable.
# In our case, a pairwise matrix of phylogenetic distances between communities

# bioFormat: Inform the function the format of biodata

# preData : A site by predictor matrix containing all the GDM predictors in column. 

# In the matrices must be present the site ID and the geographical coordinates. 

# Without ACP
GDM_Data <- formatsitepair(bioData = GDM_dist,  # The biological data, a site by site distance matrix
                           bioFormat = 3,       # Inform the type of biological data
                           dist = "bray",
                           abundance = FALSE,
                           siteColumn = "ID",
                           XColumn = "x",
                           YColumn = "y",
                           predData = GDM_Env_Data,
                           verbose = T,
                           sampleSites = 1      # The percentage of sites to work with
)

# With ACP
GDM_ACP_Data <- formatsitepair(bioData = GDM_dist,  # The biological data, a site by site distance matrix
                               bioFormat = 3,       # Inform the type of biological data
                               dist = "bray",
                               abundance = FALSE,
                               siteColumn = "ID",
                               XColumn = "x",
                               YColumn = "y",
                               predData = GDM_ACP_Data,
                               verbose = T,
                               sampleSites = 1      # The percentage of sites to work with
)

# Add the 3Ddistance to the dataframe as a predictor variable
GDM_Data <- as.data.frame(cbind(GDM_Data,XYZ_Dist_Lower[,"3DDistance"])) # On peut pas faire ça, ça casse le format issus de format site pair :/ 

# --- STEP 2: Launch the GDM --- # 

# Re-scale the response data values to be comprised between 0 and 1.
GDM_Data$distance <- rescale(GDM_Data$distance)
GDM_ACP_Data$distance <- rescale(GDM_ACP_Data$distance)

# Launch the GDM
GDM <- gdm(data = GDM_Data, geo = T)
GDM_ACP <- gdm(data = GDM_ACP_Data, geo = T)

# --- STEP 3: Explore the GDM results --- # 

# Explore the summary results of the GDM
summary(GDM)
summary(GDM_ACP)

# Plot the results of the GDM
# The maximum height of each spline indicates the magnitude of total biological change along that gradient and thereby corresponds 
# to the relative importance of that predictor in contributing to biological turnover while holding all other variables constant 
# (i.e., is a partial ecological distance). The spline’s shape indicates how the rate of biological change varies with position along that gradient. 
# Thus, the splines provide insight into the total magnitude of biological change as a function of each gradient and where along each gradient those changes are most pronounced. 
plot(GDM)

#--------------------------------------------#
##### STEP 1.b : Map of species richness #####
#--------------------------------------------#

# Load the map
Swiss_Map <- raster('Data/mnt_mean.tif')

# Plot the map using base R
raster::plot(Swiss_Map,
             col = viridis(200))

# Add the richness data 
points(
    pch = 16,
    Env_Data[,1:2],
    # size the points according to the temp value
    cex = Sp_Richness/15,
    col = alpha("black",0.7)
    )

# POINTS ARE NOT SCALING WITH THE MAP











