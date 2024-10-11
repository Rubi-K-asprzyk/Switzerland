#!/usr/bin/env Rscript

# --------------------------- #
# Swiss Patterns / K.Thibault #
# --------------------------- #

#------------------------#
##### INITIALISATION #####
#------------------------#

 
{ # !! Initialisation

  suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
  # Display a beginning message. 
  cat(rule(left = "SCRIPT SWISSPATTERNS.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 
  cat(rule(left = "INITIALISATION", line_col = "green", line = "-", col = "br_green"))
  
  # Set Working Directory
  setwd("/home/thibault/Documents/PhD-Thesis/Research/Switzerland")

  # Install/load pacman.
  suppressPackageStartupMessages(if (!require(pacman)) {
    install.packages("pacman")
    library(pacman)
  })

  # Install/load tons of packages.
  p_load(
    doParallel, # Allow parallel computation
    ape, # Multiple tools of phylogenetic analyses
    dplyr, # Allow the use of the plyr syntax
    tidyr, # Allow the use of function "pivot_longer"
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
    abdiv,
    patchwork,
    multcomp,
    multcompView,
    emmeans,
    stringr,
    gridExtra,
    ggplotify,
    data.table,
    ggh4x,
    readr
  )

  # Set the parallel backend
  registerDoParallel(cores = 10)

  # Set a global option to disable dplyr::summarise() warnings when grouping.
  options(dplyr.summarise.inform = TRUE)

  # Call the functions needed to transform the PW matrices into vector
  source(paste0(getwd(), "/Utilities/PW_to_Vector.R"))

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
      margin = margin(t = 10, b = 10)
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
    theme(
      plot.background = element_rect(fill = "grey98", color = "grey98"),
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
        margin = margin(l = 10, t = 10, b = 20),
        hjust = 0.5
      ),
      # Title and caption are going to be aligned
      plot.caption = element_text(
        color = "grey30",
        size = 8,
        lineheight = 1.2,
        hjust = 0.1,
        margin = margin(t = 10, b = 10) # Large margin on the top of the caption.
      )
    )
  } # End of arrange_theme

  # Message
  cat(rule(left = "- Theme set - ", line_col = "white", line = " ", col = "green"))

  #### . Custom Function . #####

  # Create a function to format pariwise results from Hardy Metrics
  # ARGUMENTS:
  # - Hardy.result: A pairwise result from SpacodiR
  # - Name: The name of the metric

  Hardy_format <- function(Hardy.result, Name) {
    Metric <- Hardy.result %>%
      # Change colnames and rownames to be a simple number (and not "plt.X")
      `colnames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
      `rownames<-`(sub(x = rownames(.), pattern = "plt.", replacement = "")) %>%
      # Transform the pairwise
      PW_to_Vector(Colname = Name)
    # mutate_at(c("PlotA", "PlotB"), as.numeric)
    # # Join the metadata for both plots used in the pairwise computation.
    # left_join(MetaData, c("PlotA" = "Plot")) %>% # Plot A
    # left_join(MetaData, c("PlotB" = "Plot"), suffix = c("_A", "_B")) %>% # Plot B + add suffixes "A" and "B" to distinguish the two plots.
    # # - Add the absolute difference of altitude between the two plots (delta_z) - #
    # mutate(delta_z = abs(z_A - z_B), .after = z_B) %>%
    # # - Add the complete 3D-distance between the two plots (delta_xyz) - #
    # mutate(delta_xyz = sqrt((x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2), .after = delta_z)
  }

  #### . Data Preparation . #####

  # -- Bryophytes -- #

  # Load the datasets
  # Occurence
  bryoData <- read.csv("Data/FinalMatrixLv95_SynChangedV3.csv", row.names = 1) %>%
    # Rename "Sites_notre_num_rotation_"
    rename("Site_VDP" = "Sites_notre_num_rotation_")
  # Traits
  bryoTraits <- read.csv("Data/Bryophytes_Traits_SpRichness.csv", row.names = 1)
  # Environment Variables
  bryoEnv <- read.csv("Data/EnvironmentalValues_Bryophytes.csv", row.names = 1) %>% # mnt_mean was changed directly in the dataset to "z"
    # Rename "Sites_notre_num_rotation_"
    rename("Site_VDP" = "Sites_notre_num_rotation_")
  # Load the phylotrees
  Liver_Tree <- read.tree("PhyloTree/timetree50mod-liverwortsV2.nwk")
  Mosses_Tree <- read.tree("PhyloTree/timetree50mod-mossesV2.nwk")

  # Load the Null Phylotrees
  Liver_Tree_Null <- read.tree("PhyloTree/timetree-liverworts-NM.tree")
  Mosses_Tree_Null <- read.tree("PhyloTree/timetree-mosses-NM.tree")

  # -- Split between Mosses and Liverworts -- #

  # Load the names of the mosses present in the dataframe of Occ_Data_Moss for further splitting between Mosses and Liverworts.
  Bryo_Names <- read.csv(file = "Utilities/Bryophyte_Names.csv", row.names = 1)
  Moss_Names <- dplyr::select(filter(Bryo_Names, Taxa == "Moss"), Species) %>% as.matrix()
  Liver_Names <- dplyr::select(filter(Bryo_Names, Taxa == "Liverwort"), Species) %>% as.matrix()

  # Extract the meta data names
  Bryo_Meta <- colnames(bryoData)[1:4]

  # Split the Occ_Data_Moss between Mosses and Liverworts
  liverData <- dplyr::select(bryoData, any_of(c(Bryo_Meta, Liver_Names)))
  mossData <- dplyr::select(bryoData, any_of(c(Bryo_Meta, Moss_Names)))

  # Get the names of the liverworts and mosses present
  Moss_Names <- colnames(mossData)[-c(1:4)]
  Liver_Names <- colnames(liverData)[-c(1:4)]

  # Message
  cat(rule(left = "- Data loaded - ", line_col = "white", line = " ", col = "green"))

} # !! End of the initialisation.

#--------------------------------------------------#
##### RANGE-ANALYSES: CREATION OF THE FUNCTION #####
#--------------------------------------------------#


  # . Creation of a function "range-analyses" to realize all the wanted analyses for each taxa and threshold wanted . # 

  # - ARGUMENTS - #
    # - Data: A community matrix of occurence with plot in rows and species in column (with of without metadata associated).
    # - Sp_Names: A vector of species names. 
    # - Phylo_Tree : Phylogenetic Tree associated with the dataset. 
    # - Threshold: The minimum number of species wanted in each plot. (Default = 2).
    # - Nrep: The number of NullModel Replicates to do. 
    # - Alt_limits: A vector of altitudinal limits to split the altitudinal gradient into stripes. (Default = NULL).
    # - N_stripes: If {Alt_limits == NULL}, the number of stripes to split the altitudinal gradient. (Default = 6).
    # - Weighted: If {Alt_limits == NULL}, Are the stripes weighted by the number of plots (i.e with varying size) or not (i.e varying in the number of plots). (Default = T). 

  # - OUTPUTS - #

  # - TEST - #

    # - Mosses - # 
  Data <- mossData
  Sp_Names <- Moss_Names
  Phylo_Tree <- Mosses_Tree
  Phylo_Tree_Null <- Mosses_Tree_Null
  Nrep <- 3
  Threshold <- 4 # Mosses only accept a threshold of 4
  Alt_limits <- c(0,1000,1400,1800,2000,2200,Inf)

    # - Liverworts - # 
  Data <- liverData
  Sp_Names <- Liver_Names
  Phylo_Tree <- Liver_Tree
  Phylo_Tree_Null <- Liver_Tree_Null
  Nrep <- 3
  Threshold <- 3 # Mosses only accept a threshold of 4
  Alt_limits <- c(0,1400,2000,Inf)

    # - Global - # 

  N_stripes <- 6
  Weighted <- TRUE


  # - COMPUTED METRICS - # 
      # -- ALPHA -- # 
    # - Species Richness (SR)
    # - Gini-Simpson Index (GS)
    # - Phylogenetic Diversity (PD)
    # - Mean Phylogenetic Diversity (MPD)
    # - Mean Nearest Neighbour Distance (MNTD)

      # -- BETA -- #
    # - Simpson Total (beta.sim)
    # - Simpson turnover (beta.sim)
    # - Simpson nestedness (beta.sim)


range_analyses <- function(Data, Sp_Names, Threshold,Phylo_Tree,Phylo_Tree_Null,Nrep, Alt_limits = NULL, N_stripes = 6, Weighted = TRUE) {

  #-------------------------------#
  ##### ALPHA METRICS ANALYSES ####
  #-------------------------------#

  { # !! Alpha Metric Computation 

  # Computation of the metrics and filtering of the data based on the specie richness Threshold
  Alpha.Obs <- Data %>%
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

  # Add the OBSERVED Phylogenetic metrics to the dataframe
  Alpha.Obs <- Alpha.Obs %>%
    # Compute the mean phylogenetic distance
    mutate(MPD = picante::mpd(samp = .[,colnames(Data.Diss)], dis = Data.Diss), .after = y) %>%
    # Compute the mean nearest neighbour distance
    mutate(MNTD = picante::mntd(samp = .[,colnames(Data.Diss)], dis = Data.Diss), .after = y) %>%
    # Add the "Type" of the metric (OBS vs NM)
    mutate(Type = "OBS")

    # --- Compute the NULL MODEL values --- #

      # Message
      cat(rule(left = "- BEGIN: Computation of Phylo-Alpha Metrics NULL MODEL - ", line_col = "#3C917E", line = " ", col = "#3C917E"))

      # /!\ .combine uses a custom function with "x" the previous iteration and "y" the new iteration. It is used to give the "by" argument without display message and clutter the console.

    Alpha.Null <- foreach(i = 1:Nrep, .combine = function(x,y){full_join(x,y, by = join_by(Site_VDP, Site_Suisse, x, y, Type, Metric, Value))}) %dopar% {

      # Compute and reorder the cophenetic distances
      Data.Diss <- cophenetic(Phylo_Tree_Null[[i]])
      # Only keep the species that are present in the treehttps://www.twitch.tv/solaryhs
      Sp_Names_Present <- Sp_Names[Sp_Names %in% colnames(Data.Diss)]
      # Reorder the distance matrix
      Data.Diss <- Data.Diss[Sp_Names_Present, Sp_Names_Present]

      # Compute the MPD
      Alpha.NM <- Data %>%
        # Compute the mean phylogenetic distance
        mutate(MPD = picante::mpd(samp = .[, colnames(Data.Diss)], dis = Data.Diss), .after = y) %>%
        # Compute the mean nearest neighbour distance
        mutate(MNTD = picante::mntd(samp = .[, colnames(Data.Diss)], dis = Data.Diss), .after = y) %>%
        # Add the "Type" of the metric (OBS vs NM)
        mutate(Type = paste0("NM")) %>%
        # Only select the wanted metrics
        dplyr::select(any_of(c("Site_VDP", "Site_Suisse", "x", "y", "Type", "MPD", "MNTD"))) %>%
        # Pivot longer the metrics
        pivot_longer(cols = c(MPD, MNTD), values_to = "Value", names_to = "Metric")

      # Print a message
      cat(rule(left = paste0("- Rep ",i," / ",Nrep," - "), line_col = "#73913C", line = " ", col = "#73913C"))

      # return
      return(Alpha.NM)
      
    } # End of Null Models

    # Message
    cat(rule(left = "- DONE: Computation of Phylo-Alpha Metrics NULL MODEL - ", line_col = "#3C917E", line = " ", col = "#3C917E"))

    # Combine Observed and Null Results
    Alpha.Obs <- Alpha.Obs %>%
    # Pivot longer the metrics
      pivot_longer(cols = c(SR,GS,PD,MPD,MNTD), values_to = "Value", names_to = "Metric")

    Data.filtered <- full_join(Alpha.Obs, Alpha.Null, by = join_by(Site_VDP, Site_Suisse, x, y, Type, Metric, Value)) %>%
      # Relocate the wanted columns.
      relocate(any_of(c("Type", "Metric", "Value")), .after = y) %>%
      # Remove the Number contained in the "Type" column.
      mutate(across(Type,~ gsub('[[:digit:]]+', '', .)))


  # Message
  cat(rule(left = paste0("- Data filtered based on species richness / THRESHOLD = ",Threshold," - "), line_col = "white", line = " ", col = "green"))
  cat(paste0(" - Number of initial plots: ",length(unique(Data$Site_VDP))," -\n"))
  cat(paste0(" - Number of plot left: ",length(unique(Alpha.Obs$Site_VDP))," -\n"))
  cat(paste0(" - Number of plot eliminated: ",(length(unique(Data$Site_VDP)) - length(unique(Alpha.Obs$Site_VDP)))," -\n"))



  } # !! End of Alpha Metric Computation

  # ------------------------------------------------------------------------------------------------------------------------------ #

  # ----- RANGE CHOICE ----- #

  { # !! Range Choice 

    # Creation of the ranges based on the parameters entered in the function.

  # If Alt_limits are entered, create the dataframe based on them.
  if (is.null(Alt_limits) == FALSE) {
    # Create a dataframe of intervals as well as corresponding stripe number
    Range_Intervals <- paste0("[", Alt_limits[-length(Alt_limits)], ";", Alt_limits[-1], "[") %>%
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

  } # !! End of Range Addition 

  #------------------------------#
  ##### BETA METRICS ANALYSES ####
  #------------------------------#

  { # !! Taxonomic Beta Metric Computation 

  # Add the rowname as a column to keep this info
  Beta.Data <- Data.filtered %>%
    # Remove all "NM" values
    filter(Type == "OBS") %>%
    # Pivot back the Alpha metrics computed to remove them (Or all the lines will be multiplied)
    pivot_wider(names_from = Metric, values_from = Value) %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric)

  # Save the Zdist for mantel tests before vectorization
  Zdist <- dist(Beta.Data$z)

  # Get the Metadata
  MetaData <- Beta.Data %>% 
    # Select only the Metadata (based on the species names )
    dplyr::select(!any_of(Sp_Names))

  # -- Sorensen -- # 

  Sorensen <- Beta.Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Compute the Sorensen index
    beta.pair(x = ., index.family = "sorensen") %>%
    # Transform into a matrix
    lapply(as.matrix)

  # Save it for mantel tests before vectorization
  Sorensen_Mantel <- Sorensen 

    # Transform it into a vector
  Sorensen <-
    lapply(names(Sorensen), function(x) {
      PW_to_Vector(Sorensen[[x]], Colname = x)
    }) %>%
    # Merge the results altogether
    purrr::reduce(merge)

  # -- Jaccard -- # 

  Jaccard <- Beta.Data %>%
    # Select only the occurence data
    dplyr::select(!colnames(MetaData)) %>%
    # Compute the Sorensen index
    beta.pair(x = ., index.family = "jaccard") %>%
    # Transform into a matrix
    lapply(as.matrix) 

    # Save it for mantel tests before vectorization
  Jaccard_Mantel <- Jaccard 

      # Transform it into a vector
  Jaccard <-
    lapply(names(Jaccard), function(x) {
    PW_to_Vector(Jaccard[[x]], Colname=x)}) %>%
    # Merge the results altogether
    purrr::reduce(merge)

  # Change the metadata to NOT contains the alpha metrics results
  MetaData <- dplyr::select(MetaData, !unique(Data.filtered$Metric))

  # -- Bind the metrics together -- #

  Taxo.Beta <- merge(Sorensen,Jaccard)

  } # !! End of Taxonomic Beta Computation

  # --------------------------------------------------------------------------------------------------------------------------------- #

  { # !! Phylogenetic Beta Metric Computation 

  # Add the rowname as a column to keep this info
  Beta.Data <- Data.filtered %>%
    # Remove all "NM" values
    filter(Type == "OBS") %>%
    # Pivot back the Alpha metrics computed to remove them (Or all the lines will be multiplied)
    pivot_wider(names_from = Metric, values_from = Value) %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric) %>%
    # Select only the occurence data
    dplyr::select(any_of(Sp_Names)) %>%
    # Transpose the data for the analysis
    t() %>%
    # Add the colnames as the plot numbers
    `colnames<-`(1:ncol(.))

  # Get the Metadata
  MetaData <- Data.filtered %>%
    # Create a "Plot" column
    rownames_to_column(var = "Plot") %>%
    mutate_at("Plot", as.numeric) %>% 
    # Select only the Metadata (based on the species names )
    dplyr::select(!any_of(c(Sp_Names,"Metric","Value","Type")))

# -- Compute the PIst metrics -- # 

  Hardy_Metrics <- spacodiR::spacodi.calc(
      sp.plot = Beta.Data,         # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
      phy = Phylo_Tree,             # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                   # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
      all.together = TRUE,     # whether to treat all traits together or separately
      prune = TRUE,
      pairwise = TRUE)

# -- Extract all Metrics -- #
  
  ## PST ##

Pst <- Hardy_format(Hardy.result = Hardy_Metrics$pairwise.Pst, Name = "Pst") 

  ## PIST ##

PIst <- Hardy_format(Hardy.result = Hardy_Metrics$pairwise.PIst, Name = "PIst") 
    
# --- Combine the three altogether and pivot longer the metrics --- # 

Phylo.beta.Obs <- purrr::reduce(list(Pst,PIst),function(x, y) {
  full_join(x, y, by = join_by(Sample, PlotA, PlotB))}) %>%
  # Pivot the metrics
  pivot_longer(cols = c(Pst,PIst), values_to = "Value", names_to = "Metric") %>%
  # Mutate the type of metric, "OBS vs NM"
  mutate(Type = "OBS")

# --- Compute the NULL MODEL values --- #

  # Message
  cat(rule(left = "- BEGIN: Computation of Phylo-Beta Metrics NULL MODEL - ", line_col = "#3C917E", line = " ", col = "#3C917E"))

  # /!\ .combine uses a custom function with "x" the previous iteration and "y" the new iteration. It is used to give the "by" argument without display message and clutter the console.
Phylo.beta.Null <- foreach(i = 1:Nrep, .combine = function(x, y) {
  full_join(x, y, by = join_by(Sample, PlotA, PlotB, Metric, Value, Type))}
) %dopar% {

  # -- Compute the PIst metrics -- #

  Hardy_Metrics <- spacodiR::spacodi.calc(
    sp.plot = Beta.Data, # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
    phy = Phylo_Tree_Null[[i]], # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                   # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
    all.together = TRUE, # whether to treat all traits together or separately
    prune = TRUE,
    pairwise = TRUE
  )

  # -- Extract all Metrics -- #

  ## PST ##

  Pst <- Hardy_format(Hardy.result = Hardy_Metrics$pairwise.Pst, Name = "Pst")

  ## PIST ##

  PIst <- Hardy_format(Hardy.result = Hardy_Metrics$pairwise.PIst, Name = "PIst")

  # --- Combine the three altogether and pivot longer the metrics --- #

  Phylo.beta.Null <- purrr::reduce(list(Pst, PIst), function(x, y) {
  full_join(x, y, by = join_by(Sample, PlotA, PlotB))}) %>%
    # Pivot the metrics
    pivot_longer(cols = c(Pst, PIst), values_to = "Value", names_to = "Metric") %>%
    # Mutate the type of metric, "OBS vs NM"
    mutate(Type = "NM")

  # Print a message
  cat(rule(left = paste0("- Rep ",i," / ",Nrep," - "), line_col = "#73913C", line = " ", col = "#73913C"))

  # Return the NM values
  return(Phylo.beta.Null)

} # End of Beta Null

  # Message
  cat(rule(left = "- END: Computation of Phylo-Beta Metrics NULL MODEL - ", line_col = "#3C917E", line = " ", col = "#3C917E"))

  # Combine Observed And Null Model results
  Phylo.beta <- full_join(Phylo.beta.Obs,Phylo.beta.Null,by = join_by(Sample, PlotA, PlotB, Metric, Value, Type)) 

    # ------------- #

  # Combine the taxonomic and phylogenetic beta results

  Beta.results <- full_join(Phylo.beta,Taxo.Beta,by = join_by(Sample, PlotA, PlotB)) %>%
  mutate_at(c("PlotA", "PlotB"), as.numeric) %>%
  # Join the metadata for both plots used in the pairwise computation.
  left_join(MetaData, c("PlotA" = "Plot")) %>% # Plot A
  left_join(MetaData, c("PlotB" = "Plot"), suffix = c("_A", "_B")) %>% # Plot B + add suffixes "A" and "B" to distinguish the two plots.
  # - Add the absolute difference of altitude between the two plots (delta_z) - #
  mutate(delta_z = abs(z_A - z_B), .after = z_B) %>%
  # - Add the complete 3D-distance between the two plots (delta_xyz) - #
  mutate(delta_xyz = sqrt((x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2), .after = delta_z) %>%
  # Relocate columns
  dplyr::select(order(colnames(.))) %>%
  # Mutate the type of metric, "OBS vs NM"
  relocate(any_of(c("Metric","Value","Type")),.before = 1) %>%
  relocate(any_of(c("delta_z","delta_xyz")),.after = z_B) %>%
  # Create a column that combines the two stripes
  mutate(Range_Number_AB = paste0(Range_Number_A,"-",Range_Number_B)) %>%
  # Create a column that is the stripe distance between the two plots.
  mutate(Range_Number_Distance = as.factor(abs(as.numeric(Range_Number_A)-as.numeric(Range_Number_B))))


} # !! End of Phylogenetic Beta Computation

  # Creation of a dataset with only observed values
  Beta.results.Obs <- Beta.results %>%
    # Filter the observed values
    filter(Type == "OBS")
  
  # --------------------------------------------------------------------------------------------------------------------------------- #

  #----------------------------------#
  #### ----- PLOTTING ALPHA ----- ####
  #----------------------------------#

    # Message
  cat(rule(left = "- BEGIN: Plotting - ", line_col = "white", line = " ", col = "green"))

  { # !! FIGURE 1: Alpha Table ~ Range 

# Create a global dataframe containing the results for all taxa. 
Data.filtered.summarize <- Data.filtered %>%
  # Select only the wanted data
  dplyr::select("Site_VDP":"Type") %>%
  # Group the data
  group_by(Range_Number, Metric, Type) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  mutate(Sum_Value = sum(Value), Mean_Value = mean(Value),Var_Value = var(Value), Count=n()) %>%
  # -!- Reorder the breaks -!- #
  group_by(Metric) %>% 
  # Arrange by the altitude to order in the same time the breaks
  arrange(z, .by_group = TRUE) %>%
  mutate(Break = factor(Range_Number, levels = unique(.$Range_Number)))


# Compute and create the summary table of the metric(s) splitted between stripes for each taxa
Alpha.summary <- Data.filtered %>%
  # Filter only for the Observed Values
  dplyr::filter(Type == "OBS") %>%
  # Group the data
  group_by(Range_Intervals, Metric) %>%
  # Compute the summary stat
  get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
  # Draw the summary statistics
  ggsummarytable(
    x = "Range_Intervals",                           # Split by stripes
    y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
    digits = 2,                            # Number of digits 
    size = 4,                             # Size of the text
    color = "Range_Intervals",                       # Color by stripes
    facet.by = c("Metric"),
    scales = "free_x", 
    ggtheme = arrange_theme() + theme(legend.position = "none")
    ) +
  # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  # Labs
  labs(subtitle = "Metrics summary ~ Altitudinal ranges entered.\nn = Number of plots.")
    
# Compute and create the summary table of the metric(s) unsplitted between stripes
Alpha.summary.unsplitted <- Data.filtered %>%
  # Filter only for the Observed Values
  dplyr::filter(Type == "OBS") %>%
  # Group the data
  group_by(Metric) %>%
  # Compute the summary stat
  get_summary_stats(Value, type = "common") %>% # Compute the summary statistics for each groups
  # Draw the summary statistics
  ggsummarytable(
    x = "Metric", # Split by stripes
    y = c("n", "min", "max", "mean", "sd"), # Choose the metrics we want to display
    digits = 2, # Number of digits
    size = 5, # Size of the text
    scales = "free_x",
    ggtheme = arrange_theme() + # Theme
      theme(legend.position = "none")
  ) +
  # Labs
  labs(subtitle = "Metrics summary total.\nn = Number of plots.")

# Combination of both plots together.
Figure_1 <- as.ggplot(ggarrange(
  Alpha.summary,Alpha.summary.unsplitted,common.legend = T, legend = "right",widths = c(2, 1)))

   # Print a message
  cat(rule(left = "Figure 1: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 1

  # --------------- # 

  { # !! FIGURE 2: Scatterplot of Observed Metric Values ~ Range_Intervals - # 

# Find the mean value of z for x_axis centering and Value for y-axis
Value.Mean <- Data.filtered %>%
  # Filter only for the Observed Values
  dplyr::filter(Type == "OBS") %>%
  # Group the data
  group_by(Metric,Range_Intervals) %>%
  # Compute the Sum of the Metric and the number of plots for each group. 
  summarise(z = mean(z), Mean_Value = mean(Value,  na.rm = TRUE), .groups = "keep") %>%
  # reorder z by group to have the breaks in the good order
  arrange(z, .by_group = TRUE)

# Draw the plots.
Figure_2 <- Data.filtered %>%
  # Filter only for the Observed Values
  dplyr::filter(Type == "OBS") %>%
  # Aes
  ggplot(aes(x = z, y = Value, color = Range_Intervals)) + # Reorder the factors correctly
  # Plot all the values
  geom_point() +
  # Facet the plot
  facet_grid(Metric ~ .,
    scales = "free_y"
  ) +
  # Plot the connected scatter plot of the mean values for each altitudinal stripe.
  new_scale_color() + # Define scales before initiating a new one
  # Use the new scale
  geom_point(data = Value.Mean, aes(y = Mean_Value), size = 2) +
  geom_line(data = Value.Mean, aes(y = Mean_Value), linewidth = 1, alpha = 0.7, group = 1) +
  # Guides
  guides(size = "none") +
  # Labels
  xlab("Altitude") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("ScatterPlot of Metrics ~ Altitude, colored by altitudinal ranges."),
    subtitle = paste0(
      "Black dotted line represent the mean metric value for each of the altitudinal stripes."
    )
  )

  # Print a message
  cat(rule(left = "Figure 2: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 2

  # ------ # 

  { # !! FIGURE 3: Boxplots of Metrics Values ~ Altitudinal stripes - # 

# Realization of an ANOVA of the metrics values ~  Range_Intervals * Metric
anova <- aov(Value ~ Range_Intervals * Metric, data = Data.filtered)

# Get (adjusted) means to have results groupes by metrics.
model_means <- emmeans(object = anova,
                       specs = ~ Range_Intervals | Metric) 

# Add significativity letters to each mean.
model_means_cld <- cld(object = model_means,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

# Compute the boxplots stats to have access to the upper hinge value.
  Alpha.Stats.boxplots <- Data.filtered %>% 
  # Filter only for the Observed Values
  dplyr::filter(Type == "OBS") %>%
  # Group_by
  group_by(Range_Intervals,Metric) %>%
  # Extract the values
  summarise(as_tibble_row(quantile(Value, na.rm = TRUE), .name_repair = \(x) paste0('q', parse_number(x))), .groups = "keep") 

  # Combine the model and the stats to use them altogether
  Labels <- full_join(model_means_cld,Alpha.Stats.boxplots,by = join_by(Range_Intervals, Metric))

  # Draw the Boxplots
  Figure_3 <- Data.filtered %>% 
  # Filter only for the Observed Values
  dplyr::filter(Type == "OBS") %>%
  # Aes
  ggplot(aes(x = Range_Intervals, y = Value, fill = Range_Intervals)) + # Reorder the factors correctly
  # Plot all the values
  geom_boxplot(outlier.shape = NA) +
  # Facet the plot
  ggh4x::facet_grid2(. ~ Metric, scales = "free_y", independent = "y") +
  # Add the letter
  geom_label(data = Labels, aes(label = toupper(str_trim(.group)), y = q75, fontface = "bold"), alpha = 0.7, vjust = -0.5) +
  # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  # Labs
  labs(title ="Metrics summary ~ Altitudinal ranges entered.",
    subtitle = "Separatedly per Metric, Metric means by altitudinal range followed by a common letter are not significantly different according to the Tukey-test") 

      # Print a message
  cat(rule(left = "Figure 3: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 3

    # ----- # 

  { # !! FIGURE 4: Boxplots of Null-Models ~ Altitudinal stripes - # 

  # Draw the Boxplots
  Figure_4 <- Data.filtered %>% 
  # Filter only for the phylogenetic metrics
  dplyr::filter(Metric %in% c("MPD","MNTD")) %>%
  # Aes
  ggplot(aes(x = Range_Intervals, y = Value, fill = Type)) + # Reorder the factors correctly
  # Plot all the values
  geom_boxplot(outlier.shape = FALSE) +
  # Facet the plot
  ggh4x::facet_grid2(. ~ Metric, scales = "free_y", independent = "y") +
  # 
      stat_compare_means(
        label = "p.signif",
        method = "wilcox.test",
        symnum.args = list(
          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
          symbols = c("****", "***", "**", "*", "ns")
        )
      ) + 
  # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  # Labs
  labs(title =" Observed VS NullModel Metrics ~ Altitudinal ranges entered.",
    subtitle = "Separatedly per Metric and Range Intervals, Wilcoxon tests of Observed VS NullModel Values") 

          # Print a message
  cat(rule(left = "Figure 4: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 4

  #----------------------------------#
  #### ----- PLOTTING BETA ----- ####
  #----------------------------------#

  { # !! FIGURE 5: Range Boxplot of Metrics ~ Stripe distance -- #

  # Boxplots of the metrics based on the stripe distance between the two plots used for the computation. 

# Realization of an ANOVA of the metrics values ~  Range_Number_Distance * Metric
anova <- aov(Value ~ Range_Number_Distance * Metric, data = Beta.results.Obs)

# Get (adjusted) means to have results groupes by metrics.
model_means <- emmeans(object = anova,
                       specs = ~ Range_Number_Distance | Metric) 

# Add significativity letters to each mean.
model_means_cld <- cld(object = model_means,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

# Compute the boxplots stats to have access to the upper hinge value.
  Beta.Stats.boxplots <- Beta.results.Obs %>% 
  # Group_by
  group_by(Range_Number_Distance,Metric) %>%
  # Extract the values
  summarise(as_tibble_row(quantile(Value, na.rm = TRUE), .name_repair = \(x) paste0('q', parse_number(x))), .groups = "keep") 

  # Combine the model and the stats to use them altogether
  Labels <- full_join(model_means_cld,Beta.Stats.boxplots,by = join_by(Range_Number_Distance, Metric))

  # Draw the Boxplots
  Figure_5 <- Beta.results.Obs %>% 
  # Aes
  ggplot(aes(x = Range_Number_Distance, y = Value, fill = Range_Number_Distance)) + # Reorder the factors correctly
  # Plot all the values
  geom_boxplot(outlier.shape = NA) +
  # Facet the plot
  ggh4x::facet_grid2(. ~ Metric, scales = "free_y", independent = "y") +
  # Add the letter
  geom_label(data = Labels, aes(label = toupper(str_trim(.group)), y = q75, fontface = "bold"), alpha = 0.7, vjust = -0.5) +
  # Theme
  # theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  # Labs
  labs(title ="Metrics summary ~ Range_Number_Distance.",
    subtitle = "Separatedly per Metric, Metric means computed by stripes distance between the two plots used for the computation. For each sub-plot, bars topped by a common letter are not significantly different according to the Tukey-test")

              # Print a message
  cat(rule(left = "Figure 5: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 5

    # ----- # 

  { # !! FIGURE 6 :Boxplot of Metrics ~ Intra Stripe -- #

    # Boxplots of the metrics only for the results computed between plots from the same stripe, called "Intra-Stripe".

# Create the data frame

Beta.results.intra <- Beta.results.Obs %>%
  # Filter the dataset to only keep the pairwise plots from the same stripe
  dplyr::filter(Range_Number_A == Range_Number_B) %>%
  dplyr::filter(PlotA != PlotB) %>%
  mutate_at("Range_Number_AB", as.character) %>%
  # Reorder the stripes to have 10-10 effectively last
  mutate(Range_Number_AB = fct_relevel(Range_Number_AB,mixedsort(unique(.$Range_Number_AB))))

# Compute the boxplots stats to have access to the upper hinge value.
Beta.Stats.boxplots <- Beta.results.intra %>% 
  # Group_by
  group_by(Range_Number_AB,Metric) %>%
  # Extract the values
  summarise(as_tibble_row(quantile(Value), .name_repair = \(x) paste0('q', parse_number(x))), .groups = "keep") 

# Realization of an ANOVA of the metrics values ~  Range_Number_Distance * Metric
anova <- aov(Value ~ Range_Number_AB * Metric, data = Beta.results.intra)

# Get (adjusted) means to have results groupes by metrics.
model_means <- emmeans(object = anova,
                       specs = ~ Range_Number_AB | Metric) 

# Add significativity letters to each mean.
model_means_cld <- cld(object = model_means,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

# Combine the model and the stats to use them altogether
Labels <- full_join(model_means_cld,Beta.Stats.boxplots,by = join_by(Range_Number_AB, Metric))


  # Draw the Boxplots
  Figure_6 <- Beta.results.intra %>% 
  # Aes
  ggplot(aes(x = Range_Number_AB, y = Value, fill = Range_Number_AB)) + # Reorder the factors correctly
  # Plot all the values
  geom_boxplot(outlier.shape = NA) +
  # Facet the plot
  ggh4x::facet_grid2(. ~ Metric, scales = "free_y", independent = "y") +
  # Add the letter
  geom_label(data = Labels, aes(label = toupper(str_trim(.group)), y = q75, fontface = "bold"), alpha = 0.7, vjust = -0.5) +
  # Theme
  # theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  # Labs
  labs(title ="Metrics summary ~ Range_Number.",
    subtitle = "Separatedly per Metric, Metric means computed between two plots from the same stripe. For each sub-plot, bars topped by a common letter are not significantly different according to the Tukey-test")

    cat(rule(left = "Figure 6: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 6

    # ----- # 

  { # !! FIGURE 7: Boxplots of Null-Models ~ Altitudinal stripes - # 

  Beta.results.intra <- Beta.results %>%
  # Filter the dataset to only keep the pairwise plots from the same stripe
  dplyr::filter(Range_Number_A == Range_Number_B) %>%
  dplyr::filter(PlotA != PlotB) %>%
  mutate_at("Range_Number_AB", as.character) %>%
  # Reorder the stripes to have 10-10 effectively last
  mutate(Range_Number_AB = fct_relevel(Range_Number_AB,mixedsort(unique(.$Range_Number_AB))))

  # Draw the Boxplots
  Figure_7 <- Beta.results.intra %>% 
  # Filter only for the Observed Values
  dplyr::filter(Metric %in% c("Pst","PIst")) %>%
  # Aes
  ggplot(aes(x = Range_Number_AB, y = Value, fill = Type)) + # Reorder the factors correctly
  # Plot all the values
  geom_boxplot(outlier.shape = NA) +
  # Facet the plot
  ggh4x::facet_grid2(. ~ Metric, scales = "free_y", independent = "y") +
  # 
      stat_compare_means(
        label = "p.signif",
        method = "wilcox.test",
        symnum.args = list(
          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
          symbols = c("****", "***", "**", "*", "ns")
        )
      ) + 
  # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  # Labs
  labs(title =" Observed VS NullModel Metrics ~ Altitudinal ranges entered.",
    subtitle = "Separatedly per Metric and Range Intervals, Wilcoxon tests of Observed VS NullModel Values") 

        cat(rule(left = "Figure 7: Done.", line_col = "#73913C", line = " ", col = "#73913C"))

  } # !! End of Figure 7

  #-------------------------------#
  ##### MANTEL + PEARSON TEST #####
  #-------------------------------#

  # Mantel Test are only used for beta metrics (to compare matrices against matrices)
  # Pearson Test are only used for alpha metrics

  { # !! Mantel Tests

  # -- Extract all the Metrics we want to keep -- #

      # Hardy Metrics
    PIst <- Hardy_Metrics$pairwise.PIst %>% as.dist()
    Pst <- Hardy_Metrics$pairwise.Pst %>% as.dist()
    Bst <- Hardy_Metrics$pairwise.Bst %>% as.dist()

      # Sorensen Metrics    
    beta.sim <- Sorensen_Mantel$beta.sim %>% as.dist()
    beta.sne <- Sorensen_Mantel$beta.sne %>% as.dist()
    beta.sor <- Sorensen_Mantel$beta.sor %>% as.dist()
      # Jaccard Metrics
    beta.jtu <- Jaccard_Mantel$beta.jtu %>% as.dist()
    beta.jne <- Jaccard_Mantel$beta.jne %>% as.dist()
    beta.jac <- Jaccard_Mantel$beta.jac %>% as.dist()

    # Return the results
    Mantel.Data <- (list(PIst,Pst,Bst,beta.jac,beta.jne,beta.jtu,beta.sim,beta.sne,beta.sor) %>%
              setNames(c("PIst","Pst","Bst","beta.jac","beta.jne","beta.jtu","beta.sim","beta.sne","beta.sor")))

    # PIst, Pst and Bst does not have the same dimension as Zdist. 

  #### --- Computation of the Mantel Test --- # 

  Mantel.Results <- foreach(Data = Mantel.Data, .combine = rbind) %do% {


      # Compute the test
      Test <- mantel(Data,Zdist, method = "spearman", permutations = 99, na.rm = TRUE)

      # Extract the name of the metric, the significance and the p-value
      Result <- data.frame(Statistic = Test$statistic, Pvalue = Test$signif)
      # Return the results

    } %>% mutate(Metric = names(Mantel.Data), .before = Statistic)

  } # !! End of Mantel Tests

  { # !! Pearson Correlation

  #### --- Computation of the pearson correlation --- # 

  # Pearson correlations are used with alpha metrics against a variable (z)

  Pearson.Correlation <- Data.filtered %>%
    # Filter for "OBS" values
    filter(Type == "OBS") %>%
    group_by(Metric) %>%
    summarise(Statistic = cor(Value, z))

  } # !! End of Pearson Correlation 
  
  ##### --- Combine Alpha and Beta results --- #####

  Correlation.results <- full_join(Mantel.Results,Pearson.Correlation,by = join_by(Metric, Statistic))

  # Transform it into a table to plot
  # Grob the parameters 
  Correlation.results <-tableGrob(Correlation.results)

    # Message
  cat(rule(left = "- END: Plotting - ", line_col = "white", line = " ", col = "green"))

  # Return the plots created
  return(list(Figure_1,Figure_2,Figure_3,Figure_4,Figure_5,Figure_6,Figure_7,Correlation.results))

} # End of range analyses.

#----------------------------------------------------#
##### RANGE-ANALYSES: COMPUTATION OF THE METRICS #####
#----------------------------------------------------#

Moss.results <- range_analyses(Data = mossData, Sp_Names = Moss_Names, Phylo_Tree = Mosses_Tree, Phylo_Tree_Null = Mosses_Tree_Null,  Threshold = 4, Nrep = 99, Alt_limits = c(0,1000,1400,1800,2000,2200,Inf))

# Create the pdf
Moss.pdf <- marrangeGrob(Moss.results, nrow=1, ncol = 1)

# Save the pdf
ggsave(filename = paste0("Exploratory_Metrics_Moss.pdf"),
       plot = Moss.pdf,
       device = "pdf",
       path = paste0(getwd()),
       width = 35,
       height = 25,
       units = "cm")

  # ----------------------- # 

Liver.results <- range_analyses(Data = liverData, Sp_Names = Liver_Names, Phylo_Tree = Liver_Tree,  Threshold = 3, Phylo_Tree_Null = Liver_Tree_Null,Nrep = 3, Alt_limits = c(0,1400,2000,Inf))

# Create the pdf
Liver.pdf <- marrangeGrob(Liver.results, nrow=1, ncol = 1)

# Save the pdf
ggsave(filename = paste0("Exploratory_Metrics_Liver.pdf"),
       plot = Liver.pdf,
       device = "pdf",
       path = paste0(getwd()),
       width = 35,
       height = 35,
       units = "cm")


