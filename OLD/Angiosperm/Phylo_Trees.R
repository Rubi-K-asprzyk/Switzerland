#!/usr/bin/env Rscript

# ---------- #
# K.Thibault #
# ---------- #

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT PHYLO_TREES BEGINNING", line_col = "red", line = "-", col = "br_red")) 
cat(rule(left = "INITIALISATION", line_col = "green", line = "-", col = "br_green"))

#------------------------#
##### INITIALISATION #####
#------------------------#

# Set the Working directory
setwd("~/Documents/PhD-Thesis/Research/Switzerland")

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})

# Install/load tons of packages.
p_load("doParallel", # Allow parallel computation
       "ape",        # Multiple tools of phylogenetic analyses
       "ggplot2",    # Graphical representation tools
       "dplyr",      # Load the dplyr syntax
       "stringr",     # String manipulation
       "tidyr"
)

# Set a global option to disable dplyr::summarise() warnings when grouping. 
options(dplyr.summarise.inform = FALSE)

# --- Load the data --- #

# Load the tracheophyte tree
  # OPTION 1: Baker 2022 / A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life
BakerTree <- read.tree(file = "Data/Angiosperm/PhyloTree/Baker2022/treeoflife_current.tree")
  # OPTION 2: Zanne 2014 / Three keys to the radiation of angiosperms into freezing environments
ZanneTree <- read.tree(file = "Data/Angiosperm/PhyloTree/Zanne2014/Vascular_Plants_rooted.dated.tree")

# Load the TracheoOccurences
TracheoOcc<- read.csv(file = "Data/Angiosperm/Occurence_Tracheo.csv", row.names = 1)

#------------------------------#
##### ANGIOSPERM OCCURENCE #####
#------------------------------#

# --- Modify the data --- #

# Replace column Abies.Alba by Abies.Alba1 and remove Abies.Alba1
TracheoOcc$Abies.alba <- TracheoOcc$Abies.alba.1
TracheoOcc <- select(TracheoOcc, -c(Abies.alba.1))

# Split between MetaData and OccurenceData
TracheoMeta <- TracheoOcc[,1:7]       # Verification: table(as.matrix(TracheoOcc))
TracheoOcc <- TracheoOcc[,-c(1:7)]

# Get the species names
Sp_names <- colnames(TracheoOcc)

# Replace "." by "_"
Sp_names <- gsub(".","_",Sp_names,fixed = T)
# Keep only the genuses names
Gn_names <- sub("_.*","",Sp_names) %>%
  unique() # Keep the genuses names

# --- Find the percentage of overlap between the phylotree and the species we have.

# Extract the tips from the tree
BakerTip <- BakerTree$tip.label
ZanneTip <- ZanneTree$tip.label

# Modify the BarkerTips :Only keep the species names (the 3-4 fields)
  # Split the data (Stop just after the species names and tranform into a dataframe
BakerTip <- t(as.data.frame(str_split(BakerTip, "_", n = 5))) %>%
  `rownames<-`(NULL) %>%
  `colnames<-`(c("Order","Family","Genus","Specific_epithet","Sequence_ID")) %>%
  as.data.frame() %>%
  unite("Species",Genus:Specific_epithet, remove = F)  # Combine the genus and the specific epithet into a specie name

# Create a version with only the genus to find the percentage of overlap of the tree genuses with the genus of the species we have.
Zanne_Gn_names <- sub("_.*","",ZanneTip) %>%
  unique() # Keep the genuses names

# Find the species present in the tree
Sp_Found_Zanne <- Sp_names[which(Sp_names %in% ZanneTip)] # %: 67.9  (length(Sp_Found_Zanne)/length(Sp_names))*100 
Sp_Found_Baker <- Sp_names[which(Sp_names %in% BakerTip$Species)] #  %: 13.96  (length(Sp_Found_Baker)/length(Sp_names))*100 

# Find the Genuses present in the tree
Gn_Found_Zanne <- Gn_names[which(Gn_names %in% Zanne_Gn_names)] # %: 95 (length(Gn_Found_Zanne)/length(Gn_names))*100 
Gn_Found_Baker <- Gn_names[which(Gn_names %in% BakerTip$Genus)] # %: 83 (length(Gn_Found_Baker)/length(Gn_names))*100 





















