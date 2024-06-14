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
setwd("~/Documents/PhD-Thesis/Research/Switzerland/Data")

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})

# Install/load tons of packages.
p_load("doParallel", # Allow parallel computation
       "ape",        # Multiple tools of phylogenetic analyses
       "ggplot2",    # Graphical representation tools
       "dplyr",      # Load the dplyr syntax
       "stringr",     # String manipulation
       "tidyr",
       "gdm",
       "scales",
       "FactoMineR",
       "factoextra",
       "broom",
       "spaa"
)

# Set a global option to disable dplyr::summarise() warnings when grouping. 
options(dplyr.summarise.inform = FALSE)

# --- Load the data --- #

# Load the environmental data
Env_Data <- read.csv(file = "Env_Var_Bryo.csv", row.names = 1)

# Load the occurrence data
Occ_Data <- read.csv(file = "Occurence_Bryo.csv", row.names = 1)

#-----------------------------------#
##### MAKE ENVIRONNEMTENTAL ACP #####
#-----------------------------------#

# Make a copy of the environnemental data
ACP_Data <- Env_Data

# Split the data between the different "types" of variables (the +2 allows to avoid the two first columns that are x and y coordinates)
ACP_Edaphic <- ACP_Data[,c(1:3,6:8)+2]
ACP_Clim <- ACP_Data[,c(4,21:44,68:69)+2]
ACP_Topo <- ACP_Data[,c(5,45:53,67)+2]
ACP_LC <- ACP_Data[,c(9:20,54:66)+2]

# Make the ACPs
ACP_Edaphic <- PCA(ACP_Edaphic, scale.unit = TRUE, ncp = 5, graph = TRUE)  # Auto correlated 
ACP_Clim <- PCA(ACP_Clim, scale.unit = TRUE, ncp = 5, graph = TRUE) # Auto correlated
ACP_Topo <- PCA(ACP_Topo, scale.unit = TRUE, ncp = 5, graph = TRUE) # Not auto correlated
ACP_LC <- PCA(ACP_LC, scale.unit = TRUE, ncp = 5, graph = TRUE) # Not auto correlated

# Create the dataset to make the analysis on. 
ACP_Data <- data.frame("Edaphic" = ACP_Edaphic$ind$coord[,1],
                       "Climatic" = ACP_Clim$ind$coord[,1],
                       "Topo" = ACP_Topo$ind$coord[,1],
                       "LC" = ACP_LC$ind$coord[,1])

# Add the coordinates
ACP_Data <- cbind(Env_Data[,1:2],ACP_Data)

#----------------------#
##### MAKE THE GDM #####
#----------------------#

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

  # preData : A wite by predictor matric containing all the GDM predictors in column. 

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



