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
# setwd("~/Documents/PhD-Thesis/Research/Switzerland")

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})

# Install/load tons of packages.
p_load("doParallel", # Allow parallel computation
       "ape",        # Multiple tools of phylogenetic analyses
       "ggplot2",    # Graphical representation tools
       "dplyr",      # Load the dplyr syntax
       "stringr",     # String manipulation
       "tidyr",
       "phangorn"
)

# Set a global option to disable dplyr::summarise() warnings when grouping. 
options(dplyr.summarise.inform = FALSE)

# --- Load the data --- #

# Load the tracheophyte data
TracheoOcc  <- read.csv("TracheophyteData.csv",row.names=1)

# Load the tracheophyte tree
# OPTION 1: Baker 2022 / A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life
BakerTree <- read.tree(file = "PhyloTree/Baker2022/treeoflife_current.tree")
# OPTION 2: Zanne 2014 / Three keys to the radiation of angiosperms into freezing environments
ZanneTree <- read.tree(file = "PhyloTree/Zanne2014/Vascular_Plants_rooted.dated.tree")
# OPTION 3: Zuntini 2024 / Phylogenomics and the rise of the angiosperms.
ZuntiniTree <- read.tree(file = "PhyloTree/Zuntini2024/Zuntini2024.tre")[[1]]


#------------------------------#
##### ANGIOSPERM OCCURENCE #####
#------------------------------#

# --- Modify the data --- #

# Split between MetaData and OccurenceData
TracheoMeta <- TracheoOcc[,1:2]       # Verification: table(as.matrix(TracheoOcc))
TracheoOcc <- TracheoOcc[,-c(1:2)]

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
ZuntiniTip <- ZuntiniTree$tip.label

  # -- Baker -- #
  
# Modify the BarkerTips :Only keep the species names (the 3-4 fields)
# Split the data (Stop just after the species names and tranform into a dataframe
BakerTip <- t(as.data.frame(str_split(BakerTip, "_", n = 5))) %>%
  `rownames<-`(NULL) %>%
  `colnames<-`(c("Order","Family","Genus","Specific_epithet","Sequence_ID")) %>%
  as.data.frame() %>%
  unite("Species",Genus:Specific_epithet, remove = F)  # Combine the genus and the specific epithet into a specie name

  # -- Zuntini -- # 

# Remove all the "'"
ZuntiniTip <- gsub("'"," ",ZuntiniTip,fixed = T) 
# Change "_" into " " 
ZuntiniTip <- gsub("_"," ",ZuntiniTip,fixed = T) 

# Modify the ZuntiniTip :Only keep the species names (the 3-4 fields)
# Split the data (Stop just after the species names and tranform into a dataframe
ZuntiniTip <- as.data.frame(str_split_fixed(ZuntiniTip, " ", n = 6)) %>%
  # Only keep the V4 column
  dplyr::select(4:5) %>%
  # remove the ".sp"
  dplyr::filter(V4 != "sp.")  %>%
  unite("Species",V4:V5, remove = F)  # Combine the genus and the specific epithet into a specie name

  # -- Zanne -- # 

# Create a version with only the genus to find the percentage of overlap of the tree genuses with the genus of the species we have.
Zanne_Gn_names <- sub("_.*","",ZanneTip) %>%
  unique() # Keep the genuses names

# Find the species present in the tree
Sp_Found_Zanne <- Sp_names[which(Sp_names %in% ZanneTip)] # %: 67.6  (length(Sp_Found_Zanne)/length(Sp_names))*100 
Sp_Found_Baker <- Sp_names[which(Sp_names %in% BakerTip$Species)]  # %: 14.51  (length(Sp_Found_Baker)/length(Sp_names))*100 
Sp_Found_Zuntini <- Sp_names[which(Sp_names %in% ZuntiniTip$Species)]  # %: 14.51  (length(Sp_Found_Baker)/length(Sp_names))*100 

# Find the Genuses present in the tree
Gn_Found_Zanne <- Gn_names[which(Gn_names %in% Zanne_Gn_names)] # %: 95 (length(Gn_Found_Zanne)/length(Gn_names))*100 
Gn_Found_Baker <- Gn_names[which(Gn_names %in% BakerTip$Genus)] # %: 83 (length(Gn_Found_Baker)/length(Gn_names))*100 

  # ! WE SELECT THE ZANNE TREE ! # 

# ------------------------------------------------------------------------------------------------- #

#----------------------------------#
##### PHYLOTREE MODIFICIATIONS #####
#----------------------------------#

  # We will have to change the tips of the liverworts and mosses to match with the new databases of occurence (BryophytesData OR FinalMatrixLv95_SynChangedV3)
  # See treemodif for modification done. 

# Load the names of the mosses present in the dataframe of Occ_Data_Moss for further splitting between Mosses and Liverworts.
Bryo_Names <- read.csv(file = "Utilities/Bryophyte_Names.csv", row.names = 1)
Moss_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Moss"),Species)$Species
Liver_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Liverwort"),Species)$Species

# Load the data for further testing 

  # BryophyteData is the dataframe in common with the tracheophytes for the Lee Analyses.
  # FinalMatrix is the global bryophyte dataset with all our plots. 

bryoData <- read.csv("BryophyteData.csv",row.names=1)
  # OR #
bryoData <- read.csv("FinalMatrixLv95_SynChangedV3.csv",row.names=1)

# Split the Occ_Data_Moss between Mosses and Liverworts
liverData <- dplyr::select(bryoData,any_of(Liver_Names))
mossData <- dplyr::select(bryoData,any_of(Moss_Names))

# Get the names of the liverworts and mosses present
Moss_Names <- colnames(mossData)
Liver_Names <- colnames(liverData)

# Load the phylotrees
Liver_Tree <- read.tree("PhyloTree/timetree50mod-liverworts.nwk")
Mosses_Tree <- read.tree("PhyloTree/timetree50mod-mosses.nwk")

# Get the tips clean
Liver_Tree$tip.label <- Liver_Tree$tip.label %>% gsub("'","",.)
Mosses_Tree$tip.label<- Mosses_Tree$tip.label %>% gsub("'","",.)  

  # ----- LIVERWORTS ----- # 

# Tips to delete (They are deleted because they will not be find in the DB of occurence, therefore it is easier to verify that all tips are found afterwards)
Tips_deletion <- c("Marchantia_alpestris","Leiosporoceros_dussii","Anthoceros_agrestis","Scapania_cf_helvetica")

# Tips to create a polytomy with the tips "Jungermannia"
Tips_Poly_Jungermannia <- c("Jungermannia_atrovirens","Jungermannia_borealis","Jungermannia_cf_obovata","Jungermannia_cf_polaris")

# Tips to change the name from OLD to NEW
Tips_OLD <- c("Lophozia_bicrenata","Lophozia_sudetica","Barbilophizoa_lycopodioides","Calypogeia_fissa","Jungermannia_gracillima","Leiocolea_bantriensis","Leiocolea_heterocolpos","Pellia_cf_epiphylla","Preissia_quadrata","Reboulia_hemispherica","Mannia_cf_controversa")
Tips_NEW <- c("Isopaches_bicrenatus","Barbilophozia_sudetica","Barbilophozia_lycopodioides","Calypogeia_azurea","Solenostoma_gracillimum","Mesoptychia_bantriensis","Mesoptychia_heterocolpos","Pellia_epiphylla","Marchantia_quadrata","Reboulia_hemisphaerica","Mannia_cf_androgyna")

  # -- Transform the tree -- #

New_Liver_Tree <- Liver_Tree %>%
  # Delete the tips
  drop.tip(Tips_deletion) %>%
  # Create the polytomy and remove the old
  add.tips(Tips_Poly_Jungermannia,"Jungermannia") %>%
  drop.tip("Jungermannia")
  
# Change the tips with the new name
New_Liver_Tree$tip.label[match(Tips_OLD, New_Liver_Tree$tip.label)] <- Tips_NEW

# Get the new tips
Liver_Tips <- New_Liver_Tree$tip.label

  # -- Verification -- #

# Found species from the community present in the tree. 
Int_Liver_Tips <-  base::intersect(Liver_Names,Liver_Tips)   

# Found species from the community not present in the tree
Miss_Liver_Tips <- Liver_Names[!Liver_Names %in% intersect(Liver_Tips,Liver_Names)] 

# Found Species present in the tree but not present in the communities
Liver_Missing <- Liver_Tips[!Liver_Tips %in% intersect(Liver_Tips,Liver_Names)]

  # -- Saving -- # 

# Save the new tree
write.tree(phy = New_Liver_Tree, file = "PhyloTree/timetree50mod-liverwortsV2.nwk")


# ----- MOSSES ----- # 

# Tips to delete (They are deleted because they will not be find in the DB of occurence, therefore it is easier to verify that all tips are found afterwards)
Tips_deletion <- c("Calliergon_cordifolium","Schistidium_papillosum","Schistidium_helveticum","Schistidium_elegantulum_subsp_wilsonii","Schistidium_grande","Schistidium_brunnescens_subsp_griseum","Schistidium_brunnescens_subsp_brunnescens","Schistidium_memnonium","Takakia_lepidozioides")

# Tips to create a polytomy
Tips_Poly_P_falcata <- c("Palustriella_falcata_var._sulcata","Palustriella.falcata_s.str.") # Palustriella_falcata
Tips_Poly_R_beskeanum <- c("Rhynchostegium_murale","Rhynchostegium_megapolitanum") # Rhynchostegium_beskeanum
Tips_Poly_S_lescurii <- c("Sphagnum_capillifolium","Sphagnum_fimbriatum","Sphagnum_subsecundum") # Sphagnum_lescurii

# Tips to change the name from OLD to NEW
Tips_OLD <- c("Hypnum_callichroum",
              "Hygrohypnum_eugyrium",
              "Hypnum_procerrimum",
              "Thuidium_recognitum",
              "Campylium_protensum",
              "Drepanocladus_longifolius",
              "Campylium_chrysophyllum",
              "Campylophyllum_calcareum",
              "Homomallium_adnatum",
              "Rhytidiadelphus_triquetrus",
              "Brachythecium_cirrosum",
              "Ditrichum_flexicaule",
              "Rhodobryum_ontariense",
              "Plagiobryum_zierii",
              "Bryum_elegans",
              "Plagiomnium_elatum",
              "Cyrtomium_hymenophylloides",
              "Mnium_spinulosum",
              "Bartramia",
              "Tayloria_serrata_agg",
              "Dicranoweisia_crispula",
              "Phascum_cupsidatum",
              "Pottia_truncata",
              "Schistidium_elegantulum_subsp_elegantulum")
              
Tips_NEW <- c("Stereodon_callichrous",
              "Hygrohypnum_luridum",
              "Pseudostereodon_procerrimus",
              "Thuidium_assimile",
              "Campylium_stellatum",
              "Drepanocladus_aduncus",
              "Campyliadelphus_chysophyllus",
              "Campylophyllopsis_calcarea",
              "Homallium_incurvatum",
              "Hylocomiadelphus_triquetrus",
              "Brachythecium_cirrhosum",
              "Flexitrichum_flexicaule",
              "Rhodobryum_roseum",
              "Ptychostomum_zieri",
              "Ptychostomum_elegans",
              "Plagiomnium_ellipticum",
              "Cyrtomnium_hymenophylloides",
              "Mnium_spinosum.agg.",
              "Bartramia_ithiphylla",
              "Tayloria_serrata",
              "Hymenoloma_crispulum",
              "Tortula_acaulon",
              "Tortula_truncata",
              "Schistidium_elegantulum")

# -- Transform the tree -- #

New_Moss_Tree <- Mosses_Tree %>%
  # Delete the tips
  drop.tip(Tips_deletion) %>%
  # Create the polytomy and remove the old
  add.tips(Tips_Poly_P_falcata,"Palustriella_falcata") %>% 
    drop.tip("Palustriella_falcata") %>%
  add.tips(Tips_Poly_R_beskeanum,"Rhynchostegium_beskeanum") %>%
    drop.tip("Rhynchostegium_beskeanum") %>%
  add.tips(Tips_Poly_S_lescurii,"Sphagnum_lescurii") %>%
    drop.tip("Sphagnum_lescurii") %>%
  # Add "Schistidium_brunnescens_subsp_brunnescens_alpine_morph" to "Schistidium_brunnescens"
  add.tips("Schistidium_brunnescens_subsp_brunnescens_alpine_morph","Schistidium_brunnescens") %>%
    drop.tip("Schistidium_brunnescens_subsp_brunnescens_alpine_morph")

# Change the tips with the new name
New_Moss_Tree$tip.label[match(Tips_OLD, New_Moss_Tree$tip.label)] <- Tips_NEW

# Get the new tips
Moss_Tips <- New_Moss_Tree$tip.label

# -- Verification -- #

# Found species from the community present in the tree. 
Int_Moss_Tips <-  base::intersect(Moss_Names,Moss_Tips)   

# Found species from the community not present in the tree
Miss_Moss_Tips <- Moss_Names[!Moss_Names %in% intersect(Moss_Tips,Moss_Names)] 

# Found Species present in the tree but not present in the communities
Moss_Missing <- Moss_Tips[!Moss_Tips %in% intersect(Moss_Tips,Moss_Names)]

# -- Saving -- # 

# Save the new tree
write.tree(phy = New_Moss_Tree, file = "PhyloTree/timetree50mod-mossesV2.nwk")
