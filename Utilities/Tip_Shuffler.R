#--------------------------------------#
##### CREATION OF NULL MODEL TREES #####
#--------------------------------------#

    # The goal of this function is to create a number of Null Model trees by tip-shuffling. 
    # It will create a multiphylo object contaiing all the created trees.
    # This function will work on parallel cores.

# - ARGUMENTS - #

    # Tree: A Tree of class "phylo".
    # NM: Number of Null-Model trees to be created.
    # Ncore: Number of cores to use. (Default: Half of the cores present in the computer). 


# - TEST - #

# Tree <- read.tree("PhyloTree/timetree50mod-mossesV2.nwk")
# NM <- 10
# Test <- tip.shuffler(Tree = Tree, NM = 5, Ncore =2)

# - FUNCTION - #

tip.shuffler <- function(Tree, NM, Ncore = detectCores() / 2) {
    # - Make test - #
    if (round(NM) != NM) {
        stop("The number of null model tree to realize is not an integer, please enter an integer as the number of replicates to perform.")
    }
    if (class(Tree) != "phylo") {
        stop("The tree file entered is not of class \"phylo\".")
    }

    # - Load the needed packages - #
    if (!require(foreach)) {
        install.packages("foreach")
        require(foreach)
    }

    # - Computation of the null model trees - #
    Phylo_NULL <- foreach(j = 1:NM) %dopar% {
        # Get a vector of random names from the species names
        Random_Name <- sample(Tree$tip.label, size = length(Tree$tip.label))
        # Get the good tree topology
        Tree_NULL_MODEL <- Tree
        # Create the RANDOM tree (by putting back the vector of species names shuffled)
        Tree_NULL_MODEL$tip.label <- Random_Name
        # Return the RANDOM tree.
        return(Tree_NULL_MODEL)
    }

    # Transform the list into a multiphylo_Object
    class(Phylo_NULL) <- "multiPhylo"
    # Rename it
    names(Phylo_NULL) <- paste0("Tree_NULL_MODEL_", 1:NM)
    # Return it
    return(Phylo_NULL)
}

