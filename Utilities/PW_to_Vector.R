##### ----- TEST ----- #####
# PW_matrix <- Metric_List[[1]] ; Type <- "Intra" # Intra
# PW_matrix <- Metric_List[[6]] ; Type <- "Inter" # Inter
# 
# Test <- PW_to_Vector(Intra, Colname = "Test")

##### ----- TEST ----- ##### 

PW_to_Vector <- function(PW_matrix, Colname, Diag = T){

  # Do we want to keep the intra plot values i.e the diagonal of the matrix ? 
  if (Diag == T) {Compound_new <- upper.tri(PW_matrix , diag = TRUE)}   # YES
  if (Diag == F) {Compound_new <- upper.tri(PW_matrix , diag = FALSE)}  # NO 
  
  # Transform the pairwise matrix into a vector 
  i2 <- which(Compound_new, arr.ind=TRUE)
    Compound_final <- data.frame(Sample = paste0(rownames(PW_matrix)[i2[,1]],"-",colnames(PW_matrix)[i2[,2]]),
                                 PlotA = rownames(PW_matrix)[i2[,1]], 
                                 PlotB = colnames(PW_matrix)[i2[,2]], 
                                 Value = PW_matrix[Compound_new])
  
  # Rename the value column with the name of the metric 
  colnames(Compound_final)[4] <- Colname
    
  # Return the new dataframe
  return(Compound_final)
  
}
