###Create model formulas with different predictor combinations fro the nint nurse trait models####
library(tidyverse)

###Function to check if model is valid####
is_valid_model <- function(model) {
  terms <- unlist(strsplit(model, "\\+"))
  
  # Define main effects and their corresponding interaction/squared terms
  interactions <- list("graz" = c("graz:aridity", "graz:RASE", "graz:AMT", "graz:pH", "graz:SAC", 
                                  "graz:log_nurse_meanLA", "graz:log_nurse_meanSLA", "graz:log_nurse_meanH", "graz:log_nurse_meanCNratio"),
                       "aridity" = c("graz:aridity", "RASE:aridity", "AMT:aridity", 
                                     "aridity:log_nurse_meanLA", "aridity:log_nurse_meanSLA", "aridity:log_nurse_meanH", "aridity:log_nurse_meanCNratio"),
                       "RASE" = c("graz:RASE", "RASE:AMT", "RASE:aridity", 
                                  "RASE:log_nurse_meanLA", "RASE:log_nurse_meanSLA", "RASE:log_nurse_meanH", "RASE:log_nurse_meanCNratio"),
                       "AMT" = c("graz:AMT", "RASE:AMT", "AMT:aridity", 
                                 "AMT:log_nurse_meanLA", "AMT:log_nurse_meanSLA", "AMT:log_nurse_meanH", "AMT:log_nurse_meanCNratio"),
                       "pH" = c("graz:pH", "pH:log_nurse_meanLA", "pH:log_nurse_meanSLA", "pH:log_nurse_meanH", "pH:log_nurse_meanCNratio"),
                       "SAC" = c("graz:SAC", "SAC:log_nurse_meanLA", "SAC:log_nurse_meanSLA", "SAC:log_nurse_meanH", "SAC:log_nurse_meanCNratio"), 
                       "log_nurse_meanLA" = c("graz:log_nurse_meanLA", "aridity:log_nurse_meanLA", "AMT:log_nurse_meanLA", 
                                              "RASE:log_nurse_meanLA", "pH:log_nurse_meanLA", "SAC:log_nurse_meanLA"), 
                       "log_nurse_meanSLA" = c("graz:log_nurse_meanSLA", "aridity:log_nurse_meanSLA", "AMT:log_nurse_meanSLA", 
                                               "RASE:log_nurse_meanSLA", "pH:log_nurse_meanSLA", "SAC:log_nurse_meanSLA"), 
                       "log_nurse_meanH" = c("graz:log_nurse_meanH", "aridity:log_nurse_meanH", "AMT:log_nurse_meanH", 
                                             "RASE:log_nurse_meanH", "pH:log_nurse_meanH", "SAC:log_nurse_meanH"),
                       "log_nurse_meanCNratio" = c("graz:log_nurse_meanCNratio", "aridity:log_nurse_meanCNratio", "AMT:log_nurse_meanCNratio", 
                                                   "RASE:log_nurse_meanCNratio", "pH:log_nurse_meanCNratio", "SAC:log_nurse_meanCNratio"))
  
  squared_terms <- list("aridity" = "aridity2", "AMT" = "AMT2")
  
  # Check for interaction terms without main effects
  for (main_effect in names(interactions)) {
    if (any(interactions[[main_effect]] %in% terms) && !(main_effect %in% terms)) {
      return(FALSE)
    }
  }
  
  # Check for squared terms without main effects
  for (main_effect in names(squared_terms)) {
    if ((squared_terms[[main_effect]] %in% terms) && !(main_effect %in% terms)) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}


###Loop that creates the model formulas####
predictors <- c("graz", "aridity", "aridity2", "AMT", "AMT2", "RASE", "pH", "SAC", 
                "log_nurse_meanLA", "log_nurse_meanSLA", "log_nurse_meanH", "log_nurse_meanCNratio",
                "graz:aridity", "graz:RASE", "graz:AMT", #grazing-climate interactions
                "graz:pH", "graz:SAC", #grazing-soil interactions
                "RASE:AMT", "RASE:aridity", "AMT:aridity", #climate-climate interactions
                #trait-environment interactions
                "graz:log_nurse_meanLA", "graz:log_nurse_meanSLA", "graz:log_nurse_meanH", "graz:log_nurse_meanCNratio", 
                "aridity:log_nurse_meanLA", "aridity:log_nurse_meanSLA", "aridity:log_nurse_meanH", "aridity:log_nurse_meanCNratio", 
                "AMT:log_nurse_meanLA", "AMT:log_nurse_meanSLA", "AMT:log_nurse_meanH", "AMT:log_nurse_meanCNratio",
                "RASE:log_nurse_meanLA", "RASE:log_nurse_meanSLA", "RASE:log_nurse_meanH", "RASE:log_nurse_meanCNratio",
                "pH:log_nurse_meanLA", "pH:log_nurse_meanSLA", "pH:log_nurse_meanH", "pH:log_nurse_meanCNratio",
                "SAC:log_nurse_meanLA", "SAC:log_nurse_meanSLA", "SAC:log_nurse_meanH", "SAC:log_nurse_meanCNratio")

#how many combinations are possible?
n_possible_models = 2^length(predictors) -1 #eish 1.759219e+13 = 17 592 190 000 000
#max size of an R df is 2^31 -1 = 2 147 483 648

### Loop that creates the model formulas in chunks ####
chunk_size <- 2^31-1 # Adjust this size based on available memory
output_file <- "Functional trait data\\results\\nint_nurse_trait_clim_soil_formulas.csv"

# Initialize the output file
write.csv(data.frame(formula = character()), output_file, row.names = FALSE)

for (counter1 in 1:length(predictors)) {
  # Create a matrix where each column is a combination of n = counter1 variables
  combos <- as.matrix(combn(predictors, counter1)) 
  n_combos <- ncol(combos)
  
  # Process in chunks
  for (start_idx in seq(1, n_combos, by = chunk_size)) {
    end_idx <- min(start_idx + chunk_size - 1, n_combos)
    chunk <- combos[, c(start_idx:end_idx), drop = FALSE] #isolate a chunk
    
    modlist <- data.frame(formula = character()) #create df to put the formulas in 
    l <- 1
    
    for (counter2 in 1:ncol(chunk)) {
      # Make a formula out of each column in the matrix
      mod <- paste(c(chunk[, counter2]), collapse = "+") 
      #check that the model is valid
      validity <- is_valid_model(mod)
      
      if (validity == TRUE) { # Only add it to modlist if validity is true
        modlist[l, 1] <- mod
        l <- l + 1
      }
    }
    
    if (nrow(modlist) > 0) {
      write.table(modlist, output_file, append = TRUE, sep = ",", row.names = FALSE) #append the new model to the existing file
    }
  }
}
