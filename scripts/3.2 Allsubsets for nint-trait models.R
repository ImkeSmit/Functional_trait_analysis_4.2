###Create model formulas with different predictor combinations fro the nint nurse trait models####
library(tidyverse)

###Get the model formulas####
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
n_possible_models = 2^length(predictors) -1

modlist <- data.frame(formula = character())
l = 1
for(counter1 in 1:length(predictors)) {
  combos <- as.matrix(combn(predictors, counter1))
  
  for(counter2 in 1:ncol(combos)) {
    mod <- paste(c(combos[, counter2]), collapse = "+")
    
    modlist[l, 1] <- mod
    l = l+1
  }}

# Function to check if a model is valid
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

#run modlist through the function to see which models are valid
validity = c()
for(m in 1:nrow(modlist)) {
  validity[m] <- is_valid_model(modlist[m, 1])
}
#subset modlist to keep only valid models
valid_modlist <- data.frame(predictors = modlist[c(which(validity == TRUE)), ])

write.csv(valid_modlist, "Facilitation data\\results\\nint_clim_soil_model_formulas_22Jun2024.csv")
write.csv(formula_table, "Functional trait data\\Results\\nint_nurse_trait_formulas.csv")
