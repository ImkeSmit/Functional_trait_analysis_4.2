library(MuMIn)
library(glmmTMB)
library(tidyverse)
library(tidylog)

modeldat_final <- read.csv("Functional trait data\\Clean data\\nint_nurse_traits.csv", row.names = 1) |> 
  mutate(aridity2 = aridity^2) |> 
  filter(!is.na(NIntc_richness_binom)) |> 
  filter(!is.na(NIntc_cover_binom)) |> 
  filter(!is.na(NInta_richness_binom)) |> 
  filter(!is.na(NInta_cover_binom))
modeldat_final$nurse_sp <- as.factor(modeldat_final$nurse_sp)
modeldat_final$graz <- as.factor(modeldat_final$graz)
modeldat_final$site_ID <- as.factor(modeldat_final$site_ID)
modeldat_final$ID <- as.factor(modeldat_final$ID)

##Add the other environmental covariates to modeldat final
#import siteinfo, we will use this to add ID to drypop
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  mutate(plotref = str_c(SITE, PLOT, sep = "_")) |> 
  select(ID, plotref) |> 
  distinct() |> 
  na.omit()

#import drypop, so which contains the env covariates
drypop <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Functional trait data\\Raw data\\drypop_20MAy.csv") |> 
  mutate(plotref = str_c(Site, Plot, sep = "_")) |> #create a variable to identify each plot
  select(plotref, AMT, RAI, RASE, pH.b, SAC.b) |> 
  distinct() |> 
  left_join(siteinfo, by = "plotref") |> 
  select(!plotref)
drypop$ID <- as.factor(drypop$ID)

#join the env covariates to the nurse nint data
modeldat_final <- modeldat_final |> 
  inner_join(drypop, by = "ID") |> 
  rename(pH = "pH.b", SAC = "SAC.b") |> 
  mutate(AMT2 = AMT^2)

# Define the full model formula
full_formula <- as.formula("NIntc_richness_binom ~ 
    graz + aridity + aridity2 + AMT + AMT2 + RASE + pH + SAC +
    log_nurse_meanLA + log_nurse_meanSLA + log_nurse_meanH + log_nurse_meanCNratio +
    graz:aridity + graz:RASE + graz:AMT + graz:pH + graz:SAC +
    RASE:AMT + RASE:aridity + AMT:aridity") #+
    #graz:log_nurse_meanLA + graz:log_nurse_meanSLA + graz:log_nurse_meanH + graz:log_nurse_meanCNratio +
    #aridity:log_nurse_meanLA + aridity:log_nurse_meanSLA + aridity:log_nurse_meanH + aridity:log_nurse_meanCNratio +
    #AMT:log_nurse_meanLA + AMT:log_nurse_meanSLA + AMT:log_nurse_meanH + AMT:log_nurse_meanCNratio +
    #RASE:log_nurse_meanLA + RASE:log_nurse_meanSLA + RASE:log_nurse_meanH + RASE:log_nurse_meanCNratio +
    #(1|nurse_sp) + (1|site_ID/ID)") ##haha obviously does not converge
#pH:log_nurse_meanLA + pH:log_nurse_meanSLA + pH:log_nurse_meanH + pH:log_nurse_meanCNratio +
    #SAC:log_nurse_meanLA + SAC:log_nurse_meanSLA + SAC:log_nurse_meanH + SAC:log_nurse_meanCNratio +

options(na.action = "na.omit") 

# Fit the full model
full_model <- glmmTMB(
  formula = full_formula,
  data = modeldat_final,    # Replace 'your_data' with your dataset name
  family = binomial
)

# Ensure all models maintain random effects by excluding them from being dropped
options(na.action = "na.fail") # Required for dredge function

# Perform stepwise model selection using dredge
model_selection <- dredge(
  full_model, 
  fixed = c("(1|nurse_sp)", "(1|site_ID/ID)"),  # Keep random effects in all models
  rank = "AIC"                                # Use AIC for model ranking
) #took 50 mins to run

# View model selection table
print(model_selection)

# Extract the best model based on AICc
best_model <- get.models(model_selection, subset = 1)[[1]] #subset = 1 gets the model with the absolute lowest AIC.

#extract models with AIC difference less than 2
eq_model <- get.models(model_selection, subset = delta < 2)
#average these models
avg_models <- model.avg(eq_model)


# Display summary of the best model
summary(best_model)

# Save model selection results to a CSV file
write.csv(as.data.frame(model_selection), "mumin_model_selection_results.csv", row.names = FALSE)

# Optional: Refit the best model if further validation is required
refit_best_model <- update(best_model, data = your_data)
summary(refit_best_model)
