###Script to run stepwise model selection procedure to determine how traits, the environment, 
###and their interaction affect interaction outcomes
library(glmmTMB)
library(MuMIn)
library(tidyverse)

#full model: nint ~ nurse_maxh + nurse_ldmc + graz + climate + soil
                  # + graz:climate + graz:soil + soil:climate + climate:climate +
                  # + maxh:climate + maxh:soil + maxh:graz
                  #ldmc:climate + ldmc:soil + ldmc:graz

modeldat <- read.csv("Functional trait data\\Clean data\\nint_nurse_traits.csv", row.names = 1) 
modeldat$nurse_sp <- as.factor(modeldat$nurse_sp)
modeldat$graz <- as.factor(modeldat$graz)
modeldat$site_ID <- as.factor(modeldat$site_ID)
modeldat$ID <- as.factor(modeldat$ID)

##Add the other environmental covariates to modeldat final
#import siteinfo, we will use this to add ID to drypop
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  mutate(plotref = str_c(SITE, PLOT, sep = "_")) |> 
  select(ID, plotref, Lat_decimal, Long_decimal) |> 
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
modeldat_final <- modeldat |> 
  inner_join(drypop, by = "ID") |> 
  rename(pH = "pH.b", SAC = "SAC.b") |> 
  mutate(AMT2 = AMT^2, 
         aridity2 = aridity^2) |> 
  #remove all rows which have NA values in any of our modelling variables
  drop_na(NIntc_richness_binom, NIntc_cover_binom,NInta_richness_binom, NInta_cover_binom, 
          log_nurse_meanLDMC, log_nurse_meanH, aridity, AMT, RASE, SAC, pH, graz)

####define the full model formula ####

full_formula <- as.formula("NIntc_richness_binom ~ 
    graz + aridity + AMT + RASE + pH + SAC + log_nurse_meanH + log_nurse_meanLDMC +
    graz:aridity + graz:RASE + graz:AMT + graz:pH + graz:SAC +
    #pH:aridity + pH:AMT + pH:RASE + SAC:aridity + SAC:AMT + SAC:RASE +
    RASE:AMT + RASE:aridity + AMT:aridity +
graz:log_nurse_meanH + graz:log_nurse_meanLDMC +
aridity:log_nurse_meanH + aridity:log_nurse_meanLDMC +
AMT:log_nurse_meanH + AMT:log_nurse_meanLDMC +
RASE:log_nurse_meanH + RASE:log_nurse_meanLDMC +
pH:log_nurse_meanH + pH:log_nurse_meanLDMC +
SAC:log_nurse_meanH + SAC:log_nurse_meanLDMC +
Lat_decimal +Long_decimal + #add these to account for spatial structure instead of (1|site_ID/ID)
(1|nurse_sp)") 



# Fit the full model
options(na.action = "na.omit")
full_model <- glmmTMB(
  formula = full_formula,
  data = modeldat_final,    
  family = binomial  #had to remove sq terms and soil:climate interactions to make model converge
)

# Ensure all models maintain random effects by excluding them from being dropped
options(na.action = "na.fail") # Required for dredge function

# Perform stepwise model selection using dredge
model_selection <- dredge(
  full_model, 
  fixed = c("(1|nurse_sp)", "Long_decimal", "Lat_decimal"),  # Keep random effects in all models
  rank = "AIC"                                # Use AIC for model ranking
) #took 50 mins to run

write.csv(model_selection, "Functional trait data\\paper_results\\stepwise_results.csv")

# View model selection table
print(model_selection)


