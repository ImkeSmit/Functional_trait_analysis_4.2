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

full_formula2 <- as.formula("NIntc_richness_binom ~ 
                            graz*aridity + graz*RASE + graz*AMT +
                            graz*pH + graz*SAC +
                            aridity*RASE + aridity*AMT +
                            graz*log_nurse_meanH + graz*log_nurse_meanLDMC +
                            aridity*log_nurse_meanH + aridity*log_nurse_meanLDMC +
                            AMT*log_nurse_meanH + AMT*log_nurse_meanLDMC +
                            RASE*log_nurse_meanH + RASE*log_nurse_meanLDMC +
                            #pH*log_nurse_meanH + pH*log_nurse_meanLDMC +
                            #SAC*log_nurse_meanH + SAC*log_nurse_meanLDMC +
                            Lat_decimal +Long_decimal + #add these to account for spatial structure instead of (1|site_ID/ID)
                            (1|nurse_sp)")



# Fit the full model
options(na.action = "na.omit")
full_model <- glmmTMB(
  formula = full_formula2,
  data = modeldat_final,    
  family = binomial  #had to remove sq terms and soil:climate and soil:trait interactions to make model converge
)

test_model <- glmmTMB(NIntc_richness_binom ~ aridity*log_nurse_meanH +  aridity*log_nurse_meanLDMC + 
                        aridity2 + Lat_decimal + Long_decimal + (1|nurse_sp), data = modeldat_final, family = binomial)

# Ensure all models maintain random effects by excluding them from being dropped
options(na.action = "na.fail") # Required for dredge function

# Perform stepwise model selection using dredge
model_selection <- dredge(
  full_model,
  fixed = c("cond(Lat_decimal)","cond(Long_decimal)"), #random effects are automatically included in all models due to the structure of tMB
  rank = "AIC"                                # Use AIC for model ranking
) 
#maybe we can find a way to exclude the model with only lat and long as predictors in the dredge output. 
#but maybe that is not correct?

#save model selection results
write.csv(as.data.frame(model_selection), row.names = FALSE ,"Functional trait data\\paper results\\stepwise_results.csv")

# View model selection table
print(model_selection)
#in output: cond(Int) refers to the intercept of the conditional model (test_model in this case).
#disp(Int) is the intercept of the dispersion model. Since you did not specify a dispersion structure in TMB, this is default settings. do not worry about this column
#if a column has a number or + in the output, it means the variable was included in that model. 


# Extract the best model based on AICc
best_model <- get.models(model_selection, subset = 1)[[1]] #subset = 1 gets the model with the absolute lowest AIC.

#extract models with AIC difference less than 2
eq_model <- get.models(model_selection, subset = delta < 2)
#average these models
avg_models <- model.avg(eq_model)
#from vignette:
#The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. 
#An alternative, the ‘full’ average assumes that a variable is included in every model, but in some models 
#the corresponding coefficient (and its respective variance) is set to zero. Unlike the ‘subset average’, 
#it does not have a tendency of biasing the value away from zero. The ‘full’ average is a type of shrinkage 
#estimator, and for variables with a weak relationship to the response it is smaller than ‘subset’ estimators.
formula(avg_models) #get formula of the averaged model
summary(avg_models) #get summary of the averaged model
