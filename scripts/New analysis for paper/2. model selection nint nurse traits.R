###Script to run model selection procedure to determine how traits, the environment, 
###and their interaction affect interaction outcomes
library(glmmTMB)
library(MuMIn)
library(MASS)
library(tidyverse)
library(parallel)

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
  dplyr::select(ID, plotref, Lat_decimal, Long_decimal) |> 
  distinct() |> 
  na.omit()

#import drypop, so which contains the env covariates
drypop <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Functional trait data\\Raw data\\drypop_20MAy.csv") |> 
  mutate(plotref = str_c(Site, Plot, sep = "_")) |> #create a variable to identify each plot
  dplyr::select(plotref, AMT, RAI, RASE, pH.b, SAC.b) |> 
  distinct() |> 
  left_join(siteinfo, by = "plotref") |> 
  dplyr::select(!plotref)
drypop$ID <- as.factor(drypop$ID)

#join the env covariates to the nurse nint data
modeldat_final <- modeldat |> 
  inner_join(drypop, by = "ID") |> 
  rename(pH = "pH.b", SAC = "SAC.b") |> 
  mutate(AMT2 = AMT^2, 
         aridity2 = aridity^2, 
         sin_lat = sin(Lat_decimal), #transform from a circular to a linear variable
         sin_long = sin(Long_decimal)) |> 
  #remove all rows which have NA values in any of our modelling variables
  drop_na(NIntc_richness_binom, NIntc_cover_binom,NInta_richness_binom, NInta_cover_binom, 
          log_nurse_meanLDMC, log_nurse_meanH, aridity, AMT, RASE, SAC, pH, graz)

#how many replicates included?
n_reps <- modeldat_final |> 
  distinct(ID, replicate_no) |> 
  summarise(nreps = n())
percent_reps_included <- n_reps$nreps/3735 *100 #66.96, we get the 3735 from the nint_trait_models script

#how many dominant sp
n_doms <- modeldat_final |> 
  distinct(nurse_sp) |> 
  summarise(n_nurses = n())


####define the full model formula ###
#formula without climate*climate interactions
#with trait*soil interactions
#This model converges!
full_formula3 <- as.formula("NIntc_richness_binom ~ 
                            graz*aridity + graz*RASE + graz*AMT +
                            graz*pH + graz*SAC +
                            graz*log_nurse_meanH + graz*log_nurse_meanLDMC +
                            aridity*log_nurse_meanH + aridity*log_nurse_meanLDMC +
                            AMT*log_nurse_meanH + AMT*log_nurse_meanLDMC +
                            RASE*log_nurse_meanH + RASE*log_nurse_meanLDMC +
                            pH*log_nurse_meanH + pH*log_nurse_meanLDMC +
                            SAC*log_nurse_meanH + SAC*log_nurse_meanLDMC +
                            sin_lat + sin_long + #add these to account for spatial structure instead of (1|site_ID/ID)
                            (1|nurse_sp)")

####MODEL SELECTION FOR NINTC RICHNESS####

# Fit the full model
options(na.action = "na.omit")
full_model <- glmmTMB(
  formula = full_formula3,
  data = modeldat_final,    
  family = binomial  #had to remove sq terms and soil:climate and soil:trait interactions to make model converge
)


# Ensure all models maintain random effects by excluding them from being dropped
options(na.action = "na.fail") # Required for dredge function

#create a cluster obeject to run the function over 8 cores
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK" 
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))

# Export necessary objects and functions to the cluster
clusterExport(clust, varlist = c("modeldat_final", "full_formula3"), envir = environment())
clusterEvalQ(clust, library(glmmTMB))
clusterEvalQ(clust, library(MuMIn))

#peform model selection using 8 cores
# Perform model selection using dredge
model_selection_par <- dredge(
  full_model,
  fixed = c("cond(sin_lat)","cond(sin_long)"), #random effects are automatically included in all models due to the structure of tMB
  rank = "AIC", # Use AIC for model ranking
cluster = clust) #start Sunday 16:33, end Tuesday 6:30

# Stop the cluster after use
stopCluster(clust)

#save model selection results
saveRDS(model_selection_par, "Functional trait data\\paper results\\nint_richness_nurse_trait_dredge_result.rds")

# View model selection table
model_selection_par <- readRDS("Functional trait data\\paper results\\nint_richness_nurse_trait_dredge_result.rds")
print(model_selection_par)
#in output: cond(Int) refers to the intercept of the conditional model (test_model in this case).
#disp(Int) is the intercept of the dispersion model. Since you did not specify a dispersion structure in TMB, this is default settings. do not worry about this column
#if a column has a number or + in the output, it means the variable was included in that model. 


# Extract the best model based on AICc
best_model <- get.models(model_selection_par, subset = 1)[[1]] #subset = 1 gets the model with the absolute lowest AIC.


#extract models with AIC difference less than 2
eq_model <- get.models(
  model_selection_par, 
  subset = delta < 2 )

#average these models
avg_models <- model.avg(eq_model)
#from vignette:
#The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. 
#An alternative, the ‘full’ average assumes that a variable is included in every model, but in some models 
#the corresponding coefficient (and its respective variance) is set to zero. Unlike the ‘subset average’, 
#it does not have a tendency of biasing the value away from zero. The ‘full’ average is a type of shrinkage 
#estimator, and for variables with a weak relationship to the response it is smaller than ‘subset’ estimators.

avg_model_formula <-formula(avg_models) #get formula of the averaged model
#NIntc_richness_binom ~ 0 + AMT + aridity + graz + log_nurse_meanH + 
#log_nurse_meanLDMC + pH + RASE + SAC + sin_lat + sin_long + 
#  AMT:log_nurse_meanH + AMT:log_nurse_meanLDMC + graz:log_nurse_meanLDMC + 
#  graz:RASE + graz:SAC + log_nurse_meanH:pH + log_nurse_meanH:SAC + 
#  log_nurse_meanLDMC:pH + log_nurse_meanLDMC:RASE + log_nurse_meanLDMC:SAC

summary(avg_models) #get summary of the averaged model

options(na.action = "na.omit")
final_model <-glmmTMB(NIntc_richness_binom ~ AMT + aridity + graz + log_nurse_meanH + 
                        log_nurse_meanLDMC + pH + RASE + SAC + sin_lat + sin_long + 
                        AMT:log_nurse_meanH + AMT:log_nurse_meanLDMC + graz:log_nurse_meanLDMC + 
                        graz:RASE + graz:SAC + log_nurse_meanH:pH + log_nurse_meanH:SAC + 
                        log_nurse_meanLDMC:pH + log_nurse_meanLDMC:RASE + log_nurse_meanLDMC:SAC+
                        sin_lat + sin_long + (1|nurse_sp), data = modeldat_final, family = "binomial", weights = NULL)
summary(final_model)
r.squaredGLMM(final_model) #r squared is not working

# Alternative method because Mumin doesn't work
library(performance)
performance::r2(final_model)

#get the variable importance
importance <- sw(avg_models)
#save to csv file
write.csv(as.data.frame(importance), "Functional trait data\\paper results\\nint_rich_model_weights.csv")


####MODEL SELECTION FOR NINTC COVER####
####define the full model formula ####
#formula without climate*climate interactions
#with trait*soil interactions
#This model converges!
full_formula4 <- as.formula("NIntc_cover_binom ~ 
                            graz*aridity + graz*RASE + graz*AMT +
                            graz*pH + graz*SAC +
                            graz*log_nurse_meanH + graz*log_nurse_meanLDMC +
                            aridity*log_nurse_meanH + aridity*log_nurse_meanLDMC +
                            AMT*log_nurse_meanH + AMT*log_nurse_meanLDMC +
                            RASE*log_nurse_meanH + RASE*log_nurse_meanLDMC +
                            pH*log_nurse_meanH + pH*log_nurse_meanLDMC +
                            SAC*log_nurse_meanH + SAC*log_nurse_meanLDMC +
                            sin_lat + sin_long + #add these to account for spatial structure instead of (1|site_ID/ID)
                            (1|nurse_sp)")


# Fit the full model
options(na.action = "na.omit")
full_model_cov <- glmmTMB(
  formula = full_formula4,
  data = modeldat_final,    
  family = binomial  #had to remove sq terms and soil:climate and soil:trait interactions to make model converge
)


#do stepAIC first on the full models to eliminate some variables
stepwise_results <- stepAIC(full_model_cov, direction = "both")
#stepAIC says the best option is to remove no variables!!!

# Ensure all models maintain random effects by excluding them from being dropped
options(na.action = "na.fail") # Required for dredge function

#create a cluster obeject to run the function over 8 cores
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK" 
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))

# Export necessary objects and functions to the cluster
clusterExport(clust, varlist = c("modeldat_final", "full_formula4"), envir = environment())
clusterEvalQ(clust, library(glmmTMB))
clusterEvalQ(clust, library(MuMIn))

#peform model selection using 8 cores
# Perform model selection using dredge
cov_model_selection_par <- dredge(
  full_model_cov,
  fixed = c("cond(sin_lat)","cond(sin_long)"), #random effects are automatically included in all models due to the structure of tMB
  rank = "AIC", # Use AIC for model ranking
  cluster = clust) #start 18:43 on Thursday, finished 10:00 on Sunday (5 hour break without running)

# Stop the cluster after use
stopCluster(clust)

#save model selection results
saveRDS(cov_model_selection_par, "Functional trait data\\paper results\\nint_cover_nurse_trait_dredge_result.rds")

#read in model selection results
cov_model_selection_par <- readRDS("Functional trait data\\paper results\\nint_cover_nurse_trait_dredge_result.rds")

# Extract the best model based on AICc
cov_best_model <- get.models(cov_model_selection_par, subset = 1)[[1]] #subset = 1 gets the model with the absolute lowest AIC.

#extract models with AIC difference less than 2
cov_eq_model <- get.models(
  cov_model_selection_par, 
  subset = delta < 2 )

#average these models
cov_avg_models <- model.avg(cov_eq_model)

cov_avg_model_formula <-formula(cov_avg_models) #get formula of the averaged model

summary(cov_avg_models) #get summary of the averaged model

options(na.action = "na.omit")
cov_final_model <-glmmTMB(NIntc_cover_binom ~ AMT + aridity + graz + log_nurse_meanH + 
                            log_nurse_meanLDMC + pH + RASE + SAC + sin_lat + sin_long + 
                            AMT:log_nurse_meanH + AMT:log_nurse_meanLDMC + aridity:graz + 
                            graz:log_nurse_meanH + graz:log_nurse_meanLDMC + graz:RASE + 
                            graz:SAC + log_nurse_meanH:pH + log_nurse_meanH:RASE + log_nurse_meanH:SAC + 
                            log_nurse_meanLDMC:pH + log_nurse_meanLDMC:RASE + log_nurse_meanLDMC:SAC + (1|nurse_sp), 
                          data = modeldat_final, family = "binomial")
cov_null_model <- glmmTMB(NIntc_cover_binom ~ 1+ (1|nurse_sp), data = modeldat_final, family = "binomial")
r.squaredGLMM(cov_final_model) #r squared is not working

#use performance package to calculate R2 since Mumin isn't working
library(performance)
performance::r2(cov_final_model)

#get the variable importance
cov_importance <- sw(cov_avg_models)

#save to csv file
write.csv(as.data.frame(cov_importance), "Functional trait data\\paper results\\nint_cov_model_weights.csv")
