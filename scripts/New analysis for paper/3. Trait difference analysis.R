##This script is to run the second analysis to determine whther the functional distance between dominant and 
#target species are affected by target species' association and environmental factors
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(car)

#Import pairwise differences between traits
trait_fdist <- read.csv("Functional trait data\\results\\trait_differences_between_2sp_traits_vary.csv", row.names = 1) |> 
  filter(trait %in% c("MaxH", "MeanLDMC"))
trait_fdist$SITE_ID <- as.factor(trait_fdist$SITE_ID)
trait_fdist$ID <- as.factor(trait_fdist$ID)
##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
ass$ID <- as.factor(ass$ID)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

#import siteinfo which has lat and long for each plot
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  mutate(plotref = str_c(SITE, PLOT, sep = "_")) |>
  select(ID, plotref, Lat_decimal, Long_decimal) |> 
  mutate(sin_lat = sin(Lat_decimal), 
         sin_long = sin(Long_decimal)) |> 
  select(!c(Lat_decimal, Long_decimal))

#import drypop, so which contains the env covariates
drypop <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Functional trait data\\Raw data\\drypop_20MAy.csv") |> 
  mutate(plotref = str_c(Site, Plot, sep = "_")) |> #create a variable to identify each plot
  dplyr::select(plotref, AMT, RAI, RASE, pH.b, SAC.b) |> 
  distinct() |> 
  left_join(siteinfo, by = "plotref") |> 
  dplyr::select(!plotref) |> 
  rename(pH = pH.b, SAC = SAC.b)
drypop$ID <- as.factor(drypop$ID)

#join the associations and the coordinates to the trait differences
trait_ass_join <- trait_fdist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) |> #only work with these associations
  left_join(drypop, by = "ID") |> 
  rename(nurse_sp = nurse)
trait_ass_join$association <- as.factor(trait_ass_join$association)
trait_ass_join$nurse <- as.factor(trait_ass_join$nurse)
trait_ass_join$SITE_ID <- as.factor(trait_ass_join$SITE_ID)
trait_ass_join$ID <- as.factor(trait_ass_join$ID)

####MODEL SELECTION FOR MAXH####
#MaxH model#
maxh_data <- trait_ass_join |> 
  filter(trait == "MaxH") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(maxh_data$neginv_trait_difference)
hist(maxh_data$trait_difference)

###full formula for model selection
full_formula <- as.formula("trait_difference ~ association*GRAZ + 
                           association*AMT + association*RASE + association*ARIDITY.v3 + 
                           association*SAC + association*pH +
                           sin_lat + sin_long + (1|nurse_sp)")

maxh_full_model <- glmmTMB(formula = full_formula, data = maxh_data)

#rank models with dredge
options(na.action = "na.fail")
maxh_model_selection <- dredge(maxh_full_model, 
                               fixed = c("cond(sin_lat)","cond(sin_long)"), #random effects are automatically included in all models due to the structure of tMB
                               rank = "AIC") #start16:08, end 16:12 wow

#get the bets models
maxh_eq_models <- get.models(maxh_model_selection, subset= delta <2) # the best model is the full model. There is no equivalent model


##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
maxh_test_nurse <- t.test(maxh_data[which(maxh_data$association == "nurse") , ]$trait_difference, 
                          mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0

maxh_test_bare <- t.test(maxh_data[which(maxh_data$association == "bare") , ]$trait_difference, 
                         mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0


####MODEL SELECTION FOR MeanLDMC####
ldmc_data <- trait_ass_join |> 
  filter(trait == "MeanLDMC") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(ldmc_data$trait_difference)

###full formula for model selection
full_formula <- as.formula("trait_difference ~ association*GRAZ + 
                           association*AMT + association*RASE + association*ARIDITY.v3 + 
                           association*SAC + association*pH + 
                           sin_lat + sin_long + (1|nurse_sp)")

ldmc_full_model <- glmmTMB(formula = full_formula, data = ldmc_data)

#rank models with dredge
options(na.action = "na.fail")
ldmc_model_selection <- dredge(ldmc_full_model, 
                               fixed = c("cond(sin_lat)","cond(sin_long)"), #random effects are automatically included in all models due to the structure of tMB
                               rank = "AIC") 

#get the bets models
ldmc_eq_models <- get.models(ldmc_model_selection, subset= delta <2) # the best model is the full model. There is no equivalent model

#average the equivalent models
avg_ldmc_models <- model.avg(ldmc_eq_models)
summary(avg_ldmc_models)

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
meanldmc_test_nurse <- t.test(ldmc_data[which(ldmc_data$association == "nurse") , ]$trait_difference, 
                              mu = 0, alternative = "two.sided")
#p <0.001, mean greater than 0

meanldmc_test_bare <- t.test(ldmc_data[which(ldmc_data$association == "bare") , ]$trait_difference, 
                             mu = 0, alternative = "two.sided")
#p <0.001, mean greater than 0