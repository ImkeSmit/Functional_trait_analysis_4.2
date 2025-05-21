##Package and export the data for submission to figshare####
library(tidyverse)
library(tidylog)
library(openxlsx)

####1. model selection for nint ~ nurse traits + env####
#exact code used in the model selection nint nurse traits script: 
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
  dplyr::select(plotref, Country, AMT, RAI, RASE, pH.b, SAC.b) |> 
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
#we lose many plots because of the na dropping
#if we do this again we should create separate datasets for the nintc richness and nintc cover becayuse the nintc cover is more gappy

###remove unnesecary variables:
nint_nurse_trait_data_for_export <- modeldat_final |> 
  select(c(ID, site_ID, replicate_no, nurse_sp, NIntc_richness, NIntc_cover, NIntc_richness_binom, NIntc_cover_binom, 
           aridity, graz, RASE, pH, AMT, SAC,  nurse_mean_H, nurse_meanLDMC, log_nurse_meanH, log_nurse_meanLDMC, sin_lat, sin_long))

#now export
write.xlsx(nint_nurse_trait_data_for_export,
            "C:\\Users\\imke6\\Documents\\Msc Projek\\Interaction-environment manuscript (mock paper)\\submission to GEB\\data submission\\nint_nurse_trait_env_data.xlsx" 
            )


####2. Trait difference analysis####
#copy code from script for trait difference analysis:
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
  dplyr::select(plotref, Country, AMT, RAI, RASE, pH.b, SAC.b) |> 
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
trait_ass_join$GRAZ <- as.factor(trait_ass_join$GRAZ)


###remove unnecessary variables
trait_difference_data_for_export <- trait_ass_join |> 
  select(SITE_ID, ID, replicate,trait, trait_difference, nurse_sp, target,association, GRAZ, ARIDITY.v3,  
         AMT, RASE, pH, SAC, sin_lat, sin_long) |> 
  rename(site_ID = SITE_ID, 
         replicate_no = replicate, 
         aridity = ARIDITY.v3, 
         graz = GRAZ)
  