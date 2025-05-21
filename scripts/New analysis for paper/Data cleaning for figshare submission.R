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
