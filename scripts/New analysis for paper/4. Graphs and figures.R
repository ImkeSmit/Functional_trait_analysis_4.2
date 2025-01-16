###Script to make graphs and figures for the paper####
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)

####AKAIKE WEIGHTS FIGURE####
#import akaike weights of nint richness model
rich_sw <- read.csv("Functional trait data\\paper results\\nint_rich_model_weights.csv") |> 
  rename(var = X) |> 
  mutate(var2 = c("AMT", "graz", "H", "LDMC", "pH", "RASE", "SAC", "latitude", "longitude", 
                  "graz x LDMC", "graz x RASE", "graz x SAC", "H x pH", "LDMC x pH", "aridity", "LDMC x RASE", 
                  "AMT x LDMC", "AMT x H", "H x RASE", "LDMC x SAC", "H x SAC"))
rich_sw$var2 <- reorder(rich_sw$var2, rich_sw$importance)

#import akaike weights of nint cover model
cov_sw <- read.csv("Functional trait data\\paper results\\nint_cov_model_weights.csv") |> 
  rename(var = X, importance = cov_importance) |> 
  mutate(var2 = c("graz", "LDMC", "pH", "RASE", "SAC", "latitude", "longitude", 
                  "graz x LDMC", "graz x RASE", "graz x SAC", "LDMC x pH", "H", "H x pH", "AMT", "AMT x H", 
                  "aridity", "aridity x graz", "LDMC x SAC", "graz x H", "AMT x LDMC",  "LDMC x RASE", 
                  "H x RASE", "H x SAC"))
cov_sw$var2 <- reorder(cov_sw$var2, cov_sw$importance)


rich_sw_bar <- ggplot(aes(x = var2, y = importance), data = rich_sw) +
  geom_bar(stat = "identity") +
  xlab("Variable") +
  ylab("Importance") +
  coord_flip()+
  theme_classic()

cov_sw_bar <- ggplot(aes(x = var2, y = importance), data = cov_sw) +
  geom_bar(stat = "identity") +
  xlab("Variable") +
  ylab("Importance") +
  coord_flip()+
  theme_classic()

sw_combo <- ggarrange(rich_sw_bar, cov_sw_bar, nrow = 1, ncol = 2, labels = c("a", "b"))
ggsave("sw_bar.png", plot = sw_combo, path = "Figures", height = 1000, width = 2000, units = 'px')

####GRID OF VARIABLE EFFECTS ON NINTC####
#import modelling data
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
