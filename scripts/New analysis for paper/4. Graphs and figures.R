###Script to make graphs and figures for the paper####
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(RColorBrewer)

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


###Get core df to make predictions over
temp_core0 <- modeldat_final |> 
  select(ID, nurse_sp, graz, RASE, SAC, aridity, pH, AMT,
         log_nurse_meanH, log_nurse_meanLDMC, sin_lat, sin_long) |> 
  distinct(ID, nurse_sp, .keep_all = T) |> 
  select(!ID) |> 
  mutate(graz = 0)
temp_core1 <- temp_core0 |> 
  mutate(graz = 1)
temp_core2 <- temp_core0 |> 
  mutate(graz = 2)
temp_core3 <- temp_core0 |> 
  mutate(graz = 3)

pred_dat_core <- bind_rows(temp_core0, temp_core1, temp_core2, temp_core3)
pred_dat_core$graz <- as.factor(pred_dat_core$graz) #this has every possible variable value for every grazing level

#best models (averaged model formula)
nintc_richness_bestmod <- glmmTMB(NIntc_richness_binom ~ AMT + aridity + graz + log_nurse_meanH + 
                                    log_nurse_meanLDMC + pH + RASE + SAC + sin_lat + sin_long + 
                                    AMT:log_nurse_meanH + AMT:log_nurse_meanLDMC + graz:log_nurse_meanLDMC + 
                                    graz:RASE + graz:SAC + log_nurse_meanH:pH + log_nurse_meanH:SAC + 
                                    log_nurse_meanLDMC:pH + log_nurse_meanLDMC:RASE + log_nurse_meanLDMC:SAC +
                                    sin_lat + sin_long,
                                  family = binomial, data = modeldat_final) #no RE for plotting

nintc_cover_bestmod <- glmmTMB(NIntc_cover_binom ~ AMT + aridity + graz + log_nurse_meanH + 
                                 log_nurse_meanLDMC + pH + RASE + SAC + sin_lat + sin_long + 
                                 AMT:log_nurse_meanH + AMT:log_nurse_meanLDMC + aridity:graz + 
                                 graz:log_nurse_meanH + graz:log_nurse_meanLDMC + graz:RASE + 
                                 graz:SAC + log_nurse_meanH:pH + log_nurse_meanH:RASE + log_nurse_meanH:SAC + 
                                 log_nurse_meanLDMC:pH + log_nurse_meanLDMC:RASE + log_nurse_meanLDMC:SAC + 
                                 sin_lat + sin_long, 
                               data = modeldat_final, family = binomial)

###NINtc richness~ RASE*graz
pred_dat1 <- pred_dat_core |> 
  mutate(sin_lat = mean(sin_lat), sin_long = mean(sin_long),
         AMT = mean(AMT), aridity = mean(aridity),
         pH = mean(pH), SAC = mean(SAC),
         log_nurse_meanLDMC = mean(log_nurse_meanLDMC), 
         log_nurse_meanH = mean(log_nurse_meanH)) 

pred_dat1$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, pred_dat1, type = "response")
pred_dat1$nintc_richness_true_prediction <- 2*pred_dat1$nintc_richness_binom_prediction -1 #backtransform from binomial

rich_RASE_graz <- ggplot(modeldat_final, aes(y = NIntc_richness, x = RASE)) +
  geom_jitter(height = 0.01, width = 2, color = "azure3", alpha = 0.4, size = 1.5) +
  geom_line(data = pred_dat1, aes(x = RASE, y = nintc_richness_true_prediction, color = graz), lwd = 1.5) +
  scale_color_manual(labels = c("ungrazed", "low", "medium", "high"),
                     values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" ))+
  labs(color = "Grazing pressure", y = expression(NInt[C]~richness), x = "RASE") +
  theme_classic() +
  theme(legend.position = "right")

###NIntc richness ~ SAC*graz
pred_dat2 <- pred_dat_core |> 
  mutate(sin_lat = mean(sin_lat), sin_long = mean(sin_long),
         AMT = mean(AMT), aridity = mean(aridity),
         pH = mean(pH), RASE = mean(RASE),
         log_nurse_meanLDMC = mean(log_nurse_meanLDMC), 
         log_nurse_meanH = mean(log_nurse_meanH)) 

pred_dat2$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, pred_dat2, type = "response")
pred_dat2$nintc_richness_true_prediction <- 2*pred_dat2$nintc_richness_binom_prediction -1 #backtransform from binomial

rich_SAC_graz <- ggplot(modeldat_final, aes(y = NIntc_richness, x = SAC)) +
  geom_jitter(height = 0.01, width = 2, color = "azure3", alpha = 0.4, size = 1.5) +
  geom_line(data = pred_dat2, aes(x = SAC, y = nintc_richness_true_prediction, color = graz), lwd = 1.5) +
  scale_color_manual(labels = c("ungrazed", "low", "medium", "high"),
                     values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" ))+
  labs(color = "Grazing pressure", y = expression(NInt[C]~richness), x = "SAC") +
  theme_classic() +
  theme(legend.position = "right")

###NIntc richness ~ LDMC*graz
pred_dat3 <- pred_dat_core |> 
  mutate(sin_lat = mean(sin_lat), sin_long = mean(sin_long),
         AMT = mean(AMT), aridity = mean(aridity),
         pH = mean(pH), RASE = mean(RASE), SAC = mean(SAC),
         log_nurse_meanH = mean(log_nurse_meanH)) 

pred_dat3$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, pred_dat3, type = "response")
pred_dat3$nintc_richness_true_prediction <- 2*pred_dat3$nintc_richness_binom_prediction -1 #backtransform from binomial

rich_LDMC_graz <- ggplot(modeldat_final, aes(y = NIntc_richness, x = log_nurse_meanLDMC)) +
  geom_jitter(height = 0.01, width = 0.01, color = "azure3", alpha = 0.4, size = 1.5) +
  geom_line(data = pred_dat3, aes(x = log_nurse_meanLDMC, y = nintc_richness_true_prediction, color = graz), lwd = 1.5) +
  scale_color_manual(labels = c("ungrazed", "low", "medium", "high"),
                     values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" ))+
  labs(color = "Grazing pressure", y = expression(NInt[C]~richness), x = "log(LDMC)") +
  theme_classic() +
  theme(legend.position = "right")

###NIntc richness ~ AMT*H
#to show the effect of H dependent on AMT, we should create categories of AMT
temp1 <- pred_dat_core |> 
  filter(graz == 1) |>  #we only need one level of graz
  mutate(sin_lat = mean(sin_lat), sin_long = mean(sin_long),
          aridity = mean(aridity),
         pH = mean(pH), RASE = mean(RASE), SAC = mean(SAC), #set other variables to their mean
         log_nurse_meanLDMC = mean(log_nurse_meanLDMC))
temp2 <- temp1 |> #create dataframes for 3 different values of AMT, then rbind tem
  mutate(AMT = mean(AMT), label = "mean")
temp3 <- temp1 |> 
  mutate(AMT = mean(AMT) + sd(AMT), label = "mean + sd")
temp4 <- temp1 |> 
  mutate(AMT = mean(AMT) - sd(AMT), label = "mean - sd")

temp2$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, temp2, type = "response")
temp2$nintc_richness_true_prediction <- 2*temp2$nintc_richness_binom_prediction -1 #backtransform from binomial

temp3$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, temp3, type = "response")
temp3$nintc_richness_true_prediction <- 2*temp3$nintc_richness_binom_prediction -1 #backtransform from binomial

temp4$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, temp4, type = "response")
temp4$nintc_richness_true_prediction <- 2*temp4$nintc_richness_binom_prediction -1 #backtransform from binomial

rich_H_AMT <- ggplot(modeldat_final, aes(y = NIntc_richness, x = log_nurse_meanH)) +
  geom_jitter(height = 0.01, width = 0.01, color = "azure3", alpha = 0.4, size = 1.5) +
  geom_line(data = temp2, aes(x = log_nurse_meanH, y = nintc_richness_true_prediction, color = "mean"), lwd = 1.5) +
  geom_line(data = temp3, aes(x = log_nurse_meanH, y = nintc_richness_true_prediction, color = "mean + sd"), lwd = 1.5) +
  geom_line(data = temp4, aes(x = log_nurse_meanH, y = nintc_richness_true_prediction, color = "mean - sd"), lwd = 1.5) +
  scale_color_manual(name = "Value of AMT", 
                     breaks = c("mean", "mean + sd", "mean - sd"), 
                     values = c("mean" = brewer.pal(8, "YlOrRd")[6], "mean + sd" = brewer.pal(8, "YlOrRd")[8], 
                                "mean - sd" = brewer.pal(8, "YlOrRd")[4])) +
  labs(y = expression(NInt[C]~richness), x = "log(H)") +
  theme_classic() 


###NIntc richness ~ LDMC*H
#to show the effect of H dependent on LDMC, we should create categories of AMT
temp1 <- pred_dat_core |> 
  filter(graz == 1) |>  #we only need one level of graz
  mutate(sin_lat = mean(sin_lat), sin_long = mean(sin_long),
         aridity = mean(aridity),
         pH = mean(pH), RASE = mean(RASE), SAC = mean(SAC), #set other variables to their mean
         log_nurse_meanH = mean(log_nurse_meanH))
temp2 <- temp1 |> #create dataframes for 3 different values of AMT, then rbind tem
  mutate(AMT = mean(AMT), label = "mean")
temp3 <- temp1 |> 
  mutate(AMT = mean(AMT) + sd(AMT), label = "mean + sd")
temp4 <- temp1 |> 
  mutate(AMT = mean(AMT) - sd(AMT), label = "mean - sd")

temp2$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, temp2, type = "response")
temp2$nintc_richness_true_prediction <- 2*temp2$nintc_richness_binom_prediction -1 #backtransform from binomial

temp3$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, temp3, type = "response")
temp3$nintc_richness_true_prediction <- 2*temp3$nintc_richness_binom_prediction -1 #backtransform from binomial

temp4$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, temp4, type = "response")
temp4$nintc_richness_true_prediction <- 2*temp4$nintc_richness_binom_prediction -1 #backtransform from binomial

rich_LDMC_AMT <- ggplot(modeldat_final, aes(y = NIntc_richness, x = log_nurse_meanLDMC)) +
  geom_jitter(height = 0.01, width = 0.01, color = "azure3", alpha = 0.4, size = 1.5) +
  geom_line(data = temp2, aes(x = log_nurse_meanLDMC, y = nintc_richness_true_prediction, color = "mean"), lwd = 1.5) +
  geom_line(data = temp3, aes(x = log_nurse_meanLDMC, y = nintc_richness_true_prediction, color = "mean + sd"), lwd = 1.5) +
  geom_line(data = temp4, aes(x = log_nurse_meanLDMC, y = nintc_richness_true_prediction, color = "mean - sd"), lwd = 1.5) +
  scale_color_manual(name = "Value of AMT", 
                     breaks = c("mean", "mean + sd", "mean - sd"), 
                     values = c("mean" = brewer.pal(8, "YlOrRd")[6], "mean + sd" = brewer.pal(8, "YlOrRd")[8], 
                                "mean - sd" = brewer.pal(8, "YlOrRd")[4])) +
  labs(y = expression(NInt[C]~richness), x = "log(LDMC)") +
  theme_classic() 
