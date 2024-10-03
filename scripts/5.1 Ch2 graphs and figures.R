###CHAPTER 2 GRAPHS AND FIGURES###
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(multcomp)
library(multcompView)
library(glmmTMB)
library(car)
library(tidyverse)
library(tidylog)

####RaoQ accross aridity and grazing####
#import FD results
FD_results <- read.csv("Functional trait data\\results\\FD_results_20May2024.csv", row.names = 1)
FD_results$ID <- as.factor(FD_results$ID)
#add grazing and aridity
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(c(ID, GRAZ, ARIDITY.v3))
siteinfo$ID <-as.factor(siteinfo$ID)
#do the join
long_FD_results <- FD_results |> 
  inner_join(siteinfo, by = "ID") |> 
  select(!c(qual.FRic, nsp)) |> 
  pivot_longer(cols = c(FRic, FEve, FDiv, RaoQ), names_to = "FD_metrics", values_to = "value" )
long_FD_results$GRAZ <- as.factor(long_FD_results$GRAZ)

###RAOQ###
Raoq_aridity <- long_FD_results |> 
  filter(FD_metrics == "RaoQ") |> 
  ggplot(aes(y = log(value), x = ARIDITY.v3)) +
  geom_point(alpha = 0.6, colour = "darkslategrey")+
  theme_classic() +
  xlab("Aridity") +
  ylab("ln(RaoQ)")

Raoq_graz <- long_FD_results |> 
  filter(FD_metrics == "RaoQ") |> 
  ggplot(aes(y = log(value), x = GRAZ, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6)+
  stat_summary(fun = mean, geom="point", shape = 23, size = 2, fill = "white", color = "black") +
  scale_fill_manual(values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4")) +
  theme_classic() +
  ylab("") +
  xlab("Grazing pressure") +
  scale_x_discrete(labels = c("Ungrazed", "Low", "Medium", "High")) +
  theme(legend.position = "none")
  
Raoq_graz_arid <- ggarrange(Raoq_aridity, Raoq_graz, ncol = 2, nrow= 1, labels = c("a", "b"))
ggsave("Raoq_aridity_grazing.png", Raoq_graz_arid, height = 800, width = 1500, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")


##Nint ~ RaoQ####
#import nint results
#Add the NIntc
nint_result <- 
  read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  filter(ID %in% c(FD_results$ID))
nint_result$ID <- as.factor(nint_result$ID)


##summarise the nintc by plot
nint_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NIntc_richness, NIntc_cover, NInta_richness, NInta_cover) |> 
  group_by(ID) |> 
  #calculate mean and sd of NIntc richness
  mutate(mean_NIntc_rich = mean(NIntc_richness, na.rm = TRUE), 
         se_NIntc_rich = sd(NIntc_richness, na.rm = TRUE)/sqrt(sum(!is.na(NIntc_richness)))) |> 
  #calculate the mean and sd of NIntc cover
  mutate(mean_NIntc_cov = mean(NIntc_cover, na.rm = TRUE), 
         se_NIntc_cov = sd(NIntc_cover, na.rm = T)/sqrt(sum(!is.na(NIntc_cover)))) |> 
  #calculate mean and sd of NInta richness
  mutate(mean_NInta_rich = mean(NInta_richness, na.rm = T), 
         se_NInta_rich = sd(NInta_richness, na.rm = T)/sqrt(sum(!is.na(NInta_richness)))) |> 
  #calculate the mean and sd of NInta cover
  mutate(mean_NInta_cov = mean(NInta_cover, na.rm = T), 
         se_NInta_cov = sd(NInta_cover, na.rm = T)/sqrt(sum(!is.na(NInta_cover)))) |>
  ungroup() |> 
  select(!c(NIntc_richness, NIntc_cover, NInta_richness, NInta_cover)) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  filter(!is.na(FEve))


##NIntc against the RaoQ metrics##
rich_raoq <- ggplot(nint_sum, aes(x = log(RaoQ), y = mean_NIntc_rich)) +
  geom_point(alpha = 0.6, colour = "darkslategrey") +
  ylab(expression(NInt[C]~richness))+
  xlab("ln(RaoQ)") +
  geom_errorbar(aes(ymin = mean_NIntc_rich - se_NIntc_rich, ymax =  mean_NIntc_rich + se_NIntc_rich), 
                width = 0.01, color = "darkslategrey") +
  theme_classic()

cov_raoq <- ggplot(nint_sum, aes(x = log(RaoQ), y = mean_NIntc_cov)) +
  geom_point(alpha = 0.6, colour = "darkslategrey") +
  ylab(expression(NInt[C]~cover))+
  xlab("ln(RaoQ)") +
  geom_errorbar(aes(ymin = mean_NIntc_cov - se_NIntc_cov, ymax =  mean_NIntc_cov + se_NIntc_cov), 
                width = 0.01, color = "darkslategrey") +
  theme_classic()

nintc_Raoq <- ggarrange(rich_raoq, cov_raoq, ncol = 2, nrow = 1, labels = c("a", "b"))
ggsave("NIntc_RaoQ_scatterplots.png", nintc_Raoq, height = 800, width = 1500, units = "px", 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")


##NInta against the RaoQ metrics##
ninta_rich_raoq <- ggplot(nint_sum, aes(x = log(RaoQ), y = mean_NInta_rich)) +
  geom_point(alpha = 0.6, colour = "darkslategrey") +
  ylab(expression(NInt[A]~richness))+
  xlab("ln(RaoQ)") +
  geom_errorbar(aes(ymin = mean_NInta_rich - se_NInta_rich, ymax =  mean_NInta_rich + se_NInta_rich), 
                width = 0.01, color = "darkslategrey") +
  theme_classic()

ninta_cov_raoq <- ggplot(nint_sum, aes(x = log(RaoQ), y = mean_NInta_cov)) +
  geom_point(alpha = 0.6, colour = "darkslategrey") +
  ylab(expression(NInt[A]~cover))+
  xlab("ln(RaoQ)") +
  geom_errorbar(aes(ymin = mean_NInta_cov - se_NInta_cov, ymax =  mean_NInta_cov + se_NInta_cov), 
                width = 0.01, color = "darkslategrey") +
  theme_classic()

ninta_Raoq <- ggarrange(ninta_rich_raoq, ninta_cov_raoq, ncol = 2, nrow = 1, labels = c("a", "b"))
ggsave("NInta_RaoQ_scatterplots.png", ninta_Raoq, height = 800, width = 1500, units = "px", 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")


####NIntc ~ trait graphs####
#import the nint nurse traits
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
##we will make scatterplots with model predictions overlaid

#choose a colour for the line;
chosen_col <- brewer.pal(6, "Dark2")[6]

###Get core df to make predictions over
pred_dat_core <- modeldat_final |> 
  select(ID, site_ID, nurse_sp, graz, RASE, SAC, aridity, pH,
         log_nurse_meanLA, log_nurse_meanH, log_nurse_meanSLA, log_nurse_meanCNratio) |> 
  distinct(ID, nurse_sp, .keep_all = T)

#best models
nintc_richness_bestmod <- glmmTMB(NIntc_richness_binom ~  graz+RASE+SAC+log_nurse_meanLA+log_nurse_meanH+
                                    log_nurse_meanCNratio+graz:SAC, 
                                  family = binomial, data = modeldat_final)

nintc_cover_bestmod <- glmmTMB(NIntc_cover_binom ~  graz+aridity+RASE+pH+SAC+
                                 log_nurse_meanLA+log_nurse_meanSLA+log_nurse_meanH+log_nurse_meanCNratio+
                                 graz:RASE+graz:pH+graz:SAC , 
                               family = binomial, data = modeldat_final)

###NINtc richness~ C:N
pred_dat1 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanCNratio)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except CN to their mean

pred_dat1$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, pred_dat1, type = "response")
pred_dat1$nintc_richness_true_prediction <- 2*pred_dat1$nintc_richness_binom_prediction -1 #backtransform from binomial

#how many points on graph?
modeldat_final |> 
  filter(!is.na(nurse_mean_C_N_ratio) & !is.na(NIntc_richness)) |> 
  summarise(n = n()) #2659

nintc_richness_CN <- ggplot(modeldat_final, aes(x = log_nurse_meanCNratio, y = NIntc_richness)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[C]~richness)) +
  #xlab("log(C:N)") +
  xlab(" ")+
  geom_line(data = pred_dat1, 
            aes(x = log_nurse_meanCNratio, y = nintc_richness_true_prediction), color = chosen_col, lwd = 1)

###NINtc richness~ LA
pred_dat3 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanLA)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except LA to their mean

pred_dat3$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, pred_dat3, type = "response")
pred_dat3$nintc_richness_true_prediction <- 2*pred_dat3$nintc_richness_binom_prediction -1 #backtransform from binomial

nintc_richness_LA <- ggplot(modeldat_final, aes(x = log_nurse_meanLA, y = NIntc_richness)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  xlab(" ")+
  #ylab(expression(NInt[C]~richness)) +
  #xlab("log(LA)") +
  geom_line(data = pred_dat3, 
            aes(x = log_nurse_meanLA, y = nintc_richness_true_prediction), color = chosen_col, lwd = 1)

###NINtc richness~ H
pred_dat4 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanH)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), 
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T)) #set all variables except H to their mean

pred_dat4$nintc_richness_binom_prediction <- predict(nintc_richness_bestmod, pred_dat4, type = "response")
pred_dat4$nintc_richness_true_prediction <- 2*pred_dat4$nintc_richness_binom_prediction -1 #backtransform from binomial

nintc_richness_H <- ggplot(modeldat_final, aes(x = log_nurse_meanH, y = NIntc_richness)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  xlab(" ")+
  #ylab(expression(NInt[C]~richness)) +
  #xlab("log(H)") +
  geom_line(data = pred_dat4, 
            aes(x = log_nurse_meanH, y = nintc_richness_true_prediction), color = chosen_col, lwd = 1)


###NIntc cover ~ C:N
pred_dat2 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanCNratio)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T), log_nurse_meanSLA = mean(log_nurse_meanSLA, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except CN to their mean

pred_dat2$nintc_cover_binom_prediction <- predict(nintc_cover_bestmod, pred_dat2, type = "response")
pred_dat2$nintc_cover_true_prediction <- 2*pred_dat2$nintc_cover_binom_prediction -1 #backtransform from binomial

#how many points on graph?
modeldat_final |> 
  filter(!is.na(nurse_mean_C_N_ratio) & !is.na(NIntc_cover)) |> 
  summarise(n = n()) #2659

nintc_cover_CN <- ggplot(modeldat_final, aes(x = log_nurse_meanCNratio, y = NIntc_cover)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[C]~cover)) +
  xlab("log(C:N)") +
  geom_line(data = pred_dat2, 
            aes(x = log_nurse_meanCNratio, y = nintc_cover_true_prediction), color = chosen_col, lwd = 1)

###NINtc cover~ LA
pred_dat5 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanLA)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), log_nurse_meanSLA = mean(log_nurse_meanSLA, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except CN to their mean

pred_dat5$nintc_cover_binom_prediction <- predict(nintc_cover_bestmod, pred_dat5, type = "response")
pred_dat5$nintc_cover_true_prediction <- 2*pred_dat5$nintc_cover_binom_prediction -1 #backtransform from binomial

nintc_cover_LA <- ggplot(modeldat_final, aes(x = log_nurse_meanLA, y = NIntc_cover)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  #ylab(expression(NInt[C]~cover)) +
  xlab("log(LA)") +
  geom_line(data = pred_dat5, 
            aes(x = log_nurse_meanLA, y = nintc_cover_true_prediction), color = chosen_col, lwd = 1)

###NINtc cover~ H
pred_dat6 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanH)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), log_nurse_meanSLA = mean(log_nurse_meanSLA, na.rm = T), 
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T)) #set all variables except CN to their mean

pred_dat6$nintc_cover_binom_prediction <- predict(nintc_cover_bestmod, pred_dat6, type = "response")
pred_dat6$nintc_cover_true_prediction <- 2*pred_dat6$nintc_cover_binom_prediction -1 #backtransform from binomial

nintc_cover_H <- ggplot(modeldat_final, aes(x = log_nurse_meanH, y = NIntc_cover)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  #ylab(expression(NInt[C]~cover)) +
  xlab("log(H)") +
  geom_line(data = pred_dat6, 
            aes(x = log_nurse_meanH, y = nintc_cover_true_prediction), color = chosen_col, lwd = 1)


###NINtc cover~ SLA
pred_dat7 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanSLA)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), log_nurse_meanH = mean(log_nurse_meanH, na.rm = T), 
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T)) #set all variables except CN to their mean

pred_dat7$nintc_cover_binom_prediction <- predict(nintc_cover_bestmod, pred_dat7, type = "response")
pred_dat7$nintc_cover_true_prediction <- 2*pred_dat7$nintc_cover_binom_prediction -1 #backtransform from binomial

nintc_cover_SLA <- ggplot(modeldat_final, aes(x = log_nurse_meanSLA, y = NIntc_cover)) +
  geom_jitter(height = 0.05,width = 0.01, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  #ylab(expression(NInt[C]~cover)) +
  xlab("log(SLA)") +
  geom_line(data = pred_dat7, 
            aes(x = log_nurse_meanSLA, y = nintc_cover_true_prediction), color = chosen_col, lwd = 1)

#arrange the above seven figures on the same plot
#create a blank plot to act as a placeholder
blank_plot <- ggplot() + theme_void()

nintc_nurse_traits <- ggarrange(ggarrange(nintc_richness_CN, nintc_richness_H,nintc_richness_LA, blank_plot, nrow = 1, ncol = 4, labels = c("a", "b", "c")),
                                ggarrange(nintc_cover_CN, nintc_cover_H,nintc_cover_LA, nintc_cover_SLA, nrow = 1, ncol = 4, labels = c("d", "e", "f", "g")),
                                nrow = 2, ncol = 1)
ggsave("nintc_trait_scatterplots.png", nintc_nurse_traits, height = 1500, width = 2500, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")



####NInta ~ trait graphs####
#import the nint nurse traits
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
##we will make scatterplots with model predictions overlaid

#choose a colour for the line;
chosen_col <- brewer.pal(6, "Dark2")[6]

###Get core df to make predictions over
pred_dat_core <- modeldat_final |> 
  select(ID, site_ID, nurse_sp, graz, RASE, SAC, aridity, pH,
         log_nurse_meanLA, log_nurse_meanH, log_nurse_meanSLA, log_nurse_meanCNratio) |> 
  distinct(ID, nurse_sp, .keep_all = T)

#best models
ninta_richness_bestmod <- glmmTMB(NInta_richness_binom ~ graz+pH+SAC+log_nurse_meanLA+log_nurse_meanH+
                                    log_nurse_meanCNratio+graz:pH+graz:SAC, 
                                  family = binomial, data = modeldat_final)


ninta_cover_bestmod <- glmmTMB(NInta_cover_binom ~ graz+aridity+pH+SAC+log_nurse_meanLA+log_nurse_meanSLA+log_nurse_meanH+
                                 log_nurse_meanCNratio+graz:pH+graz:SAC, 
                               family = binomial, data = modeldat_final)

###NINtc richness~ C:N
pred_dat1 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanCNratio)) |> 
  mutate(graz = 1, SAC = mean(SAC), pH = mean(pH), log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except CN to their mean

pred_dat1$ninta_richness_binom_prediction <- predict(ninta_richness_bestmod, pred_dat1, type = "response")
pred_dat1$ninta_richness_true_prediction <- 3*pred_dat1$ninta_richness_binom_prediction -1 #backtransform from binomial

#how many points on graph?
modeldat_final |> 
  filter(!is.na(nurse_mean_C_N_ratio) & !is.na(NInta_richness)) |> 
  summarise(n = n()) #2659

ninta_richness_CN <- ggplot(modeldat_final, aes(x = log_nurse_meanCNratio, y = NInta_richness)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[A]~richness)) +
  #xlab("log(C:N)") +
  xlab(" ")+
  geom_line(data = pred_dat1, 
            aes(x = log_nurse_meanCNratio, y = ninta_richness_true_prediction), color = chosen_col, lwd = 1)

###NINtc richness~ LA
pred_dat3 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanLA)) |> 
  mutate(graz = 1, SAC = mean(SAC), pH = mean(pH), log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except LA to their mean

pred_dat3$ninta_richness_binom_prediction <- predict(ninta_richness_bestmod, pred_dat3, type = "response")
pred_dat3$ninta_richness_true_prediction <- 3*pred_dat3$ninta_richness_binom_prediction -1 #backtransform from binomial

ninta_richness_LA <- ggplot(modeldat_final, aes(x = log_nurse_meanLA, y = NInta_richness)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  xlab(" ")+
  #ylab(expression(NInt[A]~richness)) +
  #xlab("log(LA)") +
  geom_line(data = pred_dat3, 
            aes(x = log_nurse_meanLA, y = ninta_richness_true_prediction), color = chosen_col, lwd = 1)

###NINtc richness~ H
pred_dat4 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanH)) |> 
  mutate(graz = 1, SAC = mean(SAC), pH = mean(pH), log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), 
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T)) #set all variables except H to their mean

pred_dat4$ninta_richness_binom_prediction <- predict(ninta_richness_bestmod, pred_dat4, type = "response")
pred_dat4$ninta_richness_true_prediction <- 3*pred_dat4$ninta_richness_binom_prediction -1 #backtransform from binomial

ninta_richness_H <- ggplot(modeldat_final, aes(x = log_nurse_meanH, y = NInta_richness)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  xlab(" ")+
  #ylab(expression(NInt[A]~richness)) +
  #xlab("log(H)") +
  geom_line(data = pred_dat4, 
            aes(x = log_nurse_meanH, y = ninta_richness_true_prediction), color = chosen_col, lwd = 1)


###NIntc cover ~ C:N
pred_dat2 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanCNratio)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T), log_nurse_meanSLA = mean(log_nurse_meanSLA, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except CN to their mean

pred_dat2$ninta_cover_binom_prediction <- predict(ninta_cover_bestmod, pred_dat2, type = "response")
pred_dat2$ninta_cover_true_prediction <- 3*pred_dat2$ninta_cover_binom_prediction -1 #backtransform from binomial

#how many points on graph?
modeldat_final |> 
  filter(!is.na(nurse_mean_C_N_ratio) & !is.na(NInta_cover)) |> 
  summarise(n = n()) #2659

ninta_cover_CN <- ggplot(modeldat_final, aes(x = log_nurse_meanCNratio, y = NInta_cover)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[A]~cover)) +
  xlab("log(C:N)") +
  geom_line(data = pred_dat2, 
            aes(x = log_nurse_meanCNratio, y = ninta_cover_true_prediction), color = chosen_col, lwd = 1)

###NINtc cover~ LA
pred_dat5 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanLA)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), log_nurse_meanSLA = mean(log_nurse_meanSLA, na.rm = T), 
         log_nurse_meanH = mean(log_nurse_meanH, na.rm = T)) #set all variables except CN to their mean

pred_dat5$ninta_cover_binom_prediction <- predict(ninta_cover_bestmod, pred_dat5, type = "response")
pred_dat5$ninta_cover_true_prediction <- 3*pred_dat5$ninta_cover_binom_prediction -1 #backtransform from binomial

ninta_cover_LA <- ggplot(modeldat_final, aes(x = log_nurse_meanLA, y = NInta_cover)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  #ylab(expression(NInt[A]~cover)) +
  xlab("log(LA)") +
  geom_line(data = pred_dat5, 
            aes(x = log_nurse_meanLA, y = ninta_cover_true_prediction), color = chosen_col, lwd = 1)

###NINtc cover~ H
pred_dat6 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanH)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), log_nurse_meanSLA = mean(log_nurse_meanSLA, na.rm = T), 
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T)) #set all variables except CN to their mean

pred_dat6$ninta_cover_binom_prediction <- predict(ninta_cover_bestmod, pred_dat6, type = "response")
pred_dat6$ninta_cover_true_prediction <- 3*pred_dat6$ninta_cover_binom_prediction -1 #backtransform from binomial

ninta_cover_H <- ggplot(modeldat_final, aes(x = log_nurse_meanH, y = NInta_cover)) +
  geom_jitter(height = 0.05, width = 0.05, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  #ylab(expression(NInt[A]~cover)) +
  xlab("log(H)") +
  geom_line(data = pred_dat6, 
            aes(x = log_nurse_meanH, y = ninta_cover_true_prediction), color = chosen_col, lwd = 1)


###NINtc cover~ SLA
pred_dat7 <- pred_dat_core |> 
  filter(!is.na(log_nurse_meanSLA)) |> 
  mutate(graz = 1, SAC = mean(SAC), RASE = mean(RASE), pH = mean(pH), aridity = mean(aridity),
         log_nurse_meanCNratio = mean(log_nurse_meanCNratio, na.rm = T), log_nurse_meanH = mean(log_nurse_meanH, na.rm = T), 
         log_nurse_meanLA = mean(log_nurse_meanLA, na.rm = T)) #set all variables except CN to their mean

pred_dat7$ninta_cover_binom_prediction <- predict(ninta_cover_bestmod, pred_dat7, type = "response")
pred_dat7$ninta_cover_true_prediction <- 3*pred_dat7$ninta_cover_binom_prediction -1 #backtransform from binomial

ninta_cover_SLA <- ggplot(modeldat_final, aes(x = log_nurse_meanSLA, y = NInta_cover)) +
  geom_jitter(height = 0.05,width = 0.01, alpha = 0.6, size = 1, colour = "darkslategrey") +
  theme_classic() +
  ylab(" ")+
  #ylab(expression(NInt[A]~cover)) +
  xlab("log(SLA)") +
  geom_line(data = pred_dat7, 
            aes(x = log_nurse_meanSLA, y = ninta_cover_true_prediction), color = chosen_col, lwd = 1)

#arrange the above seven figures on the same plot
#create a blank plot to act as a placeholder
blank_plot <- ggplot() + theme_void()

ninta_nurse_traits <- ggarrange(ggarrange(ninta_richness_CN, ninta_richness_H,ninta_richness_LA, blank_plot, nrow = 1, ncol = 4, labels = c("a", "b", "c")),
                                ggarrange(ninta_cover_CN, ninta_cover_H,ninta_cover_LA, ninta_cover_SLA, nrow = 1, ncol = 4, labels = c("d", "e", "f", "g")),
                                nrow = 2, ncol = 1)
ggsave("ninta_trait_scatterplots.png", ninta_nurse_traits, height = 1500, width = 2500, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")

###Fdist ~ association####
#import functional distances
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp_traits_varying.csv", row.names = 1)

###Lets join the results of the CHi2 tests to twosp-dist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

dist_ass_join <- twosp_dist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) #only work with these associations
dist_ass_join$association <- as.factor(dist_ass_join$association)

stars <- data.frame(symbols = c("*", "*"), ycoord = c(7.9, 7.9), association = c("bare", "nurse"))
#letters <- data.frame(symbols = c("a", "b"), ycoord = c(7.9, 7.9), association = c("bare", "nurse"))

dist_ass <- ggplot(dist_ass_join, aes(x = association, y = euclidean_dist, fill = association)) +
            geom_boxplot(alpha = 0.6) +
            scale_fill_manual(values = c(brewer.pal(8, "Dark2")[7], brewer.pal(8, "Dark2")[1]))+
            stat_summary(fun = mean, geom="point", shape = 23, size = 2, fill = "white", color = "black") +
            geom_text(data = stars, aes(x = association, y = ycoord, label = symbols), size = 7, color = "brown3") +
            #geom_text(data = letters, aes(x = association, y = ycoord, label = symbols), size = 5, color = "brown3") +
            ylim(0,8) +
            ylab("Distance") +
            xlab("Target species association") +
            scale_x_discrete(labels = c("bare", "dominant")) +
            theme_classic() +
            theme(legend.position = "none")
ggsave("fdist_association_boxplot.png", dist_ass, height = 1000, width = 800, units = "px", 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")


###one dimensional Fdist~association with raw data####
#import 1D trait difference data
trait_diff <- read.csv("Functional trait data\\results\\trait_differences_between_2sp_traits_vary.csv", row.names = 1)
##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  dplyr::select(ID, species, association) |> 
  rename(target = species)

#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites
trait_ass_join <- trait_diff |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) |>  #only work with these associations
  mutate(trait = case_when(trait == "MaxH" ~ "H~(cm)", #rename the traits so that they are labelled nicely in the plot
                   trait == "MaxLS" ~ "LS~(cm^2)",
                   trait == "MeanLA" ~ "LA~(cm^2)",
                   trait == "MeanLDMC" ~ "LDMC~('%')",
                   trait == "MeanLL" ~ "LL~(cm)",
                   trait == "MeanSLA" ~ "SLA~(cm^2/g)", 
                   trait == "C_N_ratio" ~ "C:N"))
trait_ass_join$association <- as.factor(trait_ass_join$association)
trait_ass_join$trait_difference <- as.numeric(trait_ass_join$trait_difference)

annotations <- data.frame(trait = rep(c(unique(trait_ass_join$trait)),2), 
                             association = c(rep("bare", 7), rep("nurse", 7)))
annotations<- annotations[order(annotations$trait) , ]
annotations$t_test_significance <- c("  *", "  *", #C:N
                                        "  *", "  *", #H
                                        "  *", " ", #LA
                                        "  *", "  *", #LDMC
                                        "  *", "  *", #LL
                                        "  *", "  *", #LS
                                        "  *", "  *" )#SLA
annotations$anova_significance <- c("    a", "    b", #C:N
                                     "    a", "    b", #H
                                    "    a", "    b", #LA
                                     "    a", "    b", #LDMC
                                     " ", " ", #LL
                                    "    a", "    b", #LS
                                     "    a", "    b" )#SLA
annotations$ycoord_t <- c(15, 8, #C:N
                          100, 150, #H
                          2, 2, #LA
                          0.2, 0.15, #LDMC
                          8, 5, #LL
                          50000, 90000, #LS
                          20, 20)#SLA

annotations$ycoord_anova <- c(32, 32, #C:N
                          640, 640, #H
                          23, 23, #LA
                          0.7, 0.7, #LDMC
                          0, 0, #LL
                          700000, 700000, #LS
                          95, 95)#SLA

trait_distances <- ggplot(trait_ass_join, aes(x = association, y = trait_difference, fill = association)) +
  geom_boxplot(alpha = 0.6)+
  stat_summary(fun = mean, geom="point", shape = 23, size = 2, fill = "white", color = "black") +
  scale_fill_manual(values = c(brewer.pal(8, "Dark2")[7], brewer.pal(8, "Dark2")[1])) + 
  scale_x_discrete(labels = c(expression("∆"["Db"]), expression("∆"["Dd"]))) +
  facet_wrap(~trait, scales = "free_y", labeller = label_parsed) +
  ylab("Difference") +
  xlab("Target species association") +
  geom_text(data = annotations, aes(x = annotations$association, y = ycoord_t, label = t_test_significance), size = 8, color = "brown3")+
  geom_text(data = annotations, aes(x = annotations$association, y = ycoord_anova, label = anova_significance), size = 5, color = "brown3")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = (element_text(size = 11)))

ggsave("one_dimensional_trait_distances.png", trait_distances, path = "Figures", 
       height = 1700, width = 2000, units = "px")


###one dimensional Fdist~association with model estimates####
#import the data we model trait differences with
trait_fdist <- read.csv("Functional trait data\\results\\trait_differences_between_2sp_traits_vary.csv", row.names = 1)
trait_fdist$SITE_ID <- as.factor(trait_fdist$SITE_ID)
##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

trait_ass_join <- trait_fdist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) #only work with these associations
trait_ass_join$association <- as.factor(trait_ass_join$association)
trait_ass_join$nurse <- as.factor(trait_ass_join$nurse)
trait_ass_join$SITE_ID <- as.factor(trait_ass_join$SITE_ID)
trait_ass_join$ID <- as.factor(trait_ass_join$ID)


###get the predictions for each trait model

traits <- c("MaxH", "MaxLS", "MeanLA", "MeanLDMC", "MeanLL", "MeanSLA", "C_N_ratio")

for(t in 1:length(traits)) {
  one_trait_data <- trait_ass_join |> 
    filter(trait == traits[t]) 
  
  #create data to predict over
  pred_df_nurse <- one_trait_data |> 
    dplyr::select(SITE_ID, ID, nurse) |> 
    distinct() |> 
    mutate(association = "nurse")
  
  pred_df_bare <- one_trait_data |> 
    dplyr::select(SITE_ID, ID, nurse) |> 
    distinct() |> 
    mutate(association = "bare")
  
  pred_df_final <- rbind(pred_df_nurse, pred_df_bare)
  
  if(t %in% c(2,6)) { #for MaxLS and SLA we need to use lmer, not TMB
    #select the model to use for prediction
    one_trait_model <- lmer(formula = as.formula("trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID)"), data = one_trait_data)
  } else {
  #select the model to use for prediction
  one_trait_model <- glmmTMB(formula = as.formula("trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID)"), data = one_trait_data) 
  }
  
  pred_df_final$trait_diff_predictions <- predict(one_trait_model, pred_df_final)
  pred_df_final$trait <- traits[t]
  
  if(t == 1) {
    predictions_combo <- pred_df_final
  } else {predictions_combo <- rbind(predictions_combo, pred_df_final) }
  
  rm(pred_df_final)
}


###Now we can make the graph
#Change the labels of the traits
predictions_combo2 <- predictions_combo |> 
  mutate(trait = case_when(trait == "MaxH" ~ "H~(cm)", #rename the traits so that they are labelled nicely in the plot
                           trait == "MaxLS" ~ "LS~(cm^2)",
                           trait == "MeanLA" ~ "LA~(cm^2)",
                           trait == "MeanLDMC" ~ "LDMC~('%')",
                           trait == "MeanLL" ~ "LL~(cm)",
                           trait == "MeanSLA" ~ "SLA~(cm^2/g)", 
                           trait == "C_N_ratio" ~ "C:N"))


#df with the annotations
annotations <- data.frame(trait = rep(c(unique(predictions_combo2$trait)),2), 
                          association = c(rep("bare", 7), rep("nurse", 7)))
annotations<- annotations[order(annotations$trait) , ]
annotations$t_test_significance <- c("  *", "  *", #C:N
                                     "  *", "  *", #H
                                     "  *", " ", #LA
                                     "  *", "  *", #LDMC
                                     "  *", "  *", #LL
                                     "  *", "  *", #LS
                                     "  *", "  *" )#SLA
annotations$anova_significance <- c("    a", "    b", #C:N
                                    "    a", "    b", #H
                                    "    a", "    b", #LA
                                    "    a", "    b", #LDMC
                                    " ", " ", #LL
                                    "    a", "    b", #LS
                                    "    a", "    b" )#SLA
annotations$ycoord_t <- c(6, 8, #C:N
                          160, 150, #H
                          3, 3, #LA
                          0.05, 0.2, #LDMC
                          5, 5, #LL
                          80000, 80000, #LS
                          25, 20)#SLA

annotations$ycoord_anova <- c(30, 30, #C:N
                              630, 630, #H
                              30, 30, #LA
                              0.60, 0.60, #LDMC
                              0, 0, #LL
                              800000, 800000, #LS
                              105, 105)#SLA


predicted_trait_distances <- ggplot(predictions_combo2, aes(x = association, y = trait_diff_predictions, fill = association)) +
  geom_boxplot(alpha = 0.6)+
  stat_summary(fun = mean, geom="point", shape = 23, size = 2, fill = "white", color = "black") +
  scale_fill_manual(values = c(brewer.pal(8, "Dark2")[7], brewer.pal(8, "Dark2")[1])) + 
  scale_x_discrete(labels = c(expression("∆"["Db"]), expression("∆"["Dd"]))) +
  facet_wrap(~trait, scales = "free_y", labeller = label_parsed) +
  ylab("Predicted trait difference") +
  xlab("Target species association") +
  geom_text(data = annotations, aes(x = association, y = ycoord_t, label = t_test_significance), size = 8, color = "brown3")+
  geom_text(data = annotations, aes(x = association, y = ycoord_anova, label = anova_significance), size = 5, color = "brown3")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = (element_text(size = 11)))

ggsave("model_predictions_one_dimensional_trait_distances.png", predicted_trait_distances, path = "Figures", 
       height = 1700, width = 2000, units = "px")  


###Functional distance ~ GRAZ, ARIDITY, microsite affinty####
#import data
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)

#the model
twosp_dist$arid_sq <- (twosp_dist$ARIDITY.v3)^2
twosp2 <- glmmTMB(euclidean_dist ~ grouping*GRAZ +grouping*ARIDITY.v3 + 
                    ARIDITY.v3*GRAZ + grouping*arid_sq + arid_sq*GRAZ + (1|SITE_ID), data = twosp_dist)
Anova(twosp2, type = 2) #type 2 is seq SS
summary(twosp2)

cld(glht(model = twosp2, mcp(GRAZ = "Tukey")))
emmeans(twosp2, specs = "GRAZ")

##distance~grouping, colour by graz
#we want to know which levels of graz are significantly different from each other, within gspecific groupings
#concateante graz and grouping to create treatemnt ID
twosp_dist <- twosp_dist |> 
  mutate(treatmentID = str_c(grouping, GRAZ, sep = "_"))
twosp_dist$treatmentID <- as.factor(twosp_dist$treatmentID)

treatment_mod <- lm(euclidean_dist ~ treatmentID, data = twosp_dist)
Anova(treatment_mod)
letters <- cld(glht(model = treatment_mod, mcp(treatmentID = "Tukey")))
#look at letters and make a vector to of them to add to the graph
#put them in the same order as the boxplots on the graph
# Create a data frame with the letters and corresponding x positions
letters_data <- data.frame(grouping = rep(unique(twosp_dist$grouping), each = length(unique(twosp_dist$GRAZ))),
                           GRAZ = rep(as.factor(c(0,1,2,3))), times = length(unique(twosp_dist$grouping)),
                           letter = c("a", "bc", "cde", "ef", "ag", "bc", "ef", "ce", "af", "b", "bd", "efg"), 
                           ycoord = c(rep(11.5, 12)))

dist_bygrouping_graz <- ggplot(twosp_dist, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(labels = c("ungrazed", "low", "medium", "high"), 
                    values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  ylim(0, 12) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Microsite affinity of target species", 
       fill = "") +
  scale_x_discrete(labels = c("bare", "dominant", "both")) +
  geom_text(data = letters_data, aes(x = grouping, y = ycoord, label = letter, group = GRAZ),
            position = position_dodge(width = 0.75), vjust = -0.5, hjust = 0.5, size = 4, color = "black") +
  theme_classic() +
  theme(legend.position = "bottom") 

ggsave("dist_grouping_graz.png", dist_bygrouping_graz, height = 1500, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#distance ~graz
twosp1 <- glmmTMB(euclidean_dist ~ grouping +GRAZ + ARIDITY.v3 + arid_sq 
                  + (1|SITE_ID), data = twosp_dist)
summary(twosp1)

letters <- cld(glht(model = twosp1, mcp(GRAZ = "Tukey")))
letters_df <- data.frame(GRAZ = c("0", "1", "2", "3"), letters = c("a", "b", "c", "d"), ycoord = c(12, 12, 12, 12))


dist_graz<- ggplot(twosp_dist, aes(x = GRAZ, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  ylim(0, 12) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Grazing pressure", 
       fill = "") +
  scale_x_discrete(labels = c("ungrazed", "low", "medium", "high")) +
  geom_text(data = letters_df, aes(x = GRAZ, y = ycoord, label = letters)) +
  theme_classic() +
  theme(legend.position = "none") 
ggsave("dist_graz.png", dist_graz, 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


##distance~aridity
#Get model predictions
model_forgraph <- glmmTMB(euclidean_dist ~ grouping*ARIDITY.v3*arid_sq, data = twosp_dist)
Anova(model_forgraph, type = 2)

pred_df <- data.frame(ARIDITY.v3 = c(rep(unique(twosp_dist$ARIDITY.v3), 12)), 
                      arid_sq = c(rep(unique(twosp_dist$arid_sq), 12)),
                      grouping = as.factor(c(rep("nurse_bare_only", length(unique(twosp_dist$ARIDITY.v3))*4), 
                                   rep("nurse_nurse_only", length(unique(twosp_dist$ARIDITY.v3))*4), 
                                   rep("nurse_both", length(unique(twosp_dist$ARIDITY.v3))*4)))) 


twosp2_predictions <- predict(model_forgraph, pred_df)
twosp2_predictions <- pred_df |> 
  add_column(predictions = c(twosp2_predictions))

#dist~aridity
line_color <- brewer.pal(8, "Dark2")[7]
dist_aridity <- ggplot(twosp_dist, aes(x = ARIDITY.v3, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  #geom_smooth(method = "lm", color = line_color)+
  #add model predictions to the plot
  geom_line(data = twosp2_predictions, aes(x = ARIDITY.v3, y = predictions), lty = 1, lwd = 1, color = line_color) +
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity") +
  theme_classic()

ggsave("dist_aridity.png", dist_aridity, height = 1400, width = 1100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#dist~aridity squared
line_color <- brewer.pal(8, "Dark2")[7]
dist_aridsq <- ggplot(twosp_dist, aes(x = arid_sq, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  #geom_smooth(method = "lm", color = line_color)+
  #add model predictions to the plot
  geom_line(data = twosp2_predictions, aes(x = arid_sq, y = predictions), lty = 1, lwd = 1, color = line_color) +
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity^2") +
  theme_classic()

ggsave("dist_aridity.png", dist_aridity, height = 1400, width = 900, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


###ONLY FOR FACILITATIVE REPS: Functional distance ~ GRAZ, ARIDITY, microsite affinty####
#import data
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1) |> 
  mutate(aridsq = ARIDITY.v3^2)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)

nint_result <- 
  read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  filter(ID %in% c(twosp_dist$ID)) |>
  filter(!is.na(NIntc_richness)) |> 
  select(!c(NIntc_cover, NInta_richness, NInta_cover, NIntc_shannon, NInta_shannon)) |> 
  mutate(ID_rep = str_c(ID, replicate_no, sep = "_")) |> 
  mutate(interaction = case_when(NIntc_richness > 0 ~ "facilitation", 
                                 NIntc_richness < 0 ~ "competition", 
                                 NIntc_richness == 0 ~ "neutral")) |> 
  select(ID_rep, interaction)

#now isolate facilitative reps in twosp_dist
twosp_dist_fac <- twosp_dist |> 
  mutate(ID_rep = str_c(ID, replicate, sep = "_")) |> 
  left_join(nint_result, by = "ID_rep")  |> 
  filter(interaction == "facilitation")


#we want to know which levels of graz are significantly different from each other, within specific groupings
#concateante graz and grouping to create treatemnt ID
twosp_dist_fac <- twosp_dist_fac |> 
  mutate(treatmentID = str_c(grouping, GRAZ, sep = "_"))
twosp_dist_fac$treatmentID <- as.factor(twosp_dist_fac$treatmentID)

fac_treatment_mod <- lm(euclidean_dist ~ treatmentID, data = twosp_dist_fac)
Anova(fac_treatment_mod)
fac_letters <- cld(glht(model = fac_treatment_mod, mcp(treatmentID = "Tukey")))
#look at letters and make a vector to of them to add to the graph
#put them in the same order as the boxplots on the graph
# Create a data frame with the letters and corresponding x positions
fac_letters_data <- data.frame(grouping = rep(unique(twosp_dist_fac$grouping), each = length(unique(twosp_dist_fac$GRAZ))),
                           GRAZ = rep(as.factor(c(0,1,2,3))), times = length(unique(twosp_dist_fac$grouping)),
                           letter = c("a", "bc", "cd", "de", "ef", "bc", "bd", "cd", "af", "c", "bc", "df"), 
                           ycoord = c(rep(11.5, 12)))


#dist~grouping, colour by graz
FAC_dist_bygrouping_graz <- ggplot(twosp_dist_fac, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(labels = c("ungrazed", "low", "medium", "high"), 
                    values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  ylim(0, 12) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Microsite affinity of target species", 
       fill = "") +
  scale_x_discrete(labels = c("bare", "dominant", "both")) +
  geom_text(data = fac_letters_data, aes(x = grouping, y = ycoord, label = letter, group = GRAZ),
            position = position_dodge(width = 0.75), vjust = -0.5, hjust = 0.5, size = 4, color = "black") +
  theme_classic() +
  theme(legend.position = "bottom") 
ggsave("facilitative_reps_dist_grouping_graz.png", FAC_dist_bygrouping_graz, height = 1500, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")

##distance~aridity
#Get model predictions
fac_model_forgraph <- glmmTMB(euclidean_dist ~ grouping*ARIDITY.v3*aridsq, data = twosp_dist_fac)
Anova(fac_model_forgraph, type = 2)

pred_df <- data.frame(ARIDITY.v3 = c(rep(unique(twosp_dist_fac$ARIDITY.v3), 12)), 
                      aridsq = c(rep(unique(twosp_dist_fac$aridsq), 12)),
                      grouping = as.factor(c(rep("nurse_bare_only", length(unique(twosp_dist_fac$ARIDITY.v3))*4), 
                                             rep("nurse_nurse_only", length(unique(twosp_dist_fac$ARIDITY.v3))*4), 
                                             rep("nurse_both", length(unique(twosp_dist_fac$ARIDITY.v3))*4)))) 


twosp_fac_predictions <- predict(fac_model_forgraph, pred_df)
twosp_fac_predictions <- pred_df |> 
  add_column(predictions = c(twosp_fac_predictions))


#dist~aridity
line_color <- brewer.pal(8, "Dark2")[7]
fac_dist_aridity <- ggplot(twosp_dist_fac, aes(x = ARIDITY.v3, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  #geom_smooth(method = "lm", color = line_color)+
  #add model predictions to the plot
  geom_line(data = twosp_fac_predictions, aes(x = ARIDITY.v3, y = predictions), lty = 1, lwd = 1, color = line_color) +
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity") +
  theme_classic()

ggsave("facilitative_reps_dist_aridity.png", fac_dist_aridity, height = 1400, width = 1100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")




###Biplot####
#import standardised FT data to use in pca
FT_wide_FGR <- read.csv("Functional trait data//results//standard_FT_FGR.csv", row.names = 1)
FT_wide_FGR$Functional_group <- as.factor(FT_wide_FGR$Functional_group) 

base_pca <- princomp(FT_wide_FGR[ , -c(8,9)], scores = T)

pal <- brewer.pal(3, "Dark2")
#PC 1 and 2 with FGR
pc1_2 <- ggbiplot(base_pca, choices = c(1,2), 
                  groups = FT_wide_FGR$Functional_group, 
                  ellipse = FALSE, 
                  var.axes = T) +
  scale_color_manual(labels = c("1 - Large plants", "2 - Tough plants", "3 - Palatable plants"), 
                     values = pal)+
  labs(x = "PC 1 (29.5%)", y = "PC 2 (22.2%)", color = "Functional groups") +
  xlim(-3, 5)+
  ylim(-3, 5 )+
  theme_classic()

ggsave("pc1_2_biplot.png", pc1_2, height = 1500, width = 2100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#PC 1 and 2 without FGR
pc1_2_nogroups <- ggbiplot(base_pca, choices = c(1,2), var.axes = T, alpha = 0.5, 
                  varname.color = "darkolivegreen4", varname.adjust = 1.2, varname.size = 6) +
  labs(x = "PC 1 (29.5%)", y = "PC 2 (22.2%)") +
  xlim(-3, 5)+
  ylim(-3, 5 )+
  theme_classic()

ggsave("pc1_2_biplot_nogroups.png", pc1_2_nogroups, height = 1500, width = 2100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#PC 1 and 3
pc1_3 <- ggbiplot(base_pca, choices = c(1,3), 
                  groups = FT_wide_FGR$Functional_group, 
                  ellipse = FALSE) +
  scale_color_manual(labels = c("1 - Large plants", "2 - Tough plants", "3 - Palatable plants"), 
                     values = pal)+
  labs(x = "PC 1 (29.5%)", y = "PC 3 (17.6%)", color = "Functional groups") +
  xlim(-3,4)+
  theme_classic()

ggsave("pc1_3_biplot.png", pc1_3,height = 1500, width = 2100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")



####Trait differences between species associations####
#import data 
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv", row.names = 1) |> 
  pivot_wider(names_from = trait, values_from = value) |> 
  #calculate the C:N ratio
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN)) |> 
  pivot_longer(cols = c("MeanLL","MeanSLA","MeanLDMC","MeanLA","MaxH","MaxLS","C_N_ratio"),
               names_to = "trait", values_to = "value")

##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(taxon= species)

FT_ass_join <- FT |> 
  left_join(ass, by = c("ID", "taxon")) |>  #ass is missing species in FT, because ass is missing the dominant species
  filter(!association %in% c("neutral", "too_rare")) |> 
  mutate(association = case_when(association == "nurse" ~ "nurse_associated", 
                                 association == "bare" ~ "bare_associated", 
                                 .default = association)) 

##Now we need to add values indicitaing which species are the dominant species 
#import the species position data
sp_positions <- read.csv("Functional trait data//Clean data//sp_positions.csv", row.names = 1) |> 
  select(!replicate) |> 
  distinct(ID, taxon, position) |> 
  filter(position == "nurse_species")

#overwrite the association column in FT_ass_join with "nurse_species" if it is a nurse in that plot
IDlist <- c(unique(sp_positions$ID))
FT_3_ass <- FT_ass_join

for(i in 1:length(IDlist)) {
  
  sp_positions_plot <- sp_positions[which(sp_positions$ID == IDlist[i]) , ]
  nurse_taxa <- c(sp_positions_plot$taxon)
  
  for(t in 1:length(nurse_taxa)) {
    FT_3_ass[which(FT_ass_join$ID == IDlist[i] & FT_ass_join$taxon == nurse_taxa[t]) 
                , which(colnames(FT_ass_join) == "association")] <- "nurse_species"
  }
}
#final cleaning of FT_ass_join
FT_3_ass <- FT_3_ass |> 
  filter(!is.na(value)) |> 
  mutate(trait = case_when(trait == "MaxH" ~ "H~(cm)", #rename the traits so that they are labelled nicely in the plot
                           trait == "MaxLS" ~ "LS~(cm^2)",
                           trait == "MeanLA" ~ "LA~(cm^2)",
                           trait == "MeanLDMC" ~ "LDMC~('%')",
                           trait == "MeanLL" ~ "LL~(cm)",
                           trait == "MeanSLA" ~ "SLA~(cm^2/g)", 
                           trait == "C_N_ratio" ~ "C:N"))
FT_3_ass$SITE_ID <- as.factor(FT_3_ass$SITE_ID)  
FT_3_ass$association <- as.factor(FT_3_ass$association)

##NOw we can make the figure###

#dataframe containing significance letters
annotations <- data.frame(trait = c(rep(unique(FT_3_ass$trait), 3)), 
                          association = c(rep("bare_associated", 7), rep("nurse_associated", 7), rep("nurse_species", 7)))
annotations <- annotations[order(annotations$trait) , ]
#add significance letters from analysis script
annotations$letters <- c("", "", "",  #cn ratio
                         "a", "b", "c", #maxh
                         "", "", "", #meanla
                         "", "", "", #meansldmc
                         "", "", "", #meanll
                         "a", "a", "b", #maxls
                         "a", "ab", "b" #meansla
                         )
annotations$ycoord <- c(rep(0,3), rep(560,3), rep(0,3), rep(0,3), rep(0, 3), rep(600000,3),rep(280,3))

#dataframe containing p values
p_vals <- data.frame(trait = c(unique(FT_3_ass$trait)), 
                     p_value = c( "", #LL
                                  "               *", #SLA
                                  "", #LDMC
                                  "", #LA
                                  "            ***", #H
                                  "            ***", #LS
                                  ""),         #C:N
                     ycoord = c(0,300,0,0,585,670000, 0))


trait_differences <- ggplot(FT_3_ass, aes(x = association, y = value)) +
  geom_boxplot(fill = "darkslategrey", alpha = 0.6) +
  facet_wrap(~trait, scale = "free_y", labeller = "label_parsed", ncol = 2) +
  ylab("Trait value") +
  xlab("") +
  scale_x_discrete(labels = c("bare-associated", "dominant-associated", "dominant species")) +
  geom_text(data = annotations, aes(y = ycoord, label = letters), color = "brown3", size = 8) +
  geom_text(data = p_vals, aes(x = "nurse_species",y = ycoord, label = p_value), color = "brown3", size = 8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18), 
        strip.text = element_text(size = 16))

ggsave("trait_differences.png", trait_differences, path = "Figures",  height = 3500, width = 3400, units = "px")
  
