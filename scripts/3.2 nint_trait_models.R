###Do nurse traits affect interaction outcomes??###
###Nintc ~ nurse traits###
library(tidyverse)
library(tidylog)
library(vegan)
#library(DescTools)
library(glmmTMB)
library(car)
library(ggplot2)
library(DHARMa)
library(MuMIn)
library(corrplot)
library(ggpubr)

###Get the traits of the dominant sp in each replicate####
#import nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names = 1) 

##import trait data
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv",
               row.names = 1) 

##For each replicate, we need to get the NIntc value and the traits of the nurse###
#names of the traits collectd
traits_collected <- c(unique(FT$trait))
#plot Id's
IDlist <- c(unique(nint_result$ID))
#empty table with explanatory variables
modeldat <- data.frame(ID = NA, site_ID = NA, replicate_no = NA, nurse_sp = NA, 
                       NIntc_richness = NA, NIntc_cover = NA, NInta_richness = NA, NInta_cover = NA, aridity = NA, graz = NA, 
                       nurse_mean_LL = NA, nurse_meanSLA = NA, nurse_meanLDMC = NA,nurse_meanLA = NA, 
                       nurse_mean_H = NA,nurse_mean_LS = NA, nurse_mean_percentN = NA, nurse_mean_percentC = NA)
#columns with trait values in modeldat
trait_mean_names <- c(colnames(modeldat)[11:18])

###Loop starts here
l = 1
for (i in 1:length(IDlist)) {
  
  #isolate a plot in the facilitation and FT data
  fac_plot <- nint_result[which(nint_result$ID == IDlist[i]) , ]
  FT_plot <- FT[which(FT$ID == IDlist[i]) , ]
  
  #list of reps in this plot
  replist <- c(unique(fac_plot$replicate_no))
  
  #for each rep:
  for (r in 1:length(replist)) {
    one_rep <- fac_plot[which(fac_plot$replicate_no == replist[r]) , ]
    #get the nurse species
    nurse_sp <- one_rep$nurse
    
    #get the trait values of the nurse species
    for (t in 1:length(traits_collected)) {
      val <- FT_plot |> 
        filter(taxon == nurse_sp, 
               trait == traits_collected[t]) |> 
        select(value)
      
      #if there are no values for this trait, put NA in the matrix
      if(nrow(val) == 0) {modeldat[l, which(colnames(modeldat) == trait_mean_names[t])] <- NA} 
      #if there are 2 or more trait measurements, get the mean and put that in the matrix
      else if (nrow(val) > 1) {mean_val <- mean(val$value)
      modeldat[l, which(colnames(modeldat) == trait_mean_names[t])] <- mean_val} 
      #if there is only on trait value, put that in the matrix
      else {modeldat[l, which(colnames(modeldat) == trait_mean_names[t])] <- val$value}
    }#loop through traits end
    
    #fill the rest of the table
    modeldat[l,1] <- IDlist[i]
    modeldat[l,2] <- one_rep$site_ID
    modeldat[l,3] <- replist[r]
    modeldat[l,4] <- nurse_sp
    modeldat[l,5] <- one_rep$NIntc_richness
    modeldat[l,6] <- one_rep$NIntc_cover
    modeldat[l,7] <- one_rep$NInta_richness
    modeldat[l,8] <- one_rep$NInta_cover
    modeldat[l,9] <- one_rep$aridity
    modeldat[l,10] <- one_rep$graz

    l = l+1
  }#loop through reps end
}#loop through plots end


#transform the nint to binomial 
modeldat_final <- modeldat |> 
  mutate(NIntc_richness_binom = (NIntc_richness + 1)/2, 
       NIntc_cover_binom = (NIntc_cover + 1)/2, 
       NInta_richness_binom = (NInta_richness - (-1)) / (2 - (-1)), 
       NInta_cover_binom = (NInta_cover - (-1)) / (2 - (-1)), 
       nurse_mean_C_N_ratio = nurse_mean_percentC/nurse_mean_percentN) |> 
  filter(!is.na(NIntc_richness_binom)) |> 
  #and log transform the traits
  mutate(log_nurse_meanLL = log10(nurse_mean_LL), 
         log_nurse_meanSLA = log10(nurse_meanSLA), 
         log_nurse_meanLDMC = log10(nurse_meanLDMC), 
         log_nurse_meanLA = log10(nurse_meanLA), 
         log_nurse_meanH = log10(nurse_mean_H), 
         log_nurse_meanLS = log10(nurse_mean_LS), 
         log_nurse_meanCNratio = log10(nurse_mean_C_N_ratio))

write.csv(modeldat_final, "Functional trait data//Clean data//nint_nurse_traits.csv")

###descriptive statistics of nint_nurse_traits####
modeldat_final <- read.csv("Functional trait data//Clean data//nint_nurse_traits.csv", row.names = 1)
modeldat_final$nurse_sp <- as.factor(modeldat_final$nurse_sp)
modeldat_final$graz <- as.factor(modeldat_final$graz)
modeldat_final$site_ID <- as.factor(modeldat_final$site_ID)

##examine raw distributions
LL_hist <- ggplot(modeldat_final, aes(x = nurse_mean_LL)) +geom_histogram()

SLA_hist <- ggplot(modeldat_final, aes(x = nurse_meanSLA)) +geom_histogram()

LDMC_hist <- ggplot(modeldat_final, aes(x = nurse_meanLDMC)) +geom_histogram()

LA_hist <- ggplot(modeldat_final, aes(x = nurse_meanLA)) +geom_histogram()

H_hist <- ggplot(modeldat_final, aes(x = nurse_mean_H)) +geom_histogram()

LS_hist <- ggplot(modeldat_final, aes(x = nurse_mean_LS)) +geom_histogram()

CN_hist <- ggplot(modeldat_final, aes(x = nurse_mean_C_N_ratio)) +geom_histogram()

raw_histograms <- ggarrange(LL_hist, SLA_hist, LDMC_hist, LA_hist, H_hist, LS_hist, CN_hist)
ggsave("raw_nurse_trait_histograms.png", raw_histograms, path = "Figures", width = 1800, height = 1400, units = "px")


##examine logged distributions
##examine raw distributions
LL_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanLL)) +geom_histogram()

SLA_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanSLA)) +geom_histogram()

LDMC_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanLDMC)) +geom_histogram()

LA_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanLA)) +geom_histogram()

H_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanH)) +geom_histogram()

LS_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanLS)) +geom_histogram()

CN_hist <- ggplot(modeldat_final, aes(x = log_nurse_meanCNratio)) +geom_histogram()

log_histograms <- ggarrange(LL_hist, SLA_hist, LDMC_hist, LA_hist, H_hist, LS_hist, CN_hist)
ggsave("log_nurse_trait_histograms.png", log_histograms, path = "Figures", width = 1800, height = 1400, units = "px")



##do correlations 
cordata <- modeldat_final |> 
  distinct(ID, nurse_sp, .keep_all = T) |> 
  select(contains("mean")) |>
  select(!contains("percent")) |> 
  select(!contains("log")) |> 
  na.omit() #remove all rows that have an NA in any column

png("Figures\\nurse_trait_correlation.png")

cormat <- cor(cordata, method = "pearson")
corrplot(cormat, method = "number", type = "lower")
dev.off()

##LS is highly correlated with H!

###Now we can run the allsubsets models####
#import the model formulas
formula_table <- read.table("Functional trait data\\results\\nint_nurse_trait_clim_soil_formulas.csv", sep = ",", header = T) |> 
  rename(predictors = formula) |> 
  mutate(predictors = paste(predictors, "(1|nurse_sp)+(1|site_ID/ID)", sep = "+")) |> #add the random effect to all formulas
  add_row(predictors = "1+(1|nurse_sp)+(1|site_ID/ID)")   #add the null model

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

hist(modeldat_final$AMT)
hist(modeldat_final$RASE)
hist(modeldat_final$aridity)
hist(modeldat_final$pH)
hist(modeldat_final$SAC)
hist(modeldat_final$RAI)
range(modeldat_final$RAI)

####descriptive stats####
#number of dominant species
length(unique(modeldat_final$nurse_sp)) #87

#number of replicates included in analysis
modeldat_final |> 
  filter(!is.na(nurse_meanSLA)) |> 
  filter(!is.na(nurse_meanLA)) |>
  filter(!is.na(nurse_mean_H)) |>
  filter(!is.na(nurse_mean_C_N_ratio)) |> 
  summarise(n = n()) #2625 complete cases

#%of replicates included in this analysis
#complete cases in nint_result
nint_result |> 
  filter(!is.na(NIntc_richness)) |> 
  filter(!is.na(NIntc_cover)) |> 
  summarise(n = n()) #3735

2625/3735 *100 #70.28112

#histograms of the env variables
##examine raw distributions
aridity_hist <- ggplot(modeldat_final, aes(x = aridity)) +geom_histogram()

AMT_hist <- ggplot(modeldat_final, aes(x = AMT)) +geom_histogram()

RASE_hist <- ggplot(modeldat_final, aes(x = RASE)) +geom_histogram()

SAC_hist <- ggplot(modeldat_final, aes(x = SAC)) +geom_histogram()

pH_hist <- ggplot(modeldat_final, aes(x = pH)) +geom_histogram()

raw_env_histograms <- ggarrange(aridity_hist, AMT_hist, RASE_hist, SAC_hist, pH_hist)
ggsave("raw_env_histograms.png", raw_env_histograms, path = "Figures", width = 1800, height = 1400, units = "px")


#examine logged distributions
a <- ggplot(modeldat_final, aes(x = log10(aridity))) +geom_histogram()

b <- ggplot(modeldat_final, aes(x = log10(AMT))) +geom_histogram()

c <- ggplot(modeldat_final, aes(x = log10(RASE))) +geom_histogram()

d <- ggplot(modeldat_final, aes(x = log10(SAC))) +geom_histogram()

e <- ggplot(modeldat_final, aes(x = log10(pH))) +geom_histogram()

log_env_histograms <- ggarrange(a,b,c,d,e)
ggsave("log_env_histograms.png", log_env_histograms, path = "Figures", width = 1800, height = 1400, units = "px")

###Loop through the formulas for NIntc ~ nurse traits####
#Initialise output file for results
output_file <- "Functional trait data\\results\\nintA_nurse_traits_clim_soil_nestedRE_model_results_31Aug2024.csv"

# Initialize the output file
write.csv(data.frame(Response = character(), Model = character(), AIC = numeric(), BIC = numeric(), 
                     Warnings = character()), output_file, row.names = FALSE)

# Initialize warning_msg outside the loop
warning_msg <- ""

##Also loop through response variables
#loop through Nintc first
#response_list <- c("NIntc_richness_binom", "NIntc_cover_binom", "NInta_richness_binom", "NInta_cover_binom")
#datalist = c("modeldat_final", "modeldat_final", "modeldat_final", "modeldat_final")

response_list <- c("NInta_cover_binom")
datalist = c("modeldat_final")

# Initialize a list to store result rows
result_list <- list()

# Specify chunk size
chunk_size <- 150

##LOOP THROUGH MODELS STARTS HERE##
#Loop through response variables
for(r in 1:length(response_list)) {
  
  response_var <- response_list[r]  
  data = get(datalist[r])
  
  #Loop through response variables
  for (f in 1:nrow(formula_table)) {
    
    predictors <- as.character(formula_table[f, ])
    formula <- as.formula(paste(response_var, "~",  predictors))
    
    # Clear existing warning messages
    warnings()
    
    # Initialize AIC_model outside the tryCatch block
    AIC_model <- NULL
    BIC_model <- NULL
    
    tryCatch( #tryCatch looks for errors and warinngs in the expression
      expr = {
        model <- glmmTMB(formula, family = binomial, data = data)
        
        # Get AIC
        AIC_model <- AIC(model)
        BIC_model <- BIC(model)
        
        warning_messages <- warnings()
        
        ##Do nothing if the warinng is about non integer successes
        # Check for the non-integer #successes warning
        #if ("non-integer #successes" %in% warning_messages) {
          # Handle non-integer #successes warning (e.g., print a message)
         # message("Ignoring non-integer #successes warning")
        #}
        
        #Print the warning message if it is about model fit
        # Check for other warnings, excluding the non-integer #successes warning
        other_warnings <- setdiff(warning_messages, "non-integer #successes")
        if (length(other_warnings) > 0) {
          warning_msg <- paste("warning :", as.character(other_warnings), collapse = "; ")
          message(paste("WARNING_", "r =" , response_var, "f =", f, warning_msg))
        }
      }, 
      
      #Also show me errors
      error = function(e) {
        message(paste("ERROR_", "r =" , response_var, "f =", f, conditionMessage(e)))
        print(e)
      }
    )
    
    # Extract relevant information
    result_row <- data.frame(Response = response_var,
                             Model = paste(response_var, "~",  predictors), 
                             AIC = ifelse(!is.null(AIC_model), AIC_model, NA),
                             BIC = ifelse(!is.null(BIC_model), BIC_model, NA),
                             Warnings = warning_msg)
    
    # Add the result row to the list
    result_list[[length(result_list) + 1]] <- result_row
    
    # If chunk size is reached, write to the file and reset the list
    if (length(result_list) == chunk_size) {
      write.table(do.call(rbind, result_list), output_file, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
      result_list <- list()  # Reset the list
    }
  }
} 
# Write any remaining rows to the file after the loop
  if (length(result_list) > 0) {
    write.table(do.call(rbind, result_list), output_file, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
}

##LOOK AT WARNINGS
##The warnings are sort of cumulative. Once the is one warning of fitTMB, all warning_msg will contain it.
##I don't know how to fix this yet
##Look at the warning meassages printed. Remmber that every model will have the non integer warning.
##!!MOdels with NA in the CHisq have convergence problems. 

##Output is saved in two files, nintA_nurse_traits_nestedRE_clim_soil_model_results_24Aug2024, and nintC_nurse_traits_nestedRE_clim_soil_model_results_24Aug2024

##new run with nested RE (chunk size = 200)
#start nintC richness 14:03 23 Aug
#end nintc cover 20:00 25 Aug

##run nintA with nested RE (chink size = 150)
#only load nintA richness onto the loop
#start nintA richness at 16:33 on 16 Aug
#forgot to note end time

#start nintA cover on 19:46 on 3 Sept


###Post hoc tests on best models####
###select the model with the lowest AIC###
nintC_results <- read.csv("Functional trait data\\results\\nintC_nurse_traits_clim_soil_nestedRE_model_results_23Aug2024.csv", 
                          col.names = c("Response", "Model","AIC", "BIC", "Warnings")) |> 
  filter(Response %in% c("NIntc_richness_binom", "NIntc_cover_binom")) #remove some nintA results that snuck in before I stopped the loop
nintA_results <- read.csv("Functional trait data\\results\\nintA_nurse_traits_clim_soil_nestedRE_model_results_31Aug2024.csv")

best_subset_models <- nintC_results |> 
  bind_rows(nintA_results) |> 
  filter(!is.na(AIC)) |> 
  group_by(Response) |>
  filter(AIC == min(AIC))

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

###Nintc richness

#nintc richness null model:
nintc_richness_null <- glmmTMB(NIntc_richness_binom ~ 1 + (1|nurse_sp) +(1|site_ID/ID), family = binomial, data = modeldat_final)

#NIntc richness best model
nintc_richness_bestmod <- glmmTMB(NIntc_richness_binom ~  graz+RASE+SAC+log_nurse_meanLA+log_nurse_meanH+log_nurse_meanCNratio+graz:SAC+(1|nurse_sp)+(1|site_ID/ID), 
                                 family = binomial, data = modeldat_final)


summary(nintc_richness_bestmod)
anova(nintc_richness_null, nintc_richness_bestmod) #p = 0.05172
Anova(nintc_richness_bestmod)
r.squaredGLMM(nintc_richness_bestmod)
#model diagnostics
nintc_richness_bestmod_simres <- simulateResiduals(nintc_richness_bestmod)
plot(nintc_richness_bestmod_simres)#HOV violated
#residuals underdispersed
testDispersion(simulateResiduals(nintc_richness_bestmod, re.form = NULL)) #dispersion test significant
testZeroInflation(nintc_richness_bestmod_simres) #more zeroes than expected

emmeans(nintc_richness_bestmod, specs = "graz")

#how many reps included?
modeldat_final |> 
  filter(!is.na(NIntc_richness_binom)) |> 
  filter(!is.na(log_nurse_meanLA)) |>
  filter(!is.na(log_nurse_meanH)) |>
  filter(!is.na(log_nurse_meanCNratio)) |>
  summarise(n = n()) #2625

###Nintc cover

#nintc cover null model:
nintc_cover_null <- glmmTMB(NIntc_cover_binom ~ 1 + (1|nurse_sp) +(1|site_ID/ID), family = binomial, data = modeldat_final)

#NIntc cover best model
nintc_cover_bestmod <- glmmTMB(NIntc_cover_binom ~  graz+aridity+RASE+pH+SAC+
                                 log_nurse_meanLA+log_nurse_meanSLA+log_nurse_meanH+log_nurse_meanCNratio+
                                 graz:RASE+graz:pH+graz:SAC + (1|nurse_sp) +(1|site_ID/ID), 
                                  family = binomial, data = modeldat_final)

summary(nintc_cover_bestmod)
anova(nintc_cover_null, nintc_cover_bestmod) #p = 0.003207
Anova(nintc_cover_bestmod)
r.squaredGLMM(nintc_cover_bestmod)
#model diagnostics
nintc_cover_bestmod_simres <- simulateResiduals(nintc_cover_bestmod)
plot(nintc_cover_bestmod_simres)#HOV looks ok
#residuals underdiespersed
testDispersion(nintc_cover_bestmod_simres) #dispersion test significant
testZeroInflation(nintc_cover_bestmod_simres) #more zeroes than expected

emmeans(nintc_cover_bestmod, specs = "graz")

#how many reps included?
modeldat_final |> 
  filter(!is.na(NIntc_cover_binom)) |> 
  filter(!is.na(log_nurse_meanLA)) |>
  filter(!is.na(log_nurse_meanSLA)) |>
  filter(!is.na(log_nurse_meanH)) |>
  filter(!is.na(log_nurse_meanCNratio)) |>
  summarise(n = n()) #2625


###NIntA richness

#ninta richness null model:
ninta_richness_null <- glmmTMB(NInta_richness_binom ~ 1 + (1|nurse_sp) +(1|site_ID/ID), family = binomial, data = modeldat_final)

#NInta richness best model
ninta_richness_bestmod <- glmmTMB(NInta_richness_binom ~ graz+pH+SAC+log_nurse_meanLA+log_nurse_meanH+log_nurse_meanCNratio+graz:pH+graz:SAC+(1|nurse_sp)+(1|site_ID/ID), 
                                  family = binomial, data = modeldat_final)

summary(ninta_richness_bestmod)
anova(ninta_richness_null, ninta_richness_bestmod) #p = 0.04243
Anova(ninta_richness_bestmod)
r.squaredGLMM(ninta_richness_bestmod)
#model diagnostics
ninta_richness_bestmod_simres <- simulateResiduals(ninta_richness_bestmod)
plot(ninta_richness_bestmod_simres)#HOV looks ok
#residuals underdispersed
testDispersion(ninta_richness_bestmod_simres) #jup underdispersion
testZeroInflation(ninta_richness_bestmod_simres) #less zeroes than expected
emmeans(ninta_richness_bestmod, specs = "graz")

###NintA cover

#ninta cover null model:
ninta_cover_null <- glmmTMB(NInta_cover_binom ~ 1 + (1|nurse_sp) +(1|site_ID/ID), family = binomial, data = modeldat_final)

#NInta cover best model
ninta_cover_bestmod <- glmmTMB(NInta_cover_binom ~ graz+aridity+pH+SAC+log_nurse_meanLA+log_nurse_meanSLA+log_nurse_meanH+log_nurse_meanCNratio+graz:pH+graz:SAC+(1|nurse_sp)+(1|site_ID/ID), 
                               family = binomial, data = modeldat_final)

summary(ninta_cover_bestmod)
anova(ninta_cover_null, ninta_cover_bestmod) #p = 0.004049
Anova(ninta_cover_bestmod)
r.squaredGLMM(ninta_cover_bestmod)
#model diagnostics
ninta_cover_bestmod_simres <- simulateResiduals(ninta_cover_bestmod)
plot(ninta_cover_bestmod_simres)#HOV looks ok
#residuals underdispersed
testDispersion(ninta_cover_bestmod_simres) #jup underdispersed
testZeroInflation(ninta_cover_bestmod_simres) #less zeroes than expected

emmeans(ninta_cover_bestmod, specs = "graz")
