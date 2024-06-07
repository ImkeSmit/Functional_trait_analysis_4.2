###Do nurse traits affect interaction outcomes??###
###Nintc ~ nurse traits###
library(tidyverse)
library(tidylog)
library(vegan)
library(DescTools)
library(glmmTMB)
library(car)
library(ggplot2)
library(DHARMa)
library(MuMIn)
library(corrplot)
library(ggpubr)

#import nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names = 1) 

##import trait data
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv",
               row.names = 1) 

##For each replicate, we need to get the NIntc value and the traits of the nurse
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
  filter(!is.na(NIntc_richness_binom), 
         !is.na(nurse_meanLA), 
         !is.na(nurse_mean_LS), 
         !is.na(nurse_mean_C_N_ratio), 
         !is.na(nurse_meanSLA), 
         !is.na(nurse_mean_H))

write.csv(modeldat_final, "Functional trait data//Clean data//nint_nurse_traits.csv")

#import the data that will be used for modelling
modeldat_final <- read.csv("Functional trait data//Clean data//nint_nurse_traits.csv", row.names = 1)
modeldat_final$nurse_sp <- as.factor(modeldat_final$nurse_sp)
modeldat_final$graz <- as.factor(modeldat_final$graz)
modeldat_final$site_ID <- as.factor(modeldat_final$site_ID)

cormat <- cor(modeldat_final[, which(colnames(modeldat_final) %like% "%mean%")], method = "pearson")
corrplot(cormat, method = "number")
#nothing is highly correlated except for N and CN ratio

##All variables are very left skewed except for CN ratio which is very right skewed.
hist(modeldat_final$nurse_mean_H)
hist(modeldat_final$nurse_meanLA)
hist(modeldat_final$nurse_meanSLA)
hist(modeldat_final$nurse_mean_LS)
hist(modeldat_final$nurse_mean_C_N_ratio)


###Now we can run the allsubsets models####
#import the model formulas
formula_table <- read.csv("Functional trait data\\results\\nint_nurse_trait_formulas.csv", row.names = 1) |> 
  separate_wider_delim(formula, delim = "~", names = c("response", "predictors")) |> 
  select(predictors) |> 
  distinct(predictors) |> 
  add_row(predictors = "1+(1|nurse_sp)+(1|site_ID)")

#import the data
modeldat_final <- read.csv("Functional trait data\\Clean data\\nint_nurse_traits.csv", row.names = 1) |> 
  mutate(aridity2 = aridity^2) |> 
  filter(!is.na(NIntc_richness_binom)) |> 
  filter(!is.na(NIntc_cover_binom)) |> 
  filter(!is.na(NInta_richness_binom)) |> 
  filter(!is.na(NInta_cover_binom))
modeldat_final$nurse_sp <- as.factor(modeldat_final$nurse_sp)
modeldat_final$graz <- as.factor(modeldat_final$graz)
modeldat_final$site_ID <- as.factor(modeldat_final$site_ID)

####descriptive stats####
#number of dominant species
length(unique(modeldat_final$nurse_sp)) #66

#number of replicates
nrow(modeldat_final) #2625

#%of replicates included in this analysis
nrow(modeldat_final)/nrow(nint_result) *100 #57.5


###Loop through the formulas for NIntc ~ nurse traits####
#Create a table for results
results_table <- data.frame(Response = character(), Model = character(), Chisq = numeric(), 
                            Df = integer(), Pr_value = numeric(), AIC = numeric(), 
                            Warnings = character(), row.names = NULL)

# Initialize warning_msg outside the loop
warning_msg <- ""

##Also loop through response variables
#loop through Nintc first
response_list <- c("NIntc_richness_binom", "NIntc_cover_binom", "NInta_richness_binom", "NInta_cover_binom")
datalist = c("modeldat_final", "modeldat_final", "modeldat_final", "modeldat_final")

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
    
    # Initialize anova_result and AIC_model outside the tryCatch block
    anova_result <- NULL
    AIC_model <- NULL
    
    tryCatch( #tryCatch looks for errors and warinngs in the expression
      expr = {
        model <- glmmTMB(formula, family = binomial, data = data)
        
        # Perform Anova 
        anova_result <- Anova(model, type = 2)
        # Get AIC
        AIC_model <- AIC(model)
        
        warning_messages <- warnings()
        
        ##Do nothing if the warinng is about non integer successes
        # Check for the non-integer #successes warning
        if ("non-integer #successes" %in% warning_messages) {
          # Handle non-integer #successes warning (e.g., print a message)
          message("Ignoring non-integer #successes warning")
        }
        
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
                             Chisq = ifelse(!is.null(anova_result), anova_result$Chisq[1], NA), 
                             Df = ifelse(!is.null(anova_result), anova_result$"Df"[1], NA), 
                             Pr_value = ifelse(!is.null(anova_result), anova_result$"Pr(>Chisq)"[1], NA), 
                             AIC = ifelse(!is.null(AIC_model), AIC_model, NA),
                             Warnings = warning_msg)
    
    
    results_table <- rbind(results_table, result_row)
  }
}

results_table
##LOOK AT WARNINGS
##The warnings are sort of cumulative. Once the is one warning of fitTMB, all warning_msg will contain it.
##I don't know how to fix this yet
##Look at the warning meassages printed. Remmber that every model will have the non integer warning.
##!!MOdels with NA in the CHisq have convergence problems. 

write.csv(results_table, "Functional trait data\\results\\nint_nurse_traits_model_results_20May2024.csv")


###select the model with the lowest AIC####
results_table <- read.csv("Functional trait data\\results\\nint_nurse_traits_model_results_17Apr2024.csv", row.names = 1)

best_subset_models <- results_table |> 
  filter(!is.na(AIC)) |> 
  group_by(Response) |> 
  top_n(-3, AIC) # get the 3 lowest AIC values
#the models selcted for nintc and ninta cover contain aridity2 but not aridity, so I guess we select the next lowest model


##Get the info of the best subset models##
#nintc richness null model:
nintc_richness_null <- glmmTMB(NIntc_richness_binom ~ 1 + (1|nurse_sp) +(1|site_ID), family = binomial, data = modeldat_final)

#NIntc richness best model
nintc_richness_bestmod <- glmmTMB(NIntc_richness_binom ~ aridity + nurse_meanSLA + nurse_mean_C_N_ratio + (1|nurse_sp) +(1|site_ID), 
                                 family = binomial, data = modeldat_final)

summary(nintc_richness_bestmod)
anova(nintc_richness_null, nintc_richness_bestmod) #p = 0.05172
Anova(nintc_richness_bestmod)
r.squaredGLMM(nintc_richness_bestmod)
#model diagnostics
nintc_richness_bestmod_simres <- simulateResiduals(nintc_richness_bestmod)
plot(nintc_richness_bestmod_simres)#HOV looks ok
#residuals underdispersed
testDispersion(simulateResiduals(nintc_richness_bestmod, re.form = NULL)) #dispersion test significant
testZeroInflation(nintc_richness_bestmod_simres) #less zeroes than expected


###

#nintc cover null model:
nintc_cover_null <- glmmTMB(NIntc_cover_binom ~ 1 + (1|nurse_sp) +(1|site_ID), family = binomial, data = modeldat_final)

#NIntc cover best model
nintc_cover_bestmod <- glmmTMB(NIntc_cover_binom ~ aridity+ nurse_meanSLA+ nurse_meanLA+ nurse_mean_C_N_ratio + (1|nurse_sp) +(1|site_ID), 
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
testZeroInflation(nintc_cover_bestmod_simres) #less zeroes than expected

###

#ninta richness null model:
ninta_richness_null <- glmmTMB(NInta_richness_binom ~ 1 + (1|nurse_sp) +(1|site_ID), family = binomial, data = modeldat_final)

#NInta richness best model
ninta_richness_bestmod <- glmmTMB(NInta_richness_binom ~ nurse_meanSLA+ nurse_mean_C_N_ratio + (1|nurse_sp) +(1|site_ID), 
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

###

#ninta cover null model:
ninta_cover_null <- glmmTMB(NInta_cover_binom ~ 1 + (1|nurse_sp) +(1|site_ID), family = binomial, data = modeldat_final)

#NInta cover best model
ninta_cover_bestmod <- glmmTMB(NInta_cover_binom ~  aridity+ nurse_meanSLA+ nurse_mean_C_N_ratio+ aridity2 + (1|nurse_sp) +(1|site_ID), 
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


####Let's make some figures####

###Nintc richness ~ SLA
nintc_richness_SLA <- ggplot(modeldat_final, aes(x = nurse_meanSLA, y = NIntc_richness)) +
  geom_jitter(width = 2, height = 0.1, alpha = 0.6, size = 2, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[C]~richness)) +
  xlab("mean SLA of dominant plant")

###NIntc_richness ~ C:N ratio
nintc_richness_CN <- ggplot(modeldat_final, aes(x = nurse_mean_C_N_ratio, y = NIntc_richness)) +
  geom_jitter(width = 2, height = 0.1, alpha = 0.6, size = 2, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[C]~richness)) +
  xlab("mean C:N ratio of dominant plant")

###Combine above figures
nintc_richness_var <- ggarrange(nintc_richness_SLA, nintc_richness_CN, labels = c("a", "b"))


###Nintc cover ~ SLA
#first get the model predictions of nintc cover ~ SLA
nintc_cover_SLA_mod <- glmmTMB(NIntc_cover_binom ~ nurse_meanSLA, data = modeldat_final, family = binomial)
pred_data <- data.frame(nurse_meanSLA = c(unique(modeldat_final$nurse_meanSLA)))
pred_data$binom_prediction <- predict(nintc_cover_SLA_mod, pred_data, type = "response")
pred_data$true_prediction <- 2*pred_data$binom_prediction -1

nintc_cover_SLA <- ggplot(modeldat_final, aes(x = nurse_meanSLA, y = NIntc_cover)) +
  geom_jitter(width = 2, height = 0.1, alpha = 0.6, size = 2, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[C]~cover)) +
  xlab("mean SLA of dominant plant") +
  geom_line(data = pred_data, aes(x = nurse_meanSLA, y = true_prediction), color = "turquoise", lwd = 1.5)


###NIntc cover ~ C:N ratio
#first get the model predictions of nintc cover ~ C:N
nintc_cover_CN_mod <- glmmTMB(NIntc_cover_binom ~ nurse_mean_C_N_ratio, data = modeldat_final, family = binomial)
pred_data <- data.frame(nurse_mean_C_N_ratio = c(unique(modeldat_final$nurse_mean_C_N_ratio)))
pred_data$binom_prediction <- predict(nintc_cover_CN_mod, pred_data, type = "response")
pred_data$true_prediction <- 2*pred_data$binom_prediction -1

nintc_cover_CN <- ggplot(modeldat_final, aes(x = nurse_mean_C_N_ratio, y = NIntc_cover)) +
  geom_jitter(width = 5, height = 0.1, alpha = 0.6, size = 2, colour = "darkslategrey") +
  theme_classic() +
  ylab(expression(NInt[C]~cover)) +
  xlab("mean C:N ratio of dominant plant") +
  geom_line(data = pred_data, aes(x = nurse_mean_C_N_ratio, y = true_prediction), color = "green", lwd = 1.5)




