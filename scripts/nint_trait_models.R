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
library(corrplot)

#import nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names = 1) 

##import trait data
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv",
               row.names = 1) #|> 
  #standardise trait values
  #group_by(trait) |> 
  #mutate(sd_value = sd(value), 
  #       mean_value = mean(value)) |> 
  #ungroup() |> 
  #mutate(value_std = (value - mean_value)/sd_value)

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

modeldat_final$nurse_sp <- as.factor(modeldat_final$nurse_sp)
modeldat_final$graz <- as.factor(modeldat_final$graz)
modeldat_final$site_ID <- as.factor(modeldat_final$site_ID)

cormat <- cor(modeldat_final[, which(colnames(modeldat_final) %like% "%mean%")])
corrplot(cormat, method = "number")
#nothing is gighly correlated except for N and CN ratio


##Now we can do the model building
#import the model formulas
formula_table <- read.csv("Functional trait data\\results\\nint_nurse_trait_formulas.csv") |> 
  separate_wider_delim(formula, delim = "~", names = c("response", "predictors")) |> 
  select(predictors) |> 
  distinct(predictors) |> 
  add_row(predictors = "1+(1|nurse_sp)+(1|site_ID)")
#import the data
modeldat_final <- read.csv("Functional trait data\\Clean data\\nint_nurse_traits.csv")

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
datalist = c("modeldat_final")

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

write.csv(results_table, "Functional trait data\\results\\nint_nurse_traits_model_results_17Apr2024.csv")



##Now we can run the model
nint_richness_null <- glmmTMB(NIntc_richness_binom ~ 1 + (1|nurse_sp) +(1|site_ID), family = binomial, data = modeldat_final)

nintc_richness_model <- 
  glmmTMB(NIntc_richness_binom ~ nurse_meanLA + nurse_mean_LS + nurse_mean_H + nurse_mean_C_N_ratio + nurse_meanSLA + graz + aridity + (1|nurse_sp) +(1|site_ID), 
          family = binomial, data = modeldat_final)

summary(nintc_richness_model)
Anova(nintc_richness_model)
anova(nint_richness_null, nintc_richness_model)

ggplot(modeldat_final, aes(x = nurse_meanLA, y = NIntc_richness)) +
  geom_point() +
  theme_classic()

ggplot(modeldat_final, aes(x = nurse_meanSLA, y = NIntc_richness)) +
  geom_point() +
  theme_classic()

ggplot(modeldat_final, aes(x = nurse_mean_LS, y = NIntc_richness)) +
  geom_point() +
  theme_classic()

ggplot(modeldat_final, aes(x = nurse_mean_H, y = NIntc_richness)) +
  geom_point() +
  theme_classic()

ggplot(modeldat_final, aes(x = nurse_mean_C_N_ratio, y = NIntc_richness)) +
  geom_point() +
  theme_classic()


##All variables are very left skewed except for CN ratio which is very right skewed.
hist(modeldat_final$nurse_mean_H)
hist(modeldat_final$nurse_meanLA)
hist(modeldat_final$nurse_meanSLA)
hist(modeldat_final$nurse_mean_LS)
hist(modeldat_final$nurse_mean_C_N_ratio)

plotResiduals(nintc_richness_model)
