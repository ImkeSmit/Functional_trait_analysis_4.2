###MAKE MODELS WITH NINT ~ RaoQ + graz + aridity####
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(ggplot2)
library(car)
library(corrplot)

FD_results <- read.csv("Functional trait data\\results\\FD_results_4Mar2024.csv", row.names = 1) 

FD_results$FRic <- as.numeric(FD_results$FRic)
FD_results$qual.FRic <- as.numeric(FD_results$qual.FRic)
FD_results$FEve <- as.numeric(FD_results$FEve)
FD_results$FDiv <- as.numeric(FD_results$FDiv)
FD_results$RaoQ <- as.numeric(FD_results$RaoQ) 
FD_results$ID <- as.factor(FD_results$ID)

##Import NIntc results
nint_result <- 
  read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  filter(ID %in% c(FD_results$ID))

##Histograms of predictors and responses
ggplot(nint_result, aes(x = NIntc_richness)) +
  geom_histogram()

ggplot(FD_results, aes(x = RaoQ)) + ##very very left skewed
  geom_histogram()

ggplot(FD_results, aes(x = log(RaoQ))) + ##log is better
  geom_histogram()

ggplot(FD_results, aes(x = sqrt(RaoQ))) + ##sqrt is slightly better
  geom_histogram()


##Correlations
cormat <- cor(FD_results[which(!is.na(FD_results$FEve)) , which(colnames(FD_results) %in% c("FRic", "FEve", "FDiv", "RaoQ"))],  method = "spearman")
corrplot(cormat, method = "number") ##correlations are all very low


#NIntc is bounded beween -1 and 1, so binomial family is appropriate
#However the function requires that the response be bounded between 0 and 1, so rescale NIntc
#x-min/max- min (here the formula is just already simplified)
nint_result$NIntc_richness_binom <- (nint_result$NIntc_richness + 1)/2
nint_result$NIntc_cover_binom <- (nint_result$NIntc_cover + 1)/2
nint_result$NIntc_shannon_binom <- (nint_result$NIntc_shannon + 1)/2

#x-min/max- min
nint_result$NInta_richness_binom <- (nint_result$NInta_richness - (-1)) / (2 - (-1))
nint_result$NInta_cover_binom <- (nint_result$NInta_cover - (-1)) / (2 - (-1))
nint_result$NInta_shannon_binom <- (nint_result$NInta_shannon - (-1)) / (2 - (-1))

nint_result$site_ID <- as.factor(nint_result$site_ID)
nint_result$graz <- as.factor(nint_result$graz)
nint_result$ID <- as.factor(nint_result$ID)


##summarise the nintc richness by plot
nintc_rich_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NIntc_richness_binom) |> 
  filter(!is.na(NIntc_richness_binom)) |> #remove NA's
  group_by(ID) |> #calculate mean and sd of NIntc richness binom in each plot
  mutate(mean_NIntc_rich_binom = mean(NIntc_richness_binom), 
         sd_NIntc_rich_binom = sd(NIntc_richness_binom), 
         n_obs = n()) |> 
  ungroup() |> 
  select(!NIntc_richness_binom) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  mutate(aridity2 = aridity^2) |> 
  filter(!is.na(RaoQ))

##summarise the ninta richness by plot
ninta_rich_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NInta_richness_binom) |> 
  filter(!is.na(NInta_richness_binom)) |> #remove NA's
  group_by(ID) |> #calculate mean and sd of NIntc richness binom in each plot
  mutate(mean_NInta_rich_binom = mean(NInta_richness_binom), 
         sd_NInta_rich_binom = sd(NInta_richness_binom), 
         n_obs = n()) |> 
  ungroup() |> 
  select(!NInta_richness_binom) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  mutate(aridity2 = aridity^2) |> 
  filter(!is.na(RaoQ))


##summarise nintc cover by plot
nintc_cov_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NIntc_cover_binom) |> 
  filter(!is.na(NIntc_cover_binom)) |> #remove NA's
  group_by(ID) |> #calculate mean and sd of NIntc richness binom in each plot
  mutate(mean_NIntc_cov_binom = mean(NIntc_cover_binom), 
         sd_NIntc_cov_binom = sd(NIntc_cover_binom), 
         n_obs = n()) |> 
  ungroup() |> 
  select(!NIntc_cover_binom) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  mutate(aridity2 = aridity^2) |> 
  filter(!is.na(RaoQ))

##summarise ninta cover by plot
ninta_cov_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NInta_cover_binom) |> 
  filter(!is.na(NInta_cover_binom)) |> #remove NA's
  group_by(ID) |> #calculate mean and sd of NIntc richness binom in each plot
  mutate(mean_NInta_cov_binom = mean(NInta_cover_binom), 
         sd_NInta_cov_binom = sd(NInta_cover_binom), 
         n_obs = n()) |> 
  ungroup() |> 
  select(!NInta_cover_binom) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  mutate(aridity2 = aridity^2) |> 
  filter(!is.na(RaoQ))


##Import formulas and get them ready to run models with
#lets isolate the response combinations
formula_table <- read.csv("Functional trait data\\Results\\model_formulas_with_RaoQ.csv", row.names = 1) |> 
  separate_wider_delim(formula, delim = "~", names = c("response", "predictors")) |> 
  select(predictors) |> 
  distinct(predictors) |> 
  add_row(predictors = "1+(1|site_ID)")  #add the null model
##!formula_table contains some nonsense models that contain aridity2 but not aridity!##
#These are too difficult to remove using partial string matches
#Just run all the nionsense models as well. 
#If they are selected with AIC, we just disregard them then



###Loop through the formulas####
#Create a table for results
results_table <- data.frame(Response = character(), Model = character(), Chisq = numeric(), 
                            Df = integer(), Pr_value = numeric(), AIC = numeric(), 
                            Warnings = character(), row.names = NULL)

# Initialize warning_msg outside the loop
warning_msg <- ""

##Also loop through response variables
#loop through Nintc first
response_list <- c("mean_NIntc_rich_binom", "mean_NIntc_cov_binom", "mean_NInta_rich_binom", "mean_NInta_cov_binom")
datalist = c("nintc_rich_sum", "nintc_cov_sum", "ninta_rich_sum", "ninta_cov_sum")

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

write.csv(results_table, "Functional trait data\\results\\nint_raoq_model_results_22Apr2024.csv")


###Interpret model results####
mod_results <- read.csv("Functional trait data\\results\\nint_raoq_model_results_22Apr2024.csv", row.names = 1)

##Which model had the lowest AIC?
mod_results |> 
  group_by(Response) |> 
  filter(!is.na(AIC)) |> #filter out records with model convergence problems
  filter(AIC == min(AIC))
##Null models have lowest AIC

##What are the 10 lowest AIC values
low_10_AIC <- mod_results |> 
  group_by(Response) |> 
  slice_min(AIC, n = 10) |> 
  arrange(AIC, .by_group = TRUE)

##There are models with similar AIC vals


###Are results different if we run it with log(RaoQ)####
log_nintc_rich_sum <- nintc_rich_sum |> 
  mutate(RaoQ = log(RaoQ))

log_nintc_cov_sum <- nintc_cov_sum |> 
  mutate(RaoQ = log(RaoQ))

log_ninta_rich_sum <- ninta_rich_sum |> 
  mutate(RaoQ = log(RaoQ))

log_ninta_cov_sum <- ninta_cov_sum |> 
  mutate(RaoQ = log(RaoQ))

###Loop through the  LOG formulas####
#Create a table for results
log_results_table <- data.frame(Response = character(), Model = character(), Chisq = numeric(), 
                            Df = integer(), Pr_value = numeric(), AIC = numeric(), 
                            Warnings = character(), row.names = NULL)

# Initialize warning_msg outside the loop
warning_msg <- ""

##Also loop through response variables
#loop through Nintc first
response_list <- c("mean_NIntc_rich_binom", "mean_NIntc_cov_binom", "mean_NInta_rich_binom", "mean_NInta_cov_binom")
datalist = c("log_nintc_rich_sum", "log_nintc_cov_sum", "log_ninta_rich_sum", "log_ninta_cov_sum")

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
    
    
    log_results_table <- rbind(log_results_table, result_row)
  }
}

#look for the best models
log_results_table |> 
  group_by(Response) |> 
  filter(!is.na(AIC)) |> #filter out records with model convergence problems
  filter(AIC == min(AIC))
###Results are the same!###


####Some graphs####
ggplot(nintc_rich_sum, aes(x = RaoQ, y = mean_NIntc_rich_binom)) +
  geom_point()+
  theme_classic()


ggplot(nintc_cov_sum, aes(x = RaoQ, y = mean_NIntc_cov_binom)) +
  geom_point()+
  theme_classic()


