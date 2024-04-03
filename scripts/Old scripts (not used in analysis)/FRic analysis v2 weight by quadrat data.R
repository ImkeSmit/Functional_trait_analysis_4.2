###FUNCTIONAL RICHNESS ANALYSIS###
#Get the functional richness of each plot with FD
#Then Build models for Nint ~ Fric

library(tidyverse)
library(tidylog)
library(traitstrap)
library(FD)
library(glmmTMB)
library(car)
library(DescTools)
library(corrplot)

##Import filled FT data for the facilitation plots
FT_raw <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots.csv", row.names = 1) 

FT <- FT_raw|> 
  mutate(taxon = str_replace(taxon, " ", "_")) |> 
  pivot_wider(names_from = trait, values_from = value) |> 
  #calculate the C:N ratio
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN)) |> 
  pivot_longer(cols = c("MeanLL","MeanSLA","MeanLDMC","MeanLA","MaxH","MaxLS","C_N_ratio"),
               names_to = "trait", values_to = "value")

#get species in the FT data
FT_sp <- FT |> 
  distinct(taxon)


##Import quadrat survey for the faciliation plots
quad <- read.csv("Functional trait data\\quadrat_survey_for_facilitation_plots.csv", row.names = 1) |> 
  mutate(taxon = str_replace(taxon, " ", "_"))
length(c(unique(quad$taxon)))
#there are more sp in quad than in FT. Remove unmatched species from quad because they do not have trait values.

#FD also requires a community species by site matrix
#get the poltelevel cover of each species (this is the same as percent_cover_perplot in FT)
quad_wide <- quad |> 
  #remove species that are not in the FT data
  filter(taxon %in% c(FT_sp$taxon)) |> 
  mutate(quadrat = as.numeric(quadrat)) |> 
  #summarise the percent cover for each sp per plot, not per quadrat
  group_by(ID, taxon) |> 
  mutate(sum_plot_cover = sum(percent_cover)) |>  
  ungroup(taxon) |> 
  mutate(max_quadrats = max(quadrat)) |> 
  ungroup () |> 
  mutate(percent_cover_perplot = sum_plot_cover/max_quadrats) |> 
  select(!c(quadrat, percent_cover, sum_plot_cover, max_quadrats, site_plot)) |> 
  distinct() |> 
  #change to wide
  pivot_wider(names_from = taxon, values_from = percent_cover_perplot) |> 
  column_to_rownames(var = "ID")

sum(quad_wide[which(!is.na(quad_wide$Ajuga_iva)) , 1]) #0.29
ncol(quad_wide)

#replace NA's in quad_wide with 0
for (r in 1:nrow(quad_wide)) {
  for(c in 1:ncol(quad_wide)) {
    
    record <- quad_wide[r , c]
    
    if(is.na(record)) {
      quad_wide[r , c] <- 0
    }
  }
}
sum(quad_wide[ , 1]) #0.29


###Create separate species x trait matrices for each plot, put them in a list####
#get the mean trait value where there are multiple sp entries in a plot
FT_mean <- FT |> 
  filter(!is.na(value)) |> 
  group_by(ID, taxon, trait) |> 
  summarise(mean_value = mean(value)) |> 
  ungroup()

IDlist <- c(unique(FT_mean$ID))

xlist <- vector(mode='list', length = length(IDlist))
names(xlist) <- IDlist

alist <- vector(mode='list', length = length(IDlist))
names(alist) <- IDlist

for(i in 1:length(IDlist)) {
  plot_wide <- FT_mean |> 
    filter(ID == IDlist[i]) |> 
    select(!ID) |> 
    pivot_wider(names_from = trait, values_from = mean_value) |> 
    column_to_rownames(var = "taxon")
  
  FT_plot_names <- c(row.names(plot_wide))
  
  xlist[[i]] <- plot_wide
  
  #subset quad_wide for the  ID and the sp in plot_wide
  #because dbfd needs the sp between the traits and teh community to correspond exactly
  quad_sub <- quad_wide |> 
    filter(row.names(quad_wide) == IDlist[i]) |> 
    select(any_of(FT_plot_names))
  
  alist[[i]] <- quad_sub
}


###Quality assesment - do we have enough species?####
IDlist <- c(row.names(quad_wide))

nsp_dat <- data.frame(ID = IDlist , trait_nsp = c(rep(NA, length(IDlist))) , 
                      quad_nsp = c(rep(NA, length(IDlist))), 
                      trait_percent_cover = c(rep(NA, length(IDlist))), 
                      quad_percent_cover = c(rep(NA, length(IDlist))))
l = 1
for (t in 1:length(IDlist)) {
  
  ##number of species in FT data
  trait_plot <- xlist[[name = IDlist[t]]]
  trait_nsp <- nrow(trait_plot)
  
  ##number of sp in quad data
  quad_nsp <- quad |> #use quad before we made it wide and removed extra species
    filter(ID == IDlist[t]) |> 
    distinct(taxon) |> 
    nrow()
  
  ##added cover of all species in FT data
  trait_cov <- FT |> 
    filter(ID == IDlist[t]) |> 
    distinct(taxon, percent_cover_perplot) |> 
    summarise(trait_percent_cover = sum(percent_cover_perplot))
  
  ##added cover of all species in quad data
  quad_cov <- quad |> #use quad before we made it wide and removed extra species
    filter(ID == IDlist[t]) |> 
    group_by(quadrat) |> 
    mutate(sum_quadrat_cover = sum(percent_cover)) |> 
    ungroup() |>
    distinct(quadrat, sum_quadrat_cover) |> 
    summarise(sum_plot_cover = sum(sum_quadrat_cover), max_quadrats = max(quadrat), 
              plot_percent_cover = sum_plot_cover/max_quadrats)
  
  nsp_dat[l, 2] <- trait_nsp
  nsp_dat[l, 3] <- quad_nsp
  nsp_dat[l, 4] <- trait_cov
  nsp_dat[l, 5] <- quad_cov[, 3]
  
  l = l+1
}

nsp_dat <- nsp_dat |> 
  mutate(percent_sp_in_FT = (trait_nsp/quad_nsp)*100, 
         trait_cover_div_quad_cover = (trait_percent_cover/quad_percent_cover)*100)


#plots where we have data for more than 80% cover
plots_for_FD <- nsp_dat |> 
  filter(trait_cover_div_quad_cover >= 80) |> 
  distinct(ID)

##Remove the plots with less than 80% cover from xlist and alist
remove_plots <- nsp_dat |> 
  filter(trait_cover_div_quad_cover <= 80) |> 
  distinct(ID)

xlist_red <- xlist[-which(names(xlist) %in% c(as.character(remove_plots$ID)))]

alist_red <- alist[-which(names(alist) %in% c(as.character(remove_plots$ID)))]


###Now calculate FD for each plot
#Table to put results in 
column_headers <- c("ID", "nsp", "FRic", "qual.FRic", "FEve", "FDiv", "FDis", "RaoQ")
FD_results <- data.frame(matrix(nrow = length(IDlist), ncol = length(column_headers)))
colnames(FD_results) <- column_headers

IDlist <- c(plots_for_FD$ID)

##plot 115 breaks because Ferula sink has many NA's. Let's remove this species.
new_115 <- xlist_red[[which(names(xlist_red) == "115")]] |> 
  filter(!row.names(xlist_red[[which(names(xlist_red) == "115")]]) == "Ferula_sinkiangensis")
xlist_red[[which(names(xlist_red) == "115")]] <- new_115

new_cov_115 <- alist_red[[which(names(alist_red) == "115")]] |> 
  select(!Ferula_sinkiangensis)
alist_red[[which(names(alist_red) == "115")]]  <- new_cov_115

for (k in 1:length(xlist_red)) {
  
  nsp <- nrow(xlist_red[[k]])
  
  if(nsp > 1) {
    
    result <- dbFD(x = xlist_red[[k]], #trait matrix
                   a = alist_red[[k]],
                   w.abun = T,     #weight metrics by relative abundance of species
                   stand.x = T,    #standardise trait to mean 0 and unit variance
                   corr = "cailliez",  #if distance matrix cannot be represented in euclidean space, take the sqrt of the distance
                   calc.FRic = T, stand.FRic = F, #do not standardise FRic to the global FRic (which includes all sp)
                   m = 2,          #the number of PCoA axes to use when calculating FRic
                   calc.FDiv = T, 
                   print.pco = T)  #return eigenvalues and PCoA axes
    
    FD_results[k , ]$ID <- IDlist[k]
    FD_results[k , ]$nsp <- result$nbsp
    FD_results[k , ]$FRic <- result$FRic
    FD_results[k , ]$qual.FRic <- result$qual.FRic
    FD_results[k , ]$FEve <- result$FEve
    FD_results[k , ]$FDiv <- result$FDiv
    FD_results[k , ]$FDis <- result$FDis
    FD_results[k , ]$RaoQ <- result$RaoQ
    
  }else {
    FD_results[k , ]$ID <- IDlist[k]
    FD_results[k , c(2:8)] <- "only one sp"
  }
}


write.csv(FD_results, "Functional trait data\\results\\FD_results_15Feb2024.csv")



###NInt ~ FRic models####
FD_results <- read.csv("Functional trait data\\results\\FD_results_15Feb2024.csv", row.names = 1) 

##Import NIntc results
nint_result <- 
  read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  filter(ID %in% c(FD_results$ID))

FD_results$FRic <- as.numeric(FD_results$FRic)
FD_results$qual.FRic <- as.numeric(FD_results$qual.FRic)
FD_results$FEve <- as.numeric(FD_results$FEve)
FD_results$FDiv <- as.numeric(FD_results$FDiv)
FD_results$RaoQ <- as.numeric(FD_results$RaoQ) 
FD_results$ID <- as.factor(FD_results$ID)

##Histograms of predictors and responses
ggplot(nint_result, aes(x = NIntc_richness)) +
  geom_histogram()

ggplot(FD_results, aes(x = FRic)) +
  geom_histogram()

ggplot(FD_results, aes(x = nsp, y = FRic)) +
  geom_point()

ggplot(FD_results, aes(x = FEve)) +
  geom_histogram()

ggplot(FD_results, aes(x = FDiv)) +
  geom_histogram()

ggplot(FD_results, aes(x = RaoQ)) +
  geom_histogram()


##Correlations
cormat <- cor(FD_results[which(!is.na(FD_results$FEve)) , which(colnames(FD_results) %in% c("FRic", "FEve", "FDiv"))],  method = "spearman")
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
  mutate(arid_sq = aridity^2) |> 
  filter(!is.na(FEve))

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
  mutate(arid_sq = aridity^2) |> 
  filter(!is.na(FEve))


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
  mutate(arid_sq = aridity^2) |> 
  filter(!is.na(FEve))

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
  mutate(arid_sq = aridity^2) |> 
  filter(!is.na(FEve))

is.numeric(nintc_cov_sum$FEve)
is.numeric(nintc_cov_sum$aridity)



###Create model formulas with different predictor combinations####
##Ceate a function called AllSubsets to do this

###CREATE substrRight FUNCTION TO USE IN AllSubsets####
substrRight <- function(x, n){substr(x, nchar(x) - n + 1, nchar(x))}

###CREATE AllSubsets FUNCTION####
AllSubsets <- function(ResponseVariableColumn, PredictorsColumns, data.source = data1, Add.PolynomialTerms = FALSE, Polynom.exclude = NA, 
                       Polynom.order = NA, Do.PredictorInteractions = FALSE, Interaction.Level = NA, ModelProportion = NA, 
                       Do.Random.effect = FALSE, random.effect = NA, scale.poly = TRUE, exclude_from_interactions = NA) {
  
  PredictorCombinations <- list()
  PredictorCombinations2 <- list()
  InteractionCombinations <- list()
  InteractionCombinations2 <- list()
  
  AllPredictorsColumns <- c(PredictorsColumns)
  
  # convert column numbers to variable names
  AllPredictorsNames <- colnames(data.source)[AllPredictorsColumns]
  ResponseVariable.name <- colnames(data.source)[ResponseVariableColumn]
  
  # count number of explanatory variables
  NumberExplanatoryVariables <- length(AllPredictorsNames)
  
  # add new response variable if using proportional odds model
  if (!is.na(ModelProportion == TRUE)) {ResponseVariable.name <- ModelProportion}
  
  # only allow polynomial terms or interaction terms
  if (Add.PolynomialTerms == TRUE) {Do.PredictorInteractions <- FALSE}
  if (Add.PolynomialTerms == TRUE & Polynom.order < 2 | Add.PolynomialTerms == TRUE & is.na(Polynom.order)) {Polynom.order <- 2}
  
  # add polynomial terms
  if (Add.PolynomialTerms == TRUE) {
    if (any(!is.na(Polynom.exclude))) {Polynom.predictors <- setdiff(PredictorsColumns, Polynom.exclude)} else {Polynom.predictors <- PredictorsColumns}
    Polynom.predictors.names <- colnames(data.source)[Polynom.predictors]
    Polynom.ColNames <- paste(Polynom.predictors.names, 2, sep = "")
    Polynom.ColNames2 <- Polynom.ColNames
    for (i in 1 : length(Polynom.predictors)) {
      if (scale.poly == TRUE) {
        temp.matrix <- scale(data.source[, Polynom.predictors[i]]^2)
      } else {
        temp.matrix <- data.source[, Polynom.predictors[i]]^2
      }
      data.source <- cbind(data.source, temp.matrix)
      colnames(data.source)[dim(data.source)[2]] <- Polynom.ColNames[i]
    }
    if (Polynom.order > 2) {
      for (j in 3 : Polynom.order) {
        Polynom.ColNames <- paste(Polynom.predictors.names, j, sep = "")
        Polynom.ColNames2 <- c(Polynom.ColNames2, Polynom.ColNames)
        for (i in 1 : length(Polynom.predictors)) {
          if (scale.poly == TRUE) {
            temp.matrix <- scale(data.source[, Polynom.predictors[i]]^j)
          } else {                    
            temp.matrix <- data.source[, Polynom.predictors[i]]^j
          }
          data.source <- cbind(data.source, temp.matrix)
          colnames(data.source)[dim(data.source)[2]] <- Polynom.ColNames[i]
        }
      }
    }
    AllPredictorsNames <- c(AllPredictorsNames, Polynom.ColNames2)
    NumberExplanatoryVariables <- length(AllPredictorsNames)
  }
  
  # add interaction terms if required
  if (Do.PredictorInteractions) {
    for (counter in 1 : Interaction.Level) {
      if(is.na(exclude_from_interactions) == TRUE) {
        temp.combn <- t(combn(AllPredictorsNames, counter))
      }else {
        temp.combn <- t(combn(AllPredictorsNames[-exclude_from_interactions], counter)) 
      }
      for (counter2 in 1 : length(temp.combn[,1])) {InteractionCombinations[[length(InteractionCombinations) + 1]] <- temp.combn[counter2, ]}}
    if(is.na(exclude_from_interactions) == FALSE) {
      InteractionCombinations <- InteractionCombinations[- c(1 : (NumberExplanatoryVariables - length(exclude_from_interactions)))]
    }
    # InteractionCombinations <- InteractionCombinations[- length(InteractionCombinations)]
    for (counter3 in 1 : length(InteractionCombinations)) {
      temp.number.predictors <- length(InteractionCombinations[[counter3]])
      for (counter4 in 2 : temp.number.predictors) {
        if (counter4 == 2) {temp.form <- paste(InteractionCombinations[[counter3]][1], InteractionCombinations[[counter3]][2], sep = ":")} else {temp.form <- paste(temp.form, InteractionCombinations[[counter3]][counter4], sep = ":")}
        InteractionCombinations2[[counter3]] <- temp.form
      }
    }
    AllPredictorsNames <- c(AllPredictorsNames, InteractionCombinations2)
    NumberExplanatoryVariables <- length(AllPredictorsNames)
  }
  
  
  # loop through counter of 1 : NumberExplanatoryVariables, creating each combination of variables, then placing in list
  for (counter in c(1 : NumberExplanatoryVariables)) {
    
    temp.combn <- t(combn(AllPredictorsNames, counter))
    #remove unessecary combos
    #Interactions where functional div indices interact, and where the interaction has an NA in it
    temp.combn <- t(t(temp.combn[-which(temp.combn[ , c(1:ncol(temp.combn))] %in% c("FRic:FEve", "FEve:FRic", "FRic:FDiv", "FDiv:FRic", "FEve:FDiv", "FDiv:FEve")) , ]))
    temp.combn <- t(t(temp.combn[-which(temp.combn[ , c(1:ncol(temp.combn))] %like% "%:NA%") , ]))
    
    for (counter2 in 1 : length(temp.combn[,1])) {
      PredictorCombinations[[length(PredictorCombinations) + 1]] <- temp.combn[counter2, ]
      
      #keep track of how many combinations there are in the console
      print(paste("counter= ", counter, "out of", NumberExplanatoryVariables, "counter2= ", counter2, "out of", length(temp.combn[,1])))
    }
  }
  
  # eliminate any potential models that contain interaction or polynomial, but not base or linear terms
  models.to.drop <- NULL
  
  # check for interaction terms
  for (check1 in 1 : length(PredictorCombinations)) {
    if(length(grep(":", PredictorCombinations[[check1]])) > 0) {
      interaction.terms.temp <- NULL
      for (check2 in grep(":", PredictorCombinations[[check1]])) {
        interaction.terms.temp <- c(interaction.terms.temp, strsplit(as.character(PredictorCombinations[[check1]][check2]), ":"))}
      interaction.terms.temp <- unique(unlist(interaction.terms.temp))
      Univariate.models <- setdiff(c(1 : length(PredictorCombinations[[check1]])), grep(":", PredictorCombinations[[check1]]))
      Univariate.terms <- NULL
      for (check3 in 1 : length(Univariate.models)) {
        Univariate.terms <- c(Univariate.terms, PredictorCombinations[[check1]][check3])}
      Univariate.terms <- unlist(Univariate.terms)
      if (any(is.na(match(interaction.terms.temp, Univariate.terms)))) {
        models.to.drop <- c(models.to.drop, check1)
      }
    }
  }
  PredictorCombinations[models.to.drop] <- NULL
  
  # check for polynomial terms
  if (is.na(Polynom.order)) { 
    Polynom.order = 2
    models.to.drop <- NULL
    for (check1 in 1 : length(PredictorCombinations)) {                 # loop through all models
      pred.terms <- unlist(PredictorCombinations[check1])
      for (check2 in seq(from = Polynom.order, to = 2)) {                             # loop through all possible polynomial terms
        for (check3 in 1 : length(pred.terms)) {                    # loop through each term
          if(substrRight(pred.terms[check3], 1) == check2) {
            base.var = substr(pred.terms[check3], 1, nchar(pred.terms[check3]) - 1)
            required.vars = c(base.var, paste(base.var, 2 : check2, sep = ""))
            # check if all required variables are present in model
            if (!all(required.vars %in% pred.terms) == TRUE) models.to.drop <- c(models.to.drop, check1)
          }}}}}
  
  
  # check that only unique model numbers given for models to drop
  models.to.drop <- unique(models.to.drop)
  PredictorCombinations[models.to.drop] <- NULL
  
  #loop through each combination in "PredictorCombinations"
  for (counter3 in 1 : length(PredictorCombinations)) {
    temp.number.predictors <- length(PredictorCombinations[[counter3]])
    for (counter4 in 1 : temp.number.predictors) {
      if (counter4 == 1) {temp.form <- paste(ResponseVariable.name, PredictorCombinations[[counter3]][1], sep = "~")} else {temp.form <- paste(temp.form, PredictorCombinations[[counter3]][counter4], sep = "+")}
      if (Do.Random.effect == TRUE) {PredictorCombinations2[[counter3]] <- paste(temp.form, random.effect, sep = "+")} else {PredictorCombinations2[[counter3]] <- temp.form}
    }
  }
  data1 <<- data.source	
  PredictorCombinations2
}  # end of function


####RUN AllSubsets####
formulas <- AllSubsets(ResponseVariableColumn = which(colnames(nintc_rich_sum) == "mean_NIntc_rich_binom"), 
                     PredictorsColumns = c(which(colnames(nintc_rich_sum) %in% c("FRic", "FEve", "FDiv", "graz", "aridity"))), 
                     data.source = nintc_rich_sum, 
                     Add.PolynomialTerms = FALSE, 
                     Do.PredictorInteractions = TRUE, 
                     Interaction.Level = 2, #interaction level = 3 takes wayyy too long
                     Do.Random.effect = TRUE, 
                     random.effect = "(1|site_ID)") 
#Save the formulas output because it takes so long to run

formula_table <- data.frame(formula  = character(length = length(formulas)))
for(f in 1:length(formulas)) {
  one_formula <- formulas[[f]]
  formula_table[f, 1] <- one_formula
}
write.csv(formula_table, "Functional trait data\\Results\\model_formulas.csv")


##Get formulas ready to run models with
#lets isolate the response combinations
formula_table <- read.csv("Functional trait data\\Results\\model_formulas.csv", row.names = 1) |> 
  separate_wider_delim(formula, delim = "~", names = c("response", "predictors")) |> 
  select(predictors) |> 
  distinct(predictors) |> 
  filter(!predictors %like% "%FRic:FEve%", #remove these interactions because we dont expect them to be at play
         !predictors %like% "%FRic:FDiv%", 
         !predictors %like% "%FEve:FDiv%", 
         !predictors %like% "%FEve:FRic%", 
         !predictors %like% "%FDiv:FRic%",
         !predictors %like% "%FDiv:FEve%") |> 
  add_row(predictors = "1+(1|site_ID)") #add the null model


###Loop through the formulas ####
#Create a table for results
results_table <- data.frame(Response = character(), Model = character(), Chisq = numeric(), 
                            Df = integer(), Pr_value = numeric(), AIC = numeric(), 
                            Warnings = character(), row.names = NULL)

# Initialize warning_msg outside the loop
warning_msg <- ""

##Also loop through response variables
response_list <- c("mean_NIntc_rich_binom", "mean_NIntc_cov_binom", "mean_NInta_rich_binom", "mean_NInta_cov_binom")
datalist = c("nintc_rich_sum", "nintc_cov_sum", "ninta_rich_sum", "ninta_cov_sum")

##LOOP TO CREATE MODELS STARTS HERE##
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
          anova_result <- Anova(model)
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

write.csv(results_table, "Functional trait data\\results\\model_results_21Feb2024.csv")


###Interpret model results####
mod_results <- read.csv("Functional trait data\\results\\model_results_21Feb2024.csv", row.names = 1)

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


#run the null models to double check for convergence problems
bmodel <- glmmTMB(mean_NIntc_rich_binom ~ (1)+ (1|site_ID), 
                  family = "binomial", data = nintc_rich_sum)
summary(bmodel)



