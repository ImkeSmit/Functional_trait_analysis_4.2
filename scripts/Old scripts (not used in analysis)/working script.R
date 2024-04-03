substrRight <- function(x, n){substr(x, nchar(x) - n + 1, nchar(x))}


ResponseVariableColumn = which(colnames(nintc_rich_sum) == "mean_NIntc_rich_binom")
PredictorsColumns = c(which(colnames(nintc_rich_sum) %in% c("FRic" ,"graz", "aridity")))
data.source = nintc_rich_sum
Add.PolynomialTerms = TRUE
Polynom.order = 2
Polynom.exclude = c(which(colnames(nintc_rich_sum) %in% c("FRic" ,"graz")))
Do.PredictorInteractions = TRUE 
Interaction.Level = 2
Do.Random.effect = TRUE
random.effect = "(1|site_ID)"
exclude_from_interactions = NA


ModelProportion = NA
scale.poly = FALSE



####CREATE FUNCTION####
AllSubsets <- function(ResponseVariableColumn, PredictorsColumns, data.source = data1, Add.PolynomialTerms = FALSE, Polynom.exclude = NA, 
                       Polynom.order = NA, Do.PredictorInteractions = FALSE, Interaction.Level = NA, ModelProportion = NA, 
                       Do.Random.effect = FALSE, random.effect = NA, scale.poly = TRUE, exclude_from_interactions = NA) {
  
  PredictorCombinations <- list()
  PredictorCombinations2 <- list()
  InteractionCombinations <- list()
  InteractionCombinations2 <- list()
  InteractionCombinations3 <- list()
  
  AllPredictorsColumns <- c(PredictorsColumns)
  
  # convert column numbers to variable names
  AllPredictorsNames <- colnames(data.source)[AllPredictorsColumns]
  ResponseVariable.name <- colnames(data.source)[ResponseVariableColumn]
  
  # count number of explanatory variables
  NumberExplanatoryVariables <- length(AllPredictorsNames)
  
  # add new response variable if using proportional odds model
  if (!is.na(ModelProportion == TRUE)) {ResponseVariable.name <- ModelProportion}
  
  # only allow polynomial terms or interaction terms
  #if (Add.PolynomialTerms == TRUE) {Do.PredictorInteractions <- FALSE}
  #if (Add.PolynomialTerms == TRUE & Polynom.order < 2 | Add.PolynomialTerms == TRUE & is.na(Polynom.order)) {Polynom.order <- 2}
  
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
      for (counter2 in 1 : length(temp.combn[,1])) {
        InteractionCombinations[[length(InteractionCombinations) + 1]] <- temp.combn[counter2, ]
        }}
    if(is.na(exclude_from_interactions) == FALSE) {
      InteractionCombinations <- InteractionCombinations[- c(1 : (NumberExplanatoryVariables - length(exclude_from_interactions)))]
        }
    # InteractionCombinations <- InteractionCombinations[- length(InteractionCombinations)]
    
    #remove intercations with NA in them
    for(i in 1:length(InteractionCombinations)) {
      if(length(InteractionCombinations[[i]]) > 1) {
        InteractionCombinations2[[length(InteractionCombinations2) +1]] <- InteractionCombinations[[i]]
      }
    } 
    
    for (counter3 in 1 : length(InteractionCombinations2)) {
      temp.number.predictors <- length(InteractionCombinations2[[counter3]])
      for (counter4 in 2 : temp.number.predictors) {
        if (counter4 == 2) {
          temp.form <- paste(InteractionCombinations2[[counter3]][1], InteractionCombinations2[[counter3]][2], sep = ":")
          } else {temp.form <- paste(temp.form, InteractionCombinations2[[counter3]][counter4], sep = ":")}
        InteractionCombinations3[[counter3]] <- temp.form
      }
    }
    #remove the aridity:aridity2 interaction
    InteractionCombinations3 <- InteractionCombinations3[- which(InteractionCombinations3 == "aridity:aridity2")]
    AllPredictorsNames <- c(AllPredictorsNames, InteractionCombinations3)
    NumberExplanatoryVariables <- length(AllPredictorsNames)
  }
  
  
  # loop through counter of 1 : NumberExplanatoryVariables, creating each combination of variables, then placing in list
  for (counter in c(1 : NumberExplanatoryVariables)) {
    
    temp.combn <- t(combn(AllPredictorsNames, counter))
    #remove unessecary combos
    #Interactions where functional div indices interact, and where the interaction has an NA in it
    #temp.combn <- t(t(temp.combn[-which(temp.combn[ , c(1:ncol(temp.combn))] %in% c("FRic:FEve", "FEve:FRic", "FRic:FDiv", "FDiv:FRic", "FEve:FDiv", "FDiv:FEve")) , ]))
    #temp.combn <- t(t(temp.combn[-which(temp.combn[ , c(1:ncol(temp.combn))] %like% "%:NA%") , ]))

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
    if(length(grep(":", PredictorCombinations[[check1]])) > 0) { #if there are interactions in the model
      interaction.terms.temp <- NULL
      for (check2 in grep(":", PredictorCombinations[[check1]])) {
        interaction.terms.temp <- c(interaction.terms.temp, strsplit(as.character(PredictorCombinations[[check1]][check2]), ":"))
        }
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
              if (!all(required.vars %in% pred.terms) == TRUE) {
                models.to.drop <- c(models.to.drop, check1)
                }
            }}}}
  
  
  # check that only unique model numbers given for models to drop
  models.to.drop <- unique(models.to.drop)
  PredictorCombinations[models.to.drop] <- NULL
  }
  
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



####RUN FUNCTION####
result <- AllSubsets(ResponseVariableColumn = which(colnames(nintc_rich_sum) == "mean_NIntc_rich_binom"),
                     PredictorsColumns = c(which(colnames(nintc_rich_sum) %in% c("FRic" ,"graz", "aridity"))),
                     data.source = nintc_rich_sum,
                     Add.PolynomialTerms = TRUE,
                     Polynom.order = 2,
                     Polynom.exclude = c(which(colnames(nintc_rich_sum) %in% c("FRic" ,"graz"))),
                     Do.PredictorInteractions = TRUE, 
                     Interaction.Level = 2,
                     Do.Random.effect = TRUE,
                     random.effect = "(1|site_ID)",
                     exclude_from_interactions = NA,
                     ModelProportion = NA,
                     scale.poly = FALSE) 
