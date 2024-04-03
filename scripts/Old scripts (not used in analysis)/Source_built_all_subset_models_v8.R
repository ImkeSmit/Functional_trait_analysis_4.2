##Received from Pete: 20 Feb 2024

## all subsets model building code - 24 May 2023

library(lme4)
library(MuMIn)

# load data as CSV
data1 <- read.csv(file.choose())

# clean up data
data1$Regime <- as.factor(data1$Regime)
data1$Marsh_cat <- as.factor(data1$Marsh_cat)
data1$Patch_ID <- as.factor(data1$Patch_ID)
data1$Position_in_patch <- as.factor(data1$Position_in_patch)
data1$Ecotype <- as.factor(data1$Ecotype)
data1$Survival <- as.factor(data1$Survival)

data1 <- data1[!is.na(data1$Survival),]

data1$Survival_score <- as.numeric(data1$Survival)
data1$Survival_score[data1$Survival_score == 2] <- 0       # setting alive = 1, dead = 0

# set parameters
OutputFileName = "MidSeason_survival"
setwd("D:\\UP_research\\Isabelle_Buyens\\Model_building")
ResponseVariableDistribution <- "binomial"
ResponseVariableColumn <- 14          
number.random.effects = 1
random.effect = "(1|Patch_ID)"
final_env_cols <- c(2, 3, 5, 8, 10, 12)                # predictor columns
exclude_from_interactions <- c(2, 3, 4)

# generate list of all models to be tested
list.of.models <- AllSubsets(ResponseVariableColumn = ResponseVariableColumn, PredictorsColumns = final_env_cols, Do.Random.effect = TRUE, random.effect = random.effect, Do.PredictorInteractions = TRUE, Interaction.Level = 2) 



# temp, testing code
data.source = data1
PredictorsColumns <- final_env_cols
Interaction.Level = 2
Do.PredictorInteractions = TRUE


##########################
# run all models as LMMs #
##########################

library(lme4)
library(MuMIn)
library(lmerTest)

# run null model first, for comparison
model.counter <<- 0
ResponseVariable.name <- colnames(data1)[ResponseVariableColumn]
assign(paste(OutputFileName, "_Mod_0", sep = ""), lmer(eval(parse(text = paste(ResponseVariable.name, " ~ 1 + ", random.effect, sep = ""))), data = data1))
ExtractRegResLMM(eval(parse(text = paste(OutputFileName, "_Mod_0", sep = ""))), MaxVIF = 0, max.terms = 1)

for (model.counter in 1 : length(list.of.models)) {
    model.counter <<- model.counter
    assign(paste(OutputFileName, "_Mod_", model.counter, sep = ""), lmer(eval(parse(text = list.of.models[[model.counter]])), data = data1))
    temp.vif = 0
    ExtractRegResLMM(eval(parse(text = paste(OutputFileName, "_Mod_", model.counter, sep = ""))), MaxVIF = NA, max.terms = 1, null.mod.name = eval(parse(text = paste(OutputFileName, "_Mod_0", sep = ""))))
}




###########################
# run all models as GLMMs #      Still needs to be completed
###########################
library(lme4)
library(MuMIn)
library(lmerTest)

# run null model first, for comparison
model.counter <<- 0
ResponseVariable.name <- colnames(data1)[ResponseVariableColumn]
assign(paste(OutputFileName, "_Mod_0", sep = ""), glmer(eval(parse(text = paste(ResponseVariable.name, " ~ 1 + ", random.effect, sep = ""))), family = ResponseVariableDistribution, data = data1))
ExtractRegResGLMM(ModelOutput = eval(parse(text = paste(OutputFileName, "_Mod_0", sep = ""))), max.terms = 1)

for (model.counter in 1 : length(list.of.models)) {
    model.counter <<- model.counter
    assign(paste(OutputFileName, "_Mod_", model.counter, sep = ""), glmer(eval(parse(text = list.of.models[[model.counter]])), family = ResponseVariableDistribution, data = data1))
    ExtractRegResGLMM(ModelOutput = eval(parse(text = paste(OutputFileName, "_Mod_", model.counter, sep = ""))), max.terms = 1, null.mod.name = eval(parse(text = paste(OutputFileName, "_Mod_0", sep = ""))))
    print(model.counter)
}





###################
# Short functions #
###################

substrRight <- function(x, n){substr(x, nchar(x) - n + 1, nchar(x))}
D2 <- function(glm.output) {(glm.output$null.deviance - glm.output$deviance) / glm.output$null.deviance}

StandardizeCoeffs <- function(glm.output, predictor.number = 1, dataset = "data1") {
    raw.predictor <- eval(parse(text = paste(dataset, "$", row.names(summary(glm.output)$coeff)[predictor.number + 1], sep = "")))
    raw.response <-  eval(parse(text = paste(dataset, "$", glm.output$terms[[2]], sep = "")))
    # stand. coeff = raw coeff * (predictor SD / response SD)
    summary(glm.output)$coeff[predictor.number, 1] * (sd(raw.predictor) / sd(raw.response))
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

####################################################################
# Function for creating all possible models from set of predictors #
####################################################################

AllSubsets <- function(ResponseVariableColumn, PredictorsColumns, data.source = data1, Add.PolynomialTerms = FALSE, Polynom.exclude = NA, Polynom.order = NA, Do.PredictorInteractions = FALSE, Interaction.Level = NA, ModelProportion = NA, Do.Random.effect = FALSE, random.effect = NA, scale.poly = TRUE) {

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
        if (Add.PolynomialTerms == TRUE) Do.PredictorInteractions <- FALSE
        if (Add.PolynomialTerms == TRUE & Polynom.order < 2 | Add.PolynomialTerms == TRUE & is.na(Polynom.order)) Polynom.order <- 2
        
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
            		temp.combn <- t(combn(AllPredictorsNames[-exclude_from_interactions], counter))
            		for (counter2 in 1 : length(temp.combn[,1])) {InteractionCombinations[[length(InteractionCombinations) + 1]] <- temp.combn[counter2, ]}}
           InteractionCombinations <- InteractionCombinations[- c(1 : (NumberExplanatoryVariables - length(exclude_from_interactions)))]
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
       	for (counter in 1 : NumberExplanatoryVariables) {
        		temp.combn <- t(combn(AllPredictorsNames, counter))
        		for (counter2 in 1 : length(temp.combn[,1])) {PredictorCombinations[[length(PredictorCombinations) + 1]] <- temp.combn[counter2, ]}}
        		
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
        if (is.na(Polynom.order)) Polynom.order = 2
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
        }}}}
           

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



#################################################
# Create function to extract regression results #  FOR GLM
#################################################

ExtractRegRes <- function (ModelOutput, OutputFile = OutputFileName, MaxVIF = NA, max.terms) {
	ModelOutputSummary <- summary(ModelOutput)
	TempResult <- data.frame("Result" = NA)
		TempResult[1] <- ResponseVariable.name
		TempResult[2] <- ifelse(model.counter == 0, "null", noquote(as.character(eval(parse(text = ModelOutput$call[2]))[3])))
		TempResult[3] <- length(data1[, 1])
		TempResult[4] <- ModelOutput$df.residual
		TempResult[5] <- round(ModelOutput$null.deviance - ModelOutput$deviance, 2)	#model stat - a.k.a. likelihood ratio statistic
		TempResult[6] <- round(1 - pchisq(ModelOutput$null.deviance - ModelOutput$deviance, ModelOutput$df.null - ModelOutput$df.residual), 3)
		TempResult[7] <- ModelOutput$aic
		TempResult[8] <- ModelOutput$deviance				# G2 / Deviance = measure of goodness of fit
		TempResult[9] <- 1 - pchisq(ModelOutput$deviance, ModelOutput$df.residual)
		TempResult[10] <- length(ModelOutput$coeff)
		TempResult[11] <- ModelOutput$null.deviance
		TempResult[12] <- ModelOutput$deviance
		TempResult[13] <- round((ModelOutput$null.deviance - ModelOutput$deviance)/ModelOutput$null.deviance, 4)
		TempResult[14] <- round(1 - (((length(data1[, 1]) - 1)/(length(data1[, 1]) - length(ModelOutput$coeff))) * (1 - ((ModelOutput$null.deviance - ModelOutput$deviance) / ModelOutput$null.deviance))), 4)
		TempResult[15] <- MaxVIF
		for (ParameterNumber in 1 : length(ModelOutput$coeff)) {
			if (dim(ModelOutputSummary$coeff)[1] < ParameterNumber) {break}
			TempResult[16 + ((ParameterNumber - 1) * 5)] <- noquote(names(ModelOutput$coeff)[ParameterNumber])
			TempResult[17 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber,1], 4)
			TempResult[18 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber,2], 4)
			TempResult[19 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber,3], 4)
			TempResult[20 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber,4], 3)}
    # add extra headers
    if (model.counter == 0) {TempResult <- cbind(TempResult, t(data.frame("empty" = rep("", max.terms * 5))))}

	if (file.exists(paste(OutputFile, ".txt", sep = ""))) {
		write.table(TempResult, file = paste(OutputFile, ".txt", sep = ""), append = ifelse(model.counter == 0, FALSE, TRUE), sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = FALSE, row.names = FALSE)
		if (model.counter == 0) {
        print("Old output file overwritten")
        }
		} else {
		Results.lab1 <- c("ResponseVariable", "Model", "n", "Residual d.f.", "Model stat.", "Model p", "AIC", "G2 g.o.f.", "g2 p", "Number parameters", "Null dev.", "Residual dev.", "%dev explained", "adj. D2")
		Results.lab3 <- "Max_VIF"
		Results.lab2 <- c("Variable", "Coefficient", "SE", "Stat.", "p-value")
    counter2 <- length(ModelOutput$coeff) + 1
			for (ParameterNumber in 1:counter2) {
			for (counter in 1 : 5) {
			Results.lab3 <- c(Results.lab3, paste(Results.lab2[counter], "_", ParameterNumber, sep = ""))}}
		names(TempResult) <- c(Results.lab1, Results.lab3)
		write.table(TempResult, file = paste(OutputFile, ".txt", sep = ""), append = FALSE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = TRUE, row.names = FALSE)}

} # ends function definition


#################################################
# Create function to extract regression results #  FOR LMM
#################################################

ExtractRegResLMM <- function (ModelOutput, OutputFile = OutputFileName, MaxVIF = NA, max.terms, null.mod.name) {
	ModelOutputSummary <- summary(ModelOutput)
	TempResult <- data.frame("Result" = NA)
		TempResult[1] <- ResponseVariable.name
		TempResult[2] <- ifelse(model.counter == 0, "null", paste(attr(ModelOutput@frame, "formula")[3]))
		TempResult[3] <- length(data1[, 1])
		TempResult[4] <- NA
		TempResult[5] <- NA
		TempResult[6] <- ifelse(model.counter == 0, "null", anova(null.mod.name, ModelOutput)$Pr[2])
		TempResult[7] <- AIC(ModelOutput)
		TempResult[8] <- NA
		TempResult[9] <- NA
		TempResult[10] <- length(coef(ModelOutputSummary)[,1]) - 1
		TempResult[11] <- NA
		TempResult[12] <- NA
		TempResult[13] <- r.squaredGLMM(ModelOutput)[2]
		TempResult[14] <- r.squaredGLMM(ModelOutput)[1]
		TempResult[15] <- NA
		for (ParameterNumber in 1 : dim(coef(ModelOutputSummary))[1]) {
			if (dim(coef(ModelOutputSummary))[1] < ParameterNumber) {break}
			TempResult[16 + ((ParameterNumber - 1) * 5)] <- noquote(row.names(coef(ModelOutputSummary))[ParameterNumber])
			TempResult[17 + ((ParameterNumber - 1) * 5)] <- round(coef(ModelOutputSummary)[ParameterNumber, 1], 4)
			TempResult[18 + ((ParameterNumber - 1) * 5)] <- round(coef(ModelOutputSummary)[ParameterNumber, 2], 4)
			TempResult[19 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber, 4], 4)
			TempResult[20 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber, 5], 4)}       # If needed use lmerTest: https://stackoverflow.com/questions/37336510/how-to-extract-fixed-effects-part-of-summary-from-lme4
    # add extra headers
    if (model.counter == 0) {TempResult <- cbind(TempResult, t(data.frame("empty" = rep("", max.terms * 5))))}

	if (file.exists(paste(OutputFile, ".txt", sep = ""))) {
		write.table(TempResult, file = paste(OutputFile, ".txt", sep = ""), append = ifelse(model.counter == 0, FALSE, TRUE), sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = FALSE, row.names = FALSE)
		if (model.counter == 0) {
        print("Old output file overwritten")
        }
		} else {
		Results.lab1 <- c("ResponseVariable", "Model", "n", "Residual d.f.", "Model stat.", "Model p", "AIC", "G2 g.o.f.", "g2 p", "Number parameters", "Null dev.", "Residual dev.", "r2_cond", "r2_marg")
		Results.lab3 <- "Max_VIF"
		Results.lab2 <- c("Variable", "Coefficient", "SE", "Stat.", "p-value")
    counter2 <- dim(coef(ModelOutputSummary))[1]
			for (ParameterNumber in 1:counter2) {
			for (counter in 1 : 5) {
			Results.lab3 <- c(Results.lab3, paste(Results.lab2[counter], "_", ParameterNumber, sep = ""))}}
		names(TempResult) <- c(Results.lab1, Results.lab3)
		write.table(TempResult, file = paste(OutputFile, ".txt", sep = ""), append = FALSE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = TRUE, row.names = FALSE)}

} # ends function definition





#################################################
# Create function to extract regression results #  FOR GLMM
#################################################

ExtractRegResGLMM <- function (ModelOutput, OutputFile = OutputFileName, max.terms, null.mod.name) {
	ModelOutputSummary <- summary(ModelOutput)
	TempResult <- data.frame("Result" = NA)
		TempResult[1] <- ResponseVariable.name
		TempResult[2] <- ifelse(model.counter == 0, "null", paste(attr(ModelOutput@frame, "formula")[3]))
		TempResult[3] <- length(data1[, 1])
		TempResult[4] <- NA
		TempResult[5] <- NA
		TempResult[6] <- ifelse(model.counter == 0, "null", anova(null.mod.name, ModelOutput)$Pr[2])
		TempResult[7] <- AIC(ModelOutput)
		TempResult[8] <- NA
		TempResult[9] <- NA
		TempResult[10] <- length(coef(ModelOutputSummary)[,1]) - 1
		TempResult[11] <- NA
		TempResult[12] <- NA
		TempResult[13] <- r.squaredGLMM(ModelOutput)[2]
		TempResult[14] <- r.squaredGLMM(ModelOutput)[1]
		TempResult[15] <- NA
		for (ParameterNumber in 1 : dim(coef(ModelOutputSummary))[1]) {
			if (dim(coef(ModelOutputSummary))[1] < ParameterNumber) {break}
			TempResult[16 + ((ParameterNumber - 1) * 5)] <- noquote(row.names(coef(ModelOutputSummary))[ParameterNumber])
			TempResult[17 + ((ParameterNumber - 1) * 5)] <- round(coef(ModelOutputSummary)[ParameterNumber, 1], 4)
			TempResult[18 + ((ParameterNumber - 1) * 5)] <- round(coef(ModelOutputSummary)[ParameterNumber, 2], 4)
			TempResult[19 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber, 3], 4)
			TempResult[20 + ((ParameterNumber - 1) * 5)] <- round(ModelOutputSummary$coeff[ParameterNumber, 4], 4)}       # If needed use lmerTest: https://stackoverflow.com/questions/37336510/how-to-extract-fixed-effects-part-of-summary-from-lme4
    # add extra headers
    if (model.counter == 0) {TempResult <- cbind(TempResult, t(data.frame("empty" = rep("", max.terms * 5))))}

	if (file.exists(paste(OutputFile, ".txt", sep = ""))) {
		write.table(TempResult, file = paste(OutputFile, ".txt", sep = ""), append = ifelse(model.counter == 0, FALSE, TRUE), sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = FALSE, row.names = FALSE)
		if (model.counter == 0) {
        print("Old output file overwritten")
        }
		} else {
		Results.lab1 <- c("ResponseVariable", "Model", "n", "Residual d.f.", "Model stat.", "Model p", "AIC", "G2 g.o.f.", "g2 p", "Number parameters", "Null dev.", "Residual dev.", "r2_cond", "r2_marg")
		Results.lab3 <- "Max_VIF"
		Results.lab2 <- c("Variable", "Coefficient", "SE", "Stat.", "p-value")
    counter2 <- dim(coef(ModelOutputSummary))[1]
			for (ParameterNumber in 1:counter2) {
			for (counter in 1 : 5) {
			Results.lab3 <- c(Results.lab3, paste(Results.lab2[counter], "_", ParameterNumber, sep = ""))}}
		names(TempResult) <- c(Results.lab1, Results.lab3)
		write.table(TempResult, file = paste(OutputFile, ".txt", sep = ""), append = FALSE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = TRUE, row.names = FALSE)}

} # ends function definition





#############################################################
# Function which combines the above functions into one call # 
#############################################################


RunSubsets <- function (ResponseVariableColumn, PredictorsColumns, data1 = data1) {                    
    if (!is.numeric(ResponseVariableColumn)) ResponseVariableColumn <- which(colnames(data1) == ResponseVariableColumn)
    if (any(!is.numeric(PredictorsColumns))) {
        for (i in 1 : length(PredictorsColumns)) {
            if (!is.numeric(PredictorsColumns[i])) {PredictorsColumns[i] <- which(colnames(data1) == PredictorsColumns[i])}}}
    PredictorsColumns <- sort(PredictorsColumns)
    
    ResponseVariable.name <<- colnames(data1)[ResponseVariableColumn]
    Predictors.names <- colnames(data1)[PredictorsColumns]
    
    # create all subsets of predictor variables
    list.of.models <<- AllSubsets(ResponseVariableColumn = ResponseVariableColumn, PredictorsColumns = PredictorsColumns, data.source = data1, Add.PolynomialTerms = c.Add.PolynomialTerms,
                        Polynom.exclude = c.Polynom.exclude, Polynom.order = c.Polynom.order, Do.PredictorInteractions = c.Do.PredictorInteractions, Interaction.Level = c.Interaction.Level)
    
    # run null model first, for comparison
    model.counter <<- 0
    assign(paste(OutputFileName, "_Mod_0", sep = ""), glm(eval(parse(text = ResponseVariable.name)) ~ 1, family = ResponseVariableDistribution, data = data1))
    ExtractRegRes(eval(parse(text = paste(OutputFileName, "_Mod_0", sep = ""))), MaxVIF = 0, max.terms = length(unlist(list.of.models[length(list.of.models)])))
    
    # run all models as GLMs
    for (model.counter in 1 : length(list.of.models)) {
        model.counter <<- model.counter
        assign(paste(OutputFileName, "_Mod_", model.counter, sep = ""), glm(eval(parse(text = list.of.models[[model.counter]])), family = ResponseVariableDistribution, data = data1))
        if (eval(parse(text = paste("length(attr(", OutputFileName, "_Mod_", model.counter ,"$terms, 'term.labels'))", sep = ""))) > 1) {temp.vif <- round(max(vif(eval(parse(text = paste(OutputFileName, "_Mod_", model.counter, sep = ""))))), 3) } else {temp.vif = 0}
        ExtractRegRes(eval(parse(text = paste(OutputFileName, "_Mod_", model.counter, sep = ""))), MaxVIF = temp.vif, max.terms = length(unlist(list.of.models[length(list.of.models)])))
    }
    
    # open the output file, sort by AIC and extract best model
    max.columns <- dim(read.delim(paste(OutputFileName, ".txt", sep = ""), header = FALSE, skip = model.counter + 1))[2]
    GLM.results <- read.delim(paste(OutputFileName, ".txt", sep = ""), skip = 1, header = FALSE, col.names = paste("V", 1 : max.columns, sep = ""))
    BestModel.row <- which(GLM.results[, 7] == min(GLM.results[, 7]))
    GLM.results <- GLM.results[sort.list(GLM.results[, 7]), ]
    write.table(GLM.results, file = paste(OutputFileName, "_sorted.txt", sep = ""), append = FALSE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = FALSE, row.names = FALSE)
    assign(paste("Best_model_", OutputFileName, sep = ""), eval(parse(text = paste(OutputFileName, "_Mod_", BestModel.row - 1, sep = ""))))		
    # remove(list = ls(pattern = paste("^", OutputFileName, sep = "")))

    eval(parse(text = paste("Best_model_", OutputFileName, sep = "")))    
}
