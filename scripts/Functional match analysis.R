###FUNCTIONAL MATCH ANALYSIS###
library(tidyverse)
library(tidylog)
library(vegan)
library(ggplot2)
library(ggbiplot)
library(traitstrap)

#From the filled trait data for plotspecific species
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv", row.names = 1)

##We need to standardise each trait value to mean 0 and variance of 1 with Z = (x-mean)/sd
##We need only one trait value per species

#Get the mean of each trait for each sp
FT_mean <- FT |> 
  group_by(taxon,trait) |> 
  summarise(mean_value = mean(value))

#Get it into wide format
FT_wide <- FT_mean |>
  pivot_wider(names_from = trait, values_from = mean_value) |> 
  column_to_rownames(var = "taxon") |> 
  filter(!is.na(MaxH),
         !is.na(MaxLS), 
         !is.na(MeanLA),
         !is.na(MeanLDMC),
         !is.na(MeanLL),
         !is.na(MeanSLA), 
         !is.na(percentC),
         !is.na(percentN)) |> 
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN))



#standardise trait values
std_FT_wide <- FT_wide

traitlist <- c("MaxH","MaxLS","MeanLA","MeanLDMC","MeanLL","MeanSLA","C_N_ratio")

for(t in 1:length(traitlist)) {
  
  #get the grand mean and sd for each trait
  m <- FT_wide |> 
    summarise(m = mean(FT_wide[ , which(colnames(FT_wide) == traitlist[t])]))
  
  sd <- FT_wide |> 
    summarise(m = sd(FT_wide[ , which(colnames(FT_wide) == traitlist[t])]))
  
  for (i in 1:nrow(FT_wide)) {
    
    raw_value <- FT_wide[i , which(colnames(FT_wide) == traitlist[t])]
    #replace the raw value with the standardised value using m and sd
    std_FT_wide[i , which(colnames(std_FT_wide) == traitlist[t])] <- (raw_value - m)/sd
  
    }
}



###CLUSTERING####
clusters <- hclust(dist(std_FT_wide, method = "euclidean"), method = "ward.D")
summary(clusters)
plot(clusters)
abline(h = 60, col = "red")

FGR <- as.data.frame(cutree(clusters, h = 60))
FGR[2] <- rownames(FGR)
rownames(FGR) <- c(1:nrow(FGR))
colnames(FGR) <- c("Functional_group", "taxon")

write.csv(FGR, "Functional trait data\\results\\Functional_groups3.csv")
FGR <- read.csv("Functional trait data\\results\\Functional_groups3.csv", row.names = 1)

###Now we need to see in which FGR nurses and target plants fall on the aridity gradient####
#We require raw country data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}


#Get the species that are nurses, growing with nurses, or in the bare microsite in each country
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  if(t == 1) {
  cou_ass <- cou |> 
    filter(!is.na(Species.within.quadrat)) |> 
    group_by(ID, ID_Microsite) |>
    summarise(target_sp = unique(Species.within.quadrat))
  } 
  else{
    temp_cou_ass <- cou |> 
      filter(!is.na(Species.within.quadrat)) |> 
      group_by(ID, ID_Microsite) |>
      summarise(target_sp = unique(Species.within.quadrat))
    
    cou_ass <- bind_rows(cou_ass, temp_cou_ass)
  }
}


##Now we have to merge the FGR on to cou_ass
FGR_cou_ass <- cou_ass |> 
  left_join(FGR, by = c("ID_Microsite" = "taxon")) |> 
  rename(nurse = ID_Microsite, nurse_functional_group = Functional_group) |> 
  left_join(FGR, by = c("target_sp" = "taxon")) |> 
  rename(target_functional_group = Functional_group) |> 
  mutate(nurse_target_match = str_c(nurse_functional_group, target_functional_group, sep = "_")) |> 
  ungroup()
#The functional group columns have NA's where that sp is not in the FGR data. That means we don't have all the trait data for that sp


##Merge the plotinfo to FGR_cou_ass
#Import the information about the biodesert sites
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  filter(!is.na(ID))

FGR_cou_ass <- FGR_cou_ass |> 
  left_join(siteinfo, by = "ID")

#save the file
write.csv(FGR_cou_ass, "Functional trait data//results//FGR3_nurse_target_species.csv")


###What traits do each FGR have - PCA####
FGR_cou_ass <- read.csv("Functional trait data//results//FGR3_nurse_target_species.csv", row.names = 1)

FT_wide_FGR <- std_FT_wide |> 
  mutate(taxon = rownames(std_FT_wide)) |> 
  left_join(FGR, by = "taxon")
FT_wide_FGR$Functional_group <- as.factor(FT_wide_FGR$Functional_group) 
#write to csv to use in graphing
write.csv(FT_wide_FGR, "Functional trait data//results//standard_FT_FGR.csv")


base_pca <- princomp(FT_wide_FGR[ , -c(8,9)], scores = T)
summary(base_pca) #proportion of variance is the variance explained by the PC
base_pca$scores #
base_pca$loadings #How much each var contributed to building the PC
base_pca$scale #scaling applied to each variable
base_pca$center #means
biplot(base_pca, choices = c("Comp.1", "Comp.2"))
biplot(base_pca, choices = c("Comp.1", "Comp.3"))
biplot(base_pca$scores)

pc1_2 <- ggbiplot(base_pca, varname.color = "black", 
         groups = FT_wide_FGR$Functional_group, 
         ellipse = FALSE) +
  theme_classic()
ggsave("PCA_1_2_3FGR.png", pc1_2, path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


pc1_3 <- ggbiplot(base_pca, choices = c(1,3),
         groups = FT_wide_FGR$Functional_group, 
         ellipse = FALSE)+
  theme_classic()
ggsave("PCA_1_3_3FGR.png", pc1_3, path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")



###How does the functional match change along aridity?
##Get unique nurse and target FGR's tied to ID and aridity
Fmatch <- FGR_cou_ass |> 
  select(ID, nurse_functional_group, target_functional_group, nurse_target_match, ARIDITY.v3, GRAZ) |> 
  filter(!is.na(nurse_target_match)) |> 
  distinct() |> 
  arrange(ARIDITY.v3)

#Get counts of each FGR per ID
Fmatch_counts <- Fmatch |> 
  group_by(ID, ARIDITY.v3, nurse_functional_group) |> 
  mutate(nurse_FGR_count= n()) |> 
  ungroup() |> 
  group_by(ID, ARIDITY.v3, target_functional_group) |> 
  mutate(target_FGR_count= n()) |> 
  ungroup() |> 
  group_by(ID, ARIDITY.v3, nurse_target_match) |> 
  mutate(match_FGR_count= n()) |> 
  ungroup() |> 
  distinct()
  
  
Fmatch_counts$ID <- as.factor(Fmatch_counts$ID)
Fmatch_counts$nurse_functional_group <- as.factor(Fmatch_counts$nurse_functional_group)
Fmatch_counts$target_functional_group <- as.factor(Fmatch_counts$target_functional_group)

##Let's plot:
#order by aridity
Fmatch_counts$ID <- with(Fmatch_counts, reorder(ID, ARIDITY.v3, mean)) #order ID by aridity

#counts of nurse FGR
Fmatch_counts |> 
  ggplot(aes(x = ID, y = nurse_FGR_count, fill = nurse_functional_group)) +
    geom_bar(stat = "identity")

#targetFGR
Fmatch_counts |> 
  ggplot(aes(x = ID, y = target_FGR_count, fill = target_functional_group)) +
  geom_bar(stat = "identity")

#nurse-targetmatch
Fmatch_counts |> 
  ggplot(aes(x = ID, y = match_FGR_count, fill = nurse_target_match)) +
  geom_bar(stat = "identity")


#How does the functional match change with grazing pressure?
Fmatch_counts_bygraz <- Fmatch |> 
  group_by(GRAZ, nurse_functional_group) |> 
  mutate(nurse_FGR_count= n()) |> 
  ungroup() |> 
  group_by(GRAZ, target_functional_group) |> 
  mutate(target_FGR_count= n()) |> 
  ungroup() |> 
  group_by(GRAZ, nurse_target_match) |> 
  mutate(match_FGR_count= n()) |> 
  ungroup() |> 
  select(!c(ID, ARIDITY.v3)) |> 
  distinct() |> 
  arrange(GRAZ)

Fmatch_counts_bygraz$GRAZ <- factor(Fmatch_counts_bygraz$GRAZ, levels = c("0", "1", "2", "3"))
Fmatch_counts_bygraz$nurse_functional_group <- as.factor(Fmatch_counts_bygraz$nurse_functional_group)
Fmatch_counts_bygraz$target_functional_group <- as.factor(Fmatch_counts_bygraz$target_functional_group)

#nurse FGR  
Fmatch_counts_bygraz |> 
  ggplot(aes(x = GRAZ, y = nurse_FGR_count, fill = nurse_functional_group)) +
  geom_bar(stat = "identity")

#target FGR
Fmatch_counts_bygraz |> 
  ggplot(aes(x = GRAZ, y = target_FGR_count, fill = target_functional_group)) +
  geom_bar(stat = "identity")

#functional match
Fmatch_counts_bygraz |> 
  ggplot(aes(x = GRAZ, y = match_FGR_count, fill = nurse_target_match)) +
  geom_bar(stat = "identity")


###Chisq tests to look for association between nurse-target match and interaction direction####
#rbind raw country data, use only variables we need
vars <- c("ID","Number.of.replicate", "Microsite", "ID_Microsite", "Species.within.quadrat" )
raw_fac_survey <- rbind(algeria[ , which(colnames(algeria) %in% vars)], argentina[ , which(colnames(argentina) %in% vars)], australia[ , which(colnames(australia) %in% vars)], chile[ , which(colnames(chile) %in% vars)], 
              chinachong[ , which(colnames(chinachong) %in% vars)], chinaxin[ , which(colnames(chinaxin) %in% vars)], iranabedi[ , which(colnames(iranabedi) %in% vars)], iranfarzam[ , which(colnames(iranfarzam) %in% vars)], 
              israel[ , which(colnames(israel) %in% vars)], namibiablaum[ , which(colnames(namibiablaum) %in% vars)], namibiawang[ , which(colnames(namibiawang) %in% vars)], southafrica[ , which(colnames(southafrica) %in% vars)],  
              spainmaestre[ , which(colnames(spainmaestre) %in% vars)], spainrey[ , which(colnames(spainrey) %in% vars)])

#import nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1)   
#classify intercations based on nintc richness
nint_result_direction <- nint_result |>
  mutate(interaction = case_when(NIntc_richness > 0 ~  "facilitation", 
                                 NIntc_richness < 0 ~  "competition", 
                                 NIntc_richness == 0 ~  "neutral")) |> 
  mutate(ID_rep = str_c(ID, replicate_no, sep = "_")) #create an identifier from ID and rep



#Get the nurse and target sp associated with these fac rep
sum_fac_survey <- raw_fac_survey |> 
  #remove bare microsites
  filter(Microsite == 2) |> 
  mutate(ID_rep = str_c(ID, Number.of.replicate, sep = "_")) |>  #create an identifier from ID and rep
  #add FGR of nurse sp
  left_join(FGR, by = c("ID_Microsite" = "taxon")) |> #inner join will remove sp that do not have functional group (ie for which we dont have complete FT data)
  rename(nurse_functional_group = Functional_group) |> 
  #add FGR of target sp |>
  left_join(FGR, by = c("Species.within.quadrat" = "taxon")) |> #inner join will remove sp that do not have functional group (ie for which we dont have complete FT data)
  rename(target_functional_group = Functional_group) |> 
  mutate(nurse_target_match = str_c(nurse_functional_group, target_functional_group, sep = "_")) |> 
  #join the direction of the interactions
  left_join(nint_result_direction[, which(colnames(nint_result_direction) %in% c("ID_rep", "interaction"))], by = "ID_rep")

#Do a chi test over all plots
complete_con_table <- sum_fac_survey |> 
  select(nurse_target_match, interaction) |> 
  filter(!is.na(nurse_target_match), 
         !nurse_target_match == "3_1") #remove this row because it has expected vals lower than 5
complete_con_table <- table(complete_con_table)

complete_chitest <- chisq.test(complete_con_table)
complete_chitest$stdres
#do we need to do this for each ID??
mosaicplot(complete_con_table, shade= T)


#create contingency tables for each ID
IDlist <- c(unique(sum_fac_survey$ID))
plot_table_list <- vector("list", length = length(IDlist))
names(plot_table_list) <- IDlist

for(n in 1:length(IDlist)) {
  plot <- sum_fac_survey |> 
    filter(ID == IDlist[n])
  
  plot_table <- plot |> 
    select(nurse_target_match, interaction) |> 
    filter(!is.na(nurse_target_match)) 
  plot_table <- table(plot_table)
  
  plot_table_list[[n]] <- plot_table
}


###Now we need to do the Chisquared tests on each contingency table
# Create an empty data frame to store Chi-squared test results
Chisq_results <- data.frame(ID = character(0), p_value = numeric(0), exclusion = character(0))
#Create list to put the residuals in 
Chisq_residuals <- vector("list", length = length(IDlist))
names(Chisq_residuals) <- IDlist
#create a list to put the expected values in 
Chisq_expected <- vector("list", length = length(IDlist))
names(Chisq_expected) <- IDlist

exclusion_criteria = " "
# Iterate through each contingency table
for (i in 1:length(plot_table_list)) {
  focus_table <- plot_table_list[[i]]
  
  if(length(focus_table) >1) { #only do test if there are at least two rows and two columns in the table
  
  # Perform the Chi-squared test
  chi_squared_result <- chisq.test(focus_table)
  
  # Extract the p-value from the test result
  p_value <- chi_squared_result$p.value
  
  #extract the expected results 
  exp <- chi_squared_result$expected
  #tell me if there is an expected value less than 5
  if(TRUE %in% c(c(1:5) %in% exp)) {
    exclusion_criteria <- paste(exclusion_criteria, "expected val less than 5")
  }
  
  #extract the standard residuals
  res <- chi_squared_result$residuals
  
  # Add the ID, p-value, 
  Chisq_results[i, 1] <- names(plot_table_list)[i]
  Chisq_results[i, 2] <- p_value
  Chisq_results[i, 3] <- exclusion_criteria
  
  #put the expected vals and residuals in a list
  Chisq_expected[[i]] <- exp
  Chisq_residuals[[i]] <- res
  
  }
  
  else { #if the table is empty it means we have no FGR for nurse or target species, do not do test
    
    exclusion_criteria <- paste(exclusion_criteria, "table length < 1", sep = "_")
    
    Chisq_results[i, 1] <- names(plot_table_list)[i]
    Chisq_results[i, 2] <- NA
    Chisq_results[i, 3] <- exclusion_criteria
    
    #do not put anything in the expected or residual lists. 
  }
  
  exclusion_criteria = " "                       
}##There are warnings that it the results may not be accurate??
#In chisq.test(focus_table) : Chi-squared approximation may be incorrect


##In which direction is the association for each ID?
sign_IDs <- Chisq_results |> 
  filter(p_value < 0.05,      #work with significant tests
         exclusion == " ") |> #ignore tests that meet exclusion criteria
  select(ID)
#only 5, can we do anything with that?

##Look at residuals for these 5 sign plots
Chisq_residuals[which(names(Chisq_residuals) %in% c(sign_IDs$ID))]




##Let's get CWM for each group####

##You have to use the output from trait_fill for CWM
#so now we have to do trait filling all over again

###
##COPIED FROM TRAIT FILLING:
#trait data 
FT_all_sites_raw <- read.csv("Functional trait data\\FT_all_sites.csv", row.names = 1) 

FT_all_sites_long <- FT_all_sites_raw|> 
  #calculate the C:N ratio
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN, Genus, Species)) |> 
  inner_join(FGR, by = "taxon") |>
  #make long format
  pivot_longer(cols = c("MeanLL","MeanSLA","MeanLDMC","MeanLA","MaxH","MaxLS","C_N_ratio"),
               names_to = "trait", values_to = "value")

FT_all_sites_long$ID <- as.factor(FT_all_sites_long$ID)
FT_all_sites_long$SITE_ID <- as.factor(FT_all_sites_long$SITE_ID)
FT_all_sites_long$GRAZ <- factor(FT_all_sites_long$GRAZ, levels = c(0,1,2,3))


#quadrat data
quad <- read.csv("Functional trait data\\quadrat_survey_all_plots.csv", row.names = 1)

quad_summary <- quad |> 
  mutate(quadrat = as.numeric(quadrat)) |> 
  #summarise the percent cover for each sp per plot, not per quadrat
  group_by(ID, taxon) |> 
  mutate(sum_plot_cover = sum(percent_cover)) |>  
  ungroup(taxon) |> 
  mutate(max_quadrats = max(quadrat)) |> 
  ungroup () |> 
  mutate(percent_cover_perplot = sum_plot_cover/max_quadrats) |> 
  select(!c(quadrat, percent_cover, sum_plot_cover, max_quadrats)) |> 
  distinct() |> 
  inner_join(FGR, by = "taxon") #only keep the species that are also in FGR

quad_summary$ID <- as.factor(quad_summary$ID)
quad_summary$SITE_ID <- as.factor(quad_summary$SITE_ID)
quad_summary$GRAZ <- factor(quad_summary$GRAZ, levels = c(0,1,2,3))


##Trait filling:
FT_filled <- trait_fill(
  comm = quad_summary, #community data with abundance values
  traits = FT_all_sites_long, #trait data
  abundance_col = "percent_cover_perplot", #the column with the abundance data
  taxon_col = "taxon", #column with the speciesnames (must be the same in both datasets)
  value_col = "value",
  trait_col = "trait",
  
  global = FALSE, #do not calculate traits at the global scale
  
  keep_all = FALSE, #do not keep trait data at all availible levels, only on the finest scale availible
  
  # specifies sampling hierarchy
  scale_hierarchy = c("COU", "SITE_ID", "ID"),
  
  other_col = c("Functional_group"),
  
  # min number of samples
  min_n_in_sample = 1 #if there is one trait value in the sample, do not search for values higher up in the hierarchy
)



#Get bootstrapped distributions
FT_bs <- trait_np_bootstrap(FT_filled, sample_size = 100, nrep = 50)

#Look at BS moments
FT_moments <- trait_summarise_boot_moments(
  FT_bs, 
  parametric = TRUE) #uses moment +- 1sd

#Now subset for the facilitation plots
#loop to get the IDs of plot in the fac data
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")

fac_IDs <- c(rep(NA, 3))

for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  if(t ==1) {
    
    fac_IDs <- c(unique(cou$ID))
    fac_aridities <- c(unique(cou$ARIDITY.v3))
    
  }else {
    
    temp_fac_IDs <- c(unique(cou$ID))
    fac_IDs <-c(fac_IDs, temp_fac_IDs)
  }
}

##First add siteinfo of each ID
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(c(ID, GRAZ, ARIDITY.v3))
siteinfo$ID <-as.factor(siteinfo$ID)

FT_moments_facplots <- FT_moments |> 
  filter(ID %in% fac_IDs) |> 
  inner_join(siteinfo, by = "ID")

write.csv(FT_moments_facplots, "Functional trait data\\results\\FT_moments_facplots_26Feb2024.csv")


#Get the CWm for each FGR, disregarding any other grouping variables
FGR_means <- FT_bs |> 
  ungroup() |> 
  group_by(Functional_group, trait) |> 
  summarise(mean_trait  = mean(mean))

FGR_means$Functional_group <- as.factor(FGR_means$Functional_group)
FGR_means |> 
  ggplot(aes(x = Functional_group, y = mean_trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(facets = "trait")
##Plot is difficult to see because the y axis is so wide. Standardise var to the same scale?


###Now plot how CWM of each trait changes with aridity
#C_N_ratio
FT_moments_facplots|> 
  filter(trait == "C_N_ratio") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 

#MAxH
FT_moments_facplots|> 
  filter(trait == "MaxH") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 

#MaxLS
FT_moments_facplots|> 
  filter(trait == "MaxLS") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 

#MeanLA
FT_moments_facplots|> 
  filter(trait == "MeanLA") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 

#MeanLDMC
FT_moments_facplots|> 
  filter(trait == "MeanLDMC") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 

#MeanLL
FT_moments_facplots|> 
  filter(trait == "MeanLL") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 

#MeanSLA
FT_moments_facplots|> 
  filter(trait == "MeanSLA") |>
  ggplot(aes(y = mean, x = ARIDITY.v3)) +
  geom_point() 




