##Functional distance analysis
library(tidyverse)
library(tidylog)
library(ggplot2)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(MuMIn)
library(multcomp)
library(multcompView)

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


##Get the species that are nurses, growing with nurses, or in the bare microsite in each country####
#We require raw country data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}

countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")

##LOOP starts here
l = 1
for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  IDlist <- unique(cou$ID)
  
    for(i in 1:length(IDlist)) {
    plot <- cou |> 
      filter(ID == IDlist[i])
    
    replist <- unique(plot$Number.of.replicate)
  
      for(r in 1:length(replist)) {
        one_rep <- plot |> 
          filter(Number.of.replicate == replist[r])
        
        NURSE <- one_rep |> 
          filter(Microsite == 2) |> 
          distinct(ID_Microsite)
        
        bare_sp <- one_rep |> #sp in bare microsite
          filter(Microsite == 1) |> 
          select(Species.within.quadrat)
        
        nurse_sp <- one_rep |> #sp in nurse microsite
          filter(Microsite == 2) |> 
          select(Species.within.quadrat)
        
        #get the species that are in both nurse and bare microsites
        both_sp <- c(bare_sp[match(bare_sp$Species.within.quadrat, nurse_sp$Species.within.quadrat) , ])
        both_sp <- both_sp[!is.na(both_sp)]
        #if there are no species in both microsites, assign NA
        if(length(both_sp) == 0) {
          both_sp <- NA
        }
        
        #only if both_sp is not NA:
        if(FALSE %in% is.na(both_sp)) {
        
        #species ONLY in nurse microsite:
        nurse_only <- nurse_sp[-which(nurse_sp$Species.within.quadrat %in% both_sp), ]
        #if there are no species, assign NA
        if(length(nurse_only) == 0) {
          nurse_only <- NA
        }
        
        #species ONLY in bare microsites
        bare_only <- bare_sp[-which(bare_sp$Species.within.quadrat %in% both_sp) , ]
        #if there are no species, assign NA
        if(length(bare_only) == 0) {
          bare_only <- NA
        }
        
        #if both sp is NA:
        }else {
          nurse_only <- nurse_sp$Species.within.quadrat
          bare_only <- bare_sp$Species.within.quadrat
        }
        
        #only for the first run of the loop:
        if(l == 1){
        #put species and their classifications in the dataframe
        nurse_identity_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = NURSE$ID_Microsite, position = "nurse_species")
        nurse_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = nurse_only, position = "nurse_only")
        bare_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = bare_only, position = "bare_only")
        both_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = both_sp, position = "both")
        
        sp_positions <- rbind(nurse_identity_df, nurse_df, bare_df, both_df)
        
        } else { #for subsequent runs we rbind to the df we made whn l was 1
          nurse_identity_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = NURSE$ID_Microsite, position = "nurse_species")
          nurse_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = nurse_only, position = "nurse_only")
          bare_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = bare_only, position = "bare_only")
          both_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = both_sp, position = "both")
          
          sp_positions <- rbind(sp_positions, nurse_identity_df, nurse_df, bare_df, both_df)
        }
        
        l = l+1
      }#close loop through replicates
  }#close loop through ID's
}#close loop through countries

sp_positions <- as.data.frame(sp_positions)
#save the table 
write.csv(sp_positions, "Functional trait data\\results\\sp_positions.csv")



###Retreive distances between the nurse and each species also growing in the rep####
#then classify that distance as D_bare_only, D_nurse_only or D_both
#This distance is only between 2 sp
#write the pairwise_fdist function to do this
pairwise_fdist <- function(distmat, sp_positions) {
  
  IDlist <- unique(sp_positions$ID)
  
  twosp_dist <- cbind("ID" = character(), "replicate" = character(), 
                      "euclidean_dist" = numeric(), "grouping" = character(), "nurse" = character(), "target" = character())
  
  for(i in 1:length(IDlist)) {
    plot <- sp_positions |> 
      filter(ID == IDlist[i])
    
    replist <- unique(plot$replicate) 
    
    for(r in 1:length(replist)) {
      #Get the names of nurse, bare only, nurse only and both species
      NURSE <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "nurse_species") |> 
        select(taxon)
      NURSE <- NURSE$taxon
      
      bare_only <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "bare_only") |> 
        select(taxon)
      bare_only <- c(bare_only$taxon)
      
      nurse_only <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "nurse_only") |> 
        select(taxon)
      nurse_only <- c(nurse_only$taxon)
      
      both <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "both") |> 
        select(taxon)
      both <- c(both$taxon)
      
      ###
      for(b in 1:length(bare_only)) {
        D_bare_only <- 
          mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == bare_only[b])])
        
        if(b == 1){
          twosp_bare_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                   "euclidean_dist" = D_bare_only, "grouping" = "nurse_bare_only", 
                                   "nurse" = NURSE, "target" = bare_only[b])
        } else {
          temp_twosp_bare_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                        "euclidean_dist" = D_bare_only, "grouping" = "nurse_bare_only", 
                                        "nurse" = NURSE, "target" = bare_only[b])
          twosp_bare_only <- rbind(twosp_bare_only, temp_twosp_bare_only)
        }}
      
      ###
      for(n in 1:length(nurse_only)) {
        D_nurse_only <- 
          mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == nurse_only[n])])
        
        if(n == 1){
          twosp_nurse_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                    "euclidean_dist" = D_nurse_only, "grouping" = "nurse_nurse_only", 
                                    "nurse" = NURSE, "target" = nurse_only[n])
        } else {
          temp_twosp_nurse_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                         "euclidean_dist" = D_nurse_only, "grouping" = "nurse_nurse_only", 
                                         "nurse" = NURSE, "target" = nurse_only[n])
          twosp_nurse_only <- rbind(twosp_nurse_only, temp_twosp_nurse_only)
        }}
      
      ###
      for(z in 1:length(both)) {
        D_both <- 
          mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == both[z])])
        
        if(z == 1){
          twosp_both <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                              "euclidean_dist" = D_both, "grouping" = "nurse_both", 
                              "nurse" = NURSE, "target" = both[z])
        } else {
          temp_twosp_both <- twosp_both <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                                 "euclidean_dist" = D_both, "grouping" = "nurse_both", 
                                                 "nurse" = NURSE, "target" = both[z])
          twosp_both <- rbind(twosp_both, temp_twosp_both)
        }}
      
      twosp_dist <- rbind(twosp_dist, twosp_bare_only,twosp_nurse_only, twosp_both)
    }
  }
  #The distance is NaN if nurse_only, bare_only or both is NA
  
  return(twosp_dist)
}


####run the function
#Get the euclidean distance between each pair of species
distmat <- as.matrix(dist(std_FT_wide, method = "euclidean"))

#import sp_positions
sp_positions <- read.csv("Functional trait data\\results\\sp_positions.csv", row.names = 1) 

twosp_dist <- pairwise_fdist(distmat = distmat, sp_positions = sp_positions)

#Add ardidty and graz
twosp_dist <- as.data.frame(twosp_dist)
twosp_dist$ID <- as.numeric(twosp_dist$ID)
#import siteinfo
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  filter(!is.na(ID))

#do the join
twosp_dist <- twosp_dist |> 
  filter(!is.na(euclidean_dist), !euclidean_dist == "NaN") |> 
  inner_join(siteinfo, by = "ID") 

write.csv(twosp_dist, "Functional trait data\\results\\Functional_distances_between_2sp.csv")


####Models of dist ~ association####
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)
twosp_dist$arid_sq <- (twosp_dist$ARIDITY.v3)^2


###Lets join the results of the CHi2 tests to twosp-dist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

dist_ass_join <- twosp_dist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) #only work with these associations
dist_ass_join$association <- as.factor(dist_ass_join$association)
dist_ass_join$nurse <- as.factor(dist_ass_join$nurse)
dist_ass_join$SITE_ID <- as.factor(dist_ass_join$SITE_ID)
dist_ass_join$euclidean_dist <- as.numeric(dist_ass_join$euclidean_dist)

##does the distance between nurses and bare associated species differ from the distance between nurses and nurse associated species?
dist_ass_null <- glmmTMB(euclidean_dist ~ 1 + (1|nurse) + (1|SITE_ID), data = dist_ass_join)

dist_ass_mod <- glmmTMB(euclidean_dist ~ association + (1|nurse) + (1|SITE_ID), data = dist_ass_join)
summary(dist_ass_mod)
Anova(dist_ass_mod)
anova(dist_ass_null, dist_ass_mod) #p = 0.1645 

emmeans(dist_ass_mod, specs = "association")
cld(glht(model = dist_ass_mod, mcp(association = "Tukey")))

#model diagnostics
simres <- simulateResiduals(dist_ass_mod)
plot(simres)
#a little underdispersed
#HOV violated


ggplot(dist_ass_join, aes(x = association, y = euclidean_dist)) +
  geom_boxplot() +
  ylab("Euclidean distance between dominant and target species") +
  xlab("target species association") +
  annotate(geom = "text", x = unique(dist_ass_join$association), y = c(12,12,12,12) , 
           label = c("b", "b", "a", "ab")) + #put letters in order of dist$association
  theme_classic()

##Nurse and bare are not significantly different, so the difference in the nurse and target traits do not matter for facilitation


###One-dimensional (trait) distance between species####
#get the functional distance between species in terms of one trait
##SLA
sla_distmat <- as.matrix(dist(std_FT_wide[, which(colnames(std_FT_wide) == "MeanSLA")], method = "euclidean"))




###Old code and models below####

#model without interactions
twosp1 <- glmmTMB(euclidean_dist ~ grouping +GRAZ + ARIDITY.v3 + arid_sq 
                     + (1|SITE_ID), data = twosp_dist)
summary(twosp1)
Anova(twosp1) #only grazing is significant
AIC(twosp1) #36140.13



#model with interactions 
twosp2 <- glmmTMB(euclidean_dist ~ grouping*GRAZ +grouping*ARIDITY.v3 + 
                    ARIDITY.v3*GRAZ + grouping*arid_sq + arid_sq*GRAZ + (1|SITE_ID), data = twosp_dist)
Anova(twosp2)
summary(twosp2)
AIC(twosp2) #35883.66
emmeans(twosp2, specs = "grouping")
emmeans(twosp2, specs = "GRAZ")

twosp2_simres <- simulateResiduals(twosp2)
plot(twosp2_simres)

twosp_null <- glmmTMB(euclidean_dist ~ 1 + (1|SITE_ID), data = twosp_dist)
anova(twosp_null, twosp2) #whole model p = < 2.2e-16

#get R squared
r.squaredGLMM(twosp2) #take the theoretical


##Let's run the model only for replicates with Nintc_richness > 0
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

#model for only facilitative interactions
twosp_fac <- glmmTMB(euclidean_dist ~ grouping*GRAZ +grouping*ARIDITY.v3 + 
                       ARIDITY.v3*GRAZ + grouping*arid_sq + arid_sq*GRAZ + (1|SITE_ID), data = twosp_dist)

Anova(twosp_fac)
summary(twosp_fac)
  


#PLOTS#
#dist ~ graz
ggplot(twosp_dist, aes(x = GRAZ, y = euclidean_dist)) +
  geom_boxplot()

#dist~graz, colour by grouping
ggplot(twosp_dist, aes(x = GRAZ, y = euclidean_dist, fill = grouping)) +
  geom_boxplot()

#dist~grouping, colour by graz
ggplot(twosp_dist, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot()


#dist~aridity
ggplot(twosp_dist , aes(x = ARIDITY.v3, y = euclidean_dist, color = grouping)) +
  geom_jitter(width = 0.1, height = 0.1)
