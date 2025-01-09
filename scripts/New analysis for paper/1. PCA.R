#Script to do a PCA of trait data and get the independent traits with the highest loadings on the axes
library(tidyverse)
library(tidylog)
library(ggplot2)
library(ggbiplot)

#import trait data
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv",
               row.names = 1)

####Transform FT data to wide format and standardise trait vlaues####

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
  dplyr::select(!c(percentC, percentN))

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


####Do the principal component analysis####
all_FT_pca <- princomp(std_FT_wide, scores = T)
summary(all_FT_pca) #proportion of variance is the variance explained by the PC
all_FT_pca$scores #
all_FT_pca$loadings #How much each var contributed to building the PC
all_FT_pca$scale #scaling applied to each variable
all_FT_pca$center #means

#make biplot
biplot(all_FT_pca, choices = c("Comp.1", "Comp.2"))
plot(all_FT_pca$scores[, 1], all_FT_pca$scores[, 2])

##LDMC has highest loading on PC1, correlation is positive
##maxH has highest loading on PC2, correlation is positive

#So the two traits chosen are LDMC and maxH!


###GGplot biplot###
trait_pca <- ggbiplot(all_FT_pca, choices = c(1,2), 
         varname.size = 4, varname.color = "black") +
  geom_point(colour = "azure3", alpha = 0.8)+
  theme_classic() 
trait_pca$layers <- c(trait_pca$layers, trait_pca$layers[[2]]) #move the arrows to plot in the foreground
ggsave("trait_biplot.png", trait_pca, path = "Figures")
