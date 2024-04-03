library(tidyverse)
library(tidylog)
library(ggplot2)
library(ggpubr)

#Do PCA with all species
#From the filled trait data for plotspecific species
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv", row.names = 1)
##!remember that these species were not necessarily filled from the same graz level
#!redo it with that filling??

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

all_FT_pca <- princomp(std_FT_wide, scores = T)
summary(all_FT_pca) #proportion of variance is the variance explained by the PC
all_FT_pca$scores #
all_FT_pca$loadings #How much each var contributed to building the PC
all_FT_pca$scale #scaling applied to each variable
all_FT_pca$center #means
biplot(all_FT_pca, choices = c("Comp.1", "Comp.2"))

plot(all_FT_pca$scores[, 1], all_FT_pca$scores[, 2])


#Now let's get the species that occur in each grazing level, and which of those are nurses and targets
#all the species in the trait data
sp_graz <- FT |> 
  group_by(GRAZ) |> 
  distinct(taxon) |> 
  ungroup()

#get the species that were nurses in each grazing level
#import the nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1)   

graz_nurses <- nint_result |> 
  group_by(graz) |> 
  distinct(nurse) |> 
  ungroup()


#now get whether each species in scores is present at each grazing level, and whether they are nurses or not
scores <- data.frame(all_FT_pca$scores)
scores$taxon <- row.names(scores)
scores$graz0presence <- NA
scores$graz1presence <- NA
scores$graz2presence <- NA
scores$graz3presence <- NA

grazlevels <- c(0,1,2,3)
graz_presence_columns <- c("graz0presence", "graz1presence", "graz2presence", "graz3presence")

for(g in 1:length(grazlevels)) {
  nurses <- graz_nurses |> 
    filter(graz == grazlevels[g])
  
  allsp <- sp_graz |> 
    filter(GRAZ == grazlevels[g])
  
  targets <- allsp[-which(allsp$taxon %in% nurses$nurse) , ]

for(i in 1:nrow(scores)) {
  if(scores$taxon[i] %in% nurses$nurse) {
    scores[i, which(colnames(scores) == graz_presence_columns[g])] <- "nurse" }

  if(scores$taxon[i] %in% targets$taxon) {
    scores[i, which(colnames(scores) == graz_presence_columns[g])] <- "target" }
}}



#plot the nurses and targets in graz0:
pal <- brewer.pal(3, "Dark2")[1:2]

graz0_biplot <- scores |> 
  filter(!is.na(graz0presence)) |> 
  ggplot(aes(x = Comp.1, y = Comp.2, color = graz0presence)) +
  geom_point(alpha = 0.6)+
  scale_color_manual(values = pal) +
  labs(color = " ") +
  theme_classic() 

#plot the nurses and targets in graz1:
graz1_biplot <- scores |> 
  filter(!is.na(graz1presence)) |> 
  ggplot(aes(x = Comp.1, y = Comp.2, color = graz1presence)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = pal) +
  labs(color = " ") +
  theme_classic() 

graz2_biplot <- scores |> 
  filter(!is.na(graz2presence)) |> 
  ggplot(aes(x = Comp.1, y = Comp.2, color = graz2presence)) +
  geom_point(alpha = 0.6)+
  scale_color_manual(values = pal) +
  labs(color = " ") +
  theme_classic() 

graz3_biplot <- scores |> 
  filter(!is.na(graz3presence)) |> 
  ggplot(aes(x = Comp.1, y = Comp.2, color = graz3presence)) +
  geom_point(alpha = 0.6)+
  scale_color_manual(values = pal) +
  labs(color = " ") +
  theme_classic()

pca_bygraz <- ggarrange(graz0_biplot, graz1_biplot, graz2_biplot, graz3_biplot, 
          common.legend = TRUE,
          nrow = 2, ncol = 2, labels = c("ungrazed", "low", "medium", "high"), 
          hjust = c(-0.7, -1.7, -0.8, -1.4), 
          font.label = list(size = 12), legend= "bottom")

ggsave("PCA_bygraz.png", pca_bygraz, height = 1500, width = 1500, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


