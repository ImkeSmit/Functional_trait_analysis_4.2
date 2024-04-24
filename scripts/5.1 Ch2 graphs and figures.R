###CHAPTER 2 GRAPHS AND FIGURES###
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(multcomp)
library(multcompView)
library(glmmTMB)
library(car)
library(tidyverse)
library(tidylog)

####FD metrics accross aridity and grazing####
#import FD results
FD_results <- read.csv("Functional trait data\\results\\FD_results_4Mar2024.csv", row.names = 1)
FD_results$ID <- as.factor(FD_results$ID)
#add grazing and aridity
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(c(ID, GRAZ, ARIDITY.v3))
siteinfo$ID <-as.factor(siteinfo$ID)
#do the join
long_FD_results <- FD_results |> 
  inner_join(siteinfo, by = "ID") |> 
  select(!c(qual.FRic, nsp)) |> 
  pivot_longer(cols = c(FRic, FEve, FDiv, RaoQ), names_to = "FD_metrics", values_to = "value" )
long_FD_results$GRAZ <- as.factor(long_FD_results$GRAZ)

###RAOQ###
Raoq_aridity <- long_FD_results |> 
  filter(FD_metrics == "RaoQ") |> 
  ggplot(aes(y = log(value), x = ARIDITY.v3)) +
  geom_point()+
  theme_classic() +
  xlab("Aridity") +
  ylab("ln(RaoQ)")

Raoq_graz <- long_FD_results |> 
  filter(FD_metrics == "RaoQ") |> 
  ggplot(aes(y = log(value), x = GRAZ, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6)+
  scale_fill_manual(values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4")) +
  theme_classic() +
  ylab("") +
  xlab("Grazing pressure") +
  scale_x_discrete(labels = c("Ungrazed", "Low", "Medium", "High")) +
  theme(legend.position = "none")
  
Raoq_graz_arid <- ggarrange(Raoq_aridity, Raoq_graz, ncol = 2, nrow= 1, labels = c("a", "b"))
ggsave("Raoq_aridity_grazing.png", Raoq_graz_arid, height = 900, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")


###FDiv, FEve, FRic###
FD_arid_grad <- ggplot(long_FD_results, aes(y = value, x = ARIDITY.v3)) +
  geom_point() +
  theme_classic() +
  xlab("Aridity") +
  facet_wrap(~ factor(FD_metrics, levels = c("FRic", "FEve", "FDiv")), scales = "free_y")

FD_graz_grad <- ggplot(long_FD_results, aes(y = value, x = GRAZ, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  xlab("") +
  scale_x_discrete(labels = c("Ungrazed", "Low grazing", "Medium grazing", "High grazing")) +
  facet_wrap(~ factor(FD_metrics, levels = c("FRic", "FEve", "FDiv")), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank(), legend.position = "none")

FD_grad_plots <- ggarrange(FD_arid_grad, FD_graz_grad, ncol = 1, nrow= 2, labels = c("a", "b"))
ggsave("FD_grad_plots.png", FD_grad_plots, height = 1600, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")



##Nint ~ FD metrics####
#import nint results
#Add the NIntc
nint_result <- 
  read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1) |> 
  filter(ID %in% c(FD_results$ID))
nint_result$ID <- as.factor(nint_result$ID)


##summarise the nintc by plot
nint_sum <- nint_result |> 
  select(country, site_ID, ID, graz, aridity, NIntc_richness, NIntc_cover, NInta_richness, NInta_cover) |> 
  group_by(ID) |> 
  #calculate mean and sd of NIntc richness
  mutate(mean_NIntc_rich = mean(NIntc_richness, na.rm = TRUE), 
         se_NIntc_rich = sd(NIntc_richness, na.rm = TRUE)/sqrt(sum(!is.na(NIntc_richness)))) |> 
  #calculate the mean and sd of NIntc cover
  mutate(mean_NIntc_cov = mean(NIntc_cover, na.rm = TRUE), 
         se_NIntc_cov = sd(NIntc_cover, na.rm = T)/sqrt(sum(!is.na(NIntc_cover)))) |> 
  #calculate mean and sd of NInta richness
  mutate(mean_NInta_rich = mean(NInta_richness, na.rm = T), 
         se_NInta_rich = sd(NInta_richness, na.rm = T)/sqrt(sum(!is.na(NInta_richness)))) |> 
  #calculate the mean and sd of NInta cover
  mutate(mean_NInta_cov = mean(NInta_cover, na.rm = T), 
         se_NInta_cov = sd(NInta_cover, na.rm = T)/sqrt(sum(!is.na(NInta_cover)))) |>
  ungroup() |> 
  select(!c(NIntc_richness, NIntc_cover, NInta_richness, NInta_cover)) |> 
  distinct() |> #remove duplicate rows, only need one eman per plot
  left_join(FD_results, by = "ID") |>  #join to the FD_results
  filter(!is.na(FEve))


##NIntc against the RaoQ metrics##
rich_raoq <- ggplot(nint_sum, aes(x = log(RaoQ), y = mean_NIntc_rich)) +
  geom_point() +
  ylab(expression(NInt[C]~richness))+
  xlab("ln(RaoQ)") +
  geom_errorbar(aes(ymin = mean_NIntc_rich - se_NIntc_rich, ymax =  mean_NIntc_rich + se_NIntc_rich), width = 0.01) +
  theme_classic()

cov_raoq <- ggplot(nint_sum, aes(x = log(RaoQ), y = mean_NIntc_cov)) +
  geom_point() +
  ylab(expression(NInt[C]~cover))+
  xlab("ln(RaoQ)") +
  geom_errorbar(aes(ymin = mean_NIntc_cov - se_NIntc_cov, ymax =  mean_NIntc_cov + se_NIntc_cov), width = 0.01) +
  theme_classic()

nintc_Raoq <- ggarrange(rich_raoq, cov_raoq, ncol = 2, nrow = 1, labels = c("a", "b"))
ggsave("NIntc_RaoQ_scatterplots.png", nintc_Raoq, height = 900, width = 1300, units = "px", 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis clone\\Figures")


###Functional distance ~ GRAZ, ARIDITY, microsite affinty####
#import data
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)

#the model
twosp_dist$arid_sq <- (twosp_dist$ARIDITY.v3)^2
twosp2 <- glmmTMB(euclidean_dist ~ grouping*GRAZ +grouping*ARIDITY.v3 + 
                    ARIDITY.v3*GRAZ + grouping*arid_sq + arid_sq*GRAZ + (1|SITE_ID), data = twosp_dist)
Anova(twosp2, type = 2) #type 2 is seq SS
summary(twosp2)

cld(glht(model = twosp2, mcp(GRAZ = "Tukey")))
emmeans(twosp2, specs = "GRAZ")

##distance~grouping, colour by graz
#we want to know which levels of graz are significantly different from each other, within gspecific groupings
#concateante graz and grouping to create treatemnt ID
twosp_dist <- twosp_dist |> 
  mutate(treatmentID = str_c(grouping, GRAZ, sep = "_"))
twosp_dist$treatmentID <- as.factor(twosp_dist$treatmentID)

treatment_mod <- lm(euclidean_dist ~ treatmentID, data = twosp_dist)
Anova(treatment_mod)
letters <- cld(glht(model = treatment_mod, mcp(treatmentID = "Tukey")))
#look at letters and make a vector to of them to add to the graph
#put them in the same order as the boxplots on the graph
# Create a data frame with the letters and corresponding x positions
letters_data <- data.frame(grouping = rep(unique(twosp_dist$grouping), each = length(unique(twosp_dist$GRAZ))),
                           GRAZ = rep(as.factor(c(0,1,2,3))), times = length(unique(twosp_dist$grouping)),
                           letter = c("a", "bc", "cde", "ef", "ag", "bc", "ef", "ce", "af", "b", "bd", "efg"), 
                           ycoord = c(rep(11.5, 12)))

dist_bygrouping_graz <- ggplot(twosp_dist, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(labels = c("ungrazed", "low", "medium", "high"), 
                    values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  ylim(0, 12) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Microsite affinity of target species", 
       fill = "") +
  scale_x_discrete(labels = c("bare", "dominant", "both")) +
  geom_text(data = letters_data, aes(x = grouping, y = ycoord, label = letter, group = GRAZ),
            position = position_dodge(width = 0.75), vjust = -0.5, hjust = 0.5, size = 4, color = "black") +
  theme_classic() +
  theme(legend.position = "bottom") 

ggsave("dist_grouping_graz.png", dist_bygrouping_graz, height = 1500, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#distance ~graz
twosp1 <- glmmTMB(euclidean_dist ~ grouping +GRAZ + ARIDITY.v3 + arid_sq 
                  + (1|SITE_ID), data = twosp_dist)
summary(twosp1)

letters <- cld(glht(model = twosp1, mcp(GRAZ = "Tukey")))
letters_df <- data.frame(GRAZ = c("0", "1", "2", "3"), letters = c("a", "b", "c", "d"), ycoord = c(12, 12, 12, 12))


dist_graz<- ggplot(twosp_dist, aes(x = GRAZ, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  ylim(0, 12) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Grazing pressure", 
       fill = "") +
  scale_x_discrete(labels = c("ungrazed", "low", "medium", "high")) +
  geom_text(data = letters_df, aes(x = GRAZ, y = ycoord, label = letters)) +
  theme_classic() +
  theme(legend.position = "none") 
ggsave("dist_graz.png", dist_graz, 
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


##distance~aridity
#Get model predictions
model_forgraph <- glmmTMB(euclidean_dist ~ grouping*ARIDITY.v3*arid_sq, data = twosp_dist)
Anova(model_forgraph, type = 2)

pred_df <- data.frame(ARIDITY.v3 = c(rep(unique(twosp_dist$ARIDITY.v3), 12)), 
                      arid_sq = c(rep(unique(twosp_dist$arid_sq), 12)),
                      grouping = as.factor(c(rep("nurse_bare_only", length(unique(twosp_dist$ARIDITY.v3))*4), 
                                   rep("nurse_nurse_only", length(unique(twosp_dist$ARIDITY.v3))*4), 
                                   rep("nurse_both", length(unique(twosp_dist$ARIDITY.v3))*4)))) 


twosp2_predictions <- predict(model_forgraph, pred_df)
twosp2_predictions <- pred_df |> 
  add_column(predictions = c(twosp2_predictions))

#dist~aridity
line_color <- brewer.pal(8, "Dark2")[7]
dist_aridity <- ggplot(twosp_dist, aes(x = ARIDITY.v3, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  #geom_smooth(method = "lm", color = line_color)+
  #add model predictions to the plot
  geom_line(data = twosp2_predictions, aes(x = ARIDITY.v3, y = predictions), lty = 1, lwd = 1, color = line_color) +
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity") +
  theme_classic()

ggsave("dist_aridity.png", dist_aridity, height = 1400, width = 1100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#dist~aridity squared
line_color <- brewer.pal(8, "Dark2")[7]
dist_aridsq <- ggplot(twosp_dist, aes(x = arid_sq, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  #geom_smooth(method = "lm", color = line_color)+
  #add model predictions to the plot
  geom_line(data = twosp2_predictions, aes(x = arid_sq, y = predictions), lty = 1, lwd = 1, color = line_color) +
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity^2") +
  theme_classic()

ggsave("dist_aridity.png", dist_aridity, height = 1400, width = 900, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")



###ONLY FOR FACILITATIVE REPS: Functional distance ~ GRAZ, ARIDITY, microsite affinty####
#import data
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp.csv", row.names = 1) |> 
  mutate(aridsq = ARIDITY.v3^2)
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)

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


#we want to know which levels of graz are significantly different from each other, within specific groupings
#concateante graz and grouping to create treatemnt ID
twosp_dist_fac <- twosp_dist_fac |> 
  mutate(treatmentID = str_c(grouping, GRAZ, sep = "_"))
twosp_dist_fac$treatmentID <- as.factor(twosp_dist_fac$treatmentID)

fac_treatment_mod <- lm(euclidean_dist ~ treatmentID, data = twosp_dist_fac)
Anova(fac_treatment_mod)
fac_letters <- cld(glht(model = fac_treatment_mod, mcp(treatmentID = "Tukey")))
#look at letters and make a vector to of them to add to the graph
#put them in the same order as the boxplots on the graph
# Create a data frame with the letters and corresponding x positions
fac_letters_data <- data.frame(grouping = rep(unique(twosp_dist_fac$grouping), each = length(unique(twosp_dist_fac$GRAZ))),
                           GRAZ = rep(as.factor(c(0,1,2,3))), times = length(unique(twosp_dist_fac$grouping)),
                           letter = c("a", "bc", "cd", "de", "ef", "bc", "bd", "cd", "af", "c", "bc", "df"), 
                           ycoord = c(rep(11.5, 12)))


#dist~grouping, colour by graz
FAC_dist_bygrouping_graz <- ggplot(twosp_dist_fac, aes(x = grouping, y = euclidean_dist, fill = GRAZ)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(labels = c("ungrazed", "low", "medium", "high"), 
                    values = c("darkgreen", "chartreuse2" , "darkolivegreen3", "darkgoldenrod4", "azure4" )) +
  ylim(0, 12) +
  labs(y = "Euclidean distance between dominant and target species", 
       x = "Microsite affinity of target species", 
       fill = "") +
  scale_x_discrete(labels = c("bare", "dominant", "both")) +
  geom_text(data = fac_letters_data, aes(x = grouping, y = ycoord, label = letter, group = GRAZ),
            position = position_dodge(width = 0.75), vjust = -0.5, hjust = 0.5, size = 4, color = "black") +
  theme_classic() +
  theme(legend.position = "bottom") 
ggsave("facilitative_reps_dist_grouping_graz.png", FAC_dist_bygrouping_graz, height = 1500, width = 1300, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")

##distance~aridity
#Get model predictions
fac_model_forgraph <- glmmTMB(euclidean_dist ~ grouping*ARIDITY.v3*aridsq, data = twosp_dist_fac)
Anova(fac_model_forgraph, type = 2)

pred_df <- data.frame(ARIDITY.v3 = c(rep(unique(twosp_dist_fac$ARIDITY.v3), 12)), 
                      aridsq = c(rep(unique(twosp_dist_fac$aridsq), 12)),
                      grouping = as.factor(c(rep("nurse_bare_only", length(unique(twosp_dist_fac$ARIDITY.v3))*4), 
                                             rep("nurse_nurse_only", length(unique(twosp_dist_fac$ARIDITY.v3))*4), 
                                             rep("nurse_both", length(unique(twosp_dist_fac$ARIDITY.v3))*4)))) 


twosp_fac_predictions <- predict(fac_model_forgraph, pred_df)
twosp_fac_predictions <- pred_df |> 
  add_column(predictions = c(twosp_fac_predictions))


#dist~aridity
line_color <- brewer.pal(8, "Dark2")[7]
fac_dist_aridity <- ggplot(twosp_dist_fac, aes(x = ARIDITY.v3, y = euclidean_dist)) +
  geom_jitter(width = 0.08, height = 0.2, alpha = 0.6) +
  facet_wrap(~grouping, ncol = 1,
             labeller = as_labeller(c("nurse_bare_only" = "bare", "nurse_both" = "both", "nurse_nurse_only" = "nurse")))+
  #geom_smooth(method = "lm", color = line_color)+
  #add model predictions to the plot
  geom_line(data = twosp_fac_predictions, aes(x = ARIDITY.v3, y = predictions), lty = 1, lwd = 1, color = line_color) +
  labs(y = "Euclidean distance between dominant and target species", x = "Aridity") +
  theme_classic()

ggsave("facilitative_reps_dist_aridity.png", fac_dist_aridity, height = 1400, width = 1100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")




###Biplot####
#import standardised FT data to use in pca
FT_wide_FGR <- read.csv("Functional trait data//results//standard_FT_FGR.csv", row.names = 1)
FT_wide_FGR$Functional_group <- as.factor(FT_wide_FGR$Functional_group) 

base_pca <- princomp(FT_wide_FGR[ , -c(8,9)], scores = T)

pal <- brewer.pal(3, "Dark2")
#PC 1 and 2 with FGR
pc1_2 <- ggbiplot(base_pca, choices = c(1,2), 
                  groups = FT_wide_FGR$Functional_group, 
                  ellipse = FALSE, 
                  var.axes = T) +
  scale_color_manual(labels = c("1 - Large plants", "2 - Tough plants", "3 - Palatable plants"), 
                     values = pal)+
  labs(x = "PC 1 (29.5%)", y = "PC 2 (22.2%)", color = "Functional groups") +
  xlim(-3, 5)+
  ylim(-3, 5 )+
  theme_classic()

ggsave("pc1_2_biplot.png", pc1_2, height = 1500, width = 2100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#PC 1 and 2 without FGR
pc1_2_nogroups <- ggbiplot(base_pca, choices = c(1,2), var.axes = T, alpha = 0.5, 
                  varname.color = "darkolivegreen4", varname.adjust = 1.2, varname.size = 6) +
  labs(x = "PC 1 (29.5%)", y = "PC 2 (22.2%)") +
  xlim(-3, 5)+
  ylim(-3, 5 )+
  theme_classic()

ggsave("pc1_2_biplot_nogroups.png", pc1_2_nogroups, height = 1500, width = 2100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


#PC 1 and 3
pc1_3 <- ggbiplot(base_pca, choices = c(1,3), 
                  groups = FT_wide_FGR$Functional_group, 
                  ellipse = FALSE) +
  scale_color_manual(labels = c("1 - Large plants", "2 - Tough plants", "3 - Palatable plants"), 
                     values = pal)+
  labs(x = "PC 1 (29.5%)", y = "PC 3 (17.6%)", color = "Functional groups") +
  xlim(-3,4)+
  theme_classic()

ggsave("pc1_3_biplot.png", pc1_3,height = 1500, width = 2100, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")


