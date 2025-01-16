###Script to make graphs and figures for the paper####
library(ggplot2)
library(ggpubr)

####AKAIKE WEIGHTS FIGURE####
#import akaike weights of nint richness model
rich_sw <- read.csv("Functional trait data\\paper results\\nint_rich_model_weights.csv") |> 
  rename(var = X) |> 
  mutate(var2 = c("AMT", "graz", "H", "LDMC", "pH", "RASE", "SAC", "latitude", "longitude", 
                  "graz x LDMC", "graz x RASE", "graz x SAC", "H x pH", "LDMC x pH", "aridity", "LDMC x RASE", 
                  "AMT x LDMC", "AMT x H", "H x RASE", "LDMC x SAC", "H x SAC"))
rich_sw$var2 <- reorder(rich_sw$var2, rich_sw$importance)

#import akaike weights of nint cover model
cov_sw <- read.csv("Functional trait data\\paper results\\nint_cov_model_weights.csv") |> 
  rename(var = X, importance = cov_importance) |> 
  mutate(var2 = c("graz", "LDMC", "pH", "RASE", "SAC", "latitude", "longitude", 
                  "graz x LDMC", "graz x RASE", "graz x SAC", "LDMC x pH", "H", "H x pH", "AMT", "AMT x H", 
                  "aridity", "aridity x graz", "LDMC x SAC", "graz x H", "AMT x LDMC",  "LDMC x RASE", 
                  "H x RASE", "H x SAC"))
cov_sw$var2 <- reorder(cov_sw$var2, cov_sw$importance)


rich_sw_bar <- ggplot(aes(x = var2, y = importance), data = rich_sw) +
  geom_bar(stat = "identity") +
  xlab("Variable") +
  ylab("Importance") +
  coord_flip()+
  theme_classic()


cov_sw_bar <- ggplot(aes(x = var2, y = importance), data = cov_sw) +
  geom_bar(stat = "identity") +
  xlab("Variable") +
  ylab("Importance") +
  coord_flip()+
  theme_classic()

sw_combo <- ggarrange(rich_sw_bar, cov_sw_bar, nrow = 1, ncol = 2, labels = c("a", "b"))
ggsave("sw_bar.png", plot = sw_combo, path = "Figures", height = 1000, width = 2000, units = 'px')


