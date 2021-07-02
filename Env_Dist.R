#Quantitative Biology 
#Measures of Biodiversity 
#Author: Benjamin Gazeau
#Date: 01 July 2021


library(vegan)
library(ggplot2)
library(geodist)
library(ggpubr) 
library(readr)

xyz <- read.csv("~/Desktop/Quantitative_Ecology-main/exercises/diversity/Euclidian_distance_demo_data_xyz.csv")
dim(xyz)
xyz 

#the vegdist() function in vegan package to calculate euclidean distances 
#it applies pythagoras theorem 
#'upper = FALSE' eliminates the display of upper triangle in square matrix 
#'diag = TRUE' shows the calculation over the same sites causing a zero diagonal over the matrix 
xyz_euc <- round(vegdist(xyz[, 2:4], method = "euclidian", upper = FALSE, diag = TRUE), 4) # select only cols 2, 3 and 4
xyz_euc 
#[, 2:4] - take note how ratio is after the comma - refers to column selection - sites column removed 

xyz_df <- as.data.frame(round(as.matrix(xyz_euc), 4))
xyz_df #presenting matrix as a data frame 

#Analyzing the seaweed environmental data 

load("~/Desktop/Quantitative_Ecology-main/exercises/diversity/SeaweedEnv.RData")
dim(env)

round(env[1:5, 1:5], 4)
round(env[(nrow(env) - 5):nrow(env), (ncol(env) - 5):ncol(env)], 4)

colnames(env)



#selecting the environmental that we are interested in 
env1 <- dplyr::select(env, febMean, febRange, febSD, augMean,
                      augRange, augSD, annMean, annRange, annSD)


E1 <- round(decostand(env1, method = "standardize"), 4)
E1[1:5, 1:5]


E1_euc <- round(vegdist(E1, method = "euclidian", upper = TRUE), 4)
E1_df <- as.data.frame(as.matrix(E1_euc))
E1_df[1:10, 1:10] # the first 10 rows and columns

ggplot(data = E1_df, (aes(x = 1:58, y = `1`))) +
  ggtitle("Graph showing environmental dissimilarity") +
  theme(plot.title = element_text( colour = "forestgreen", hjust = 0.5, face = "bold")) +
  geom_line(colour = "forestgreen") + xlab("Coastal section, west to east") + ylab("Environmental distance")

#Euclidean distances of geographical data 

geo <- read.csv("~/Desktop/Quantitative_Ecology-main/exercises/diversity/sites.csv")
dim(geo) 
head(geo)

dists <- geodist(geo, paired = TRUE, measure = "geodesic")
dists_df <- as.data.frame(as.matrix(dists))
colnames(dists_df) <- seq(1:58) #dataframe to allow for graphic displays 
dists_df[1:5, 1:5]

plt1 <- ggplot(data = dists_df, (aes(x = 1:58, y = `1`/1000))) +
  ggtitle("Graph showing actual distance across geographic points") +
  theme(plot.title = element_text( colour = "red", hjust = 0.5, face = "bold")) +
  geom_line( colour = "red") + xlab("Coastal section, west to east") + ylab("Distance (km)")
plt1

dists_euc <- vegdist(geo, method = "euclidian")
dists_euc_df <- round(as.data.frame(as.matrix(dists_euc)), 4)
dists_euc_df[1:5, 1:5]
library(ggthemes)
plt2 <- ggplot(data = dists_euc_df, (aes(x = 1:58, y = `1`))) +
  ggtitle("Graph showing euclidian distance across geographic points") +
  theme(plot.title = element_text( colour = "black", hjust = 0.5, face = "bold")) +
  geom_line(colour = "black") + xlab("Coastal section, west to east(km)") + ylab("Euclidian distance")
plt2
ggarrange(plt1, plt2, nrow = 2) #This shows both of the graphs  
