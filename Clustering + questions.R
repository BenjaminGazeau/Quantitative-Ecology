#Quantitative Biology 
#Clustering
#Author: Benjamin Gazeau
#Date: 08 August 2021

# Load Packages -----------------------------------------------------------
  
library(tidyverse) 
library(cluster)
library(ggcorrplot)
library(factoextra)
library(vegan)
library(ggpubr)
library(Rcpp)
library(ggplot2)
install.packages('Rcpp') 


# Load Data ---------------------------------------------------------------

SDGs <- read_csv("SDG_complete.csv")
SDGs[1:5, 1:8] # View first 5 rows and first 8 columns 

# The parent locations:
unique(SDGs$ParentLocation)

# The number of countries:
length(SDGs$Location)

# Correlation analysis of explanatory variables 
corr <- round(cor(SDGs[3:ncol(SDGs)]), 1)
ggcorrplot(corr, type = 'upper', outline.col = "white", 
           colors = c("navy", "white", "#FC4E07"), 
           lab = TRUE)
# there are highly correlated variables, but let's continue (not ideal)

# Standardize variables
SDGs_std <- decostand(SDGs[3:ncol(SDGs)], method = "standardize")
rownames(SDGs_std) <- SDGs$Location # carry location names into output


# Use the factoextra package to try and decide how many clusters to use

# using silhouette analysis
plt1 <- fviz_nbclust(SDGs_std, cluster::pam, method = "silhouette") + theme_grey()

# total within cluster sum of square / elbow analysis
plt2 <- fviz_nbclust(SDGs_std, cluster::pam, method = "wss") + theme_grey()

# gap statistics
plt3 <- fviz_nbclust(SDGs_std, cluster::pam, method = "gap_stat") + theme_grey()

ggarrange(plt1, plt2, plt3, nrow = 3)


# Let's proceed with 3 clusters (as 2 seem insufficient)
SDGs_pam <- pam(SDGs_std, metric = "euclidean", k = 3)

fviz_cluster(SDGs_pam, geom = "point", ellipse.type = "convex", palette = c("#FC4E07", "violetred3", "deepskyblue3"), ellipse.alpha = 0.05) +
  geom_text(aes(label = SDGs$Location), size = 2.5)


# We cannot clearly see SA, so let's edit the plot
# scale SA bigger for plotting
SDGs <- SDGs |> 
  mutate(col_vec = ifelse(Location == "South Africa", "black", "grey50"),
         scale_vec = ifelse(Location == "South Africa", 3.5, 2.5))

fviz_cluster(SDGs_pam, geom = "point", ellipse.type = "convex",
             palette = c("#FC4E07", "violetred3", "deepskyblue3"),
             ellipse.alpha = 0.05, pointsize = 2.0) +
  geom_text(aes(label = SDGs$Location), size = SDGs$scale_vec, col = SDGs$col_vec)

# Let's make it a star plot
fviz_cluster(SDGs_pam, palette = c("#FC4E07", "violetred3", "deepskyblue3"), ellipse.type = "euclid", 
             star.plot = TRUE, repel = TRUE, pointsize = SDGs$scale_vec * 0.8) + # SA, no 147, plotted slightly bigger
  theme_grey()


# Do silhouette analysis to check cluster fidelity
fviz_silhouette(SDGs_pam, palette = c("#FC4E07", "violetred3", "deepskyblue3"), ggtheme = theme_grey())

# Once happy with the number of clusters, find the median value for each cluster:
SDGs_centroids <- SDGs |> 
  mutate(cluster = SDGs_pam$clustering) |> 
  group_by(cluster) |> 
  summarise_at(vars(other_1:SDG3.b_5), median, na.rm = TRUE)
SDGs_centroids

# pam() can also provide the most representative example countries of each cluster. Note: medoids report the standardised data
SDGs_pam$medoids

# We can do a coloured pairwise scatterplot to check data details. We limit it here to the pairs of the first 7 columns because of the large number of possible combinations
pairs(SDGs[, 3:10], col = c("#FC4E07", "violetred3", "deepskyblue3")[SDGs_pam$clustering])

# Questions:

# 1) As the number of clusters grows, so do the groupings of related countries, which become smaller 
# and more specialized. However, we observe some significant overlaps in the clusters. 

# 2) I strongly believe that 4 clusters would be the most ideal as its simmplicity as well as it's 
# accuracy its perfectly balanced. 4 clusters also does not make for untidyness and the data can still
# be viewed with ease. m