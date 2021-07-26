#Quantitative Biology 
#Measures of Biodiversity 
#Author: Benjamin Gazeau
#Date: 26 July 2021

# 1. Using 2 unconstrained ordination techniques, analyse the mite data
# SECTION A: nMDS -> mite environmental data
# SECTION B: CA -> mite species data

# 2. Using 2 other unconstrained ordination techniques, analyse the dune data
# SECTION C: PCA -> dune species data (tb-PCA)
# SECTION D: PCoA -> dune environmental data


#Loading the packages

library(tidyverse)
library(vegan)
library(cluster)
source("cleanplot.pca.R")


#Reading in the Data 

data("mite")
data("mite.env") 
data("dune") 
data("dune.env")


#Section A: nMDS 

glimpse(mite.env) 

#Calculating distance matrices with Gower's distance

mite.env.matrix <- as.matrix(daisy(mite.env, metric = "gower"))

#Doing the nMDS 

mite.env.nMDS <- metaMDS(mite.env.matrix) 
mite.env.nMDS # Stress = 0.18 (in between fair and suspect)

stressplot(mite.env.nMDS)

#Ploting sites

pl1 <- ordiplot(mite.env.nMDS, type = "none", main = "nMDS Mite Environments")
points(pl1, "sites", pch = 20, col = "black", cex = 1.25)
text(pl1, "sites", pos = 2, cex = 0.75)
abline(h = 0, v = 0, lty = 3)

#envfit() and ploting the environmental variables:

mite.env.fit <- envfit(mite.env.nMDS, mite.env)
plot(mite.env.fit, col = "red", cex = 0.75)

# Findings:

# The ordering relationships between items are represented by nMDS. In terms 
# of environmental characteristics, tight clusters of sites will be quite similar. 
# The environmental variable that has the greatest impact on a site might be used to categorize it.


#Section B: CA 

head(mite)

mite.ca <- cca(mite)
summary(mite.ca) 

# CA Ordination Diagrams
CA_plot_1 <- ordiplot(mite.ca, type = "none", scaling = 1, main = "CA Mite Abundances - Scaling 1")
text(CA_plot_1, "sites", cex = 0.6, col = "navy")
text(CA_plot_1, "species", cex = 0.5, col = "darkred")

CA_plot_2 <- ordiplot(mite.ca, type = "none", scaling = 2, main = "CA Mite Abundances - Scaling 2")
text(CA_plot_2, "sites", cex = 0.5, col = "navy")
text(CA_plot_2, "species", cex = 0.6, col = "darkred")

# FINDINGS 
# On the left, a dense cluster of sites with numbers ranging from 5 to 29. 
# Similar species are most likely to be found at those locations (and relative 
# frequencies of those species). Those sites primarily feature Brachy and ONOV. 
# The large cluster of species in the upper left makes a larger contribution to 
# those locations than the others. The remaining locations are usually grouped 
# together on the right (31-70).
# The species in the bottom right contribute mostly to these sites. 

# Scaling 2 - ordination of species:

# The top left cluster of species will have similar frequencies across all sites.
# The others are a little more evenly distributed.
# The tight cluster of species is linked to the above-mentioned tight cluster of places (5-29)
# These locations are more likely to have these species.

# Conclusion:
# There are certain sites (those found in the 5-29 cluster) which have a higher 
# mite diversity and abundance than the rest.  

#Section C: PCA

head(dune) 

dune_pa <- decostand(dune, method = "pa")
dune_hel <- decostand(dune_pa, method = "hellinger") 
dune.pca <- rda(dune_hel)
summary(dune.pca)

#Creating ordination diagrams:


cleanplot.pca(dune.pca, scaling = 1)
cleanplot.pca(dune.pca, scaling = 2)

# Scaling 1 - the euclidean distances are represented by the distance between the objects 
# There is a linear gradient between the sites 
# Outside the circle of equilibrium, species contribute more 
# to where the sites are located on the plot.
# This represents ten of the thirty species.
# Scaling 2 - Angles between species represent their connections 


#Section D:PCoA 
# ----------------------------------------------------------

glimpse(dune.env)


dune.env.matrix <- as.matrix(daisy(dune.env, metric = "gower")) #calculating the distance matrix with Gower's distance

#PCoa:

dune.env.pcoa <- capscale(dune.env.matrix ~ 1) 
summary(dune.env.pcoa)

# Plot sites
pl2 <- ordiplot(dune.env.pcoa, type = "none", main = "PCoA Dune Environments")
points(pl2, "sites", pch = 20, col = "black", cex = 1.25)
text(pl2, "sites", pos = 2, cex = 0.75)
dune.env.fit <- envfit(dune.env.pcoa, dune.env)
plot(dune.env.fit, col = "red", cex = 0.75)

# PCoA creates a Euclidean representation of object associations based on Gower's index.
# There is no clear linear gradient between the sites.
# Environmental variables aid in the explanation of site correlations.
# The amounts of moisture and land use are plotted as linear gradients.
# The relationship between sites 14 and 15 is best described by soil thickness (A1), as they plot away from the other sites. 
# "Nature Conservation Management" and "0 Manure" plot together with sites 19 and 20,
# which differentiate them from the other sites. 
# Biological farming is concerned with level one manure and moisture, 
# while hobby farming is concerned with level two manure and moisture.
# Standard Farming plots with level 3 Manure and Pasture Use, explaining 
# the relationships between the sites in the plot's bottom left corner.