#Quantitative Biology 
#Unconstrained Ordination: Re-cap
#Author: Benjamin Gazeau
#Date: 08 August 2021

# Load Packages

library(tidyverse)
library(vegan)
library(cluster)
source("cleanplot.pca.R")


# Read in Data 

data("mite") # mite species abundance data
data("mite.env") # mite environmental data 
data("dune") # dune species abundance data
data("dune.env") # dune environmental data


# SECTION A: nMDS 

glimpse(mite.env) # Variables are not homogeneous

# Calculate distance matrix with Gower's distance
mite.env.matrix <- as.matrix(daisy(mite.env, metric = "gower"))

# Do nMDS 
mite.env.nMDS <- metaMDS(mite.env.matrix) 
mite.env.nMDS # Stress = 0.18 (in between fair and suspect)

stressplot(mite.env.nMDS)

# Plot sites
pl1 <- ordiplot(mite.env.nMDS, type = "none", main = "nMDS Mite Environments")
points(pl1, "sites", pch = 20, col = "black", cex = 1.25)
text(pl1, "sites", pos = 2, cex = 0.75)
abline(h = 0, v = 0, lty = 3)

# Do envfit() and plot environmental variables
mite.env.fit <- envfit(mite.env.nMDS, mite.env)
plot(mite.env.fit, col = "red", cex = 0.75)

# SECTION B: CA

head(mite) # Abundance data, appropriate for CA

mite.ca <- cca(mite)
summary(mite.ca) # Total Inertia = 1.696, 75% explained by first 7 components 

# CA Ordination Diagrams
CA_plot_1 <- ordiplot(mite.ca, type = "none", scaling = 1, main = "CA Mite Abundances - Scaling 1")
text(CA_plot_1, "sites", cex = 0.6, col = "navy")
text(CA_plot_1, "species", cex = 0.5, col = "darkred")

CA_plot_2 <- ordiplot(mite.ca, type = "none", scaling = 2, main = "CA Mite Abundances - Scaling 2")
text(CA_plot_2, "sites", cex = 0.5, col = "navy")
text(CA_plot_2, "species", cex = 0.6, col = "darkred")


# SECTION C: PCA 

head(dune) # Abundance data

dune_pa <- decostand(dune, method = "pa") # convert to presence-absence data
dune_hel <- decostand(dune_pa, method = "hellinger") # transform data to Hellinger distance (immune to double-zero problem)

# Do PCA
dune.pca <- rda(dune_hel)
summary(dune.pca)

# Use cleanplot.pca() to create ordination diagrams
cleanplot.pca(dune.pca, scaling = 1)
cleanplot.pca(dune.pca, scaling = 2)

# SECTION D: PCoA 

glimpse(dune.env) 

# calculate distance matrix with Gower's distance
dune.env.matrix <- as.matrix(daisy(dune.env, metric = "gower"))

# Do the PCoa
dune.env.pcoa <- capscale(dune.env.matrix ~ 1) 
summary(dune.env.pcoa)

# Plot sites
pl2 <- ordiplot(dune.env.pcoa, type = "none", main = "PCoA Dune Environments")
points(pl2, "sites", pch = 20, col = "black", cex = 1.25)
text(pl2, "sites", pos = 2, cex = 0.75)

# Do envfit() and plot environmental variables
dune.env.fit <- envfit(dune.env.pcoa, dune.env)
plot(dune.env.fit, col = "red", cex = 0.75)
