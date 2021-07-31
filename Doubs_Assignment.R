# Load required packages

library(vegan)
library(tidyverse)
library(corrplot)
library(ggcorrplot)


spe <- as.tibble(read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1))
spe <- spe[-8,]
spe

# calculate the species dissimilarities based on the Jaccard Index

spe_dist <- round(vegdist(spe, method = "jaccard", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist))

# Environmental data frame
env <- as.tibble(read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1))
env <- env[-8,]
env

(env_std <- as.tibble(round(decostand(env, method = "standardize"), 2)))
env_dist_all <- round(vegdist(env_std, method = "euclidian", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist_all))

# we can also calculate distances based on only one (or a few) variables
env_dist_alt <- round(vegdist(env_std[, "alt"], method = "euclidian", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist_alt))

# PCA ---------------------------------------------------------------------

library(tidyverse)
library(vegan)


# The Doubs River data ----------------------------------------------------

env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)
head(env)


# Do the PCA --------------------------------------------------------------

env_pca <- rda(env, scale = TRUE)
str(env_pca)
env_pca

env_pca$CA$eig

round(env_pca$CA$eig[1], 3)

sum(env_pca$CA$eig)

round(env_pca$CA$eig[1] / sum(env_pca$CA$eig) * 100, 1) # result in %

summary(env_pca)


# Graphical represenations of ordinations 

par(mfrow = c(1, 2))
biplot(env_pca, scaling = 1, main = "PCA scaling 1", choices = c(1, 2))
biplot(env_pca, scaling = 2, main = "PCA scaling 2", choices = c(1, 2))


source("cleanplot.pca.R")
cleanplot.pca(env_pca, scaling = 1)
cleanplot.pca(env_pca, scaling = 2)


biplot(env_pca, type = c("text", "points"), col = c("black", "black"))
ordisurf(env_pca ~ bod, env, add = TRUE, col = "turquoise", knots = 1)
ordisurf(env_pca ~ alt, env, add = TRUE, col = "salmon", knots = 1)

# Doing the CA ------------------------------------------------------

# The Doubs River species data --------------------------------------------

spe <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpe.csv")
spe <- dplyr::select(spe, -1)
head(spe, 8)


# Do the CA ---------------------------------------------------------------

spe_ca <- cca(spe)
apply(spe, 1, sum)
spe <- spe[rowSums(spe) > 0, ]
head(spe, 8)


spe_ca <- cca(spe)
spe_ca

summary(spe_ca)

# sum of all eigenvalues
# is the total inertia

round(sum(spe_ca$CA$eig), 5)

# the first eigenvalue
round(spe_ca$CA$eig[1], 5)

# sum of CA1 and CA2
round(sum(spe_ca$CA$eig[1:2]), 5)

# proportion explained by CA1 and CA2
round(sum(spe_ca$CA$eig[1:2]) / sum(spe_ca$CA$eig) * 100, 2)


# Ordination diagrams 

par(mfrow = c(1, 2))
plot(spe_ca, scaling = 1, main = "CA fish abundances \n - biplot scaling 1")
plot(spe_ca, scaling = 2, main = "CA fish abundances \n - biplot scaling 2")

# Scaling 1: This scaling emphasises relationships between rows accurately in low-dimensional ordination space. Distances among objects (samples or sites) in the biplot are approximations of their 𝜒2  distances in multidimensional space. Objects found near a point representing a species are likely to contain a high contribution of that species.

# Scaling 2: This scaling emphasises relationships between columns accurately in low-dimensional ordination space. Distances among objects (samples or sites) in the biplot are not approximations of their 𝜒2  distances in multidimensional space, but the distances among species are. Species positioned close to the point representing an object (a sample or site) are more likely to be found in that object or to have higher frequency there.

require('viridis')
palette(viridis(8))
par(mar = c(4, 4, 0.9, 0.5) + .1, mfrow = c(2, 2))
with(spe, tmp <- ordisurf(spe_ca ~ Satr, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Satr"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_ca ~ Scer, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Scer"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_ca ~ Teso, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Teso"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_ca ~ Cogo, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Cogo"))
abline(h = 0, v = 0, lty = 3)

# a posteriori projection of environmental variables in a CA
env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)


env <- dplyr::slice(env, -8)

(spe_ca_env <- envfit(spe_ca, env, scaling = 2))
plot(spe_ca_env, col = "grey40")
plot(spe_ca_env, p.max = 0.05, col = "red") 


# Doing the NMDS plots ------------------------------------------------------

spe <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpe.csv")
spe <- dplyr::select(spe, -1)
spe <- dplyr::slice(spe, -8)


spe_nmds <- metaMDS(spe, distance = "bray")
spe_nmds

summary(spe_nmds) 


# Ordination diagrams 

par(mfrow = c(2, 2))
stressplot(spe_nmds, main = "Shepard plot")
ordiplot(spe_nmds, type = "t", cex = 1.5, main = paste0("nMDS stress = ", round(spe_nmds$stress, 2)))
gof = goodness(spe_nmds)
plot(spe_nmds, type = "t", main = "Goodness of fit")
points(spe_nmds, display = "sites", cex = gof * 200) # bigger bubbles indicate a worse fit

pl <- ordiplot(spe_nmds, type = "none", main = "nMDS fish abundances ")
points(pl, "sites", pch = 21, cex = 1.75, col = "grey80", bg = "grey80")
points(pl, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl, "species", col = "blue4", cex = 0.9)
text(pl, "sites", col = "red4", cex = 0.9)

env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)
env <- dplyr::slice(env, -8)

(spe_nmds_env <- envfit(spe_nmds, env))
plot(spe_nmds_env, col = "grey40")
plot(spe_nmds_env, p.max = 0.05, col = "red")

require('viridis')
palette(viridis(8))
par(mar = c(4, 4, 0.9, 0.5) + .1, mfrow = c(2, 2))
with(spe, tmp <- ordisurf(spe_nmds ~ Satr, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Satr"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_nmds ~ Scer, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Scer"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_nmds ~ Teso, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Teso"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_nmds ~ Cogo, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Cogo"))
abline(h = 0, v = 0, lty = 3)

env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)
env <- dplyr::slice(env, -8)

(spe_nmds_env <- envfit(spe_nmds, env))
plot(spe_nmds_env, col = "grey40")
plot(spe_nmds_env, p.max = 0.05, col = "red")

sites <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpa.csv", row.names = 1)

spe <- slice(spe, -8) # No observations in row 8, therefore it is removed
env <- slice(env, -8) 
sites <- slice(sites, -8) 

# Env Variables (Correlations) 

corr <- cor(env) # matrix of correlations amongst env variables
pmat <- cor_pmat(env) # matrix of p-values for correlations

# Plot of pairwise correlations
corpl <- ggcorrplot(corr,
                    hc.order = TRUE, # orders plot by hierarchical clustering
                    type = "lower", # lower triangle only
                    lab = TRUE, # displays values
                    p.mat = pmat, # only shows significant correlations
                    insig = "blank", # blocks out the rest
                    legend.title = "") 
corpl
