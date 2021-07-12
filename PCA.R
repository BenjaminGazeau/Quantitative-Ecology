#Quantitative Biology 
#PCA
#Author: Benjamin Gazeau
#Date: 12 July 2021



# Load Packages 

library(tidyverse)
library(vegan)

# Reading in the Data

env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1)
head(env)


# Doing the PCA: 

# The rda() function does the PCA.
# The PCA usually automatically standardizes data. This will be necessary as we 
# are using the environmental data 
# the PCA preserves Euclidean distances and the
# relationships detected are linear 
env.pca <- rda(env, scale = TRUE)
env.pca

# Eigenvalues for unconstrained axes indicate the relative importance of the resulting reduced axes 
# and can be used to calculate the percentage of total inertia captured by one of the axis.

# To extract the first eigenvalue we must:
round(env.pca$CA$eig[1], 3)
# The total inertia is:
sum(env.pca$CA$eig)
# So the proportion of variation explained by the first PC is:
round(env.pca$CA$eig[1] / sum(env.pca$CA$eig) * 100, 1)
# The result is giving in a percentage 


# QUESTIONS (A) -----------------------------------------------------------

# Why can a PCA not explain all of the variation in a dataset? In other words, why is it best to only use the first few Principal Components? What is explained by the remaining PC axes?

# PCA is a heuristic approach that allows a user to express important properties
# of data along a smaller number of axes than a statistical test.
# As a result, the user only needs to interpret the number of axes required to 
# represent +-75 percent of the data variation (or however many axes represents
# variance that is of interest). This is usually sufficient, and the first two 
# or three Principal Components are frequently used to indicate this.

# Species and Site Scores

summary(env.pca)

# The STRENGTH OF CONTRIBUTION of original environmental factors to the 
# Principal Components is indicated by species scores. Don't let the name 
# deceive you! We're working with environmental data, which the software refers 
# to as'species scores.' The larger (more +) and smaller (more -) numbers, but 
# in opposite directions, indicate greater contribution.

# The positions of the sites are plotted in 2D or 3D ordination space using 
# site scores. In terms of environmental circumstances, sites that are further 
# apart differ significantly.

# The distances between locations are preserved when calling SCALING 1 
# (to interpret species), allowing for a clearer comprehension of how sites 
# relate to one another. The angles between variables are kept when calling 
# SCALING 2 (to interpret variables), resulting in narrower angles between 
# variable vectors reflecting greater correlations.


# Graphical Representations of Ordinations

biplot(env.pca, scaling = 1, main = "PCA Scaling 1", choices = c(1,2))
biplot(env.pca, sclaing = 2, main = "PCA Scaling 2", choices = c(1,2))


source("C:\Users\benji\Desktop\Stats (2021)\Quantitative_Ecology-main\Num_Ecol_R_book_ed1\cleanplot.pca.R")
cleanplot.pca(env.pca, scaling = 1)
cleanplot.pca(env.pca, scaling = 2)

# We can also plot the underlying environmental gradients using odisurf() 


biplot(env.pca, type = c("text", "points"), col = c("black", "black"))
ordisurf(env.pca ~ bod, env, add = TRUE, col = "turquoise", knots = 1)
ordisurf(env.pca ~ alt, env, add = TRUE, col = "salmon", knots = 1)


#QUESTIONS(B) 

#1. Replicate analysis on environmental data included in:

#A) Bird Communities in Yushan Mountain 


ybirds.env <- read.delim('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/ybirds_env.txt', row.names = 1)
ybirds.env <- select(ybirds.env, c(-'Veg.', -'Veg_ext'))
birds.pca <- rda(ybirds.env, scale = TRUE)
summary(birds.pca)
par(mfrow = c(1,2)) 
cleanplot.pca(birds.pca, scaling = 1)
cleanplot.pca(birds.pca, scaling = 2)

#B) Alpine Plant Communities in Aravo 


library (ade4)
data (aravo)
aravo.env <- aravo$env
aravo.env <- select(aravo.env, c(-'Form', -'ZoogD'))
aravo.pca <- rda(aravo.env, scale = TRUE)
summary(aravo.pca)
par(mfrow = c(1,2))
cleanplot.pca(aravo.pca, scaling = 1)
cleanplot.pca(aravo.pca, scaling = 2)


# 2) Discuss the patterns observed:
# {a} Explain the ordination diagram with particular reference to the pattern shown


# The sites in the Bird's PCA biplot exhibit a gradient from left to
# right, essentially in the shape of a 'W'. The highest values of tree density 
# and variety, as well as secondary tree cover, are found in the first cluster 
# of variables. After then, the ecosystem shifts to one with high leaf height 
# diversity, tree circumference, and canopy cover values.
# The final environmental cjhange is from high ground cover to low ground cover and 
# The next large cluster has high elevation and conifers. Overall, 
# there is a natural evolution of the ecology along the sites of:
# Lots of trees with a lot of variety -> Tall trees -> 
# Conifers as the height rises -> High ground cover %