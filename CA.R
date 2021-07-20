#Quantitative Biology 
#Correspondence Analysis 
#Author: Benjamin Gazeau
#Date: 16 July 2021


#Loading the Packages 

library(tidyverse)
library(vegan)
library(ade4)


#Reading in the Data 

spe <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1)
head(spe, 8)
summary(spe)


#Doing the Correspondence Analysis 

spe_ca <- cca(spe)
apply(spe, 1, sum) 
#The '1' ensures the sum function is applied to the rows and not the columns 
#We can see that row 8 sums to 0, and therefore contains 0 species 
spe <- spe[rowSums(spe) > 0,]
spe_ca <- cca(spe)
summary(spe_ca)

#Total inertia
round(sum(spe_ca$CA$eig), 5)

#Inertia of the first axis
round(spe_ca$CA$eig[1], 5)

#Inertia of CA1 and CA2
round(sum(spe_ca$CA$eig[1:2]), 5)

#Fraction of the variance explained by CA1 and CA2
round(sum(spe_ca$CA$eig[1:2]) / sum(spe_ca$CA$eig) * 100, 2)

#The ordination Diagrams
plot(spe_ca, scaling = 1, main = "CA fish abundances - biplot scaling 1") # Scaling 1 emphasizes relationships between ROWS (Sites)
plot(spe_ca, scaling = 2, main = "CA fish abundances - biplot scaling 2") # Scaling 2 emphasizes relationships between COLUMNS (Species)

#Making biplots focusing on 4 specific species:
library(viridis)
library(viridisLite)

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

#a posteriori projection of environmental variables in a CA

env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1)

#Removing the 8th row in spe

env <- dplyr::slice(env, -8)

# the last plot produced (CA scaling 2) must be active
# scaling 2 is default
(spe_ca_env <- envfit(spe_ca, env, scaling = 2))
plot(spe_ca_env, col = "grey40")
plot(spe_ca_env, p.max = 0.05, col = "red") # plot significant variables with a different colour



# Questions:

# 1. How would you explain the patterns seen in the four panels of the above figure?


# Each plot focuses on a certain species and its relationship to the several locales.
# A bubble is plotted for each site, with the size of the bubble corresponding to 
# the species score of the species at that site.
# As a result, the pattern of the bubbles reveals the relative impact a species 
# has on the locations based on its abundance. The contour lines divide sites into 
# groups based on species abundances, allowing us to estimate the gradient of species 
# impact across the sites. The Satr species is abundant in the head and tail ends of 
# the river, but not in the center, according to the Doubs data. Scer is only seen in
# large numbers in a few locations in the center of the river. In the areas along the 
# tail end of the river, Teso and Cogo are abundant.

# 2. Apply these approaches to the [A] Birds and [B] Aravo datasets

#[A] Bird Communities
ybirds.spe <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/ybirds_spe.txt', row.names = 1)

birds_ca <- cca(ybirds.spe)
summary(birds_ca)

par(mfrow = c(1,1)) # View 1 plot on full pane 
plot(spe_ca, scaling = 1, main = "CA Bird Communtites - Biplot Scaling 1")
plot(spe_ca, scaling = 2, main = "CA Bird Communtites - Biplot Scaling 2")

# Posteriori projection of env variables
ybirds.env <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/ybirds_env.txt', row.names = 1)
birds_ca_env <- envfit(birds_ca, ybirds.env, scaling = 2) # Scaling 2 is used and must be active plot
plot(birds_ca_env, p.max = 0.05, col = "grey40") # significant variables are plotted


# [B] Alpine Plant Communities in Aravo 
data(aravo)
aravo.spe <- aravo$spe

aravo_ca <- cca(aravo.spe)
summary(aravo_ca)

par(mfrow = c(1,1)) 
plot(aravo_ca, scaling = 1, main = "CA Alpine Plant Communities - Biplot Scaling 1")
plot(aravo_ca, scaling = 2, main = "CA Alpine Plant Communities - Biplot Scaling 2")

# Posteriori projection of env variables
aravo.env <- aravo$env
aravo_ca_env <- envfit(aravo_ca, aravo.env, scaling = 2) 
plot(aravo_ca_env, p.max = 0.05, col = "grey40") 
