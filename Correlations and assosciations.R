#Quantitative Biology 
#Correlations and assosciations
#Author: Benjamin Gazeau
#Date: 06 July 2021



#Load Packages 

library(tidyverse)
library(vegan)
library(Hmisc) 
library(ggcorrplot)


#Reading in the Data 

env <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1) # row.names ensures that random X col doesn't appear
head(env)


#Correlations

#Environmental factors are related to one another across sample sites via correlations. 
#There is no need to standardize, but data transformations may be necessary in some circumstances.

round(cor(env), 2)

rcorr(as.matrix(env))

#Questions(A) 

# 1. Creating a plot of pairwise correlations


corr <- cor(env) 
p.mat <- cor_pmat(env) # computes correlation matrix of p-values

ggcorrplot(corr,  
           hc.order = TRUE,  #reorders matrix using hierarchical clustering 
           type = "lower",  #displays lower triangle only 
           lab = TRUE,  
           p.mat = p.mat,  #bars non-significant correlations
           insig = "blank") 

# 2) Name two top positive and two top negative statistically-significant correlations

# POSITIVE: pho-amm and dfs-flo & NEGATIVE: alt-dfs and alt-flo 

# 3) For each, discuss the mechanism behind the relationships. Why do these relationships exist?

# Phosphate concentration and Ammonium concentration (+):
# Distance from source and Mean minimum discharge (+):
# Altitude and Distance from source (-): 
# The lower the altitude, the further away from the source. The source of  the river is located at 
# a high altitude, which flows downwards. 
# Altitude and Mean minimum discharge (-): 


#Associations

# Doubs River fish species dataset contains abundance data

spp <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1)
head(spp)

#In order to calculate an association matrix, we must first transpose the data

spp_t <- t(spp)



#Questions(B)

# 1. Why do we need to transpose the data? 
# The species are required to be in rows in order to calculate an associations matrix


# 2. What are the properties of the transposed species table?

# The transposed table has the same features as the original matrix, but the columns and rows are reversed.
# As a result, the species are regarded as observations, whereas the locations are regarded as variables.


#Association Matrix

spp_assoc1 <- vegdist(spp_t, method = "jaccard")
as.matrix((spp_assoc1))[1:10, 1:10]

spp_assoc2 <- vegdist(spp_t, method = "jaccard", binary = TRUE)
as.matrix((spp_assoc2))[1:10, 1:10]


#Questions(C)

# 1. What are the properties of the association matrix? How do these properties differ from a i) species dissimilarity matrix and ii) correlation matrix?

# Jaccard's coefficients are the attributes of the association matrix. A similarity matrix measures the 
# similarity between species (or the percentage overlap at each site), whereas a dissimilarity matrix does 
# the converse, and a correlation matrix calculates the dependency of one variable on another.

# 2. What is the difference between spp_assoc 1 and 2? Is the information markedly different?

# spp assoc1 uses abundance data to determine association, while spp assoc2 employs presence-absence data. 
# In the case of presence-absence data, the information isn't significantly different; it's just less 
# severe (variable).

# 3. Explain the kind of insight we can glean from a species association matrix.

# We may deduce which species appear more frequently at various sites using the Jaccard similarity, 
# as a higher coefficient indicates a greater amount of overlap.