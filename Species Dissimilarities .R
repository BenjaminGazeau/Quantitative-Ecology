#Quantitative Biology 
#Measures of Biodiversity 
#Author: Benjamin Gazeau
#Date: 02 July 2021

# Load Packages

library(vegan)
library(ggplot2)
library(BiodiversityR)
library(betapart)
library(dplyr)
# Read Data 

spp <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/exercises/diversity/SeaweedsSpp.csv")
spp <- select(spp, -1) #remove X column


# Whittaker's Beta-diversity 
# 1) True Beta-diversity
# Gamma-diversity / Alpha-diversity
true_beta <- data.frame(
  beta = ncol(spp) / specnumber(spp, MARGIN = 1),
  section_no = c(1:58)
)
# Plot true beta
ggplot(data = true_beta, (aes(x = section_no, y = beta))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("True beta-diversity")


# 2) Absolute Species Turnover
# Gamma-diversity - Alpha-diversity 
abs_beta <- data.frame(
  beta = ncol(spp) - specnumber(spp, MARGIN = 1),
  section_no = c(1:58)
)
# Plot absolute species turnover
ggplot(data = abs_beta, (aes(x = section_no, y = beta))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("Absolute beta-diversity")



# Contemporary definitions of Beta-diversity ------------------------------

# 1. Species Turnover (beta-sim)
# Same alpha-diversity between sections, but different species makeup 
# Hence, this beta-diversity refers to processes that cause communities to differ
# due to species being lost/gained from section to section without corresponding  
# changes in alpha-diversity 

# 2. Nestedness-Resultant Beta-Diversity (beta-sne)
# Two places share same species, number of species can differ amongst quadrants 
# Nestedness-resultant beta-diversity refers to processes that cause species to be 
# gained or lost, and the community with the lowest alpha-diversity is a subset of
# the richer community


# Calculating Turnover and Nestedness-Resultant  --------------------------

# Decompose total Sørensen dissimilarity into turnover and nestedness-resultant components:
Y.core <- betapart.core(spp)
Y.pair <- beta.pair(Y.core, index.family = "sor")

# Let Y1 be the turnover component (beta-sim):
Y1 <- as.matrix(Y.pair$beta.sim)

# Let Y2 be the nestedness-resultant component (beta-sne):
Y2 <- as.matrix(Y.pair$beta.sne)

round(Y1[1:10, 1:10], 4) # Matrix of beta-sim
round(Y2[1:10, 1:10], 4) # Matrix of beta-sne


# QUESTIONS ---------------------------------------------------------------

# 1. Plot species turnover as a function of Section number, and provide a mechanistic exaplanation for the pattern observed.
y1_plot <- as.data.frame(Y1)
ggplot(y1_plot, aes(x = 1:58, y = y1_plot[,1])) + 
  geom_line() +
  labs(x = "Coastal Section, West to East", y = "Species Turnover (beta-sim)")

# The graph shows an increasing trend as we move from 1 to 58, meaning there are
# different species compositions as we move from site 1 on the west coast to site
# 58 on the east coast. The species composition along the shore is constantly changing, 
# and the species composition of site 1 is completely different from that of site 58. 
# As we calculated beta-sim, while the makeup of species changes, the species richness does not.


# 2. Based on an assessment of literature on the topic, provide a discussion of nestedness-resultant
# β-diversity. Use either a marine or terrestrial example to explain this mode of structuring biodiversity.
# Biotas from sites with lower species richness are subsets of biotas from sites with higher species richness,
# resulting in nestedness of species assemblages. This represents the non-random process of species extinction
# at different locations (Baselga 2010). In the simplest terms, this is demonstrated by the presence of 12 
# species at Site A, 4 species (of the original 12) at Site B, and just 2 species at Site C. There is no 
# turnover in species, but a difference in richness, where the composition of Site C is a subset of Site A. 