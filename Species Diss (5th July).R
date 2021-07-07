#Quantitative Biology 
#Measures of Biodiversity 
#Author: Benjamin Gazeau
#Date: 05 July 2021

# Quantitative Ecology 2021
# Species dissimilarities 
# 04/07/2021
# Jesse Phillips 


# Loading packages

library(dplyr)
library(vegan)

# Reading in the data 

doubs <- read.csv("C:/Users/benji/Desktop/Stats (2021)/Quantitative_Ecology-main/Num_Ecol_R_book_ed1/DoubsSpe.csv")
doubs <- select(doubs, -X) # remove redundant X col


# QUESTIONS ---------------------------------------------------------------

# 1. Look at the dataset and explain its structure in words

# The collection includes statistics on the abundance of 27 different species at 30 different locations.

# 2. Would we use Bray-Curtis or Jaccard dissimilarities? 
# Bray-Curtis, because we are working with abundance data. 

# 3. Applying the calculation:

BC_doubs <- vegdist(doubs, method = "bray", diag = TRUE)

# 4. Explain the meaning of the data in broad terms

# The Bray-Curtis index measures how different species compositions are at different locations. 
# Every site is compared pairwise, and indices near to 0 indicate that the sites are quite similar,
# while indices close to 1 indicate that the sites are very distinct. 

# 5. Examine it more closely: what general pattern comes out?

# As you progress from site 1 to site 30, the sites grow increasingly distinct from one another. 
# The sites appear to be equally distinct from site 1 from site 19 forward. This is taken to mean 
# that the species composition of site 1 differs significantly from that of the other locations.

# 6. Plotting this pattern 

BC_doubs_plot <- as.data.frame(as.matrix(BC_doubs))

ggplot(BC_doubs_plot, aes(x = 1:30, y = BC_doubs_plot[,1])) +
  geom_line(col = "red4", size = 1) + 
  labs(x = "Sites", y = "Bray-Curtis Dissimilarity Index") +
  theme_bw()

# 7. What explanations can you offer for this pattern?

# Site 1 has considerably different environmental circumstances than site 5, and there appears to 
# be a sharp environmental gradient across the first five locations. There are likely differences 
# in environmental conditions between sites 5 and 20, as well as between those locations.
# The landscape appears to be consistent from site 20 to site 30, although it is vastly different
# from that of site 1. As all those sites are equally dissimilar from site 1. 

# 8. Using the decostand() function, create presence/absence data, and apply the appropriate vegdist() function to obtain a suitable dissimilarity matrix

PA_doubs <- decostand(doubs, method = "pa")

Sor_doubs <- vegdist(PA_doubs, binary = TRUE, diag = TRUE)

# 9. Creating another plot and examining the pattern

Sor_doubs_plot <- as.data.frame(as.matrix(Sor_doubs))

ggplot(Sor_doubs_plot, aes(x = 1:30, y = Sor_doubs_plot[,1])) +
  geom_line(col = "blue3", size = 1) + 
  labs(x = "Sites", y = "Sorenson Dissimilarity Index") + 
  theme_bw()

# The curve's shape is extremely similar to, but not identical to, the curve constructed with abundance data.
# Sites 5 to 15 have less dramatic peaks and troughs, indicating that identical species occur there, although
# their abundances range considerably.