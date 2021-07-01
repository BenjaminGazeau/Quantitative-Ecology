#Quamtitative Biology 
#Measures of Biodiversity 
#Author: Benjamin Gazeau
#Date: 29 June 2021

#Load the packages 

library(vegan)
library(ggplot2)

#Reading in the seaweed data:
#to read in csv files on mac you click on the "session tab",
#choose "set working directory" and choose file,
#in the console there will be the directory/pathway that can be copied 

spp <-  read.csv("~/Desktop/Quantitative_Ecology-main/exercises/diversity/SeaweedsSpp.csv")

#Understand the data 
#how many rows and columns? use dim function 
dim(spp)

spp[1:10, 1:10]

spp <- dplyr::select(spp, -1)

#Looking at species richness 

library(tcltk)
library(BiodiversityR)

spp_richness <- diversityresult(spp, index = 'richness', method = 'each site')

specnumber(spp, MARGIN= 1)

ggplot(data= spp_richness, (aes(x = 1:58, y= richness)))+
  geom_line(col="red")+ xlab('Coastal selection, west to east')+
  ylab("species Richness")

#univariate diversity indices
Light <- read.csv("~/Desktop/Quantitative_Ecology-main/exercises/diversity/light_levels.csv")

#preview "light" dataset
Light

#Calculating the Jaccard dissimilarity 

Light_div <- data.frame(
  site = c("Low_light", "mid_light", "high_light"),
  richness = specnumber(Light[,2:7], MARGIN = 1),
  shannon = round(diversity(Light[,2:7], MARGIN = 1, index = "shannon"), 2),
  simpson = round(diversity(Light[,2:7], MARGIN = 1, index = "shannon"), 2)
)

#Preview "Light_div
Light_div

# Dissimilarity indices 
# Back to the seaweed data 
# Presence/absence data, therefore use sorenson dissimilarity 
# by setting binary = 2

sor <- round(vegdist(spp, binary = TRUE, diag = TRUE), 4)

#What are the class and dimensions?

class(sor)
dim(sor)  

sor_df <- as.data.frame(as.matrix(sor)) 
class(sor_df)  
dim(sor_df)  
sor_df[1:10, 1:10]  

#can write to csv: write.csv(sor,df, file = "~/Desktop/Quantitative_Ecology-main/exercises/diversity/SeaweedSpp_dis_matrix.csv")

#Gamma diversity 

#the number of columns gives the total number of species in this example 

ncol(spp)

#we may also use:

diversityresult(spp, index = 'richness', method = 'pooled')
#why this difference?

#Beta Diversity 

#Whittaker's concept of Beta Diversity 
#true beta 

true_beta <- data.frame(
  beta =  specnumber(spp, MARGIN = 1) / ncol(spp),
  section_no  = c(1:58)
)
ggplot(data = true_beta, (aes(x= section_no, y = beta))) +
  geom_line() + xlab("Coastal section, west to east")+ ylab("True beta-diversity")

#absolute species turnover 

abs_beta <- data.frame(
  beta = ncol(spp)- specnumber(spp, MARGIN = 1),
  section_no = c(1:58)
)

ggplot(data = abs_beta, (aes(x = section_no, y = beta)))+
  geom_line()+ xlab("Coastal section, west to east")+ ylab("Absolute beta-diversity")

# Contemporary difinitons beta-diversity

library(betapart)

#decompose total Sorensen Dissimilarity into turnover and nestedness-resultant components 

#Compute the basic quantities needed for computing  the multiple-site beta diversity 
#measures and pairwise dissimilarity matrices

Y.core <- betapart.core(spp)
str(Y.core)

# Compute 3 distances matrices accounting for the 
# (i) turnover (replacement)
# (ii) nestedness-resultant component, and
# (iii) total dissimilarity (i.e. the sum of both components)
Y.pair <-  beta.pair(Y.core, index.family = "sor")
str(Y.pair)

# let Y1 be the turnover component 9beta-sim):
Y1 <- as.matrix(Y.pair$beta.sim)

# let y2 be the nestedness-resultant component (beta-sne):
Y2 <- as.matrix(Y.pair$beta.sne)

#Questions:

# 1) Why is the matrix square, and what determines the number of rows/columns? 

# the dissimilarity values of every site is compared to the values of every other site
# (including itself). Therefore, the 58 rows are compared to the other 58 columns. This 
# results in all of the values being squared. 

# 2) What is the meaning of the diagonal? 

# The diagonal values will always equal to zero due to the fact that each site is being 
# compared to itself. This means that these sites have a dissimilarity of zero.

# 3) What is the meaning of the non-diagonal elements?

# The non-diagonal elements represent the dissimilarity values which are a result of one 
# site being compared to another. The dissimilarity values fall between 0-1. The closer
# the value to 1 the more dissimilar the sites are to each other. 

# 4) Take the data in row 1 and create a line graph that shows these values as a function 
#    of section number. 

ggplot(data = sor_df, (aes(x= 1:58, y = sor_df[,1])))+
  geom_line()+ xlab("Site")+ ylab("Site")


# 5) Provide a mechanistic (ecological) explanation for why this figure takes the shape that it does.

# The graph starts at zero as it is the value of site one being compared with itself. The dissimilarity
# values increase as you move away from site one due to the fact that the sites further away from site 1 will
# tend to be more dissimilar. 
#

# Questions for section 2:

# 1) There is a difference because the function ncol physically counts all of the columns
# including the first column which is usually contains the heading. The ncol function returns 
# the number of columns of a matrix or data frame.The diversityresult 
# function only calculates the diversity using the columns with values.

# 2) The diversityresult function is the correct method.

