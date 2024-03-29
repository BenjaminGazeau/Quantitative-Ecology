# Packages

library(tidyverse)
library(vegan)
library(missMDA) 
library(ggcorrplot) 


# SDG 1.a Domestic general government health expenditure
SDG1.a <- read.csv("WHO_SDG1.a_domestic_health_expenditure.csv")%>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG1.a")

# SDG 3.1 Maternal mortality ratio
SDG3.1_1 <- read.csv("WHO_SDG3.1_maternal_mort.csv") %>%
  filter(Period == 2016,
         Indicator == "Maternal mortality ratio (per 100 000 live births)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.1_1")

# SDG 3.1 Births attended by skilled health personnel 
SDG3.1_2 <- read.csv("WHO_SDG3.1_skilled_births.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.1_2")

# SDG 3.2 Number of neonatal deaths
SDG3.2_1 <- read.csv("WHO_SDG3.2_neonatal_deaths.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.2_1")

# SDG 3.2 Number of under-five deaths
SDG3.2_2 <- read.csv("WHO_SDG3.2_under_5_deaths.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.2_2")

# SDG 3.2 Number of infant deaths
SDG3.2_3 <- read.csv("WHO_SDG3.2_infant_deaths.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.2_3")

# SDG 3.3 New HIV infections 
SDG3.3_1 <- read.csv("WHO_SDG3.3_new_HIV_infections.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_1")

# SDG 3.3 Incidence of tuberculosis
SDG3.3_2 <- read.csv("WHO_SDG3.3_TB.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_2")

# SDG 3.3 Malaria incidence
SDG3.3_3 <- read.csv("WHO_SDG3.3_malaria.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_3")

# SDG 3.3 Hepatitis B surface antigen prevalence among children under 5
SDG3.3_4 <- read.csv("WHO_SDG3.3_hepatitis_B.csv") %>%
  filter(Period == 2015) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_4")

# SDG 3.3 Reported number of people requiring interventions against NTDs
SDG3.3_5 <- read.csv("WHO_SDG3.3_NCD_interventions.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_5")

# SDG 3.4 Adult mortality rate 
SDG3.4_1 <- read.csv("WHO_SDG3.4_adult_death_prob.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_1")

# SDG 3.4 Number of deaths attributed to non-communicable diseases, by type of disease and sex
SDG3.4_2 <- read.csv("WHO_SDG3.4_NCD_by_cause.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes",
         Dim2 == "Diabetes mellitus") %>%
  mutate(Indicator = Dim2) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_2")

SDG3.4_3 <- read.csv("WHO_SDG3.4_NCD_by_cause.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes",
         Dim2 == "Cardiovascular diseases") %>%
  mutate(Indicator = Dim2) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_3")

SDG3.4_4 <- read.csv("WHO_SDG3.4_NCD_by_cause.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes",
         Dim2 == "Respiratory diseases") %>%
  mutate(Indicator = Dim2) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_4")

# SDG 3.4 Crude suicide rates
SDG3.4_5 <- read.csv("WHO_SDG3.4_suicides.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_5")

# SDG 3.4 Total NCD deaths
SDG3.4_6 <- read.csv("WHO_SDG3.4_NCD_data_total.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_6")

# SDG 3.5 Alcohol, total per capita consumption
SDG3.5 <- read.csv("WHO_SDG3.5_alcohol_consumption.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.5")

# SDG 3.6 Estimated road traffic death rate
SDG3.6 <- read.csv("WHO_SDG3.6_traffic_deaths_prop.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.6")

# SDG 3.7 Adolescent birth rate
SDG3.7 <- read.csv("WHO_SDG3.7_adolescent_births.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.7")

# SDG 3.8 UHC Index of service coverage 
SDG3.8_1 <- read.csv("WHO_SDG3.8_UHC_data_availability.csv") %>%
  filter(Period == "2013-2017") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.8_1")

# SDG 3.8 Data availability for UHC index of essential service coverage 
SDG3.8_2 <- read.csv("WHO_SDG3.8_UHC_index_of_service_coverage.csv") %>%
  filter(Period == 2017) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.8_2")

# SDG 3.9 Poison control and unintentional poisoning 
SDG3.9_1 <- read.csv("WHO_SDG3.9_unintentional_poisoning_prop.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.9_1")

# SDG 3.9 Mortality attributed to exposure of unsafe WASH services
SDG3.9_3 <- read.csv("WHO_SDG3.9_WASH_mortalities.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.9_3")

# SDG 16.1 Estimates of rate of homicides
SDG16.1 <- read.csv("WHO_SDG16.1_homicides.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG16.1")

# SDG 3.a Prevalence of current tobacco use 
SDG3.a <- read.csv("WHO_SDG3.a_tobacco_control.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.a")

# SDG 3.b Total net official development assistance to medical research and basic health
SDG3.b_1 <- read.csv("WHO_SDG3.b_dev_assistence_for_med_research.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_1")

# SDG 3.b Measles-containing-vaccine second-dose immunization coverage 
SDG3.b_2 <- read.csv("WHO_SDG3.b_measles_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_2")

# SDG 3.b Diphtheria tetanus toxoid and pertussis immunization coverage 
SDG3.b_3 <- read.csv("WHO_SDG3.b_diphtheria_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_3")

# SDG 3.b Pneumococcal conjugate vaccines immunization coverage
SDG3.b_4 <- read.csv("WHO_SDG3.b_pneumococcal_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_4")

# SDG 3.b Girls aged 15 years old that received the recommended doses of HPV vaccine
SDG3.b_5 <- read.csv("WHO_SDG3.b_HPV_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_5")

# SDG 3.c Health Workforce
SDG3.c_1 <- read.csv("WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Medical doctors (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_1")

SDG3.c_2 <- read.csv("WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Nursing and midwifery personnel (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_2")

SDG3.c_3 <- read.csv("WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Dentists (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_3")

SDG3.c_4 <- read.csv("WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Pharmacists  (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_4")

# SDG 3.d Average of 13 International Health Regulations core capacity scores
SDG3.d_1 <- read.csv("WHO_SDG3.d_health_risks.csv")  %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.d_1")

# Other Life expectancy at birth
other_1 <- read.csv("WHO_Other_life_expectancy.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes",
         Indicator == "Life expectancy at birth (years)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "other_1")

# Other Life expectancy at age 60 
other_2 <- read.csv("WHO_Other_life_expectancy.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes",
         Indicator == "Life expectancy at age 60 (years)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "other_2")


# Discussion 

# Each of the applicable indicators for SDG 3 - "Good Health and Well-Being" - 
# is read in and allocated to its own object in the R environment. Only the 
# time period we wish to use (primarily 2016) is selected, and only the variables 
# that are required are selected, 
# those being the description of the indicator and a numerical measurement of 
# indicator, as well as the measuring site The goal of reading in the data 
# one at a time like this is to make sure we're only working with the 
# time periods we need for each indicator. We have complete control on the data 
# we're entering, which allows us to be more confident in our PCA.


# rbind the Data 

health <- do.call("rbind", lapply(ls(),get))
head(health)

# Discussion 

# To merge all of the items currently in the R 
# environment, use the rbind() method. Each data frame is stacked on top of the one before it, 
# resulting in a full, long-format data frame with all of the data.


# Create List of SDGs Used 

unique(health[, c(5, 1)])
write_csv(unique(health[, c(5, 1)]), file = "SDG_description.csv")

# Discussion 

# Each explanation of the indicators is extracted using the unique() function, 
# and this is written to a new csv. This is done so that we can keep track of 
# what each indicator indicates on a reference page.

# Pivot Wider 

health_wide <- health %>%
  arrange(Location) %>%
  select(-Indicator) %>%
  pivot_wider(names_from = SDG, values_from = FactValueNumeric) %>%
  as_tibble()

# Discussion

# The full data frame in long format is converted to wide format. This is done to 
# ensure that each country appears only once as an observation and that each SDG is a variable.


# Add World Population Data 

popl <- read_csv("WHO_population.csv") %>%
  filter(Year == 2016) %>%
  rename(popl_size = `Population (in thousands) total`,
         Location = Country) %>%
  select(Location, popl_size) %>%
  mutate(popl_size = as.numeric(gsub("[[:space:]]", "", popl_size)) * 1000)

health_wide <- health_wide %>%
  left_join(popl)


# Express Some Variables to Unit of Population Size 

health_wide <- health_wide %>%
  mutate(SDG3.4_4 = SDG3.4_4 / popl_size * 100000,
         SDG3.4_3 = SDG3.4_3 / popl_size * 100000,
         SDG3.4_2 = SDG3.4_2 / popl_size * 100000,
         SDG3.4_6 = SDG3.4_6 / 100,
         SDG3.2_2 = SDG3.2_2 / popl_size * 100000,
         SDG3.2_3 = SDG3.2_3 / popl_size * 100000,
         SDG3.2_1 = SDG3.2_1 / popl_size * 100000)

# Discussion

# Variables formerly expressed as values (per 100,000 population) have been converted to 
# be directly comparable to the population variable from the previous phase.
# We want all of our variables to be directly comparable before we run our 
# PCA so that we can test for collinearity.


# Histograms of Missing Values, and Correlations 
# calculate histograms

health_wide$na_count <- apply(health_wide[, 3:(ncol(health_wide) - 1)], 1, function(x) sum(is.na(x)))
hist(health_wide$na_count, breaks = 14, plot = TRUE)

# remove rows where there are more than 10 NAs
health_wide <- health_wide %>%
  filter(na_count <= 10) %>%
  select(-na_count)

# calculate pairwise correlations
corr <- round(cor(health_wide[, 3:(ncol(health_wide) - 1)]), 1)

# visualization of the correlation matrix
ggcorrplot(corr, type = 'upper', outline.col = "grey60",
           colors = c("#1679a1", "white", "#f8766d"),
           lab = TRUE)

# Discussion 

# We make a histogram that shows how many NAs we have in our dataframe, and then we eliminate rows 
# (countries) where we have more than 10 NAs, indicating that we are lacking more than 10 indications.
# The pairwise correlations between all of the indicators are then calculated and displayed 
# as a correlation plot. This is done to see which indicators have a high correlation with one another.


# Impute Remaining NAs 

health_wide_complete <- imputePCA(health_wide[, 3:(ncol(health_wide) - 1)])$completeObs

# save for later use
SGD_data <- cbind(health_wide[, 1:2], health_wide_complete)
write_csv(SGD_data, file = "data/WHO/SDG_complete.csv")

# We impute (assign) the NA values with predicted values (such as the mean of the   
# variable). This is done before the PCA so that we don't have NAs in the data frame.
# We return only the imputed dataset from the imputePCA() function.


# Scale and Center the Data and do the PCA --------------------------------

health_wide_complete_std <- decostand(health_wide_complete, method = "standardize")
health_pca <- rda(health_wide_complete_std)
health_pca
summary(health_pca)

# We standardize our data using the decostand() function and perform the PCA using
# the rda() function. The total inertia of the analysis is 37 and the first two 
# principal components explain about 55% of the inertia. 


# Graphical Displays ------------------------------------------------------

# vegan Biplots
biplot(health_pca, scaling = 1, main = "PCA scaling 1", choices = c(1, 2))
biplot(health_pca, scaling = 2, main = "PCA scaling 2", choices = c(1, 2))

# Assembled using vegan component functions
pl1 <- ordiplot(health_pca, type = "none", scaling = 1, main = "PCA WHO/SDG")
points(pl1, "sites", pch = 21, cex = 1.0, col = "grey20", bg = "grey80")
points(pl1, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl1, "species", col = "blue4", cex = 0.9)
text(pl1, "sites", col = "red4", cex = 0.9)

pl2 <- ordiplot(health_pca, type = "none", scaling = 2, main = "PCA WHO/SDG")
points(pl2, "sites", pch = 21, cex = 1.75, col = "grey80", bg = "grey80")
points(pl2, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl2, "species", col = "blue4", cex = 0.9)
text(pl2, "sites", col = "red4", cex = 0.9)

# Using ggplot()
site_scores <- tibble(ParentLocation = health_wide$ParentLocation,
                      Location = health_wide$Location)
site_scores <- tibble(cbind(site_scores, scores(health_pca, display = "sites", choices = c(1:7))))
species_scores <- data.frame(scores(health_pca, display = "species", choices = c(1:7)))
species_scores$species <- rownames(species_scores)
species_scores <- tibble(species_scores)

ggplot(data = site_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = ParentLocation)) +
  geom_segment(data = species_scores, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.4, "cm"), type = "closed"), 
               color = "lightseagreen", alpha = 1, size = 0.3) +
  geom_text(data = species_scores, 
            aes(x = PC1, y = PC2, label = species),
            color = "black") +
  xlab("PC1") + ylab("PC2") + 
  ggtitle("WHO SDGs, Scaling 2")



# Questions:

# 1) THE CODE HAS BEEN DISCUSSED AFTER EACH SECTION.

# 2. Discuss and explain the patterns observed. How does South Africa fare in 
# terms of attaining SDGs? Contrast with some key countries of your choice to 
# make your points. Label the key countries that you refer to in your text by 
# updating the code accordingly.

# Patterns:

# There is a clear distinct separation of African countries.
# African countries show high incidences of infectious diseases such a malaria and Hepatitis B. 
# There is also a large percentage of indivuals living in African countires that suffer from HIV.
# These medical issues are also exacerbated by poor medical facilities and a lack of medical staff.
# European countries show the opposite trends. Most of the European countires are third world countries
# and people living in these countries tend to have a higher standard of living. 
# Government health spending is also higher in European and Asian countries.
# Commonalities between Africa and Europe: diabetes prevalence and suicide rate; the 
# Americas have the highest rates of both.In terms of SDG 3, "excellent health and well-being," 
# African countries appear to be doing the poorest. This is due to a lack of access to quality healthcare.
# The Americas, as well as countries in the Western Pacific and Eastern Mediterranean, 
# are in the middle of the range, with high suicide rates despite being completely vaccinated.

# South Africa:

# South Africa is the best country in the whole of Africa with regards to achieving our SDG's. 
# SA falls into the middle category in our PCA. Our healthcare sector also struggles as there is a 
# deficit of healhcare workers aswell as a deficit in healthcare facilities. Our country struggles with
# diabetes and HIV and we also have high suicide rates. Crime rates in SA are also at an alltime high
# which exacerbates many issues we struggle with in out country. 
# unfortuntaely 

