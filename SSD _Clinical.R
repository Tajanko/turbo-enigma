# Load necessary libraries
if(!require(readxl)) install.packages("readxl", dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require(tidyr)) install.packages("tidyr", dependencies=TRUE)
if(!require(broom)) install.packages("broom", dependencies=TRUE)
if(!require(purrr)) install.packages("purrr", dependencies=TRUE)
if(!require(multcomp)) install.packages("multcomp", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(rstatix)) install.packages("rstatix", dependencies=TRUE)

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(multcomp)
library(ggplot2)
library(rstatix)

# Read the data and rename columns
data <- read_excel("SepticShockDataR.xlsx", sheet = "Clinical Data") %>%
  rename_all(make.names)

# Convert to tidy format
tidy_data<-data %>% 
  pivot_longer(cols = Mass.1..g.:SBC..mmol.L., names_to = "Observation", values_to = "Values")

# Perform Tukey's HSD test
tidy_tukey<-tidy_data %>% 
  group_by(Observation) %>% 
  tukey_hsd(Values ~ Sex * Treatment)

# Filter to retain sets where at least one p-value is less than 0.05
significant_observations <- tidy_tukey %>%
  group_by(Observation) %>%
  filter(any(p.adj < 0.05)) %>%
  ungroup()

# Extract unique names of significant observations
unique_significant_observations <- significant_observations %>%
  distinct(Observation)

# Save the summary table to a CSV file
write.csv(unique_significant_observations, "SSD_Clinical_sig_observations2.csv", row.names = FALSE)

# Save significant Tukey results to a CSV file
write.csv(significant_observations, "SSD_Clinical_tukey_sig_results2.csv", row.names = FALSE)
