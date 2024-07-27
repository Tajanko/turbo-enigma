# Load necessary libraries
if(!require(readxl)) install.packages("readxl", dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require(tidyr)) install.packages("tidyr", dependencies=TRUE)
if(!require(broom)) install.packages("broom", dependencies=TRUE)
if(!require(purrr)) install.packages("purrr", dependencies=TRUE)
if(!require(multcomp)) install.packages("multcomp", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(multcomp)
library(ggplot2)

# Read the data from the first sheet "Clinical Data"
data <- read_excel("SepticShockDataR.xlsx", sheet = "Core Data") %>%
  rename_all(make.names)

# Select relevant columns starting from column 5
core_data <- data[,-1:-4]

# Extract column names for clinical observations
observation_columns <- names(core_data)

# Function to perform two-way ANOVA for each observation
perform_anova <- function(column_name) {
  formula <- as.formula(paste(column_name, "~ Sex * Treatment", sep = " "))
  result <- aov(formula, data = data)
  tidy_result <- tidy(result)
  tidy_result <- tidy_result %>% mutate(Observation = column_name)
  return(tidy_result)
}

# Apply the function to all clinical observations
anova_results <- map_dfr(observation_columns, perform_anova)

# Filter to retain sets where at least one p-value is less than 0.05
significant_observations <- anova_results %>%
  group_by(Observation) %>%
  filter(any(p.value < 0.05)) %>%
  ungroup()

# Extract unique names of significant observations
unique_significant_observations <- significant_observations %>%
  distinct(Observation)

# Create a summary table of unique significant observations
summary_table <- unique_significant_observations

# Perform Tukey's HSD test for each significant observation
perform_tukey <- function(column_name) {
  formula <- as.formula(paste(column_name, "~ Sex * Treatment", sep = " "))
  aov_result <- aov(formula, data = data)
  tukey_result <- TukeyHSD(aov_result)
  tidy_tukey <- tidy(tukey_result) %>% mutate(Observation = column_name)
  return(tidy_tukey)
}

# Apply the function to all significant observations
tukey_results <- map_dfr(unique_significant_observations$Observation, perform_tukey)

# Filter Tukey results to retain sets where at least one p-value is less than 0.05
tukey_significant_results <- tukey_results %>%
  filter(adj.p.value < 0.05)

# Save the summary table to a CSV file
write.csv(summary_table, "SSD_Core_sig_observations.csv", row.names = FALSE)

# Save significant ANOVA results to a CSV file
write.csv(significant_observations, "SSD_Core_anova_sig_results.csv", row.names = FALSE)

# Save significant Tukey results to a CSV file
write.csv(tukey_significant_results, "SSD_Core_tukey_sig_results.csv", row.names = FALSE)