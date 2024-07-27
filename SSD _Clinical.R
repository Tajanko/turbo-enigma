# Load necessary libraries
if(!require(readxl)) install.packages("readxl", dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(viridis)) install.packages("viridis", dependencies=TRUE)
if(!require(rstatix)) install.packages("rstatix", dependencies=TRUE)


library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)
library(rstatix)

# Read the data and rename columns
data <- read_excel("SepticShockDataR.xlsx", sheet = "Clinical Data") %>%
  rename_all(make.names)

# Convert to tidy format
tidy_data<-data %>% 
  pivot_longer(cols = Mass.1..g.:SBC..mmol.L., names_to = "Observation", values_to = "Values") %>%
  mutate("Sex:Treatment" = paste0(Sex, ":", Treatment)) %>%
  select(!c("Phase.1", "Phase.2")) %>%
  drop_na(Values)

tidy_data$`Sex:Treatment` <- factor(tidy_data$`Sex:Treatment`)
tidy_data$Patient.ID <- factor(tidy_data$Patient.ID)
tidy_data$Species <- factor(tidy_data$Species)
tidy_data$Sex <- factor(tidy_data$Sex)
tidy_data$Treatment <- factor(tidy_data$Treatment)
tidy_data$Observation <- factor(tidy_data$Observation)

# Read the summary table of unique significant observations
#tidy_tukey <- read.csv("SSD_Clinical_tukey_sig_results2.csv")

# Perform Tukey's HSD test
tidy_tukey<-tidy_data %>% 
  group_by(Observation) %>% 
  tukey_hsd(Values ~ Sex * Treatment) %>%
  filter(term == "Sex:Treatment") %>%
  ungroup()

# Filter to retain sets where at least one p-value is less than 0.05
significant_observations <- tidy_tukey %>%
  group_by(Observation) %>%
  filter(any(p.adj < 0.05)) %>%
  filter(term == "Sex:Treatment") %>%
  ungroup()

ungroup(tidy_data)
tidy_tukey <- tidy_tukey %>%
  # add_y_position(test = ., data = tidy_data, formula = Values ~ group, fun = "max")
add_y_position(test = ., data = group_by(tidy_data, Observation), formula = Values ~ Sex:Treatment, fun = "max")

