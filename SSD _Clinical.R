# Load necessary libraries
if(!require(readxl)) install.packages("readxl", dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(viridis)) install.packages("viridis", dependencies=TRUE)
if(!require(rstatix)) install.packages("rstatix", dependencies=TRUE)
if(!require(tidyr)) install.packages("tidyr", dependencies=TRUE)
if(!require(ggpubr)) install.packages("ggpubr", dependencies=TRUE)


library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)
library(rstatix)
library(tidyr)
library(ggpubr)

#Read the original data from the first sheet "Clinical Data" to get the original column names
original_data <- read_excel("SepticShockDataR.xlsx", sheet = "Clinical Data")

# Create a mapping of original names to modified names
original_names <- names(original_data)[-1:-6]
modified_names <- make.names(original_names)
name_mapping <- setNames(original_names, modified_names)


# Read the data and rename columns
data <- original_data %>%
  rename_all(make.names)

# Convert to tidy format
tidy_data<-data %>% 
  pivot_longer(cols = Mass.1..g.:SBC..mmol.L., names_to = "Observation", values_to = "Values") %>%
  mutate("Sex:Treatment" = paste0(Sex, ":", Treatment)) %>%
  select(!c("Phase.1", "Phase.2")) %>%
  drop_na(Values)

tidy_data$`Sex:Treatment` <- factor(tidy_data$`Sex:Treatment`, levels = c("Female:Control", "Male:Control", "Female:LPS", "Male:LPS"))
tidy_data$Patient.ID <- factor(tidy_data$Patient.ID)
tidy_data$Species <- factor(tidy_data$Species)
tidy_data$Sex <- factor(tidy_data$Sex)
tidy_data$Treatment <- factor(tidy_data$Treatment)
tidy_data$Observation <- factor(tidy_data$Observation)

# Perform Tukey's HSD test
tidy_tukey<-tidy_data %>% 
  group_by(Observation) %>% 
  tukey_hsd(Values ~ Sex * Treatment) %>%
  filter(term == "Sex:Treatment") %>%
  ungroup()

ungroup(tidy_data)
tidy_tukey <- tidy_tukey %>%
  add_y_position(test = ., data = group_by(tidy_data, Observation), formula = Values ~ Sex:Treatment, fun = "max", scales = "free") %>%
  add_x_position(test = ., x = 'Sex:Treatment')

# Filter to retain sets where at least one p-value is less than 0.05
significant_observations <- tidy_tukey %>%
  group_by(Observation) %>%
  filter(any(p.adj < 0.05)) %>%
  filter(term == "Sex:Treatment") %>%
  ungroup()


# Extract unique names of significant observations
summary_table <- significant_observations %>%
  distinct(Observation)


############
for (observation in summary_table$Observation) {
  original_observation <- name_mapping[[observation]]
  
  plot_data <- tidy_data %>%
    filter(Observation == observation) %>%
    dplyr::select(Sex, Treatment, Values)
  plot_tukey <- tidy_tukey %>%
    filter(Observation == observation)
  
  plot_tukey <- mutate(plot_tukey, xmin = case_when(xmin == 1 ~ 0.8,
                                                    xmin == 2 ~ 1.2,
                                                    xmin == 3 ~ 1.8))
  plot_tukey <- mutate(plot_tukey, xmax = case_when(xmax == 2 ~ 1.2,
                                                    xmax == 3 ~ 1.8,
                                                    xmax == 4 ~ 2.2))

  p <- ggplot(plot_data, aes(x = Treatment, y = !!rlang::sym("Values"))) +
    geom_boxplot(outlier.shape = NA, aes(fill = Sex), key_glyph = draw_key_rect) +  # Remove outliers from boxplot to avoid duplication with jitter
    geom_point(aes(fill = Sex), size = 1.5, position=position_jitterdodge(jitter.width = 0.2), show.legend = F) +
    labs(title = paste("Results for", original_observation),
         y = original_observation) +
    stat_pvalue_manual(plot_tukey, label = "p.adj.signif", y.position = "y.position", xmin = "xmin", xmax = "xmax", hide.ns = TRUE, tip.length = 0)+
    scale_fill_viridis(discrete = TRUE, alpha = 0.5, option="viridis") + theme_bw() +
    theme(legend.key = element_rect(color = "black"), legend.key.spacing.y = unit(0.2, "cm"))
  
  print(p)
  
  # Save each plot as a PNG file
  #ggsave(filename = paste0("boxplot_", observation, ".png"), plot = p, width = 8, height = 6)
}

