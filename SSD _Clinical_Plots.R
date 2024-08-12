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

#Set your working directory
setwd("~/R/SepticShock")

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

# Calculate means and standard deviations
mean_sd_table <- tidy_data %>%
  group_by(Observation, `Sex:Treatment`) %>%
  summarize(Mean = mean(Values, na.rm = TRUE),
            Std_Dev = sd(Values, na.rm = TRUE)) %>%
  mutate(Formatted = paste0(round(Mean, 2), " ± ", round(Std_Dev, 2))) %>%
  select(Observation, `Sex:Treatment`, Formatted) %>%
  pivot_wider(names_from = `Sex:Treatment`, 
              values_from = Formatted, 
              names_glue = "{gsub(':', ' ', `Sex:Treatment`)}") %>%
  group_by(Observation) %>%
  summarize(across(everything(), ~ paste(na.omit(.), collapse = ""))) %>%
  ungroup() %>%
  # Replace observation names with original names
  mutate(Observation = name_mapping[Observation])

# Save the Stats table to a CSV file
write.csv(mean_sd_table, "SSD_Clinical_StatsTable.csv", row.names = FALSE, fileEncoding = "UTF-8")

# Perform Tukey's HSD test
tidy_tukey<-tidy_data %>% 
  group_by(Observation) %>% 
  tukey_hsd(Values ~ Sex * Treatment) %>%
  ungroup()


write.csv(tidy_tukey, "SSD_Clinical_TukeyHSD.csv", row.names = FALSE)

# Filter to retain only Group Differences between reasonable groups
tidy_tukey<-tidy_tukey %>%
  filter(term == "Sex:Treatment") %>%
  filter(group1 != "Female:Control" | group2 != "Male:LPS") %>%
  filter(group1 != "Male:Control" | group2 != "Female:LPS") %>%
  filter(p.adj.signif != 'ns') 

# Filter to retain sets where at least one p-value is less than 0.05
sig_observations <- tidy_tukey %>%
  group_by(Observation) %>%
  filter(any(p.adj < 0.05)) %>%
  ungroup()

# Extract unique names of significant observations
sig_observations <- sig_observations %>%
  distinct(Observation)

write.csv(sig_observations, "SSD_Clinical_Sig_Observations.csv", row.names = FALSE)

############
for (observation in sig_observations$Observation) {
  original_observation <- name_mapping[[observation]]
  
  plot_data <- tidy_data %>%
    filter(Observation == observation) %>%
    dplyr::select(Sex, Treatment, 'Sex:Treatment', Values)
  plot_tukey <- tidy_tukey %>%
    filter(Observation == observation) %>%
    add_y_position(test = ., 
                   data = plot_data, 
                   formula = Values ~ Sex:Treatment, 
                   fun = "max", 
                   scales = "free", 
                   comparisons = list(c("Female:LPS", "Male:LPS"),
                                      c("Male:Control", "Male:LPS"), 
                                      c("Female:Control", "Female:LPS"),
                                      c("Female:Control", "Male:Control")), 
                   step.increase = 0.07) %>%
    add_x_position(test = ., x = 'Sex:Treatment')
  
  plot_tukey <- mutate(plot_tukey, xmin = case_when(xmin == 1 ~ 0.8,
                                                    xmin == 2 ~ 1.2,
                                                    xmin == 3 ~ 1.8))
  plot_tukey <- mutate(plot_tukey, xmax = case_when(xmax == 2 ~ 1.2,
                                                    xmax == 3 ~ 1.8,
                                                    xmax == 4 ~ 2.2))
  
  p <- ggplot(plot_data, 
              aes(x = Treatment, y = !!rlang::sym("Values"))) +
    geom_boxplot(outlier.shape = NA, 
                 aes(fill = Sex), 
                 key_glyph = draw_key_rect) +  # Remove outliers from boxplot to avoid duplication with jitter
    geom_point(aes(fill = Sex), 
               size = 1.5, 
               position=position_jitterdodge(jitter.width = 0.2), 
               show.legend = F) +
    labs(title = paste("Results for", original_observation),
         y = original_observation) +
    stat_pvalue_manual(plot_tukey, 
                       label = "p.adj.signif", 
                       y.position = "y.position", 
                       xmin = "xmin", 
                       xmax = "xmax", 
                       hide.ns = TRUE, 
                       tip.length = 0.02, 
                       bracket.nudge.y = 0) +
    scale_fill_manual(values = alpha(c("#21918c","#fde725"), 0.5))+ 
    # scale_fill_viridis(discrete = TRUE, 
    #                    alpha = 0.5, 
    #                    option="viridis") + 
    theme_bw() +
    theme(legend.key = element_rect(color = "black"), 
          legend.key.spacing.y = unit(0.2, "cm"))
  
  print(p)
  
  # Save each plot as a PNG file
  ggsave(filename = paste0("boxplot_", 
                           observation, 
                           ".png"), 
         path = '~/R/SepticShock/Plots/Clinical',
         plot = p, 
         width = 8, 
         height = 6)
}

