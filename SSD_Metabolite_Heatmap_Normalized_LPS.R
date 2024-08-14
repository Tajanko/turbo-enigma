# Load necessary libraries
if(!require(readxl)) install.packages("readxl", dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(viridis)) install.packages("viridis", dependencies=TRUE)
if(!require(rstatix)) install.packages("rstatix", dependencies=TRUE)
if(!require(tidyr)) install.packages("tidyr", dependencies=TRUE)
if(!require(ggpubr)) install.packages("ggpubr", dependencies=TRUE)
if(!require(stringr)) install.packages("stringr", dependencies=TRUE)
if(!require(forcats)) install.packages("forcats", dependencies=TRUE)

library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)
library(rstatix)
library(tidyr)
library(ggpubr)
library(stringr)
library(forcats)

#Set your working directory
setwd("~/R/SepticShock")

# Read the original data from the first sheet "Clinical Data" to get the original column names
original_data <- read_excel("SepticShockDataR.xlsx", sheet = "Normalized Metabolite Data")

# Create a mapping of original names to modified names
original_names <- names(original_data)[-1:-5]
modified_names <- make.names(original_names)
name_mapping <- setNames(original_names, modified_names)

# Read the data and rename columns
data <- original_data %>%
  rename_all(make.names)
data[is.na(data)] <- 1

# Convert to tidy format
tidy_data<-data %>% 
  pivot_longer(cols = L.alanine:X.5Z.8Z.11Z.14Z.17Z..Icosapentaenoic.acid, names_to = "Observation", values_to = "Values") %>%
  mutate("Sex:Treatment" = paste0(Sex, ":", Treatment)) %>%
  drop_na(Values)

tidy_data$`Sex:Treatment` <- factor(tidy_data$`Sex:Treatment`, 
                                    levels = c("Female:Control", 
                                               "Male:Control", 
                                               "Estradiol:LPS", 
                                               "Female:LPS", 
                                               "Male:LPS"))
tidy_data$Patient.ID <- factor(tidy_data$Patient.ID)
tidy_data$Species <- factor(tidy_data$Species)
tidy_data$Sex <- factor(tidy_data$Sex)
tidy_data$Treatment <- factor(tidy_data$Treatment)
tidy_data$Observation <- factor(tidy_data$Observation)


# Calculate means
mean_table <- tidy_data %>%
  group_by(Observation, `Sex:Treatment`) %>%
  summarize(Mean = mean(Values, na.rm = TRUE)) %>%
  mutate(Formatted = paste0(round(Mean, 2))) %>%
  select(Observation, `Sex:Treatment`, Formatted) %>%
  pivot_wider(names_from = `Sex:Treatment`, 
              values_from = Formatted, 
              names_glue = "{gsub(':', ' ', `Sex:Treatment`)}") %>%
  group_by(Observation) %>%
  summarize(across(everything(), ~ paste(na.omit(.), collapse = ""))) %>%
  ungroup() %>%
  # Replace observation names with original names
  mutate(Observation = recode(Observation, !!!name_mapping)) %>%
# Add a factor to reorder according to the name mapping
  mutate(Observation = factor(Observation, levels = original_names)) %>%
  arrange(Observation)  # Reorder based on the original names

tidy_table<-mean_table %>% 
  pivot_longer(cols = 'Estradiol LPS':'Male LPS', names_to = "Group", values_to = "Mean") %>%
  select(Observation, Group, Mean)

tidy_table <- mutate(tidy_table, Group=as.character(Group)) %>%
  mutate(tidy_table,Group=sapply(strsplit(tidy_table$Group, split=' ', fixed=TRUE), function(x) (x[1])))

tidy_table$Group <- factor(tidy_table$Group, 
                           levels = c("Male", 
                                      "Female", 
                                      "Estradiol"))
tidy_table$Mean <-  as.numeric(tidy_table$Mean)

physio_list <- c("Amino Acids",'Nucleotides','Glycolysis and TCA cycle','Misc Processes','Lipid Metabolism and Biosynthesis','Fatty Acids')
physio_list<-as.data.frame(physio_list)
colnames(physio_list) <- "Physio"
tidy_table$Physio <- 0
tidy_table$Physio[1:60] <- physio_list[1,1]
tidy_table$Physio[61:111] <- physio_list[2,1]
tidy_table$Physio[112:147] <- physio_list[3,1]
tidy_table$Physio[148:225] <- physio_list[4,1]
tidy_table$Physio[226:300] <- physio_list[5,1]
tidy_table$Physio[301:351] <- physio_list[6,1]




for (physio in physio_list$Physio) {

plot_table <- tidy_table %>%
  filter(Physio == physio)

p <- ggplot(plot_table, 
             aes(x = Group, y = fct_rev(as_factor(Observation)), fill = Mean)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_viridis(alpha = 0.5,
                     option="viridis") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = paste(physio, "Relative to Baseline"),
       y = '',
       x = '')

print(p)

# Save each plot as a PNG file
ggsave(filename = paste0("Heatmap_",physio,".png"),
       path = '~/R/SepticShock/Heatmaps/Metabolite',
       plot = p, 
       width = 8, 
       height = 6)
}
