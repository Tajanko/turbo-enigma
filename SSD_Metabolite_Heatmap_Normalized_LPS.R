####### PREAMBLE ###############################################################

requiredPackages <- c("readxl", "dplyr", "ggplot2", "viridis", "rstatix", "tidyr", "ggpubr", "stringr", "forcats", "egg")


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(requiredPackages)

#Set your working directory
setwd("~/R/SepticShock")


####### DATA PREP ##############################################################

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


####### ANALYSIS ###############################################################

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

# Convert to tidy format
tidy_table<-mean_table %>% 
  pivot_longer(cols = 'Estradiol LPS':'Male LPS', names_to = "Group", values_to = "Mean") %>%
  select(Observation, Group, Mean)

# Remove LPS tag
tidy_table <- mutate(tidy_table, Group=as.character(Group)) %>%
  mutate(tidy_table,Group=sapply(strsplit(tidy_table$Group, split=' ', fixed=TRUE), function(x) (x[1])))

# Reorder Group levels
tidy_table$Group <- factor(tidy_table$Group, 
                           levels = c("Male", 
                                      "Female", 
                                      "Estradiol"))
tidy_table$Mean <-  as.numeric(tidy_table$Mean)

tidy_table <- tidy_table %>%
  group_by(Observation) %>%
  mutate(Male_Mean = Mean[Group == "Male"]) %>%
  mutate(Ratio = Mean / Male_Mean) %>%
  select(Observation, Group, Mean, Ratio) %>%
  ungroup()

# Assign metabolite Groupings
physio_list <- c("Amino Acids",'Nucleotides','Glycolysis and TCA cycle','Misc Processes','Lipid Biosynthesis','Fatty Acids')
physio_list<-as.data.frame(physio_list)
colnames(physio_list) <- "Physio"
tidy_table$Physio <- 0
tidy_table$Physio[1:60] <- physio_list[1,1]
tidy_table$Physio[61:111] <- physio_list[2,1]
tidy_table$Physio[112:147] <- physio_list[3,1]
tidy_table$Physio[148:225] <- physio_list[4,1]
tidy_table$Physio[226:300] <- physio_list[5,1]
tidy_table$Physio[301:351] <- physio_list[6,1]


####### PLOTTING ###############################################################

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
        plot.title = element_text(hjust = 0.5),
        legend.key.height = 0.2*unit(n_distinct(plot_table$Observation), 'cm'),
        legend.key.width = 0.25*unit(n_distinct(plot_table$Group), 'cm'),
        legend.title=element_blank()) +
  labs(title = paste(physio, "Relative to Baseline"),
       y = '',
       x = '')

print(p)

q <- egg::set_panel_size(p, height=unit(n_distinct(plot_table$Observation), "cm"),
                         width=4*unit(n_distinct(plot_table$Group), "cm") )

# Save each plot as a PNG file
ggsave(filename = paste0("Heatmap_",physio,".png"),
       path = '~/R/SepticShock/Heatmaps/Metabolite',
       plot = q, 
       width = 3*unit(n_distinct(plot_table$Group), "cm"), 
       height = 0.45*unit(n_distinct(plot_table$Observation), "cm"))


p <- ggplot(plot_table, 
            aes(x = Group, y = fct_rev(as_factor(Observation)), fill = Ratio)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(low = '#440154', mid = "white", high = "#fde725", midpoint = 1, 
                       limits = c(0.25, 
                                  1.25),
                       name = "Ratio to Male") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = 0.2*unit(n_distinct(plot_table$Observation), 'cm'),
        legend.key.width = 0.25*unit(n_distinct(plot_table$Group), 'cm'),
        legend.title=element_blank()) +
  labs(title = paste(physio, "Relative to Baseline as a Ratio of Male"),
       y = '',
       x = '')

print(p)

q <- egg::set_panel_size(p, height=unit(n_distinct(plot_table$Observation), "cm"),
                         width=4*unit(n_distinct(plot_table$Group), "cm") )

# Save each plot as a PNG file
ggsave(filename = paste0("HeatmapRatio_",physio,".png"),
       path = '~/R/SepticShock/Heatmaps/Metabolite',
       plot = q, 
       width = 3*unit(n_distinct(plot_table$Group), "cm"), 
       height = 0.45*unit(n_distinct(plot_table$Observation), "cm"))

}
