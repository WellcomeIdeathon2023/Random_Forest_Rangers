# Import tidyverse
library(tidyverse)
install.packages("readxl") # Install the package
library(readxl) 
install.packages("readr") # Install the package
library(readr) # Load the package into the current R session

sdy180_data_Serology_by_Cohort <- read_excel('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy180/resultfiles/neutralizing_antibody_titer_result/Serology_by_Cohort.423737.xlsx')
sdy180_data_mbaa_result <- read_csv('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy180/resultfiles/mbaa_result.csv')
sdy180_data_SDY180_CBC_Results_and_Dictionary <- read_excel('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy180/studyfiles/SDY180_CBC_Results_and_Dictionary.xlsx')

sdy296_data_hai_result <- read_csv('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy296/resultfiles/hai_result.csv')
sdy296_data_neut_ab_titer_result <- read_csv('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy296/resultfiles/neut_ab_titer_result.csv')

sdy301_data_hai_result <- read_csv('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy301/resultfiles/hai_result.csv')
sdy301_data_neut_ab_titer_result <- read_csv('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy301/resultfiles/neut_ab_titer_result.csv')


# Load the datasets

# Extract the subject_accession numbers from each dataset
s1 <- unique(sdy180_data_Serology_by_Cohort$`Sub Org Accession`)
s2 <- unique(sdy296_data_hai_result$SUBJECT_ACCESSION)
s3 <- unique(sdy296_data_neut_ab_titer_result$SUBJECT_ACCESSION)
s4 <- unique(sdy301_data_hai_result$SUBJECT_ACCESSION)
s5 <- unique(sdy301_data_neut_ab_titer_result$SUBJECT_ACCESSION)

# Find the common subject_accession numbers
common_subjects <- s1[s1 %in% s2 & s1 %in% s3 & s1 %in% s4 & s1 %in% s5]

# Check if there are any common subjects
if (length(common_subjects) > 0) {
  cat("Common subject_accession numbers found in the datasets:\n")
  print(common_subjects)
} else {
  cat("No common subject_accession numbers found in the datasets.\n")
}
#All subject accession umbers are unique across the 3 studies 

library(ggplot2)
library(broom)

generate_plots <- function(data, data_name) {
  lapply(unique(data$VIRUS_STRAIN_REPORTED), function(virus_strain) {
    filtered_data <- data[data$VIRUS_STRAIN_REPORTED == virus_strain, ]
    t_test <- t.test(VALUE_REPORTED ~ STUDY_TIME_COLLECTED, data = filtered_data)
    p_value <- tidy(t_test)$p.value
    
    plot <- ggplot(filtered_data, aes(x = factor(STUDY_TIME_COLLECTED), y = VALUE_REPORTED, group = STUDY_TIME_COLLECTED)) +
      geom_boxplot() +
      labs(title = paste("Reported Value Responses for", data_name, virus_strain),
           x = "Study Time Collected",
           y = "Reported Value",
           subtitle = paste("p-value =", formatC(p_value, digits = 4, format = "f"))) +
      theme_bw() +
      ylim(0, 3000)
    
    # Save the plot
    ggsave(paste("/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/results/", data_name, "_", gsub("/", "_", virus_strain), ".png", sep = ""), plot = plot)
    
    return(plot)
  })
}

plots_hai_result <- generate_plots(sdy301_data_hai_result, "sdy301_data_hai_result")
plots_neut_ab_titer_result <- generate_plots(sdy301_data_neut_ab_titer_result, "sdy301_data_neut_ab_titer_result")

for (i in seq_along(plots_hai_result)) {
  print(plots_hai_result[[i]])
}

for (i in seq_along(plots_neut_ab_titer_result)) {
  print(plots_neut_ab_titer_result[[i]])
}

library(ggplot2)
library(broom)

generate_plots <- function(data, data_name) {
  lapply(unique(data$VIRUS_STRAIN_REPORTED), function(virus_strain) {
    filtered_data <- data[data$VIRUS_STRAIN_REPORTED == virus_strain, ]
    t_test <- t.test(VALUE_REPORTED ~ STUDY_TIME_COLLECTED, data = filtered_data)
    p_value <- tidy(t_test)$p.value
    
    plot <- ggplot(filtered_data, aes(x = factor(STUDY_TIME_COLLECTED), y = VALUE_REPORTED, group = STUDY_TIME_COLLECTED)) +
      geom_boxplot() +
      labs(title = paste("Reported Value Responses for", data_name, virus_strain),
           x = "Study Time Collected",
           y = "Reported Value",
           subtitle = paste("p-value =", formatC(p_value, digits = 4, format = "f"))) +
      theme_bw() +
      ylim(0, 3000)
    
    # Save the plot
    ggsave(paste("/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/results/", data_name, "_", gsub("/", "_", virus_strain), ".png", sep = ""), plot = plot)
    
    return(plot)
  })
}

plots_hai_result <- generate_plots(sdy296_data_hai_result, "sdy296_data_hai_result")
plots_neut_ab_titer_result <- generate_plots(sdy296_data_neut_ab_titer_result, "sdy296_data_neut_ab_titer_result")

for (i in seq_along(plots_hai_result)) {
  print(plots_hai_result[[i]])
}

for (i in seq_along(plots_neut_ab_titer_result)) {
  print(plots_neut_ab_titer_result[[i]])
}




# install packages if not already installed
if (!require(tidyverse)) install.packages('tidyverse')

# load the package
library(tidyverse)

# reshape the data to long format
df_long <- sdy180_data_Serology_by_Cohort %>%
  select(-matches("triplicate|FC")) %>%
  select(Timepoint, matches("^Brisbane")) %>%
  pivot_longer(-Timepoint, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("location", "virus", "type"), sep = " ") %>%
  unite("virus_type", virus, type, sep = "_")

# compare HAI to mean VN for each virus at each time point using boxplots
p <- ggplot(df_long, aes(x = Timepoint, y = value, fill = virus_type)) +
  geom_boxplot() +
  facet_wrap(~ location, scales = "free") +
  theme_bw() +
  labs(x = "Time", y = "Value", title = "sdy180 Comparing HAI to mean VN for each virus over time")

# print the plot
print(p)

# save the plot to a file
ggsave(filename = "/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/results/sdy180_Comparing_HAI_to_mean_VN_for_each_virus_over_time.png", plot = p)




# Load required package
library(tidyverse)


# Filter rows where VIRUS_STRAIN_REPORTED equals "Brisbane"
sdy296_filtered <- sdy296_data_hai_result %>% 
  filter(VIRUS_STRAIN_REPORTED == "B/Brisbane/60/2008") %>% 
  select(SUBJECT_ACCESSION, STUDY_TIME_COLLECTED, VIRUS_STRAIN_REPORTED, VALUE_REPORTED) %>%
  mutate(STUDY_TIME_COLLECTED = as.character(STUDY_TIME_COLLECTED)) 

# Transform sdy180 to long format
sdy180_long <- sdy180_data_Serology_by_Cohort %>% 
  gather("key", "value", starts_with("Brisbane")) %>%
  select('Sub Org Accession', 'Timepoint', key, value) 



# Load required packages
library(tidyverse)

# Add a new column to indicate the study
sdy180_long$study <- 'sdy180'
sdy296_filtered$study <- 'sdy296'

# Combine both datasets
combined_data <- rbind(
  sdy180_long %>% select(study, Timepoint, value = value),
  sdy296_filtered %>% select(study, Timepoint = STUDY_TIME_COLLECTED, value = VALUE_REPORTED)
)

# Convert the Timepoint and value variables to numeric
combined_data <- combined_data %>%
  mutate(
    Timepoint = as.numeric(Timepoint),
    value = as.numeric(value)
  )

# Remove NA values
combined_data <- na.omit(combined_data)

# Create a boxplot
ggplot(combined_data, aes(x = Timepoint, y = value, fill = study)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Responses from SDY180 and SDY296 Studies",
    x = "Timepoint",
    y = "Value",
    fill = "Study"
  )
ggsave(filename = "/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/results/Responses from SDY180 and SDY296 Studies.png")


# Load required package
library(dplyr)

# Remove rows with NA in value column
sdy180_long <- sdy180_long %>% 
  filter(!is.na(value))

# Check the result
head(sdy180_long)


# Load required package
library(dplyr)

# Remove rows with 'triplicate' or 'FC' in key column
sdy180_long <- sdy180_long %>% 
  filter(!grepl("triplicate|FC", key))

# Check the result
head(sdy180_long)

# Checking the distribution of 'Timepoint' in each dataset
table(sdy180_long$Timepoint)
table(sdy296_filtered$STUDY_TIME_COLLECTED)

# Load required package
library(dplyr)

# Remove 'd' and convert to numeric
sdy180_long <- sdy180_long %>%
  mutate(Timepoint = as.numeric(gsub("d", "", Timepoint)))

# Check the result
head(sdy180_long)

# Check the distribution of 'Timepoint' again
table(sdy180_long$Timepoint)




