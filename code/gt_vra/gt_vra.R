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

# Combine the hai_result and neut_ab_titre_result datasets
combined_data <- rbind(sdy296_data_hai_result, sdy296_data_neut_ab_titer_result)

# Create separate plots for each virus strain reported
plots <- lapply(unique(combined_data$VIRUS_STRAIN_REPORTED), function(virus_strain) {
  # Filter the data for the current virus strain
  filtered_data <- combined_data[combined_data$VIRUS_STRAIN_REPORTED == virus_strain, ]
  
  # Create the plot
  plot <- ggplot(filtered_data, aes(x = STUDY_TIME_COLLECTED, y = VALUE_REPORTED, group = STUDY_TIME_COLLECTED )) +
    geom_boxplot() +
    labs(title = paste("Reported Value Responses for sdy296", virus_strain),
         x = "Study Time Collected",
         y = "Reported Value") +
    theme_bw()
  
  return(plot)
})

# Print the plots
for (i in seq_along(plots)) {
  print(plots[[i]])
}




