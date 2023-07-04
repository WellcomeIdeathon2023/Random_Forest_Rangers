# Install the necessary packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")

# Load the necessary libraries
library(ggplot2)
library(reshape2)

# Create a data frame with your tasks and timeline
tasks <- data.frame(
  name = factor(c("UX Improvements", 
                  "Implementation of UX Improvements",
                  "Incorporation of Advanced ML & Statistical Methods",
                  "Testing & Iteration I",
                  "Integration of Diverse Data Types & Structures",
                  "Testing & Iteration II",
                  "Finalizing the Platform & Preparing for Deployment",
                  "Deployment & User Training"), 
                levels = c("Deployment & User Training",
                           "Finalizing the Platform & Preparing for Deployment",
                           "Testing & Iteration II",
                           "Integration of Diverse Data Types & Structures",
                           "Testing & Iteration I",
                           "Incorporation of Advanced ML & Statistical Methods",
                           "Implementation of UX Improvements",
                           "UX Improvements")),
  start = c(1, 2, 4, 6, 7, 9, 10, 11),
  end = c(3, 4, 6, 8, 9, 11, 12, 13)
)

# Melt the data to have one row per time point per task
tasks_melt <- melt(tasks, id.vars='name')

# Plot the Gantt chart
ggplot(tasks_melt, aes(value, name)) +
  geom_line(size=10, color="lightgreen") +
  xlab('Month') + ylab('Task') +
  theme_minimal() +
  ggtitle('Project Outline')

ggsave(filename = "/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/results/Gantt Chart.png")
