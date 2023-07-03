
#Analysis of mbaa data file - cytokines responses 
mbaa <- read_csv('/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy180/resultfiles/mbaa_result.csv')

# Load the necessary library
library(ggplot2)
library(dplyr)

# Convert the necessary columns to numeric
mbaa$STUDY_TIME_COLLECTED <- as.numeric(mbaa$STUDY_TIME_COLLECTED)

# Handle "OOR" values in CONCENTRATION_VALUE_REPORTED column
mbaa$CONCENTRATION_VALUE_REPORTED[mbaa$CONCENTRATION_VALUE_REPORTED == "OOR"] <- NA
mbaa$CONCENTRATION_VALUE_REPORTED <- as.numeric(mbaa$CONCENTRATION_VALUE_REPORTED)

ggplot(mbaa, aes(x = STUDY_TIME_COLLECTED, y = CONCENTRATION_VALUE_REPORTED)) +
  geom_line() +
  facet_wrap(~ ANALYTE_REPORTED) +
  scale_y_continuous(trans = 'log10', limits = c(1, 10000000)) +
  labs(x = "Study Time Collected", y = "Concentration Value Reported") +
  theme_minimal()






library(dplyr)

# Filter out the "standard curve" rows and set STUDY_TIME_COLLECTED to 1 for "control" rows
mbaa_filtered <- mbaa %>% 
  filter(SOURCE_TYPE != "STANDARD CURVE") %>%
  mutate(STUDY_TIME_COLLECTED = ifelse(SOURCE_TYPE == "CONTROL SAMPLE", 1, STUDY_TIME_COLLECTED))

ggplot() +
  geom_line(data = mbaa_filtered %>% filter(SOURCE_TYPE == "EXPSAMPLE"), 
            aes(x = STUDY_TIME_COLLECTED, y = CONCENTRATION_VALUE_REPORTED, color = "EXPSAMPLE")) +
  geom_line(data = mbaa_filtered %>% filter(SOURCE_TYPE == "CONTROL SAMPLE"), 
            aes(x = STUDY_TIME_COLLECTED, y = CONCENTRATION_VALUE_REPORTED, color = "+ve CONTROL SAMPLE")) +
  facet_wrap(~ ANALYTE_REPORTED) +
  scale_y_continuous(trans = 'log10', limits = c(0.01, 100000)) +
  labs(x = "Study Time Collected", y = "Concentration Value Reported", color = "Source Type") +
  theme_minimal()
ggsave(filename = "/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/results/sdyCytokineResponses.png")
