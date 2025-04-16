library(data.table)
library(Matrix)
library(dplyr)
library(SKAT)
library(SPAtest)

# Load the modified MetaSAIGE.R file
source('R/MetaSAIGE.R')

# Print the function signature of load_cohort to confirm it's defined
cat("Load cohort function defined:", exists("load_cohort"), "\n")
print(args(load_cohort))

cat("Testing complete!") 