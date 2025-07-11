################################################################################
#
# 4) Resample data and estimate age- and sex-specific prevalence x 2000 
#
################################################################################

# NOTE: Code partially redacted due to data protection and privacy reasons.

# Load package
library(srvyr)

# Read in the cleaned up NAKO data 

# Resample and repeat the analysis 2000 times
n_iterations <- 2000 

# Placeholder for the results
resample_prev_a <- list() # anxiety
resample_prev_d <- list() # depression

# For reproducibility
set.seed(123)

# Perform resampling and analysis
for (i in 1:n_iterations) {
  n <- nrow(dat)
  
  # Resample the dataset with replacement
  resample <- dat[sample.int(n, n, replace = TRUE), ]
  
  # Compute the total weight for each study center in the ith sample

  # Merge the original and new study center weights

  # Create a survey design object for the ith sample

  # Estimate the age- and sex-specific anxiety prevalence for the ith sample

  # Estimate the age- and sex-specific depression prevalence for the ith sample
  
  # Store the result for the current iteration
  resample_prev_a[[i]] <- prev_a
  resample_prev_d[[i]] <- prev_d
}

# Combine all the results into one dataframe
merged_prev_a <- bind_rows(resample_prev_a) # anxiety
merged_prev_d <- bind_rows(resample_prev_d) # depression
