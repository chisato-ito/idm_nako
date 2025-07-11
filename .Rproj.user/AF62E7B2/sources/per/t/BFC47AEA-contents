################################################################################
#
# 3) Read in mortality rate ratio from the Danish registry and model by age, 
#    anxiety and depression, respectively
#
################################################################################

# Mortality rate ratios (MRR) from the Danish Atlas of Disease Mortality are
# available at: https://csievert.shinyapps.io/.
# [Plana-Ripoll O, Dreier JW, Momen NC, Prior A, Weye N, Mortensen PB, et al. 
# Analysis of mortality metrics associated with a comprehensive range of 
# disorders in Denmark, 2000 to 2018: A population-based cohort study. PLoS Med 
# 2022;19:e1004023. https://doi.org/10.1371/journal.pmed.1004023.]

# Load packages
library(here)
library(tidyverse)

# Anxiety ----------------------------------------------------------------------
# Danish data on mortality rates: anxiety (F41)
fct_R_a <- function(a){ 
  return(exp(approx(c(17.5, 27.5, 37.5, 72.5, 97.5), 
                    log(c(1.5, 3.9, 4.4, 2.8, 1.2)), xout = a, rule = 2)$y)) 
  }

# Danish data on mortality rates: depression (F41)
m.f41 <- read.csv(here("data","df_m_F41.csv"))
mrr.df <- m.f41 %>% 
  select(ageVals, MRR) %>% 
  filter(!is.na(MRR)) %>% 
  rename(a = ageVals)

# Visual check
ages <- 19:74
par(mfrow = c(1, 1), las = 1)
matplot(ages, fct_R_a(ages), type = "l", lty = 1, col = "green", 
        xlab = "Age (years)", ylab = "MRR")
matplot(mrr.df$a, mrr.df$MRR, type = "p", pch = 3, add = TRUE)


# Depression -------------------------------------------------------------------
# Danish data on mortality rates: depression (F32)
fct_R_d <- function(a){ 
  return(exp(approx(c(12.5, 32.5, 92.5), c(log(2.7), log(5.3), log(1.5)), 
                    xout = a, rule = 2)$y)) 
  }

# Danish data on mortality rates: depression (F32)
m.f32 <- read.csv(here("data","df_m_F32.csv"))
mrr.df <- m.f32 %>% 
  select(ageVals, MRR) %>% 
  filter(!is.na(MRR)) %>% 
  rename(a = ageVals)

# Visual check
par(mfrow = c(1, 1), las = 1)
matplot(ages, fct_R_d(ages), type = "l", lty = 1, col = "green", 
        xlab = "Age (years)", ylab = "MRR")
matplot(mrr.df$a, mrr.df$MRR, type = "p", pch = 3, add = TRUE)
