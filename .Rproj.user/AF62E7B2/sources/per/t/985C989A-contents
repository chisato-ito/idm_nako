################################################################################
#
# 4a) re-estimate age- and sex-specific prevalence incorporating sensitivity &
#     specificity
#
################################################################################

# NOTE: Code partially redacted due to data protection and privacy reasons.

# Sensitivity and specificity of questionnaires
# Anxiety
# sensitivity: 89%; specificity: 82% (Spitzer et al. 2006)
se.a <- .89
sp.a <- .82

# Depression
# sensitivity: 78%; specificity: 87% (Moriarty et al. 2015)
se.d <- .78
sp.d <- .87

# The prevalence incorporating sensitivity and specificity of the tool can be 
# obtained from the observed prevalence by the following equation:
# p = (p_obs - 1 + sp)/(se + sp - 1) 

merged_prev_a_ <- merged_prev_a %>% 
  mutate(prev.sesp = (prev - 1 + sp.a)/(se.a + sp.a - 1))

merged_prev_d_ <- merged_prev_d %>% 
  mutate(prev.sesp = (prev - 1 + sp.d)/(se.d + sp.d - 1))
