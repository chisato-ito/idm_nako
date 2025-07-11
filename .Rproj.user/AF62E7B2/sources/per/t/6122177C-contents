################################################################################
#
# 5) Model prevalence & derivative then model incidence & remission and optimize
# ANXIETY & DEPRESSION
#
################################################################################

# Here we present the anxiety model as a sample, using Augsburg (center 11) as 
# a reference.

# Auxiliary functions
logit <- function(x){log(x/(1-x))}
expit <- function(x){y<-exp(x); y/(1+y)}

# Fit a prevalence model
# Dataframe p.df is the weighted NAKO prevalence data, but some of the variables
# are renamed (a = age, t = time, sx = sex, sc = study center)
mod.p_a  <- lm(logit(p) ~ (I(a^2) + a) + t + sx + as.factor(sc), data = p.df, 
               subset = (p > 0 & p < 1))
summary(mod.p_a)

# The following is hard coded with Augsburg (center 11) as a reference for code 
# sharing purposes
fct_p_a <- function(t, a, sx){
  logitP <- -10.84  -1.034e-03*a^2+ 6.991e-02*a + 3.193e-03*t + 4.393e-01*sx
  return(expit(logitP))
}

fct_dp_a <- function(t, a, sx){
  return((fct_p_a(t+.001, a+.001, sx) - fct_p_a(t-.001, a-.001, sx))/.002)
}

# Check
sex <- 1 # 1: Men; 2: Women
year <- 2015 # 2018

matplot(ages, 1e3*fct_p_a(year, ages, sex, c), type = "l", lty = 2, 
        col = "#A01A9CFF",
        ylim = c(0, 200), main = paste(sex,"-",year,": Center", c), 
        xlab = "Age (years)", ylab = "Prevalence (per 1000)")
with(p.df, matplot(a[t == year & sx == sex & sc == c],
                   1e3*p[t == year & sx == sex & sc == c], type = "p", 
                   col = "#B52F8CFF", pch = 4, add = TRUE))
for (ii in 1:length(unique(p.df$a))) {
  aa <- (p.df$a    [p.df$t == year & p.df$sx == sex & p.df$sc == c]) [ii]
  up <- (p.df$p_upp[p.df$t == year & p.df$sx == sex & p.df$sc == c]) [ii]
  lw <- (p.df$p_low[p.df$t == year & p.df$sx == sex & p.df$sc == c]) [ii]
  lines(rep(aa,2), 1e3*c(lw,up), col = "#B52F8CFF")
}

# Model incidence and remission as two Gaussian --------------------------------
# Women
fct_i_a_w <- function(beta1 = 30e-4, beta2 = 15, beta3 = 15){
  return(beta1*exp(-0.5*((age-beta2)/beta3)^2))
}
fct_r_a_w <- function(beta1 = 10e-4, beta2 = 58, beta3 = 7){
  return(beta1*exp(-0.5*((age-beta2)/beta3)^2))
}
# Men
fct_i_a_m <- function(beta1 = 30e-4, beta2 = 10, beta3 = 20){
  return(beta1*exp(-0.5*((age-beta2)/beta3)^2))
}
fct_r_a_m <- function(beta1 = 7e-4, beta2 = 58, beta3 = 7){
  return(beta1*exp(-0.5*((age-beta2)/beta3)^2))
}

# Optimize incidence and remission parameters ----------------------------------
sex <- 2 # 1: Men; 2: Women
age <- 19:74
centers <- c(11) # 11: Augsburg
MRR <- fct_R_a(age)


dp_ <- fct_dp_a(2016.5, age, sex)
p_  <- fct_p_a (2016.5, age, sex)
m_  <- exp(fct_logm_w(age))
  
fct_target_a <- function(x){
    rhs_ <- (1-p_)*(fct_i_a_w(x[1],x[2],x[3]) - m_*p_*(MRR-1)/(1+p_*(MRR-1))) - fct_r_a_w(x[4],x[5],x[6])*p_
    return(sum( (dp_ - (rhs_))^2))
}
  
x0 <- c(0.003,15,15,0.001,58,7)
optim_res <- optim(x0, fct_target_a, control = list(maxit = 5000), hessian = TRUE)
print(optim_res$par)
print(optim_res$convergence)
  
par(mfrow = c(1, 2), las = 1)
# Prevalence derivative
matplot(age, 1e4*dp_, type = "l", col = "blue", las = 1, 
        ylab = "", xlab = "Age (years)",
        main = paste("Women: Center", c))
x <- optim_res$par
rhs_ <- (1-p_)*(fct_i_a_w(x[1],x[2],x[3]) - m_*p_*(MRR-1)/(1+p_*(MRR-1))) - fct_r_a_w(x[4],x[5],x[6])*p_
matplot(age, 1e4*rhs_, type = "l", lty = 2, lwd = 2, add = TRUE)
legend("topright", legend = c("Derivative", "Right hand side"), lwd = c(1, 2), 
       col = c("blue", "black"), lty = c(1,2), bty="n")
  
# Incidence & remission
matplot(age, 1e2*fct_r_a_w(x[4],x[5],x[6]), type = "l", col = "blue", las = 1, 
        main = paste("Women: Center", c), ylab = "Rate (per 100 py)", 
        xlab = "Age (years)", ylim = c(0, 20))
matplot(age, 1e3*fct_i_a_w(x[1],x[2],x[3]), type = "l", lty = 2, lwd = 2, 
        add = TRUE)
legend("topleft", legend = c("Remission", "Incidence (x 10)"), lwd = c(1, 2), 
       col = c("blue", "black"), lty = c(1,2), bty="n")


# Repeat with the resampled datasets -------------------------------------------
