#--------------------------------------------------------------------------------------#
### This file contains simulation code for deriving properties of the national poverty index when:
### 1. k = 1 
### 2. the true distributions have equal means and equal variances (population param are equal)
### 3. sample size is increasing, i.e. n = 2,4,...

library(tidyverse)
library(ggplot2)

main <- function(n_1, z_1){
    mainResults <- matrix(NA, nrow=1, ncol=18)
    currentRow <- 1
    
    mu = 11.0302
    sigma = sqrt(0.6688)
    
    y_deviates = rlnorm(n_1, meanlog = mu , sdlog = sigma) # generating lognormal deviates
    
    mu_hat <- function(y_deviates, n_1){
        mu_hat <- 0
        for(j in 1:n_1){
            mu_hat <- mu_hat + log(y_deviates[j])/n_1
        }
        return(mu_hat)
    } # mle - mean of the lognormal
    
    sigma_sq_hat <- function(y_deviates, n_1){
        sigma_sq_hat <- 0
        for(j in 1:n_1){
            sigma_sq_hat <- sigma_sq_hat + (log(y_deviates[j])-mu_hat_val)^2/(n_1)
        }
        return(sigma_sq_hat)
    } # mle - st deviation of the lognormal
    
    mu_hat_val = mu_hat(y_deviates, n_1) # calculate the sample mean
    sigma_sq_hat_val = sigma_sq_hat(y_deviates, n_1) # calculate the sample variance
    sigma_hat_val = sqrt(sigma_sq_hat_val) # calculate the sample sd
    
    # poverty index
    sq_pov_index <- function(z_1, mu_1, sigma_1) {
        g <- function(y) ((z_1-y)/z_1)^2*dlnorm(y, meanlog = mu_1, sdlog = sigma_1)
        integrate(g, lower=0, upper = z_1)$value
    }
    
    sq_pov_index_true = sq_pov_index(z_1, mu, sigma) # calculate true pov index
    sq_pov_index_est = sq_pov_index(z_1, mu_hat_val, sigma_hat_val) # calculate pov estimate
    
    #---------------------------------------------------------------------------------------#
    
    # decomposed poverty index
    
    decomp_pov_index <- function(z_1, mu_1, sigma_1){
        p_2 <- pnorm((log(z_1)-mu_1)/sigma_1) -
            (2/z_1)*exp(mu_1+sigma_1^2/2)*(pnorm((log(z_1)-mu_1-sigma_1^2)/sigma_1)) +
            (1/z_1^2)*exp(2*mu_1+2*sigma_1^2)*((pnorm((log(z_1)-mu_1-2*sigma_1^2)/sigma_1)))
        return(p_2)
    }
    
    true_decomp_pov_index = decomp_pov_index(z_1, mu, sigma) 
    est_decomp_pov_index = decomp_pov_index(z_1, mu_hat_val, sigma_hat_val)
    
    mainResults[currentRow,] <- c(n_1, z_1, sq_pov_index_true, sq_pov_index_est, 
                                  true_decomp_pov_index, est_decomp_pov_index,
                                  mu_hat_val, sigma_hat_val)
    currentRow <- currentRow + 1
    
    return(mainResults)
    
}

df_sims <- data.frame()

# function to create poverty line
create_z <- function(hc_ratio){
    z = vector()
    for(i in 1:length(hc_ratio)){
        z[i] <- qlnorm(hc_ratio[i], meanlog = 11.0302, sdlog = sqrt(0.6688))
    }
    return(z)
}

n = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50) # vector of sample sizes
z = create_z(c(0.1, 0.2)) # vector of predefined poverty lines 

# 10,000 iterations for each combination of n and z
for(k in 1:length(n)){
    for(m in 1:length(z)){
        for (l in 1:10000){
            currentResults <- main(n[k], z[m])
            df_sims <- rbind(df_sims, currentResults) 
        } 
    }
}
