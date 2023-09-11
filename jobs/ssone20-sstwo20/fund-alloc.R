# Load libraries
library(tidyverse)
library(ggplot2)


# Two subgroup distributions have same means and same variances  
mu_1 = 11.0302
sigma_1 = sqrt(0.6688)

mu_2 = 11.0302
sigma_2 = sqrt(0.6688)

ssone = 20
sstwo = 20
taskid <- commandArgs(trailingOnly=TRUE)

#--------------------------------------------------------------------#

# Decomposed poverty index function
sq_pov_index <- function(z_1, b_1, n, mu, sigma){
    p_2 <- (1-(b_1/(n*z_1)))^2*pnorm((log(z_1-(b_1/n))-mu)/sigma) -
        (2/z_1)*(1-(b_1/(n*z_1)))*exp(mu+sigma^2/2)*(pnorm((log(z_1-(b_1/n))-mu-sigma^2)/sigma)) +
        (1/z_1^2)*exp(2*mu+2*sigma^2)*((pnorm((log(z_1-(b_1/n))-mu-2*sigma^2)/sigma)))
    return(p_2)
}

pov_index <- function(z_1, b_1, n, mu, sigma){
    p_1 <- (1-(b_1/(n*z_1)))*pnorm((log(z_1-(b_1/n))-mu)/sigma) -
        (1/z_1)*exp(mu+sigma^2/2)*(pnorm((log(z_1-(b_1/n))-mu-sigma^2)/sigma)) 
    return(p_1)
}

# Gradient function
gradient_f <- function(z_1, b_1, n, n_sample, mu, sigma){
    g <- (1-(b_1/(n*z_1)))*pnorm((log(z_1-(b_1/n))-mu)/sigma) -
        (1/z_1)*(1+((sigma^2)/(2*n_sample))+((sigma^4)/(4*(n_sample-1))))*exp(mu+((sigma^2)/2))*
        (pnorm((log(z_1-(b_1/n))-mu-(sigma^2))/sigma)) +
        z_1*(1-(b_1/(n*z_1)))^2*dlnorm(z_1-(b_1/n), meanlog = mu, sdlog = sigma)*
        (((sigma^4)/(4*(n_sample-1)))+((sigma^2)/(2*n_sample))*(((n_sample*(log(z_1-(b_1/n))-mu))/(2*(n_sample-1)))-(n_sample/(2*(n_sample-1)))+1)+
             (((log(z_1-(b_1/n))-mu)^2)/(4*(n_sample-1))))
    return(g)
}

# national expected poverty index function for two groups
exp_pov_index <- function(z_1, b_1, n, n_sample, mu, sigma){
    f_1 <- n*(1-(b_1/(n*z_1)))^2*pnorm((log(z_1-(b_1/n))-mu)/sigma)-
        (n/z_1)*(1-(b_1/(n*z_1)))*(2+((sigma^2)/n_sample)+((sigma^4)/(2*(n_sample-1))))*exp(mu+(sigma^2/2))*(pnorm((log(z_1-(b_1/n))-mu-(sigma^2))/sigma))+
        (n/(z_1^2))*(1+(2*((sigma^2)/n_sample))+((4*(sigma^4))/(n_sample-1)))*exp(2*mu+2*(sigma^2))*(pnorm((log(z_1-(b_1/n))-mu-2*(sigma^2))/sigma))+
        n*((sigma^4)/(n_sample-1))*(1-(b_1/(n*z)))^2*((-log(z-(b_1/n))+mu-3*(sigma^2))/(2*(sigma^2)))*(z_1-(b_1/n))*dlnorm(z-(b_1/n), meanlog=mu, sdlog = sigma)
    return(f_1)
}


main <- function(n_1, n_2, z_1, b_1, mu_hat_1, mu_hat_2, sigma_hat_1, sigma_hat_2){
    mainResults <- matrix(NA, nrow = 1, ncol = 32)
    currentRow <- 1
    
    #---------------------------------------------------------------------------------------#
    
    true_pov_index_1 = sq_pov_index(z_1, b_1, N/2, mu_1, sigma_1)   
    est_pov_index_1 = sq_pov_index(z_1, b_1, N/2, mu_hat_1, sigma_hat_1)
    
    true_pov_index_2 = sq_pov_index(z_1, fixed_budget - b_1, N/2, mu_2, sigma_2)   
    est_pov_index_2 = sq_pov_index(z_1, fixed_budget - b_1, N/2, mu_hat_2, sigma_hat_2)
    
    # Compute national poverty index
    national_p2 = true_pov_index_1 # true national poverty index
    est_national_p2 = 0.5*est_pov_index_1+0.5*est_pov_index_2 # equal shares of pop in two districts
    
    national_p2_initial = sq_pov_index(z_1, 0, N/2, mu_1, sigma_1)   
    est_national_p2_initial = 0.5*sq_pov_index(z_1, 0, N/2, mu_hat_1, sigma_hat_1)+
        0.5*sq_pov_index(z_1, 0, N/2, mu_hat_2, sigma_hat_2)
    
    true_kanbur_needs_gr1 = pov_index(z_1, b_1, N/2, mu_1, sigma_1)
    true_kanbur_needs_gr2 = pov_index(z_1, fixed_budget - b_1, N/2, mu_2, sigma_2)
    
    est_kanbur_needs_gr1 = pov_index(z_1, b_1, N/2, mu_hat_1, sigma_hat_1)
    est_kanbur_needs_gr2 = pov_index(z_1, fixed_budget - b_1, N/2, mu_hat_2, sigma_hat_2)
    
    #---------------------------------------------------------------------------------------#
    
    # Compute true gradients in groups 1 and 2 post-transfer
    true_needs_gr1 <- gradient_f(z_1, b_1, N/2, n_1, mu_1, sigma_1)
    true_needs_gr2 <- gradient_f(z_1, fixed_budget - b_1, N/2, n_2, mu_2, sigma_2)
    
    # Compute estimated gradients in groups 1 and 2 post-transfer
    est_needs_gr1 <- gradient_f(z_1, b_1, N/2, n_1, mu_hat_1, sigma_hat_1)
    est_needs_gr2 <- gradient_f(z_1, fixed_budget-b_1, N/2, n_2, mu_hat_2, sigma_hat_2)
    
    # Compute true gradients in groups 1 and 2 pre-transfer
    true_needs_gr1_initial <- gradient_f(z_1, 0, N/2, n_1, mu_1, sigma_1)
    true_needs_gr2_initial <- gradient_f(z_1, 0, N/2, n_2, mu_2, sigma_2)
    
    # Compute estimated gradients in groups 1 and 2 pre-transfer
    est_needs_gr1_initial <- gradient_f(z_1, 0, N/2, n_1, mu_hat_1, sigma_hat_1)
    est_needs_gr2_initial <- gradient_f(z_1, 0, N/2, n_2, mu_hat_2, sigma_hat_2)
    
    # Compute true and estimated national expected squared poverty gap
    true_national_exp_pov_index = (exp_pov_index(z_1, b_1, N/2, n_1, mu_1, sigma_1) +
        exp_pov_index(z_1, fixed_budget-b_1, N/2, n_2, mu_2, sigma_2))/N
    est_national_exp_pov_index = (exp_pov_index(z_1, b_1, N/2, n_1, mu_hat_1, sigma_hat_1) +
        exp_pov_index(z_1, fixed_budget-b_1, N/2, n_2, mu_hat_2, sigma_hat_2))/N
    
    # Compute true and estimated national expected squared poverty gap pre-transfer
    true_national_exp_pov_index_initial = (exp_pov_index(z_1, 0, N/2, n_1, mu_1, sigma_1) +
        exp_pov_index(z_1, 0, N/2, n_2, mu_2, sigma_2))/N
    est_national_exp_pov_index_initial = (exp_pov_index(z_1, 0, N/2, n_1, mu_hat_1, sigma_hat_1) +
        exp_pov_index(z_1, 0, N/2, n_2, mu_hat_2, sigma_hat_2))/N
    
    
    mainResults[currentRow,] <- c(n_1, n_2, z_1, b_1, mu_hat_1, mu_hat_2, sigma_hat_1, sigma_hat_2, 
                                  true_pov_index_1, est_pov_index_1, true_pov_index_2, est_pov_index_2,
                                  national_p2, est_national_p2, national_p2_initial, est_national_p2_initial,
                                  true_kanbur_needs_gr1, est_kanbur_needs_gr1, true_kanbur_needs_gr2, est_kanbur_needs_gr2,
                                  true_needs_gr1, true_needs_gr2, est_needs_gr1, est_needs_gr2,
                                  true_needs_gr1_initial, true_needs_gr2_initial, est_needs_gr1_initial, est_needs_gr2_initial,
                                  true_national_exp_pov_index, est_national_exp_pov_index, true_national_exp_pov_index_initial,
                                  est_national_exp_pov_index_initial)
    currentRow <- currentRow + 1
    
    return(mainResults)
    
}

# Function to create poverty line
create_z <- function(hc_ratio){
    z = vector()
    for(i in 1:length(hc_ratio)){
        z[i] <- qlnorm(hc_ratio[i], meanlog = 11.0302, sdlog = sqrt(0.6688))
    }
    return(z)
}

z = create_z(0.1) # Vector of predefined poverty lines 

# The budget is set equal to 5% of the per-capita income value of the 25^th percentile 
# of the income distribution, scaled by the total population in a country. 

# Function to create the total fixed budget available for allocation
create_fixed_budget <- function(percentile_inc_dst){
    fixed_budget_1 = vector()
    for(i in 1:length(percentile_inc_dst)){
        fixed_budget_1[i] <- qlnorm(percentile_inc_dst[i], meanlog = 11.0302, sdlog = sqrt(0.6688))
    }
    return(fixed_budget_1)
}

fixed_budget_1 = create_fixed_budget(percentile_inc_dst = 0.25)

N = 100 # total number of people in nation
fixed_budget <- c(0.05*N*fixed_budget_1) # Fixed budget externally set

B_1 = c(seq(0, fixed_budget, 1)) 

generate_mu_hat_val_1 <- function(ssone){
    return(rnorm(1, mean = 11.0302, sd = sqrt(0.6688)/ssone))
}

generate_mu_hat_val_2 <- function(sstwo){
    return(rnorm(1, mean = 11.0302, sd = sqrt(0.6688)/sstwo))
}

generate_sigma_hat_val_1 <- function(ssone){
    return(rchisq(1, ssone-1))
}

generate_sigma_hat_val_2 <- function(sstwo){
    return(rchisq(1, sstwo-1))
}

df_sims <- data.frame()

for(g in 1:1){
    mu_hat_val_1 = generate_mu_hat_val_1(ssone)
    mu_hat_val_2 = generate_mu_hat_val_2(sstwo)
    sigma_hat_val_1 = sqrt(sigma_1^2*generate_sigma_hat_val_1(ssone)/(ssone-1))
    sigma_hat_val_2 = sqrt(sigma_2^2*generate_sigma_hat_val_2(sstwo)/(sstwo-1))
    
    for(a in 1:length(B_1)){
        for(m in 1:length(z)){
                currentResults <- main(ssone, sstwo, z[m], B_1[a], mu_hat_val_1, mu_hat_val_2, 
                                       sigma_hat_val_1, sigma_hat_val_2)
                df_sims <- rbind(df_sims, currentResults) 
        }
    }
}

saveRDS(df_sims, file=paste0("pov-ssone",ssone,"-sstwo",sstwo,"-task",taskid,".RDS"))



