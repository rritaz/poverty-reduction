library(matlib)
# Run algorithm with true means and variances in the two groups
# Note that the population moments are the same, the groups are at an optimum 
# (i.e. their poverty measures are equal)

# At iteration 0 (i.e. t-1 = 1-1), the cash transfer is $0
# Stop the algorithm when the difference between lambda at iteration t and t-1 is within some small threshold

# Two subgroup distributions have same means and same variances  
# first term is true mean, first quartile, mean of mean, third quartile
mu_1 = 11.0302
sigma_1 = sqrt(0.6687999)

mu_2 = 11.0302
sigma_2 = sqrt(0.6687999)


#------------------------------------------------------------------------#

# Average post-transfer income of individuals below the poverty line for j-th group, where j=1,2

avg_pti1_prev_iter <- function(z_1, mu_1, sigma_1, initial_alloc_vec1){
    f <- function(y) (y*dlnorm(y, meanlog = mu_1, sdlog = sigma_1))
    h <- integrate(f, lower = 0, upper = z_1[1] - initial_alloc_vec1)$value
    
    g <- function(y) (dlnorm(y, meanlog = mu_1, sdlog = sigma_1))
    m <- integrate(g, lower = 0, upper = z_1[1] - initial_alloc_vec1)$value
    
    return(h/m)
} # Group 1

avg_pti2_prev_iter <- function(z_1, mu_2, sigma_2, initial_alloc_vec2){
    f <- function(y) (y*dlnorm(y, meanlog = mu_2, sdlog = sigma_2))
    h <- integrate(f, lower = 0, upper = z_1[1] - initial_alloc_vec2)$value
    
    g <- function(y) (dlnorm(y, meanlog = mu_2, sdlog = sigma_2))
    m <- integrate(g, lower = 0, upper = z_1[1] - initial_alloc_vec2)$value
    
    return(h/m)
} # Group 2

#------------------------------------------------------------------------#
# P0 - headcount ratio

P01_prev_iter <- function(z_1, mu_1, sigma_1, initial_alloc_vec1){
    g <- function(y) (dlnorm(y, meanlog = mu_1, sdlog = sigma_1))
    h <- integrate(g, lower = 0, upper = z_1[1] - initial_alloc_vec1)$value
    
    return(h)
} # Group 1

P02_prev_iter <- function(z_1, mu_2, sigma_2, initial_alloc_vec2){
    g <- function(y) (dlnorm(y, meanlog = mu_2, sdlog = sigma_2))
    h <- integrate(g, lower = 0, upper = z_1[1] - initial_alloc_vec2)$value
    
    return(h)
} # Group 2

#------------------------------------------------------------------------#

P0_prev_iter_mat <- function(){ 
    mats <- matrix(0, nrow = 2, ncol = 2)
    mats[1,1] <- P01_prev_iter(z_1, mu_1, sigma_1, initial_alloc_vec1)
    mats[2, 2] <- P02_prev_iter(z_1, mu_2, sigma_2, initial_alloc_vec2)
    return(mats)
} # Matrix A11 

#------------------------------------------------------------------------#

# Function to create poverty line
create_z <- function(hc_ratio){
    z = vector()
    for(i in 1:length(hc_ratio)){
        z[i] <- qlnorm(hc_ratio[i], meanlog = 11.0302, sdlog = sqrt(0.6688))
    }
    return(z)
}

z = create_z(hc_ratio = 0.2)
z_1 = c(z, z) # Vector of a repeated pre-defined national poverty line

#------------------------------------------------------------------------#

n = c(50, 50) # Vector of sample sizes for the two groups

w_prev_iter <- function(){
    C11 <- matrix(0, nrow = 2, ncol = 1)
    C11[1] <- (z_1[1]-avg_pti1_prev_iter(z_1[1], mu_1, sigma_1, initial_alloc_vec1))*
        P01_prev_iter(z_1[1], mu_1, sigma_1, initial_alloc_vec1)
    C11[2] <- (z_1[1]-avg_pti2_prev_iter(z_1[1], mu_2, sigma_2, initial_alloc_vec2))*
        P02_prev_iter(z_1[1], mu_2, sigma_2, initial_alloc_vec2)
    return(C11)
} # Matrix C11

zero <- c(0) # A single value - zero

#------------------------------------------------------------------------#

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

N = 100
fixed_budget <- c(0.05*N*fixed_budget_1) # Value of matrix C21 - set value

#------------------------------------------------------------------------#

# Initial values for allocation vector
initial_alloc_vec1 <- 0
initial_alloc_vec2 <- 0

criteria_met_right_side <- FALSE
previous_right_side <- c(10000000, 1000000, 10000000)
m <- 1

while(m<=100){
    column1 <- rbind(P0_prev_iter_mat(), n)
    column2 <- t(cbind(t(z_1), zero))
        
    A <- cbind(column1, column2) # Matrix A
        
    C <- rbind(w_prev_iter(), t(fixed_budget)) # Matrix C
        
    right_side <- MASS::ginv(A) %*% C

# Stopping criteria:
    if( abs(right_side[1] - previous_right_side[1]) < 0.0000001 & 
        abs(right_side[2] - previous_right_side[2]) < 0.0000001 &
        abs(right_side[3] - previous_right_side[3]) < 0.0000001) {
        criteria_met_right_side <- TRUE
        print(mu_1)
        print(sigma_1)
        print(mu_2)
        print(sigma_2)
        print("---------------")
        print(right_side)
        print("---------------")
        print("---------------")
        break;
        }
        initial_alloc_vec1 <- right_side[1]
        initial_alloc_vec2 <- right_side[2]
        previous_right_side = right_side
        m <- m+1
        print(mu_1)
        print(sigma_1)
        print(mu_2)
        print(sigma_2)
        print("---------------")
        print(right_side)
        print("---------------")
        print("---------------")
}



