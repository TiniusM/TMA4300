library(boot)
library(ggplot2)
library(invgamma)

#######################################
# Problem A1
# Here we we start with the first year with 0 accidents and count one more accident per row in the data.frame
# The last year in the data do not represent another accident. Hence the two last elements in cumulative number of accidents are equal
cum_num = c(c(0:189), 189)
cumulative_plot = ggplot(data = coal, aes(x = date, y = cum_num)) + 
  geom_path() 
cumulative_plot
# The plot shows a steep linear trend until around 1988, then the growth of the cumulative decreased. 
# i.e. time between accidents increased, it became safer.



####################################################################################
####################################################################################
get_y = function(x, t){
  # the number of observations from t_k to t_{k+1}, but we only have two intervals, hence
  y_0 = sum(x <= t)
  y_1 = length(x) - y_0
  return(list("y_0" = y_0, "y_1" = y_1))
}

log_fc_t_1 = function(x, t){
  return(get_y(x, t)$y_0*log(lambda_0) + get_y(x, t)$y_1*log(lambda_1) + t*(lambda_1-lambda_0))
}

single_site_mcmc = function(t_1, lambda_0, lambda_1, beta, x, sigma, n = 10000){
  # We initialise the data frame where we will store the values
  data_mcmc =  data.frame(matrix(nrow = n, ncol = 4))
  colnames(data_mcmc) = c("t_1", "lambda_0", "lambda_1", "beta")
  data_mcmc[1, ] = c(t_1, lambda_0, lambda_1, beta)
  
  # We calculate the log full conditional
  log_fc_new = log_fc_t_1(x, t_1)
  
  success_count = 0
  for (i in 2:n){
    # Generate proposals
    t_1_new = rnorm(1, mean = t_1, sd = sigma)
    
    # if the proposal is outside the interval [t_0, t_2]
    if((t_1_new < t_0) | (t_1_new > t_2)){
      log_accept = -Inf # corresponds to alpha = 0 on non-log scale
    }
    else{
      # We store the log full conditional of the last iteration
      log_fc = log_fc_new
      
      # We get the y values for the new t_1
      y_new = get_y(x, t_1_new)
      
      # We calculate the log of the full conditional for t_1_new
      log_fc_new = log_fc_t_1(x, t_1_new)
      # Calculate the acceptance probability on log scale
      log_accept = log_fc_new - log_fc
    }
    
    # Accept or reject
    u <- runif(1)
    # If u < alpha, we accept the proposal
    if(u < exp(log_accept)){
      t_1 = t_1_new
      success_count = success_count + 1
    }
    # else{ do nothing }
    
    # We add the selected values in the data frame
    data_mcmc[i, ] = c(t_1, lambda_0, lambda_1, beta)
    
    # Full conditional for the parameters
    beta = rinvgamma(1, shape = 4, rate = (lambda_0 + lambda_1 + 1)) 
    lambda_0 = rgamma(1, shape = get_y(x, t_1)$y_0 + 2, rate = (t_1 - t_0 + 1/beta))
    lambda_1 = rgamma(1, shape = get_y(x, t_1)$y_1 + 2, rate = (t_2 - t_1 + 1/beta))
    
  }
  return(list("chain" = data_mcmc, "success_rate" = success_count/n))
}


# These are global variables because we assume n = 1
t_0 = coal[1,]
t_2 = coal[dim(coal)[1], ]
# Remove the first and last element of the coal data frame
x = coal[-1, ]
x = x[-length(x)]

####################################
# Initialising values
t_init = (t_0 + t_2)/2
beta = 1

lambda_0 = rgamma(1, shape = get_y(x, t = t_init)$y_0 + 2, rate = (t_init - t_0 + 1/beta))
lambda_1 = rgamma(1, shape = get_y(x, t = t_init)$y_1 + 2, rate = (t_2 - t_init + 1/beta))
lambda_0;lambda_1

# Now we test for different sigmas
# first for sigma = 1
sigma1 = 1
mcmc_relisation1 = single_site_mcmc(t_init, lambda_0, lambda_1, beta, x, sigma = sigma1, n = 5000)

# sigma = 4
sigma2 = 4
mcmc_relisation2 = single_site_mcmc(t_init, lambda_0, lambda_1, beta, x, sigma = sigma2, n = 5000)

# sigma = 15
sigma3 = 15
mcmc_relisation3 = single_site_mcmc(t_init, lambda_0, lambda_1, beta, x, sigma = sigma3, n = 5000)


par(mfrow = c(3, 1))
plot(mcmc_relisation1$chain$t_1, type = 'l', ylab = expression(t[1]), xlab = "Iterations")
plot(mcmc_relisation2$chain$t_1, type = 'l', ylab = expression(t[1]), xlab = "Iterations")
plot(mcmc_relisation3$chain$t_1, type = 'l', ylab = expression(t[1]), xlab = "Iterations")

mcmc_relisation1$success_rate; mcmc_relisation2$success_rate; mcmc_relisation3$success_rate

par(mfrow = c(3, 1))
hist(mcmc_relisation1$chain$t_1)
hist(mcmc_relisation2$chain$t_1)
hist(mcmc_relisation3$chain$t_1)

slope_0 = mean(mcmc_relisation2$chain$lambda_0)
slope_1 = mean(mcmc_relisation2$chain$lambda_1)
t_1_mean = mean(mcmc_relisation2$chain$t_1)

# To make the plots
slope_0_df = data.frame(list('x' = seq(t_0, t_1_mean), 'y' = slope_0*(seq(t_0,t_1_mean)-t_0)))
slope_1_df = data.frame(list('x' = seq(t_1_mean, t_2), 'y' = slope_1*(seq(t_1_mean,t_2)-t_1_mean) + get_y(x, t_1_mean)$y_0))
cumulative_plot + geom_vline(xintercept = t_1_mean) +
  geom_path(data = slope_0_df, aes(x = x, y = y, col = 'slope = lambda_0')) +
  geom_path(data = slope_1_df, aes(x = x, y = y, col = 'slope = lambda_1'))

get_y(x, t_2)
length(x)
  
