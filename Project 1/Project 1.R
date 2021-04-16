library(ggplot2)

#####################################
# Problem A 1)
exp_generator = function(n=1, lambda){
  u = runif(n = n) # u = a vector of n uniformly distributed numbers between 0 and 1
  return(-1/lambda*log(u))
}
lambda = 2
n = 10000
exp_vector = data.frame(x = exp_generator(n, lambda))


ggplot(exp_vector, aes(x)) + 
  geom_histogram(aes(y = stat(density))) + 
  stat_function(fun = dexp, args = lambda)



#####################################
# 2) 

rg = function(n, alpha){
  u = runif(n = n)
  boundary = exp(1) / (alpha + exp(1))
  left_boundary = u < boundary
  right_boundary = !left_boundary
  
  u[left_boundary] = (u[left_boundary] / boundary) ** (1 / alpha)
  u[right_boundary] = -log((1-u[right_boundary]) / (boundary * alpha))
  return(u)
}

dg = function(x, alpha = 1) {
  c = alpha * exp(1) / (alpha + exp(1))
  d <- rep(0, length(x))
  left_indices <- 0 < x & x < 1
  right_indices <- 1 <= x
  d[left_indices] <- c * (x[left_indices] ** (alpha - 1))
  d[right_indices] <- c * exp(-x[right_indices])
  return(d)
}


#####################################
# Problem B
#####################################

# 1

gamma_samples = function(n = 1, alpha = 0.5, beta = 1){
  
  # define the acceptance probability 
  a = (1/alpha + 1/exp(1))/gamma(alpha)
  # Generate n samples from g
  #samp_g = rg(n, alpha, beta = beta)
  # Generating n uniform samples on (0, 1)
  u = runif(n)
  # initialise vector of samples
  samples = c()
  # number of succesful samples
  num_of_samp = 0
  repeat{
    # Generating n uniform samples on (0, 1)
    u = runif(n)
    # genereate samples from g as in A)
    x = rg(n = n, alpha =  alpha)
    # probability of accepting
    prob = dgamma(x, alpha, beta)/(a*dg(x, alpha))
    # appending the x-s where u < prob
    samples = c(samples, x[u < prob])
    # increas num_of_samp by the number of samples accepted
    num_of_samp = num_of_samp + length(x[u < prob])
    # Here we remove every accepted samples in addition to the n we want
    if(length(samples) > n){
      samples = samples[1:n]
      break
    }
  }
  return(samples)
}
####################
# We genereate many samples each iteration to save number of iterations 
# in the loop. We do this because vector operations are faster than loops 
# in R. The way we implemented it, the number of iterations in the for 
# loop is proportional to alpha.
####################

alpha = 0.5
smpls = gamma_samples(n = 100000, alpha = alpha)  
# to be able to use ggplot we transform the samples to a data frame
smpls = data.frame(samples = smpls) 
# making a nice plot for comparison
ggplot(data = smpls, aes(x = samples)) + 
  geom_histogram(data = smpls, 
                 aes(x=samples, y=..density..), 
                 binwidth = 0.02) + 
  stat_function(fun = dgamma, 
                geom = "line", 
                size = 1, 
                args = list(shape=alpha, rate=1), 
                aes(x=samples)) +
  ylim(0, 2) + 
  xlim(0, 4)

####################################################
# 2
####################################################

rum_gamma = function(n = 1, alpha = 2, beta = 1){
  # Calculating the limits of the sample region using the log
  log_a  = 1/2*(alpha-1)*(log(alpha-1) - 1)
  b1 = 0
  log_b2 = 1/2*(alpha+1)*(log(alpha+1) - 1)
  # creat an empty vector to store samples
  samples = c()
  tries = 0
  while(length(samples) < n){
    tries = tries + 1
    # Sampling from the uniform distribution on (0, 1) and 
    # transforming the samples to appropriate regions
    log_y1 = log_a + log(runif(1))
    log_y2 = log_b2 + log(runif(1))
    if((alpha+1)*log_y1 <= (alpha-1)*log_y2 - exp(log_y2 - log_y1)){
      y = exp(log_y2 - log_y1)
      samples = c(samples, y)
    }
  }
  return(c(tries, samples))
}
####################################################
# Plotting success rates
success_rates = function(alphas){
  result = c()
  n = 1000
  for(i in 1:length(alphas)){
    # Sample from the algorithm
    smpls_rum = rum_gamma(n = n, alpha = alphas[i])
    vect = c(alphas[i], smpls_rum[1]/n)
    result = cbind(result, vect)
  }
  result = data.frame(rates = t(result))
  colnames(result) = c("Alphas", "Success_rates")
  return(result)
}

alphas = c(2, 5, 30, 100, 200, 500, 1000, 1500, 2000)
df = success_rates(alphas)
ggplot(data = df) + 
  geom_col(data = df, 
           aes(x = Alphas, y = Success_rates), 
           width = 30)


###############################################
# plotting one realization
alpha = 50
n = 20000
# Sample from the algorithm
smpls_rum = rum_gamma(n = n, alpha = alpha)
tries = smpls_rum[1]
smpls_rum = smpls_rum[-1]
# to be able to use ggplot we transform the samples to a data frame
smpls_rum = data.frame(samples = smpls_rum) 

# making a nice plot for comparison
ggplot(data = smpls_rum, aes(x = samples)) + 
  geom_histogram(data = smpls_rum, 
                 aes(x=samples, y=..density..), 
                 binwidth = 1) + 
  stat_function(fun = dgamma, 
                geom = "line", 
                size = 1, 
                args = list(shape=alpha, rate=1), 
                aes(x=samples)) 

####################################################
# 3
####################################################
# The reccursion demands to much memory for large n
gam_from_exp = function(n = 1, alpha = 2){
  # sum of alpha exp(1) distributed variables using function from A 1)
  if(n>0){
    # concatinating current result with reccursive call
    result = c(sum(exp_generator(n = alpha, lambda = 1)), gam_from_exp(n = n-1, alpha))
    return(result)
  }
}
# Need to fix this function if large n is needed
# gam_from_exp = function(n = 1, alpha = 2){
#   result = c()
#   for(i in 1:n){
#     result = c(result, exp_generator(n = alpha, lambda = 1))
#   }
#   return(result)
# }
n = 10000
alpha = 5000
gam_realizations = gam_from_exp(n, alpha = alpha)
gam_realizations = data.frame(samples = gam_realizations)

# making a nice plot for comparison
ggplot(data = gam_realizations, aes(x = samples)) + 
  geom_histogram(data = gam_realizations, 
                 aes(x=samples, y=..density..)) + 
  stat_function(fun = dgamma, 
                geom = "line", 
                size = 1, 
                args = list(shape=alpha, rate=1), 
                aes(x=samples)) 


# ## 3 No longer included
# ### a)
# _Proof._
# $$
#   X_1 \sim \text{Gamma}(\alpha_1, \beta = 1)
# $$
#   and 
# $$
#   X_2 \sim \text{Gamma}(\alpha_2, \beta = 1)
# $$
#   are independent. By defining
# $$
#   X = X_1 + X_2
# $$
#   we get that
# $$
#   M_X(t) \overset{ind.}{=} M_{X_1}(t)M_{X_2}(t) \overset{\beta = 1}{=} \bigg{(}\frac{1}{1-t}\bigg{)}^{\alpha_1} \bigg{(}\frac{1}{1-t} \bigg{)}^{\alpha_2} = \bigg{(}\frac{1}{1-t} \bigg{)}^{\alpha_1+\alpha_2} \; \blacksquare
# $$
#   
#   
#   ### b)
#   ```{r}
# gam_from_exp = function(n = 1, alpha = 2){
#   # sum of alpha exp(1) distributed variables using function from A 1)
#   result = c()
#   for(i in n:1){
#     result = c(result, sum(exp_generator(n = alpha, lambda = 1)))
#   }
#   return(result)
# }
# 
# ```
# 
# #### Plotting
# We plot this to illustrate that it works
# ```{r}
# n = 10000
# alpha = 5000
# gam_realizations = gam_from_exp(n, alpha = alpha)
# gam_realizations = data.frame(samples = gam_realizations)
# 
# # making a nice plot for comparison
# ggplot(data = gam_realizations, aes(x = samples)) + 
#   geom_histogram(data = gam_realizations, 
#                  aes(x=samples, y=..density..)) + 
#   stat_function(fun = dgamma, 
#                 geom = "line", 
#                 size = 1, 
#                 args = list(shape=alpha, rate=1), 
#                 aes(x=samples)) 
# ```
# 


####################################################
# 4
####################################################

gam4 = function(n = 1, alpha = 2, beta = 1){
  beta*gam_from_exp(n, alpha)
}

n = 2000
alpha = 100
beta = 5

gam4_realizations = gam4(n, alpha, beta)
gam4_realizations = data.frame(samples = gam4_realizations)

# making a nice plot for comparison
ggplot(data = gam4_realizations, aes(x = samples)) + 
  geom_histogram(data = gam4_realizations, 
                 aes(x=samples, y=..density..)) + 
  stat_function(fun = dgamma, 
                geom = "line", 
                size = 1, 
                args = list(shape=alpha, rate=1/beta), 
                aes(x=samples)) 

####################################################
# 5
####################################################

beta_gen = function(n, alpha, beta){
  x = gam_from_exp(n, alpha)
  y = gam_from_exp(n, beta)
  return(x/(x+y))
}

n = 20000
alpha = 2
beta = 5
beta_vect = beta_gen(n, alpha, beta)
beta_vect = data.frame(samples = beta_vect)
# Plot
ggplot(data = beta_vect, aes(x = samples)) + 
  geom_histogram(data = beta_vect, 
                 aes(x=samples, y=..density..)) + 
  stat_function(fun = dbeta, 
                geom = "line", 
                size = 1, 
                args = list(shape1=alpha, shape2=beta), 
                aes(x=samples)) 

####################################################
# Probelm D
####################################################
# 1












