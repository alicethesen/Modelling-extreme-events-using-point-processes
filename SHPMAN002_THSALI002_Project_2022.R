## Load required packages

library(ggplot2)
library(evmix)
library(POT)
library(gridExtra)
library(ggpubr)
library(moments)


### Data  =================================================================================
# read in data
sav <- read.csv("BenSav_L1.csv", sep = ";")

# Making date-time variable
sav$date_time <- as.POSIXct(sav$date, format = "%m/%d/%Y %H:%M")

# Make just date variable 
sav$just_date <- as.POSIXct(sav$date, format = "%m/%d/%Y")

# Make just time variable
sav$just_time <- format(sav$date_time, format = "%H:%M")

###### Important varibles included in functions =================================================================================

# pars: optimized parameters for a particular model
# pars consists of lambda, alpha and beta for the Hawkes models
# pars consists of lambda, alpha, beta and gamma for the marked Hawkes models

# H_t: contains the indicies of the histroy 
# t: a sequence from 1 to the last time index


###### Hawkes functions =================================================================================


hawkes_intensity <- function(pars, H_t, t){
  ## Hawkes intensity with Exponential Decay
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  
  # Actual Hawkes intensity vector
  hawkes_vec <- c()
  
  # Vectors that will be used to plot the intensity over time
  hawkes_plot <- c()
  plot_index <- c()
  

  for (s in  t){
    # iterate over the entire time interval
    
    if (any(H_t <= t[s])){
      # If an event has occurred before time s
      
      ind <- which(H_t <= t[s]) # Indices of events that have occurred before time s
      
      lambda_x <-   lambda + alpha*sum(exp(-beta*(t[s] - H_t[ind]) ) )
      hawkes_vec <- append(hawkes_vec, lambda_x)
      
      if (any(H_t == t[s])){
        # determine if an event occurred at time s
        
        event <- which(H_t == t[s])
        
        if(event == 1){
          # If the event arrival at time s is the first event
          
          hawkes_plot <- append(hawkes_plot, lambda)
          plot_index <- append(plot_index, s)          
          
        } else {
          # If the event arrival at time s is not the first event
          
          ind_new <- ind[1:(length(ind)-1)]
          
          plot_x <-   lambda + alpha*sum(exp(-beta*(t[s] -H_t[ind_new]) ) )
          hawkes_plot <- append(hawkes_plot, plot_x)
          plot_index <- append(plot_index, s)
        }
      } 
      
      hawkes_plot <- append(hawkes_plot, lambda_x)
      plot_index <- append(plot_index, s)
      
    } else {
      # If no event has occurred yet 
      
      lambda_x <- lambda 
      hawkes_vec <- append(hawkes_vec, lambda_x)
      hawkes_plot <- append(hawkes_plot, lambda_x)
      plot_index <- append(plot_index, s)
    }
  }
  return(list(hawkes_vec = hawkes_vec, hawkes_plot = hawkes_plot, plot_index = plot_index))
}

hawkes_inst_cond_intensity <- function(pars, H_t, time){
  # Hawkes instantaneous conditional intensity 
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  
  # Get indices of events that have happened previously to the current time index
  H_tt = H_t[H_t < time]
  
  if(length(H_tt) == 0){ 
    # If no event has happened 
    ICI = lambda
  }
  
  else{ 
    # If there has been at least one event iin the history
    ICI = lambda + sum( alpha * exp( - beta * ( time - H_tt)))
  }
  return(ICI)
}

hawkes_compensator <- function(pars, H_t, t){
  # Hawkes compensator 
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  
  comp_vec <- c()
  
  for (s  in  t){
    # iterate over the entire time interval
    
    if (any(H_t <= t[s])){
      # If an event has occurred before time s
      
      comp <- lambda * s
      
      ind <- which(H_t <= t[s]) # Indices of events that have occurred before time s
      
      comp <- comp + (alpha/beta)*sum((1 - exp(-beta *(t[s] - H_t[ind] ) ) ) )
      
      comp_vec <- append(comp_vec, comp)
      
    } else {
      # If no event has occurred before time s
      comp_vec <- append(comp_vec, 0)
    }
  }
  comp_df <- data.frame(comp_vec, t)
  return(comp_df)
}

hawkes_log_like <- function(pars, H_t, t){
  # Log-Likelihood of Hawkes with Exponential Decay
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  
  LT <- t[length(t)] # last time index 
  N <- length(H_t) # number of events in the history 
  A <- rep(0, N) # A(1) = 0
  
  hawkes_vec <- c()
  
  for (i in 2:N){
    # Recursive term 
    
    A[i] <- exp(-beta * (H_t[i] - H_t[i - 1])) * (1 + A[i - 1])
  }
  
  log_like <- -lambda * LT
  
  for (i in 1:N){
    # sum over all the events in the history 
    
    log_like <- log_like + log(lambda + alpha * A[i]) - (alpha/beta) * (1 - exp(-beta*(LT - H_t[i])))
    
  }
  return(-log_like)  
}

###### Marked Hawkes functions =================================================================================


marked_intensity <- function(pars, H_t, t, m){
  # Marked Hawkes with Exponential Decay and exponential impact function
  # pars, history, time_index, marks
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  gamma <- pars[4]
  
  # Actual Hawkes intensity vector
  hawkes_vec <- c()
  
  # Vectors that will be used to plot the intensity over time
  hawkes_plot <- c()
  plot_index <- c()
  
  for (s in t){
    # iterate over the entire time interval
    
    if (any(H_t <= t[s])){
      
      ind <- which(H_t <= t[s])
      
      lambda_x <-   lambda + alpha*sum(exp(gamma*m[ind]-beta*(t[s] - H_t[ind]) ) )
      hawkes_vec <- append(hawkes_vec, lambda_x)
      
      if (any(H_t == t[s])){
        # determine if an event occurred at time s
        
        event <- which(H_t== t[s])
        
        if(event == 1){
          # If the event arrival at time s is the first event
          
          hawkes_plot <- append(hawkes_plot, lambda)
          plot_index <- append(plot_index, s)          
          
        } else {
          # If the event arrival at time s is not the first event
          
          ind_new <- ind[1:(length(ind)-1)]
          
          plot_x <-   lambda + alpha*sum(exp(gamma*m[ind_new]-beta*(t[s] - H_t[ind_new]) ) )
          hawkes_plot <- append(hawkes_plot, plot_x)
          plot_index <- append(plot_index, s)
        }
      }
      
      hawkes_plot <- append(hawkes_plot, lambda_x)
      plot_index <- append(plot_index, s)
      
    } else {
      # If no event has occurred yet 
      
      lambda_x <- lambda 
      hawkes_vec <- append(hawkes_vec, lambda_x)
      hawkes_plot <- append(hawkes_plot, lambda_x)
      plot_index <- append(plot_index, s)
    }
    
  }
  return(list(hawkes_vec = hawkes_vec, hawkes_plot = hawkes_plot, plot_index = plot_index))
  
}

marked_inst_cond_intensity <- function(pars, H_t, m, time){
  # Marked Hawkes instantaneous conditional intensity 
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  gamma <- pars[4]
  
  # Get indices of events that have happened previously to the current time index
  H_tt = H_t[H_t < time]
  
  if(length(H_tt) == 0){ 
    # If no event has happened 
    ICI = lambda
  }
  
  else{ 
    # If there has been at least one event iin the history
    ICI = lambda + sum( alpha * exp(gamma * m - beta * ( time - H_tt)))
  }
  return(ICI)
}

marked_compensator <- function(pars, H_t, t, m){
  # marked Hawkes compensator 
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  gamma <- pars[4]
  
  comp_vec <- c()
  
  for (s  in  t){
    # iterate over the entire time interval
    
    if (any(H_t <= t[s])){
      # If an event has occurred before time s
      
      comp <- lambda * s
      
      ind <- which(H_t <= t[s]) # Indices of events that have occurred before time s
      
      comp <- comp + (alpha/beta)*sum(exp(gamma*m[ind])*(1 - exp(-beta *(t[s] - H_t[ind] ) ) ) )
      
      comp_vec <- append(comp_vec, comp)
      
    } else {
      # If no event has occurred before time s
      comp_vec <- append(comp_vec, 0)
    }
  }  
  comp_df <- data.frame(comp_vec, t)
  return(comp_df)
}

marked_log_like <- function(pars, H_t, t, m){
  # Log-Likelihood of Marked Hawkes with Exponential Decay and exponential impact function
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  gamma <- pars[4]
  
  # Rate of marks distribution
  rate <- 1/mean(m)
  
  LT <- t[length(t)] # last time index 
  N <- length(H_t) # number of events in the history 
  A <- rep(lambda, N) # A(1) = lambda
  
  hawkes_vec <- c()
  
  for (i in 2:N){
    # Recursive term 
    
    A[i] <- lambda + exp(-beta * (H_t[i] - H_t[i - 1])) * (A[i - 1]- lambda) + alpha * exp(-beta * (H_t[i] - H_t[i - 1]))*exp(gamma*m[i-1])
  }
  
  log_like <- -lambda * LT + N * log(rate)
  
  for (i in 1:N){
    # sum over all the events in the history 
    
    log_like <- log_like + log(A[i]) - (alpha/beta) *exp(gamma*m[i]) *(1 - exp(-beta*(LT - H_t[i]))) - rate*m[i]
    
  }
  
  return(-log_like)  
}

###### AIC/BIC function =================================================================================

AIC.BIC <- function(pars, H_t, log_likelihood){
  # Find the AIC and BIC for model comparison 
  
  AIC <- 2 * log_likelihood + 2 * length(pars)
  BIC <- 2 * log_likelihood + length(pars) * length(H_t)
  return(list(log_likelihood = log_likelihood, AIC = AIC, BIC = BIC))
}

###### Counting process functions =================================================================================

counting_process <- function(H_t, t){
  # Counting procss N(t)
  
  count_vec <- c()
  count_num <- 0 
  
  for (s  in  t){
    # iterate over the entire time interval
    
    if (any(H_t <= t[s])){
      # If an event has occurred before time s
      
      ind <- which(H_t <= t[s]) # Indices of events that have occurred before time s
      
      if (any((t[s]- 0.5) < H_t & H_t <= (t[s] + 0.5))){ 
        # If an events is 0.5 on either side of the current time index
        # 0.5 is used as the time interval is continuous even though t is discrete
        
        number <- length(which((t[s]- 0.5) < H_t & H_t <= (t[s] + 0.5)))
        count_num <- count_num + number
      }
      
      count_vec <- append(count_vec, count_num)
      
    } else {
      # If no event has occurred yet then count will be equal to zero 
      count_vec <- append(count_vec, count_num)
    }
  }
  count_df <- data.frame(count_vec, t)
  return(count_df)
}

#### Moments =================================================================================

hawkes_intensity_expect <- function(pars, t, intensity_vec){
  # Hawkes expectation conditional intensity given initial intensity  
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  
  lambda_inital <- intensity_vec[1] # lambda at time 1
  K <- beta - alpha 
  
  m1 <- c()
  
  for (s in 1:length(t)){
    # iterate over the entire time interval
    
    m1[s] <- (lambda*beta)/K + (lambda_inital - ((lambda*beta)/K)) * exp(-K * t[s])
  }
  return(m1)
}

hawkes_count_expect <- function(pars, t, intensity_vec){
  # Hawkes expectation of counting process 
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  
  lambda_inital <- intensity_vec[1] # lambda at time 1
  K <- beta - alpha
  
  # get the first moment of the conditional intensity 
  m1 <- hawkes_intensity_expect(pars, t, intensity_vec)
  w1 <- c()
  
  for (s in 1:length(t)){
    # iterate over the entire time interval
    w1[s] <- m1[s]* t[s] + (lambda_inital - m1[s])* (1/K) * (1 - exp(-K * t[s]))
  }
  return(w1)
}


marked_intensity_expect <- function(pars, t, intensity_vec, m){
  # marked Hawkes expectation conditional intensity given initial intensity
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  gamma <- pars[4]
  
  lambda_inital <- intensity_vec[1] # lambda at time 1
  
  rate <- 1/mean(m)
  
  par_marks <- rate/(rate - gamma)
  
  K <- beta - alpha*par_marks
  
  m1 <- c()
  for (s in 1:length(t)){
    # iterate over the entire time interval
    m1[s] <- (lambda*beta)/K + (lambda_inital - ((lambda*beta)/K)) * exp(-K * t[s])
  }
  return(m1)
}



marked_count_expect <- function(pars, t, intensity_vec, m){
  # marked Hawkes expectation of counting process
  
  lambda <- pars[1]
  alpha <- pars[2]
  beta <- pars[3]
  gamma <- pars[4]
  
  lambda_inital <- intensity_vec[1] # lambda at time 1
  
  rate <- 1/mean(m)
  
  # Moment generating function 
  par_marks <- rate/(rate - gamma)
  
  K <- beta - alpha*par_marks
  
  # get the first moment of the conditional intensity 
  m1 <- marked_intensity_expect(pars, t, intensity_vec, m)
  w1 <- c()
  
  for (s in 1:length(t)){
    # iterate over the entire time interval
    w1[s] <- m1[s]* t[s] + (lambda_inital - m1[s])* (1/K) * (1 - exp(-K * t[s]))
  }
  return(w1)
}



###### Simulation by thinning functions =================================================================================

hawkes_sim_thinning = function(pars, t, kappa){
  # Simulation Hawkes with Exponential Decay By Thinning
  
  LT <- t[length(t)]
  
  time = 0
  H_t = c()
  
  while(time < LT){
    # until time is equal to the last time index 
    
    # Get Hawkes instantaneous conditional intensity at time t + 1e-04
    Mt = hawkes_inst_cond_intensity(pars, H_t, time + 1e-04)
    Lt = kappa * Mt
    R = rexp(1, Mt)
    
    if( R > Lt){ 
      # no event has occurred 
      time = time + Lt
    }
    
    else{
      
      # Get Hawkes instantaneous conditional intensity at time t + R
      cond = hawkes_inst_cond_intensity(pars, H_t, (time + R))/Mt
      
      # Random uniform for whether event has occurred or not 
      U = runif(1, min=0, max=1)
      
      if( U[1] > cond[1]) {
        # no event has occurred 
        time = time + R
        
      } 
      else {
        # event has occurred 
        time = time + R
        H_t = c(H_t, time)
        
      }
    }
  }
  
  H_t = H_t[ H_t < LT ]
  
  return(list(H_t = H_t))
}

marked_sim_thinning = function(pars, t, mark_rate, kappa){
  # Simulation Marked Hawkes with exponential decay and impact By Thinning
  
  LT <- t[length(t)]
  
  time = 0
  H_t = c()
  X_t = c()
  
  while(time < LT){
    # until time is equal to the last time index 
    
    # Get marked Hawkes instantaneous conditional intensity at time t + 1e-04
    Mt = marked_inst_cond_intensity(pars, H_t, X_t, time + 1e-04)
    Lt = kappa * Mt
    R = rexp(1,  Mt)
    
    if (R > Lt){ 
      # no event has occurred 
      time = time + Lt
    }
    
    else{
      # Get marked Hawkes instantaneous conditional intensity at time t + R
      cond = marked_inst_cond_intensity(pars, H_t, X_t, (time + R))/Mt
      
      # Random uniform for whether event has occurred or not 
      U = runif( 1, min=0, max=1)
      
      if( U[1] > cond[1]) {
        # no event has occurred 
        time = time + R
        
      } 
      else {
        # event has occurred - obtain a mark for that event 
        M = rexp(1, rate = mark_rate)
        time = time + R
        X_t = c( X_t, M)
        H_t = c( H_t, time)
        
      }
    }
  }
  
  H_t = H_t[ H_t < LT ]
  X_t = X_t[1:length(H_t)]
  
  return(list(H_t = H_t, X_t = X_t))
}

### Hawkes QQ-plots from simulating functions =================================================================================

Hawkes_QQ <- function(N, IAT, pars, t){
  
  pp <- ppoints(100)
  obs_quant <- quantile(sort(IAT), pp)
  
  # Sample Quantile Matrix
  quant_mat <- matrix(0, nrow = N, ncol = 100)
  sim_history <- list()
  
  for (i in 1:N){
    
    print("Simulations")
    print(i)
    
    # Use simulation by thinning to get simulated event history 
    simulation <-  hawkes_sim_thinning(pars = pars, t, kappa = 1)
    
    # Get simulated event history inter-arrival times 
    sim_IAT <- (simulation$H_t[2:length(simulation$H_t)] - simulation$H_t[1:(length(simulation$H_t)-1)])+1
    
    quant_mat[i,] <- quantile(sim_IAT, pp)
    
    # Add histroy 
    sim_history <- append(sim_history, simulation)
    
  }
  
  # Aggregating Quantiles Across Samples
  sim_2.5  <- rep(0, 100)
  sim_mean <- rep(0, 100)
  sim_97.5 <- rep(0, 100)
  
  for (k in 1:100){
    # Find mean and 95% confidence interval 
    
    sim_2.5[k] <- quantile(quant_mat[,k], 0.025)
    sim_mean[k] <- mean(quant_mat[,k])
    sim_97.5[k] <- quantile(quant_mat[,k], 0.975)
  }
  
  # Making Data Frames
  
  qq_sim_2.5 <- as.data.frame(cbind(sim_mean, sim_2.5))
  qq_sim_mean <- as.data.frame(cbind(obs_quant, sim_mean))
  qq_sim_97.5 <- as.data.frame(cbind(sim_mean, sim_97.5))
  
  return(list(Sim_2.5 = qq_sim_2.5, Sim_mean = qq_sim_mean, Sim_97.5 = qq_sim_97.5, History = sim_history))
  
}


### Marked Hawkes QQ-plots from simulating =================================================================================

Marked_QQ <- function(N, IAT, pars, t, m){
  
  pp <- ppoints(100)
  
  obs_quant <- quantile(sort(IAT), pp)
  
  rate <- 1/mean(m)
  
  # Making qqline:
  
  # Sample Quantile Matrix
  quant_mat <- matrix(0, nrow = N, ncol = 100)
  sim_history <- list()
  
  for (i in 1:N){
    
    print("Simulation")
    print(i)
    
    # Use simulation by thinning to get simulated event history 
    simulation <-  marked_sim_thinning(pars = pars, t =  t, rate, kappa = 1 )
    
    # Get simulated event history inter-arrival times 
    sim_IAT <- (simulation$H_t[2:length(simulation$H_t)] - simulation$H_t[1:(length(simulation$H_t)-1)])+1
    quant_mat[i,]  <- quantile(sim_IAT, pp)
    
    # Add histroy
    sim_history <- append(sim_history, simulation)
  }
  
  # Aggregating Quantiles Across Samples
  sim_2.5  <- rep(0, 100)
  sim_mean <- rep(0, 100)
  sim_97.5 <- rep(0, 100)
  
  for (k in 1:100){
    # Find mean and 95% confidence interval 
    
    sim_2.5[k] <- quantile(quant_mat[,k], 0.025)
    sim_mean[k] <- mean(quant_mat[,k])
    sim_97.5[k] <- quantile(quant_mat[,k], 0.975)
    
  }
  
  # Making Data Frames
  
  qq_sim_2.5 <- as.data.frame(cbind(sim_mean, sim_2.5))
  qq_sim_mean <- as.data.frame(cbind(obs_quant, sim_mean))
  qq_sim_97.5 <- as.data.frame(cbind(sim_mean, sim_97.5))
  
  return(list(Sim_2.5 = qq_sim_2.5, Sim_mean = qq_sim_mean, Sim_97.5 = qq_sim_97.5, History = sim_history))
}


###### Optim grid search =================================================================================

hawkes_optim <- function(H_t, time_index){
  # Hawkes grid search for finding max log-likelihood 
  
  # Create sets of search values 
  hawkes_grid_search <- expand.grid(seq(0.01, 0.1, length.out = 15), 
                                    seq(0.01, 0.3, length.out = 15), seq(0.01, 0.3, length.out = 15))
  
  lambda_seq <- hawkes_grid_search$Var1
  alpha_seq <- hawkes_grid_search$Var2
  beta_seq <- hawkes_grid_search$Var3
  
  search_mat <- matrix(0, nrow = length(hawkes_grid_search$Var1), ncol = 4)
  
  for (i in 1:length(hawkes_grid_search$Var1)){
    # Iterate through all of the parameter initial value sets
    
    search_optim <- optim(par = c(lambda_seq[i], alpha_seq[i], beta_seq[i]), fn = hawkes_log_like, 
                          method = "L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001),
                          upper = c(Inf, Inf, Inf), H_t = H_t, t = time_index)
    
    # store optimal paramater values and corresponding log-likelihood 
    search_out <- c(search_optim$par, search_optim$value)
    
    search_mat[i,] <- search_out
    
    print(i)
  }
  
  return(search_mat)
}

marked_optim <- function(H_t, time_index, m){
  # marked Hawkes grid search for finding max log-likelihood 
  
  # Create sets of search values 
  marked_grid_search  <- expand.grid(seq(0.001, 0.01, length.out = 8), seq(0.1, 0.3, length.out = 8), 
                                     seq(0.2, 0.4, length.out = 8), seq(0.001, 0.01, length.out = 7))
  
  lambda_seq <- marked_grid_search$Var1
  alpha_seq <- marked_grid_search$Var2
  beta_seq <- marked_grid_search$Var3
  gamma_seq <- marked_grid_search$Var4
  
  search_mat <- matrix(0, nrow = length(marked_grid_search$Var1), ncol = 5)
  
  for (i in 1:length(marked_grid_search$Var1)){
    # Iterate through all of the parameter initial value sets
    
    search_optim <- try(optim(par = c(lambda_seq[i], alpha_seq[i], beta_seq[i], gamma_seq[i]), fn = marked_log_like, 
                              method = "L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001),
                              upper = c(Inf, Inf, Inf), H_t = H_t, t = time_index, m = m))
    
    if (class(search_optim) != "try-error"){
      # If there is an error, still continue running function 
      
      # store optimal paramater values and corresponding log-likelihood 
      search_out <- c(search_optim$par, search_optim$value)
      
      search_mat[i,] <- search_out
    }
    
    print(i)
  }
  
  return(search_mat)
  
}


### Exploring CO2 variable =================================================================================

# remove NAs from CO2 data 
CO2_NA_omit <- na.omit(sav$CO2)

# 5 number summary and other numerical summaries
round(fivenum(CO2_NA_omit), 3) # 351.904   403.015   409.863   419.528   549.300
mean(CO2_NA_omit) # 412.6652
var(CO2_NA_omit) # 252.5086
sd(CO2_NA_omit) # 15.89052
skewness(CO2_NA_omit) # 2.010641 
kurtosis(CO2_NA_omit) # 11.76855 

# CO2 histogram 
ggplot(sav, aes(x = CO2)) + 
  geom_histogram(color = "black", fill = "seagreen", bins = 50, na.rm= TRUE, lwd = 0.4) + 
  ylab("Count") + xlab(expression(paste("CO"[2],"  (mmol per mol)",sep='')))


### Splitting into wet and dry season =================================================================================

# first period of wet (incomplete)
start_wet1 <- which(sav$date_time == as.POSIXct("1/14/2020 14:30", format = "%m/%d/%Y %H:%M"))
end_wet1 <- which(sav$date_time == as.POSIXct("04/15/2020 23:30", format = "%m/%d/%Y %H:%M"))

# index of first period of wet 
wet_index1 <- seq(start_wet1, end_wet1, 1)


# First dry (complete)
start_dry1 <- which(sav$date_time == as.POSIXct("04/16/2020 0:00", format = "%m/%d/%Y %H:%M"))
end_dry1 <- which(sav$date_time == as.POSIXct("10/19/2020 23:30", format = "%m/%d/%Y %H:%M") )

# index of first dry
dry_index1 <- seq(start_dry1, end_dry1, 1)


# Second wet (complete)
start_wet2 <- which(sav$date_time == as.POSIXct("10/20/2020 0:00", format = "%m/%d/%Y %H:%M"))
end_wet2 <- which(sav$date_time == as.POSIXct("04/15/2021 23:30", format = "%m/%d/%Y %H:%M"))

# index of second wet 
wet_index2 <- seq(start_wet2, end_wet2, 1)


# Second dry (complete)
start_dry2 <- which(sav$date_time == as.POSIXct("04/16/2021 0:00", format = "%m/%d/%Y %H:%M"))
end_dry2 <- which(sav$date_time == as.POSIXct("10/19/2021 23:30", format = "%m/%d/%Y %H:%M") )

# index of second dry 
dry_index2 <- seq(start_dry2, end_dry2, 1)


#  Third wet (incomplete)
start_wet3 <- which(sav$date_time == as.POSIXct("10/20/2021 0:00", format = "%m/%d/%Y %H:%M"))
end_wet3 <- which(sav$date_time == as.POSIXct("1/2/2022 18:00", format = "%m/%d/%Y %H:%M"))

# index of third wet
wet_index3 <- seq(start_wet3, end_wet3, 1)


# Create sub set for each season
wet_season1 <- sav[wet_index1,]
dry_season1 <- sav[dry_index1,]
wet_season2 <- sav[wet_index2,]
dry_season2 <- sav[dry_index2,]
wet_season3 <- sav[wet_index3,]


# Making Categorical Season Variable
# Only Wet vs Dry season 
season_overall <- vector("character", 34520)
season_overall[wet_index1] <- "Wet"
season_overall[wet_index2] <- "Wet"
season_overall[wet_index3] <- "Wet"
season_overall[dry_index1] <- "Dry"
season_overall[dry_index2] <- "Dry"

# Specific Wet vs Dry season 
season_specific <- vector("character", 34520)
season_specific[wet_index1] <- "Wet1"
season_specific[wet_index2] <- "Wet2"
season_specific[wet_index3] <- "Wet3"
season_specific[dry_index1] <- "Dry1"
season_specific[dry_index2] <- "Dry2"

# appending categorical season to dataset
sav$season_overall <- season_overall
sav$season_specific <- season_specific


### Exploring seasons =================================================================================

# Average CO2 for each season
round(mean(wet_season1$CO2, na.rm = T),3)
round(mean(dry_season1$CO2, na.rm = T),3)
round(mean(wet_season2$CO2, na.rm = T),3)
round(mean(dry_season2$CO2, na.rm = T),3)
round(mean(wet_season2$CO2, na.rm = T),3)


# CO2 separated by seasons 
# create dataframe
sav.df <- sav
sav.df$season_overall_wet <- sav.df$CO2
sav.df$season_overall_dry <- sav.df$CO2

# Give NA values when it is the dry season 
sav.df$season_overall_wet[which(sav.df$season_overall == "Dry")] <- NA

# Give NA values when it is the wet season 
sav.df$season_overall_dry[which(sav.df$season_overall == "Wet")] <- NA

# Create CO2 separated by seasons plot
ggplot(data=sav.df, aes(x= date_time)) +
  geom_line(aes(y = season_overall_wet, colour = "Wet"), na.rm = TRUE)+
  geom_line(aes(y = season_overall_dry, colour = "Dry"), na.rm = TRUE)+
  geom_smooth(aes(y = CO2),na.rm = TRUE, alpha = 0.4, col = "black")+
  labs(x ="Date", y = expression(paste("CO"[2],"  (mmol per mol)",sep=''))) +
  scale_colour_manual(name = 'Season', values=c(Wet = "lightseagreen", Dry = "firebrick3"))


# Plotting CO2 Historgrams for each season 

## Dry 1

# 5 number summary and other numerical summaries 
round(fivenum(dry_season1$CO2),3) # 351.904   401.918   407.977   413.380   546.168
mean(na.omit(dry_season1$CO2)) # 407.9969
sd(na.omit(dry_season1$CO2)) # 9.831803
skewness(na.omit(dry_season1$CO2))  # 0.7974247
kurtosis(na.omit(dry_season1$CO2)) # 9.596339 

# histogram 
dry1 <- ggplot(dry_season1, aes(x = CO2)) + 
  geom_histogram(color = "black", fill = "firebrick3", bins = 40, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + ylab("Count") + 
  ggtitle("Dry Season 1")+ xlab(expression(paste("CO"[2],"  (mmol per mol)",sep='')))


## Wet 2 (complete)

# 5 number summary and other numerical summaries 
round(fivenum(wet_season2$CO2),3) # 361.511   411.670   417.370   426.892   549.300
mean(na.omit(wet_season2$CO2)) # 422.3847
sd(na.omit(wet_season2$CO2)) # 18.08004
skewness(na.omit(wet_season2$CO2))  #  2.631455
kurtosis(na.omit(wet_season2$CO2)) # 13.23738

# histogram
wet2 <- ggplot(wet_season2, aes(x = CO2)) + 
  geom_histogram(color = "black", fill = "lightseagreen", bins = 40, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + ylab("") + 
  ggtitle("Wet Season")+ xlab(expression(paste("CO"[2],"  (mmol per mol)",sep='')))


# Dry 2

# 5 number summary and other numerical summaries 
round(fivenum(dry_season2$CO2),3) # 357.615   398.590   404.054   409.041   545.562
mean(na.omit(dry_season2$CO2)) # 404.1417
sd(na.omit(dry_season2$CO2)) # 9.453532
skewness(na.omit(dry_season2$CO2))  # 3.144173
kurtosis(na.omit(dry_season2$CO2)) # 36.33641

# histogram
dry2 <- ggplot(dry_season2, aes(x = CO2)) + 
  geom_histogram(color = "black", fill = "firebrick3", bins = 40, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + ylab("") + 
  ggtitle("Dry Season 2")+ xlab(expression(paste("CO"[2],"  (mmol per mol)",sep='')))

# All histogram plots together 
ggarrange(dry1, wet2, dry2, ncol=3, common.legend = F)


# boxplot by season (dry 1, wet, dry 2):
complete_seasons_index <- seq(start_dry1, end_dry2, 1)
ggplot(sav[complete_seasons_index,], aes(y = CO2, x = season_specific, fill = season_specific)) + 
  geom_boxplot(na.rm= TRUE, alpha = 0.8)+
  scale_fill_manual(values=c("firebrick3","firebrick3","lightseagreen"))+  xlab("Season")+
  scale_x_discrete(limits=c("Dry1", "Wet2", "Dry2"), labels = c("Dry 1", "Wet", "Dry 2"))+
  guides(fill = "none")+ ylab(expression(paste("CO"[2],"  (mmol per mol)",sep='')))


### Wet season threshold to get wet extremes =================================================================================

# find indices of all wet season observations 
all_wet_ind <- which(sav$season_overall == "Wet")

# make wet season data set 
wet_all <- sav[all_wet_ind,]
wet_NA_omit <- na.omit(wet_all$CO2)

# finding 90th percentile
wet_quant_90 <- quantile(wet_NA_omit, 0.9)

# Finding Threshold, mean excess, number of exceedances and confidence interval
mrl_wet <- evmix::mrlplot(wet_NA_omit)
num_wet <- mrl_wet$nu[c(1,18,32,47,60,78,100)]

# Plot threshold
plot.new()
grid(nx = NULL, ny = NULL,
     lty = 1, col = "lightgray", lwd = 0.6)
par(new = TRUE)
mrlplot(wet_all$CO2, nt = 200, col = c("seagreen", "black", "seagreen"), 
        lwd = c(1,2,1), xlim = c(420,545), ylim = c(0,40), 
        xlab = expression(paste("Threshold " , "(", italic(u) ,")",sep='')))
axis(3, at=mrl_wet$u[c(1,18,32,47,60,78,100)], lab=num_wet)
mtext("Number of Exceedances", side = 3, line = 2)
abline(v= c(wet_quant_90, 480), 
       col = c("darkorange", "yellowgreen"), 
       lwd = 1.6)
legend("topright", c("90th percentile", "Chosen threshold"), lty = 1, lwd = 2.4,
       col = c("darkorange", "yellowgreen"),
       bty = "n", cex = 0.7)


# RULE OF THUMB
wet_k1 <- round(sqrt(length(wet_NA_omit)))
wet_k2 <- round(((length(wet_NA_omit))^(2/3))/log(log(length(wet_NA_omit))))

wet_X1 <- length(wet_NA_omit) - wet_k1 
wet_X2 <- length(wet_NA_omit) - wet_k2 

wet_sort <- sort(wet_NA_omit)

# k values from rule (threshold values)
wet_rule1 <-  wet_sort[wet_X1] # says to use 488.4273
wet_rule2 <-  wet_sort[wet_X2] # says to use 472.2358

# number of event greater than the threshold rule
length(which(wet_NA_omit > wet_rule1))
length(which(wet_NA_omit > wet_rule2))

# setting wet threshold 
wet_thresh <- 480

### Wet season using above threshold  =================================================================================

# wet season 2 data set: wet_season2 

# number of events greater than the chosen threshold
length(which(wet_season2$CO2 > wet_thresh))

# indicator vector for wet whether CO2 is greater than the threshold
indicator_wet <- vector()

for (i in 1:length(wet_season2$CO2)){
  if (is.na(wet_season2$CO2[i])) {
    indicator_wet[i] = wet_season2$CO2[i]
  } else if (wet_season2$CO2[i] > wet_thresh) {
    indicator_wet[i] = "T"
  } else if (wet_season2$CO2[i] <= wet_thresh) {
    indicator_wet[i] = "F"
  }
}

# Add indicator vector to the wet season dataset 
wet_season2$indicator_wet <- indicator_wet

# Renaming to Indicator
names(wet_season2)[names(wet_season2) == "indicator_wet"] <- "Indicator"
wet_season2$indicator_wet <- indicator_wet


# Making magnitude above threshold vector (marks)
mag_wet <- vector()

for (i in 1:length(wet_season2$CO2)){
  if (is.na(wet_season2$CO2[i])) {
    mag_wet[i] = wet_season2$CO2[i]
  } else if (wet_season2$CO2[i] > wet_thresh) {
    mag_wet[i] = wet_season2$CO2[i] - wet_thresh
  } else if (wet_season2$CO2[i] <= wet_thresh) {
    mag_wet[i] = NA
  }
}

# Add magnitude vector to the wet season dataset 
wet_season2$mag_wet <- mag_wet

# making time index variable and adding to the wet season dataset 
wet_index <- seq(1:length(wet_index2))
wet_season2$wet_index <- wet_index

# get subset which only contains extremes 
wet_extremes <- wet_season2[which(wet_season2$indicator_wet == "T"),]


### Dry season threshold to get dry extremes =================================================================================

# find indices of all dry season observations 
all_dry_ind <- which(sav$season_overall == "Dry")

# make dry season data set 
dry_all <- sav[all_dry_ind,]
dry_NA_omit <- na.omit(dry_all$CO2)

# finding 90th percentile
dry_quant_90 <- quantile(dry_NA_omit, 0.9)

# Finding Threshold, mean excess, number of exceedances and confidence interval
mrl_dry <- evmix::mrlplot(na.omit(dry_all$CO2), legend.loc = NULL, bty = "n")
num_dry <- mrl_dry$nu[c(1,17,32,47,60,78,100)]

# Plot threshold
plot.new()
grid(nx = NULL, ny = NULL,
     lty = 1, col = "lightgray", lwd = 0.6)
par(new = TRUE)
mrlplot(dry_all$CO2, nt = 200, col = c("seagreen", "black", "seagreen"), 
        lwd = c(1,2,1), xlim = c(400,525), 
        xlab = expression(paste("Threshold " , "(", italic(u) ,")",sep='')))
axis(3, at=mrl_dry$u[c(1,17,32,47,60,78,100)], lab=num_dry)
mtext("Number of Exceedances", side = 3, line = 2)
abline(v= c(dry_quant_90, 430), 
       col = c("darkorange", "yellowgreen"), 
       lwd = 1.6)
legend("topright", c("90th percentile", "Chosen threshold"), lty = 1, lwd = 2.4,
       col = c("darkorange", "yellowgreen"),
       bty = "n", cex = 0.7)


# RULE OF THUMB
dry_k1 <- round(sqrt(length(dry_NA_omit)))
dry_k2 <- round(((length(dry_NA_omit))^(2/3))/log(log(length(dry_NA_omit))))

dry_X1 <- length(dry_NA_omit) - dry_k1 
dry_X2 <- length(dry_NA_omit) - dry_k2 

dry_sort <- sort(dry_NA_omit)

# k values from rule (threshold values)
dry_rule1 <-  dry_sort[dry_X1] # says to use  436.0575
dry_rule2 <-  dry_sort[dry_X2] # says to use 428.326

# number of event greater than the threshold rule
length(which(dry_NA_omit > dry_rule1))
length(which(dry_NA_omit > dry_rule2))

# setting dry threshold
dry_thresh <- 430

### Dry season 1 using above threshold  =================================================================================

# number of events greater than the chosen threshold
length(which(dry_season1$CO2 > dry_thresh))

# indicator vector for dry whether CO2 is greater than the threshold
indicator_dry1 <- vector()

for (i in 1:length(dry_season1$CO2)){
  if (is.na(dry_season1$CO2[i])) {
    indicator_dry1[i] = dry_season1$CO2[i]
  } else if (dry_season1$CO2[i] > dry_thresh) {
    indicator_dry1[i] = "T"
  } else if (dry_season1$CO2[i] <= dry_thresh) {
    indicator_dry1[i] = "F"
  }
}

# Add indicator vector to the dry season 1 dataset 
dry_season1$indicator_dry1 <- indicator_dry1

# Renaming to Indicator
names(dry_season1)[names(dry_season1) == "indicator_dry1"] <- "Indicator"
dry_season1$indicator_dry1 <- indicator_dry1


# Making magnitude above threshold vector (marks)
mag_dry1 <- vector()

for (i in 1:length(dry_season1$CO2)){
  if (is.na(dry_season1$CO2[i])) {
    mag_dry1[i] = dry_season1$CO2[i]
  } else if (dry_season1$CO2[i] > dry_thresh) {
    mag_dry1[i] = dry_season1$CO2[i] - dry_thresh
  } else if (dry_season1$CO2[i] <= dry_thresh) {
    mag_dry1[i] = NA
  }
}

# Add magnitude vector to the dry season 1 dataset 
dry_season1$mag_dry1 <- mag_dry1

# making time index variable and adding to the dry season 1 dataset 
dry_index1 <- seq(1:length(dry_index1))
dry_season1 <- cbind(dry_season1, dry_index1)

# get subset which only contains extremes 
dry_extremes1 <- dry_season1[which(dry_season1$indicator_dry1 == "T"),]


### Dry season 2 using above threshold  =================================================================================

# number of events greater than the chosen threshold
length(which(dry_season2$CO2 > dry_thresh))

# indicator vector for dry whether CO2 is greater than the threshold
indicator_dry2 <- vector()

for (i in 1:length(dry_season2$CO2)){
  if (is.na(dry_season2$CO2[i])) {
    indicator_dry2[i] = dry_season2$CO2[i]
  } else if (dry_season2$CO2[i] > dry_thresh) {
    indicator_dry2[i] = "T"
  } else if (dry_season2$CO2[i] <= dry_thresh) {
    indicator_dry2[i] = "F"
  }
}
# Add indicator vector to the dry season 2 dataset 
dry_season2$indicator_dry2 <- indicator_dry2

# Renaming to Indicator
names(dry_season2)[names(dry_season2) == "indicator_dry2"] <- "Indicator"
dry_season2$indicator_dry2 <- indicator_dry2


# Making magnitude above threshold vector (marks)
mag_dry2 <- vector()

for (i in 1:length(dry_season2$CO2)){
  if (is.na(dry_season2$CO2[i])) {
    mag_dry2[i] = dry_season2$CO2[i]
  } else if (dry_season2$CO2[i] > dry_thresh) {
    mag_dry2[i] = dry_season2$CO2[i] - dry_thresh
  } else if (dry_season2$CO2[i] <= dry_thresh) {
    mag_dry2[i] = NA
  }
}

# Add magnitude vector to the dry season 2 dataset 
dry_season2$mag_dry2 <- mag_dry2

# making time index variable and adding to the dry season 2 dataset 
dry_index2 <- seq(1:length(dry_index2))
dry_season2 <- cbind(dry_season2, dry_index2)

# get subset which only contains extremes 
dry_extremes2 <- dry_season2[which(dry_season2$indicator_dry2 == "T"),]


### Exploring the threshold =================================================================================

par(mfrow= c(1,2))
plot.new()
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "#f7f7f7") 
grid(nx = NULL, ny = NULL, lty = 1, col = "white", lwd = 2)
par(new = TRUE)
mrlplot(wet_all$CO2, nt = 200, col = c("seagreen", "black", "seagreen"), 
        lwd = c(1.2,2,1.2), xlim = c(420,545), ylim = c(0,40), lty = c(2,1,2),
        xlab = expression(paste("Threshold " , "(", italic(u) ,")",sep='')),
        main = "")
axis(3, at=mrl_wet$u[c(1,18,32,47,60,78,100)], lab=num_wet, cex.axis=0.9)
mtext("Number of Exceedances", side = 3, line = 2, cex = 0.9)
title("Wet season", line = 3)
abline(v= c(wet_quant_90, 480), 
       col = c("purple", "orange"), 
       lwd = 1.8)

plot.new()
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "#f7f7f7") 
grid(nx = NULL, ny = NULL, lty = 1, col = "white", lwd = 2)
par(new = TRUE)
mrlplot(dry_all$CO2, nt = 200, col = c("seagreen", "black", "seagreen"), 
        lwd = c(1.2,2,1.2), xlim = c(400,525), lty = c(2,1,2),
        xlab = expression(paste("Threshold " , "(", italic(u) ,")",sep='')),
        main = "", ylab = "")
axis(3, at=mrl_dry$u[c(1,17,32,47,60,78,100)], lab=num_dry, cex.axis=0.9)
mtext("Number of Exceedances", side = 3, line = 2, cex = 0.9)
abline(v= c(dry_quant_90, 430),  col = c("purple", "orange"),  lwd = 1.8)
title("Dry season", line = 3)
legend("topright", c("95% CI","90th percentile", "Chosen threshold"), lty = c(1,1,1), lwd = 2.4,
       col = c("seagreen","purple", "orange"),
       bty = "n", cex = 0.7)


## Plotting above threshold for dry season 1
# Create data frame 
dry1.df <- dry_season1
dry1.df$above_threshold <- dry1.df$CO2
dry1.df$above_threshold[which(dry_season1$CO2 < dry_thresh)] <- NA

# Plot 
dry1_above <- ggplot(data=dry1.df, aes(x= date_time, y = CO2)) +
  geom_line(aes(y =  CO2), na.rm = TRUE)+ ylim(350,550)+
  geom_line(aes(y = dry_thresh), colour = "firebrick3", lwd = 1.2)+
  geom_point(aes(x = date_time, y = above_threshold), na.rm = TRUE, colour = "firebrick3", cex = 1)+
  labs(title = "Dry Season 1", x ="Date", y = expression(paste("CO"[2],"  (mmol per mol)",sep=''))) 

## Plotting above threshold wet season
# Create data frame 
wet.df <- wet_season2
wet.df$above_threshold <- wet.df$CO2
wet.df$above_threshold[which(wet_season2$CO2 < wet_thresh)] <- NA

# Plot
wet_above <- ggplot(data=wet.df, aes(x= date_time, y = CO2)) +
  geom_line(aes(y =  CO2), na.rm = TRUE)+ ylim(350,550)+
  geom_line(aes(y = wet_thresh), colour = "lightseagreen", lwd = 1.2)+
  geom_point(aes(x = date_time, y = above_threshold), na.rm = TRUE, colour = "lightseagreen", cex = 1)+
  labs(title ="Wet Season", x ="Date", y = NULL) 

## Plotting above threshold dry season 2
# Create data frame 
dry2.df <- dry_season2
dry2.df$above_threshold <- dry2.df$CO2
dry2.df$above_threshold[which(dry_season2$CO2 < dry_thresh)] <- NA

# Plot
dry2_above <- ggplot(data=dry2.df, aes(x= date_time, y = CO2)) +
  geom_line(aes(y =  CO2), na.rm = TRUE)+ ylim(350,550) +
  geom_line(aes(y = dry_thresh), colour = "firebrick3", lwd = 1.2)+
  geom_point(aes(x = date_time, y = above_threshold), na.rm = TRUE, colour = "firebrick3", cex = 1)+
  labs(title = "Dry Season 2", x ="Date", y = NULL) 

# Above seasonal threshold
ggarrange(dry1_above, wet_above, dry2_above, ncol=3, common.legend = F)


# Comparing extremes 
par(mfrow = c(1,3))
plot(dry_extremes1$CO2 ~ dry_extremes1$date_time, type = "h", xlab = "Date", 
     ylab = expression(paste("CO"[2])), main = "Dry Season 1")
plot(wet_extremes$CO2 ~ wet_extremes$date_time, type = "h", xlab = "Date", 
     ylab = expression(paste("CO"[2])), main = "Wet Season")
plot(dry_extremes2$CO2 ~ dry_extremes2$date_time, type = "h", xlab = "Date", 
     ylab = expression(paste("CO")), main = "Dry Season 2")


### Making Variables required for hawkes & marked =================================================================================

# Time index: Entire time interval of that season 
# History: Index of all extreme events that are above the threshold
# Marks: Mark value for all the history events 
# IAT: Inter-arrival times between the histroy events 

# Making Variables for wet 
time_index_wet <- wet_season2$wet_index
history_wet <- which(wet_season2$indicator_wet == "T")
marks_wet <- wet_season2$mag_wet[history_wet]
IAT_wet <- history_wet[2:length(history_wet)] - history_wet[1:(length(history_wet)-1)]

# Making Variables for dry 1
time_index_dry1 <- dry_season1$dry_index1
history_dry1 <- which(dry_season1$indicator_dry == "T")
marks_dry1 <- dry_season1$mag_dry1[history_dry1]
IAT_dry1 <- history_dry1[2:length(history_dry1)] - history_dry1[1:(length(history_dry1)-1)]

# Making Variables for dry 2
time_index_dry2 <- dry_season2$dry_index2
history_dry2 <- which(dry_season2$indicator_dry == "T")
marks_dry2 <- dry_season2$mag_dry2[history_dry2]
IAT_dry2 <- history_dry2[2:length(history_dry2)] - history_dry2[1:(length(history_dry2)-1)]



### Exploring inter arrival times =================================================================================

## DRY 1 Inter-Arrival Times

# 5 number summary and other numerical summaries
round(fivenum(IAT_dry1),3) #  1    1    1    3 5342
mean(na.omit(IAT_dry1)) # 53.96319
sd(na.omit(IAT_dry1)) # 436.0286

# create dataframe
IAT_dry1_df <- as.data.frame(IAT_dry1)

# histogram
IAT_dry1_hist <- ggplot(IAT_dry1_df , aes(x = IAT_dry1)) + 
  geom_histogram(color = "black", fill = "firebrick3", na.rm= TRUE, bins = 60, alpha = 0.8, lwd = 0.4) +
  ggtitle("Dry extreme 1")+ xlab("Inter-arrival time (30-minutes)")+ylab("Count")

# boxplot
IAT_dry1_box <- ggplot(IAT_dry1_df, aes(y = IAT_dry1)) + 
  geom_boxplot(na.rm= TRUE, fill = "firebrick3")+ ylab("Inter-arrival time (30-minutes)")+  xlab("")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())


## WET Inter-Arrival Times

# 5 number summary and other numerical summaries
round(fivenum(IAT_wet),3) # 1    1    4   48    568
mean(na.omit(IAT_wet)) # 53.69784
sd(na.omit(IAT_wet)) # 107.2937

# create dataframe
IAT_wet_df <- as.data.frame(IAT_wet)

# histogram
IAT_wet_hist <- ggplot(IAT_wet_df , aes(x = IAT_wet)) + 
  geom_histogram(color = "black", fill = "lightseagreen", bins = 60, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + 
  ggtitle("Wet extreme")+ xlab("Inter-arrival time (30-minutes)")+ xlim(NA,600)+ylab("")

# boxplot
IAT_wet_box <- ggplot(IAT_wet_df, aes(y = IAT_wet)) + 
  geom_boxplot(na.rm= TRUE, fill = "lightseagreen", alpha = 0.8)+ ylab("")+  xlab("")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())



## DRY 2 Inter-Arrival Times

# 5 number summary and other numerical summaries
round(fivenum(IAT_dry2),3) #  1    1    1    4 3171
mean(na.omit(IAT_dry2)) # 95.96512
sd(na.omit(IAT_dry2)) # 440.5908

# create dataframe
IAT_dry2_df <- as.data.frame(IAT_dry2)

# histogram
IAT_dry2_hist <- ggplot(IAT_dry2_df , aes(x = IAT_dry2)) + 
  geom_histogram(color = "black", fill = "firebrick3", na.rm= TRUE, bins = 60, alpha = 0.8, lwd = 0.4) +
  ggtitle("Dry extreme 2")+ xlab("Inter-arrival time (30-minutes)")+ylab("")

# boxplot
IAT_dry2_box <- ggplot(IAT_dry2_df, aes(y = IAT_dry2)) + 
  geom_boxplot(na.rm= TRUE, fill = "firebrick3")+ ylab("")+  xlab("")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 


# all seasons Inter-Arrival times
ggarrange(IAT_dry1_hist, IAT_wet_hist, IAT_dry2_hist, IAT_dry1_box, IAT_wet_box, IAT_dry2_box, ncol=3, nrow = 2, common.legend = F)


### Exploring extreme marks =================================================================================

## DRY 1 marks

# 5 number summary and other numerical summaries
round(fivenum(marks_dry1),3) #  0.046   2.985   6.009  11.174 116.168
mean(na.omit(marks_dry1)) # 9.166
sd(na.omit(marks_dry1)) # 12.178

# create dataframe
marks_dry1_df <- data.frame(marks_dry1)

# histogram
marks_dry1_hist <- ggplot(marks_dry1_df , aes(x = marks_dry1)) + 
  geom_histogram(color = "black", fill = "firebrick3", bins = 60, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + 
  ggtitle("Dry extreme 1")+ xlab("Marks (mmol per mol)") +ylab("Count")

# boxplot
marks_dry1_box <- ggplot(marks_dry1_df, aes(y = marks_dry1)) + 
  geom_boxplot(na.rm= TRUE, fill = "firebrick3", alpha = 0.8)+ ylab("Marks (mmol per mol)")+  xlab("")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 



## WET marks

# 5 number summary and other numerical summaries
round(fivenum(marks_wet),3) #  0.097  8.906 24.713 40.657 69.300
mean(na.omit(marks_wet)) # 26.743
sd(na.omit(marks_wet)) # 19.664

# create dataframe
marks_wet_df <- data.frame(marks_wet)

# histogram
marks_wet_hist <- ggplot(marks_wet_df , aes(x = marks_wet)) + 
  geom_histogram(color = "black", fill = "lightseagreen", bins = 60, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + 
  ggtitle("Wet extreme")+ xlab("Marks (mmol per mol)") +ylab("")

# boxplot
marks_wet_box <- ggplot(marks_wet_df, aes(y = marks_wet)) + 
  geom_boxplot(na.rm= TRUE, fill = "lightseagreen", alpha = 0.8)+ ylab("")+  xlab("")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 


## DRY 2

# 5 number summary and other numerical summaries
round(fivenum(marks_dry2),3) #  0.148   3.339   9.496  37.097 115.562
mean(na.omit(marks_dry2)) # 23.324
sd(na.omit(marks_dry2)) # 28.374

# create dataframe
marks_dry2_df <- data.frame(marks_dry2)

# histogram
marks_dry2_hist <- ggplot(marks_dry2_df , aes(x = marks_dry2)) + 
  geom_histogram(color = "black", fill = "firebrick3", bins = 60, na.rm= TRUE, alpha = 0.8, lwd = 0.4) + 
  ggtitle("Dry extreme 2")+ xlab("Inter-arrival times")+ xlab("Marks (mmol per mol)") +ylab("")

# boxplot
marks_dry2_box <- ggplot(marks_dry2_df, aes(y = marks_dry2)) + 
  geom_boxplot(na.rm= TRUE, fill = "firebrick3", alpha = 0.8)+ ylab("")+  xlab("")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())


# both seasons marks
ggarrange(marks_dry1_hist, marks_wet_hist, marks_dry2_hist, marks_dry1_box, marks_wet_box, marks_dry2_box, ncol=3, nrow = 2, common.legend = F)


### QQ plots for marks 

# set seed for reproducible results
set.seed(10)

Mark_dist_QQ <- function(mark_dist){
  
  # Generating Probability values to find Quantiles
  p <- ppoints(100)
  
  # Observed Quantiles
  obs_mark <- quantile(mark_dist, p)
  
  # Theoretical mean
  mu <- 1/mean(mark_dist)
  
  
  # Finding QQ-line by simulation:
  
  # Making Matrix to store samples
  samp_matrix <- matrix(0, nrow = 1000, ncol = length(mark_dist))
  
  # Sample stored as row in matrix
  for (i in 1:1000){
    sample <- rexp(length(mark_dist), rate = mu)
    samp_matrix[i,] <- sample
  }
  
  # Making Quantile Matrix
  quant_mat <- matrix(0, nrow = 1000, ncol = 100)
  
  # Storing quantiles of samples in quantile matrix
  for (k in 1:1000){
    quant_mat[k,] <- quantile(samp_matrix[k,], p)
  }
  
  # Making aggregate quantile vectors
  mean <- rep(0, 100)
  upper <- rep(0, 100)
  lower <- rep(0, 100)
  
  # Aggregating Quantiles Across Samples
  for (s in 1:100){
    mean[s] <- mean(quant_mat[,s])
    upper[s] <- quantile(quant_mat[,s], 0.975)
    lower[s] <- quantile(quant_mat[,s], 0.025)
  }
  
  # Making Data Frame
  qline <- as.data.frame(cbind(obs_mark, mean))
  q_upper <- as.data.frame(cbind(mean, upper))
  q_lower <- as.data.frame(cbind(mean, lower))
  
  return(list(qline = qline, q_upper = q_upper, q_lower = q_lower))
}

## Dry 1 
QQ_mark_dry1 <- Mark_dist_QQ(marks_dry1)

dry1_markdist_qq <- ggplot(QQ_mark_dry1$qline) + geom_point(aes(x=mean, y = obs_mark)) +
  geom_line(data = QQ_mark_dry1$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = QQ_mark_dry1$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Dry extreme 1")


## Wet 
QQ_mark_wet <- Mark_dist_QQ(marks_wet)

wet_markdist_qq <- ggplot(QQ_mark_wet$qline) + geom_point(aes(x=mean, y = obs_mark)) +
  geom_line(data = QQ_mark_wet$q_lower, aes(x = mean, y = lower), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_line(data = QQ_mark_wet$q_upper, aes(x = mean, y = upper), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("") + 
  labs(title = "Wet extreme")


## Dry 2 
QQ_mark_dry2 <- Mark_dist_QQ(marks_dry2)

dry2_markdist_qq <- ggplot(QQ_mark_dry2$qline) + geom_point(aes(x=mean, y = obs_mark)) +
  geom_line(data = QQ_mark_dry2$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = QQ_mark_dry2$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("") + 
  labs(title = "Dry extreme 2")


ggarrange(dry1_markdist_qq, wet_markdist_qq, dry2_markdist_qq, ncol = 3)


### QQ Poisson plots =================================================================================

# Function For generations Poisson QQ plots:

Poisson_QQ <- function(obs_IAT){
  
  # Generating Probability values to find Quantiles 
  p <- ppoints(100)
  
  # Observed Quantiles
  obs_quant <- quantile(obs_IAT, p)
  
  # Theoretical Lambda
  lambda_test <- 1/mean(obs_IAT)
  
  # Finding QQ-line by simulation:
  
  # Making Matrix to store samples
  samp_matrix <- matrix(0, nrow = 1000, ncol = length(obs_IAT))
  
  # Sample stored as row in matrix
  for (i in 1:1000){
    sample <- rexp(length(obs_IAT), rate = lambda_test)
    samp_matrix[i,] <- sample
  }
  
  # Making Quantile Matrix
  quant_mat <- matrix(0, nrow = 1000, ncol = 100)
  
  # Storing quantiles of samples in quantile matrix
  for (k in 1:1000){
    quant_mat[k,] <- quantile(samp_matrix[k,], p)
  }
  
  # Making aggregate quantile vectors
  mean <- rep(0, 100)
  upper <- rep(0, 100)
  lower <- rep(0, 100)
  
  # Aggregating Quantiles Across Samples
  for (s in 1:100){
    mean[s] <- mean(quant_mat[,s])
    upper[s] <- quantile(quant_mat[,s], 0.975)
    lower[s] <- quantile(quant_mat[,s], 0.025)
  }
  
  # Making Data Frame
  qline <- as.data.frame(cbind(obs_quant, mean))
  q_upper <- as.data.frame(cbind(mean, upper))
  q_lower <- as.data.frame(cbind(mean, lower))
  
  return(list(qline = qline, q_upper = q_upper, q_lower = q_lower))
}

## DRY 1 POISSON QQ-PLOT
1/mean(IAT_dry1)

dry1_mat <- Poisson_QQ(IAT_dry1)

dry1_pois_qq <- ggplot(dry1_mat$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = dry1_mat$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = dry1_mat$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") 


# WET POISSON QQ-PLOT
1/mean(IAT_wet)

wet_mat <- Poisson_QQ(IAT_wet)

wet_pois_qq <- ggplot(wet_mat$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = wet_mat$q_lower, aes(x = mean, y = lower), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_line(data = wet_mat$q_upper, aes(x = mean, y = upper), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles")


# DRY 2 POISSON QQ-PLOT
1/mean(IAT_dry2)

dry2_mat <- Poisson_QQ(IAT_dry2)

dry2_pois_qq <- ggplot(dry2_mat$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = dry2_mat$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = dry2_mat$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles")



### Histograms for comparison =================================================================================

# DRY 1
samp_dry1_hist <- rexp(n=1000, rate = 1/(mean(IAT_dry1)))
samp_dry1_hist <- as.data.frame(samp_dry1_hist)

poisson_hist_dry1 <- ggplot() + geom_histogram(data = IAT_dry1_df, aes( x = IAT_dry1, y = ..density..),
                                               fill = "red", alpha = 1, col = "black", binwidth = 50, lwd = 0.3) +
  geom_histogram(data = samp_dry1_hist, aes(x=samp_dry1_hist,y = ..density..), 
                 fill = "lightblue", alpha = 0.7, col = "black", binwidth = 50, lwd = 0.3) +
  labs(x= "Inter-Arrival time (30-minutes)", y = "Density")


# WET 
samp_wet_hist <- rexp(n=1000, rate = 1/(mean(IAT_wet)))
samp_wet_hist <- as.data.frame(samp_wet_hist)

poisson_hist_wet <- ggplot() + geom_histogram(data = IAT_wet_df, aes( x = IAT_wet, y = ..density..),
                                              fill = "seagreen", alpha = 1, col = "black", binwidth = 20, lwd = 0.3) +
  geom_histogram(data = samp_wet_hist, aes(x=samp_wet_hist,y = ..density..), 
                 fill = "lightblue", alpha = 0.7, col = "black", binwidth = 20, lwd = 0.3) +
  labs(x= "Inter-Arrival time (30-minutes)", y = "Density")


# DRY 2
samp_dry2_hist <- rexp(n=1000, rate = 1/(mean(IAT_dry2)))
samp_dry2_hist <- as.data.frame(samp_dry2_hist)

poisson_hist_dry2 <- ggplot() + geom_histogram(data = IAT_dry2_df, aes( x = IAT_dry2, y = ..density..),
                                               fill = "red", alpha = 1, col = "black", binwidth = 35, lwd = 0.3) +
  geom_histogram(data = samp_dry2_hist, aes(x=samp_dry2_hist,y = ..density..), 
                 fill = "lightblue", alpha = 0.7, col = "black", binwidth = 35, lwd = 0.3) +
  labs(x= "Inter-Arrival time (30-minutes)", y = "Density")



ggarrange(dry1_pois_qq, poisson_hist_dry1, ncol=2, common.legend = F)
ggarrange(wet_pois_qq, poisson_hist_wet, ncol=2, common.legend = F)
ggarrange(dry2_pois_qq, poisson_hist_dry2, ncol=2, common.legend = F)


### Simulating Simple Point Process =================================================================================

# Function to simulate Point Plots:

pointplot <- function(history, IAT, t){
  
  actual <- history
  
  sample <- rexp(n=1, rate = 1/mean(IAT))
  
  i = 0
  while (i < length(t)){
    samp <- rexp(n=1, rate = 1/mean(IAT))
    i = i + samp
    sample <- append(sample, samp)
  }
  # generating sample history from sampled IAT
  
  sims <- sample[1]
  
  for (i in 2:length(sample)){
    val = sample[i]
    sims <- append(sims, sims[i-1]+val)
  }
  
  n <- length(sims)
  
  df_actual <- as.data.frame(cbind(actual, 0.1))
  colnames(df_actual) <- c("Observed", "Plot")
  df_sims <- as.data.frame(cbind(sims, 0))
  colnames(df_sims) <- c("Sims", "Plot")
  
  return(list(actual = df_actual, sims = df_sims, count = n))
  
}


## Dry 1

point_dry1 <- pointplot(history = history_dry1, IAT = IAT_dry1, t = time_index_dry1)

ggplot() + 
  geom_point(data = point_dry1$actual,aes (x=point_dry1$actual$Observed, y=point_dry1$actual$Plot,
                                           colour = "Observed"), alpha = 0.4, cex = 2) +
  geom_point(data = point_dry1$sims, aes(x=point_dry1$sims$Sims, y =point_dry1$sims$Plot, 
                                         colour = "Simulated"), alpha = 0.4, cex = 2) +
  labs(x = "Time Index", y = "") + 
  ylim(-0.05, 0.15) +  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  scale_colour_manual(name = 'Data Type', values=c(Observed = "firebrick3", 
                                                   Simulated = "seagreen")) 


## Wet 

point_wet <- pointplot(history = history_wet, IAT = IAT_wet, t = time_index_wet)

ggplot() + 
  geom_point(data = point_wet$actual, aes (x=point_wet$actual$Observed, y=point_wet$actual$Plot,
                                           colour = "Observed"), alpha = 0.4, cex = 2) +
  geom_point(data = point_wet$sims, aes(x=point_wet$sims$Sims, y =point_wet$sims$Plot, 
                                        colour = "Simulated"), alpha = 0.4, cex = 2) +
  labs(x = "Time Index", y = "") +
  ylim(-0.05, 0.15) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  scale_colour_manual(name = 'Data Type', values=c(Observed = "lightseagreen", 
                                                   Simulated = "seagreen"))

## Dry 2 

point_dry2 <- pointplot(history = history_dry2, IAT = IAT_dry2, t = time_index_dry2)

ggplot() + 
  geom_point(data = point_dry2$actual, aes (x=point_dry2$actual$Observed, y=point_dry2$actual$Plot,
                                            colour = "Observed"), alpha = 0.4, cex = 2) +
  geom_point(data = point_dry2$sims, aes(x=point_dry2$sims$Sims, y =point_dry2$sims$Plot, 
                                         colour = "Simulated"), alpha = 0.4, cex = 2) +
  labs(x = "Time Index", y = "") +
  ylim(-0.05, 0.15) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  scale_colour_manual(name = 'Data Type', values=c(Observed = "firebrick3", 
                                                   Simulated = "seagreen"))




### Finding hawkes optim values =================================================================================

# DRY SEASON 1

# Run optim grid search for dry season 1
hawkes_dry1_optim <- hawkes_optim(H_t = history_dry1, time_index = time_index_dry1)

# find minimum log-likelihood
hawkes_dry1_min <- which(hawkes_dry1_optim[,4] == min(hawkes_dry1_optim[,4]))

# store parameters and log-likelihood
pars_hawkes_dry1 <- hawkes_dry1_optim[hawkes_dry1_min,][1:3] # 0.003864217 0.323006370 0.409284526
ll_hawkes_dry1 <- hawkes_dry1_optim[hawkes_dry1_min,][4] # 486.6221

# Finding variance-covaraince matrix

hawkes_cov_dry1 <- cov(hawkes_dry1_optim[,1:3])

# Confidence interval 
hawkes_lam_dry1 <- c(pars_hawkes_dry1[1] - 1.96*sqrt(hawkes_cov_dry1[1,1]), pars_hawkes_dry1[1] + 1.96*sqrt(hawkes_cov_dry1[1,1]))
# 0.003630743 0.004097692

hawkes_alpha_dry1 <- c(pars_hawkes_dry1[2] - 1.96*sqrt(hawkes_cov_dry1[2,2]), pars_hawkes_dry1[2] + 1.96*sqrt(hawkes_cov_dry1[2,2]))
# 0.2665376 0.3794751

hawkes_beta_dry1 <- c(pars_hawkes_dry1[3] - 1.96*sqrt(hawkes_cov_dry1[3,3]), pars_hawkes_dry1[3] + 1.96*sqrt(hawkes_cov_dry1[3,3]))
# 0.3508647 0.4677043


# WET SEASON

# Run optim grid search for wet season
hawkes_wet_optim <- hawkes_optim(H_t = history_wet, time_index = time_index_wet)

# find minimum log-likelihood
hawkes_wet_min <- which(hawkes_wet_optim[,4] == min(hawkes_wet_optim[,4]))

# store parameters and log-likelihood
pars_hawkes_wet <- hawkes_wet_optim[hawkes_wet_min,][1:3] #  0.006515119 0.141358593 0.234531759
ll_hawkes_wet <- hawkes_wet_optim[hawkes_wet_min,][4] # 586.659

# Finding variance-covaraince matrix
hawkes_cov_wet <- cov(hawkes_wet_optim[,1:3])

# Confidence interval 
hawkes_lam_wet <- c(pars_hawkes_wet[1] - 1.96*sqrt(hawkes_cov_wet[1,1]), pars_hawkes_wet[1] + 1.96*sqrt(hawkes_cov_wet[1,1]))
# 0.006486809 0.006543429

hawkes_alpha_wet <- c(pars_hawkes_wet[2] - 1.96*sqrt(hawkes_cov_wet[2,2]), pars_hawkes_wet[2] + 1.96*sqrt(hawkes_cov_wet[2,2]))
# 0.1382549 0.1444623

hawkes_beta_wet <- c(pars_hawkes_wet[3] - 1.96*sqrt(hawkes_cov_wet[3,3]), pars_hawkes_wet[3] + 1.96*sqrt(hawkes_cov_wet[3,3]))
# 0.2302329 0.2388306



# DRY SEASON 2

# Run optim grid search for dry season 2
hawkes_dry2_optim <- hawkes_optim(H_t = history_dry2, time_index = time_index_dry2)

# find minimum log-likelihood
hawkes_dry2_min <- which(hawkes_dry2_optim[,4] == min(hawkes_dry2_optim[,4]))

# store parameters and log-likelihood
pars_hawkes_dry2 <- hawkes_dry2_optim[hawkes_dry2_min,][1:3] #  0.001512975 0.183655968 0.218108370
ll_hawkes_dry2 <- hawkes_dry2_optim[hawkes_dry2_min,][4] # 274.6863


# Finding variance-covaraince matrix

hawkes_cov_dry2 <- cov(hawkes_dry2_optim[,1:3])

# Confidence interval 
hawkes_lam_dry2 <- c(pars_hawkes_dry2[1] - 1.96*sqrt(hawkes_cov_dry2[1,1]), pars_hawkes_dry2[1] + 1.96*sqrt(hawkes_cov_dry2[1,1]))
# 0.001156483 0.001869466

hawkes_alpha_dry2 <- c(pars_hawkes_dry2[2] - 1.96*sqrt(hawkes_cov_dry2[2,2]), pars_hawkes_dry2[2] + 1.96*sqrt(hawkes_cov_dry2[2,2]))
#  -0.005531765  0.372843702

hawkes_beta_dry2 <- c(pars_hawkes_dry2[3] - 1.96*sqrt(hawkes_cov_dry2[3,3]), pars_hawkes_dry2[3] + 1.96*sqrt(hawkes_cov_dry2[3,3]))
# 0.0009164673 0.4353002736

### Finding marked hawkes optim values ===============================================================================

# DRY SEASON 1

# Run optim grid search for dry season 1
marked_dry1_optim <- marked_optim(H_t = history_dry1, time_index = time_index_dry1,   m = marks_dry1)

# Remove any rows with zeros
opt_dry1_zero <- which(marked_dry1_optim[,1] == 0)
marked_dry1_optim <- marked_dry1_optim[-(opt_dry1_zero),]

# find minimum log-likelihood
marked_dry1_min <- which(marked_dry1_optim[,5] == min(marked_dry1_optim[,5]))

# store parameters and log-likelihood
pars_marked_dry1 <- marked_dry1_optim[marked_dry1_min ,][1:4] # 0.003862897 0.322926403 0.409425733 0.000000100
ll_marked_dry1 <- marked_dry1_optim[marked_dry1_min ,][5] # 1013.961


# Finding variance-covaraince matrix


marked_cov_dry1 <- cov(marked_dry1_optim[,1:4])

# Confidence interval 
marked_lam_dry1 <- c(pars_marked_dry1[1] - 1.96*sqrt(marked_cov_dry1[1,1]), pars_marked_dry1[1] + 1.96*sqrt(marked_cov_dry1[1,1]))
# 0.003777911 0.003947884

marked_alpha_dry1 <- c(pars_marked_dry1[2] - 1.96*sqrt(marked_cov_dry1[2,2]), pars_marked_dry1[2] + 1.96*sqrt(marked_cov_dry1[2,2]))
# 0.2830382 0.3628146

marked_beta_dry1 <- c(pars_marked_dry1[3] - 1.96*sqrt(marked_cov_dry1[3,3]), pars_marked_dry1[3] + 1.96*sqrt(marked_cov_dry1[3,3]))
# 0.3677955 0.4510559

marked_gamma_dry1 <- c(pars_marked_dry1[4] - 1.96*sqrt(marked_cov_dry1[4,4]), pars_marked_dry1[4] + 1.96*sqrt(marked_cov_dry1[4,4]))
# (0 , 0) 



# WET SEASON 

# Run optim grid search for wet season
marked_wet_optim <- marked_optim(H_t = history_wet, time_index = time_index_wet,   m = marks_wet)

# Remove any rows with zeros
opt_wet_zero <- which(marked_wet_optim[,1] == 0)
marked_wet_optim <- marked_wet_optim[-(opt_wet_zero),]

# find minimum log-likelihood
marked_wet_min <- which(marked_wet_optim[,5] == min(marked_wet_optim[,5]))

# store parameters and log-likelihood
pars_marked_wet <- marked_wet_optim[marked_wet_min,][1:4] # 0.006499944 0.104583821 0.229674052 0.009820809
ll_marked_wet <- marked_wet_optim[marked_wet_min,][5] # 1185.693


# Finding variance-covaraince matrix

marked_cov_wet <- cov(marked_wet_optim[,1:4])

# Confidence interval 
marked_lam_wet <- c(pars_marked_wet[1] - 1.96*sqrt(marked_cov_wet[1,1]), pars_marked_wet[1] + 1.96*sqrt(marked_cov_wet[1,1]))
# 0.006467890 0.006531998

marked_alpha_wet <- c(pars_marked_wet[2] - 1.96*sqrt(marked_cov_wet[2,2]), pars_marked_wet[2] + 1.96*sqrt(marked_cov_wet[2,2]))
# 0.006467890 0.006531998

marked_beta_wet <- c(pars_marked_wet[3] - 1.96*sqrt(marked_cov_wet[3,3]), pars_marked_wet[3] + 1.96*sqrt(marked_cov_wet[3,3]))
# 0.2228150 0.2365331

marked_gamma_wet <- c(pars_marked_wet[4] - 1.96*sqrt(marked_cov_wet[4,4]), pars_marked_wet[4] + 1.96*sqrt(marked_cov_wet[4,4]))
# 0.009179949 0.010461669




# DRY SEASON 2 

# Run optim grid search for dry season 2
marked_dry2_optim <- marked_optim(H_t = history_dry2, time_index = time_index_dry2, m = marks_dry2)

# Remove any rows with zeros
opt_dry2_zero <- which(marked_dry2_optim[,1] == 0)
marked_dry2_optim <- marked_dry2_optim[-(opt_dry2_zero),]

# find minimum log-likelihood
marked_dry2_min <- which(marked_dry2_optim[,5] == min(marked_dry2_optim[,5]))

# store parameters and log-likelihood
pars_marked_dry2 <- marked_dry2_optim[marked_dry2_min ,][1:4] # 0.001566468 0.168955241 0.222774947 0.004594195
ll_marked_dry2 <- marked_dry2_optim[marked_dry2_min ,][5] #  635.1447


# Finding variance-covaraince matrix

marked_cov_dry2 <- cov(marked_dry2_optim[,1:4])

# Confidence interval 
marked_lam_dry2 <- c(pars_marked_dry2[1] - 1.96*sqrt(marked_cov_dry2[1,1]), pars_marked_dry2[1] + 1.96*sqrt(marked_cov_dry2[1,1]))
#  0.001327266 0.001805671

marked_alpha_dry2 <- c(pars_marked_dry2[2] - 1.96*sqrt(marked_cov_dry2[2,2]), pars_marked_dry2[2] + 1.96*sqrt(marked_cov_dry2[2,2]))
# 0.05681922 0.28109126

marked_beta_dry2 <- c(pars_marked_dry2[3] - 1.96*sqrt(marked_cov_dry2[3,3]), pars_marked_dry2[3] + 1.96*sqrt(marked_cov_dry2[3,3]))
# 0.1188935 0.3266564

marked_gamma_dry2 <- c(pars_marked_dry2[4] - 1.96*sqrt(marked_cov_dry2[4,4]), pars_marked_dry2[4] + 1.96*sqrt(marked_cov_dry2[4,4]))
# -0.004265844  0.013454234



###### Dry season 1 intensity estimation =================================================================================

# Hawkes
# Use the optimized parameters to get the conditional intensity over time 
hawkes_intensity_dry1  <- hawkes_intensity(pars = pars_hawkes_dry1, H_t = history_dry1, t = time_index_dry1)

# First moment 
hawkes_intens_exp_dry1 <- hawkes_intensity_expect(pars = pars_hawkes_dry1, t = time_index_dry1, intensity_vec = hawkes_intensity_dry1$hawkes_vec)

# Marked
# Use the optimized parameters to get the conditional intensity over time 
marked_intensity_dry1 <- marked_intensity(pars = pars_marked_dry1, H_t = history_dry1, t = time_index_dry1, m = marks_dry1)

# First moment 
marked_intens_exp_dry1 <- marked_intensity_expect(pars = pars_marked_dry1, t = time_index_dry1, intensity_vec = marked_intensity_dry1$hawkes_vec, m = marks_dry1)


# Create dots where events are for plotting
dot_dry1 <- vector()

for (i in 1:length(hawkes_intensity_dry1$plot_index)){
  if (any(hawkes_intensity_dry1$plot_index[i]== history_dry1)){
    dot_dry1 <- append(dot_dry1, -0.02)
  } 
  else {
    dot_dry1 <- append(dot_dry1, NA)
  }
}

# Create data frame of intensities for plotting
intensity_dry1.df <- as.data.frame(cbind(hawkes_intensity_dry1$hawkes_plot,  marked_intensity_dry1$hawkes_plot, hawkes_intensity_dry1$plot_index, dot_dry1))
colnames(intensity_dry1.df) <- c("Hawkes", "Marked", "Time Index", "Dot")

# Create data frame of expectation for plotting
expect_dry1.df <- data.frame(cbind(hawkes_intens_exp_dry1, marked_intens_exp_dry1, time_index_dry1))
colnames(expect_dry1.df) <- c("Hawkes_exp", "Marked_exp","Time Index")


# Plot of Hawkes 
ggplot(intensity_dry1.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line()+ geom_point(aes(y = Dot), col = "firebrick3", cex = 2, alpha = 0.4)+ ylim(NA, 1)+
  ylab("Conditional Intensity") +
  geom_line(data = expect_dry1.df, aes(y = `Hawkes_exp`), col = "orange", lty = 2, lwd = 0.8)

# Plot of Hawkes sub interval
ggplot(intensity_dry1.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line() + coord_cartesian(xlim = c(50, 160)) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4) + ylim(NA, 1)+
  ylab("Conditional Intensity") +
  geom_line(data = expect_dry1.df, aes(y = `Hawkes_exp`), col = "orange", lty = 2, lwd = 1)


# Plot of Marked
ggplot(intensity_dry1.df, aes(x= `Time Index`, y = `Marked`))+
  geom_line()+ geom_point(aes(y = Dot), col = "firebrick3", cex = 2, alpha = 0.4) + ylim(NA, 1)+
  ylab("Conditional Intensity") +
  geom_line(data = expect_dry1.df, aes(y = `Marked_exp`), col = "orange", lty = 2, lwd = 0.8)


# Plot of Marked sub interval
ggplot(intensity_dry1.df, aes(x= `Time Index`, y = `Marked`))+
  geom_line() + coord_cartesian(xlim = c(50, 160)) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4) + ylim(NA, 1)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_dry1.df, aes(y = `Marked_exp`), col = "orange", lty = 2, lwd = 1)


# Compare Hawkes vs Marked 
ggplot(intensity_dry1.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line() + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4)+ ylim(NA, 1)+
  geom_line(aes(y = Marked), lty =1, lwd = 0.6, col = "seagreen", alpha = 0.5)+
  ylab("Conditional Intensity")

# Compare Hawkes vs Marked sub interval
ggplot(intensity_dry1.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line(lwd = 0.6, col = "purple", alpha = 0.7) + coord_cartesian(xlim = c(50, 160)) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4)+ ylim(NA, 1)+
  geom_line(aes(y = Marked), lty =1, lwd = 0.6, col = "orange", alpha = 0.9)+
  ylab("Conditional Intensity")


### Model Comparing - Log Likelihood, AIC, BIC
# Hawkes
AIC.BIC(pars_hawkes_dry1, history_dry1, log_likelihood = ll_hawkes_dry1)

# Marked 
AIC.BIC(pars_marked_dry1, history_dry1, log_likelihood = ll_marked_dry1)



###### Wet season intensity estimation =================================================================================

# Hawkes 
# Use the optimized parameters to get the conditional intensity over time 
hawkes_intensity_wet <- hawkes_intensity(pars = pars_hawkes_wet, H_t = history_wet, t = time_index_wet)

# First moment 
hawkes_intens_exp_wet <- hawkes_intensity_expect(pars = pars_hawkes_wet, t = time_index_wet, intensity_vec = hawkes_intensity_wet$hawkes_vec)

# Marked
# Use the optimized parameters to get the conditional intensity over time 
marked_intensity_wet <- marked_intensity(pars = pars_marked_wet, H_t = history_wet, t = time_index_wet, m = marks_wet)

# First moment 
marked_intens_exp_wet <- marked_intensity_expect(pars = pars_marked_wet, t = time_index_wet, intensity_vec = marked_intensity_wet$hawkes_vec, m = marks_wet)


# Create dots where events are for plotting
dot_wet <- vector()

for (i in 1:length(hawkes_intensity_wet$plot_index)){
  if (any(hawkes_intensity_wet$plot_index[i]== history_wet)){
    dot_wet <- append(dot_wet, -0.01)
  } 
  else {
    dot_wet <- append(dot_wet, NA)
  }
}

# Create data frame of intensities for plotting
intensity_wet.df <- as.data.frame(cbind(hawkes_intensity_wet$hawkes_plot,  marked_intensity_wet$hawkes_plot, 
                                        hawkes_intensity_wet$plot_index, dot_wet))
colnames(intensity_wet.df) <- c("Hawkes", "Marked", "Time Index", "Dot")

# Create data frame of expectation for plotting
expect_wet.df <- data.frame(cbind(hawkes_intens_exp_wet, marked_intens_exp_wet, time_index_wet))
colnames(expect_wet.df) <- c("Hawkes_exp", "Marked_exp","Time Index")


# Plot of Hawkes 
ggplot(intensity_wet.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line() + 
  geom_point(aes(y = Dot), col = "lightseagreen", cex = 2, alpha = 0.4) + ylim(NA, 0.6)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_wet.df, aes(y = `Hawkes_exp`), col = "orange", lty = 2, lwd = 0.8)

# Plot of Hawkes sub interval
ggplot(intensity_wet.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line() + coord_cartesian(xlim = c(3590,3700)) + 
  geom_point(aes(y = Dot), col = "lightseagreen", cex = 2.4, alpha = 0.4) + ylim(NA, 0.6)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_wet.df, aes(y = `Hawkes_exp`), col = "orange", lty = 2, lwd = 1)


# Plot of Marked
ggplot(intensity_wet.df, aes(x= `Time Index`, y = `Marked`))+
  geom_line()+ geom_point(aes(y = Dot), col = "lightseagreen", cex = 2, alpha = 0.4) + ylim(NA, 0.6)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_wet.df, aes(y = `Marked_exp`), col = "orange", lty = 2, lwd = 0.8)

# Plot of Marked sub interval 
ggplot(intensity_wet.df, aes(x= `Time Index`, y = `Marked`))+
  geom_line() + coord_cartesian(xlim = c(3590,3700)) + 
  geom_point(aes(y = Dot), col = "lightseagreen", cex = 2.4, alpha = 0.4) + ylim(NA, 0.6)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_wet.df, aes(y = `Marked_exp`), col = "orange", lty = 2, lwd = 1)


# Compare Hawkes vs Marked 
ggplot(intensity_wet.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line(lwd = 0.6) + 
  geom_point(aes(y = Dot), col = "lightseagreen", cex = 2.4, alpha = 0.4)+ ylim(NA, 0.6)+
  geom_line(aes(y = Marked), lty =1, lwd = 0.6, col = "purple", alpha = 0.5)+
  ylab("Conditional Intensity")

# Compare Hawkes vs Marked sub interval
ggplot(intensity_wet.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line(lwd = 0.6, col = "purple", alpha = 0.7) + coord_cartesian(xlim = c(3590,3700)) + 
  geom_point(aes(y = Dot), col = "lightseagreen", cex = 2.4, alpha = 0.4)+ ylim(NA, 0.6)+
  geom_line(aes(y = Marked), lty =1, lwd = 0.6, col = "orange", alpha = 0.7)+
  ylab("Conditional Intensity")


### Model Comparing - Log Likelihood, AIC, BIC
# Hawkes
AIC.BIC(pars_hawkes_wet, history_wet, ll_hawkes_wet)

# Marked 
AIC.BIC(pars_marked_wet, history_wet, ll_marked_wet)



###### Dry season 2 intensity estimation =================================================================================

# Hawkes
# Use the optimized parameters to get the conditional intensity over time 
hawkes_intensity_dry2  <- hawkes_intensity(pars = pars_hawkes_dry2, H_t = history_dry2, t = time_index_dry2)

# First moment 
hawkes_intens_exp_dry2 <- hawkes_intensity_expect(pars = pars_hawkes_dry2, t = time_index_dry2, intensity_vec = hawkes_intensity_dry2$hawkes_vec)

# Marked
# Use the optimized parameters to get the conditional intensity over time  
marked_intensity_dry2 <- marked_intensity(pars = pars_marked_dry2, H_t = history_dry2, t = time_index_dry2, m = marks_dry2)

# First moment 
marked_intens_exp_dry2 <- marked_intensity_expect(pars = pars_marked_dry2, t = time_index_dry2, intensity_vec = marked_intensity_dry2$hawkes_vec, m = marks_dry2)


# Create dots where events are for plotting
dot_dry2 <- vector()

for (i in 1:length(hawkes_intensity_dry2$plot_index)){
  if (any(hawkes_intensity_dry2$plot_index[i]== history_dry2)){
    dot_dry2 <- append(dot_dry2, -0.02)
  } 
  else {
    dot_dry2 <- append(dot_dry2, NA)
  }
}

# Create data frame of intensities for plotting
intensity_dry2.df <- as.data.frame(cbind(hawkes_intensity_dry2$hawkes_plot,  marked_intensity_dry2$hawkes_plot, hawkes_intensity_dry2$plot_index, dot_dry2))
colnames(intensity_dry2.df) <- c("Hawkes", "Marked", "Time Index", "Dot")

# Create data frame of expectation for plotting
expect_dry2.df <- data.frame(cbind(hawkes_intens_exp_dry2, marked_intens_exp_dry2, time_index_dry2))
colnames(expect_dry2.df) <- c("Hawkes_exp", "Marked_exp","Time Index")


# Plot of Hawkes 
ggplot(intensity_dry2.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line()+ geom_point(aes(y = Dot), col = "firebrick3", cex = 2, alpha = 0.4)+ ylim(NA, 0.9)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_dry2.df, aes(y = `Hawkes_exp`), col = "orange", lty = 2, lwd = 0.8)


# Plot of Hawkes sub interval
ggplot(intensity_dry2.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line() + coord_cartesian(xlim = c(590,700)) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4) + ylim(NA, 0.9)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_dry2.df, aes(y = `Hawkes_exp`), col = "orange", lty = 2, lwd = 1)


# Plot of Marked
ggplot(intensity_dry2.df, aes(x= `Time Index`, y = `Marked`))+
  geom_line()+ geom_point(aes(y = Dot), col = "firebrick3", cex = 2, alpha = 0.4) + ylim(NA, 0.9)+
  ylab("Conditional Intensity") +
  geom_line(data = expect_dry2.df, aes(y = `Marked_exp`), col = "orange", lty = 2, lwd = 0.8)

# Plot of Marked sub interval
ggplot(intensity_dry2.df, aes(x= `Time Index`, y = `Marked`))+
  geom_line() + coord_cartesian(xlim = c(590,700)) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4) + ylim(NA, 0.9)+
  ylab("Conditional Intensity")+
  geom_line(data = expect_dry2.df, aes(y = `Marked_exp`), col = "orange", lty = 2, lwd = 1)


# Compare Hawkes vs Marked 
ggplot(intensity_dry2.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line(lwd = 0.6) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4)+ ylim(NA, 0.9)+
  geom_line(aes(y = Marked), lty =1, lwd = 0.6, col = "seagreen", alpha = 0.5)+
  ylab("Conditional Intensity")

# Compare Hawkes vs Marked sub interval
ggplot(intensity_dry2.df, aes(x= `Time Index`, y = `Hawkes`))+
  geom_line(lwd = 0.6, col = "purple", alpha = 0.7) + coord_cartesian(xlim = c(590,700)) + 
  geom_point(aes(y = Dot), col = "firebrick3", cex = 2.4, alpha = 0.4)+ ylim(NA, 0.9)+
  geom_line(aes(y = Marked), lty =1, lwd = 0.6, col = "orange", alpha = 0.7)+
  ylab("Conditional Intensity")


### Model Comparing - Log Likelihood, AIC, BIC
# Hawkes
AIC.BIC(pars_hawkes_dry2, history_dry2, log_likelihood = ll_hawkes_dry2)

# Marked 
AIC.BIC(pars_marked_dry2, history_dry2, log_likelihood = ll_marked_dry2)


###### Produce simulated data for QQ plots using optim values  =================================================================================

## Hawkes

# Wet 
set.seed(20)

hawkes_sim_wet_1000 <- Hawkes_QQ(N = 1000, IAT = IAT_wet, pars = pars_hawkes_wet, t = time_index_wet)
save(hawkes_sim_wet_1000, file = "wet_Hawkes_Sims_1000.Rdata")



# Dry 1
set.seed(30)

hawkes_sim_dry1_1000 <- Hawkes_QQ(N = 1000, IAT = IAT_dry1, pars = pars_hawkes_dry1, t = time_index_dry1)
save(hawkes_sim_dry1_1000, file = "dry1_Hawkes_Sims_1000.Rdata")


# Dry 2
set.seed(40)

hawkes_sim_dry2_1000 <- Hawkes_QQ(N = 1000, IAT = IAT_dry2, pars = pars_hawkes_dry2, t = time_index_dry2)
save(hawkes_sim_dry2_1000, file = "dry2_Hawkes_Sims_1000.Rdata")




## Marked Hawkes

# Wet 
set.seed(50)

marked_sim_wet_1000 <- Marked_QQ(N = 1000, IAT = IAT_wet, pars = pars_marked_wet, t = time_index_wet, m = marks_wet)
save(marked_sim_wet_1000, file = "wet_Marked_Sims_1000.Rdata")

# Dry 1
set.seed(60)

marked_sim_dry1_1000 <- Marked_QQ(N = 1000, IAT = IAT_dry1, pars = pars_marked_dry1, t = time_index_dry1, m = marks_dry1)
save(marked_sim_dry1_1000, file = "dry1_Marked_Sims_1000.Rdata")


# Dry 2
set.seed(70)

marked_sim_dry2_1000 <- Marked_QQ(N = 1000, IAT = IAT_dry2, pars = pars_marked_dry2, t = time_index_dry2, m = marks_dry2)
save(marked_sim_dry2_1000, file = "dry2_Marked_Sims_1000.Rdata")



###### Model checking - QQ plots  =================================================================================


# Loading Objects:

# Hakwes
load("wet_Hawkes_Sims_1000.Rdata")
load("dry1_Hawkes_Sims_1000.Rdata")
load("dry2_Hawkes_Sims_1000.Rdata")

# Marked
load("wet_Marked_Sims_1000.Rdata")
load("dry1_Marked_Sims_1000.Rdata")
load("dry2_Marked_Sims_1000.Rdata")


### Hakwes QQ Plot:

# Dry 1

qqhawkes_dry1 <- ggplot() + geom_point(data = hawkes_sim_dry1_1000$Sim_mean, aes(x = sim_mean, y = obs_quant)) +
  geom_line(data = hawkes_sim_dry1_1000$Sim_2.5, aes(x=sim_mean, y=sim_2.5), lty =2, col = "red", lwd = 0.6) +
  geom_line(data = hawkes_sim_dry1_1000$Sim_97.5, aes(x=sim_mean, y=sim_97.5), lty =2, col = "red", lwd = 0.6) +
  geom_abline(slope = 1, intercept = 0, lty =2, size = 0.6) +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") +
  labs(title = "Hawkes process")

# Wet

qqhawkes_wet <- ggplot() + geom_point(data = hawkes_sim_wet_1000$Sim_mean, aes(x = sim_mean, y = obs_quant)) +
  geom_line(data = hawkes_sim_wet_1000$Sim_2.5, aes(x=sim_mean, y=sim_2.5), lty =2, col = "seagreen", lwd = 0.6) +
  geom_line(data = hawkes_sim_wet_1000$Sim_97.5, aes(x=sim_mean, y=sim_97.5), lty =2, col = "seagreen", lwd = 0.6) +
  geom_abline(slope = 1, intercept = 0, lty =2, size = 0.6) +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") +
  labs(title = "Hawkes process")+ ylim(NA, 900)


# Dry 2

qqhawkes_dry2 <- ggplot() + geom_point(data = hawkes_sim_dry2_1000$Sim_mean, aes(x = sim_mean, y = obs_quant)) +
  geom_line(data = hawkes_sim_dry2_1000$Sim_2.5, aes(x=sim_mean, y=sim_2.5), lty =2, col = "red", lwd = 0.6) +
  geom_line(data = hawkes_sim_dry2_1000$Sim_97.5, aes(x=sim_mean, y=sim_97.5), lty =2, col = "red", lwd = 0.6) +
  geom_abline(slope = 1, intercept = 0, lty =2, size = 0.6) +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") +
  labs(title = "Hawkes process")+ ylim(NA, 3200)



### Marked QQ Plots:

# Dry 1:

qqmarked_dry1 <- ggplot() + geom_point(data = marked_sim_dry1_1000$Sim_mean, aes(x = sim_mean, y = obs_quant)) +
  geom_line(data = marked_sim_dry1_1000$Sim_2.5, aes(x=sim_mean, y=sim_2.5), lty =2, col = "red", lwd = 0.6) +
  geom_line(data = marked_sim_dry1_1000$Sim_97.5, aes(x=sim_mean, y=sim_97.5), lty =2, col = "red", lwd = 0.6) +
  geom_abline(slope = 1, intercept = 0, lty =2, size = 0.6) +
  xlab("Theoretical Quantiles") +
  ylab("") +
  labs(title = "Marked Hawkes process")

# Wet

qqmarked_wet <- ggplot() + geom_point(data = marked_sim_wet_1000$Sim_mean, aes(x = sim_mean, y = obs_quant)) +
  geom_line(data = marked_sim_wet_1000$Sim_2.5, aes(x=sim_mean, y=sim_2.5), lty =2, col = "seagreen", lwd = 0.6) +
  geom_line(data = marked_sim_wet_1000$Sim_97.5, aes(x=sim_mean, y=sim_97.5), lty =2, col = "seagreen", lwd = 0.6) +
  geom_abline(slope = 1, intercept = 0, lty =2, size = 0.6) +
  xlab("Theoretical Quantiles") +
  ylab("") +
  labs(title = "Marked Hawkes process")+ ylim(NA, 900)


# Dry 2

qqmarked_dry2 <- ggplot() + geom_point(data = marked_sim_dry2_1000$Sim_mean, aes(x = sim_mean, y = obs_quant)) +
  geom_line(data = marked_sim_dry2_1000$Sim_2.5, aes(x=sim_mean, y=sim_2.5), lty =2, col = "red", lwd = 0.6) +
  geom_line(data = marked_sim_dry2_1000$Sim_97.5, aes(x=sim_mean, y=sim_97.5), lty =2, col = "red", lwd = 0.6) +
  geom_abline(slope = 1, intercept = 0, lty =2, size = 0.6) +
  xlab("Theoretical Quantiles") +
  ylab("") +
  labs(title = "Marked Hawkes process")+ ylim(NA, 3200)


ggarrange(qqhawkes_dry1, qqmarked_dry1, ncol = 2)
ggarrange(qqhawkes_wet, qqmarked_wet, ncol = 2)
ggarrange(qqhawkes_dry2, qqmarked_dry2, ncol = 2)

### Model checking - counting process with moments =================================================================================

## Wet 
# Get counting process
counting_process_wet <- counting_process(H_t = history_wet, t = time_index_wet)

# Find hawkes and marked first moment 
hawkes_count_expect_wet <- hawkes_count_expect(pars_hawkes_wet, time_index_wet, hawkes_intensity_wet$hawkes_vec)
marked_count_expect_wet <- marked_count_expect(pars_marked_wet, time_index_wet, marked_intensity_wet$hawkes_vec, marks_wet)

# Add the first moments to dataframe
counting_process_wet$hawkes_exp <- hawkes_count_expect_wet
counting_process_wet$marked_exp <- marked_count_expect_wet

# plot
hawkes_count_wet.gg <- ggplot(counting_process_wet, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "(A) Counting process")+
  geom_line(aes(y = `hawkes_exp`), col = "lightseagreen", lwd = 0.8)

marked_count_wet.gg <- ggplot(counting_process_wet, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "(A) Counting process")+
  geom_line(aes(y = `marked_exp`), col = "lightseagreen", lwd = 0.8) 

both_count_wet.gg <- ggplot(counting_process_wet, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "Counting process")+
  geom_line(aes(y = `hawkes_exp`, col = "Hawkes"), lwd = 0.9)+
  geom_line(aes(y = `marked_exp`, col = "Marked"), lwd = 0.9)+
  scale_colour_manual(name = 'Expectation', values=c(Hawkes = "purple", Marked = "orange"))+
  theme(legend.position = "bottom")




## Dry 1
# Get counting process
counting_process_dry1 <- counting_process(H_t = history_dry1, t = time_index_dry1)

# Find hawkes and marked first moment 
hawkes_count_expect_dry1 <- hawkes_count_expect(pars_hawkes_dry1, time_index_dry1, hawkes_intensity_dry1$hawkes_vec)
marked_count_expect_dry1 <- marked_count_expect(pars_marked_dry1, time_index_dry1, marked_intensity_dry1$hawkes_vec, marks_dry1)

# Add the first moments to dataframe
counting_process_dry1$hawkes_exp <- hawkes_count_expect_dry1
counting_process_dry1$marked_exp <- marked_count_expect_dry1

# plot
hawkes_count_dry1.gg <- ggplot(counting_process_dry1, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "(A) Counting process")+
  geom_line(aes(y = `hawkes_exp`), col = "firebrick3", lwd = 0.8)

marked_count_dry1.gg <- ggplot(counting_process_dry1, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "(A) Counting process")+
  geom_line(aes(y = `marked_exp`), col = "firebrick3", lwd = 0.8) +
  geom_line(data = counting_process_dry1, aes(x = `t`,y = `hawkes_exp`), col = "red")

both_count_dry1.gg <- ggplot(counting_process_dry1, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "Counting process")+
  geom_line(aes(y = `hawkes_exp`, col = "Hawkes"), lwd = 0.9)+
  geom_line(aes(y = `marked_exp`, col = "Marked"), lwd = 0.9)+
  scale_colour_manual(name = 'Expectation', values=c(Hawkes = "purple", Marked = "orange"))+
  theme(legend.position = "bottom")



## Dry 2
# Get counting process
counting_process_dry2 <- counting_process(H_t = history_dry2, t = time_index_dry2)

# Find hawkes and marked first moment 
hawkes_count_expect_dry2 <- hawkes_count_expect(pars_hawkes_dry2, time_index_dry2, hawkes_intensity_dry2$hawkes_vec)
marked_count_expect_dry2 <- marked_count_expect(pars_marked_dry2, time_index_dry2, marked_intensity_dry2$hawkes_vec, marks_dry2)

# Add the first moments to dataframe
counting_process_dry2$hawkes_exp <- hawkes_count_expect_dry2
counting_process_dry2$marked_exp <- marked_count_expect_dry2

# plot
hawkes_count_dry2.gg <- ggplot(counting_process_dry2, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "(A) Counting process")+
  geom_line(aes(y = `hawkes_exp`), col = "firebrick3", lwd = 0.8)

marked_count_dry2.gg <- ggplot(counting_process_dry2, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "(A) Counting process")+
  geom_line(aes(y = `marked_exp`), col = "firebrick3", lwd = 0.8) 

both_count_dry2.gg <- ggplot(counting_process_dry2, aes(x = `t`, y = `count_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Count", title = "Counting process")+
  geom_line(aes(y = `hawkes_exp`, col = "Hawkes"), lwd = 0.9)+
  geom_line(aes(y = `marked_exp`, col = "Marked"), lwd = 0.9)+
  scale_colour_manual(name = 'Expectation', values=c(Hawkes = "purple", Marked = "orange"))+
  theme(legend.position = "bottom")


### Model checking - compensator =================================================================================

### Hawkes & Marked

## Wet 

# Obtain Hawkes and marked compensator vectors 
hawkes_compensator_wet <- hawkes_compensator(pars = pars_hawkes_wet, H_t = history_wet, t = time_index_wet)

marked_compensator_wet <- marked_compensator(pars = pars_marked_wet, H_t = history_wet, t = time_index_wet, m = marks_wet)

# plot
hawkes_comp_wet.gg <- ggplot(hawkes_compensator_wet, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")

marked_comp_wet.gg <- ggplot(marked_compensator_wet, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")

both_comp_wet.gg <- ggplot(hawkes_compensator_wet, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")+
  geom_line(data = marked_compensator_wet, aes(x = `t`, y = `comp_vec`), col = "orange", lwd = 0.7)



## Dry 1
# Obtain Hawkes and marked compensator vectors 
hawkes_compensator_dry1 <- hawkes_compensator(pars = pars_hawkes_dry1, H_t = history_dry1, t = time_index_dry1)

marked_compensator_dry1 <- marked_compensator(pars = pars_marked_dry1, H_t = history_dry1, t = time_index_dry1, m = marks_dry1)

# plot
hawkes_comp_dry1.gg <- ggplot(hawkes_compensator_dry1, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")

marked_comp_dry1.gg <- ggplot(marked_compensator_dry1, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")

both_comp_dry1.gg <- ggplot(hawkes_compensator_dry1, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")+
  geom_line(data = marked_compensator_dry1, aes(x = `t`, y = `comp_vec`), col = "orange", lwd = 0.7)




## Dry 2
hawkes_compensator_dry2 <- hawkes_compensator(pars = pars_hawkes_dry2, H_t = history_dry2, t = time_index_dry2)

marked_compensator_dry2 <- marked_compensator(pars = pars_marked_dry2, H_t = history_dry2, t = time_index_dry2, m = marks_dry2)

hawkes_comp_dry2.gg <- ggplot(hawkes_compensator_dry2, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")

marked_comp_dry2.gg <- ggplot(marked_compensator_dry2, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")

both_comp_dry2.gg <- ggplot(hawkes_compensator_dry2, aes(x = `t`, y = `comp_vec`))+
  geom_line(lwd = 0.7) + labs(x = "Time Index", y = "Expected count", title = "(B) Compensator")+
  geom_line(data = marked_compensator_dry2, aes(x = `t`, y = `comp_vec`), col = "orange", lwd = 0.7)



#### Random time change of event compensators ==================================================

## WET 

# Hawkes event compensator
hawkes_comp_time_wet <- data.frame(cbind(hawkes_compensator_wet$comp_vec[history_wet], seq(1,length(history_wet))))
colnames(hawkes_comp_time_wet) <- c("Expected count", "Event index")

# Hawkes event counting process 
hawkes_comp_counting_wet <- counting_process(H_t = hawkes_comp_time_wet$`Expected count`, t = seq(1:length(history_wet)))
colnames(hawkes_comp_counting_wet) <- c("Count", "Event index")

hawkes_count_time_wet.gg <- ggplot(hawkes_comp_counting_wet, aes(x = `Event index`, y = `Count`))+
  geom_line(lwd = 0.7)+geom_line(aes(x = `Event index`, y = `Event index`), col = "lightseagreen", lwd = 0.8)+ 
  labs(title = "(C) Transformed counting process")


# marked event compensator
marked_comp_time_wet <- data.frame(cbind(marked_compensator_wet$comp_vec[history_wet], seq(1,length(history_wet))))
colnames(marked_comp_time_wet) <- c("Expected count", "Event index")

# marked event counting process 
marked_comp_counting_wet <- counting_process(H_t = marked_comp_time_wet$`Expected count`, t = seq(1:length(history_wet)))
colnames(marked_comp_counting_wet) <- c("Count", "Event index")

marked_count_time_wet.gg <- ggplot(marked_comp_counting_wet, aes(x = `Event index`, y = `Count`))+
  geom_line(lwd = 0.7)+geom_line(aes(x = `Event index`, y = `Event index`), col = "lightseagreen", lwd = 0.8)+ 
  labs(title = "(C) Transformed counting process")

# Compare 
both_count_time_wet.gg <- ggplot(hawkes_comp_counting_wet, aes(x = `Event index`, y = `Count`))+
  geom_line(aes(col = "Hawkes"), lwd = 0.8)+ geom_line(aes(x = `Event index`, y = `Event index`), lwd = 0.6)+ 
  labs(title = "Transformed counting process")+
  geom_line(data = marked_comp_counting_wet, aes(x = `Event index`, y = `Count`, col = "Marked"), lwd = 0.8)+
  scale_colour_manual(name = 'Expectation', values=c(Hawkes = "purple", Marked = "orange"))+
  theme(legend.position = "bottom")




## DRY 1
# Hawkes event compensator
hawkes_comp_time_dry1 <- data.frame(cbind(hawkes_compensator_dry1$comp_vec[history_dry1], seq(1,length(history_dry1))))
colnames(hawkes_comp_time_dry1) <- c("Expected count", "Event index")

# Hawkes event counting process 
hawkes_comp_counting_dry1 <- counting_process(H_t = hawkes_comp_time_dry1$`Expected count`, t = seq(1:length(history_dry1)))
colnames(hawkes_comp_counting_dry1) <- c("Count", "Event index")

hawkes_count_time_dry1.gg <- ggplot(hawkes_comp_counting_dry1, aes(x = `Event index`, y = `Count`))+
  geom_line(lwd = 0.7)+geom_line(aes(x = `Event index`, y = `Event index`), col = "firebrick3", lwd = 0.8) + 
  labs(title = "(C) Transformed counting process")


# marked event compensator
marked_comp_time_dry1 <- data.frame(cbind(marked_compensator_dry1$comp_vec[history_dry1], seq(1,length(history_dry1))))
colnames(marked_comp_time_dry1) <- c("Expected count", "Event index")

# marked event counting process 
marked_comp_counting_dry1 <- counting_process(H_t = marked_comp_time_dry1$`Expected count`, t = seq(1:length(history_dry1)))
colnames(marked_comp_counting_dry1) <- c("Count", "Event index")

marked_count_time_dry1.gg <- ggplot(marked_comp_counting_dry1, aes(x = `Event index`, y = `Count`))+
  geom_line(lwd = 0.7)+geom_line(aes(x = `Event index`, y = `Event index`), col = "firebrick3", lwd = 0.8)+ 
  labs(title = "(C) Transformed counting process")

# Compare 
both_count_time_dry1.gg <- ggplot(hawkes_comp_counting_dry1, aes(x = `Event index`, y = `Count`))+
  geom_line(aes(col = "Hawkes"), lwd = 0.8)+ geom_line(aes(x = `Event index`, y = `Event index`), lwd = 0.6)+ 
  labs(title = "Transformed counting process")+
  geom_line(data = marked_comp_counting_dry1, aes(x = `Event index`, y = `Count`, col = "Marked"), lwd = 0.8)+
  scale_colour_manual(name = 'Expectation', values=c(Hawkes = "purple", Marked = "orange"))+
  theme(legend.position = "bottom")



## DRY 2
# Hawkes event compensator
hawkes_comp_time_dry2 <- data.frame(cbind(hawkes_compensator_dry2$comp_vec[history_dry2], seq(1,length(history_dry2))))
colnames(hawkes_comp_time_dry2) <- c("Expected count", "Event index")

# Hawkes event counting process 
hawkes_comp_counting_dry2 <- counting_process(H_t = hawkes_comp_time_dry2$`Expected count`, t = seq(1:length(history_dry2)))
colnames(hawkes_comp_counting_dry2) <- c("Count", "Event index")

hawkes_count_time_dry2.gg <- ggplot(hawkes_comp_counting_dry2, aes(x = `Event index`, y = `Count`))+
  geom_line(lwd = 0.7)+geom_line(aes(x = `Event index`, y = `Event index`), col = "firebrick3", lwd = 0.8)+ 
  labs(title = "(C) Transformed counting process")


# marked event compensator
marked_comp_time_dry2 <- data.frame(cbind(marked_compensator_dry2$comp_vec[history_dry2], seq(1,length(history_dry2))))
colnames(marked_comp_time_dry2) <- c("Expected count", "Event index")


# marked event counting process 
marked_comp_counting_dry2 <- counting_process(H_t = marked_comp_time_dry2$`Expected count`, t = seq(1:length(history_dry2)))
colnames(marked_comp_counting_dry2) <- c("Count", "Event index")

marked_count_time_dry2.gg <- ggplot(marked_comp_counting_dry2, aes(x = `Event index`, y = `Count`))+
  geom_line(lwd = 0.7)+geom_line(aes(x = `Event index`, y = `Event index`), col = "firebrick3", lwd = 0.8) + 
  labs(title = "(C) Transformed counting process")

# Compare 
both_count_time_dry2.gg <- ggplot(hawkes_comp_counting_dry2, aes(x = `Event index`, y = `Count`))+
  geom_line(aes(col = "Hawkes"), lwd = 0.8)+ geom_line(aes(x = `Event index`, y = `Event index`), lwd = 0.6)+ 
  labs(title = "Transformed counting process")+
  geom_line(data = marked_comp_counting_dry2, aes(x = `Event index`, y = `Count`, col = "Marked"), lwd = 0.8)+
  scale_colour_manual(name = 'Expectation', values=c(Hawkes = "purple", Marked = "orange"))+
  theme(legend.position = "bottom")



#### IAT of event compensators  ==================================================

# WET
hawkes_comp_time_wet_IAT <- hawkes_comp_time_wet$`Expected count`[2:length(hawkes_comp_time_wet$`Event index`)] - hawkes_comp_time_wet$`Expected count`[1:(length(hawkes_comp_time_wet$`Event index`)-1)]

marked_comp_time_wet_IAT <- marked_comp_time_wet$`Expected count`[2:length(marked_comp_time_wet$`Event index`)] - marked_comp_time_wet$`Expected count`[1:(length(marked_comp_time_wet$`Event index`)-1)]


# DRY 1
hawkes_comp_time_dry1_IAT <- hawkes_comp_time_dry1$`Expected count`[2:length(hawkes_comp_time_dry1$`Event index`)] - hawkes_comp_time_dry1$`Expected count`[1:(length(hawkes_comp_time_dry1$`Event index`)-1)]

marked_comp_time_dry1_IAT <- marked_comp_time_dry1$`Expected count`[2:length(marked_comp_time_dry1$`Event index`)] - marked_comp_time_dry1$`Expected count`[1:(length(marked_comp_time_dry1$`Event index`)-1)]


# DRY 2
hawkes_comp_time_dry2_IAT <- hawkes_comp_time_dry2$`Expected count`[2:length(hawkes_comp_time_dry2$`Event index`)] - hawkes_comp_time_dry2$`Expected count`[1:(length(hawkes_comp_time_dry2$`Event index`)-1)]

marked_comp_time_dry2_IAT <- marked_comp_time_dry2$`Expected count`[2:length(marked_comp_time_dry2$`Event index`)] - marked_comp_time_dry2$`Expected count`[1:(length(marked_comp_time_dry2$`Event index`)-1)]


#### Basic test - unit rate Poisson   ==================================================

set.seed(80)

# Function For generations unit rate Poisson QQ plots:
Poisson_unit_QQ <- function(obs_IAT){
  
  # Generating Probability values to find Quantiles at
  p <- ppoints(100)
  # Observed Quantiles
  obs_quant <- quantile(obs_IAT, p)
  
  # Finding QQ-line by simulation:
  
  # Making Matrix to store samples
  samp_matrix <- matrix(0, nrow = 1000, ncol = length(obs_IAT))
  
  # Sample stored as row in matrix
  for (i in 1:1000){
    sample <- rexp(length(obs_IAT), rate = 1)
    samp_matrix[i,] <- sample
  }
  
  # Making Quantile Matrix
  quant_mat <- matrix(0, nrow = 1000, ncol = 100)
  
  # Storing quantiles of samples in quantile matrix
  for (k in 1:1000){
    quant_mat[k,] <- quantile(samp_matrix[k,], p)
  }
  
  # Making aggregate quantile vectors
  mean <- rep(0, 100)
  upper <- rep(0, 100)
  lower <- rep(0, 100)
  
  # Aggregating Quantiles Across Samples
  for (s in 1:100){
    mean[s] <- mean(quant_mat[,s])
    upper[s] <- quantile(quant_mat[,s], 0.975)
    lower[s] <- quantile(quant_mat[,s], 0.025)
  }
  
  # Making Data Frame
  qline <- as.data.frame(cbind(obs_quant, mean))
  q_upper <- as.data.frame(cbind(mean, upper))
  q_lower <- as.data.frame(cbind(mean, lower))
  
  
  return(list(qline = qline, q_upper = q_upper, q_lower = q_lower))
  
}



# DRY 1 
hawkes_unit_QQ_dry1 <- Poisson_unit_QQ(obs_IAT = hawkes_comp_time_dry1_IAT)

marked_unit_QQ_dry1 <- Poisson_unit_QQ(obs_IAT = marked_comp_time_dry1_IAT)


qq_unit_hawkes_dry1 <- ggplot(hawkes_unit_QQ_dry1$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = hawkes_unit_QQ_dry1$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = hawkes_unit_QQ_dry1$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Hawkes process")

qq_unit_marked_dry1 <- ggplot(marked_unit_QQ_dry1$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = marked_unit_QQ_dry1$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = marked_unit_QQ_dry1$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Marked Hawkes process")

ggarrange(qq_unit_hawkes_dry1, qq_unit_marked_dry1, ncol = 2)


# WET 
hawkes_unit_QQ_wet <- Poisson_unit_QQ(obs_IAT = hawkes_comp_time_wet_IAT)

marked_unit_QQ_wet <- Poisson_unit_QQ(obs_IAT = marked_comp_time_wet_IAT)


qq_unit_hawkes_wet <- ggplot(hawkes_unit_QQ_wet$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = hawkes_unit_QQ_wet$q_lower, aes(x = mean, y = lower), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_line(data = hawkes_unit_QQ_wet$q_upper, aes(x = mean, y = upper), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Hawkes process")

qq_unit_marked_wet <- ggplot(marked_unit_QQ_wet$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = marked_unit_QQ_wet$q_lower, aes(x = mean, y = lower), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_line(data = marked_unit_QQ_wet$q_upper, aes(x = mean, y = upper), colour = "seagreen", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Marked Hawkes process")

ggarrange(qq_unit_hawkes_wet, qq_unit_marked_wet, ncol = 2)


# DRY 2
hawkes_unit_QQ_dry2 <- Poisson_unit_QQ(obs_IAT = hawkes_comp_time_dry2_IAT)

marked_unit_QQ_dry2 <- Poisson_unit_QQ(obs_IAT = marked_comp_time_dry2_IAT)


qq_unit_hawkes_dry2 <- ggplot(hawkes_unit_QQ_dry2$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = hawkes_unit_QQ_dry2$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = hawkes_unit_QQ_dry2$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Hawkes process")

qq_unit_marked_dry2 <- ggplot(marked_unit_QQ_dry2$qline) + geom_point(aes(x=mean, y = obs_quant)) +
  geom_line(data = marked_unit_QQ_dry2$q_lower, aes(x = mean, y = lower), colour = "red", lty = 2, lwd = 0.6) +
  geom_line(data = marked_unit_QQ_dry2$q_upper, aes(x = mean, y = upper), colour = "red", lty = 2, lwd = 0.6) +
  geom_abline(intercept = 0, slope = 1, size = 0.6, lty =2)  +
  xlab("Theoretical Quantiles") +
  ylab("Observed Quantiles") + 
  labs(title = "Marked Hawkes process")

ggarrange(qq_unit_hawkes_dry2, qq_unit_marked_dry2, ncol = 2)


### Basic test - Independence    ==================================================

# WET 
# Hawkes 
hawkes_Uk_wet <- 1 - exp(-hawkes_comp_time_wet_IAT)

hawkes_Uk_wet.df <- data.frame(cbind(hawkes_Uk_wet[1:(length(hawkes_Uk_wet)-1)], hawkes_Uk_wet[2:length(hawkes_Uk_wet)]))

hawkes_Uk_wet.gg <- ggplot(hawkes_Uk_wet.df, aes(x = X1, y = X2))+
  geom_point()+labs(x = expression(paste("U"[k])), y = expression(paste("U"[k+1])), title = "(D) Point process residuals")

# Marked
marked_Uk_wet <- 1 - exp(-marked_comp_time_wet_IAT)

marked_Uk_wet.df <- data.frame(cbind(marked_Uk_wet[1:(length(marked_Uk_wet)-1)], marked_Uk_wet[2:length(marked_Uk_wet)]))

marked_Uk_wet.gg <- ggplot(marked_Uk_wet.df, aes(x = X1, y = X2))+
  geom_point()+labs(x = expression(paste("U"[k])), y = expression(paste("U"[k+1])), title = "(D) Point process residuals")


# DRY 1
# Hawkes 
hawkes_Uk_dry1 <- 1 - exp(-hawkes_comp_time_dry1_IAT)

hawkes_Uk_dry1.df <- data.frame(cbind(hawkes_Uk_dry1[1:(length(hawkes_Uk_dry1)-1)], hawkes_Uk_dry1[2:length(hawkes_Uk_dry1)]))

hawkes_Uk_dry1.gg <- ggplot(hawkes_Uk_dry1.df, aes(x = X1, y = X2))+
  geom_point()+labs(x = expression(paste("U"[k])), y = expression(paste("U"[k+1])), title = "(D) Point process residuals")

# Marked
marked_Uk_dry1 <- 1 - exp(-marked_comp_time_dry1_IAT)

marked_Uk_dry1.df <- data.frame(cbind(marked_Uk_dry1[1:(length(marked_Uk_dry1)-1)], marked_Uk_dry1[2:length(marked_Uk_dry1)]))

marked_Uk_dry1.gg <- ggplot(marked_Uk_dry1.df, aes(x = X1, y = X2))+
  geom_point()+labs(x = expression(paste("U"[k])), y = expression(paste("U"[k+1])), title = "(D) Point process residuals")



# Hawkes 
hawkes_Uk_dry2 <- 1 - exp(-hawkes_comp_time_dry2_IAT)

hawkes_Uk_dry2.df <- data.frame(cbind(hawkes_Uk_dry2[1:(length(hawkes_Uk_dry2)-1)], hawkes_Uk_dry2[2:length(hawkes_Uk_dry2)]))

hawkes_Uk_dry2.gg <- ggplot(hawkes_Uk_dry2.df, aes(x = X1, y = X2))+
  geom_point()+labs(x = expression(paste("U"[k])), y = expression(paste("U"[k+1])), title = "(D) Point process residuals")

# Marked
marked_Uk_dry2 <- 1 - exp(-marked_comp_time_dry2_IAT)

marked_Uk_dry2.df <- data.frame(cbind(marked_Uk_dry2[1:(length(marked_Uk_dry2)-1)], marked_Uk_dry2[2:length(marked_Uk_dry2)]))

marked_Uk_dry2.gg <- ggplot(marked_Uk_dry2.df, aes(x = X1, y = X2))+
  geom_point()+labs(x = expression(paste("U"[k])), y = expression(paste("U"[k+1])), title = "(D) Point process residuals")


### ALL PLOTS

# DRY 1
ggarrange(hawkes_count_dry1.gg, hawkes_comp_dry1.gg, hawkes_count_time_dry1.gg, hawkes_Uk_dry1.gg, ncol = 2, nrow = 2)
ggarrange(marked_count_dry1.gg, marked_comp_dry1.gg, marked_count_time_dry1.gg, marked_Uk_dry1.gg, ncol = 2, nrow = 2)

ggarrange(both_count_dry1.gg, both_count_time_dry1.gg, ncol = 2, common.legend = T, legend = "bottom")

# WET
ggarrange(hawkes_count_wet.gg, hawkes_comp_wet.gg, hawkes_count_time_wet.gg, hawkes_Uk_wet.gg, ncol = 2, nrow = 2)
ggarrange(marked_count_wet.gg, marked_comp_wet.gg, marked_count_time_wet.gg, marked_Uk_wet.gg, ncol = 2, nrow = 2)

ggarrange(both_count_wet.gg, both_count_time_wet.gg, ncol = 2, common.legend = T, legend = "bottom")

# DRY 2
ggarrange(hawkes_count_dry2.gg, hawkes_comp_dry2.gg, hawkes_count_time_dry2.gg, hawkes_Uk_dry2.gg, ncol = 2, nrow = 2)
ggarrange(marked_count_dry2.gg, marked_comp_dry2.gg, marked_count_time_dry2.gg, marked_Uk_dry2.gg, ncol = 2, nrow = 2)

ggarrange(both_count_dry2.gg, both_count_time_dry2.gg, ncol = 2, common.legend = T, legend = "bottom")
