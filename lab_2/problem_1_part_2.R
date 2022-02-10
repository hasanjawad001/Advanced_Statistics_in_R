################################################################################
# How many times on average would I need to roll a loaded die to be 99.999% sure 
# that it was loaded? 
################################################################################

## clean
rm(list=ls())

## likelihood
prob_like_fair <- c(1/6,1/6,1/6,1/6,1/6,1/6) ## P(Y|X=fair)
prob_like_loaded <- c(1/10,1/10,1/10,1/10,1/10,1/2) ## P(Y|X=loaded)


## Avg times needed to roll a loaded die to be 99.999% sure.
data_x <- c(1,2,3,4,5,6)
data_prob <- prob_like_loaded

limit_break <- 99.999/100 
num_trials <- 100000 ## considering 100k trials to take the avg
trial_no_vs_required_rolls <- vector()

for (trial_no in 1:num_trials){
  
  ## prior and posterior 
  prob_x_fair_loaded <- c(0.99, 0.01) ## fair, loaded
  
  i <- 0
  while (prob_x_fair_loaded[2] < limit_break){
    
    data <- sample(data_x, 1, replace=TRUE, prob=data_prob)
    
    denom <- (prob_x_fair_loaded[1] * prob_like_fair[data]) + (prob_x_fair_loaded[2] * prob_like_loaded[data])
    prob_x_fair_loaded[1] = (prob_like_fair[data] * prob_x_fair_loaded[1])/denom ## fair
    prob_x_fair_loaded[2] = (prob_like_loaded[data] * prob_x_fair_loaded[2])/denom ## loaded
    
    i <- i + 1
  }
  trial_no_vs_required_rolls[trial_no] <- i

}

paste('Avg. : ', mean(trial_no_vs_required_rolls), ' for ', num_trials, ' trials!!!')
