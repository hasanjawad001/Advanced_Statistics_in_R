################################################################################
# Repeat the simulations for a patient without the disease. About how many times on avg.
# must the test be repeated to achieve the hospital's requirements?
################################################################################

## clean
rm(list=ls())

## likelihood
prob_like_sick <- c(0.91, 0.09)  ## (+ve, -ve) given sick patient
prob_like_healthy <- c(0.16, 0.84) ## (+ve, -ve) given healthy patient

limit_break <- 0.99999

## Avg. times must be repeated to be sure
num_trials <- 100000 ## considering 100k trials to take the avg
trial_no_vs_required_tests <- vector()

for (trial_no in 1:num_trials){
  
  ## prior and posterior
  prob_x_sick_healthy <- c(0.001, 0.999) ## sick, healthy
  
  i <- 0
  while (prob_x_sick_healthy[2] < limit_break){ ## condition to be sure healthy
    
    data <- vector()
    ## simulations considering patient without disease    
    if ( runif(1) <= prob_like_healthy[1]){
      data <- 1 ## test result +ve
    }else{
      data <- 2 ## test result -ve 
    }
    
    denom <- (prob_x_sick_healthy[1] * prob_like_sick[data]) + (prob_x_sick_healthy[2] * prob_like_healthy[data])
    prob_x_sick_healthy[1] = (prob_like_sick[data] * prob_x_sick_healthy[1])/denom ## sick
    prob_x_sick_healthy[2] = (prob_like_healthy[data] * prob_x_sick_healthy[2])/denom ## healthy
    
    i <- i + 1
  }
  trial_no_vs_required_tests[trial_no] <- i
}

paste('Avg. : ', mean(trial_no_vs_required_tests), ' for ', num_trials, ' trials!!!')