################################################################################
# How much should the hospital budget to run these tests?
################################################################################

## clean
rm(list=ls())

##
## of all the patients, we expect 0.1% with disease
prob_sample_sick_healthy <- c(0.001, 0.999) 

## likelihood
prob_like_sick <- c(0.91, 0.09)  ## (+ve, -ve) given sick
prob_like_healthy <- c(0.16, 0.84) ## (+ve, -ve) given healthy

limit_break <- 0.99999

## For a million patients, how many tests can the hospital anticipate running?

num_trials <- 10 ## considering 10 repeated trials with 1 million patients to take avg.
trial_no_vs_required_tests <- vector() ## this will keep track of total required test for 1 million at each trial

for (trial_no in 1:num_trials){
  
  num_patients <- 1000000
  patient_no_vs_required_tests <- vector() ## this will keep track of required test for each patient
  
  for (patient_no in 1:num_patients){
    
    data_patient <- vector()
    if ( runif(1) <= prob_sample_sick_healthy[1]){
      data_patient <- 1 ## current patient is sick
    }else{
      data_patient <- 2 ## current patient is healthy
    }    
    
    ## prior and posterior
    prob_x_sick_healthy <- c(0.001, 0.999) # initial belief about current patient  
    
    if (data_patient == 1){ ## patient is sick, data_patient==1
      i <- 0
      while (prob_x_sick_healthy[1] < limit_break){
        
        data_test <- vector()
        if ( runif(1) <= prob_like_sick[1]){
          data_test <- 1 ## test result +ve
        }else{
          data_test <- 2 ## test result -ve 
        }
        
        denom <- (prob_x_sick_healthy[1] * prob_like_sick[data_test]) + (prob_x_sick_healthy[2] * prob_like_healthy[data_test])
        prob_x_sick_healthy[1] = (prob_like_sick[data_test] * prob_x_sick_healthy[1])/denom ## sick
        prob_x_sick_healthy[2] = (prob_like_healthy[data_test] * prob_x_sick_healthy[2])/denom ## healthy
        
        i <- i + 1
      }
      patient_no_vs_required_tests[patient_no] <- i 
      
    }else{ ## patient is healthy, data_patient==2
      
      i <- 0
      while (prob_x_sick_healthy[2] < limit_break){
        
        data_test <- vector()
        if ( runif(1) <= prob_like_healthy[1]){
          data_test <- 1 ## test result +ve
        }else{
          data_test <- 2 ## test result -ve 
        }
        
        denom <- (prob_x_sick_healthy[1] * prob_like_sick[data_test]) + (prob_x_sick_healthy[2] * prob_like_healthy[data_test])
        prob_x_sick_healthy[1] = (prob_like_sick[data_test] * prob_x_sick_healthy[1])/denom ## sick
        prob_x_sick_healthy[2] = (prob_like_healthy[data_test] * prob_x_sick_healthy[2])/denom ## healthy
        
        i <- i + 1
      }
      patient_no_vs_required_tests[patient_no] <- i 
    } 
  }
  
  sum_test_for_million <- sum(patient_no_vs_required_tests)
  trial_no_vs_required_tests[trial_no] <- sum_test_for_million
  
}
print(trial_no_vs_required_tests)
paste('Avg. : ', mean(trial_no_vs_required_tests), ' for ', num_trials, ' trials!!!')
