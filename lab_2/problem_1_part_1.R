################################################################################
# make a graph of the posterior probability of picking up a loaded die
################################################################################

##
rm(list=ls())

## prior and posterior
prob_x_fair_loaded <- c(0.99, 0.01) ## fair, loaded

## likelihood
prob_like_fair <- c(1/6,1/6,1/6,1/6,1/6,1/6) ## P(Y|X=fair)
prob_like_loaded <- c(1/10,1/10,1/10,1/10,1/10,1/2) ## P(Y|X=loaded)

## data
data <- c(2,3,2,6,3,5,6,2,6,6,2,6,6,2,3,6,6,6,5,6,6,5,6,6,6,6,6,4,6,3,3,3,6,6,5,6,6)

## graph of posterior P(loaded_die | data_number_of_times)
prob_x_fair = vector()
prob_x_loaded = vector()
title_str <- ''
for (i in 1:length(data)){
  
  prob_x_fair[i] <- prob_x_fair_loaded[1]  ## prior/posterior for fair
  prob_x_loaded[i] <- prob_x_fair_loaded[2]  ## prior/posterior for loaded 
  
  denom <- (prob_x_fair_loaded[1] * prob_like_fair[data[i]]) + (prob_x_fair_loaded[2] * prob_like_loaded[data[i]])
  prob_x_fair_loaded[1] = (prob_like_fair[data[i]] * prob_x_fair_loaded[1])/denom ## fair
  prob_x_fair_loaded[2] = (prob_like_loaded[data[i]] * prob_x_fair_loaded[2])/denom ## loaded
  
  print(prob_x_fair_loaded)
  title_str <- paste(title_str, data[i], sep='')
  plot(1:i,prob_x_loaded, main = title_str, ylim = c(0,1), xlim = c(1, length(data)+1) )
  Sys.sleep(1)
}

