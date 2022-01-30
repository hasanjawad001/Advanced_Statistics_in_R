

## data of x
data_x <- c(1,2,3,4,5,6)
print(data_x)


## probabilities of x
data_prob <- c(.1,.1,.1,.1,.1,.5)
print(data_prob)


## calculate mean
mean_x <- 0
for (i in 1: length(data_x)){
  mean_x <- mean_x + (data_prob[i] * data_x[i])
}
print(mean_x)


## calculate variance
var_x <- 0
for (i in 1: length(data_x)){
  diff <- (data_x[i]-mean_x)
  var_x <- var_x + (data_prob[i] * diff * diff)
}
print(var_x)


## function 'rollLoadedDie'
rollLoadedDie <- function(numRolls) {
  data_sample <- sample(data_x, numRolls, replace=TRUE, prob=data_prob)
  return (data_sample)
}


## make a histogram of some large number of rolls
data_rolls <- rollLoadedDie(1000000)
hist(data_rolls)
print(length(data_rolls))


## modify the code on Slide #58 of lecture #2
trial_sizes <- c(
  5,10,15,20,25,30,40,50,100,200,300,400,500,
  1000,2000,3000,4000,5000,10000,20000,30000,
  100000, 1000000, 10000000
)
means <- vector(mode='double', length=length(trial_sizes))
variances <- vector(mode='double', length=length(trial_sizes))

for (i in 1:length(trial_sizes)){
  rolls <- rollLoadedDie(trial_sizes[i])
  means[i] <- mean(rolls)
  variances[i] <- var(rolls)
}

plot(log10(trial_sizes), means)
lines(log10(trial_sizes), rep(mean_x, length(trial_sizes)))
windows()
plot(log10(trial_sizes), variances)
lines(log10(trial_sizes), rep(var_x, length(trial_sizes)))




