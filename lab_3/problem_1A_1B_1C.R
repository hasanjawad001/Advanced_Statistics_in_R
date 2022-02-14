################################################################################
# lab 03
################################################################################

rm(list=ls())

############################################################
## (1A) Plot Prior #########################################

numBreaks=1000;
xVals <- seq(0,1,1/numBreaks);
prob_x <- dexp(xVals, rate = 5)/0.9932621 
plot(xVals, prob_x, main = '(1A) exp prior') 

#############################################################
## (1B) Metropolis algorithm ################################

piOld <- 0.5
numIterations <- 500000 ## 500000
posteiorDist_metropolis <- vector(length=numIterations)
step_size <- 0.01
for( i in 1:numIterations )
{
  # our prior is exp distribution, our new data with 14 heads and 10 tails
  # pOld <- dbeta( piOld, 10,10 ) * dbinom( 14, 24, piOld ) ## update here
  pOld <- prob_x * dbinom( 14, 24, piOld )  
  piNew <- piOld + rnorm(1, 0, sd=step_size);
  if( piNew > 1) 
    piNew = 1;
  if( piNew < 0 ) 
    piNew =0;
  # pNew <- dbeta( piNew, 10,10 ) * dbinom( 14, 24, piNew ) ## update here
  pNew <- prob_x * dbinom( 14, 24, piNew )
  ratio <- pNew / pOld
  if( ratio > 1 || ratio >= runif(1) ) 
    piOld = piNew;
  posteiorDist_metropolis[i] = piOld;	
}
hist_metropolis <- hist(posteiorDist_metropolis,breaks=200,plot=FALSE)

#############################################################
## (1B) grid approximation ################################## 

posteriorDist_grid <- vector(length=length(hist_metropolis$mids))
i <- 1;
sum <- 0;
for( x in hist_metropolis$mids )
{
  # our prior is exp distribution, our new data with 14 heads and 10 tails
  # posteriorDist_grid[i] <- dbeta( x, 10,10 ) * dbinom( 14, 24, x) ## update here
  posteriorDist_grid[i] <- prob_x * dbinom( 14, 24, x)  
  sum = sum + posteriorDist_grid[i];
  i <- i + 1;	
}

ylim_max = max(hist_metropolis$counts/numIterations, posteriorDist_grid / sum)
ylim_max = max(ylim_max, dbeta(hist_metropolis$mids, 40+14, 40+10)/ sum(dbeta(hist_metropolis$mids, 40+14, 40+10)))

plot(
  hist_metropolis$mids, 
  hist_metropolis$counts/numIterations, 
  ylim=c(0,ylim_max), col='red',
  main = paste('(1B) numIterations: ', numIterations)
) 
lines(   
  hist_metropolis$mids, 
  posteriorDist_grid / sum, 
  ylim=c(0,ylim_max), col='blue' 
)
lines(
  hist_metropolis$mids, 
  dbeta(hist_metropolis$mids, 40+14, 40+10)/sum(dbeta(hist_metropolis$mids, 40+14, 40+10)), 
  ylim=c(0,ylim_max), col='green' 
)
legend(
  x='topright', 
  legend=c('metropolis', 'grid', 'beta prior'),
  col=c('red', 'blue', 'green'), lty=1, lwd=1, cex=0.8
)

#############################################################
## (1C) Metropolis algorithm ################################

piOld <- 0.5
numIterations <- 500000 ## 500000
posteiorDist_metropolis <- vector(length=numIterations)
step_size <- 0.01 ## 0.01
for( i in 1:numIterations )
{
  # our prior is exp distribution, our new data with 583 heads and 417 tails
  # pOld <- dbeta( piOld, 10,10 ) * dbinom( 14, 24, piOld ) ## update here
  pOld <- prob_x * dbinom( 583, 1000, piOld )  
  piNew <- piOld + rnorm(1, 0, sd=step_size);
  if( piNew > 1) 
    piNew = 1;
  if( piNew < 0 ) 
    piNew =0;
  # pNew <- dbeta( piNew, 10,10 ) * dbinom( 14, 24, piNew ) ## update here
  pNew <- prob_x * dbinom( 583, 1000, piNew )
  ratio <- pNew / pOld
  if( ratio > 1 || ratio >= runif(1) ) 
    piOld = piNew;
  posteiorDist_metropolis[i] = piOld;	
}
hist_metropolis <- hist(posteiorDist_metropolis,breaks=200,plot=FALSE)

#############################################################
## (1C) grid approximation ################################## 

posteriorDist_grid <- vector(length=length(hist_metropolis$mids))
i <- 1;
sum <- 0;
for( x in hist_metropolis$mids )
{
  # our prior is exp distribution, our new data with 583 heads and 417 tails
  # posteriorDist_grid[i] <- dbeta( x, 10,10 ) * dbinom( 14, 24, x) ## update here
  posteriorDist_grid[i] <- prob_x * dbinom( 583, 1000, x)  
  sum = sum + posteriorDist_grid[i];
  i <- i + 1;	
}

ylim_max = max(hist_metropolis$counts/numIterations, posteriorDist_grid / sum)
ylim_max = max(ylim_max, dbeta(hist_metropolis$mids, 40+583, 40+417)/ sum(dbeta(hist_metropolis$mids, 40+583, 40+417)))

plot(
  hist_metropolis$mids, 
  hist_metropolis$counts/numIterations, 
  ylim=c(0,ylim_max), col='red',
  main = paste('(1C) numIterations: ', numIterations)
) 
lines(   
  hist_metropolis$mids, 
  posteriorDist_grid / sum, 
  ylim=c(0,ylim_max), col='blue' 
)
lines(
  hist_metropolis$mids, 
  dbeta(hist_metropolis$mids, 40+583, 40+417)/sum(dbeta(hist_metropolis$mids, 40+583, 40+417)), 
  ylim=c(0,ylim_max), col='green' 
)
legend(
  x='topright', 
  legend=c('metropolis', 'grid', 'beta prior'),
  col=c('red', 'blue', 'green'), lty=1, lwd=1, cex=0.8
)

