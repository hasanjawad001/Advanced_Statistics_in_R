################################################################################
# lab 06
################################################################################

library("gap")
rm(list=ls())

setwd("D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\lab_6")
myT <- read.table(
  "longitdunalRNASeqData\\nc101_scaff_dataCounts.txt",
  sep="\t",
  header=TRUE,
  row.names=1
)
print(dim(myT))

## remove rare genes
myT <- myT[apply(myT,1,median) > 5,]
print(dim(myT))

## normalization
myTNorm <- myT
for(i in 1:ncol(myT))
{
  colSum = sum(myT[,i])
  myTNorm[,i] =myTNorm[,i]/colSum
}
print(dim(myTNorm))

pValuesOneWayAnova <- vector()
pValuesRegression <- vector()
pValueModelDiff <- vector()

index <- vector()
cats <- factor( c( rep("day3",3),rep("week12",3),rep("week20",5)  ))
time <- c( rep(2,3),rep(86,3),rep(128,5)  )

## loop through the data
for(i in 1:nrow(myTNorm))
{
  index[i] <- i
  
  ## populate p-values
  ## build your linear models with myData as the y-variable 
  myData <- as.numeric( myTNorm[i,] ) 
  
  ## categories, 3 param model 
  ## (B0-background cat0-day2, B1-cat1-week12, B2-cat2-week20)
  myLm_cats <- lm( myData ~ cats, x=TRUE)
  myAnova_cats <- anova(myLm_cats)
  p_value_cats <- myAnova_cats$ "Pr(>F)"[1]
  pValuesOneWayAnova[i] <- p_value_cats
  
  ## regression, 2 param model (B0, B1)
  myLm_time <- lm( myData ~ time, x=TRUE)
  myAnova_time <- anova(myLm_time)
  p_value_time <- myAnova_time$ "Pr(>F)"[1]
  pValuesRegression[i] <- p_value_time  
  
  ## model diff
  fullresiduals <- sum(residuals(myLm_cats)^2)
  reducedresiduals <- sum(residuals(myLm_time)^2)
  F = ((reducedresiduals - fullresiduals)/(9-8))/(fullresiduals/8)
  p_value_diff <- pf(F, 1, 8, lower.tail = FALSE)
  pValueModelDiff[i] <- p_value_diff
}


################################ part (A)
print(length(pValuesOneWayAnova))
hist(pValuesOneWayAnova, plot = TRUE, breaks = 20)
qqunif(pValuesOneWayAnova, logscale=FALSE)

print(sum(pValuesOneWayAnova < 0.05)) ## 1442

pValuesOneWayAnova_bh <- p.adjust (pValuesOneWayAnova, method='BH')
print(sum(pValuesOneWayAnova_bh < 0.05)) ## 612 

pValuesOneWayAnova_bonf <- p.adjust (pValuesOneWayAnova, method='bonf')
print(sum(pValuesOneWayAnova_bonf < 0.05)) ## 7
################################ part (A) end


################################ part (B) start
print(length(pValuesRegression))
hist(pValuesRegression, plot = TRUE, breaks = 20)
qqunif(pValuesRegression, logscale=FALSE)

print(sum(pValuesRegression < 0.05)) ## 1297

pValuesRegression_bh <- p.adjust (pValuesRegression, method='BH')
print(sum(pValuesRegression_bh < 0.05)) ## 448

pValuesRegression_bonf <- p.adjust (pValuesRegression, method='bonf')
print(sum(pValuesRegression_bonf < 0.05)) ## 9
################################ part (B) end


################################ part (C) start
print(length(pValueModelDiff))
hist(pValueModelDiff, plot = TRUE, breaks = 20)
qqunif(pValueModelDiff, logscale=FALSE)

print(sum(pValueModelDiff < 0.05)) ## 732

pValueModelDiff_bh <- p.adjust (pValueModelDiff, method='BH')
print(sum(pValueModelDiff_bh < 0.05)) ## 51

pValueModelDiff_bonf <- p.adjust (pValueModelDiff, method='bonf')
print(sum(pValueModelDiff_bonf < 0.05)) ## 4
################################ part (C) end


################################ part (D) start

myFrame <- data.frame( 
  index, 
  pValuesOneWayAnova,
  pValuesRegression,
  pValueModelDiff
)

## OneWayAnova
myFrame_cats <- myFrame[ order(myFrame$pValuesOneWayAnova), ] 
print(myFrame_cats$index[1]) ## 2896
boxplot( 
  as.numeric( myTNorm[ myFrame_cats$index[1],]) ~ cats, 
  main="model (A)- most significant gene",
  xlab='categories',
  ylab='relative abundance'
)

## Regression
myFrame_time <- myFrame[ order(myFrame$pValuesRegression), ] 
print(myFrame_time$index[1]) ## 2915
boxplot(
  as.numeric( myTNorm[ myFrame_time$index[1],]) ~ time,
  main="model (B)- most significant gene",
  xlab='time (in days)',
  ylab='relative abundance'  
)
## For the graph of the top hit from (B), 
## include the regression line for the plot from (B).
relative_abundance <- as.numeric( myTNorm[ myFrame_time$index[1],])
myLm_time <- lm( relative_abundance ~ time)
plot(
  time, relative_abundance,
  main='model (B)- regression line')
abline(myLm_time)

## Model Diff
myFrame_diff <- myFrame[ order(myFrame$pValueModelDiff), ] 
print(myFrame_diff$index[1]) ## 2896
boxplot( 
  as.numeric( myTNorm[ myFrame_diff$index[1],]) ~ cats, 
  main="model (C)- most significant gene",
  xlab='categories',
  ylab='relative abundance'  
)

################################ part (D) end


################################ part (E) start
# I think the three param model (A) is more appropriate for these data, 
# cause the p values from ANOVA we see more number of significant 
# values/genes for model (A) than the model(B), which I assume shows 
# that with the three param model for most of the genes, we are 
# learning additional information.
################################ part (E) end
