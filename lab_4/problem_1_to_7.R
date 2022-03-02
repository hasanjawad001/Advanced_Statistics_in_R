################################################################################
# lab 04
################################################################################

rm(list=ls())


## part 1: reading the dataset
setwd("D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\lab_4")
myT <- read.table(
  "longitdunalRNASeqData\\nc101_scaff_dataCounts.txt",
  header=TRUE,
  row.names=1
)
nrows = dim(myT)[1]
ncols = dim(myT)[2]
paste('nrows: ', nrows, ', ncols: ', ncols)


## part 2
col_D2_01 = myT[, c('D2_01')]
col_D2_02 = myT[, c('D2_02')]
plot(log10(col_D2_01), log10(col_D2_02), main= "D2_02 vs D2_01")


## part 3
plot(log10(apply(myT,1,mean)), log10(apply(myT,1,var)), main= "var vs mean") 
lines(log10(seq(0,nrows,0.1)), log10(seq(0,nrows,0.1)), col="red")


## part 4
m <- matrix(c(col_D2_01[1],col_D2_02[1],sum(col_D2_01)-col_D2_01[1],sum(col_D2_02)-col_D2_02[1]), nrow=2)
p_value = fisher.test(m, alternative = 'two.sided')$p.value
paste('p_value using first gene: ', p_value)


## part 5
# for these 2 samples (D2_01 and D2_02) run Fisher test considering every gene.
library(gap)
data_p_value <- vector(length = nrows)
sum_D2_01 = sum(col_D2_01)
sum_D2_02 = sum(col_D2_02)
for(i in 1: nrows){
  m <- matrix(c(col_D2_01[i], col_D2_02[i], sum_D2_01-col_D2_01[i], sum_D2_02-col_D2_02[i]), nrow=2)
  p_value = fisher.test(m, alternative = 'two.sided')$p.value
  data_p_value[i] = p_value
}
hist(data_p_value, plot = TRUE, breaks = 10)
qqunif(data_p_value, logscale=FALSE)
##################################################################
v2_myT <- myT[(myT$D2_01 + myT$D2_02 > 50),]
v2_nrows = dim(v2_myT)[1]
v2_ncols = dim(v2_myT)[2]
paste('v2_nrows: ', v2_nrows, ', v2_ncols: ', v2_ncols)

v2_col_D2_01 = v2_myT[, c('D2_01')]
v2_col_D2_02 = v2_myT[, c('D2_02')]
v2_data_p_value <- vector(length = v2_nrows)
v2_sum_D2_01 = sum(v2_col_D2_01)
v2_sum_D2_02 = sum(v2_col_D2_02)
for(i in 1: v2_nrows){
  m <- matrix(c(v2_col_D2_01[i], v2_col_D2_02[i], v2_sum_D2_01-v2_col_D2_01[i], v2_sum_D2_02-v2_col_D2_02[i]), nrow=2)
  p_value = fisher.test(m, alternative = 'two.sided')$p.value
  v2_data_p_value[i] = p_value
}
hist(v2_data_p_value, plot = TRUE, breaks = 10)
qqunif(v2_data_p_value, logscale=FALSE)


## part 6
v3_myT = myT + 1
v3_col_D2_01 = v3_myT[, c('D2_01')]
v3_col_D2_02 = v3_myT[, c('D2_02')]
exp_freq_p = v3_col_D2_01[1]/sum(v3_col_D2_01)
(exp_freq_p)
#
v3_pois_p_value = dpois(v3_col_D2_02[1], sum(v3_col_D2_02) * exp_freq_p)
paste('poisson p-value: ', v3_pois_p_value)


## part 7
v3_nrows = dim(v3_myT)[1]
v3_ncols = dim(v3_myT)[2]
v3_data_p_value <- vector(length = v3_nrows)
v3_sum_D2_01 = sum(v3_col_D2_01)
v3_sum_D2_02 = sum(v3_col_D2_02)
for(i in 1: v3_nrows){
  exp_freq_p = v3_col_D2_01[i]/v3_sum_D2_01
  v3_pois_p_value = dpois(v3_col_D2_02[i], v3_sum_D2_02 * exp_freq_p)
  v3_data_p_value[i] = v3_pois_p_value
}
hist(v3_data_p_value, plot = TRUE, breaks = 10)
qqunif(v3_data_p_value, logscale=FALSE)
#
plot(log10(data_p_value), log10(v3_data_p_value), main='(Y = myT+1) vs (X = myT)')

