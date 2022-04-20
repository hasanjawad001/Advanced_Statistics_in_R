################################################################################
# project final work
################################################################################

## check correlation matrix
library("gap")
setwd(
  "D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\project_final_work"
)
##
df_variable <- read.table(
  "variables.csv",
  header=TRUE,
  sep=','
  #row.names=1
)
print(dim(df_variable))
print(head(df_variable))
str(df_variable)
##
data_colname = colnames(df_variable)
print(data_colname)
##
data_colname_29 = data_colname[
  ! data_colname %in% c(data_colname[1], 'region_sim')
]
data_colname_16 = data_colname[
  data_colname %in% c(
    'Pop_bw20_34', 'logPop_tot',                              ## demographic
    'logMort_resp', 'logMort_infect',                          ## disease
    'GINI_index', 'Business',                                  ## economic
    'Temp', 'logPercip', 'logPollution',                       ## environment
    'Pop_Urban', 'Urban_pop',                                   ## habitat
    'GHS', 'logNurses',                                        ## health
    'Social_media', 'logInternet_filtering', 'log_Air_trans'   ## demographic    
  )
]
print(data_colname_29)
print(data_colname_16)
##
df_variable_29 = df_variable[, data_colname_29]
df_variable_16 = df_variable[, data_colname_16]
print(dim(df_variable))
print(dim(df_variable_29))
print(dim(df_variable_16))
##
## 1. plot correlation matrix
cr_29 = cor(df_variable_29)
cr_16 = cor(df_variable_16)
print(dim(cr_29))
print(dim(cr_16))
##
#install.packages("corrplot")
library(corrplot)
#corrplot(cr, method="circle")
#corrplot(cr, method="color")
#corrplot(cr, method="number")
##
pdf(file = "cr_29.pdf")
corrplot(
  cr_29, 
  method = "number", 
  #  type = "lower", 
  title = "correlation matrix over 29 covariates", 
  mar = c(0,0,1,0), 
  number.cex = 0.5, 
  number.digits = 2
)
dev.off()
##
#print(cr_29[,c('Log10GDP')])
##
pdf(file = "cr_16.pdf")
corrplot(
  cr_16, 
  method = "number", 
  #  type = "lower", 
  title = "correlation matrix over 16 covariates", 
  mar = c(0,0,1,0), 
  number.cex = 0.5, 
  number.digits = 2
)
dev.off()
## max correlation -0.65 between business and GHS
##
## 2. Analyse based on PCA, how many principle component for 95%? out of 16
print(dim(df_variable_29))
print(dim(df_variable_16))
pca_29 = prcomp(df_variable_29, center = TRUE,scale. = TRUE) ## princomp
pca_16 = prcomp(df_variable_16, center = TRUE,scale. = TRUE)
##
##
summary(pca_29)
summary(pca_16)
##
## 3. pair wise plotting within each category
#  'Pop_bw20_34', 'logPop_tot',                              ## demographic
#  'logMort_resp', 'logMort_infect',                          ## disease
#  'GINI_index', 'Business',                                  ## economic
#  'Temp', 'logPercip', 'logPollution',                       ## environment
#  'Pop_Urban', 'Urban_pop',                                   ## habitat
#  'GHS', 'logNurses',                                        ## health
#  'Social_media', 'logInternet_filtering', 'log_Air_trans'   ## social media    
##
print(data_colname_16)
##
val_col = c(5,12,13) ## social media (5,12,13) == (air_trans, social_media, internet)
val_title = '(Social Media)'
pairs(df_variable_16[val_col],
      col = 'blue', #modify color
      labels = data_colname_16[val_col], #modify labels
      main = paste('Pair-wise', val_title) ) #modify title
##
## group 16 into 7 categories
##
##
## 4. plot pca_1 vs pca_2 and see if it forms clusters by regions
##
print(dim(df_variable)) ## 58 X 31
print(dim(df_variable_29)) ## 58 X 29
print(dim(df_variable_16)) ## 58 X 16
##
pca_29 = prcomp(df_variable_29, center = TRUE,scale. = TRUE) ## princomp
pca_16 = prcomp(df_variable_16, center = TRUE,scale. = TRUE)
##
col_region = df_variable['region_sim'][,1]
col_color = vector(length = length(col_region))
unique_region = unique(col_region)
print(unique_region)
for(i in 1:length(col_region)){
  if ( col_region[i]=='AF') {
    col_color[i]='red'
  } else if (col_region[i]=='SA') {
    col_color[i]='blue'
  } else if (col_region[i]=='AS') {
    col_color[i]='green'
  } else if (col_region[i]=='EU') {
    col_color[i]='magenta'
  } else if (col_region[i]=='EUAS') {
    col_color[i]='cyan'
  } else if (col_region[i]=='NOR') {
    col_color[i]='orange'
  } else if (col_region[i]=='ME') {
    col_color[i]='black'
  } else{
    print('ERROR!!!')
  }
}
print(col_color)
##
png(file="pca1_vs_pca2_covariates_29.png",width=600,height=300)
plot(pca_29$x[,1],pca_29$x[,2], col=col_color)
legend(
  x='topright', 
  legend=c('AF', 'SA', 'AS', 'EU', 'EUAS', 'NOR', 'ME'), 
  col=c('red', 'green', 'blue', 'magenta', 'cyan', 'orange', 'black'), 
  lty=1:2, cex=0.8
)
dev.off()
##
png(file="pca1_vs_pca2_covariates_16.png",width=600,height=300)
plot(pca_16$x[,1],pca_16$x[,2], col=col_color)
legend(
  x='topright', 
  legend=c('AF', 'SA', 'AS', 'EU', 'EUAS', 'NOR', 'ME'), 
  col=c('red', 'green', 'blue', 'magenta', 'cyan', 'orange', 'black'), 
  lty=1:2, cex=0.8
)
dev.off()
##
## 5. 7C2=21 paired t-test to figure out how different each region from other 
list_r0 <- df_variable['R0'][,1]
col_r0 <- vector(length= length(list_r0))
for(i in 1:length(list_r0)){
  col_r0[i] = list_r0[i]
}
col_region <- factor(col_region)
#myLm <- lm(col_r0 ~ col_region, x=TRUE)
#anova(myLm)
#print(myLm$x)
pairwise.t.test(col_r0, col_region, p.adjust.method='BH')
##
## 6. replicating plot for baseline paper (nonlinear model) vs mine (linear model)
#install.packages("Hmisc")
library("Hmisc")
png(file="Ro_vs_covariates_16.png",width=1000,height=1000)
par(mfrow = c(4, 4))  
list_r0 <- df_variable['R0'][,1]
r0_max = max(list_r0)
r0_min = min(list_r0)
#print(r0_min)
#print(r0_max)
data_readable = c(
  'Youth', 'Total Pop', 'Mort Resp', 'Mort Infect',
  'GINI', 'Business', 'Temparature', 'Precipitation',
  'Pollution', 'City', 'Urbanization', 'GHS',
  'Nurses', 'Social Media', 'Internet Filtering', 'Air Transport'
)
print(data_colname_16)
for( i in 1: length(data_colname_16)){
  var <- data_colname_16[i]
  list_var <- df_variable[var][,1]
  myLm <- lm(list_r0 ~ list_var, x=TRUE)
  myAnova <- anova(myLm)
  p_value <- myAnova$"Pr(>F)"[1]
  p_value <- format(round(p_value, 2), nsmall = 2)
  plot(
    list_var, list_r0, 
    ylab = 'R0', xlab = data_readable[i],
    main = paste('(linear) p = ', p_value),
    ylim = c(r0_min-1, r0_max+1),
    cex.lab=1.5, cex.main=2
  )
  var_max = max(list_var)
  var_min = min(list_var)
  var_by = (var_max - var_min)/20
  xRange <- seq(from=var_min, to=var_max, by=var_by)
  linearMeans <- coef(myLm)[1] + coef(myLm)[2] * xRange
  meanR <- mean(residuals(myLm))
  standardError = sqrt(
    sum((residuals(myLm)-meanR)^2/(length(residuals(myLm))-2))
  )
  #abline(myLm)
  #lines(xRange, linearMeans)
  errbar(
    xRange, linearMeans, 
    linearMeans + (1*standardError), linearMeans - (1*standardError),
    add=TRUE
  )  
}
dev.off()
##
##