################################################################################
# project preliminary work
################################################################################

## check correlation matrix
library("gap")
setwd(
  "D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\project_preliminary_work"
)
df_variable <- read.table(
  "variables.csv",
  header=TRUE,
  sep=','
  #row.names=1
)
print(dim(df_variable))
print(head(df_variable))
str(df_variable)
#cov(df_variable)
#is.numeric(df_variable)
#is.logical(df_variable)
data_colname = colnames(df_variable)
#for (i in 1:length(data_colname)){
#  print(is.numeric(df_variable[, data_colname[i]]))
#}
print(data_colname)
##
data_colname_num = data_colname[
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
print(data_colname_num)
print(data_colname_16)


df_variable_num = df_variable[, data_colname_num]
df_variable_16 = df_variable[, data_colname_16]
print(dim(df_variable_num))
print(dim(df_variable_16))

## plot correlation matrix
#cv = cov(df_variable_num)
#print(dim(cv))
cr_num = cor(df_variable_num)
cr_16 = cor(df_variable_16)
print(dim(cr_num))
print(dim(cr_16))
install.packages("corrplot")
library(corrplot)
#corrplot(cr, method="circle")
#corrplot(cr, method="color")
#corrplot(cr, method="number")
##
pdf(file = "cr_num.pdf")
corrplot(
  cr_num, 
  method = "number", 
  #  type = "lower", 
  title = "correlation matrix over 29 covariates", 
  mar = c(0,0,1,0), 
  number.cex = 0.5, 
  number.digits = 2
)
dev.off()
##
#print(cr_num[,c('Log10GDP')])
##
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
##
## analyse based on PCA, how many principle component 
print(dim(df_variable_num))
print(dim(df_variable_16))
pca_29 = prcomp(df_variable_num, center = TRUE,scale. = TRUE)
pca_16 = prcomp(df_variable_16, center = TRUE,scale. = TRUE)
##
##
summary(pca_29)
summary(pca_16)
##
#  'Pop_bw20_34', 'logPop_tot',                              ## demographic
#  'logMort_resp', 'logMort_infect',                          ## disease
#  'GINI_index', 'Business',                                  ## economic
#  'Temp', 'logPercip', 'logPollution',                       ## environment
#  'Pop_Urban', 'Urban_pop',                                   ## habitat
#  'GHS', 'logNurses',                                        ## health
#  'Social_media', 'logInternet_filtering', 'log_Air_trans'   ## social media    
##
colnames(df_variable_16)
##
val_col = c(5,12,13)
val_title = '(Social Media)'
pairs(df_variable_16[val_col],
      col = 'blue', #modify color
      labels = data_colname_16[val_col], #modify labels
      main = paste('Pair-wise', val_title) ) #modify title
##