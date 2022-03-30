################################################################################
# lab 05
################################################################################

## problem (1A)
setwd("D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\lab_5")
df_cancer <- read.table(
  "cancerRisk.txt",header=TRUE,sep='\t'
  #row.names=1
)
#print(dim(df_cancer))
#print(head(df_cancer))
#str(df_cancer)
col_type = df_cancer[, c('Cancer_type')]
col_risk = df_cancer[, c('Lifetime_cancer_risk')]
col_div = df_cancer[, c('CumulativeCellDivisions')]
plot(log10(col_div), log10(col_risk), main= "risk vs cell div")
#plot(col_div, col_risk, main= "risk vs cell div")

## problem (1B)
model_lm <- lm(col_risk ~ col_div)
abline(model_lm, col = "red")

## problem (1C)
summary_lm = summary(model_lm)
#print(summary_lm)
print(summary_lm$coefficients[2,4])
print(summary_lm$r.squared)

## (1D)
resid_lm = resid(model_lm)
#print(resid_lm)
plot(model_lm)
#qqnorm(resid_lm)
## there is some deviation in the plot for residuals vs fitted values but for only a few samples

#----------------------------------------------------------------------------------------------#

## problem (2)
library("gap")

setwd("D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\lab_5")
df_case_control <- read.table(
  "caseControlData.txt",
  header=TRUE,
  sep='\t'
  #row.names=1
)
#print(dim(df_case_control))
#print(head(df_case_control))
#str(df_case_control)

setwd("D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\lab_5")
df_bmi <- read.table(
  "BMI_Data.txt",
  header=TRUE,
  sep='\t'
  #row.names=1
)
#print(dim(df_bmi))
#print(head(df_bmi))
#str(df_bmi)


## using similar column name for both dataframes
colnames(df_bmi) <- c('sample','bmi')
#str(df_bmi)

# cleaning un-necessary prefix-suffix from the case-control dataframe
convert_key <- function(key){
  # remove case and control
  key <- sub("case", "", key[1])
  key <- sub("control", "", key)
  
  # remove extraneous information from the suffix
  key <- strsplit( key, "_")[[1]][1]
  return (key)
}
newsample <- apply(df_case_control, 1, convert_key)
df_case_control$sample <- newsample

# verification of clean-up process
str(df_case_control)
str(df_bmi)

## way 1: by merging
#df_merge <- merge(df_case_control,df_bmi,by="sample")
#str(df_merge)

## way 2: by matching through 'which' condition
df_case_control$bmi <- df_bmi$bmi[which(df_bmi$sample %in% df_case_control$sample)]
str(df_case_control)

## skip column name 'sample', 'bmi' => remaining would be all 420 OTU columns
data_colname_otu = colnames(df_case_control)
data_colname_otu = data_colname_otu[! data_colname_otu %in% c('sample', 'bmi')]
length(data_colname_otu)
print(data_colname_otu)

## figure out all 420 p-values and check on that
data_p_value = vector(mode='double', length=length(data_colname_otu))
col_bmi = df_case_control[, c('bmi')]
for (col_no in 1: length(data_colname_otu)){
  colname_otu = data_colname_otu[col_no]
  #col_bmi = df_case_control[, c('bmi')] ## this would be same, no need to repeat
  col_otu = df_case_control[, colname_otu]
  
  model_lm <- lm(col_bmi ~ col_otu)
  #summary_lm = summary(model_lm)
  #p_value = summary_lm$coefficients[2,4] ## p-value from summary stat
  p_value = anova(model_lm)$"Pr(>F)"[1] ## p-value from anova stat
  data_p_value[col_no] = p_value
}
#hist(data_p_value)
hist(data_p_value, breaks=50)
qqunif(data_p_value,main='plot: qqunif')

## Are any of these associations significant at a 10% false discovery rate? : BH
print(paste(sum(data_p_value <= 0.10), 'out of',length(data_p_value))) ## 52 out of 420
data_p_value_BH = p.adjust(data_p_value, method = 'BH')
print(paste(sum(data_p_value_BH <= 0.10), 'out of',length(data_p_value))) ## 0 out of 420



