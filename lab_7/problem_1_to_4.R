################################################################################
# lab 07
################################################################################
##
## 01
library("nlme")
rm(list=ls())
setwd(
  "D:\\2022_spring\\8310_bioinformatics_advanced_stats\\R_working_dir\\advanced_statistics\\lab_7"
)
inFileName <- paste("prePostPhylum.txt", sep ="")
myT <-read.table(inFileName,header=TRUE,sep="\t")
print(dim(myT)) ## 81 X 10
print(head(myT))
print(str(myT))
numCols <- ncol(myT)
myColClasses <- c(rep("character",4), rep("numeric", numCols-4))
print(myColClasses)
myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
myTData<-myT[,5:10]
myPCOA <- princomp(myTData)
print(summary(myPCOA))
##
## 02(A) graph by genotype
print(summary(myPCOA))
print(colnames(myT))
col_gt = myT['genotype'][,1]
unique_col_gt = unique(col_gt)
print(unique_col_gt)
col_gt_color = vector(length = length(col_gt))
for(i in 1:length(col_gt)){
  if ( col_gt[i]=='WT') {
    col_gt_color[i]='red'
  } else if (col_gt[i]=='10-/-') {
    col_gt_color[i]='black'
  } else{
    print('ERROR!!!')
  }
}
print(col_gt_color)
##
png(file="plot_2A_pca1_vs_pca2_genotype.png",width=600,height=300)
plot(myPCOA$scores[,1],myPCOA$scores[,2], col=col_gt_color,
     main = paste('pca 1 vs pca 2 by genotype')
)
legend(
  x='topright',
  legend=unique_col_gt,
  col=c('red', 'black'),
  lty=1:2, cex=0.8
)
dev.off()
##
## 02(B) graph by cage
print(summary(myPCOA))
print(colnames(myT))
col_cage = myT['cage'][,1]
unique_col_cage = unique(col_cage)
print(unique_col_cage)
col_cage_color = vector(length = length(col_cage))
for(i in 1:length(col_cage)){
  if ( col_cage[i]=='Cage3_WT') {
    col_cage_color[i]='red'
  } else if (col_cage[i]=='Cage2_WT') {
    col_cage_color[i]='blue'
  } else if (col_cage[i]=='Cage1_WT') {
    col_cage_color[i]='green'
  } else if (col_cage[i]=='Cage6_WT') {
    col_cage_color[i]='magenta'
  } else if (col_cage[i]=='Cage4_10-/-') {
    col_cage_color[i]='cyan'
  } else if (col_cage[i]=='Cage3_10-/-') {
    col_cage_color[i]='orange'
  } else if (col_cage[i]=='Cage5_10-/-') {
    col_cage_color[i]='black'
  } else if (col_cage[i]=='Cage5_WT') {
    col_cage_color[i]='white'
  } else if (col_cage[i]=='Cage4_WT') {
    col_cage_color[i]='burlywood'
  } else if (col_cage[i]=='Cage2_10-/-') {
    col_cage_color[i]='violet'
  } else if (col_cage[i]=='Cage1_10-/-') {
    col_cage_color[i]='yellow'
  } else{
    print('ERROR!!!')
  }
}
print(col_cage_color)
##
png(file="plot_2B_pca1_vs_pca2_cage.png",width=600,height=300)
plot(myPCOA$scores[,1],myPCOA$scores[,2], col=col_cage_color, 
     main = paste('pca 1 vs pca 2 by cage'))
legend(
  x='topright',
  legend=unique_col_cage,
  col=c('red', 'green', 'blue', 'magenta', 'cyan', 'orange', 'black',
        'white', 'burlywood', 'violet', 'yellow'),
  lty=1:2, cex=0.8
)
dev.off()
##
## 02(C) graph by time
print(summary(myPCOA))
print(colnames(myT))
col_time = myT['time'][,1]
unique_col_time = unique(col_time)
print(unique_col_time)
col_time_color = vector(length = length(col_time))
for(i in 1:length(col_time)){
  if ( col_time[i]=='PRE') {
    col_time_color[i]='red'
  } else if (col_time[i]=='POST') {
    col_time_color[i]='black'
  } else{
    print('ERROR!!!')
  }
}
print(col_time_color)
##
png(file="plot_2C_pca1_vs_pca2_time.png",width=600,height=300)
plot(myPCOA$scores[,1],myPCOA$scores[,2], col=col_time_color, 
     main = paste('pca 1 vs pca 2 by time'))
legend(
  x='topright',
  legend=unique_col_time,
  col=c('red', 'black'),
  lty=1:2, cex=0.8
)
dev.off()
##
##
## 03
y_pca1 <- myPCOA$scores[,1]
y_pca2 <- myPCOA$scores[,2]
## 03 (A) for cage, use a way one-ANOVA

x_cage <- factor(col_cage) 

myLm_cage_pca1 <- lm(y_pca1 ~ x_cage, x=TRUE)
myLm_cage_pca2 <- lm(y_pca2 ~ x_cage, x=TRUE)

myAnova_cage_pca1 <- anova(myLm_cage_pca1)
myAnova_cage_pca2 <- anova(myLm_cage_pca2)

p_value_cage_pca1 <- myAnova_cage_pca1$ "Pr(>F)"[1]
p_value_cage_pca2 <- myAnova_cage_pca2$ "Pr(>F)"[1]

print(p_value_cage_pca1) ## 0.992
print(p_value_cage_pca2) ## 1.63e-07

## 03 (B) for genotype, use a t-test
genotype <- myT$genotype
timepoint <- myT$time
myFrame2 <- data.frame(y_pca1, y_pca2, genotype, timepoint)

a_pca1=myFrame2[myFrame2$genotype=="WT",]$y_pca1
b_pca1=myFrame2[myFrame2$genotype=="10-/-",]$y_pca1
p_value_genotype_pca1 = t.test(a_pca1, b_pca1)$p.value
print(p_value_genotype_pca1) ## 0.929701

a_pca2=myFrame2[myFrame2$genotype=="WT",]$y_pca2
b_pca2=myFrame2[myFrame2$genotype=="10-/-",]$y_pca2
p_value_genotype_pca2 = t.test(a_pca2, b_pca2)$p.value
print(p_value_genotype_pca2) ## 1.274344e-10

## 03 (C) for timepoint, use a t-test
genotype <- myT$genotype
timepoint <- myT$time
myFrame2 <- data.frame(y_pca1, y_pca2, genotype, timepoint)

a_pca1=myFrame2[myFrame2$timepoint=="PRE",]$y_pca1
b_pca1=myFrame2[myFrame2$timepoint=="POST",]$y_pca1
p_value_timepoint_pca1 = t.test(a_pca1, b_pca1)$p.value
print(p_value_timepoint_pca1) ## 2.519974e-29

a_pca2=myFrame2[myFrame2$timepoint=="PRE",]$y_pca2
b_pca2=myFrame2[myFrame2$timepoint=="POST",]$y_pca2
p_value_timepoint_pca2 = t.test(a_pca2, b_pca2)$p.value
print(p_value_timepoint_pca2) ## 0.4268188

##
##
## Which variable seems to be most associated with the first PCA axis?
##    : Time because with cage and genotype the PCA 1 does not seem to vary 
##      or does not seem to have significant slope.
## Which variable is most associated with the second PCA axis?
##    : Genotype seems to be more associated to the second PCA axis 
##      than the other two.
## Does cage seem to be having an effect on these data?
##    : not on the first PCA axis but Yes, cage seems to be having 
##      an effect on the second PCA axis.
##
##
## 04
## 04 (A) for each phyla, graph relative abundance of that phyla vs cage
##
myT_post <- myT[myT$time=='POST',] 
dim(myT_post) ## 39 X 10
x=myT_post['cage'][,1]
data_col_y = c(
  'Tenericutes', 'Verrucomicrobia', 'Bacteroidetes',
  'Actinobacteria', 'Firmicutes', 'Proteobacteria'
)
png(file="plot_4A_relab_vs_cage_all_6.png",width=900,height=600, res = 90)
par( mfrow = c( 2, 3 ))
for (i in 1:length(data_col_y)){
  y = myT_post[data_col_y[i]][,1]
  boxplot( 
    as.numeric(y) ~ x, 
    main=paste('Phyla:', data_col_y[i]),
    xlab='',
    ylab='relative abundance',
    las=2
  )
}
dev.off()
##
##
## 04 (B) 
myT_post <- myT[myT$time=='POST',] 
dim(myT_post) ## 39 X 10
png(file="plot_4B_relab_vs_cage_all_6_p_rho.png",width=900,height=600, res = 90)
par( mfrow = c( 2, 3 ))
data_col_y = c(
  'Tenericutes', 'Verrucomicrobia', 'Bacteroidetes',
  'Actinobacteria', 'Firmicutes', 'Proteobacteria'
)
rho_gls = vector(mode='double', length=length(data_col_y))
p_value_lme = vector(mode='double', length=length(data_col_y))
for (i in 1:length(data_col_y)){
  
  bug <- myT_post[data_col_y[i]][,1]
  cage <- myT_post$cage
  genotype <- myT_post$genotype
  myFrame <- data.frame(bug, cage, genotype)
  
  M.gls <- gls( 
    bug ~ genotype , method = "REML",
    correlation = corCompSymm( form = ~ 1 | cage),data=myFrame
  )
  rho_gls[i] <- coef(M.gls$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]  
  
  M.mixed <- lme( 
    bug ~ genotype, method= "REML", 
    random = ~1 | cage, data = myFrame
  )
  p_value_lme[i] <- unclass(summary(M.mixed))$tTable[2,5] # H0:= no effect of cage 
  
  boxplot(
    myFrame$bug ~myFrame$cage,
    main=paste(
      'p_val, rho:', 
      format(round(p_value_lme[i], 2), nsmall = 2), 
      ', ', 
      format(round(rho_gls[i], 2), nsmall = 2)
    ),
    xlab='',
    ylab='',
    las=2
  )
  stripchart(bug ~ cage, data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
}
dev.off()
##
##
print(length(p_value_lme))
print(sum(p_value_lme < 0.10))
p_value_lme_adj <- p.adjust (p_value_lme, method='BH')
print(sum(p_value_lme_adj < 0.10))
##
##