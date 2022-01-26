################################################################################
#############_______________ CORRELATION ANALYSIS _______________###############
################################################################################

library(readxl)
library(corrplot)

#### IMPORT OF DATA ----

main <- read_excel("Datasets/Database_Humanitas_finale.xlsx")
pc <-read_excel("Datasets/Database_Humanitas_finale.xlsx", sheet = "PORTAL CORE")
pm <-read_excel("Datasets/Database_Humanitas_finale.xlsx", sheet = "PORTAL MARGIN")

## last row of portal margin is empty
pm <- pm[-261,]

## converting categorical variables in factors
main$PATTERN <- as.factor(main$PATTERN)
main$`T VIII ed` <- as.factor(main$`T VIII ed`)
main$N <- as.factor(main$N)
main$GRADING <- as.factor(main$GRADING)
main$SEX <- as.factor(main$SEX)
main$`NUMERO LINFONODI METASTATICI` <- as.factor(main$`NUMERO LINFONODI METASTATICI`)

## Removing some columns containing text that we will not consider
main$`Risposta radiologica`<-NULL
main$INTERVENTO <- NULL
main$Detail <- NULL
main$`MORTALITA' POSTOP` <- NULL
main$M <- NULL #unbalanced 

#### CORRELATION ANALYISIS OF PORTAL CORE AT 90% ----

# in order to reduce the dimensionality we remove the highly correlated columns
pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)]

R <- cor(pc)

# visualizing the correlation matrix
x11(width=10,height=10)
corrplot(R, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.01,main="Portal Core correlation")

## we build  matrix C which has a one in position (i,j) if columns i and j have 
## a correlation above 0.9
p = dim(R)[1]
C = matrix(0,p,p)
# I: vector first indeces
# J: vector second indeces
I <- NULL
J <- NULL

for (i in (2:p-1)) {
  for (j in (i+1):p) {
    if (abs(R[i,j]) >= 0.9) {
      I = c(I, i)
      J = c(J, j)
      C[i,j] = 1
    }
  }
}

corCorrelated_pc_0.9 = rbind(I,J)

# we visualize the names of the highly correlated columns: for each group, we will
# keep only one of these columns in our dataset
colnames(pc[,corCorrelated_pc_0.9[1,]])
colnames(pc[,corCorrelated_pc_0.9[2,]])


#### CORRELATION ANALYISIS OF PORTAL MARGIN AT 90% ----

# in order to reduce the dimensionality we remove the highly correlated columns
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)]

R2 <- cor(pm)

# visualizing the correlation matrix
x11(width=10,height=10)
corrplot(R2, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.01,main="Portal Margin correlation")


## we build  matrix C which has a one in position (i,j) if columns i and j have 
## a correlation above 0.9
p2 = dim(R2)[1]
C2 = matrix(0,p2,p2)
# I: vector first indeces
# J: vector second indeces
I <- NULL
J <- NULL

for (i in (2:p2-1)) {
  for (j in (i+1):p2) {
    if (abs(R2[i,j]) >= 0.9) {
      I = c(I, i)
      J = c(J, j)
      C2[i,j] = 1
    }
  }
}

Correlated_pm_0.9 = rbind(I,J)

# we visualize the names of the highly correlated columns: for each group, we will
# keep only one of these columns in our dataset
colnames(pm[,Correlated_pm_0.9[1,]])
colnames(pm[,Correlated_pm_0.9[2,]])

#Plot of pairs to understand dependencies between groups in PORTAL CORE
x11()
pairs(pc[,c(2,5,6,7,29)])
x11()
pairs(pc[,c(4,10)])
x11()
pairs(pc[,c(11,12,13,20,23,24)])
x11()
pairs(pc[,c(14,15,17,34,35)])
x11()
pairs(pc[,c(28,43)])
x11()
pairs(pc[,c(41,46,47)])
x11()
pairs(pc[,c(42,44)])
x11()
pairs(pc[,c(19,25,26,27,32,33,36)])


#Plot of pairs to understand dependencies between groups in PORTAL MARGIN
x11()
pairs(pm[,c(2,5,6,7,29)])
x11()
pairs(pm[,c(4,10)])
x11()
pairs(pm[,c(11,12,13,20,23,24)])
x11()
pairs(pm[,c(14,15,17,34,35)])
x11()
pairs(pm[,c(28,43)])
x11()
pairs(pm[,c(41,46,47)])
x11()
pairs(pm[,c(42,44)])


#### CREATION OF DATASETS ----
write.table(main, file ="main.txt")

# saving new datasets removing the highly correlated columns
pc <- pc[,-c(2,5,7,29,10,12,13,20,23,24,15,35,19,25,26,33,36,41,47,44)]
pc.def = c(pc.label, pc)
write.table(pc.def, file ="pc_post_correlation.txt")

pm <- pm[,-c(2,5,7,29,10,12,13,20,23,24,15,35,19,25,26,33,36,41,47,44)]

for (i in 1:length(colnames(pm))){
  colnames(pm)[i] <- paste(colnames(pm)[i],"margin",sep="_")
}

pm.def = c(pm.label, pm)
write.table(pm.def, file ="pm_post_correlation.txt")


## new correlation plots, in the same order for both the data-sheets
corrplot(R, type = "upper", 
         tl.col = "black", tl.cex = 0.01,main="Portal Core correlation")
corrplot(R2, type = "upper", 
         tl.col = "black", tl.cex = 0.01,main="Portal Margin correlation")


################################################################################
########_________ CORRELATION ANALYSIS after Outliers Removal _________#########
################################################################################

# removing the outliers individuated in Applied Stat and confirmed by 
# "3_Depth_outliers": 122 204 235 241
out <- c(122,204,235,241)

pc <-read_excel("Datasets/Database_Humanitas_finale.xlsx", sheet = "PORTAL CORE")
pm <-read_excel("Datasets/Database_Humanitas_finale.xlsx", sheet = "PORTAL MARGIN")
pm <- pm[-261,]
pc <- pc[,-c(1,2)]
pm <- pm[,-c(1,2)]

pc_good <- pc[-out,]
pm_good <- pm[-out,]

R <- cor(pc)
R2 <- cor(pm)

corrplot(R, type = "upper", 
         tl.col = "black", tl.cex = 0.01,title="Portal Core correlation")
corrplot(R2, type = "upper", 
         tl.col = "black", tl.cex = 0.01,title="Portal Margin correlation")

### Same as before!!

################################################################################
#######________ CORRELATION ANALYSIS with robust estimator of R ________########
################################################################################
library(robcor)

R_rob_pc <- robcor(pc)
R_rob_pm <- robcor(pm)

corrplot(R_rob_pc, type = "upper", tl.cex = 0.01)
corrplot(R_rob_pm, type = "upper", tl.cex = 0.01)

I <- NULL
J <- NULL

for (i in (2:p-1)) {
  for (j in (i+1):p) {
    if (abs(R_rob_pc[i,j]) >= 0.9) {
      I = c(I, i)
      J = c(J, j)
    }
  }
}

rob_corr_pc_0.9 = rbind(I,J)
# compare with corCorrelated_pc_0.9: 2,5,6,7,29;   9,10;   11,12,13,20,23,24;
#                                    15,34,35;   28,43;   41,46,47;  42,44

# groups: 2,5,6,7,28,29,42,43; 9,10;  3,11,12,13,20,23,24;    14,17;  
#         15,34,35,39,48,49;  19,21,25,26,27,32,33,36,50;  41,46,47;
# keep 7, remove --> 30

# using the results from the previous tests (choose variables which are 
# significant for more outcomes)
# -> 6 -> 9 -> 11 -> 14 -> 34 -> 27 -> 46

correlated_pc <- c(2,5,7,28,29,42,43,10,3,12,13,20,23,24,17,15,35,39,
                48,49,19,21,25,26,32,33,36,50,41,47)

I <- NULL
J <- NULL

for (i in (2:p-1)) {
  for (j in (i+1):p) {
    if (abs(R_rob_pm[i,j]) >= 0.9) {
      I = c(I, i)
      J = c(J, j)
    }
  }
}

rob_corr_pm_0.9 = rbind(I,J)
# compare with Correlated_pm_0.9: 2,5,6,7,29;   9,10;   11,12,13,20,23,24;
#                                 15,34,35;   28,43;   41,46,47;  42,44
# --> remove 20

# groups: 2,5,6,7,28,29,42,43;   9,10;  11,12,13,23,24,20;   14,17; 
#         15,34,35,48,49;   19,25,21;   26,27,32,33,36,50;  41,46,47;   
# --> keep 8, remove 27

# -> 6 -> 9 -> 11 -> 14 -> 34 -> 27 -> 46
# 19/25/21 not significant --> remove all 

correlated_pm <- c(2,5,7,28,29,42,43,10,12,13,23,24,20,17,15,35,48,49,
                   19,25,21,26,32,33,36,50,41,47)
remove <- intersect(correlated_pc,correlated_pm)

length(remove) # 28


############______ SAVE NEW DATASETS POST ROBUST CORRELATION ______############

pc_rob <- pc[,-remove]
pc_rob = c(pc.label, pc_rob)
write.table(pc_rob, file ="pc_post_robcor.txt")

pm_rob <- pm[,-remove]

for (i in 1:length(colnames(pm_rob))){
  colnames(pm_rob)[i] <- paste(colnames(pm_rob)[i],"margin",sep="_")
}

pm_rob = c(pm.label, pm_rob)
write.table(pm_rob, file ="pm_post_robcor.txt")

