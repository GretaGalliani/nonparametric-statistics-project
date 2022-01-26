###############################################################################
######______ MULTIVARIATE OUTLIER DETECTION with ROBUST ESTIMATORS ______######
###############################################################################

#####---- importing libraries -----##### 
library(readxl)
library(robustbase)
library(psych)
library(MASS)
library(ellipse)
library(here)
library(DescTools)
library(knitr)
library(RobStatTM)

#####---- Dataset -----##### 
main <- read.table("Datasets/main.txt")
main.label <- main[,c(1,2)]
main <- main[,-c(1,2)]

pc_radiomics <-read.table("Datasets/pc_post_correlation.txt")
pm_radiomics <-read.table("Datasets/pm_post_correlation.txt")

pc_radiomics.label <- pc_radiomics[,c(1,2)]
pc_radiomics <- pc_radiomics[,-c(1,2)] 

pm_radiomics.label <- pm_radiomics[,c(1,2)]
pm_radiomics <- pm_radiomics[,-c(1,2)] 

# Standardize the variables
# pc_radiomics <- scale(pc_radiomics)
# pm_radiomics <- scale(pm_radiomics)

main_pc <- cbind(main,pc_radiomics)
main_pc_pm <- cbind(main[-c(258),],pc_radiomics[-c(258),], pm_radiomics)

# Remove NAs
main <- na.omit(main)
main_pc <- na.omit(main_pc)
main_pc_pm <- na.omit(main_pc_pm)

# expected outliers: 204,235,241,122 + 102 for portal margin

#####---- MCD estimator on CLINICAL VAR -----##### 
fit_MCD_main <- covMcd(x = main, alpha = .75, nsamp = "best", tolSolve = 1e-24)
ind_best_subset_main <- fit_MCD_main$best

# I believe this has some correlation problems (the warning is strange to interpret)
# since I want in the Regression part to understand also which are the outliers
# wrt this data, I will use only the variables I'd use in regression

lst <- c("AGE" ,"Ca.19.9","Degenza",
        "HBV" , "Ca19.9.55" , "CHEMIOTERAPIA.NEOADIUVANTE" , 
        "Major.Hepatectomy"  , "RESEZIONE.VIA.BILIARE" ,
        "ASSOCIATED.RESECTION" , "COMPLICANZE.SEVERE" , "CIRROSI" , 
        "SINGLE.NODULE" , "R.status" , "CHEMIOTERAPIA.ADIUVANTE" , 
        "NODULI.SATELLITI" , "GRADING" , "N" , "T.VIII.ed")

main_new <- NULL
for (i in 1:dim(main)[2]){
  if (names(main)[i] %in% lst){
    main_new <- cbind(main_new,main[,i])
  }
}
main_new <- data.frame(main_new)
# at this point it was still not working

# I try removing these two variables
#main_new <- main_new[,-c(13,14)]
# still not working

# last trial selecting some variables at hand
main_new <- data.frame(cbind(main$AGE,main$Ca.19.9,main$Degenza,main$Major.Hepatectomy,
                             main$R.status, main$CHEMIOTERAPIA.ADIUVANTE, main$NODULI.SATELLITI),
                             main$GRADING, main$N, main$T.VIII.ed)

fit_MCD_main_new <- covMcd(x = main_new, alpha = .75, nsamp = "best", tolSolve = 1e-29)
ind_best_subset_main_new <- fit_MCD_main_new$best

# IT WORKS! after a morning trying :(

plot(fit_MCD_main_new,classic=TRUE, tol = 1e-24)

# most extreme outliers from "main" are: 6, 39, 46, 48, 58, 63, 67, 77, 82, 85, 92, 98,112 141,163

#####---- MCD estimator on PORTAL CORE -----##### 

fit_MCD_pc <- covMcd(x = pc_radiomics, alpha = .75, nsamp = "best", tolSolve = 1e-24)
ind_best_subset_pc <- fit_MCD_pc$best

plot(fit_MCD_pc,classic=TRUE, tol = 1e-24)


#####---- MCD estimator on PORTAL MARGIN -----##### 

fit_MCD_pm <- covMcd(x = pm_radiomics, alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pm <- fit_MCD_pm$best

plot(fit_MCD_pm,classic=TRUE, tol = 1e-29)

# observation 102 is extremely high wrt the others --> 
# can't see what's going on with the rest, remove it

#####---- MCD estimator on PORTAL MARGIN without observation 102 -----##### 

fit_MCD_pm2<- covMcd(x = pm_radiomics[-102,], alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pm2 <- fit_MCD_pm2$best

plot(fit_MCD_pm2,classic=TRUE, tol = 1e-29)

# observation 110 is extremely high wrt the others --> 
# remove it and see what happens

#####---- MCD estimator on PORTAL MARGIN without observation 102,110 -----##### 
fit_MCD_pm3<- covMcd(x = pm_radiomics[-c(102,110),], alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pm3 <- fit_MCD_pm3$best

plot(fit_MCD_pm3,classic=TRUE, tol = 1e-29)


#####---- MCD estimator on PORTAL MARGIN without obs 102,110,109 -----##### 
fit_MCD_pm4<- covMcd(x = pm_radiomics[-c(102,109,110),], alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pm4 <- fit_MCD_pm4$best

plot(fit_MCD_pm4,classic=TRUE, tol = 1e-29)


#####---- MCD estimator on PORTAL MARGIN without obs 102,110,109,108 -----##### 
fit_MCD_pm5<- covMcd(x = pm_radiomics[-c(102,106,107,108,109,110),], alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pm5 <- fit_MCD_pm5$best

plot(fit_MCD_pm5,classic=TRUE, tol = 1e-29)

# every time there's an outlier which is extreme with respect to all others
# hard to see the trend of the plots

# ChiSquare QQ-plot difficult to interpret, but since it seems better than the 
# non-robust estimation, I'll give Robust Regression a try

### CONCLUSION: already individuated outliers from univariate perspective
# are included in the robust individuated outliers: 204,235,241,122,102 
# but there are many more!



############## RE-DO ALL ON POST ROBCOR DATASETS ##############
pc_radiomics <-read.table("Datasets/pc_post_robcor.txt")
pm_radiomics <-read.table("Datasets/pm_post_robcor.txt")

pc_radiomics.label <- pc_radiomics[,c(1,2)]
pc_radiomics <- pc_radiomics[,-c(1,2)] 

pm_radiomics.label <- pm_radiomics[,c(1,2)]
pm_radiomics <- pm_radiomics[,-c(1,2)] 


#####---- MCD estimator on PORTAL CORE -----##### 

fit_MCD_pc <- covMcd(x = pc_radiomics, 
                     alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pc <- fit_MCD_pc$best

plot(fit_MCD_pc,classic=TRUE, tol = 1e-29)

# SOME of the outlying observations, individuated by recursively
# removing the most outlying ones (plots are not readable)
# c(235,247,180,124,122,142,241,129,261,176,
# 106,77,141,145,178,138,190,224,79,173,
# 167,139,169,180,213,118,163,165,137,
# 176,208,117)


#####---- MCD estimator on PORTAL MARGIN -----##### 

fit_MCD_pm <- covMcd(x = pm_radiomics, alpha = .75, nsamp = "best", tolSolve = 1e-25)
ind_best_subset_pm <- fit_MCD_pm$best

plot(fit_MCD_pm,classic=TRUE, tol = 1e-29)




###############################################################################
###############____________ ROBUST REGRESSION for OS ____________##############
###############################################################################

main_pc$OS..Days.<- log(main_pc$OS..Days.)

data <- main_pc[,-which(names(main_pc)=='RFS..Days.')]

fit_lts <- ltsReg(OS..Days.~., alpha=.75,mcd=TRUE,data=data) 
# doesn't work because too many variables (even if #obs > 2*#covariates)

# I'll try to use the output of the gam.stepwise selection for OS.Days 
# and lasso selection of categorical variables done in "GAM_os.R"

fit_lts <- ltsReg(OS..Days.~ AGE + Ca.19.9 + Degenza + 
                    CONVENTIONAL_HUmax + GLZLM_SZHGE + GLZLM_LZLGE +
                    HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + 
                    Major.Hepatectomy  + RESEZIONE.VIA.BILIARE +
                    ASSOCIATED.RESECTION + COMPLICANZE.SEVERE + CIRROSI + 
                    SINGLE.NODULE + R.status + CHEMIOTERAPIA.ADIUVANTE + 
                    NODULI.SATELLITI + GRADING + N + T.VIII.ed, 
                    alpha=.75,mcd=TRUE,data=data) 
fit_lts <- ltsReg(OS..Days.~ AGE + Ca.19.9 + Degenza + 
                    CONVENTIONAL_HUmax + GLZLM_SZHGE + GLZLM_LZLGE +
                    Major.Hepatectomy  + R.status + CHEMIOTERAPIA.ADIUVANTE + 
                    NODULI.SATELLITI + GRADING + N + T.VIII.ed, 
                  alpha=.75,mcd=TRUE,data=data) 

plot(fit_lts) 

# covariance matrix results singular

# Normal qq-plot doesn't look so bad IF observations outside the line are those
# also identified as outliers by the robust estimators --> NOPE

# try to solve the problem with the warning






