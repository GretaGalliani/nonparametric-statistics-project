#GOAL: Find the best GAM model to predict the outcome "RFS"

#Models considered:
# - Clinical + Portal Core 
# - Clinical + Portal Core + Portal Margin
# - Clinical + Portal Core (PCA)
# - Clinical + Portal Core (PCA) + Portal Margin (PCA)

#To reduce the high starting dimensionality we use:
# - LASSO (with glmnet) to reduce categorical variables
# - STEPWISE (step.Gam) to reduce continuous variables

#After that eventually remove manually other covariates that are not significant
# in predicting the outcome

library(gam)
library(plsmselect)
library(purrr)
library(glmnet)
library(caret)
library(progress)
library(rsample)

seed <- 24091998
set.seed(seed)
#### Importing datasets ---- 
main <- read.table("Datasets/main.txt")
main.label <- main[,c(1,2)]
main <- main[,-c(1,2)]
names(main)[names(main) == "Ca19.9â..55"] <- "Ca19.9.55" #for some windows pc

# Import radiomic's dataset
pc <-read.table("Datasets/pc_post_correlation.txt")
pm <-read.table("Datasets/pm_post_correlation.txt")
pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)] 
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)] 

# Import PC's of radiomic's features
pc_pca <-read.table("Datasets/pc_post_pca.txt")
pm_pca <-read.table("Datasets/pm_post_pca.txt")
#Deleting labels
pc_pca <- pc_pca[,-c(1,2)]
pm_pca <- pm_pca[,-c(1,2)]

main$OS..Days. <- NULL
main$STATO.VIVO.MORTO <- NULL
main$RECIDIVA <- NULL

main$RFS..Days. <- log(main$RFS..Days.)

#### Feature selection on categorical variables: LASSO ----
main1 <- na.omit(main)
X <- model.matrix(~ SEX + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE
                  + PVE.preop + Major.Hepatectomy + PRIMA.RESEZIONE
                  + RESEZIONE.VIA.BILIARE + LINFOADENECTOMIA + ASSOCIATED.RESECTION
                  + COMPLICANZE.SEVERE +CIRROSI+
                  + SINGLE.NODULE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA + CHEMIOTERAPIA.ADIUVANTE 
                  + INFILTRAZIONE.PERINEURALE + NODULI.SATELLITI + PATTERN 
                  + NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI 
                  + NUMERO..NO.SATELLITI. + GRADING + N + T.VIII.ed,
                  data=main1)[,-1]

lambda.grid <- 10^seq(5,-3,length=100) #defining lambda grid
cv.lasso <- cv.glmnet(X,main1$RFS..Days.,nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0)
# HBV + Ca19.9.55 + Major.Hepatectomy + RESEZIONE.VIA.BILIARE + COMPLICANZE.SEVERE + 
# CIRROSI + SINGLE.NODULE + R.status + NODULI.SATELLITI +
# NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. +
# N + T.VIII.ed

#removing the others
main$SEX <- NULL
main$HCV <- NULL
main$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main$PVE.preop <- NULL
main$PRIMA.RESEZIONE <- NULL
main$LINFOADENECTOMIA <- NULL
main$ASSOCIATED.RESECTION <- NULL
main$INVASIONE.VASCOLARE.MACROSCOPICA <- NULL
main$INVASIONE.VASCOLARE.MICROSCOPICA <- NULL
main$CHEMIOTERAPIA.ADIUVANTE <- NULL
main$INFILTRAZIONE.PERINEURALE <- NULL
main$PATTERN <- NULL
main$GRADING <- NULL

main_pc <- data.frame(main,pc)
main_pc <- na.omit(main_pc)


######## Model 1: Clinical Variables + Portal Core Radiomics -----

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c("HBV","Ca19.9.55","Major.Hepatectomy","RESEZIONE.VIA.BILIARE","COMPLICANZE.SEVERE",
            "CIRROSI","SINGLE.NODULE","R.status","NODULI.SATELLITI","NUMERO.LINFONODI.METASTATICI",
            "NUMERO.LINFONODI.ASPORTATI","NUMERO..NO.SATELLITI.","N","T.VIII.ed")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="RFS..Days."), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(RFS..Days. ~ AGE + Degenza,
                        data = main_pc,
                        family = gaussian,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
 
#Best model:
# AGE + Ca.19.9 + s(Degenza, df = 3) + CONVENTIONAL_HUKurtosis +
# SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + SHAPE_Surface.mm2..onlyFor3DROI. +
# GLCM_Contrast..Variance. + GLRLM_SRLGE + NGLDM_Coarseness

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza, k = 3) + CONVENTIONAL_HUKurtosis +
                    SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + SHAPE_Surface.mm2..onlyFor3DROI. +
                    GLCM_Contrast..Variance. + GLRLM_SRLGE + NGLDM_Coarseness 
                  + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + COMPLICANZE.SEVERE
                  + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + NUMERO.LINFONODI.METASTATICI
                  + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

