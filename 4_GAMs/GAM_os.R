#GOAL: Find the best GAM model to predict the binary outcome "OS"
#Righe da levare 204,235,241,122
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

#Then compute via LOOCV the following metrics:
# - Accuracy
# - Sensitivity
# - Specificity
# - AUC

#Lastly some visualization tools on best model

#CONCLUSION:
#The best model is: Clinical + Portal Core + Portal Margin

source("3_Logistic_GAMs/models_functions.R")
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

main$OS..Days.<- log(main$OS..Days.)

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
cv.lasso <- cv.glmnet(X,main1$OS..Days.,nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0)  
#Categorical variable significant: 
# HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy  + RESEZIONE.VIA.BILIARE +
# ASSOCIATED.RESECTION + COMPLICANZE.SEVERE + CIRROSI + SINGLE.NODULE + R.status +
# CHEMIOTERAPIA.ADIUVANTE + NODULI.SATELLITI + GRADING + N + T.VIII.ed

#Deleting from the dataset main the categorical variables not significant
which(coef.lasso == 0)

main$SEX <- NULL
main$HCV <- NULL
main$PVE.preop <- NULL
main$PRIMA.RESEZIONE <- NULL
main$LINFOADENECTOMIA <- NULL
main$INFILTRAZIONE.PERINEURALE <- NULL
main$PATTERN <- NULL
main$NUMERO.LINFONODI.METASTATICI <- NULL
main$NUMERO.LINFONODI.ASPORTATI <-NULL
main$NUMERO..NO.SATELLITI. <- NULL


######## Model 1: Clinical Variables + Portal Core Radiomics -----

main_pc <- data.frame(main,pc)
#Leaving out other outcomes
main_pc$RECIDIVA <- NULL
main_pc$STATO.VIVO.MORTO <- NULL
main_pc$RFS..Days. <- NULL
main_pc <- na.omit(main_pc)

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c("HBV", "Ca19.9.55", "CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy","RESEZIONE.VIA.BILIARE" 
            ,"ASSOCIATED.RESECTION", "COMPLICANZE.SEVERE","CIRROSI","SINGLE.NODULE","R.status" ,
            "CHEMIOTERAPIA.ADIUVANTE" ,"NODULI.SATELLITI" ,"GRADING" ,"N" ,"T.VIII.ed")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="OS..Days."), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(OS..Days. ~ AGE ,
                        data = main_pc,
                        family = gaussian,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
#Best model:
# AGE + Ca.19.9 + s(Degenza, df = 3) + s(CONVENTIONAL_HUmax, df = 2) + GLZLM_SZHGE + GLZLM_LZLGE

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(OS..Days. ~ AGE + Ca.19.9 + s(Degenza, k = 3) + s(CONVENTIONAL_HUmax, k = 3) 
                  + GLZLM_SZHGE + GLZLM_LZLGE
                  + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy  + RESEZIONE.VIA.BILIARE +
                    ASSOCIATED.RESECTION + COMPLICANZE.SEVERE + CIRROSI + SINGLE.NODULE + R.status +
                    CHEMIOTERAPIA.ADIUVANTE + NODULI.SATELLITI + GRADING + N + T.VIII.ed,
                  data = main_pc,
                  method='REML',
                  family = gaussian,
                  select = T)
summary(fit1)

#Let us look at the residuals
hist(fit1$residuals)
qqnorm(fit1$residuals)
plot(fit1)
shapiro.test(fit1$residuals)
library(faraway)
halfnorm(fit1$residuals)
# I one at a time remove:


# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + s(Degenza, k = 3) + CONVENTIONAL_HUKurtosis +
                    SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + SHAPE_Surface.mm2..onlyFor3DROI. +
                    GLCM_Contrast..Variance. + GLRLM_SRLGE + HCV +  + Ca19.9.55 + 
                    CHEMIOTERAPIA.NEOADIUVANTE  + COMPLICANZE.SEVERE + NUMERO..NO.SATELLITI. 
                  + N + T.VIII.ed,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.23

final_model <- STATO.VIVO.MORTO ~ AGE +  + s(Degenza, k = 3) + CONVENTIONAL_HUKurtosis +
  SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + SHAPE_Surface.mm2..onlyFor3DROI. +
  GLCM_Contrast..Variance. + GLRLM_SRLGE + HCV +  + Ca19.9.55 + 
  CHEMIOTERAPIA.NEOADIUVANTE  + COMPLICANZE.SEVERE + NUMERO..NO.SATELLITI. + N + T.VIII.ed

#Cross-validation for the prediction
result <- our_loocv(main_pc, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.6952381

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.6545455

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[1,2])
specificity # 0.627451

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
# 0.7533636






######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----

main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
#Leaving out other outcomes
main_pc_pm$RECIDIVA <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL
main_pc_pm <- na.omit(main_pc_pm)

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pm
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy",
           "RESEZIONE.VIA.BILIARE","COMPLICANZE.SEVERE","R.status","INVASIONE.VASCOLARE.MACROSCOPICA",
           "INVASIONE.VASCOLARE.MICROSCOPICA", "NODULI.SATELLITI","PATTERN","N","T.VIII.ed",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="STATO.VIVO.MORTO"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(STATO.VIVO.MORTO ~ AGE + Degenza,
                        data = main_pc_pm,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm)
summary(step)
#Best model:
# AGE + Ca.19.9 + s(Degenza,df = 3) + s(CONVENTIONAL_HUmin, df = 2) + SHAPE_Volume.mL. + 
#   GLCM_Contrast..Variance. + GLRLM_SRHGE + NGLDM_Busyness + 
#   s(CONVENTIONAL_HUmax_margin, df = 3) + SHAPE_Surface.mm2..onlyFor3DROI._margin + 
#   SHAPE_Compacity.onlyFor3DROI._margin + GLRLM_LGRE_margin + 
#   GLZLM_SZE_margin + GLZLM_SZHGE_margin + GLZLM_LZLGE_margin

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza,k = 3) + s(CONVENTIONAL_HUmin, k = 3) + SHAPE_Volume.mL. + 
                    GLCM_Contrast..Variance. + GLRLM_SRHGE + NGLDM_Busyness + 
                    s(CONVENTIONAL_HUmax_margin, k = 3) + SHAPE_Surface.mm2..onlyFor3DROI._margin + 
                    SHAPE_Compacity.onlyFor3DROI._margin + GLRLM_LGRE_margin + 
                    GLZLM_SZE_margin + GLZLM_SZHGE_margin + GLZLM_LZLGE_margin
                  + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + COMPLICANZE.SEVERE
                  + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + NUMERO.LINFONODI.METASTATICI
                  + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#NODULI.SATELLITI,PATTERN,Ca.19.9,NGLDM_Busyness,SHAPE_Volume.mL.,GLZLM_LZLGE_margin
#NUMERO..NO.SATELLITI.,NUMERO.LINFONODI.METASTATICI,GLRLM_LGRE_margin,GLRLM_SRHGE
#CHEMIOTERAPIA.NEOADIUVANTE,COMPLICANZE.SEVERE,HBV,GLZLM_SZE_margin
#SHAPE_Surface.mm2..onlyFor3DROI._margin, s(CONVENTIONAL_HUmax_margin, k = 3)
#s(Degenza,k = 3)
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + s(CONVENTIONAL_HUmin, k = 3)  + 
                    GLCM_Contrast..Variance. + 
                    SHAPE_Compacity.onlyFor3DROI._margin + GLZLM_SZHGE_margin 
                  + HCV + Ca19.9.55 + Major.Hepatectomy + 
                    + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA+ NUMERO.LINFONODI.ASPORTATI  
                  + N + T.VIII.ed ,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.20

final_model <-STATO.VIVO.MORTO ~ AGE + s(CONVENTIONAL_HUmin, k = 3)  + 
  GLCM_Contrast..Variance.+ 
  SHAPE_Compacity.onlyFor3DROI._margin + GLZLM_SZHGE_margin + HCV + Ca19.9.55 + Major.Hepatectomy + 
  RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA+
  INVASIONE.VASCOLARE.MICROSCOPICA+ NUMERO.LINFONODI.ASPORTATI  +
  N + T.VIII.ed
#Cross-validation for the prediction
result <- our_loocv(main_pc_pm, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.708134

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.7181818

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[1,2])
specificity # 0.69

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.8048669



######## Model 3: Clinical Variables + PC's (Portal Core) -----

main_pc_pca <- data.frame(main,pc_pca) #Last line of Portal Margin is empty
#Leaving out other outcomes
main_pc_pca$RECIDIVA <- NULL
main_pc_pca$OS..Days. <- NULL
main_pc_pca$RFS..Days. <- NULL
main_pc_pca <- na.omit(main_pc_pca)

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pca
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy",
           "RESEZIONE.VIA.BILIARE","COMPLICANZE.SEVERE","R.status","INVASIONE.VASCOLARE.MACROSCOPICA",
           "INVASIONE.VASCOLARE.MICROSCOPICA", "NODULI.SATELLITI","PATTERN","N","T.VIII.ed",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="STATO.VIVO.MORTO"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(STATO.VIVO.MORTO ~ AGE + Degenza,
                        data = main_pc_pca,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pca)
summary(step)
#Best model:
# AGE + Ca.19.9 + Degenza + Comp.1

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + Degenza + Comp.1
                  +HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + COMPLICANZE.SEVERE
                  + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + NUMERO.LINFONODI.METASTATICI
                  + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#Ca.19.9,PATTERN,R.status,NODULI.SATELLITI,NUMERO.LINFONODI.METASTATICI
#INVASIONE.VASCOLARE.MICROSCOPICA,Major.Hepatectomy,INVASIONE.VASCOLARE.MACROSCOPICA
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE +  + Degenza + Comp.1
                  +HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE +  + COMPLICANZE.SEVERE
                  + RESEZIONE.VIA.BILIARE + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.20

final_model <-STATO.VIVO.MORTO ~ AGE +  + Degenza + Comp.1 + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE +  + COMPLICANZE.SEVERE + RESEZIONE.VIA.BILIARE + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed

#Cross-validation for the prediction
result <- our_loocv(main_pc_pca, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.6714286

# Sensitivity 
specificity <- tab[2,2]/(tab[2,2] + tab[1,2])
specificity # 0.6727273

# Sensitivity 
specificity <- tab[1,1]/(tab[1,1]+tab[1,2])
specificity # 0.6504854

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7533636




######## Model 4: Clinical Variables + PC's (Portal Core & Margin) -----

main_pc_pm_pca <- data.frame(main[-261,],pc_pca[-261,],pm_pca)
#Leaving out other outcomes
main_pc_pm_pca$RECIDIVA <- NULL
main_pc_pm_pca$OS..Days. <- NULL
main_pc_pm_pca$RFS..Days. <- NULL
main_pc_pm_pca <- na.omit(main_pc_pm_pca)

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pm_pca
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy",
           "RESEZIONE.VIA.BILIARE","COMPLICANZE.SEVERE","R.status","INVASIONE.VASCOLARE.MACROSCOPICA",
           "INVASIONE.VASCOLARE.MICROSCOPICA", "NODULI.SATELLITI","PATTERN","N","T.VIII.ed",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="STATO.VIVO.MORTO"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(STATO.VIVO.MORTO ~ AGE + Degenza,
                        data = main_pc_pm_pca,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm_pca)
summary(step)
#Best model:
# AGE + Ca.19.9 + s(Degenza, df = 3) + Comp.1 + Comp.2 + Comp.4_margin

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza, k = 3) + Comp.1 + Comp.2 + Comp.4_margin
                  +HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + COMPLICANZE.SEVERE
                  + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + NUMERO.LINFONODI.METASTATICI
                  + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#NUMERO..NO.SATELLITI.,PATTERN,Ca.19.9,Comp.2,R.status,Major.Hepatectomy
#INVASIONE.VASCOLARE.MICROSCOPICA,N,INVASIONE.VASCOLARE.MACROSCOPICA
#RESEZIONE.VIA.BILIARE,Comp.4_margin,HBV,HCV,NODULI.SATELLITI
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE +  + s(Degenza, k = 3) + Comp.1  + 
                    +Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE  + COMPLICANZE.SEVERE
                  +NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + T.VIII.ed,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.20

final_model <-AGE +  + s(Degenza, k = 3) + Comp.1 +Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE  + COMPLICANZE.SEVERE +NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + + T.VIII.ed


#Cross-validation for the prediction
result <- our_loocv(main_pc_pm_pca, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.6698565

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.6727273

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[1,2])
specificity # 0.6504854

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7533636

















