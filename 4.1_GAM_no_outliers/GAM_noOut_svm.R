#GOAL: Find the best GAM model to predict the binary outcome "STATO VIVO/MORTO"

#Models considered:
# - Clinical + Portal Core
# - Clinical + Portal Core + Portal Margin

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

source("3_Logistic_GAMs/models_functions.R")
library(gam)
library(plsmselect)
library(purrr)
library(glmnet)
library(caret)
library(ROCR)
library(progress)
library(rsample)

seed <- 24091998
set.seed(seed)
#### Importing datasets ---- 
main <- read.table("Datasets/main.txt")
main.label <- main[,c(1,2)]
main <- main[,-c(1,2)]
#names(main)[names(main) == "Ca19.9â..55"] <- "Ca19.9.55" #for some windows pc


# Import radiomic's dataset
pc <-read.table("Datasets/pc_post_robcor.txt")
pm <-read.table("Datasets/pm_post_robcor.txt")
pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)]
pc <- scale(pc)
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)] 
pm <- scale(pm)

#Remove the outliers: 204,235,241,122
main <- main[-c(204,235,241,122),]
pc <- pc[-c(204,235,241,122),]
pm <- pm[-c(204,235,241,122),]


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
cv.lasso <- cv.glmnet(X,main1$STATO.VIVO.MORTO,nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0)

#Categorical variable significant: 
# HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + COMPLICANZE.SEVERE
# + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
# + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + 
# + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed

#Deleting from the dataset main the categorical variables not significant
which(coef.lasso == 0)


######## Model 1: Clinical Variables + Portal Core Radiomics -----
main_pc <- data.frame(main,pc)
main_pc <- na.omit(main_pc)
main_pc$PVE.preop <- NULL
main_pc$SEX <- NULL
main_pc$PRIMA.RESEZIONE <- NULL
main_pc$LINFOADENECTOMIA <- NULL
main_pc$ASSOCIATED.RESECTION <- NULL
main_pc$CIRROSI <- NULL
main_pc$SINGLE.NODULE <- NULL
main_pc$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc$INFILTRAZIONE.PERINEURALE <- NULL
main_pc$GRADING <- NULL

#Leaving out other outcomes
main_pc$RECIDIVA <- NULL
main_pc$OS..Days. <- NULL
main_pc$RFS..Days. <- NULL


#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c( "HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy",
           "RESEZIONE.VIA.BILIARE","COMPLICANZE.SEVERE","R.status","INVASIONE.VASCOLARE.MACROSCOPICA",
           "INVASIONE.VASCOLARE.MICROSCOPICA", "NODULI.SATELLITI","PATTERN","N","T.VIII.ed",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="STATO.VIVO.MORTO"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(STATO.VIVO.MORTO ~ AGE + Degenza,
                        data = main_pc,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
#Best model:
# AGE + Ca.19.9 + s(Degenza, df = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + s(GLRLM_SRHGE, df = 3) + NGLDM_Coarseness


#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza, k = 3) + 
                    CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + s(GLRLM_SRHGE, k = 3) + NGLDM_Coarseness
                  + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + COMPLICANZE.SEVERE
                  + RESEZIONE.VIA.BILIARE + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + NUMERO.LINFONODI.METASTATICI
                  + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)


# I one at a time remove:
#PATTERN + COMPLICANZE.SEVERE + Ca.19.9 + Ca19.9.55 + INVASIONE.VASCOLARE.MICROSCOPICA
#NGLDM_Coarseness + RESEZIONE.VIA.BILIARE + NODULI.SATELLITI +NUMERO.LINFONODI.METASTATICI
#R.status + HBV + s(GLRLM_SRHGE, k = 3) 
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE +  + s(Degenza, k = 3) + 
                    CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + HCV 
                    + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + 
                  + INVASIONE.VASCOLARE.MACROSCOPICA + NUMERO.LINFONODI.ASPORTATI 
                  + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.5134

final_model <- STATO.VIVO.MORTO ~ AGE +  + s(Degenza, k = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + HCV + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + INVASIONE.VASCOLARE.MACROSCOPICA + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed

#Cross-validation for the prediction
result <- our_loocv(main_pc, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.6941176

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.6506024

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.7356322

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
# 0.7889489






######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----

main_pc_pm <- data.frame(main[-257,],pc[-257,],pm) #Last line of Portal Margin is empty
main_pc_pm <- na.omit(main_pc_pm)
main_pc_pm$PVE.preop <- NULL
main_pc_pm$SEX <- NULL
main_pc_pm$PRIMA.RESEZIONE <- NULL
main_pc_pm$LINFOADENECTOMIA <- NULL
main_pc_pm$ASSOCIATED.RESECTION <- NULL
main_pc_pm$CIRROSI <- NULL
main_pc_pm$SINGLE.NODULE <- NULL
main_pc_pm$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pm$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pm$GRADING <- NULL

#Leaving out other outcomes
main_pc_pm$RECIDIVA <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL


#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pm
d_cont[,c( "HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy",
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
# AGE + Ca.19.9 + s(Degenza, df = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + CONVENTIONAL_HUSkewness + 
#   s(GLRLM_SRHGE, df = 3) + s(SHAPE_Compacity.onlyFor3DROI._margin, df = 2) 
# + s(GLRLM_SRHGE_margin, df = 3) + NGLDM_Contrast_margin + GLZLM_SZHGE_margin
#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza, k = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + CONVENTIONAL_HUSkewness + 
                    s(GLRLM_SRHGE, k = 3) + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) 
                  + s(GLRLM_SRHGE_margin, k = 3) + NGLDM_Contrast_margin + GLZLM_SZHGE_margin
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
#NODULI.SATELLITI + COMPLICANZE.SEVERE + Ca.19.9 + PATTERN + RESEZIONE.VIA.BILIARE
#CONVENTIONAL_HUSkewness + HBV + INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.METASTATICI
#GLZLM_SZHGE_margin
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE  + s(Degenza, k = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd
                  + s(GLRLM_SRHGE, k = 3) + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) 
                  + s(GLRLM_SRHGE_margin, k = 3) + NGLDM_Contrast_margin + 
                  + HCV +  + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy 
                  + R.status + INVASIONE.VASCOLARE.MACROSCOPICA
                  + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #0.62

final_model <- STATO.VIVO.MORTO ~ AGE  + s(Degenza, k = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd +
 s(GLRLM_SRHGE, k = 3) + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) +
 s(GLRLM_SRHGE_margin, k = 3) + NGLDM_Contrast_margin + 
 HCV +  + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy +
 R.status + INVASIONE.VASCOLARE.MACROSCOPICA +
 NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed

#Cross-validation for the prediction
result <- our_loocv(main_pc_pm, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.7455621

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.7228916

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #  0.7674419

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.8248809



