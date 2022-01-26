#GOAL: Find the best GAM model to predict the binary outcome "INVASIONE VASCOLARE MICROSCOPICA"

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

main$GRADING <- ifelse(main$GRADING == 3,"3","1-2")


#### Feature selection on categorical variables: LASSO ----
main1 <- na.omit(main)
X <- model.matrix(~ SEX + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE
                  + PVE.preop + PRIMA.RESEZIONE+ ASSOCIATED.RESECTION
                  + CIRROSI + SINGLE.NODULE +  
                 + INFILTRAZIONE.PERINEURALE + PATTERN + NUMERO..NO.SATELLITI.+ GRADING,
                  data=main1)[,-1]

lambda.grid <- 10^seq(5,-3,length=100) #defining lambda grid
cv.lasso <- cv.glmnet(X,main1$INVASIONE.VASCOLARE.MICROSCOPICA, nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0)  
#Categorical variable significant: 
#SEX + PVE.preop + ASSOCIATED.RESECTION + INFILTRAZIONE.PERINEURALE + GRADING

#Deleting from the dataset main the categorical variables not significant
which(coef.lasso == 0)
    

######## Model 1: Clinical Variables + Portal Core Radiomics -----
main_pc <- data.frame(main,pc)
main_pc <- na.omit(main_pc)
#Delete pre surgery variables
main_pc$Degenza <- NULL
main_pc$Major.Hepatectomy <- NULL
main_pc$RESEZIONE.VIA.BILIARE <- NULL
main_pc$LINFOADENECTOMIA <- NULL
main_pc$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc$N <-NULL
main_pc$R.status <-NULL
main_pc$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc$T.VIII.ed <-NULL          
main_pc$NODULI.SATELLITI <-NULL    
main_pc$COMPLICANZE.SEVERE <-NULL
main_pc$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc$STATO.VIVO.MORTO <- NULL
main_pc$RECIDIVA <- NULL
main_pc$OS..Days. <- NULL
main_pc$RFS..Days. <- NULL

#Removing variables  discarded by Lasso
main_pc$Ca19.9.55 <- NULL
main_pc$HCV <- NULL
main_pc$HBV <- NULL
main_pc$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc$PRIMA.RESEZIONE <- NULL
main_pc$CIRROSI <- NULL
main_pc$SINGLE.NODULE <- NULL
main_pc$PATTERN <- NULL
main_pc$NUMERO..NO.SATELLITI. <- NULL

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c( "SEX","PVE.preop", "ASSOCIATED.RESECTION", "INFILTRAZIONE.PERINEURALE",
           "GRADING")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="INVASIONE.VASCOLARE.MICROSCOPICA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ AGE,
                        data = main_pc,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)

step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
#Best model:
 # Ca.19.9 + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + s(DISCRETIZED_HISTO_Entropy_log10, df = 3) + GLRLM_SRHGE

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~  Ca.19.9 + CONVENTIONAL_HUmin + 
                    CONVENTIONAL_HUstd + s(DISCRETIZED_HISTO_Entropy_log10, k = 3) + GLRLM_SRHGE+
                    SEX + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# PVE.preop + Ca.19.9 + ASSOCIATED.RESECTION + CONVENTIONAL_HUstd + SEX
# --->FINAL FIT:
fit2 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ CONVENTIONAL_HUmin + 
                     s(DISCRETIZED_HISTO_Entropy_log10, k = 3) + GLRLM_SRHGE +
                     INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") 

final_model <- INVASIONE.VASCOLARE.MICROSCOPICA ~ CONVENTIONAL_HUmin + s(DISCRETIZED_HISTO_Entropy_log10, k = 3) + GLRLM_SRHGE+INFILTRAZIONE.PERINEURALE + GRADING

#Cross-validation for the prediction
result <- our_loocv(main_pc, final_model, "INVASIONE.VASCOLARE.MICROSCOPICA")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.7294118

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.8434783

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.5087719

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
# 0.7505046



######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----
main_pc_pm <- data.frame(main[-257,],pc[-257,],pm) #Last line of Portal Margin is empty
main_pc_pm <- na.omit(main_pc_pm)
#Delete pre surgery variables
main_pc_pm$Degenza <- NULL
main_pc_pm$Major.Hepatectomy <- NULL
main_pc_pm$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm$LINFOADENECTOMIA <- NULL
main_pc_pm$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc_pm$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc_pm$N <-NULL
main_pc_pm$R.status <-NULL
main_pc_pm$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc_pm$T.VIII.ed <-NULL          
main_pc_pm$NODULI.SATELLITI <-NULL    
main_pc_pm$COMPLICANZE.SEVERE <-NULL
main_pc_pm$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc_pm$STATO.VIVO.MORTO <- NULL
main_pc_pm$RECIDIVA <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL


#Removing variables discarded by Lasso
main_pc_pm$Ca19.9.55 <- NULL
main_pc_pm$HCV <- NULL
main_pc_pm$HBV <- NULL
main_pc_pm$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pm$PRIMA.RESEZIONE <- NULL
main_pc_pm$CIRROSI <- NULL
main_pc_pm$SINGLE.NODULE <- NULL
main_pc_pm$PATTERN <- NULL
main_pc_pm$NUMERO..NO.SATELLITI. <- NULL

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pm
d_cont[,c("SEX","PVE.preop", "ASSOCIATED.RESECTION", "INFILTRAZIONE.PERINEURALE",
          "GRADING")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="INVASIONE.VASCOLARE.MICROSCOPICA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ AGE,
                        data = main_pc_pm,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm)
summary(step)
#Best model:
# Ca.19.9 + GLRLM_SRHGE + CONVENTIONAL_HUstd_margin + 
#   s(SHAPE_Sphericity.onlyFor3DROI.._margin, df = 3) + s(NGLDM_Contrast_margin, df = 2)

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ Ca.19.9 + GLRLM_SRHGE + CONVENTIONAL_HUstd_margin + 
                    s(SHAPE_Sphericity.onlyFor3DROI.._margin, k = 3) + s(NGLDM_Contrast_margin, k = 2) +
                    SEX + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#PVE.preop + Ca.19.9 + ASSOCIATED.RESECTION  +SEX +s(NGLDM_Contrast_margin, k = 3)
#s(SHAPE_Sphericity.onlyFor3DROI.._margin, k = 3) 

# --->FINAL FIT:
fit2 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~  GLRLM_SRHGE + CONVENTIONAL_HUstd_margin + 
                    INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") 

final_model <- INVASIONE.VASCOLARE.MICROSCOPICA ~  GLRLM_SRHGE + CONVENTIONAL_HUstd_margin + 
  INFILTRAZIONE.PERINEURALE + GRADING

#Cross-validation for the prediction
result <- our_loocv(main_pc_pm, final_model, "INVASIONE.VASCOLARE.MICROSCOPICA")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.739645

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.8230088

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #  0.5714286

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7594817

