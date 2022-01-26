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
#names(main)[names(main) == "Ca19.9â..55"] <- "Ca19.9.55" #for some windows pc

# Import radiomic's dataset
pc <-read.table("Datasets/pc_post_correlation.txt")
pm <-read.table("Datasets/pm_post_correlation.txt")
pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)] 
pc <- scale(pc)
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)]
pm <- scale(pm)

# Import PC's of radiomic's features
pc_pca <-read.table("Datasets/pc_post_pca.txt")
pm_pca <-read.table("Datasets/pm_post_pca.txt")
#Deleting labels
pc_pca <- pc_pca[,-c(1,2)]
pm_pca <- pm_pca[,-c(1,2)]

main$GRADING <- ifelse(main$GRADING == 3,"3","1-2") #grouping grading

#MODEL 1:
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

#MODEL 2
main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
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

#MODEL 3
main_pc_pca <- data.frame(main,pc_pca) #Last line of Portal Margin is empty
main_pc_pca <- na.omit(main_pc_pca)
#Delete pre surgery variables
main_pc_pca$Degenza <- NULL
main_pc_pca$Major.Hepatectomy <- NULL
main_pc_pca$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pca$LINFOADENECTOMIA <- NULL
main_pc_pca$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc_pca$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc_pca$N <-NULL
main_pc_pca$R.status <-NULL
main_pc_pca$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc_pca$T.VIII.ed <-NULL          
main_pc_pca$NODULI.SATELLITI <-NULL    
main_pc_pca$COMPLICANZE.SEVERE <-NULL
main_pc_pca$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc_pca$STATO.VIVO.MORTO <- NULL
main_pc_pca$RECIDIVA <- NULL
main_pc_pca$OS..Days. <- NULL
main_pc_pca$RFS..Days. <- NULL

#MODEL 4
main_pc_pm_pca <- data.frame(main[-261,],pc_pca[-261,],pm_pca)
main_pc_pm_pca <- na.omit(main_pc_pm_pca)
#Delete pre surgery variables
main_pc_pm_pca$Degenza <- NULL
main_pc_pm_pca$Major.Hepatectomy <- NULL
main_pc_pm_pca$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm_pca$LINFOADENECTOMIA <- NULL
main_pc_pm_pca$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc_pm_pca$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc_pm_pca$N <-NULL
main_pc_pm_pca$R.status <-NULL
main_pc_pm_pca$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc_pm_pca$T.VIII.ed <-NULL          
main_pc_pm_pca$NODULI.SATELLITI <-NULL    
main_pc_pm_pca$COMPLICANZE.SEVERE <-NULL
main_pc_pm_pca$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc_pm_pca$STATO.VIVO.MORTO <- NULL
main_pc_pm_pca$RECIDIVA <- NULL
main_pc_pm_pca$OS..Days. <- NULL
main_pc_pm_pca$RFS..Days. <- NULL

#### Feature selection on categorical variables: LASSO ----
main1 <- na.omit(main)
X <- model.matrix(~ SEX + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE
                  + PVE.preop + PRIMA.RESEZIONE+ ASSOCIATED.RESECTION
                  + CIRROSI + SINGLE.NODULE +  
                  + INFILTRAZIONE.PERINEURALE + PATTERN + NUMERO..NO.SATELLITI. 
                  + GRADING,
                  data=main1)[,-1]

lambda.grid <- 10^seq(5,-3,length=100) #defining lambda grid
cv.lasso <- cv.glmnet(X,main1$INVASIONE.VASCOLARE.MICROSCOPICA, nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0)  
#Categorical variable significant: 
# Ca19.9.55 + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING

#Deleting from the dataset main the categorical variables not significant
which(coef.lasso == 0)

######## Model 1: Clinical Variables + Portal Core Radiomics -----
#removing variables discarded by LASSO
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
# Ca.19.9 + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + s(DISCRETIZED_HISTO_Entropy_log10, df = 3) + 
#   GLRLM_SRHGE

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ Ca.19.9 + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + 
                    s(DISCRETIZED_HISTO_Entropy_log10, k = 3) + GLRLM_SRHGE
                  +SEX + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# PVE.preop + Ca.19.9 + ASSOCIATED.RESECTION + CONVENTIONAL_HUstd + SEX
# --->FINAL FIT:
fit2 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ CONVENTIONAL_HUmin + s(DISCRETIZED_HISTO_Entropy_log10, k = 3) + GLRLM_SRHGE
                  +INFILTRAZIONE.PERINEURALE + GRADING,
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
accuracy # 0.7325581

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
# 0.7554539



######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----
#removing variables discarded by LASSO
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
# Ca.19.9 + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd + s(DISCRETIZED_HISTO_Entropy_log10,df = 3) + GLRLM_SRHGE

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ Ca.19.9 + CONVENTIONAL_HUmin + CONVENTIONAL_HUstd 
                  + s(DISCRETIZED_HISTO_Entropy_log10,k = 3) + GLRLM_SRHGE
                  +SEX + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#PVE.preop + Ca.19.9 + ASSOCIATED.RESECTION + CONVENTIONAL_HUstd +SEX

# --->FINAL FIT:
fit2 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ CONVENTIONAL_HUmin + 
                    s(DISCRETIZED_HISTO_Entropy_log10,k = 3) + GLRLM_SRHGE
                  +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") 

final_model <- INVASIONE.VASCOLARE.MICROSCOPICA ~ CONVENTIONAL_HUmin + 
  s(DISCRETIZED_HISTO_Entropy_log10,k = 3) + GLRLM_SRHGE + INFILTRAZIONE.PERINEURALE + GRADING
  
#Cross-validation for the prediction
result <- our_loocv(main_pc_pm, final_model, "INVASIONE.VASCOLARE.MICROSCOPICA")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.7426901

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.8608696

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #  0.5

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7608696



######## Model 3: Clinical Variables + PC's (Portal Core) -----
#removing variables discarded by LASSO
main_pc_pca$Ca19.9.55 <- NULL
main_pc_pca$HCV <- NULL
main_pc_pca$HBV <- NULL
main_pc_pca$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pca$PRIMA.RESEZIONE <- NULL
main_pc_pca$CIRROSI <- NULL
main_pc_pca$SINGLE.NODULE <- NULL
main_pc_pca$PATTERN <- NULL
main_pc_pca$NUMERO..NO.SATELLITI. <- NULL

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pca
d_cont[,c("SEX","PVE.preop", "ASSOCIATED.RESECTION", "INFILTRAZIONE.PERINEURALE",
          "GRADING")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="INVASIONE.VASCOLARE.MICROSCOPICA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ AGE,
                        data = main_pc_pca,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)

step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pca)
summary(step)
#Best model:
# Ca.19.9 + s(Comp.1, df = 3)

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ Ca.19.9 + s(Comp.1, k = 3)
                  +SEX + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#PVE.preop,ASSOCIATED.RESECTION,Ca.19.9, ASSOCIATED.RESECTION, SEX

# --->FINAL FIT:
fit2 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ s(Comp.1, k = 3)+ INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") 

final_model <-INVASIONE.VASCOLARE.MICROSCOPICA ~ s(Comp.1, k = 3)+ INFILTRAZIONE.PERINEURALE + GRADING

#Cross-validation for the prediction
result <- our_loocv(main_pc_pca, final_model, "INVASIONE.VASCOLARE.MICROSCOPICA")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.744186

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.8608696

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #  0.5087719

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7491991




######## Model 4: Clinical Variables + PC's (Portal Core & Margin) -----
#removing variables discarded by LASSO
main_pc_pm_pca$Ca19.9.55 <- NULL
main_pc_pm_pca$HCV <- NULL
main_pc_pm_pca$HBV <- NULL
main_pc_pm_pca$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pm_pca$PRIMA.RESEZIONE <- NULL
main_pc_pm_pca$CIRROSI <- NULL
main_pc_pm_pca$SINGLE.NODULE <- NULL
main_pc_pm_pca$PATTERN <- NULL
main_pc_pm_pca$NUMERO..NO.SATELLITI. <- NULL


#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc_pm_pca
d_cont[,c("SEX","PVE.preop", "ASSOCIATED.RESECTION", "INFILTRAZIONE.PERINEURALE",
          "GRADING")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="INVASIONE.VASCOLARE.MICROSCOPICA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ AGE,
                        data = main_pc_pm_pca,
                        family = binomial,
                        na.action=na.omit)
summary(start_model)

step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm_pca)
summary(step)
#Best model:
# Ca.19.9 + s(Comp.1, df = 3) + Comp.1_margin

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ Ca.19.9 + s(Comp.1, k = 3) + Comp.1_margin
                  +SEX + PVE.preop + ASSOCIATED.RESECTION +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
#PVE.preop,ASSOCIATED.RESECTION,Ca.19.9
# --->FINAL FIT:
fit2 <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ s(Comp.1, k = 3) + Comp.1_margin
                  +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") 

final_model <- INVASIONE.VASCOLARE.MICROSCOPICA ~ s(Comp.1, k = 3) + Comp.1_margin+INFILTRAZIONE.PERINEURALE + GRADING

#Cross-validation for the prediction
result <- our_loocv(main_pc_pm_pca, final_model, "INVASIONE.VASCOLARE.MICROSCOPICA")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.7368421

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.8434783

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.5178571

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7464286



















######## ODDS RATIO FOR THE BEST MODEL (Model 2) ----

main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
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

#removing variables discarded by LASSO
main_pc_pm$Ca19.9.55 <- NULL
main_pc_pm$HCV <- NULL
main_pc_pm$HBV <- NULL
main_pc_pm$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pm$PRIMA.RESEZIONE <- NULL
main_pc_pm$CIRROSI <- NULL
main_pc_pm$SINGLE.NODULE <- NULL
main_pc_pm$PATTERN <- NULL
main_pc_pm$NUMERO..NO.SATELLITI. <- NULL

fit <- mgcv::gam(INVASIONE.VASCOLARE.MICROSCOPICA ~ CONVENTIONAL_HUmin + 
                    s(DISCRETIZED_HISTO_Entropy_log10,k = 3) + GLRLM_SRHGE
                  +INFILTRAZIONE.PERINEURALE + GRADING,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit)


#ODDS RATIO
alpha = 0.05
qalpha = qnorm( 1 - alpha/2 )

INF <- NULL
CENT <- NULL
SUP <- NULL
names <- NULL

#CONVENTIONAL_HUmin 
IC.sup = exp(fit$coefficients["CONVENTIONAL_HUmin"] + qalpha * summary(fit)$se["CONVENTIONAL_HUmin"])
IC.or = exp(fit$coefficients["CONVENTIONAL_HUmin"])
IC.inf = exp(fit$coefficients["CONVENTIONAL_HUmin"] - qalpha * summary(fit)$se["CONVENTIONAL_HUmin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"CONVENTIONAL_HUmin")

#GLRLM_SRHGE 
IC.sup = exp(fit$coefficients["GLRLM_SRHGE"] + qalpha * summary(fit)$se["GLRLM_SRHGE"])
IC.or = exp(fit$coefficients["GLRLM_SRHGE"])
IC.inf = exp(fit$coefficients["GLRLM_SRHGE"] - qalpha * summary(fit)$se["GLRLM_SRHGE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLRLM_SRHGE")

#INFILTRAZIONE.PERINEURALE 
IC.sup = exp(fit$coefficients["INFILTRAZIONE.PERINEURALE"] + qalpha * summary(fit)$se["INFILTRAZIONE.PERINEURALE"])
IC.or = exp(fit$coefficients["INFILTRAZIONE.PERINEURALE"])
IC.inf = exp(fit$coefficients["INFILTRAZIONE.PERINEURALE"] - qalpha * summary(fit)$se["INFILTRAZIONE.PERINEURALE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"INFILTRAZIONE.PERINEURALE")


#GRADING3
IC.sup = exp(fit$coefficients["GRADING3"] + qalpha * summary(fit)$se["GRADING3"])
IC.or = exp(fit$coefficients["GRADING3"])
IC.inf = exp(fit$coefficients["GRADING3"] - qalpha * summary(fit)$se["GRADING3"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GRADING3")

#Plot degli ODDS RATIO
df <- data.frame(x = names,
                 INF = unname(INF),
                 IC_Odds_Ratio = unname(CENT),
                 SUP = unname(SUP))

require(ggplot2)
ggplot(df, aes(x = x, y = IC_Odds_Ratio), colour =  x) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = SUP, ymin = INF, colour =  x)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(aes(yintercept = 1)) + 
  coord_cartesian(ylim=c(0, 3))

