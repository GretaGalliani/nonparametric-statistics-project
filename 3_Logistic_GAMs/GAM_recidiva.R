#GOAL: Find the best GAM model to predict the binary outcome "RECIDIVA"

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
names(main)[names(main) == "Ca19.9â..55"] <- "Ca19.9.55" #for some windows pc

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

#creating dataframes for modelling

main_pc <- data.frame(main,pc)
#Leaving out other outcomes
main_pc$STATO.VIVO.MORTO <- NULL
main_pc$OS..Days. <- NULL
main_pc$RFS..Days. <- NULL
main_pc <- na.omit(main_pc)

main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
#Leaving out other outcomes
main_pc_pm$STATO.VIVO.MORTO <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL
main_pc_pm <- na.omit(main_pc_pm)

main_pc_pca <- data.frame(main,pc_pca)
#Leaving out other outcomes
main_pc_pca$STATO.VIVO.MORTO <- NULL
main_pc_pca$OS..Days. <- NULL
main_pc_pca$RFS..Days. <- NULL
main_pc_pca <-na.omit(main_pc_pca)

main_pc_pm_pca <- data.frame(main[-261,],pc_pca[-261,],pm_pca)
#Leaving out other outcomes
main_pc_pm_pca$STATO.VIVO.MORTO  <- NULL
main_pc_pm_pca$OS..Days. <- NULL
main_pc_pm_pca$RFS..Days. <- NULL
main_pc_pm_pca <-na.omit(main_pc_pm_pca)

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
cv.lasso <- cv.glmnet(X,main1$RECIDIVA,nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0)  

# ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
# + NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI

#removing from datasets the variables discarded by lasso
#model 1
main_pc$SEX <- NULL
main_pc$HCV <- NULL
main_pc$HBV <- NULL
main_pc$Ca19.9.55 <- NULL
main_pc$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc$PVE.preop <- NULL
main_pc$Major.Hepatectomy <- NULL
main_pc$PRIMA.RESEZIONE <- NULL
main_pc$RESEZIONE.VIA.BILIARE <- NULL
main_pc$LINFOADENECTOMIA <- NULL
main_pc$COMPLICANZE.SEVERE <- NULL
main_pc$CIRROSI <- NULL
main_pc$R.status <- NULL
main_pc$INVASIONE.VASCOLARE.MACROSCOPICA <- NULL
main_pc$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc$INFILTRAZIONE.PERINEURALE <- NULL
main_pc$NODULI.SATELLITI <- NULL
main_pc$PATTERN <- NULL
main_pc$NUMERO..NO.SATELLITI. <- NULL
main_pc$GRADING <- NULL
main_pc$N <- NULL
main_pc$T.VIII.ed <- NULL

#model 2
main_pc_pm$SEX <- NULL
main_pc_pm$HCV <- NULL
main_pc_pm$HBV <- NULL
main_pc_pm$Ca19.9.55 <- NULL
main_pc_pm$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pm$PVE.preop <- NULL
main_pc_pm$Major.Hepatectomy <- NULL
main_pc_pm$PRIMA.RESEZIONE <- NULL
main_pc_pm$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm$LINFOADENECTOMIA <- NULL
main_pc_pm$COMPLICANZE.SEVERE <- NULL
main_pc_pm$CIRROSI <- NULL
main_pc_pm$R.status <- NULL
main_pc_pm$INVASIONE.VASCOLARE.MACROSCOPICA <- NULL
main_pc_pm$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pm$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pm$NODULI.SATELLITI <- NULL
main_pc_pm$PATTERN <- NULL
main_pc_pm$NUMERO..NO.SATELLITI. <- NULL
main_pc_pm$GRADING <- NULL
main_pc_pm$N <- NULL
main_pc_pm$T.VIII.ed <- NULL

#model 3
main_pc_pca$SEX <- NULL
main_pc_pca$HCV <- NULL
main_pc_pca$HBV <- NULL
main_pc_pca$Ca19.9.55 <- NULL
main_pc_pca$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pca$PVE.preop <- NULL
main_pc_pca$Major.Hepatectomy <- NULL
main_pc_pca$PRIMA.RESEZIONE <- NULL
main_pc_pca$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pca$LINFOADENECTOMIA <- NULL
main_pc_pca$COMPLICANZE.SEVERE <- NULL
main_pc_pca$CIRROSI <- NULL
main_pc_pca$R.status <- NULL
main_pc_pca$INVASIONE.VASCOLARE.MACROSCOPICA <- NULL
main_pc_pca$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pca$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pca$NODULI.SATELLITI <- NULL
main_pc_pca$PATTERN <- NULL
main_pc_pca$NUMERO..NO.SATELLITI. <- NULL
main_pc_pca$GRADING <- NULL
main_pc_pca$N <- NULL
main_pc_pca$T.VIII.ed <- NULL

#model 4
main_pc_pm_pca$SEX <- NULL
main_pc_pm_pca$HCV <- NULL
main_pc_pm_pca$HBV <- NULL
main_pc_pm_pca$Ca19.9.55 <- NULL
main_pc_pm_pca$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pm_pca$PVE.preop <- NULL
main_pc_pm_pca$Major.Hepatectomy <- NULL
main_pc_pm_pca$PRIMA.RESEZIONE <- NULL
main_pc_pm_pca$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm_pca$LINFOADENECTOMIA <- NULL
main_pc_pm_pca$COMPLICANZE.SEVERE <- NULL
main_pc_pm_pca$CIRROSI <- NULL
main_pc_pm_pca$R.status <- NULL
main_pc_pm_pca$INVASIONE.VASCOLARE.MACROSCOPICA <- NULL
main_pc_pm_pca$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pm_pca$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pm_pca$NODULI.SATELLITI <- NULL
main_pc_pm_pca$PATTERN <- NULL
main_pc_pm_pca$NUMERO..NO.SATELLITI. <- NULL
main_pc_pm_pca$GRADING <- NULL
main_pc_pm_pca$N <- NULL
main_pc_pm_pca$T.VIII.ed <- NULL

######## Model 1: Clinical Variables + Portal Core Radiomics -----

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c( "ASSOCIATED.RESECTION","SINGLE.NODULE","INVASIONE.VASCOLARE.MICROSCOPICA",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="RECIDIVA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(RECIDIVA ~ AGE + Degenza,
                        data = main_pc,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
#Best model:
# CONVENTIONAL_HUstd + CONVENTIONAL_HUKurtosis + SHAPE_Surface.mm2..onlyFor3DROI. + 
# GLRLM_LRE + GLRLM_GLNU + NGLDM_Contrast + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_GLNU + 
# GLZLM_ZP

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(RECIDIVA ~ ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + CONVENTIONAL_HUstd + 
                    CONVENTIONAL_HUKurtosis + SHAPE_Surface.mm2..onlyFor3DROI. + GLRLM_LRE + GLRLM_GLNU + 
                    NGLDM_Contrast + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_GLNU + GLZLM_ZP,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# ASSOCIATED.RESECTION + NUMERO.LINFONODI.METASTATICI + SHAPE_Surface.mm2..onlyFor3DROI. +
# CONVENTIONAL_HUKurtosis + GLRLM_LRE

# --->FINAL FIT:
fit2 <- mgcv::gam(RECIDIVA ~ SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.ASPORTATI + CONVENTIONAL_HUstd + 
                    GLRLM_GLNU + NGLDM_Contrast + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_GLNU + 
                    GLZLM_ZP,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.159

#LOOCV
f <- RECIDIVA ~ SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA + 
  NUMERO.LINFONODI.ASPORTATI + CONVENTIONAL_HUstd + GLRLM_GLNU + 
  NGLDM_Contrast + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_GLNU + 
  GLZLM_ZP

res <- our_loocv(main_pc,f,"RECIDIVA")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.7325581

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.8421053

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.5172414

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7309437





######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----

####Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM

#Dataset witch just continuous variables 
d_cont <- main_pc_pm
d_cont[,c( "ASSOCIATED.RESECTION","SINGLE.NODULE","INVASIONE.VASCOLARE.MICROSCOPICA",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="RECIDIVA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(RECIDIVA ~ AGE + Degenza,
                        data = main_pc_pm,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm)
summary(step)

#Best model:
# GLCM_Correlation + NGLDM_Contrast + GLZLM_LGZE + GLZLM_LZLGE + GLZLM_GLNU + GLZLM_ZP + 
# DISCRETIZED_HISTO_Entropy_log10_margin + GLZLM_SZE_margin + GLZLM_LZLGE_margin

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(RECIDIVA ~ ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + GLCM_Correlation + NGLDM_Contrast + 
                    GLZLM_LGZE + GLZLM_LZLGE + GLZLM_GLNU + GLZLM_ZP + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLZLM_SZE_margin + GLZLM_LZLGE_margin,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# ASSOCIATED.RESECTION + NUMERO.LINFONODI.METASTATICI + SINGLE.NODULE

# --->FINAL FIT:
fit2 <- mgcv::gam(RECIDIVA ~  INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.ASPORTATI + 
                    GLCM_Correlation + NGLDM_Contrast + GLZLM_LGZE + GLZLM_LZLGE + GLZLM_GLNU + 
                    GLZLM_ZP + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLZLM_SZE_margin + GLZLM_LZLGE_margin,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables?
anova(fit1,fit2,test = "Chisq") #ok: pvalue=0.1365

#LOOCV
f <- RECIDIVA ~  INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.ASPORTATI + 
  GLCM_Correlation + NGLDM_Contrast + GLZLM_LGZE + GLZLM_LZLGE + GLZLM_GLNU + 
  GLZLM_ZP + DISCRETIZED_HISTO_Entropy_log10_margin + 
  GLZLM_SZE_margin + GLZLM_LZLGE_margin

res <- our_loocv(main_pc_pm,f,"RECIDIVA")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.7836257

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.8947368

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.6666667

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7873192


######## Model 3: Clinical Variables + Portal Core (PCA) -----

####Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM

#Dataset witch just continuous variables 
d_cont <- main_pc_pca
d_cont[,c( "ASSOCIATED.RESECTION","SINGLE.NODULE","INVASIONE.VASCOLARE.MICROSCOPICA",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="RECIDIVA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(RECIDIVA ~ AGE + Degenza,
                        data = main_pc_pca,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pca)
summary(step)
#Best model:
# Comp.1 + Comp.3

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(RECIDIVA ~ ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + 
                    Comp.1 + Comp.3,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# ASSOCIATED.RESECTION + NUMERO.LINFONODI.METASTATICI

# --->FINAL FIT:
fit2 <- mgcv::gam(RECIDIVA ~  SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.ASPORTATI + 
                    Comp.1 + Comp.3,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables?
anova(fit1,fit2,test = "Chisq") #0.05816

#LOOCV
f <- RECIDIVA ~ SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
  NUMERO.LINFONODI.ASPORTATI + Comp.1 + Comp.3

res <- our_loocv(main_pc_pca,f,"RECIDIVA")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.7209302

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.8859649

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.3965517

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7200544

######## Model 4: Clinical Variables + Portal Core & Margin (PCA) -----

####Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM

#Dataset witch just continuous variables 
d_cont <- main_pc_pm_pca
d_cont[,c( "ASSOCIATED.RESECTION","SINGLE.NODULE","INVASIONE.VASCOLARE.MICROSCOPICA",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="RECIDIVA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(RECIDIVA ~ AGE + Degenza,
                        data = main_pc_pm_pca,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm_pca)
summary(step)
#Best model:
# s(Comp.1, df = 3) + Comp.3 + s(Comp.1_margin, df = 3)

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(RECIDIVA ~ ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI +
                    s(Comp.1, k = 3) + Comp.3 + s(Comp.1_margin, k = 3),
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# ASSOCIATED.RESECTION + NUMERO.LINFONODI.METASTATICI

# --->FINAL FIT:
fit2 <- mgcv::gam(RECIDIVA ~  SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.ASPORTATI + s(Comp.1, k = 3) + s(Comp.1_margin, k = 3),
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables?
anova(fit1,fit2,test = "Chisq") #0.09554

#LOOCV
f <- RECIDIVA ~  SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
  NUMERO.LINFONODI.ASPORTATI + s(Comp.1, k = 3) + s(Comp.1_margin, k = 3)

res <- our_loocv(main_pc_pm_pca,f,"RECIDIVA")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.7309942

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.9122807

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.3684211

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.6783626






######## ODDS RATIO FOR THE BEST MODEL (Model 2) ----
main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
#Leaving out other outcomes
main_pc_pm$STATO.VIVO.MORTO <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL
main_pc_pm <- na.omit(main_pc_pm)

main_pc_pm$SEX <- NULL
main_pc_pm$HCV <- NULL
main_pc_pm$HBV <- NULL
main_pc_pm$Ca19.9.55 <- NULL
main_pc_pm$CHEMIOTERAPIA.NEOADIUVANTE <- NULL
main_pc_pm$PVE.preop <- NULL
main_pc_pm$Major.Hepatectomy <- NULL
main_pc_pm$PRIMA.RESEZIONE <- NULL
main_pc_pm$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm$LINFOADENECTOMIA <- NULL
main_pc_pm$COMPLICANZE.SEVERE <- NULL
main_pc_pm$CIRROSI <- NULL
main_pc_pm$R.status <- NULL
main_pc_pm$INVASIONE.VASCOLARE.MACROSCOPICA <- NULL
main_pc_pm$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pm$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pm$NODULI.SATELLITI <- NULL
main_pc_pm$PATTERN <- NULL
main_pc_pm$NUMERO..NO.SATELLITI. <- NULL
main_pc_pm$GRADING <- NULL
main_pc_pm$N <- NULL
main_pc_pm$T.VIII.ed <- NULL

fit <- mgcv::gam(RECIDIVA ~  INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.ASPORTATI + 
                    GLCM_Correlation + NGLDM_Contrast + GLZLM_LGZE + GLZLM_LZLGE + GLZLM_GLNU + 
                    GLZLM_ZP + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLZLM_SZE_margin + GLZLM_LZLGE_margin,
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

#GLCM_Correlation 
IC.sup = exp(fit$coefficients["GLCM_Correlation"] + qalpha * summary(fit)$se["GLCM_Correlation"])
IC.or = exp(fit$coefficients["GLCM_Correlation"])
IC.inf = exp(fit$coefficients["GLCM_Correlation"] - qalpha * summary(fit)$se["GLCM_Correlation"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLCM_Correlation")

#NGLDM_Contrast 
IC.sup = exp(fit$coefficients["NGLDM_Contrast"] + qalpha * summary(fit)$se["NGLDM_Contrast"])
IC.or = exp(fit$coefficients["NGLDM_Contrast"])
IC.inf = exp(fit$coefficients["NGLDM_Contrast"] - qalpha * summary(fit)$se["NGLDM_Contrast"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"NGLDM_Contrast")

#GLZLM_LGZE 
IC.sup = exp(fit$coefficients["GLZLM_LGZE"] + qalpha * summary(fit)$se["GLZLM_LGZE"])
IC.or = exp(fit$coefficients["GLZLM_LGZE"])
IC.inf = exp(fit$coefficients["GLZLM_LGZE"] - qalpha * summary(fit)$se["GLZLM_LGZE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_LGZE (1e-5)")

#GLZLM_LZLGE 
IC.sup = exp(fit$coefficients["GLZLM_LZLGE"] + qalpha * summary(fit)$se["GLZLM_LZLGE"])
IC.or = exp(fit$coefficients["GLZLM_LZLGE"])
IC.inf = exp(fit$coefficients["GLZLM_LZLGE"] - qalpha * summary(fit)$se["GLZLM_LZLGE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_LZLGE")

#GLZLM_GLNU 
IC.sup = exp(fit$coefficients["GLZLM_GLNU"] + qalpha * summary(fit)$se["GLZLM_GLNU"])
IC.or = exp(fit$coefficients["GLZLM_GLNU"])
IC.inf = exp(fit$coefficients["GLZLM_GLNU"] - qalpha * summary(fit)$se["GLZLM_GLNU"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_GLNU")

#GLZLM_ZP 
IC.sup = exp(fit$coefficients["GLZLM_ZP"] + qalpha * summary(fit)$se["GLZLM_ZP"])
IC.or = exp(fit$coefficients["GLZLM_ZP"])
IC.inf = exp(fit$coefficients["GLZLM_ZP"] - qalpha * summary(fit)$se["GLZLM_ZP"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_ZP")

#DISCRETIZED_HISTO_Entropy_log10_margin 
IC.sup = exp(fit$coefficients["DISCRETIZED_HISTO_Entropy_log10_margin"] + qalpha * summary(fit)$se["DISCRETIZED_HISTO_Entropy_log10_margin"])
IC.or = exp(fit$coefficients["DISCRETIZED_HISTO_Entropy_log10_margin"])
IC.inf = exp(fit$coefficients["DISCRETIZED_HISTO_Entropy_log10_margin"] - qalpha * summary(fit)$se["DISCRETIZED_HISTO_Entropy_log10_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"DISCRETIZED_HISTO_Entropy_log10_margin")

#GLZLM_SZE_margin 
IC.sup = exp(fit$coefficients["GLZLM_SZE_margin"] + qalpha * summary(fit)$se["GLZLM_SZE_margin"])
IC.or = exp(fit$coefficients["GLZLM_SZE_margin"])
IC.inf = exp(fit$coefficients["GLZLM_SZE_margin"] - qalpha * summary(fit)$se["GLZLM_SZE_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_SZE_margin")

#GLZLM_LZLGE_margin  
IC.sup = exp(fit$coefficients["GLZLM_LZLGE_margin"] + qalpha * summary(fit)$se["GLZLM_LZLGE_margin"])
IC.or = exp(fit$coefficients["GLZLM_LZLGE_margin"])
IC.inf = exp(fit$coefficients["GLZLM_LZLGE_margin"] - qalpha * summary(fit)$se["GLZLM_LZLGE_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_LZLGE_margin ")

#INVASIONE.VASCOLARE.MICROSCOPICA 
IC.sup = exp(fit$coefficients["INVASIONE.VASCOLARE.MICROSCOPICA"] + qalpha * summary(fit)$se["INVASIONE.VASCOLARE.MICROSCOPICA"])
IC.or = exp(fit$coefficients["INVASIONE.VASCOLARE.MICROSCOPICA"])
IC.inf = exp(fit$coefficients["INVASIONE.VASCOLARE.MICROSCOPICA"] - qalpha * summary(fit)$se["INVASIONE.VASCOLARE.MICROSCOPICA"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"INVASIONE.VASCOLARE.MICROSCOPICA")


#NUMERO.LINFONODI.ASPORTATI
IC.sup = exp(fit$coefficients["NUMERO.LINFONODI.ASPORTATI"] + qalpha * summary(fit)$se["NUMERO.LINFONODI.ASPORTATI"])
IC.or = exp(fit$coefficients["NUMERO.LINFONODI.ASPORTATI"])
IC.inf = exp(fit$coefficients["NUMERO.LINFONODI.ASPORTATI"] - qalpha * summary(fit)$se["NUMERO.LINFONODI.ASPORTATI"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"NUMERO.LINFONODI.ASPORTATI")

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
