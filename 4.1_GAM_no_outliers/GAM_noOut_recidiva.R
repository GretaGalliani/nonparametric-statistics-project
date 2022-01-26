#GOAL: Find the best GAM model to predict the binary outcome "RECIDIVA"

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
names(main)[names(main) == "Ca19.9â..55"] <- "Ca19.9.55" #for some windows pc


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

#creating dataframes for modelling

main_pc <- data.frame(main,pc)
#Leaving out other outcomes
main_pc$STATO.VIVO.MORTO <- NULL
main_pc$OS..Days. <- NULL
main_pc$RFS..Days. <- NULL
main_pc <- na.omit(main_pc)

main_pc_pm <- data.frame(main[-257,],pc[-257,],pm) #Last line of Portal Margin is empty
#Leaving out other outcomes
main_pc_pm$STATO.VIVO.MORTO <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL
main_pc_pm <- na.omit(main_pc_pm)

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
#Categorical variable significant: 
# ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.METASTATICI
# + NUMERO.LINFONODI.ASPORTATI + N

#Deleting from the dataset main the categorical variables not significant
which(coef.lasso == 0)

#removing variables discarded by lasso
main_pc$SEX <- NULL
main_pc$HBV <- NULL
main_pc$HCV <- NULL
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
main_pc$GRADING <- NULL
main_pc$T.VIII.ed <- NULL

main_pc_pm$SEX <- NULL
main_pc_pm$HBV <- NULL
main_pc_pm$HCV <- NULL
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
main_pc_pm$GRADING <- NULL
main_pc_pm$T.VIII.ed <- NULL

######## Model 1: Clinical Variables + Portal Core Radiomics -----

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c( "ASSOCIATED.RESECTION","SINGLE.NODULE","INVASIONE.VASCOLARE.MICROSCOPICA",
           "NUMERO.LINFONODI.METASTATICI","NUMERO.LINFONODI.ASPORTATI","N")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="RECIDIVA"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(RECIDIVA ~ AGE + Degenza,
                        data = main_pc,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
#Best model:
# AGE + NUMERO..NO.SATELLITI. + CONVENTIONAL_HUstd + s(CONVENTIONAL_HUQ2, df = 3) + GLZLM_SZLGE

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(RECIDIVA ~ ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + AGE + NUMERO..NO.SATELLITI. + 
                    CONVENTIONAL_HUstd + s(CONVENTIONAL_HUQ2, k = 3) + GLZLM_SZLGE,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# ASSOCIATED.RESECTION + NUMERO..NO.SATELLITI. + CONVENTIONAL_HUstd + s(CONVENTIONAL_HUQ2, k = 3)

# --->FINAL FIT:
fit2 <- mgcv::gam(RECIDIVA ~ SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + AGE + 
                    GLZLM_SZLGE,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK

#LOOCV
f <- RECIDIVA ~ SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
  NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + AGE + 
  GLZLM_SZLGE

res <- our_loocv(main_pc,f,"RECIDIVA")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #.7352941

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.875

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.4655172

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7327586





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
# SHAPE_Sphericity.onlyFor3DROI.. + s(GLRLM_SRHGE, df = 3) + NGLDM_Busyness + GLZLM_LZLGE + 
# DISCRETIZED_HISTO_Entropy_log10_margin + GLZLM_SZE_margin + GLZLM_LZLGE_margin

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(RECIDIVA ~ ASSOCIATED.RESECTION + SINGLE.NODULE + INVASIONE.VASCOLARE.MICROSCOPICA +
                    NUMERO.LINFONODI.METASTATICI + NUMERO.LINFONODI.ASPORTATI + SHAPE_Sphericity.onlyFor3DROI.. + 
                    s(GLRLM_SRHGE, k = 3) + NGLDM_Busyness + GLZLM_LZLGE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLZLM_SZE_margin + GLZLM_LZLGE_margin,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# ASSOCIATED.RESECTION + SINGLE.NODULE + SHAPE_Sphericity.onlyFor3DROI.. + NUMERO.LINFONODI.METASTATICI

# --->FINAL FIT:
fit2 <- mgcv::gam(RECIDIVA ~ INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.ASPORTATI + 
                    s(GLRLM_SRHGE, k = 3) + NGLDM_Busyness + GLZLM_LZLGE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLZLM_SZE_margin + GLZLM_LZLGE_margin,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables?
anova(fit1,fit2,test = "Chisq") #ok

#LOOCV
f <- RECIDIVA ~  INVASIONE.VASCOLARE.MICROSCOPICA + NUMERO.LINFONODI.ASPORTATI + 
  s(GLRLM_SRHGE, k = 3) + NGLDM_Busyness + GLZLM_LZLGE + DISCRETIZED_HISTO_Entropy_log10_margin + 
  GLZLM_SZE_margin + GLZLM_LZLGE_margin

res <- our_loocv(main_pc_pm,f,"RECIDIVA")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.739645

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.8660714

# Specificivity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.4912281

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7586153

