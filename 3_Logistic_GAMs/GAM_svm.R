#GOAL: Find the best GAM model to predict the binary outcome "STATO VIVO/MORTO"

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
pm <-scale(pm)

# Import PC's of radiomic's features
pc_pca <-read.table("Datasets/pc_post_pca.txt")
pm_pca <-read.table("Datasets/pm_post_pca.txt")
#Deleting labels
pc_pca <- pc_pca[,-c(1,2)]
pm_pca <- pm_pca[,-c(1,2)]


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
# + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + PATTERN + NUMERO.LINFONODI.METASTATICI
# + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed

#Deleting from the dataset main the categorical variables not significant
which(coef.lasso == 0)


######## Model 1: Clinical Variables + Portal Core Radiomics -----
main_pc <- data.frame(main,pc)
main_pc <- na.omit(main_pc)
main_pc$PVE.preop <- NULL
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
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","Major.Hepatectomy",
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
# AGE + Ca.19.9 + s(Degenza,df = 3) + s(CONVENTIONAL_HUmin, df = 2) + SHAPE_Surface.mm2..onlyFor3DROI. + 
#   s(GLCM_Contrast..Variance., df = 2) + GLCM_Correlation + 
#   s(GLRLM_SRHGE, df = 3) + GLRLM_GLNU


#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza,k = 3) + s(CONVENTIONAL_HUmin, k = 3) + SHAPE_Surface.mm2..onlyFor3DROI. + 
                    s(GLCM_Contrast..Variance., k = 3) + GLCM_Correlation + 
                    s(GLRLM_SRHGE, k = 3) + GLRLM_GLNU
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
#PATTERN + NUMERO.LINFONODI.METASTATICI + s(GLRLM_SRHGE, k = 3) + Ca.19.9 + RESEZIONE.VIA.BILIARE
#COMPLICANZE.SEVERE + Ca19.9.55 + INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + R.status
#HBV
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + s(Degenza,k = 3) + s(CONVENTIONAL_HUmin, k = 3) + SHAPE_Surface.mm2..onlyFor3DROI. + 
                    s(GLCM_Contrast..Variance., k = 3) + GLCM_Correlation + GLRLM_GLNU
                    + HCV + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + N + T.VIII.ed
                    + INVASIONE.VASCOLARE.MACROSCOPICA + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI.,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.576

final_model <- STATO.VIVO.MORTO ~ AGE +  + s(Degenza,k = 3) + s(CONVENTIONAL_HUmin, k = 3) + SHAPE_Surface.mm2..onlyFor3DROI.  + s(GLCM_Contrast..Variance., k = 3) + GLCM_Correlation + GLRLM_GLNU + HCV + CHEMIOTERAPIA.NEOADIUVANTE + Major.Hepatectomy + N + T.VIII.ed + INVASIONE.VASCOLARE.MACROSCOPICA + NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI.

#Cross-validation for the prediction
result <- our_loocv(main_pc, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.7209302

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.7294118

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.7126437

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
# 0.8091954






######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----

main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
main_pc_pm <- na.omit(main_pc_pm)
main_pc_pm$PVE.preop <- NULL
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
# AGE + s(Degenza, df = 2) + 
#   CONVENTIONAL_HUmin + SHAPE_Compacity.onlyFor3DROI. + GLCM_Contrast..Variance. + 
#   GLRLM_SRHGE + GLZLM_SZE + CONVENTIONAL_HUmax_margin + CONVENTIONAL_HUKurtosis_margin + 
#   SHAPE_Volume.mL._margin + s(SHAPE_Compacity.onlyFor3DROI._margin, df = 3) + 
#   GLRLM_LGRE_margin + GLZLM_LGZE_margin + s(GLZLM_HGZE_margin, df = 3) + s(GLZLM_LZLGE_margin, df = 2)

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + s(Degenza, k = 3) + 
                    CONVENTIONAL_HUmin + SHAPE_Compacity.onlyFor3DROI. + GLCM_Contrast..Variance. + 
                    GLRLM_SRHGE + GLZLM_SZE + CONVENTIONAL_HUmax_margin + CONVENTIONAL_HUKurtosis_margin + 
                    SHAPE_Volume.mL._margin + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) + 
                    GLRLM_LGRE_margin + GLZLM_LGZE_margin + s(GLZLM_HGZE_margin, k = 3) + s(GLZLM_LZLGE_margin, k = 3)+
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
#INVASIONE.VASCOLARE.MICROSCOPICA + PATTERN + RESEZIONE.VIA.BILIARE + GLRLM_LGRE_margin
#COMPLICANZE.SEVERE + NODULI.SATELLITI + CHEMIOTERAPIA.NEOADIUVANTE + NUMERO.LINFONODI.METASTATICI
#Ca19.9.55 + GLZLM_LGZE_margin + GLZLM_SZE + HBV
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + s(Degenza, k = 3) + 
                    CONVENTIONAL_HUmin + SHAPE_Compacity.onlyFor3DROI. + GLCM_Contrast..Variance. + 
                    GLRLM_SRHGE  + CONVENTIONAL_HUmax_margin + CONVENTIONAL_HUKurtosis_margin + 
                    SHAPE_Volume.mL._margin + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) + 
                    s(GLZLM_HGZE_margin, k = 3) + s(GLZLM_LZLGE_margin, k = 3) +
                    HCV + Major.Hepatectomy + R.status + INVASIONE.VASCOLARE.MACROSCOPICA + 
                    NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") 

final_model <-STATO.VIVO.MORTO ~ AGE + s(Degenza, k = 3) + 
  CONVENTIONAL_HUmin + SHAPE_Compacity.onlyFor3DROI. + GLCM_Contrast..Variance. + 
  GLRLM_SRHGE  + CONVENTIONAL_HUmax_margin + CONVENTIONAL_HUKurtosis_margin + 
  SHAPE_Volume.mL._margin + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) + 
  s(GLZLM_HGZE_margin, k = 3) + s(GLZLM_LZLGE_margin, k = 3) +
  HCV + Major.Hepatectomy + R.status + INVASIONE.VASCOLARE.MACROSCOPICA + 
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
accuracy # 0.7426901

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.7294118

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.7093023

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.8391245

### ODDS RATIO
alpha = 0.05
qalpha = qnorm( 1 - alpha/2 )

IC.sup = exp(m4$coefficients["CONVENTIONAL_HUmin"] + qalpha * summary(m4)$coefficients["CONVENTIONAL_HUmin",2])
IC.or = exp(m4$coefficients["CONVENTIONAL_HUmin"])
IC.inf = exp(m4$coefficients["CONVENTIONAL_HUmin"] - qalpha * summary(m4)$coefficients["CONVENTIONAL_HUmin",2])
c(inf = IC.inf, or = IC.or, sup = IC.sup )


######## Model 3: Clinical Variables + PC's (Portal Core) -----

main_pc_pca <- data.frame(main,pc_pca) #Last line of Portal Margin is empty
main_pc_pca <- na.omit(main_pc_pca)
main_pc_pca$PVE.preop <- NULL
main_pc_pca$PRIMA.RESEZIONE <- NULL
main_pc_pca$LINFOADENECTOMIA <- NULL
main_pc_pca$ASSOCIATED.RESECTION <- NULL
main_pc_pca$CIRROSI <- NULL
main_pc_pca$SINGLE.NODULE <- NULL
main_pc_pca$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pca$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pca$GRADING <- NULL

#Leaving out other outcomes
main_pc_pca$RECIDIVA <- NULL
main_pc_pca$OS..Days. <- NULL
main_pc_pca$RFS..Days. <- NULL

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
# AGE + Ca.19.9 + s(Degenza, df = 3) + Comp.1 + s(Comp.2, df = 2)

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza, k = 3) + Comp.1 + s(Comp.2, k = 3)
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
#PATTERN + INVASIONE.VASCOLARE.MICROSCOPICA + R.status + NUMERO.LINFONODI.METASTATICI + Ca19.9.55
#HBV + NODULI.SATELLITI + Ca.19.9 + COMPLICANZE.SEVERE + RESEZIONE.VIA.BILIARE + HCV
#CHEMIOTERAPIA.NEOADIUVANTE + NUMERO..NO.SATELLITI.
# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE +  + s(Degenza, k = 3) + Comp.1 + s(Comp.2, k = 3)
                  + Major.Hepatectomy + INVASIONE.VASCOLARE.MACROSCOPICA + NUMERO.LINFONODI.ASPORTATI +  
                  + N + T.VIII.ed,
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.13

final_model <-STATO.VIVO.MORTO ~ AGE +  + s(Degenza, k = 3) + Comp.1 + s(Comp.2, k = 3) + 
  Major.Hepatectomy + INVASIONE.VASCOLARE.MACROSCOPICA + NUMERO.LINFONODI.ASPORTATI + N + T.VIII.ed

#Cross-validation for the prediction
result <- our_loocv(main_pc_pca, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.6802326

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.6588235

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.7011494

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
# 0.7660581




######## Model 4: Clinical Variables + PC's (Portal Core & Margin) -----

main_pc_pm_pca <- data.frame(main[-261,],pc_pca[-261,],pm_pca)
main_pc_pm_pca <- na.omit(main_pc_pm_pca)
main_pc_pm_pca$PVE.preop <- NULL
main_pc_pm_pca$PRIMA.RESEZIONE <- NULL
main_pc_pm_pca$LINFOADENECTOMIA <- NULL
main_pc_pm_pca$ASSOCIATED.RESECTION <- NULL
main_pc_pm_pca$CIRROSI <- NULL
main_pc_pm_pca$SINGLE.NODULE <- NULL
main_pc_pm_pca$CHEMIOTERAPIA.ADIUVANTE <- NULL
main_pc_pm_pca$INFILTRAZIONE.PERINEURALE <- NULL
main_pc_pm_pca$GRADING <- NULL

#Leaving out other outcomes
main_pc_pm_pca$RECIDIVA <- NULL
main_pc_pm_pca$OS..Days. <- NULL
main_pc_pm_pca$RFS..Days. <- NULL

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
# AGE + Ca.19.9 + s(Degenza, df = 3) + Comp.1 + s(Comp.2, df = 2) + Comp.3_margin + Comp.4_margin

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME

#Fitting the model with categorical variable (lasso) + continuos variable (step)

fit1 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + Ca.19.9 + s(Degenza,k = 3) + Comp.1 + s(Comp.2, k = 3) + Comp.3_margin + Comp.4_margin
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
#Ca19.9.55 +PATTERN + R.status + HBV + NODULI.SATELLITI + NUMERO.LINFONODI.METASTATICI
#RESEZIONE.VIA.BILIARE + INVASIONE.VASCOLARE.MICROSCOPICA + Comp.4_margin + Ca.19.9
#COMPLICANZE.SEVERE + CHEMIOTERAPIA.NEOADIUVANTE + HCV + INVASIONE.VASCOLARE.MACROSCOPICA
#NUMERO..NO.SATELLITI.

# --->FINAL FIT:
fit2 <- mgcv::gam(STATO.VIVO.MORTO ~ AGE +  + s(Degenza,k = 3) + Comp.1 + s(Comp.2, k = 3) + Comp.3_margin 
                  + Major.Hepatectomy + NUMERO.LINFONODI.ASPORTATI + N + T.VIII.ed,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.103

final_model <- STATO.VIVO.MORTO ~ AGE + s(Degenza,k = 3) + Comp.1 + s(Comp.2, k = 3) + Comp.3_margin + Major.Hepatectomy + NUMERO.LINFONODI.ASPORTATI + N + T.VIII.ed
  
  
#Cross-validation for the prediction
result <- our_loocv(main_pc_pm_pca, final_model, "STATO.VIVO.MORTO")
pred1 <- result$predicted
obs1 <- result$observed

# Metrics
tab <- table(predicted = round(pred1), true = obs1)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy # 0.7134503

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity # 0.6727273

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity # 0.744186

# AUC & ROC CURVE
pred <- prediction(pred1, obs1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7718194


















######## ODDS RATIO FOR THE BEST MODEL (Model 2) ----

main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
main_pc_pm <- na.omit(main_pc_pm)
main_pc_pm$PVE.preop <- NULL
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

fit <- mgcv::gam(STATO.VIVO.MORTO ~ AGE + s(Degenza, k = 3) + 
                    CONVENTIONAL_HUmin + SHAPE_Compacity.onlyFor3DROI. + GLCM_Contrast..Variance. + 
                    GLRLM_SRHGE  + CONVENTIONAL_HUmax_margin + CONVENTIONAL_HUKurtosis_margin + 
                    SHAPE_Volume.mL._margin + s(SHAPE_Compacity.onlyFor3DROI._margin, k = 3) + 
                    s(GLZLM_HGZE_margin, k = 3) + s(GLZLM_LZLGE_margin, k = 3) +
                    HCV + Major.Hepatectomy + R.status + INVASIONE.VASCOLARE.MACROSCOPICA + 
                    NUMERO.LINFONODI.ASPORTATI + NUMERO..NO.SATELLITI. + N + T.VIII.ed,
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

#AGE 
IC.sup = exp(10*fit$coefficients["AGE"] + 10*qalpha * summary(fit)$se["AGE"])
IC.or = exp(10*fit$coefficients["AGE"])
IC.inf = exp(10*fit$coefficients["AGE"] - 10*qalpha * summary(fit)$se["AGE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"AGE (10)")

#CONVENTIONAL_HUmin 
IC.sup = exp(fit$coefficients["CONVENTIONAL_HUmin"] + qalpha * summary(fit)$se["CONVENTIONAL_HUmin"])
IC.or = exp(fit$coefficients["CONVENTIONAL_HUmin"])
IC.inf = exp(fit$coefficients["CONVENTIONAL_HUmin"] - qalpha * summary(fit)$se["CONVENTIONAL_HUmin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"CONVENTIONAL_HUmin")

#SHAPE_Compacity.onlyFor3DROI. 
IC.sup = exp(fit$coefficients["SHAPE_Compacity.onlyFor3DROI."] + qalpha * summary(fit)$se["SHAPE_Compacity.onlyFor3DROI."])
IC.or = exp(fit$coefficients["SHAPE_Compacity.onlyFor3DROI."])
IC.inf = exp(fit$coefficients["SHAPE_Compacity.onlyFor3DROI."] - qalpha * summary(fit)$se["SHAPE_Compacity.onlyFor3DROI."])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"SHAPE_Compacity.onlyFor3DROI.")

#GLCM_Contrast..Variance. 
IC.sup = exp(fit$coefficients["GLCM_Contrast..Variance."] + qalpha * summary(fit)$se["GLCM_Contrast..Variance."])
IC.or = exp(fit$coefficients["GLCM_Contrast..Variance."])
IC.inf = exp(fit$coefficients["GLCM_Contrast..Variance."] - qalpha * summary(fit)$se["GLCM_Contrast..Variance."])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLCM_Contrast..Variance.")

#GLRLM_SRHGE 
IC.sup = exp(fit$coefficients["GLRLM_SRHGE"] + qalpha * summary(fit)$se["GLRLM_SRHGE"])
IC.or = exp(fit$coefficients["GLRLM_SRHGE"])
IC.inf = exp(fit$coefficients["GLRLM_SRHGE"] - qalpha * summary(fit)$se["GLRLM_SRHGE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLRLM_SRHGE")

#CONVENTIONAL_HUmax_margin 
IC.sup = exp(fit$coefficients["CONVENTIONAL_HUmax_margin"] + qalpha * summary(fit)$se["CONVENTIONAL_HUmax_margin"])
IC.or = exp(fit$coefficients["CONVENTIONAL_HUmax_margin"])
IC.inf = exp(fit$coefficients["CONVENTIONAL_HUmax_margin"] - qalpha * summary(fit)$se["CONVENTIONAL_HUmax_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"CONVENTIONAL_HUmax_margin")

#CONVENTIONAL_HUKurtosis_margin 
IC.sup = exp(fit$coefficients["CONVENTIONAL_HUKurtosis_margin"] + qalpha * summary(fit)$se["CONVENTIONAL_HUKurtosis_margin"])
IC.or = exp(fit$coefficients["CONVENTIONAL_HUKurtosis_margin"])
IC.inf = exp(fit$coefficients["CONVENTIONAL_HUKurtosis_margin"] - qalpha * summary(fit)$se["CONVENTIONAL_HUKurtosis_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"CONVENTIONAL_HUKurtosis_margin")

#SHAPE_Volume.mL._margin 
IC.sup = exp(fit$coefficients["SHAPE_Volume.mL._margin"] + qalpha * summary(fit)$se["SHAPE_Volume.mL._margin"])
IC.or = exp(fit$coefficients["SHAPE_Volume.mL._margin"])
IC.inf = exp(fit$coefficients["SHAPE_Volume.mL._margin"] - qalpha * summary(fit)$se["SHAPE_Volume.mL._margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"SHAPE_Volume.mL._margin")

#HCV
IC.sup = exp(fit$coefficients["HCV"] + qalpha * summary(fit)$se["HCV"])
IC.or = exp(fit$coefficients["HCV"])
IC.inf = exp(fit$coefficients["HCV"] - qalpha * summary(fit)$se["HCV"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"HCV")

#Major.Hepatectomy 
IC.sup = exp(fit$coefficients["Major.Hepatectomy"] + qalpha * summary(fit)$se["Major.Hepatectomy"])
IC.or = exp(fit$coefficients["Major.Hepatectomy"])
IC.inf = exp(fit$coefficients["Major.Hepatectomy"] - qalpha * summary(fit)$se["Major.Hepatectomy"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"Major.Hepatectomy")

#R.status 
IC.sup = exp(fit$coefficients["R.status"] + qalpha * summary(fit)$se["R.status"])
IC.or = exp(fit$coefficients["R.status"])
IC.inf = exp(fit$coefficients["R.status"] - qalpha * summary(fit)$se["R.status"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"R.status")

#INVASIONE.VASCOLARE.MACROSCOPICA 
IC.sup = exp(fit$coefficients["INVASIONE.VASCOLARE.MACROSCOPICA"] + qalpha * summary(fit)$se["INVASIONE.VASCOLARE.MACROSCOPICA"])
IC.or = exp(fit$coefficients["INVASIONE.VASCOLARE.MACROSCOPICA"])
IC.inf = exp(fit$coefficients["INVASIONE.VASCOLARE.MACROSCOPICA"] - qalpha * summary(fit)$se["INVASIONE.VASCOLARE.MACROSCOPICA"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"INVASIONE.VASCOLARE.MACROSCOPICA")

#NUMERO.LINFONODI.ASPORTATI 
IC.sup = exp(fit$coefficients["NUMERO.LINFONODI.ASPORTATI"] + qalpha * summary(fit)$se["NUMERO.LINFONODI.ASPORTATI"])
IC.or = exp(fit$coefficients["NUMERO.LINFONODI.ASPORTATI"])
IC.inf = exp(fit$coefficients["NUMERO.LINFONODI.ASPORTATI"] - qalpha * summary(fit)$se["NUMERO.LINFONODI.ASPORTATI"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"NUMERO.LINFONODI.ASPORTATI")

#NUMERO..NO.SATELLITI. 
IC.sup = exp(fit$coefficients["NUMERO..NO.SATELLITI."] + qalpha * summary(fit)$se["NUMERO..NO.SATELLITI."])
IC.or = exp(fit$coefficients["NUMERO..NO.SATELLITI."])
IC.inf = exp(fit$coefficients["NUMERO..NO.SATELLITI."] - qalpha * summary(fit)$se["NUMERO..NO.SATELLITI."])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"NUMERO..NO.SATELLITI.")

#N1
IC.sup = exp(fit$coefficients["N1"] + qalpha * summary(fit)$se["N1"])
IC.or = exp(fit$coefficients["N1"])
IC.inf = exp(fit$coefficients["N1"] - qalpha * summary(fit)$se["N1"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"N1")

#Nx
IC.sup = exp(fit$coefficients["Nx"] + qalpha * summary(fit)$se["Nx"])
IC.or = exp(fit$coefficients["Nx"])
IC.inf = exp(fit$coefficients["Nx"] - qalpha * summary(fit)$se["Nx"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"Nx")

#T.VIII.ed1b
IC.sup = exp(fit$coefficients["T.VIII.ed1b"] + qalpha * summary(fit)$se["T.VIII.ed1b"])
IC.or = exp(fit$coefficients["T.VIII.ed1b"])
IC.inf = exp(fit$coefficients["T.VIII.ed1b"] - qalpha * summary(fit)$se["T.VIII.ed1b"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"T.VIII.ed1b")

#T.VIII.ed2
IC.sup = exp(fit$coefficients["T.VIII.ed2"] + qalpha * summary(fit)$se["T.VIII.ed2"])
IC.or = exp(fit$coefficients["T.VIII.ed2"])
IC.inf = exp(fit$coefficients["T.VIII.ed2"] - qalpha * summary(fit)$se["T.VIII.ed2"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"T.VIII.ed2")

#T.VIII.ed3
IC.sup = exp(fit$coefficients["T.VIII.ed3"] + qalpha * summary(fit)$se["T.VIII.ed3"])
IC.or = exp(fit$coefficients["T.VIII.ed3"])
IC.inf = exp(fit$coefficients["T.VIII.ed3"] - qalpha * summary(fit)$se["T.VIII.ed3"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"T.VIII.ed3")

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


