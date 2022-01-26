#GOAL: Find the best GAM model to predict the binary outcome "GRADING"

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

#grouping grading
main$GRADING <- ifelse(main$GRADING==3,1,0) #grad=0 if 1 or 2, grad=1 if 3

#### Feature selection on categorical variables: LASSO ----
main1 <- na.omit(main)
X <- model.matrix(~ SEX + HCV + HBV + Ca19.9.55 + CHEMIOTERAPIA.NEOADIUVANTE
                  + PVE.preop + PRIMA.RESEZIONE + ASSOCIATED.RESECTION
                  + CIRROSI + SINGLE.NODULE + INFILTRAZIONE.PERINEURALE 
                  + PATTERN + NUMERO..NO.SATELLITI.,
                  data=main1)[,-1]

lambda.grid <- 10^seq(5,-3,length=100) #defining lambda grid
cv.lasso <- cv.glmnet(X,main1$GRADING,nfolds= dim(main1)[1], lambda=lambda.grid)

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

coef.lasso <- predict(cv.lasso, s=bestlam.lasso, type = 'coefficients')[1:(dim(X)[2]+1),]
coef.lasso 
which(coef.lasso != 0) #NONE


######## Model 1: Clinical Variables + Portal Core Radiomics -----

main_pc <- data.frame(main,pc)
main_pc <- na.omit(main_pc)

#removing post surgery variables
main_pc$Degenza <- NULL
main_pc$Major.Hepatectomy <- NULL
main_pc$RESEZIONE.VIA.BILIARE <- NULL
main_pc$LINFOADENECTOMIA <- NULL
main_pc$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc$N <-NULL
main_pc$R.status <-NULL
main_pc$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc$INVASIONE.VASCOLARE.MICROSCOPICA <-NULL
main_pc$T.VIII.ed <-NULL          
main_pc$NODULI.SATELLITI <-NULL    
main_pc$COMPLICANZE.SEVERE <-NULL
main_pc$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc$STATO.VIVO.MORTO <- NULL
main_pc$RECIDIVA <- NULL
main_pc$OS..Days. <- NULL
main_pc$RFS..Days. <- NULL

#### Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM
#Dataset witch just continuous variables 
d_cont <- main_pc
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","PVE.preop",
           "PRIMA.RESEZIONE","ASSOCIATED.RESECTION","CIRROSI","SINGLE.NODULE",
           "INFILTRAZIONE.PERINEURALE","PATTERN","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="GRADING"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(GRADING ~ AGE,
                        data = main_pc,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)
#Best model:
# s(DIMENSIONE.MAX.MM, df = 3) + CONVENTIONAL_HUstd + CONVENTIONAL_HUQ2 + SHAPE_Volume.mL. + 
# SHAPE_Sphericity.onlyFor3DROI.. + GLZLM_LZLGE

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME


#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUstd + 
                    CONVENTIONAL_HUQ2 + SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + 
                    GLZLM_LZLGE,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove: nothing

# --->FINAL FIT:
fit2 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUstd + 
                    CONVENTIONAL_HUQ2 + SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + 
                    GLZLM_LZLGE,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#LOOCV
f <- GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUstd + 
  CONVENTIONAL_HUQ2 + SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + 
  GLZLM_LZLGE

res <- our_loocv(main_pc,f,"GRADING")

## METRICS
tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.6529412

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.3548387

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.8240741

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.6671147



######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----

main_pc_pm <- data.frame(main[-257,],pc[-257,],pm) #Last line of Portal Margin is empty
main_pc_pm <- na.omit(main_pc_pm)

#removing post surgery variables
main_pc_pm$Degenza <- NULL
main_pc_pm$Major.Hepatectomy <- NULL
main_pc_pm$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm$LINFOADENECTOMIA <- NULL
main_pc_pm$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc_pm$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc_pm$N <-NULL
main_pc_pm$R.status <-NULL
main_pc_pm$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc_pm$INVASIONE.VASCOLARE.MICROSCOPICA <-NULL
main_pc_pm$T.VIII.ed <-NULL          
main_pc_pm$NODULI.SATELLITI <-NULL    
main_pc_pm$COMPLICANZE.SEVERE <-NULL
main_pc_pm$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc_pm$STATO.VIVO.MORTO <- NULL
main_pc_pm$RECIDIVA <- NULL
main_pc_pm$OS..Days. <- NULL
main_pc_pm$RFS..Days. <- NULL


####Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM

#Dataset witch just continuous variables 
d_cont <- main_pc_pm
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","PVE.preop",
           "PRIMA.RESEZIONE","ASSOCIATED.RESECTION","CIRROSI","SINGLE.NODULE",
           "INFILTRAZIONE.PERINEURALE","PATTERN","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="GRADING"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)


start_model <- gam::gam(GRADING ~ AGE,
                        data = main_pc_pm,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm)
summary(step)

#Best model:
# s(DIMENSIONE.MAX.MM, df = 3) + CONVENTIONAL_HUmin + SHAPE_Volume.mL. + GLRLM_SRHGE + 
# DISCRETIZED_HISTO_Entropy_log10_margin + s(NGLDM_Coarseness_margin, df = 3)

## --->FINAL FIT:
fit1 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
                    SHAPE_Volume.mL. + GLRLM_SRHGE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    s(NGLDM_Coarseness_margin, k = 3),
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove: nothing

# --->FINAL FIT:
fit2 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
                    SHAPE_Volume.mL. + GLRLM_SRHGE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    s(NGLDM_Coarseness_margin, k = 3),
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)


#LOOCV
f <- GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
  SHAPE_Volume.mL. + GLRLM_SRHGE + DISCRETIZED_HISTO_Entropy_log10_margin + 
  s(NGLDM_Coarseness_margin, k = 3)

res <- our_loocv(main_pc_pm,f,"GRADING")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.739645

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.5483871

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.8504673

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7401266

