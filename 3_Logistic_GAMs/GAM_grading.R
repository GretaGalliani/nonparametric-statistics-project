#GOAL: Find the best GAM model to predict the binary outcome "GRADING"

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

source("3_Logistic_GAMs/models_functions.R") #function for LOOCV 

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
#random start model
start_model <- gam::gam(GRADING ~ AGE,
                        data = main_pc,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc)
summary(step)

#Best model:
# s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUQ2 + DISCRETIZED_HISTO_Entropy_log10 +
# SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + GLRLM_LGRE + GLRLM_GLNU +
# GLZLM_LZLGE + GLZLM_GLNU

#AFTER DOING AN AUTOMATIC FEATURE SELECTION WE MANUALLY REMOVE SOME OTHER COVARIATES WHICH SEEM
# TO NOT BE IMPORTANT ON PREDICTING THE OUTCOME


#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUQ2 + DISCRETIZED_HISTO_Entropy_log10 +
                    SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + GLRLM_LGRE + GLRLM_GLNU +
                    GLZLM_LZLGE + GLZLM_GLNU,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# GLRLM_LGRE

# --->FINAL FIT:
fit2 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUQ2 + DISCRETIZED_HISTO_Entropy_log10 +
                    SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI..  + GLRLM_GLNU +
                    GLZLM_LZLGE + GLZLM_GLNU,
                  data = main_pc,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#Is the model losing information after we removed this variables? No
anova(fit1,fit2,test = "Chisq") #OK: 0.09743

#LOOCV
f <- GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUQ2 + 
  DISCRETIZED_HISTO_Entropy_log10 + SHAPE_Volume.mL. + SHAPE_Sphericity.onlyFor3DROI.. + 
  GLRLM_GLNU + GLZLM_LZLGE + GLZLM_GLNU

res <- our_loocv(main_pc,f2,"GRADING")

## METRICS
tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.6511628

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.3492063

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.8256881

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.6538518



######## Model 2: Clinical Variables + Portal Core & Margin Radiomics -----

main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
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
# s(DIMENSIONE.MAX.MM, df = 3) + CONVENTIONAL_HUmin + CONVENTIONAL_HUmax + SHAPE_Surface.mm2..onlyFor3DROI. + 
# GLCM_Contrast..Variance. + GLRLM_SRHGE + GLZLM_LGZE + DISCRETIZED_HISTO_Entropy_log10_margin + 
# GLRLM_LRE_margin + GLRLM_LRLGE_margin + s(NGLDM_Coarseness_margin,df = 3) + GLZLM_ZLNU_margin + 
# GLZLM_ZP_margin

## --->FINAL FIT:
fit1 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
                    CONVENTIONAL_HUmax + SHAPE_Surface.mm2..onlyFor3DROI. + GLCM_Contrast..Variance. + 
                    GLRLM_SRHGE + GLZLM_LGZE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLRLM_LRE_margin + GLRLM_LRLGE_margin + s(NGLDM_Coarseness_margin, k = 3) + 
                    GLZLM_ZLNU_margin + GLZLM_ZP_margin,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# GLRLM_LRLGE_margin

# --->FINAL FIT:
fit2 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
                    CONVENTIONAL_HUmax + SHAPE_Surface.mm2..onlyFor3DROI. + GLCM_Contrast..Variance. + 
                    GLRLM_SRHGE + GLZLM_LGZE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                    GLRLM_LRE_margin + s(NGLDM_Coarseness_margin, k = 3) + 
                    GLZLM_ZLNU_margin + GLZLM_ZP_margin,
                  data = main_pc_pm,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)


#LOOCV
f <- GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
  CONVENTIONAL_HUmax + SHAPE_Surface.mm2..onlyFor3DROI. + GLCM_Contrast..Variance. + 
  GLRLM_SRHGE + GLZLM_LGZE + DISCRETIZED_HISTO_Entropy_log10_margin + 
  GLRLM_LRE_margin + s(NGLDM_Coarseness_margin, k = 3) + 
  GLZLM_ZLNU_margin + GLZLM_ZP_margin

res <- our_loocv(main_pc_pm,f,"GRADING")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.7368421

# Sensitivity
sensiticity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.5555556

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.8425926

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7771899


######## Model 3: Clinical Variables + Portal Core (PCA) -----

main_pc_pca <- data.frame(main,pc_pca)
main_pc_pca <-na.omit(main_pc_pca)

#removing post surgery variables
main_pc_pca$Degenza <- NULL
main_pc_pca$Major.Hepatectomy <- NULL
main_pc_pca$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pca$LINFOADENECTOMIA <- NULL
main_pc_pca$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc_pca$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc_pca$N <-NULL
main_pc_pca$R.status <-NULL
main_pc_pca$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc_pca$INVASIONE.VASCOLARE.MICROSCOPICA <-NULL
main_pc_pca$T.VIII.ed <-NULL          
main_pc_pca$NODULI.SATELLITI <-NULL    
main_pc_pca$COMPLICANZE.SEVERE <-NULL
main_pc_pca$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc_pca$STATO.VIVO.MORTO <- NULL
main_pc_pca$RECIDIVA <- NULL
main_pc_pca$OS..Days. <- NULL
main_pc_pca$RFS..Days. <- NULL

####Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM

#Dataset witch just continuous variables 
d_cont <- main_pc_pca
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","PVE.preop",
           "PRIMA.RESEZIONE","ASSOCIATED.RESECTION","CIRROSI","SINGLE.NODULE",
           "INFILTRAZIONE.PERINEURALE","PATTERN","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="GRADING"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(GRADING ~ AGE,
                        data = main_pc_pca,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pca)
summary(step)
#Best model:
# s(DIMENSIONE.MAX.MM, df = 2) + s(Comp.4, df = 2)

#Fitting the model with categorical variable (lasso) + continuos variable (step)
fit1 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + s(Comp.4, k = 3),
                  data = main_pc_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)


#LOOCV
f <- GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + s(Comp.4, k = 3)

res <- our_loocv(main_pc_pca,f,"GRADING")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.6162791

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.2222222

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.8440367

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.6082714


######## Model 4: Clinical Variables + Portal Core & Margin (PCA) -----

main_pc_pm_pca <- data.frame(main[-261,],pc_pca[-261,],pm_pca)
main_pc_pm_pca <-na.omit(main_pc_pm_pca)

#removing post surgery variables
main_pc_pm_pca$Degenza <- NULL
main_pc_pm_pca$Major.Hepatectomy <- NULL
main_pc_pm_pca$RESEZIONE.VIA.BILIARE <- NULL
main_pc_pm_pca$LINFOADENECTOMIA <- NULL
main_pc_pm_pca$NUMERO.LINFONODI.ASPORTATI <-NULL
main_pc_pm_pca$NUMERO.LINFONODI.METASTATICI <-NULL
main_pc_pm_pca$N <-NULL
main_pc_pm_pca$R.status <-NULL
main_pc_pm_pca$INVASIONE.VASCOLARE.MACROSCOPICA <-NULL     
main_pc_pm_pca$INVASIONE.VASCOLARE.MICROSCOPICA <-NULL
main_pc_pm_pca$T.VIII.ed <-NULL          
main_pc_pm_pca$NODULI.SATELLITI <-NULL    
main_pc_pm_pca$COMPLICANZE.SEVERE <-NULL
main_pc_pm_pca$CHEMIOTERAPIA.ADIUVANTE <-NULL
main_pc_pm_pca$STATO.VIVO.MORTO <- NULL
main_pc_pm_pca$RECIDIVA <- NULL
main_pc_pm_pca$OS..Days. <- NULL
main_pc_pm_pca$RFS..Days. <- NULL

####Feature selection on continuous variables (on which I can fit splines): STEPWISE GAM

#Dataset witch just continuous variables 
d_cont <- main_pc_pm_pca
d_cont[,c( "SEX","HCV","HBV","Ca19.9.55","CHEMIOTERAPIA.NEOADIUVANTE","PVE.preop",
           "PRIMA.RESEZIONE","ASSOCIATED.RESECTION","CIRROSI","SINGLE.NODULE",
           "INFILTRAZIONE.PERINEURALE","PATTERN","NUMERO..NO.SATELLITI.")] <- NULL

scope_list <- gam.scope(d_cont, response=which(names(d_cont)=="GRADING"), 
                        smoother = "s", arg = c("df=2","df=3"), form = TRUE)

start_model <- gam::gam(GRADING ~ AGE,
                        data = main_pc_pm_pca,
                        family = binomial)
summary(start_model)


step <- gam::step.Gam(start_model, scope_list, direction = "both",trace=TRUE, data=main_pc_pm_pca)
summary(step)

#Best model:
# s(DIMENSIONE.MAX.MM, k = 3) + Comp.1 + Comp.3 + s(Comp.1_margin, k = 3) + Comp.2_margin

## --->FINAL FIT:
fit1 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + Comp.3 + 
                    s(Comp.1_margin, k = 3) + Comp.2_margin + Comp.4_margin,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit1)

# I one at a time remove:
# Comp.4_margin

# --->FINAL FIT:
fit2 <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + Comp.3 + 
                    s(Comp.1_margin, k = 3) + Comp.2_margin,
                  data = main_pc_pm_pca,
                  method='REML',
                  family = binomial,
                  select = T)
summary(fit2)

#LOOCV
f <- GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + Comp.3 + s(Comp.1_margin, k = 3) + 
  Comp.2_margin

res <- our_loocv(main_pc_pm_pca,f,"GRADING")

## METRICS

tab <- table(predicted = round(res$predicted), true = res$observed)
tab

# Accuracy
accuracy <- sum(diag(tab))/sum(tab)
accuracy #0.6666667

# Sensitivity 
sensitivity <- tab[2,2]/(tab[2,2] + tab[1,2])
sensitivity #0.3333333

# Specificity 
specificity <- tab[1,1]/(tab[1,1]+tab[2,1])
specificity #0.8611111

# ROC
pred <- prediction(res$predicted, res$observed)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.6512346






######## ODDS RATIO FOR THE BEST MODEL (Model 2) ----
main_pc_pm <- data.frame(main[-261,],pc[-261,],pm) #Last line of Portal Margin is empty
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

fit <- mgcv::gam(GRADING ~ s(DIMENSIONE.MAX.MM, k = 3) + CONVENTIONAL_HUmin + 
                   CONVENTIONAL_HUmax + SHAPE_Surface.mm2..onlyFor3DROI. + GLCM_Contrast..Variance. + 
                   GLRLM_SRHGE + GLZLM_LGZE + DISCRETIZED_HISTO_Entropy_log10_margin + 
                   GLRLM_LRE_margin + s(NGLDM_Coarseness_margin, k = 3) + 
                   GLZLM_ZLNU_margin + GLZLM_ZP_margin,
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

#CONVENTIONAL_HUmax 
IC.sup = exp(fit$coefficients["CONVENTIONAL_HUmax"] + qalpha * summary(fit)$se["CONVENTIONAL_HUmax"])
IC.or = exp(fit$coefficients["CONVENTIONAL_HUmax"])
IC.inf = exp(fit$coefficients["CONVENTIONAL_HUmax"] - qalpha * summary(fit)$se["CONVENTIONAL_HUmax"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"CONVENTIONAL_HUmax")

#SHAPE_Surface.mm2..onlyFor3DROI. 
IC.sup = exp(fit$coefficients["SHAPE_Surface.mm2..onlyFor3DROI."] + qalpha * summary(fit)$se["SHAPE_Surface.mm2..onlyFor3DROI."])
IC.or = exp(fit$coefficients["SHAPE_Surface.mm2..onlyFor3DROI."])
IC.inf = exp(fit$coefficients["SHAPE_Surface.mm2..onlyFor3DROI."] - qalpha * summary(fit)$se["SHAPE_Surface.mm2..onlyFor3DROI."])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"SHAPE_Surface.mm2..onlyFor3DROI.")

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

#GLZLM_LGZE 
IC.sup = exp(fit$coefficients["GLZLM_LGZE"] + qalpha * summary(fit)$se["GLZLM_LGZE"])
IC.or = exp(fit$coefficients["GLZLM_LGZE"])
IC.inf = exp(fit$coefficients["GLZLM_LGZE"] - qalpha * summary(fit)$se["GLZLM_LGZE"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_LGZE")

#DISCRETIZED_HISTO_Entropy_log10_margin 
IC.sup = exp(fit$coefficients["DISCRETIZED_HISTO_Entropy_log10_margin"] + qalpha * summary(fit)$se["DISCRETIZED_HISTO_Entropy_log10_margin"])
IC.or = exp(fit$coefficients["DISCRETIZED_HISTO_Entropy_log10_margin"])
IC.inf = exp(fit$coefficients["DISCRETIZED_HISTO_Entropy_log10_margin"] - qalpha * summary(fit)$se["DISCRETIZED_HISTO_Entropy_log10_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"DISCRETIZED_HISTO_Entropy_log10_margin")

#GLRLM_LRE_margin 
IC.sup = exp(fit$coefficients["GLRLM_LRE_margin"] + qalpha * summary(fit)$se["GLRLM_LRE_margin"])
IC.or = exp(fit$coefficients["GLRLM_LRE_margin"])
IC.inf = exp(fit$coefficients["GLRLM_LRE_margin"] - qalpha * summary(fit)$se["GLRLM_LRE_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLRLM_LRE_margin")

#GLZLM_ZLNU_margin 
IC.sup = exp(fit$coefficients["GLZLM_ZLNU_margin"] + qalpha * summary(fit)$se["GLZLM_ZLNU_margin"])
IC.or = exp(fit$coefficients["GLZLM_ZLNU_margin"])
IC.inf = exp(fit$coefficients["GLZLM_ZLNU_margin"] - qalpha * summary(fit)$se["GLZLM_ZLNU_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_ZLNU_margin")

#GLZLM_ZP_margin 
IC.sup = exp(fit$coefficients["GLZLM_ZP_margin"] + qalpha * summary(fit)$se["GLZLM_ZP_margin"])
IC.or = exp(fit$coefficients["GLZLM_ZP_margin"])
IC.inf = exp(fit$coefficients["GLZLM_ZP_margin"] - qalpha * summary(fit)$se["GLZLM_ZP_margin"])
c(inf = IC.inf, or = IC.or, sup = IC.sup )

INF <- c(INF, IC.inf)
CENT <- c(CENT, IC.or)
SUP <- c(SUP, IC.sup)
names <- c(names,"GLZLM_ZP_margin")

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

