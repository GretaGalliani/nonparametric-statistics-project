#########__________ SURVIVAL MODELS (OS DAYS)__________#########

# Loading libraries
library(readxl)
library(car)
library(caret)
library(survival)
library(survminer) 
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)
library(RColorBrewer)

# Setting a seed
set.seed(24091998)

######### Dataset preparation ######### 
# now we import the dataset after robust correlation reduction 
main <- read.table("Datasets/main.txt")
pc_radiomics <-read.table("Datasets/pc_post_robcor.txt")
pm_radiomics <-read.table("Datasets/pm_post_robcor.txt")
# discarding labels which we are not considering in these tests
pc_radiomics <- pc_radiomics[,-c(1,2)]
pm_radiomics <- pm_radiomics[,-c(1,2)]

# removing outliers found by robust analysis
pc_radiomics <- pc_radiomics[-c(102,122,204,235,241),]
pm_radiomics <- pm_radiomics[-c(102,122,204,235,241),]
main <-main[-c(102,122,204,235,241),]

# Standardize the variables
pc_radiomics <- scale(pc_radiomics)
pm_radiomics <- scale(pm_radiomics)

# Constructing the full dataset
main_pc <- cbind(main, pc_radiomics)

main_pc_pm <- cbind(main[-256,], pc_radiomics[-256,], pm_radiomics)

# Removing NAs
main_pc <- na.omit(main_pc)
main_pc_pm <- na.omit(main_pc_pm)

##### KAPLAN-MEYER ESTIMATOR #####
# Creating the survival object 
surv_obj <- Surv(main_pc$OS..Days., main_pc$STATO.VIVO.MORTO)

# Computing the Kaplan-Meier estimator
fit <- survfit(surv_obj  ~ 1, data = main_pc)

# Reporting the table for the KM estimator
kable(head(tidy(fit),20))

# The medial survival time is 
surv_median(fit) # 1187 days

##### KAPLAN-MEYER PLOTS #####
plot(fit, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='steelblue2',
     main="Kaplan-Meier Curve for Cholangiocarcinoma - OS")

# Another possible plot
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette  = brewer.pal(n = 1, name = "Set2"),
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90, # title="Kaplan-Meier Curve for Cholangiocarcinoma - OS"
           )

# Computation of CFP (cumulative probabilities of experiencing the event of interest)
cumulative_incidence <- 1 - fit$surv

# Plot this quantity
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           fun='event',
           title="Cumulative Incidence Curve for Cholangiocarcinoma - OS")

# CUMULATIVE HAZARD (it measures the total amount of risk that has been accumulated up to time t)
H <- fit$cumhaz

# Plot this quantity
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           palette  = brewer.pal(n = 1, name = "Paired")[3],
           fun='cumhaz',
           title="Cumulative Hazard Curve for Cholangiocarcinoma - OS")

##### KAPLAN-MEIER CURVES BETWEEN GROUPS ####
# Test on clinical variables to assess which of them is influencing the time-to-event outcome
# We flag as YES (significant) variables with p-values of this test <=0.05
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$SEX, data = main)  # NO
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$HCV, data = main)  # NO
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$HBV, data = main)  # NO
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$Ca19.9.55, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$CHEMIOTERAPIA.NEOADIUVANTE, data = main)  # NO
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$Major.Hepatectomy, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$RESEZIONE.VIA.BILIARE, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$LINFOADENECTOMIA, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$COMPLICANZE.SEVERE, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$CIRROSI, data = main)  # NO
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$PATTERN, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$SINGLE.NODULE, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$T.VIII.ed, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$N, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$GRADING, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$R.status, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$INVASIONE.VASCOLARE.MACROSCOPICA, data = main)  # NO
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$INVASIONE.VASCOLARE.MICROSCOPICA, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$NODULI.SATELLITI, data = main)  # YES
survdiff(Surv(main$OS..Days., main$STATO.VIVO.MORTO) ~ main$CHEMIOTERAPIA.ADIUVANTE, data = main)  # NO

# Test on CORE radiomic variables to assess which of them is influencing the time-to-event outcome
# We construct the dataset made by the time variable, the censor variable and all the core radiomics variables
main_pc_test <- data.frame(OS..Days. = main$OS..Days.,STATO.VIVO.MORTO = main$STATO.VIVO.MORTO, pc_radiomics)

n <- dim(main_pc_test)[2]
for(i in 3:n){
  limit <- c(min(main_pc_test[,i]), median(main_pc_test[,i]), max(main_pc_test[,i]))
  Dim_cut <- cut(main_pc_test[,i], breaks=limit, labels=c("low", "high"))
  print(colnames(main_pc_test[i]))
  print(survdiff( Surv(main_pc_test$OS..Days., main_pc_test$STATO.VIVO.MORTO) ~ Dim_cut , data=main_pc_test))
  
}
# Again, we consider as significant (YES) variables with p-value of this test <0.05

## CONVENTIONAL_HUmin NO
## CONVENTIONAL_HUstd NO
## CONVENTIONAL_HUmax NO
## CONVENTIONAL_HUQ2 YES
## CONVENTIONAL_HUSkewness NO
## CONVENTIONAL_HUKurtosis NO
## DISCRETIZED_HISTO_Entropy_log10 NO
## SHAPE_Volume.mL. YES
## SHAPE_Sphericity.onlyFor3DROI.. NO
## SHAPE_Compacity.onlyFor3DROI. YES
## GLCM_Correlation NO
## GLRLM_LRE NO
## GLRLM_SRLGE NO
## GLRLM_SRHGE YES
## GLRLM_GLNU YES
## NGLDM_Coarseness YES
## NGLDM_Contrast NO
## NGLDM_Busyness YES
## GLZLM_SZE NO
## GLZLM_SZLGE NO
## GLZLM_SZHGE NO
## GLZLM_LZLGE YES

# 8 out of 22 significantly different radiomic variables



# Test on MARGIN radiomic variables to assess which of them is influencing the time-to-event outcome
# We construct the dataset made by the time variable, the censor variable and all the core radiomics variables
main_test <- main[-256,]
pm_rad_test <- data.frame(OS..Days. = main_test$OS..Days.,STATO.VIVO.MORTO = main_test$STATO.VIVO.MORTO, pm_radiomics)

n <- dim(pm_rad_test)[2]
for(i in 3:n){
  limit <- c(min(pm_rad_test[,i]), median(pm_rad_test[,i]), max(pm_rad_test[,i]))
  Dim_cut <- cut(pm_rad_test[,i], breaks=limit, labels=c("low", "high"))
  print(colnames(pm_rad_test[i]))
  print(survdiff( Surv(pm_rad_test$OS..Days., pm_rad_test$STATO.VIVO.MORTO) ~ Dim_cut , data=pm_rad_test))
  
}
# Again, we consider as significant (YES) variables with p-value of this test <0.05

## CONVENTIONAL_HUmin NO
## CONVENTIONAL_HUstd NO
## CONVENTIONAL_HUmax NO
## CONVENTIONAL_HUQ2 NO
## CONVENTIONAL_HUSkewness NO
## CONVENTIONAL_HUKurtosis NO
## DISCRETIZED_HISTO_Entropy_log10 NO
## SHAPE_Volume.mL. YES
## SHAPE_Sphericity.onlyFor3DROI.. YES
## SHAPE_Compacity.onlyFor3DROI. NO
## GLCM_Correlation NO
## GLRLM_LRE NO
## GLRLM_SRLGE NO
## GLRLM_SRHGE NO
## GLRLM_GLNU YES
## NGLDM_Coarseness YES
## NGLDM_Contrast YES
## NGLDM_Busyness YES 
## GLZLM_SZE NO
## GLZLM_SZLGE NO
## GLZLM_SZHGE NO
## GLZLM_LZLGE YES (0.05)

# 7 out of 22 significantly different radiomic variables


##### COX MODEL - MAIN + CORE #####
# Fit the Cox model (starting by including only the significant variables found in the test in KAPLAN-MEIER CURVES BETWEEN GROUPS)
# for computation reasons
mod.cox <- coxph( surv_obj ~ AGE + Ca.19.9 + Major.Hepatectomy + RESEZIONE.VIA.BILIARE +
                    LINFOADENECTOMIA + COMPLICANZE.SEVERE + Degenza + PATTERN + DIMENSIONE.MAX.MM +
                    SINGLE.NODULE + NUMERO..NO.SATELLITI. + NUMERO.LINFONODI.ASPORTATI +
                    NUMERO.LINFONODI.METASTATICI + T.VIII.ed + N + GRADING + R.status + 
                    INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI +
                    CONVENTIONAL_HUQ2 + SHAPE_Compacity.onlyFor3DROI. +
                    GLRLM_SRHGE + GLRLM_GLNU + NGLDM_Coarseness + NGLDM_Busyness +
                    GLZLM_LZLGE,
                    data = main_pc)
summary(mod.cox)

## Stepwise procedure to select the best model
m_null <- coxph( surv_obj ~ 1, data = main_pc)
step(m_null, trace = F, scope = list(lower=formula(m_null), upper=formula(mod.cox)),
     direction = 'both', data =main_pc)

cox.reduced <- coxph(formula = surv_obj ~ GLRLM_SRHGE + AGE + NODULI.SATELLITI + 
                       COMPLICANZE.SEVERE + Ca.19.9 + SINGLE.NODULE + Major.Hepatectomy + 
                       GLRLM_GLNU, 
                     data = main_pc)
# We remove R.status even if the step suggest it because the test for proprotional hazard assumption
# for this covariate is not verified

summary(cox.reduced)
anova(cox.reduced, mod.cox, test = 'Chisq') # OK, the two models are not significantly different

## Some plots
x11()
ggforest(cox.reduced, data=main_pc)   # HRs and theirs CIs

x11()
ggcoxdiagnostics(cox.reduced, type = "martingale")  # Martingale Residuals have 0 mean along time
# roughly symmetrically distributed about zero with a standard deviation of 1

x11()
ggcoxdiagnostics(cox.reduced, type = "deviance")  # check outliers by visualizing the deviance residuals.
# We see symmetry around 0, this is good

x11()
# Shoenfeld residuals represent the difference 
# between the observed covariate and the expected given the risk set at that time. 
ggcoxdiagnostics(cox.reduced, type = "schoenfeld")   # They should be flat, centred about zero.

test.ph <- cox.zph(cox.reduced)
test.ph    # all p-values are high, the porportional hazard assumption is verified for all covariates



##### COX MODEL - MAIN + CORE + MARGIN #####
# Creating the survival function for the margin
surv_obj_m <- Surv(main_pc_pm$OS..Days., main_pc_pm$STATO.VIVO.MORTO)

# Fit the Cox model (starting by including only the significant variables found in the test in KAPLAN-MEIER CURVES BETWEEN GROUPS)
# for computation reasons
mod.cox_all <- coxph( surv_obj_m ~ AGE + Ca.19.9 + Major.Hepatectomy + RESEZIONE.VIA.BILIARE +
                        LINFOADENECTOMIA + COMPLICANZE.SEVERE + Degenza + PATTERN + DIMENSIONE.MAX.MM +
                        SINGLE.NODULE + NUMERO..NO.SATELLITI. + NUMERO.LINFONODI.ASPORTATI +
                        NUMERO.LINFONODI.METASTATICI + T.VIII.ed + N + GRADING + R.status + 
                        INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI +
                        CONVENTIONAL_HUQ2 + SHAPE_Compacity.onlyFor3DROI. +
                        GLRLM_SRHGE + GLRLM_GLNU + NGLDM_Coarseness + NGLDM_Busyness +
                        GLZLM_LZLGE +
                        SHAPE_Volume.mL._margin + SHAPE_Sphericity.onlyFor3DROI.._margin +
                        GLRLM_GLNU_margin + NGLDM_Coarseness_margin + NGLDM_Contrast_margin +
                        NGLDM_Busyness_margin + GLZLM_LZLGE_margin,
                      data = main_pc_pm)
summary(mod.cox)

## Stepwise procedure to select the best model
m_null <- coxph( surv_obj_m ~ 1, data = main_pc_pm)
step(m_null, trace = F, scope = list(lower=formula(m_null), upper=formula(mod.cox_all)),
     direction = 'both', data = main_pc_pm)

cox_all.reduced <- coxph(formula = surv_obj_m ~ GLRLM_SRHGE + AGE + NODULI.SATELLITI + 
                           COMPLICANZE.SEVERE + SINGLE.NODULE + Major.Hepatectomy + 
                           R.status + GLRLM_GLNU_margin + NUMERO.LINFONODI.ASPORTATI + 
                           N + SHAPE_Sphericity.onlyFor3DROI.._margin + NGLDM_Coarseness_margin + 
                           GRADING + RESEZIONE.VIA.BILIARE + GLZLM_LZLGE_margin, 
                         data = main_pc_pm)
# We remove NGLDM_Busyness_margin even if the step suggest it because the test for proprotional hazard assumption
# for this covariate is not verified
summary(cox_all.reduced)
anova(cox_all.reduced, mod.cox_all, test = 'Chisq') # OK, the two models are not significantly different


## Some plots
x11()
ggforest(cox_all.reduced, data=main_pc_pm) # HRs and theirs CIs

x11()
ggcoxdiagnostics(cox_all.reduced, type = "martingale")
# Simmetry around 0

x11()
ggcoxdiagnostics(cox_all.reduced, type = "deviance")
# Simmetry around 0

x11()
# Shoenfeld residuals represent the difference 
# between the observed covariate and the expected given the risk set at that time. 
ggcoxdiagnostics(cox_all.reduced, type = "schoenfeld")

test.ph <- cox.zph(cox_all.reduced)
test.ph   # all p-values are high, the porportional hazard assumption is verified for all covariates


# Estimated Baseline Survival Curves  
fit<-survfit(cox_all.reduced, data=main_pc_pm)
x11()
plot(fit, 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()


#### Estimated Adjusted Survival Curves ####
# We plot these curves for the reduced model with core and margin
# We plot only the curves observing the variables which are risk/protective factors in two particular cases 


#### GLRLM_SRHGE ####
# We construct the dataset where GLRLM_SRHGE has three possible values (1st quantile, median and 3rd quantile) 
# and the other variables are fixed at the median (for continuous variables)
# or at the mode (for categorical variables)
df_GLRLM_SRHGE <- with(main_pc_pm,
                       data.frame(NODULI.SATELLITI = rep(0,3),
                                  AGE = rep(median(main_pc_pm$AGE, na.rm=T),3),
                                  COMPLICANZE.SEVERE = rep(0,3),
                                  R.status = rep(0,3),
                                  NUMERO.LINFONODI.ASPORTATI = rep(median(main_pc_pm$NUMERO.LINFONODI.ASPORTATI, na.rm=T),3),
                                  SINGLE.NODULE = rep(1,3),
                                  Major.Hepatectomy = rep(0,3),
                                  N = as.factor(rep(0,3)),
                                  GRADING = rep(2,3),
                                  RESEZIONE.VIA.BILIARE = rep(0,3),
                                  GLRLM_SRHGE = c(summary(main_pc_pm$GLRLM_SRHGE)[[2]],summary(main_pc_pm$GLRLM_SRHGE)[[3]],summary(main_pc_pm$GLRLM_SRHGE)[[5]]),
                                  GLRLM_GLNU_margin = rep(median(main_pc_pm$GLRLM_GLNU_margin, na.rm=T),3),
                                  SHAPE_Sphericity.onlyFor3DROI.._margin = rep(median(main_pc_pm$SHAPE_Sphericity.onlyFor3DROI.._margin, na.rm=T),3),
                                  NGLDM_Coarseness_margin = rep(median(main_pc_pm$NGLDM_Coarseness_margin, na.rm=T),3),
                                  GLZLM_LZLGE_margin = rep(median(main_pc_pm$GLZLM_LZLGE_margin, na.rm=T),3)
                                  ))
df_GLRLM_SRHGE

# Constructing the survival object
fit.df_GLRLM_SRHGE <- survfit(cox_all.reduced, newdata = df_GLRLM_SRHGE)
fit.df_GLRLM_SRHGE

# Plotting the three curves together
ggsurvplot(fit.df_GLRLM_SRHGE, fit.df_GLRLM_SRHGE, censor=F,
           palette  = brewer.pal(n = 1, name = "Accent"))

# Not significantly different...




#### SINGLE-NODULE - SIGNIFICANTLY DIFFERENT CURVES ####
SINGLE.NODULE <- as.factor(SINGLE.NODULE)

# We construct the dataset where SINGLE.NODULE is 0 or 1 and the other variables are fixed at the median (for continuous variables)
# or at the mode (for categorical variables)
df_SN <- with(main_pc_pm,
              data.frame(NODULI.SATELLITI = rep(0,2),
                         AGE = rep(median(main_pc_pm$AGE, na.rm=T),2),
                         COMPLICANZE.SEVERE = rep(0,2),
                         R.status = rep(0,2),
                         NUMERO.LINFONODI.ASPORTATI = rep(median(main_pc_pm$NUMERO.LINFONODI.ASPORTATI, na.rm=T),2),
                         SINGLE.NODULE = c(0,1),
                         Major.Hepatectomy = rep(0,2),
                         N = as.factor(rep(0,2)),
                         GRADING = rep(2,2),
                         RESEZIONE.VIA.BILIARE = rep(0,2),
                         GLRLM_SRHGE = rep(median(main_pc_pm$GLRLM_SRHGE, na.rm=T),2),
                         GLRLM_GLNU_margin = rep(median(main_pc_pm$GLRLM_GLNU_margin, na.rm=T),2),
                         SHAPE_Sphericity.onlyFor3DROI.._margin = rep(median(main_pc_pm$SHAPE_Sphericity.onlyFor3DROI.._margin, na.rm=T),2),
                         NGLDM_Coarseness_margin = rep(median(main_pc_pm$NGLDM_Coarseness_margin, na.rm=T),2),
                         GLZLM_LZLGE_margin = rep(median(main_pc_pm$GLZLM_LZLGE_margin, na.rm=T),2)
              ))

df_SN

# Constructing the survival object
fit.df_SN <- survfit(cox_all.reduced, newdata = df_SN)
fit.df_SN

# Plotting the two curves together
ggsurvplot(fit.df_SN, fit.df_SN, censor=F,
           palette  = brewer.pal(n = 1, name = "Set2"))

# The two curves look significantly different !



