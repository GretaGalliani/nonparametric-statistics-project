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
surv_obj <- Surv(main_pc$RFS..Days., main_pc$RECIDIVA)

# Computing the Kaplan-Meier estimator
fit <- survfit(surv_obj  ~ 1, data = main_pc)

# Reporting the table for the KM estimator
kable(head(tidy(fit),20))

# The medial survival time is 
surv_median(fit) # 450 days

##### KAPLAN-MEYER PLOTS #####
plot(fit, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='steelblue2',
     main="Kaplan-Meier Curve for Cholangiocarcinoma - RFS")

# Another possible plot
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette  = brewer.pal(n = 1, name = "Set2"),
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90 # title="Kaplan-Meier Curve for Cholangiocarcinoma - RFS"
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
           title="Cumulative Incidence Curve for Cholangiocarcinoma - RFS")

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
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$AGE, data = main) # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$SEX, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$HCV, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$HBV, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$Ca19.9.55, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$CHEMIOTERAPIA.NEOADIUVANTE, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$Major.Hepatectomy, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$RESEZIONE.VIA.BILIARE, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$LINFOADENECTOMIA, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$COMPLICANZE.SEVERE, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$CIRROSI, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$PATTERN, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$SINGLE.NODULE, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$T.VIII.ed, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$N, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$GRADING, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$R.status, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$INVASIONE.VASCOLARE.MACROSCOPICA, data = main)  # NO
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$INVASIONE.VASCOLARE.MICROSCOPICA, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$NODULI.SATELLITI, data = main)  # YES
survdiff(Surv(main$RFS..Days., main$RECIDIVA) ~ main$CHEMIOTERAPIA.ADIUVANTE, data = main)  # YES


# Test on CORE radiomic variables to assess which of them is influencing the time-to-event outcome
# We construct the dataset made by the time variable, the censor variable and all the core radiomics variables
main_pc_test <- data.frame(RFS..Days. = main$RFS..Days.,RECIDIVA = main$RECIDIVA, pc_radiomics)

n <- dim(main_pc_test)[2]
for(i in 3:n){
  limit <- c(min(main_pc_test[,i]), median(main_pc_test[,i]), max(main_pc_test[,i]))
  Dim_cut <- cut(main_pc_test[,i], breaks=limit, labels=c("low", "high"))
  print(colnames(main_pc_test[i]))
  print(survdiff( Surv(main_pc_test$RFS..Days., main_pc_test$RECIDIVA) ~ Dim_cut , data=main_pc_test))
  
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


# Test on MARGIN radiomic variables to assess which of them is influencing the time-to-event outcome
# We construct the dataset made by the time variable, the censor variable and all the core radiomics variables
main_test <- main[-256,]

pm_rad_test <- data.frame(RFS..Days. = main_test$RFS..Days.,RECIDIVA = main_test$RECIDIVA, pm_radiomics)

n <- dim(pm_rad_test)[2]
for(i in 3:n){
  limit <- c(min(pm_rad_test[,i]), median(pm_rad_test[,i]), max(pm_rad_test[,i]))
  Dim_cut <- cut(pm_rad_test[,i], breaks=limit, labels=c("low", "high"))
  print(colnames(pm_rad_test[i]))
  print(survdiff( Surv(pm_rad_test$RFS..Days., pm_rad_test$RECIDIVA) ~ Dim_cut , data=pm_rad_test))
  
}
# Again, we consider as significant (YES) variables with p-value of this test <0.05

## CONVENTIONAL_HUmin YES
## CONVENTIONAL_HUstd YES
## CONVENTIONAL_HUmax NO
## CONVENTIONAL_HUQ2 NO
## CONVENTIONAL_HUSkewness YES
## CONVENTIONAL_HUKurtosis NO
## DISCRETIZED_HISTO_Entropy_log10 YES
## SHAPE_Volume.mL. YES
## SHAPE_Sphericity.onlyFor3DROI.. YES
## SHAPE_Compacity.onlyFor3DROI. NO
## GLCM_Correlation NO
## GLRLM_LRE NO
## GLRLM_SRLGE NO
## GLRLM_SRHGE NO
## GLRLM_GLNU YES
## NGLDM_Coarseness YES
## NGLDM_Contrast NO
## NGLDM_Busyness YES
## GLZLM_SZE NO
## GLZLM_SZLGE NO
## GLZLM_SZHGE NO
## GLZLM_LZLGE YES


##### COX MODEL - MAIN + CORE #####
# Fit the Cox model (starting by including only the significant variables found in the test in KAPLAN-MEIER CURVES BETWEEN GROUPS)
# for computation reasons
mod.cox <- coxph( surv_obj ~ AGE + Ca.19.9  + 
                    COMPLICANZE.SEVERE + PATTERN + SINGLE.NODULE + T.VIII.ed + N +  R.status + 
                    INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + CHEMIOTERAPIA.ADIUVANTE + 
                    CONVENTIONAL_HUQ2 + SHAPE_Volume.mL. + SHAPE_Compacity.onlyFor3DROI. + 
                    GLRLM_SRHGE + GLRLM_GLNU + NGLDM_Coarseness + NGLDM_Busyness +
                    GLZLM_LZLGE,
                  data = main_pc)
summary(mod.cox)

## Stepwise procedure to select the best model
m_null <- coxph( surv_obj ~ 1, data = main_pc)
step(m_null, trace = F, scope = list(lower=formula(m_null), upper=formula(mod.cox)),
     direction = 'both', data =main_pc)

cox.reduced <- coxph(formula = surv_obj ~ CONVENTIONAL_HUQ2 + R.status + SINGLE.NODULE + 
                       COMPLICANZE.SEVERE + GLZLM_LZLGE + N, 
                     data = main_pc)
# We remove Ca.19.9 even if the step suggest it because the test for proprotional hazard assumption
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
surv_obj_m <- Surv(main_pc_pm$RFS..Days., main_pc_pm$RECIDIVA)

# Fit the Cox model (starting by including only the significant variables found in the test in KAPLAN-MEIER CURVES BETWEEN GROUPS)
# for computation reasons
mod.cox_all <- coxph(surv_obj_m ~ AGE + Ca.19.9  + 
                        COMPLICANZE.SEVERE + PATTERN + SINGLE.NODULE + T.VIII.ed + N +  R.status + 
                        INVASIONE.VASCOLARE.MICROSCOPICA + NODULI.SATELLITI + CHEMIOTERAPIA.ADIUVANTE + 
                        CONVENTIONAL_HUQ2 + SHAPE_Volume.mL. + SHAPE_Compacity.onlyFor3DROI. + 
                        GLRLM_SRHGE + GLRLM_GLNU + NGLDM_Coarseness + NGLDM_Busyness +
                        GLZLM_LZLGE + CONVENTIONAL_HUmin_margin + CONVENTIONAL_HUstd_margin +
                        CONVENTIONAL_HUSkewness_margin + DISCRETIZED_HISTO_Entropy_log10_margin +
                        SHAPE_Volume.mL._margin + SHAPE_Sphericity.onlyFor3DROI.._margin +
                        GLRLM_GLNU_margin + NGLDM_Coarseness_margin + NGLDM_Busyness_margin +
                        GLZLM_LZLGE_margin,
                      data = main_pc_pm)
summary(mod.cox)


## Stepwise procedure to select the best model
m_null <- coxph( surv_obj_m ~ 1, data = main_pc_pm)
step(m_null, trace = F, scope = list(lower=formula(m_null), upper=formula(mod.cox_all)),
     direction = 'both', data = main_pc_pm)

cox_all.reduced <- coxph(formula = surv_obj_m ~ CONVENTIONAL_HUQ2 + R.status + SINGLE.NODULE + 
                           GLZLM_LZLGE + GLZLM_LZLGE_margin + DISCRETIZED_HISTO_Entropy_log10_margin + 
                           GLRLM_GLNU_margin + Ca.19.9 + NODULI.SATELLITI, 
                         data = main_pc_pm)
# We remove COMPLICANZE.SEVERE even if the step suggest it because the test for proprotional hazard assumption
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


