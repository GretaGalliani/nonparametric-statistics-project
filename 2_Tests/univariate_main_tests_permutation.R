# GOAL: test each clinical features to look for differences in distributions 
# among the populations induced by the grouping structures of the outcome
# PERMUTATIONAL APPROACH

##COMMENT:
#for IVM: 7 (10 11 12 21 22 24 25) variables significantly different
#for GRADING: 3 (15 23 24) variables significantly different
#for STATO.VIVO.MORTO: 12 (1  5  6 10 14 15 24 25 26 29 30 31) variables significantly different
#for RECIDIVA: 3 (21 29 31) variables significantly different
library(progress)
source("2_Tests/test_functions.R")

seed <- 24091998
alpha <- 0.05

# importing datasets
main <- read.table("Datasets/main.txt")
main <- main[-c(204,235,241,122),]

#Delete the center and codice.paziente
main$CENTRO<-NULL
main$Codice.PAZ <-NULL


#Delete this convariates since are categorical
main$T.VIII.ed <- NULL
main$N <- NULL
main$NUMERO.LINFONODI.METASTATICI <- NULL

#Convert the sex, male --> 0, female -->1
main$SEX <-ifelse(main$SEX=='F',1,0)

#Discard observation with NA
main <-na.omit(main)

p <- dim(main)[2]
p <- p-1 #since i have to delete the outcome 

# Now we perform the test for each of the 4 different outcomes (vascular invasion, grading, alive/dead status, relapse)
# For each outcome, we tests all the clinical features looking for differences
# We apply B&H correction on the p-value vector to deal with false discoveries

#### VASCULAR INVASION ----

# We divide the two populations 
main_ivm_t <- main[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1),]
main_ivm_f <- main[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0),]

#Delete the outcomes IVM
main_ivm_t$INVASIONE.VASCOLARE.MICROSCOPICA <- NULL
main_ivm_f$INVASIONE.VASCOLARE.MICROSCOPICA <- NULL


p.val.ivm.main = rep(0, p)

for (i in (1:p)){
  p.val.ivm.main[i] <- perm_mean_test(main_ivm_t[,i],main_ivm_f[,i])
}
which(p.val.ivm.main<alpha)

# Correcting the p-values
p.val.ivm.main <- p.adjust(p.val.ivm.main, method='BH')

# Showing the results
which(p.val.ivm.main<alpha) # 10 11 12 21 22 24 25 --> variables significantly different  

names(main_ivm_t[,c(10,11, 12, 21,22, 24, 25)])


#### GRADING ----
# Notice that for the case of the grading we have tested as populations the one with grading 1 or 2 versus the one with grading 3 

# We divide the two populations 
main_grad_t <- main[which(main$GRADING=="3"),]
main_grad_f <- main[which(main$GRADING=="1"|main$GRADING=="2"),]

#Delete the outcomes GRADING
main_grad_f$GRADING <- NULL
main_grad_t$GRADING <- NULL


p.val.grad.main = rep(0, p)

for (i in (1:p)){
  p.val.grad.main[i] <- perm_mean_test(main_grad_t[,i],main_grad_f[,i])
}
which(p.val.grad.main<alpha)#6

# Correcting the p-values
p.val.grad.main <- p.adjust(p.val.grad.main, method='BH')

# Showing the results
which(p.val.grad.main<alpha) # 15 23 24 -> variables significantly different  

names(main_grad_t[,c(15 ,23, 24)])


#### DEAD/ALIVE STATUS ----

# We divide the two populations 
main_da_t <- main[which(main$STATO.VIVO.MORTO==1),]
main_da_f <- main[which(main$STATO.VIVO.MORTO==0),]

#Detele the oucome STATO VIVO MORTO
main_da_t$STATO.VIVO.MORTO <- NULL
main_da_f$STATO.VIVO.MORTO <- NULL

p.val.da.main = rep(0, p)

for (i in (1:p)){
  p.val.da.main[i] <- perm_mean_test(main_da_t[,i],main_da_f[,i])
}

which(p.val.da.main<alpha)#17

# Correcting the p-values
p.val.da.main <- p.adjust(p.val.da.main, method='BH')

# Showing the results
which(p.val.da.main<alpha) # 1  5  6 10 14 15 24 25 26 29 30 31 --> variables significantly different  

names(main_da_t[,c(1 , 5  ,6, 10, 14, 15, 24 ,25, 26, 29, 30, 31)])



#### RELAPSE ----

# We divide the two populations 
main_rel_t <- main[which(main$RECIDIVA==1),]
main_rel_f <- main[which(main$RECIDIVA==0),]

#Delete the outcome RECIDIVA
main_rel_f$RECIDIVA <- NULL
main_rel_t$RECIDIVA <-NULL

# Testing the 30 features of portal core 
p.val.rel.main = rep(0, p)

for (i in (1:p)){
  p.val.rel.main[i] <- perm_mean_test(main_rel_t[,i],main_rel_f[,i])
}

which(p.val.rel.main<alpha)#7

# Correcting the p-values
p.val.rel.main <- p.adjust(p.val.rel.main, method='BH')

# Showing the results
which(p.val.rel.main<alpha) # 21 29 31 --> variables significantly different 

names(main_rel_t[,c(21, 29, 31)])



 
