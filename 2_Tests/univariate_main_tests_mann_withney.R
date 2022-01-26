# GOAL: test each radiomic features to look for differences in distributions 
# among the populations induced by the grouping structures of the outcome
# MANN-WHITNEY TEST APPROACH

###
#For IVM: 8 (10 11 12 21 22 24 25 28) variables significantly different
#For GRADING: 3 (15 23 24) variables significantly different
#For Dead/alive: 16 (1  5  6 10 11 14 15 17 18 24 25 26 27 29 30 31) variables significantly different
#For Relapse: 3 (21 29 31) variables significantly different
seed <- 24091998
alpha <- 0.05

# importing datasets
main <- read.table("Datasets/main.txt")

# removing outliers found by robust analysis
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

#Omit the NA
main <-na.omit(main)

p <- dim(main)[2]
p <- p-1 #since i have to delete the outcome 

# Now we perform the test for each of the 4 different outcomes (vascular invasion, grading, alive/dead status, relapse)
# For each outcome, we tests all the features looking for differences
# We apply B&H correction on the p-value vector to deal with false discoveries

#### VASCULAR INVASION ----

p.val.ivm.main = rep(0, p)

#Create a dataset without the outcomes of interest
main1 <- main
main1$INVASIONE.VASCOLARE.MICROSCOPICA <- NULL

pc_ivm_t <- main[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1),]
pc_ivm_f <- main[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0),]

pp = rep(0,p)
for (i in (1:p)){
  p.val.ivm.main[i] <-  wilcox.test(as.numeric(unlist(main1[,i])) ~ main$INVASIONE.VASCOLARE.MICROSCOPICA)$p.value
}
which(p.val.ivm.main<alpha)#12

# Correcting the p-values
p.val.ivm.main <- p.adjust(p.val.ivm.main, method='BH')

# Showing the results
which(p.val.ivm.main<alpha) #  10 11 12 21 22 24 25 28 significantly different variables 
names(main1[,c(10,11, 12,21, 22, 24, 25,28)])

#### GRADING ----
# Notice that for the case of the grading we have tested as populations the one with grading 1 or 2 versus the one with grading 3 
p.val.gra.main = rep(0, p)

#Create a dataset without the outcomes of interest
main1 <- main
main1$GRADING <- NULL
main1$GRADING <- ifelse(main$GRADING==1 | main$GRADING==2,0,1)

main2 <- main1
main2$GRADING <- NULL
pp = rep(0,p)
for (i in (1:p)){
  p.val.gra.main[i] <-  wilcox.test(as.numeric(unlist(main2[,i])) ~ main1$GRADING)$p.value
}
which(p.val.gra.main<alpha)#7
# Correcting the p-values
p.val.gra.main <- p.adjust(p.val.gra.main, method='BH')

# Showing the results
which(p.val.gra.main<alpha) # 15 23 24 significantly different variables 

names(main1[,c(15, 23, 24)])

#### DEAD/ALIVE STATUS ----
p.val.da.main = rep(0, p)

#Create a dataset without the outcomes of interest
main1 <- main
main1$STATO.VIVO.MORTO <- NULL

pp = rep(0,p)
for (i in (1:p)){
  p.val.da.main[i] <-  wilcox.test(as.numeric(unlist(main1[,i])) ~ main$STATO.VIVO.MORTO)$p.value
}
which(p.val.da.main<alpha)#18
# Correcting the p-values
p.val.da.main <- p.adjust(p.val.da.main, method='BH')

# Showing the results
which(p.val.da.main<alpha) # 1  5  6 10 11 14 15 17 18 24 25 26 27 29 30 31 significantly different variables 


names(main1[,c(1 , 5  ,6 ,10, 11, 14, 15, 17, 18, 24, 25, 26, 27, 29, 30, 31)])

#### RELAPSE ----
p.val.re.main = rep(0, p)

#Create a dataset without the outcomes of interest
main1 <- main
main1$RECIDIVA <- NULL

pp = rep(0,p)
for (i in (1:p)){
  p.val.re.main[i] <-  wilcox.test(as.numeric(unlist(main1[,i])) ~ main$RECIDIVA)$p.value
}
which(p.val.re.main<alpha) 
# Correcting the p-values
p.val.re.main <- p.adjust(p.val.re.main, method='BH')

# Showing the results
which(p.val.re.main<alpha) # 21 29 31 significantly different variables 
names(main1[,c(21,29,31)])
