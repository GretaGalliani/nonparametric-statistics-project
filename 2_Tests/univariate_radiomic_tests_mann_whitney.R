# GOAL: test each radiomic features to look for differences in distributions 
# among the populations induced by the grouping structures of the outcome
# MANN-WHITNEY TEST APPROACH

set.seed(24091998)
alpha <- 0.05

# importing datasets
main <- read.table("Datasets/main.txt")

# now we import the dataset after robust correlation reduction 
pc <-read.table("Datasets/pc_post_robcor.txt")
pm <-read.table("Datasets/pm_post_robcor.txt")
# discarding labels which we are not considering in these tests
pc <- pc[,-c(1,2)]
pm <- pm[,-c(1,2)]

# removing outliers found by robust analysis
pc <- pc[-c(102,122,204,235,241),]
pm <- pm[-c(102,122,204,235,241),]
main <-main[-c(102,122,204,235,241),]

p <- dim(pc)[2]

# Now we perform the test for each of the 4 different outcomes (vascular invasion, grading, alive/dead status, relapse)
# For each outcome, we tests all the features of portal core and then portal margin looking for differences
# We apply B&H correction on the p-value vector to deal with false discoveries

#### VASCULAR INVASION ----
# We divide the two populations 
pc_ivm_t <- pc[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1),]
pc_ivm_f <- pc[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0),]
pm_ivm_t <- pm[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1),]
pm_ivm_f <- pm[which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0),]

# Testing the 22 features of portal core 
p.val.ivm.pc = rep(1, p)
for (i in (1:p)){
  p.val.ivm.pc[i] <-  wilcox.test(pc_ivm_t[,i],pc_ivm_f[,i])$p.value
}

# Correcting the p-values
p.val.ivm.pc <- p.adjust(p.val.ivm.pc, method='BH')

# Showing the results
which(p.val.ivm.pc<alpha) 
length(which(p.val.ivm.pc<alpha)) # 16 significantly different variables 
                                  # they were 21 before robust analysis 



# Testing the 22 features of portal margin 
p.val.ivm.pm = rep(0, p)

for (i in (1:p)){
  p.val.ivm.pm[i] <- wilcox.test(pm_ivm_t[,i],pm_ivm_f[,i])$p.value
}

# Correcting the p-values
p.val.ivm.pm <- p.adjust(p.val.ivm.pm, method='BH')

# Showing the results
which(p.val.ivm.pm<alpha)
length(which(p.val.ivm.pm<alpha)) # 14 significantly different variables 
                                  # they were 18 before robust analysis 


#### GRADING ----
# Notice that for the case of the grading we have tested as populations the one 
# with grading 1 or 2 versus the one with grading 3 

# We divide the two populations 
pc_grad_t <- pc[which(main$GRADING=="3"),]
pc_grad_f <- pc[which(main$GRADING=="1"|main$GRADING=="2"),]
pm_grad_t <- pm[which(main$GRADING=="3"),]
pm_grad_f <- pm[which(main$GRADING=="1"|main$GRADING=="2"),]

# Testing the 22 features of portal core 
p.val.grad.pc = rep(0, p)

for (i in (1:p)){
  p.val.grad.pc[i] <- wilcox.test(pc_grad_t[,i],pc_grad_f[,i])$p.value

}

# Correcting the p-values
p.val.grad.pc <- p.adjust(p.val.grad.pc, method='BH')

# Showing the results
which(p.val.grad.pc<alpha) # 0 significantly different variables 


# Testing the 22 features of portal margin 
p.val.grad.pm = rep(0, p)

for (i in (1:p)){
  p.val.grad.pm[i] <- wilcox.test(pm_grad_t[,i],pm_grad_f[,i])$p.value
}

# Correcting the p-values
p.val.grad.pm <- p.adjust(p.val.grad.pm, method='BH')

# Showing the results
which(p.val.grad.pm<alpha) 
length(which(p.val.grad.pm<alpha)) # 7 significantly different variables 
                                  # they were 12 before robust analysis 



#### DEAD/ALIVE STATUS ----

# We divide the two populations 
pc_da_t <- pc[which(main$STATO.VIVO.MORTO==1),]
pc_da_f <- pc[which(main$STATO.VIVO.MORTO==0),]
pm_da_t <- pm[which(main$STATO.VIVO.MORTO==1),]
pm_da_f <- pm[which(main$STATO.VIVO.MORTO==0),]

# Testing the 22 features of portal core 
p.val.da.pc = rep(0, p)

for (i in (1:p)){
  p.val.da.pc[i] <- wilcox.test(pc_da_t[,i],pc_da_f[,i])$p.value
}

# Correcting the p-values
p.val.da.pc <- p.adjust(p.val.da.pc, method='BH')

# Showing the results
which(p.val.da.pc<alpha)
length(which(p.val.da.pc<alpha)) # 11 significantly different variables 
                                 # they were 17 before robust analysis 



# Testing the 22 features of portal margin 
p.val.da.pm = rep(0, p)

for (i in (1:p)){
  p.val.da.pm[i] <- wilcox.test(pm_da_t[,i],pm_da_f[,i])$p.value
}

# Correcting the p-values
p.val.da.pm <- p.adjust(p.val.da.pm, method='BH')

# Showing the results
which(p.val.da.pm<alpha) # 1 significantly different variables 
                         # they were 3 before robust analysis 



#### RELAPSE ----

# We divide the two populations 
pc_rel_t <- pc[which(main$RECIDIVA==1),]
pc_rel_f <- pc[which(main$RECIDIVA==0),]
pm_rel_t <- pm[which(main$RECIDIVA==1),]
pm_rel_f <- pm[which(main$RECIDIVA==0),]

# Testing the 22 features of portal core 
p.val.rel.pc = rep(0, p)

for (i in (1:p)){
  p.val.rel.pc[i] <- wilcox.test(pc_rel_t[,i],pc_rel_f[,i])$p.value
}

# Correcting the p-values
p.val.rel.pc <- p.adjust(p.val.rel.pc, method='BH')

# Showing the results
which(p.val.rel.pc<alpha)
length(which(p.val.rel.pc<alpha)) # 7 significantly different variables 
                                  # they were 12 before robust analysis 


# Testing the 22 features of portal margin 
p.val.rel.pm = rep(0, p)

for (i in (1:p)){
  p.val.rel.pm[i] <- wilcox.test(pm_da_t[,i],pm_da_f[,i])$p.value
}

# Correcting the p-values
p.val.rel.pm <- p.adjust(p.val.da.pm, method='BH')

# Showing the results
which(p.val.rel.pm<alpha) # 0 significantly different variables


