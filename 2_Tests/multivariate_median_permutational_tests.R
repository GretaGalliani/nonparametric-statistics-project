# GOAL: test  radiomic features to look for differences in distributions 
# among the populations induced by the grouping structures of the outcome
# This time we are going to apply a multivariate approach
# MULTIVARIATE PERMUTATIONAL TEST APPROACH

library(DepthProc)
library(progress)

# loading test functions
source("2_Tests/test_functions.R")

set.seed(13021998)

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

# We are going to perform a two-sample independent test in the multivariate setting
# the two populations are those with positive/negative outcome


# Start with Portal Core

#### IVM ####
idx_ivm_yes <- which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1)   # rows of the patients with IVM
idx_ivm_no <- which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0)

ivm_yes <- pc[idx_ivm_yes,]
ivm_no <- pc[idx_ivm_no,]

p_val_ivm_pc <- perm_mul_t_test(t1 = ivm_yes, t2 = ivm_no)  # 0.104


#### GRADING ####
idx_grad_yes <- which(main$GRADING==3)   # rows of the patients with high grading 
idx_grad_no <- which(main$GRADING==1 | main$GRADING==2)

grad_yes <- pc[idx_grad_yes,]
grad_no <- pc[idx_grad_no,]

p_val_grad_pc <- perm_mul_t_test(t1 = grad_yes, t2 = grad_no)  # 0.393


#### STATO V/M ####
idx_vm_yes <- which(main$STATO.VIVO.MORTO==1)   # dead patients
idx_vm_no <- which(main$STATO.VIVO.MORTO==0)

vm_yes <- pc[idx_vm_yes,]
vm_no <- pc[idx_vm_no,]

p_val_vm_pc <- perm_mul_t_test(t1 = vm_yes, t2 = vm_no)  # 0.238


#### RECIDIVA ####
idx_rec_yes <- which(main$RECIDIVA==1)  
idx_rec_no <- which(main$RECIDIVA==0)

rec_yes <- pc[idx_rec_yes,]
rec_no <- pc[idx_rec_no,]

p_val_rec_pc <- perm_mul_t_test(t1 = rec_yes, t2 = rec_no)  # 0.256


# Now on Portal Margin
# We remove an unit on main which is not in pm
main <- main[-256,]

#### IVM ####
idx_ivm_yes <- which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1)   # rows of the patients with IVM
idx_ivm_no <- which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0)

ivm_yes <- pm[idx_ivm_yes,]
ivm_no <- pm[idx_ivm_no,]

p_val_ivm_pm <- perm_mul_t_test(t1 = ivm_yes, t2 = ivm_no)  # 0.674


#### GRADING ####
idx_grad_yes <- which(main$GRADING==3)   # rows of the patients with high grading 
idx_grad_no <- which(main$GRADING==1 | main$GRADING==2)

grad_yes <- pm[idx_grad_yes,]
grad_no <- pm[idx_grad_no,]

p_val_grad_pm <- perm_mul_t_test(t1 = grad_yes, t2 = grad_no)  # 0.624


#### STATO V/M ####
idx_vm_yes <- which(main[1:260,]$STATO.VIVO.MORTO==1)   # dead patients
idx_vm_no <- which(main[1:260,]$STATO.VIVO.MORTO==0)

vm_yes <- pm[idx_vm_yes,]
vm_no <- pm[idx_vm_no,]

p_val_vm_pm <- perm_mul_t_test(t1 = vm_yes, t2 = vm_no)  # 0.309


#### RECIDIVA ####
idx_rec_yes <- which(main[1:260,]$RECIDIVA==1)  
idx_rec_no <- which(main[1:260,]$RECIDIVA==0)

rec_yes <- pm[idx_rec_yes,]
rec_no <- pm[idx_rec_no,]

p_val_rec_pm <- perm_mul_t_test(t1 = rec_yes, t2 = rec_no)  # 0.303



