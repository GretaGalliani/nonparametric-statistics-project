library(rgl)
library(MASS)
library(DepthProc)

####--- Differences in distributions among outcomes using DDplots ---####

# importing datasets
main <- read.table("Datasets/main.txt")

pc <-read.table("Datasets/pc_post_correlation.txt")
pm <-read.table("Datasets/pm_post_correlation.txt")
pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)]
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)]

# I start working only on Portal Core, and I do not consider the PCA dataset 

#### IVM ####
idx_ivm_yes <- which(main$INVASIONE.VASCOLARE.MICROSCOPICA==1)   # rows of the patients with IVM
idx_ivm_no <- which(main$INVASIONE.VASCOLARE.MICROSCOPICA==0)

ivm_yes <- pc[idx_ivm_yes,]
ivm_no <- pc[idx_ivm_no,]

ddPlot(x = ivm_yes,y = ivm_no,depth_params = list(method='Tukey'), title = "Vascular Invasion or not")


#### GRADING ####
idx_grad_yes <- which(main$GRADING==3)   # rows of the patients with high grading 
idx_grad_no <- which(main$GRADING==1 | main$GRADING==2)

grad_yes <- pc[idx_grad_yes,]
grad_no <- pc[idx_grad_no,]

ddPlot(x = grad_yes,y = grad_no,depth_params = list(method='Tukey'), title = "High VS Low Grading")


#### STATO V/M ####
idx_vm_yes <- which(main$STATO.VIVO.MORTO==1)   # dead patients
idx_vm_no <- which(main$STATO.VIVO.MORTO==0)

vm_yes <- pc[idx_vm_yes,]
vm_no <- pc[idx_vm_no,]

ddPlot(x = vm_yes,y = vm_no,depth_params = list(method='Tukey'), title = "Dead VS Alive")


#### RECIDIVA ####
idx_rec_yes <- which(main$RECIDIVA==1)  
idx_rec_no <- which(main$RECIDIVA==0)

rec_yes <- pc[idx_rec_yes,]
rec_no <- pc[idx_rec_no,]

ddPlot(x = rec_yes,y = rec_no,depth_params = list(method='Tukey'), title = "Relapse VS Non-relapse")


# Now on Portal Margin

#### IVM ####
idx_ivm_yes <- which(main[1:260,]$INVASIONE.VASCOLARE.MICROSCOPICA==1)   # rows of the patients with IVM
idx_ivm_no <- which(main[1:260,]$INVASIONE.VASCOLARE.MICROSCOPICA==0)

ivm_yes <- pm[idx_ivm_yes,]
ivm_no <- pm[idx_ivm_no,]

ddPlot(x = ivm_yes,y = ivm_no,depth_params = list(method='Tukey'), title = "Vascular Invasion or not")


#### GRADING ####
idx_grad_yes <- which(main[1:260,]$GRADING==3)   # rows of the patients with high grading 
idx_grad_no <- which(main[1:260,]$GRADING==1 | main[1:260,]$GRADING==2)

grad_yes <- pm[idx_grad_yes,]
grad_no <- pm[idx_grad_no,]

ddPlot(x = grad_yes,y = grad_no,depth_params = list(method='Tukey'), title = "High VS Low Gradingt")


#### STATO V/M ####
idx_vm_yes <- which(main[1:260,]$STATO.VIVO.MORTO==1)   # dead patients
idx_vm_no <- which(main[1:260,]$STATO.VIVO.MORTO==0)

vm_yes <- pm[idx_vm_yes,]
vm_no <- pm[idx_vm_no,]

ddPlot(x = vm_yes,y = vm_no,depth_params = list(method='Tukey'), title = "Dead VS Alive")


#### RECIDIVA ####
idx_rec_yes <- which(main[1:260,]$RECIDIVA==1)  
idx_rec_no <- which(main[1:260,]$RECIDIVA==0)

rec_yes <- pm[idx_rec_yes,]
rec_no <- pm[idx_rec_no,]

ddPlot(x = rec_yes,y = rec_no,depth_params = list(method='Tukey'), title = "Relapse VS Non-relapse")




####### DDPLOTS on outliers individuated by depth measures #######
out_pc <- c(77, 106 ,122 ,129 ,135, 180, 181, 195, 204, 235, 241, 247)
pc_good <- pc[-out_pc,]
pc_bad <- pc[out_pc,]

ddPlot(pc_good,pc_bad,depth_params = list(method='Tukey'))


out_pm <- c(108, 111, 118, 122, 145, 149, 183, 192, 204, 235, 241)
pm_good <- pm[-out_pm,]
pm_bad <- pm[out_pm,]

ddPlot(pm_good,pm_bad,depth_params = list(method='Tukey'))



