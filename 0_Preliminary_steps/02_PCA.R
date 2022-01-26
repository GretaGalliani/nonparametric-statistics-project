
#### IMPORT DATA ----
main <- read.table("main.txt")
pc <-read.table("pc_post_correlation.txt")
pm <-read.table("pm_post_correlation.txt")

pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)]
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)]

#### PCA ON PORTAL CORE ----

#standardize the radiomics 
pc.sd = scale(pc)

# pca 
pca.pc <- princomp(pc.sd, scores=T)
summary(pca.pc)

# scores
scores.pc <- pca.pc$scores
scores.pc

# explained variance
x11()
par(mfrow=c(1,2))

plot(pca.pc, las=2,main="PCA on portal core",ylim=c(0,8))

plot(cumsum(pca.pc$sde^2)/sum(pca.pc$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(pc.sd),labels=1:ncol(pc.sd),las=2)

# loadings
load.pc <- pca.pc$loadings
load.pc

x11()
par(mar = c(2,2,2,1), mfrow=c(2,1))
for(i in 1:2)barplot(load.pc[,i], ylim = c(-1, 1), main=paste('Loadings PC ',i,sep=''))


#### PCA ON PORTAL MARGIN ----
#standardize the radiomics 
pm.sd = scale(pm)

# pca 
pca.pm <- princomp(pm.sd, scores=T)
summary(pca.pm)

# scores
scores.pm <- pca.pm$scores
scores.pm

# explained variance
x11()
par(mfrow=c(1,2))

plot(pca.pm, las=2,main="PCA on portal margin",ylim=c(0,8))

plot(cumsum(pca.pm$sde^2)/sum(pca.pm$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(pm.sd),labels=1:ncol(pm.sd),las=2)

# loadings
load.pm <- pca.pm$loadings
load.pm

x11()
par(mar = c(2,2,2,1), mfrow=c(2,1))
for(i in 1:2)barplot(load.pm[,i], ylim = c(-1, 1), main=paste('Loadings PM',i,sep=''))


#### SAVING THE PCA DATASETS

pc_post_pca <- cbind(pc.label, scores.pc[,1:6])
write.table(pc_post_pca, "pc_post_pca.txt")

for (i in 1:length(colnames(scores.pm))){
  colnames(scores.pm)[i] <- paste(colnames(scores.pm)[i],"margin",sep="_")
}
pm_post_pca <- cbind(pm.label, scores.pm[,1:6])
write.table(pm_post_pca, "pm_post_pca.txt")
