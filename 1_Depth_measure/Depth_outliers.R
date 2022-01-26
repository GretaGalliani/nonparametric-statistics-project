library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(packagefinder)
library(aplpack)
library(robustbase)

#### IMPORT OF DATA ----

main <- read.table("Datasets/main.txt")
pc <-read.table("Datasets/pc_post_correlation.txt")
pm <-read.table("Datasets/pm_post_correlation.txt")

pc.label <- pc[,c(1,2)]
pc <- pc[,-c(1,2)]
pm.label <- pm[,c(1,2)]
pm <- pm[,-c(1,2)]

pc_pca <-read.table("Datasets/pc_post_pca.txt")
pm_pca <-read.table("Datasets/pm_post_pca.txt")

pc_pca <- pc_pca[,-c(1,2)]
pm_pca <- pm_pca[,-c(1,2)]

set.seed(1)

#### PORTAL CORE #### 

# Tukey distance
tuk_depth_core <- depth(pc,method='Tukey')
depthMedian(pc, depth_params = list(method='Tukey'))
# Mahalanobis distance
mal_depth_core <- depth(pc,method='Mahalanobis')


#### PORTAL CORE AFTER PCA #### 

# Tukey distance
tuk_depth_corepca <- depth(pc_pca[,c(1,2)],method='Tukey')
depthMedian(pc_pca[,c(1,2)], depth_params = list(method='Tukey'))

quartz()
depthContour(pc_pca[,c(1,2)],depth_params = list(method='Tukey'))

# Mahalanobis distance
mal_depth_corepca <- depth(pc_pca[,c(1,2)],method='Mahalanobis')

quartz()
depthContour(pc_pca[,c(1,2)],depth_params = list(method='Mahalanobis'))


#Bagplot
quartz()
bag_pc <- bagplot(pc_pca[,c(1,2)])

outlying_obs_pc <- bag_pc$pxy.outlier    #12 outliers
ind_outlying_obs_pc <- which(apply(pc_pca[,c(1,2)],1,function(x) all(x %in% outlying_obs_pc)))

outlier_data_pc <- cbind(main[ind_outlying_obs_pc,],pc[ind_outlying_obs_pc,])

# observing these outliers we can't find particular pattern in the clinical features
# but outliers already found by univariate procedures are contained 
# in this new subset

##--- Multivariate case

bagplot_matrix_pc <- aplpack::bagplot.pairs(pc_pca)
# this, unfortunately, doesn't allow to obtain outlying observations for every pair




#### PORTAL MARGIN AFTER PCA #### 
#Tukey distance
tuk_depth_marginpca <- depth(pm_pca[,c(1,2)],method='Tukey')
depthMedian(pm_pca[,c(1,2)], depth_params = list(method='Tukey'))

quartz()
depthContour(pm_pca[,c(1,2)],depth_params = list(method='Tukey'))

#Mahalanobis distance
mal_depth_marginpca <- depth(pm_pca[,c(1,2)],method='Mahalanobis')

quartz()
depthContour(pm_pca[,c(1,2)],depth_params = list(method='Mahalanobis'))

#Bagplot
quartz()
bag_pm <- bagplot(pm_pca[,c(1,2)])

outlying_obs_pm <- bag_pm$pxy.outlier
ind_outlying_obs_pm <- which(apply(pm_pca[,c(1,2)],1,function(x) all(x %in% outlying_obs_pm)))
#11 outliers 

outlier_data_pm <- cbind(main[ind_outlying_obs_pm,],pm[ind_outlying_obs_pm,])


##--- Multivariate case

bagplot_matrix_pm <- aplpack::bagplot.pairs(pm_pca)


### see which are the observations which are common outliers
out <- intersect(ind_outlying_obs_pm,ind_outlying_obs_pc)
outliers_intersection <- cbind(main[out,],pc[out,])


##################________________ CONCLUSIONS _______________#################
# expected outliers from Applied Stat: 204,235,241,122 + 102 for portal margin
# outliers from Portal Core: 77 106 122 129 135 180 181 195 204 235 241 247 
# outliers from Portal Margin: 108 111 118 122 145 149 183 192 204 235 241 
# common outliers: 122 204 235 241

# All outliers previously detected are also identified by depth measures

# Patterns in outliers identified by PC:
# all patients apart from 129 have IVM = 1, Single.nodule = 1 apart from 247, 
# CHEMIOTERAPIA.NEOADIUVANTE = 0 apart from 247

# Patterns in outliers identified by PM:
# all patients apart from 183 have IVM = 1, SINGLE.NODULE=1, 
# NUMERO..NO.SATELLITI =1


# next: see if all the outliers found by Depth measures are also found 
# with robust estimation techniques



