# We define here some functions to make permutational tests (taken from the scripts of the lab sessions)

# Function for univariate permutational test on the median of two independent samples
# INPUT: x, y, samples from the two univariate distributions
# OUTPUT: p-value of the permutational test H0:med(x)=med(y)

perm_median_test=function(x,y,iter=1e4){
  
  
  T0=abs(median(x)-median(y))  
  T_stat=numeric(iter)
  x_pooled=c(x,y)
  n=length(x_pooled)
  n1=length(x)
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(median(x1_perm) - median(x2_perm))
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  return(p_val)
  
}

perm_mean_test=function(x,y,iter=1e4){
  
  T0=abs(mean(x)-mean(y))  
  T_stat=numeric(iter)
  x_pooled=c(x,y)
  n=length(x_pooled)
  n1=length(x)
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  return(p_val)
  
}


## function for a two-sample independent test in the multivariate setting
## the two populations are those with positive/negative outcome

perm_mul_t_test <- function(t1,t2,B=1e3){
  # t1 and t2 are the two populations divided by outcome
  t1.median<- depthMedian(t1,depth_params = list(method='Tukey'))
  t2.median <- depthMedian(t2,depth_params = list(method='Tukey'))
  
  n1 <- dim(t1)[1]
  n2 <- dim(t2)[1]
  n  <- n1 + n2
  
  T20 <- max(abs(t1.median-t2.median))
  
  T2 = numeric(B)
  pb=progress_bar$new(total=1e3)  # create a progress object of dimension B
  pb$tick(0) 
  for(perm in 1:B){
    t_pooled = rbind(t1,t2)
    permutation = sample(n)
    t_perm = t_pooled[permutation,]    # only permute rows
    t1_perm = t_perm[1:n1,]
    t2_perm = t_perm[(n1+1):n,]
    
    t1.median_perm = depthMedian(t1_perm,depth_params = list(method='Tukey'))
    t2.median_perm = depthMedian(t2_perm,depth_params = list(method='Tukey'))
    T2[perm]  = max(abs(t1.median_perm-t2.median_perm))
    pb$tick()
  }
  p_val = sum(T2>=T20)/B
  
}
