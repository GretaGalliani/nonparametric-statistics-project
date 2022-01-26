
our_loocv = function(df, form,y){
  df_split <- loo_cv(df)
  p <- dim(df)[1]
  pred1 <- numeric(p)
  obs1 <- numeric(p)
  pb=progress_bar$new(total=p)
  pb$tick(0)
  for (i in (1:p))
  {
    fit <- mgcv::gam(formula = form,
                     data = analysis(df_split[[1]][[i]]),
                     family = binomial,
                     select = T)
    
    pred1[i] <- predict(fit,assessment(df_split[[1]][[i]]),type='response')
    obs1[i] <- assessment(df_split[[1]][[i]])[y][[1]]
    pb$tick()
  }
  result <- list("predicted"=pred1,"observed"=obs1)
  return(result)
}