lin_reg_func <- function(combs, exp1, exp2, df){
  
  ## print experiments being regressed:
  print(paste0(combs))  
  
  ## models
  reg <- lm(get(exp2) ~ get(exp1), data = df)
  reg_coef <- as.data.frame(summary(reg)$coefficients)
  reg_coef$comp <- paste0(combs)
  r <- cor(df[[exp1]], df[[exp2]],  method = "pearson", use = "complete.obs")

  return(list(r, reg_coef))
}
