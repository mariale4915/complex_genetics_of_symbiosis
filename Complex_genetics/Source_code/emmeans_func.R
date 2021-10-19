## function to calculate emmeans:
emmeans_function <- function(combs, traits, exps, df){
  
  # print comb
  print(paste0(combs))
  
  # linear mixed model
  lmm <- lmer(sqrt(get(traits)) ~ strain_ID + (1|rack), data = df %>% filter(line_exp==exps)) 
  
  # calculate emmeans
  emms <- emmeans(lmm, ~ strain_ID)
  emms_sum <- summary(emms, infer= c(TRUE,TRUE), adjust = "bon")
  emms_sel <- emms_sum[,c("strain_ID","emmean","SE")]
  emms_sel$bt_emm <- (emms_sel$emmean)^2
  emms_sel$bt_SE <- (emms_sel$SE)^2
  
  return(list(emms_sel))
}