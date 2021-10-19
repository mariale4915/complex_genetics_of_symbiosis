herit_func <- function(combs, traits, exps, df){
  
  # print combinations
  print(paste0(combs))
  
  # linear mixed models
  lmm <- lmer(sqrt(get(traits)) ~ (1|strain_ID) + (1|rack), data = df %>%
                filter(line_exp==exps))
  lmm_line <- lmer(sqrt(get(traits)) ~ (1|rack), data = df %>%
                     filter(line_exp==exps))
  
  # variance and significance
  varcor_out <- as.data.frame(VarCorr(lmm))
  aov_out <- anova(lmm, lmm_line)
  varcor_out$trait_exp <- paste(combs)
  aov_out$trait_exp <- paste(combs)
  
  return(list(varcor_out, aov_out))
}