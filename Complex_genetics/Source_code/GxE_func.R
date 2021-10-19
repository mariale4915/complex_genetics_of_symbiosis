GxE_func <- function(combs, traits, lines, df){
  
  ## print combination
  print(paste(combs))
  
  ## linear mixed models
  lmm <- lmer(sqrt(get(traits)) ~ strain_ID*exp + (1|rack), data = df %>%
                filter(line==lines))
  
  ## type 3 ANOVAs
  Anova_out <- Anova(lmm, type = 3)
  Anova_out$trait_line <- paste(combs)
  #Anova_out <- as.data.frame(Anova_out)
  lmm_rack <- lm(sqrt(get(traits)) ~ strain_ID*exp, data = df %>%
                   filter(line==lines))
  rack_sig <- anova(lmm, lmm_rack)
  rack_sig$trait_line <- paste(combs)
  #rack_sig <- as.data.frame(rack_sig)
  
  return(list(Anova_out, rack_sig))
}