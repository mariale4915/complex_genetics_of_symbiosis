corr.plot_func <- function(traits_regions, traits){
  
  ## print traits
  print(paste0(traits_regions))
  
  ## betas from ulmm and mlmm are highly correlated
  ### gather betas from UV and MV analyses
  comb_l1 <- comp.sigs_out[[traits_regions]] %>%
    select(-contains(paste0(traits))) %>% ## get rid of other betas for consistent col naming
    gather(analyses, MV_betas, beta_1:beta_4) %>%
    mutate(rs_betas = paste0(rs, "_", analyses))
  comb_l2 <- comp.sigs_out[[traits_regions]] %>%
    select(-contains("beta")) %>%
    rename(beta_1 = paste0(traits, "_DZA_1"), 
           beta_2 = paste0(traits, "_A17_2"), 
           beta_3 = paste0(traits, "_DZA_3"), 
           beta_4 = paste0(traits, "_A17_4")) %>%
    gather(analyses, UV_betas, beta_1:beta_4) %>%
    mutate(rs_betas = paste0(rs, "_", analyses))
  comb.l <- left_join(comb_l1, comb_l2[,c("rs_betas","UV_betas")], by = "rs_betas")
  ### corrplots
  corr_plot <- ggplot(comb.l, aes(x=UV_betas, y = MV_betas)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~analyses) +
    theme_bw()
  
  return(list(corr_plot))
  
}