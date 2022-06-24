corr_func <- function(traits){
  
  ## print traits
  print(paste0(traits))
  

  ## pull df for each trait
  traits_comb <- comb_out[[traits]]
  
  ## make corr plots
  corr_plot <- ggplot(traits_comb, aes(x=abs(beta), y = -log10(p_lrt))) +
    geom_point(aes(size = af), alpha = 0.5) +
    geom_smooth(method = "lm") +
    facet_wrap(~chr) +
    ggtitle(paste0(traits)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  return(corr_plot)
  
}