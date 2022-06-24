manhattan.plot_func <- function(traits_regions, traits){
  
  ## print traits
  print(paste0(traits_regions))
  
  ## load df
  comb <- comp.sigs_out[[traits_regions]]
  L <- 7965
  comb$w <- 1/L

  ## plotting
  comb$p.adj.FDR <- p.adjust(comb$p_lrt, method = "fdr")
  comb$p.adj <- comb$p_lrt/comb$w
  
  ## manhattan_plots 
  manhattan_plot <- ggplot(comb, aes(x=ps/1e6, y = -log10(p.adj))) + 
    geom_point(alpha = 0.5, color = "black") +
    geom_point(data = comb %>% filter(!is.na(RefSeq_ID), cat == "Universal"),
               color = "#9F7871", size = 4) +
    geom_point(data = comb %>% filter(!is.na(RefSeq_ID), cat == "GxG", host == "DZA only"),
              color = "#79AB39", size = 3) +
    geom_point(data = comb %>% filter(!is.na(RefSeq_ID), cat == "GxG", host == "A17 only"),
               color = "#C646A9", size = 3) +
    theme_bw() +
    ggtitle(paste0(traits_regions)) +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = 2) + ## alpha
    xlab("Genomic position (Mbp)") +
    ylab(expression(paste(-log[10], " p-value"))) +
    theme(axis.title.y = element_text(colour = "black", size = 18), 
          axis.text.y = element_text(size=14), 
          axis.title.x = element_text(colour = "black", size = 18), 
          legend.title = element_text(colour="black", size=10, face="bold"),
          legend.text = element_text(colour="black", size=8),
          strip.text = element_text(size = 12, face="bold"),
          axis.text.x = element_text(size=14),
          legend.position="right",
          legend.background = element_blank(),
          plot.title = element_text(hjust=0.5, size=18, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())

  return(list(manhattan_plot))
  
}