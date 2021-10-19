plast_func <- function(trait, df){
  
  plast_DZA <- paste0(trait, ".plast_DZA_13")
  plast_A17 <- paste0(trait, ".plast_A17_24")
  
  (plot <- ggplot(data=df, aes(x= get(plast_DZA), y = get(plast_A17))) +
      geom_smooth(method = "lm", se = FALSE) + 
      geom_text(aes(label=strain_ID)) +  
      theme_bw())
  
  (summary(lm1 <- lm(get(plast_A17) ~ get(plast_DZA), 
                     data=df))) ## NS
  
  ## rank order for plasticity
  (plot_DZA <- ggplot(df, aes(x = reorder(strain_ID, get(plast_DZA)), y = get(plast_DZA))) +
      geom_bar(stat = "identity", fill = "#79AB39") +
      xlab(NULL) + 
      ylab(NULL) +
      theme_bw() +
      theme(axis.title.y = element_text(colour = "black", size = 24), 
            axis.text.y = element_text(size=20), 
            axis.title.x = element_text(size=24), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            legend.position="none",
            legend.title = element_text(colour="black", size=16, face="bold"),
            legend.text = element_text(colour="black", size=12),
            strip.text.x = element_text(size=12, face = "bold"),
            plot.title = element_text(hjust=0.5, size=20, face = "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()))
  
  (plot_A17 <-ggplot(df, aes(y=get(plast_A17), x = reorder(strain_ID, get(plast_A17)))) +
      geom_bar(stat = "identity", fill = "#C646A9") +
      xlab(NULL) + 
      ylab(NULL) +
      theme_bw() +
      theme(axis.title.y = element_text(colour = "black", size = 24), 
            axis.text.y = element_text(size=20), 
            axis.title.x = element_text(size=24), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            legend.position="none",
            legend.title = element_text(colour="black", size=16, face="bold"),
            legend.text = element_text(colour="black", size=12),
            strip.text.x = element_text(size=12, face = "bold"),
            plot.title = element_text(hjust=0.5, size=20, face = "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()))
  
  return(list(lm1, plot_DZA, plot_A17))
  
}