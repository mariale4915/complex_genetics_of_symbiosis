rxn_func <- function(trait, df){
  
  ## print traits
  print(trait)
  
  (p <- ggplot(data = df, aes(x=exp, y=get(trait), group=strain_ID)) + 
     geom_line(position = position_dodge(0.3), color = "grey") +
     geom_point(data = df %>% filter(line == "DZA"), 
                size=3, position = position_dodge(0.3), color = "#4F820D",alpha = 0.5) + 
     geom_point(data = df %>% filter(line == "A17"), 
                size=3, position = position_dodge(0.3), color = "#AD208E",alpha = 0.5) +   
     facet_wrap(~line, scales = "free_x") +
     theme_bw() +
     xlab(NULL) + 
     ylab(paste0(trait)) +
     theme(axis.title.y = element_text(colour = "black", size = 24), 
           axis.text.y = element_text(size=20), 
           axis.title.x = element_text(size=24), 
           axis.text.x = element_text(size=20), 
           legend.position="none",
           legend.title = element_text(colour="black", size=16, face="bold"),
           legend.text = element_text(colour="black", size=12),
           strip.text.x = element_text(size=12, face = "bold"),
           plot.title = element_blank(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank()
     )
  )
  return(list(p))  
}