barplot_func <- function(traits, df){
  
  # print off trait_line combination
  print(paste(traits))
  
  bp <- ggplot(data = df %>% filter(trait == traits),
               aes(x=line, y=gene_count, fill=region)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = region_cols) +
    theme_bw() +
    xlab("Experiment") +
    ylab("Genes (no.)") +
    ggtitle(paste0(traits)) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size=14), 
          axis.title.x = element_blank(), 
          legend.title = element_text(colour="black", size=10, face="bold"),
          legend.text = element_text(colour="black", size=8),
          strip.text = element_text(size = 12, face="bold"),
          axis.text.x = element_text(size=14),
          legend.position="none",
          legend.background = element_blank(),
          plot.title = element_text(hjust=0.5, size=18, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    )
  
  
  return(bp)
}