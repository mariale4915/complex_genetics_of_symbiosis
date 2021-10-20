bar_func <- function(combs, traits, lines, df){
  
  ## print trait line comb:
  print(paste0(combs))
  
  ## designate plotting colors
  line_region_cols <- c("#DC9CCD","#CF77BA","#B32691",
                        "#B6CB9A","#99B771","#558713")
  names(line_region_cols) <- levels(df$line_region)
  colScale <- scale_fill_manual(values = line_region_cols)
  
  ## create plot
  bp <- ggplot(data = df %>% filter(trait == traits & line == lines) %>% droplevels(.),
               aes(x=exp, y=gene_count, fill=line_region)) +
    geom_bar(width = 1, stat = "identity") +
    theme_bw() +
    xlab("Experiment") +
    ylab("Genes (no.)") +
    colScale +
    scale_y_continuous(limits = c(0,550)) +
    # ggtitle(paste0(combs)) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size=14),
          legend.title = element_text(colour="black", size=10, face="bold"),
          legend.text = element_text(colour="black", size=8),
          strip.text = element_text(size = 12, face="bold"),
          legend.position="none",
          legend.background = element_blank(),
          plot.title = element_text(hjust=0.5, size=18, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    )
  
  return(list(bp))
}