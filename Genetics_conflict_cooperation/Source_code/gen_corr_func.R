gen_corr_func <- function(combs, traits1, traits2, lines, cols, df){
  
  # print off regression combinations
  print(paste0(combs))
  
  # format traits
  trait1 <- paste0(traits1,".", lines)
  trait2 <- paste0(traits2,".", lines)
  
  # linear regression model
  reg <- lm(get(trait1) ~ get(trait2), data = df)
  reg_coef <- as.data.frame(summary(reg)$coefficients)
  reg_coef$comp <- paste(trait1, "Vs", trait2, sep="_")
  
  # format significance value
  reg_coef$sig <- ifelse(reg_coef$'Pr(>|t|)' < 0.001, "***",
                         ifelse(reg_coef$'Pr(>|t|)' < 0.01, "**",
                                ifelse(reg_coef$'Pr(>|t|)' < 0.05, "*",
                                       ifelse(reg_coef$'Pr(>|t|)' < 0.1, ".", " "))))
  
  # calculate pearson r correlation coefficient
  r <- cor(df[[trait1]], df[[trait2]],  method = "pearson", use = "complete.obs")
  rr <- round(r, 3)
  
  # concatenate pearson r and sig
  r_sig <- paste0(rr, reg_coef[2,6])
  
  # plotting function
  fig <- ggplot(data=df, aes(x=get(trait2), y = get(trait1))) +
    geom_smooth(color = I(cols), method = "lm", formula = "y ~ x") +   
    geom_point(size = 2, color = I(cols)) +   
    annotate("text", x=Inf, y = Inf, vjust=-0.3, hjust=1, size = 5,
             label = r_sig) +
    coord_cartesian(clip = 'off') +
    xlab(NULL)  +
    ylab(NULL)  +  
    theme_bw() +
    theme(axis.title.y = element_text(size=16), 
          axis.text.y = element_text(size=14), 
          axis.title.x = element_text(size=16), 
          axis.text.x = element_text(size=14),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size=8),
          legend.position="none",
          legend.background = element_blank(),
          strip.text = element_text(size = 12, face="bold"),
          plot.title = element_text(hjust=0.5, size=18, face = "bold"),
          plot.margin = margin(15,5.5,5.5,5.5, "pt"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    )
  
  return(list(fig))
}
