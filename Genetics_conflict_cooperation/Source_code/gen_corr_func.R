gen_corr_func2 <- function(combs, traits_y, traits_x, cols, df){
  
  # print off regression combinations
  print(paste0(combs))
  
  # linear regression model
  reg <- lm(get(traits_y) ~ get(traits_x), data = df)
  reg_coef <- as.data.frame(summary(reg)$coefficients)
  reg_coef$comp <- paste(traits_y, "Vs", traits_x, sep="_")
  
  # format significance value
  reg_coef$sig <- ifelse(reg_coef$'Pr(>|t|)' < 0.001, "***",
                         ifelse(reg_coef$'Pr(>|t|)' < 0.01, "**",
                                ifelse(reg_coef$'Pr(>|t|)' < 0.05, "*",
                                       ifelse(reg_coef$'Pr(>|t|)' < 0.1, ".", " "))))
  
  # calculate pearson r correlation coefficient
  r <- cor(df[[traits_y]], df[[traits_x]],  method = "pearson", use = "complete.obs")
  rr <- round(r, 3)
  
  # concatenate pearson r and sig
  r_sig <- paste0(rr, reg_coef[2,6])
  
  # plotting function
  fig <- ggplot(data=df, aes(x=get(traits_x), y = get(traits_y))) +
    geom_smooth(color = I(cols), method = "lm", formula = "y ~ x") +   
    geom_point(size = 2, color = I(cols)) +   
    annotate("text", x=Inf, y = Inf, vjust=1.2, hjust=1, size = 5,
             label = r_sig) +
    xlab(NULL)  +
    ylab(NULL)  +  
    theme_bw() +
    theme(axis.title.y = element_text(size=18, face = "bold"), 
          axis.text.y = element_text(size=16), 
          axis.title.x = element_text(size=18), 
          axis.text.x = element_text(size=16),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size=8),
          legend.position="none",
          legend.background = element_blank(),
          strip.text = element_text(size = 12, face="bold"),
          plot.title = element_text(hjust=0.5, size=18, face = "bold"),
          plot.margin = margin(0, 0, 0, 0, "cm"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    )
  
  return(list(fig))
}
