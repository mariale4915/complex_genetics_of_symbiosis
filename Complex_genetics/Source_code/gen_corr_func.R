gen_corr_func <- function(traits1, traits2, cols, df){
  
  # print off regression combinations
  print(paste0(traits1,"_vs_",traits2))
  
  # linear regression model
  reg <- lm(get(traits2) ~ get(traits1), data = df)
  reg_coef <- as.data.frame(summary(reg)$coefficients)
  reg_coef$comp <- paste(traits1, "Vs", traits2, sep="_")
  
  # format significance value
  reg_coef$sig <- ifelse(reg_coef$'Pr(>|t|)' < 0.001, "***",
                         ifelse(reg_coef$'Pr(>|t|)' < 0.01, "**",
                                ifelse(reg_coef$'Pr(>|t|)' < 0.05, "*",
                                       ifelse(reg_coef$'Pr(>|t|)' < 0.1, ".", " "))))
  
  # calculate pearson r correlation coefficient
  r <- cor(df[[traits1]], df[[traits2]],  method = "pearson", use = "complete.obs")
  rr <- round(r, 3)
  
  # concatenate pearson r and sig
  r_sig <- paste0(rr, reg_coef[2,6])
  
  # plotting function
  fig <- ggplot(data=df, aes(x=get(traits1), y = get(traits2))) +
    geom_smooth(color = I(cols), method = "lm", formula = "y ~ x") +   
    geom_point(size = 3, color = I(cols)) +   
    annotate("text", x=Inf, y = -Inf, vjust=0, hjust=1, size = 6,
             label = r_sig) +
    xlab(NULL)  +
    ylab(NULL)  +  
    theme_bw() +
    theme(axis.title.y = element_text(size=24), 
          axis.text.y = element_text(size=20), 
          axis.title.x = element_text(size=24), 
          legend.title = element_text(colour="black", size=10, face="bold"),
          legend.text = element_text(colour="black", size=8),
          strip.text = element_text(size = 12, face="bold"),
          axis.text.x = element_text(size=20),
          legend.position="right",
          legend.background = element_blank(),
          plot.title = element_text(hjust=0.5, size=24),
          panel.background = element_rect(fill="transparent"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    )
  
  return(list(reg_coef, fig))
}