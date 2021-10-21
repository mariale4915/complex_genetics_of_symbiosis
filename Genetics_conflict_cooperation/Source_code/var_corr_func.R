var_corr_func <- function(combs, cols_D, cols_OD, df){
  
  ## print combinations
  print(paste0(combs)) # x vs y
  
  ## subset dataframe to specific traits
  df_sub <- df %>%
    filter(trait_comb_line == combs) %>%
    droplevels(.)
  
  ## summarize all sig vars
  df_sum <- df_sub %>%
    mutate(total = n()) %>% 
    group_by(cat, line_cat) %>%
    summarize(total = unique(total), 
              count = n(),
              prop = paste0(round(((count/total)*100),1),"%")
    )
  
  ## summarize only pleio vars
  df_pleio <- df_sub %>%
    filter(pleio == "Y") %>%
    group_by(cat, line_cat, .drop = FALSE) %>%
    summarize(count = n())
  
  ## create annotation df
  annotations <- data.frame(
    xpos = c(-Inf,Inf), ## TL, TR
    ypos =  c(Inf,Inf),
    annotateText = c(paste0(df_pleio["count"][4,1],
                            " (",df_sum["prop"][2,1],")"), ## TL: prop, OD
                     paste0(df_pleio["count"][1,1],
                            " (",df_sum["prop"][1,1],")")), ## TR: prop,D
    hjustvar = c(0,1) ,
    vjustvar = c(-0.3,-0.3),
    colors = c("#899da4","#dc863b")
    ) #<- adjust
  
  ## designate plotting colors
  line_cols <- c("#4E820C","#AD208D")
  names(line_cols) <- levels(df$line)
  colScale <- scale_color_manual(values = line_cols)
  
  ## set palette
  palette <- "Royal1" ## from wesanderson package
  
  ## plotting function
  plot <- ggplot(data = df_sub,aes(x=betas_x, y=betas_y)) +
      geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=Inf), fill="grey90", alpha=0.5) +
      geom_rect(aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=0), fill="grey90", alpha=0.5) +
      geom_point(alpha = 0.5, aes(size = af_191, color = line)) +
      geom_point(data = df_sub %>% filter(pleio == "Y"), 
                 aes(size = af_191), color = "black") +
      geom_hline(aes(yintercept = 0), linetype = 2) +
      geom_vline(aes(xintercept = 0), linetype = 2) +
      geom_text(data=annotations,
              aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),
              size = 5, color = annotations$colors) +
      coord_cartesian(clip = 'off') +
      colScale +
      xlab(NULL) +
      ylab(NULL) +
      # ggtitle(paste0(combs)) +
      ## additional aesthetics:
      theme_bw() +
      theme(axis.title.y = element_text(size=16),
            axis.title.x = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.position="none",
            legend.direction = 'vertical',
            legend.title = element_blank(),
            #legend.background = element_rect(fill=alpha('white', 0.4)),
            legend.key = element_rect(colour = "transparent", fill = "transparent"),
            legend.text = element_text(size = 12),
            plot.title = element_text(size=18, face = "bold", hjust = 0.5),
            panel.background = element_rect(fill=wes_palette(palette)[[3]]),
            strip.background = element_rect(fill=wes_palette(palette)[[4]]),
            strip.text = element_text(size = 12, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(15,15,5.5,5.5, "pt")
      )
  
  return(list(plot))
  
}