var_effect_corr_func <- function(combs, regx, regy, lines, cols, df){
  
  # print combinations
  print(paste0(combs)) # y vs x
  
  # format traits
  sig.regx <- paste0("sig.", regx,".",lines)
  beta.regx <- paste0("beta.", regx,".",lines)
  sig.regy <- paste0("sig.", regy,".",lines)
  beta.regy <- paste0("beta.", regy,".",lines)
  
  # all vars pleio vs not
  df$pleio <- ifelse(df[[sig.regx]] == "TRUE" &
                       df[[sig.regy]] == "TRUE", "Pleio","All")
  df$pleio <- factor(df$pleio, levels = c("All","Pleio"))
  
  # subset to significant SNPs (at least one trait)
  sig_snps <- df %>%
    filter(get(sig.regx) == "TRUE" | get(sig.regy) == "TRUE") %>%
    droplevels(.)
  
  ## significant SNPs, either trait
  sig_snps$traits <- ifelse(sig_snps[[sig.regx]] == "TRUE" &
                              sig_snps[[sig.regy]] == "TRUE", "Both_traits",
                            ifelse(sig_snps[[sig.regx]] == "TRUE" &
                                     sig_snps[[sig.regy]] == "FALSE", paste0(regx, " only"),
                                   paste0(regy, " only")))
  # order factor levels
  sig_snps$traits <- factor(sig_snps$traits, levels = c(paste0(regy, " only"), paste0(regx, " only"), "Both_traits"))
  # for regression lines:
  sig_snps$pleio <- ifelse(sig_snps$traits == "Both_traits","Pleio","All")
  # categorize effects
  sig_snps$effect <- ifelse(sig_snps[[beta.regx]] > 0 & 
                              sig_snps[[beta.regy]] > 0, "Top_right", 
                            ifelse(sig_snps[[beta.regx]] < 0 & 
                                     sig_snps[[beta.regy]] > 0, "Top_left", 
                                   ifelse(sig_snps[[beta.regx]] > 0 &
                                            sig_snps[[beta.regy]] < 0, "Bottom_right",
                                          "Bottom_left")))
  
  # order factor levels (clock-wise)
  sig_snps$effect <- factor(sig_snps$effect, 
                            levels = c("Top_right", "Bottom_right", "Bottom_left", "Top_left"))
  
  # add in effect direction (minor allele)
  sig_snps$trait1_dir <- ifelse(sig_snps[[beta.regx]] > 0, "+","-")
  sig_snps$trait2_dir <- ifelse(sig_snps[[beta.regy]] > 0, "+","-")
  
  # add in effect significance
  sig_snps$trait1_sig <- ifelse(sig_snps[[sig.regx]] == TRUE, "S","NS")
  sig_snps$trait2_sig <- ifelse(sig_snps[[sig.regy]] == TRUE, "S","NS")
  
  # add col for pw trait comb
  sig_snps$comb <- paste0(regx,"_vs_",regy,".",lines)
  sig_snps$comb_dir <- paste0(regx,sig_snps$trait1_dir,"_vs_",regy,sig_snps$trait2_dir)
  
  # summarize how many SNPs in each quadrant
  sig_snps.prop <- sig_snps %>%
    group_by(effect,traits, .drop = FALSE) %>%
    summarize(total_count = n()) %>%
    group_by(effect) %>%
    mutate(prop = total_count/sum(total_count)) %>%
    ungroup(.) %>%
    select(-total_count) %>%
    spread(key = traits, value = prop) %>%
    mutate(exp_pleio = get(paste0(regy, " only")) * get(paste0(regx, " only")))
  
  ## count pleio only
  sig_snps.sum.pleio <- sig_snps %>%
    filter(traits == "Both_traits") %>%
    group_by(effect, .drop = FALSE) %>%
    summarize(pleio_count = n())
  
  ## count total
  sig_snps.sum.tot <- sig_snps %>%
    group_by(effect, .drop = FALSE) %>%
    summarize(tot_count = n())
  
  ## create a list
  sig_snps.list <- list(sig_snps.prop, sig_snps.sum.pleio, sig_snps.sum.tot)
  
  ## combine lists
  sig_snps.sum <- sig_snps.list %>%
    reduce(full_join, by = "effect") %>%
    mutate(exp_count = round(tot_count * exp_pleio, 0),
           ann = paste0(pleio_count,"(",tot_count,")"))
  
  # annotation config for plot
  annotations <- data.frame(
    xpos = c(-Inf,-Inf,Inf,Inf),
    ypos =  c(-Inf, Inf,-Inf,Inf),
    annotateText = c(paste0(sig_snps.sum[3,9]), paste0(sig_snps.sum[4,9]),
                     paste0(sig_snps.sum[2,9]), paste0(sig_snps.sum[1,9])),
    hjustvar = c(0,0,1,1) ,
    vjustvar = c(0,1,0,1)) # <- adjust
  
  # use function to extract r2 and p-vals
  source("../Source_code/lmp_func.R")
  
  # conduct linear regression for all significant and pleio only vars
  lm_pleio <- lm(get(beta.regy) ~ get(beta.regx), data = sig_snps %>% filter(traits == "Both_traits"))
  lm_all <- lm(get(beta.regy) ~ get(beta.regx), data = df)
  
  # format the p and r2 vals
  pleio_p <- round(lmp_func(lm_pleio)[[1]][[1]],3)
  all_p <- round(lmp_func(lm_all)[[1]][[1]],3)
  pleio_r <- round(lmp_func(lm_pleio)[[2]][[1]],3)
  all_r <- round(lmp_func(lm_all)[[2]][[1]],3)
  
  All_stats <- bquote(italic(R)^2 == .(all_r))
  Pleio_stats <- bquote(italic(R)^2 == .(pleio_r))
  
  #All_stats <- bquote(italic(R)^2 == .(all_r)*";" ~ "p" == .(all_p))
  #Pleio_stats <- bquote(italic(R)^2 == .(pleio_r)*";" ~ "p" == .(pleio_p))
  
  # produce the plot
  palette <- "Royal1" ## from wesanderson package
  
  if (lines == "DZA"){
  plot_colors <- c("#4E820C","#A4D466","#466323")
  names(plot_colors) <- levels(sig_snps$traits)
  colScale <- scale_color_manual(values = plot_colors)
  }
  
  else if (lines == "A17"){
    plot_colors <- c("#AD208D","#DE6BC4","#69265A")
    names(plot_colors) <- levels(sig_snps$traits)
    colScale <- scale_color_manual(values = plot_colors)
  }
  
  ### universal variants to highlight
  universal_vars <- c('freebayes-snp-psyma-348072-A',
 ' freebayes-snp-psyma-395351-T',
 ' freebayes-snp-psyma-592463-T',
 ' freebayes-snp-psyma-606537-T',
 ' freebayes-snp-psymb-1642586-A',
  'freebayes-snp-psyma-359100-A',
  'freebayes-snp-psyma-359911-C',
  'freebayes-snp-psyma-360129-A',
  'freebayes-snp-psyma-360507-A',
  'freebayes-snp-psymb-812107-A',
  'freebayes-snp-psymb-812690-G',
  'freebayes-snp-psymb-812934-G',
  'freebayes-snp-psymb-812996-A',
  'freebayes-snp-psymb-813251-C',
  'freebayes-snp-psymb-813356-C',
  'freebayes-snp-psymb-813587-C',
  'freebayes-snp-psymb-813638-G',
  'freebayes-snp-psymb-813812-C')
  
  
  plot <- ggplot(data = sig_snps, aes(x=get(beta.regx), y=get(beta.regy))) +
    geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=Inf), fill="grey90", alpha=0.5) + 
    geom_rect(aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=0), fill="grey90", alpha=0.5) + 
    geom_point(data=sig_snps %>% filter(traits == paste0(regx, " only")),
               aes(color = traits, size = af), alpha = 0.5) + 
    geom_point(data=sig_snps %>% filter(traits == paste0(regy, " only")),
               aes(color = traits, size = af), alpha = 0.5) + 
    geom_smooth(data=df, method="lm", formula = "y ~ x", se=FALSE, color = I(cols),
                aes(linetype = pleio)) +   
    #geom_smooth(data=sig_snps %>% filter(traits == "Both_traits"),
    #            method="lm", formula = "y ~ x", se=FALSE, linetype = 2, color = I(cols)) +  
    geom_point(data=sig_snps %>% filter(traits == "Both_traits"),  
               aes(x=get(beta.regx), y=get(beta.regy), size = af, color = traits)) + 
    geom_point(data = sig_snps %>% filter(rs %in% universal_vars),
               color = "black", aes(size = af)) +
    geom_hline(aes(yintercept=0), linetype="dashed")+
    geom_vline(aes(xintercept=0), linetype="dashed")+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),
              color = "black", size = 5) +
    xlab(NULL) +
    ylab(NULL) +
    colScale +
    guides(size = "none", color = "none") +
    scale_linetype_discrete(breaks=c("All", "Pleio"),labels=c(All_stats,Pleio_stats)) +
    theme(axis.text.y = element_text(size=16), 
          axis.title.y = element_text(size=20, color = "#000000"), 
          axis.title.x = element_text(size=20, color = "#000000"),
          axis.text.x = element_text(size=16), legend.position=c(0.2,0.85), 
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.background = element_rect(fill=alpha('white', 0.4)),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.text = element_text(size = 12),
          plot.title = element_text(size=16, face = "bold", hjust = 0.5),
          panel.background = element_rect(fill=wes_palette(palette)[[3]]),
          strip.background = element_rect(fill=wes_palette(palette)[[4]]),
          strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.1,0.5,0.1,0, "cm")
    )
  
  return(list(plot))
  
}