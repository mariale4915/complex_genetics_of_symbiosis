tally_func <- function(groups, df){
  
  print(paste0(groups))
  
  ## summarize pleio for rhiz
  df.sum <- df %>%
    filter(pleio == "Y" & group == groups) %>%
    droplevels(.) %>%
    group_by(line, cat, group, trait_comb, pleio_cat) %>%
    summarize(total_vars = n_distinct(rs), mean_af = mean(af_191))
  
  # New facet label names for trait comb variable
  trait.labs <- c("Chlorophyll vs. \n Nodules", "Shoot biomass vs. \n Nodules", 
                  "Chlorophyll vs. \n Nodule weight", "Shoot biomass vs. \n Nodule weight",
                  "Chlorophyll vs. \n Rhizo. relative fitness", "Shoot biomass vs. \n Rhizo. relative fitness",
                  "Nodule weight vs. \n Nodules","Rhizo. relative fitness vs. \n Nodules",
                  "Rhizo. relative fitness vs. \n Nodule weight")
  names(trait.labs) <- c("chloro_vs_nod", "shoot_vs_nod", 
                         "chloro_vs_nod.weight", "shoot_vs_nod.weight",  "chloro_vs_fit",
                         "shoot_vs_fit", "nod.weight_vs_nod", "fit_vs_nod","fit_vs_nod.weight")
  
  ## need to order trait_comb
  df.sum$trait_comb <- factor(df.sum$trait_comb, 
                                         levels=c("nod.weight_vs_nod","fit_vs_nod", "fit_vs_nod.weight",
                                                  "shoot_vs_nod","shoot_vs_nod.weight", "shoot_vs_fit",
                                                  "chloro_vs_nod","chloro_vs_nod.weight", "chloro_vs_fit"))
  
  ## factorize and order:
  df.sum$cat <- factor(df.sum$cat, levels = c("OD","D"))
  
  ## count number of discordant vs concordant vars
  df.sum_counts <- df.sum %>% group_by(cat) %>%
    summarise(count = sum(total_vars))
  
  (plot <- ggplot(df.sum, aes(x=line, y = total_vars, fill = cat)) +
      geom_bar(stat = "identity") +
      facet_wrap(~trait_comb, labeller = labeller(trait_comb = trait.labs)) +
      theme_bw() +
      xlab("Plant line") +
      ylab("Proportion of SNPs") +
      scale_fill_manual(name="Fitness effect",
                             breaks=c("OD", "D"),
                             labels=c(paste0("Discordant = ", df.sum_counts[1,2]),
                                 paste0("Concordant = ", df.sum_counts[2,2])),
                             values = c("#899da4","#dc863b")) +
      # geom_text(aes(label = round((1-mean_af), 3)),
      #           position=position_stack(vjust=0.5), colour="black") +
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_text(size=16), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size=16), 
            legend.title = element_blank(),
            legend.text = element_text(colour="black", size=12),
            strip.text = element_text(size = 14, face="bold"),
            legend.position=c(0.15,0.8),
            legend.background = element_blank(),
            plot.title = element_text(hjust=0.5, size=20, face = "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = margin(0,0,0,0, "pt")
      )
  )
  
  return(plot)
}