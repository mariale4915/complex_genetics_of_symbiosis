linkage_func <- function(vals){
  
  ## print LD vals
  print(paste0(vals))
  
  LG <- read.delim(paste0("./Linkage_analyses/LG", vals, ".tsv"))
  
  ## determine number and range of vars per group
  LG_sum <- LG %>%
    group_by(group) %>%
    summarize(no_vars = n(),
              min_pos = min(pos), 
              max_pos = max(pos), 
              range = ((max_pos - min_pos) + 1)/1000,
              chroms = paste(unique(chrom), collapse = ", "))
  
  ## custom binning
  LG_sum$vars_bin <- as.factor(ifelse(LG_sum$no_vars > 15, ">15", LG_sum$no_vars))
  
  LG_sum$vars_bin <- factor(LG_sum$vars_bin, levels = c("1","2","3","4","5","6","7","8","9",
                                                        "10","11","12","13","14","15",">15"))
  
  ## group by no of variants per group
  LG_sum2 <- LG_sum %>%
    group_by(vars_bin) %>%
    summarize(no_vars_per_group = n())
  
  ## plot
  (LD_groups_plot <- ggplot(LG_sum2, aes(x=vars_bin, y = no_vars_per_group)) +
      geom_bar(stat="identity") +
      theme_bw() +
      xlab("Variants per LD group (no.)") + 
      ylab("LD groups (no.)") +
      geom_text(aes(label = no_vars_per_group), nudge_y = 100) +
      theme(axis.title.y = element_text(colour = "black", size = 18), 
            axis.text.y = element_text(size=16), 
            axis.title.x = element_text(size=18), 
            axis.text.x = element_text(size=16), 
            legend.position="none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()))
  
  # rename genomic region
  LG_sum$region <- ifelse(LG_sum$chroms == "psyma" |
                            LG_sum$chroms == "psymb", "Megaplasmids", 
                          ifelse(LG_sum$chroms == "chr", "Chromosome",
                                 "inter_rep"))
  
  ## create facets based on LD span
  LG_sum$facet <- ifelse(LG_sum$range < 10, "<10", 
                         ifelse(LG_sum$range >= 10 & LG_sum$range <= 500, "10-500", ">500"))
  
  LG_sum$facet <- factor(LG_sum$facet, levels = c("<10","10-500",">500"))
  LG_sum$region <- factor(LG_sum$region, levels = c("Megaplasmids","Chromosome"))
  
  LG_sum_intra <- filter(LG_sum, region != "inter_rep")
  
  ## plot with facetting
  (LD_groups_span <- ggplot(LG_sum_intra, aes(x=range, fill = region)) +
      geom_histogram(data = LG_sum_intra %>% filter(facet == "<10"),
                     binwidth = 1) +
      geom_histogram(data = LG_sum_intra %>% filter(facet == "10-500"), 
                     binwidth = 5) +
      geom_histogram(data = LG_sum_intra %>% filter(facet == ">500"), 
                     binwidth = 100) +
      facet_wrap(~facet, scales = "free") +
      theme_bw() +
      xlab("Genomic span (bp)") + 
      ylab("LD groups (no.)") +
      scale_fill_discrete(name = "Genomic region") +
      theme(axis.title.y = element_text(colour = "black", size = 18), 
            axis.text.y = element_text(size=16), 
            axis.title.x = element_text(size=18), 
            axis.text.x = element_text(size=16), 
            legend.position=c(0.2,0.8),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()))
  
  comb.p <- plot_grid(LD_groups_plot, LD_groups_span,
                      ncol = 1, nrow = 2, align="v", labels=c("A","B"))
  
  save_plot(paste0("./Linkage_analyses/LG", vals, ".png"), comb.p,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_width = 12, base_height = 6)
  
  return(LG)
  
}