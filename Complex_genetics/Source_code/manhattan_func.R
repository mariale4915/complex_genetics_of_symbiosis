manhattan_func <- function(traits){
  
  ## print traits
  print(paste0(traits))
                          
  ## use harmonic p package
  gwas <- comb_out[[traits]]
  L <- 7965
  gwas$w <- 1/L
  R <- 1:nrow(gwas)
  comb_p <- p.hmp(gwas$p_lrt[R],gwas$w[R],L)
  comb_p.r <- round(comb_p, 3)
  
  # Define overlapping sliding windows of 50 megabase at 10 megabase intervals
  region.list <- c("chr","psyma","psymb")
  # Calculate the combined p-values for each window
  p.region <- sapply(region.list, function(region) {
  R <- which(gwas$chr==region)
  p.hmp(gwas$p_lrt[R],gwas$w[R],L)
  })
  
  # Calculate the sums of weights for each combined test
  w.region <- sapply(region.list,function(region) {
  R <- which(gwas$chr==region)
  sum(gwas$w[R])
  })
  
  # calculate adjusted p-value for each region
  p.region.adj <- p.region/w.region
  
  ## plotting
  gwas$p.adj <- gwas$p_lrt/gwas$w
  gwas$p.adj.FDR <- p.adjust(gwas$p_lrt, method = "fdr")
  
  manhattan_plot <-  ggplot(gwas, aes(x=ps/(1e6), y = -log10(p.adj))) + 
  facet_wrap(~chr, scales = "free_x") +
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 2) + ## alpha
  geom_hline(aes(yintercept = -log10(qharmonicmeanp(0.05, L))),
             linetype = 2, color = "grey") + ## HMP threshold
  geom_hline(data = gwas %>% filter(chr == "chr"),
             aes(yintercept = -log10(p.region.adj[1])), 
             linetype = 1, size = 2, color = "grey", alpha = 0.5) + ## adjusted p value for chr
  geom_hline(data = gwas %>% filter(chr == "psyma"),
             aes(yintercept = -log10(p.region.adj[2])), 
             linetype = 1, size = 2, color = "grey", alpha = 0.5) + ## adjusted p value for psyma
  geom_hline(data = gwas %>% filter(chr == "psymb"),
             aes(yintercept = -log10(p.region.adj[3])), 
             linetype = 1, size = 2, color = "grey", alpha = 0.5) + ## adjusted p value for psymb
  geom_point(alpha = 0.5) +
  ggtitle(paste0("Genome-Wide combined P = ", comb_p.r)) +
  theme_bw() + 
  xlab("Genomic position (Mbp)") +
  ylab(expression(paste(-log[10], " p-value"))) +
  theme(axis.title.y = element_text(colour = "black", size = 18), 
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14), 
        axis.title.x = element_text(colour = "black", size = 18), 
        legend.title = element_text(colour="black", size=10, face="bold"),
        legend.text = element_text(colour="black", size=8),
        strip.text = element_text(size = 12, face="bold"),
        legend.position="right",
        legend.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=18, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
  
  return(manhattan_plot)
  
}