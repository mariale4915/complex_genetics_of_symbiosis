gene_sum_func <- function(traits, df){
  
  # print off trait
  print(paste(traits))
  
  # Summarize at the gene-level (each trait individually)
  SNPs_ann_gene <- df %>%
    filter(trait == traits & ncbi_func != "non-coding") %>%
    droplevels(.) %>%
    mutate(effect_cat = paste(traits, line_exp, sep = '-')) %>%
    group_by(region, RefSeq_ID, ncbi_func) %>%
    summarize(
      no_vars = n_distinct(ps),
      min_ps = min(ps), 
      max_ps = max(ps), 
      no_effects = n(),
      ave_score = mean(scores),
      lines = paste(unique(line), collapse = ", "),
      exps = paste(unique(exp), collapse = ", ")) %>%
    as.data.frame(.)
  
  SNPs_ann_gene$host <- ifelse(SNPs_ann_gene$lines == "A17, DZA" | SNPs_ann_gene$lines == "DZA, A17", "Both",
                               ifelse(SNPs_ann_gene$lines == "A17", "A17 only", "DZA only"))
  
  # Venn diagram summary (to region and experiment)
  SNPs_ann_gene.sum <- SNPs_ann_gene %>%
    mutate(host_exps = paste(host, exps, sep = "_")) %>%
    group_by(host_exps) %>%
    summarize(count = n())
  
  return(SNPs_ann_gene.sum)
}