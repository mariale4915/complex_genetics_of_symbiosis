gene_sum_func <- function(traits, df){
  
  # print traits
  print(traits)
  
  # Summarize at the gene-level (fitness only)
  SNPs_ann_gene <- df %>%
    filter(trait == traits & ncbi_func != "non-coding") %>%
    droplevels(.) %>%
    group_by(region, gene_ID, RefSeq_ID, ncbi_func, assoc) %>%
    mutate(
      no_effect_assoc = n(),
      assoc_no = paste0(assoc,"[",no_effect_assoc,"]")) %>%
    ungroup(.) %>%
    group_by(region, gene_ID, RefSeq_ID, ncbi_func) %>%
    summarize(
      no_vars = n_distinct(rs), ## no. of unique vars
      no_traits = n_distinct(trait), ## no. of unique traits
      no_effects = n(), ## total no. of associations
      min_ps = min(ps), 
      max_ps = max(ps), 
      min_score = min(scores), 
      max_score = max(scores), 
      ave_abs_score = mean(abs(scores)),
      ave_af = mean(af),
      assocs = paste(unique(assoc_no), collapse = ", "),
      lines = paste(unique(line), collapse = ", "),
      exps = paste(unique(exp), collapse = ", ")) %>%
    as.data.frame(.)
  
  # add in relevant cols
  SNPs_ann_gene$host <- ifelse(SNPs_ann_gene$lines == "A17, DZA" | SNPs_ann_gene$lines == "DZA, A17", "Both",
                               ifelse(SNPs_ann_gene$lines == "A17", "A17 only", "DZA only"))
  
  if (traits == "chloro" || traits == "height" || traits == "leaf" || traits == "shoot"){
    
  SNPs_ann_gene <- SNPs_ann_gene %>%
    mutate(host_exps = paste(host, exps, sep = "_"))
  
  SNPs_ann_gene$cat <- ifelse(SNPs_ann_gene$host_exps == "Both_1, 2, 3" | 
                                      SNPs_ann_gene$host_exps == "Both_1, 2, 4" |
                                      SNPs_ann_gene$host_exps == "Both_1, 3, 4" |
                                      SNPs_ann_gene$host_exps == "Both_2, 3, 4", "Partially Universal",
                                    ifelse(SNPs_ann_gene$host_exps == "Both_1, 2, 3, 4", "Universal",
                                           ifelse(SNPs_ann_gene$host_exps == "DZA only_1, 3" | 
                                                    SNPs_ann_gene$host_exps == "A17 only_2, 4", "GxG",
                                                  ifelse(SNPs_ann_gene$host_exps == "DZA only_1" | 
                                                           SNPs_ann_gene$host_exps == "DZA only_3" |
                                                           SNPs_ann_gene$host_exps == "A17 only_2" |
                                                           SNPs_ann_gene$host_exps == "A17 only_4", "GxE",
                                                         "other"))))
  
  
  ## save the summary file
  save(SNPs_ann_gene, file = paste0("./Data_output/Gene-lvl_summaries/genes_", traits, ".Rdata"))
  
  # summarize number of genes in ea. category
  SNPs_ann_gene.sum <- SNPs_ann_gene %>%
    group_by(cat, host_exps) %>%
    summarize(count = n())
  
  SNPs_ann_gene.sum$trait <- paste0(traits)
  
  ## create summary file at gene-level
  SNPs_ann_gene.out <- SNPs_ann_gene %>%
    select(-lines, -host_exps)
  
  write.csv(SNPs_ann_gene.out, paste0("./Data_output/Gene-lvl_summaries/genes_", traits, ".csv"), row.names = FALSE)
  
  return(SNPs_ann_gene.sum)
  }

  else if (traits == "chloro.plast" || traits == "height.plast" || traits == "leaf.plast" || traits == "shoot.plast"){
    
    ## save the summary file
    save(SNPs_ann_gene, file = paste0("./Data_output/Gene-lvl_summaries/genes_", traits, ".Rdata"))
    
    ## summarize number of genes, plasticity
    (genes_plast.sum <- SNPs_ann_gene %>%
        group_by(host) %>%
        summarize(count = n()))
    genes_plast.sum$trait <- paste0(traits)
    
    ## create summary file at gene-level
    SNPs_ann_gene.out <- SNPs_ann_gene %>%
      select(-lines, -exps)
    
    write.csv(SNPs_ann_gene.out, paste0("./Data_output/Gene-lvl_summaries/genes_", traits, ".csv"), row.names = FALSE)
    
    return(genes_plast.sum)
  }
  
}