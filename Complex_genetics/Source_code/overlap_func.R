overlap_func <- function(traits, df){
  
    plast <- paste0(traits,".plast")
    
    # print traits
    print(traits)
        
    # Summarize at the gene-level (shoot and plasticity only) 
    SNPs_ann_gene.plast <- SNPs_sig_long %>%
      filter(trait.c %in% c(traits, plast) & ncbi_func != "non-coding") %>%
      droplevels(.) %>%
      group_by(region, gene_ID, RefSeq_ID, ncbi_func) %>%
      summarize(
        traits = paste(unique(trait.c), collapse = ", "),
        lines = paste(unique(line), collapse = ", ")) %>%
      as.data.frame(.)
  
    # add in relevant cols
    SNPs_ann_gene.plast$overlap <- ifelse(SNPs_ann_gene.plast$traits == paste0(traits,", ", plast) |
                                             SNPs_ann_gene.plast$traits == paste0(plast,", ", traits), "both_traits",
                                           ifelse(SNPs_ann_gene.plast$traits == paste0(traits), "trait_only",
                                                  "plast_only"))
    
    # save
    write.csv(SNPs_ann_gene.plast, file = paste0("./Data_output/Gene-lvl_summaries/plast_overlap_", traits,".csv"),
              row.names = FALSE)

  (SNPs_ann_gene.plast.sum <- SNPs_ann_gene.plast %>%
      group_by(overlap) %>%
      summarize(gene_count = n()))
  
  return(SNPs_ann_gene.plast.sum)
  
}
        