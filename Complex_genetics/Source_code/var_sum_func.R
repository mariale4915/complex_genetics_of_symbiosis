var_sum_func <- function(traits, df){
  
  # print traits
  print(traits)
  
  if (traits == "all") {
    # Summarize at the ps-level
    SNPs_ann_ps <- df %>%
      group_by(region, rs, gene_name, RefSeq_ID, ncbi_func, start_pos, end_pos, effect) %>% ## summarized to rs-level
      summarize(
        no_effects = n(),
        ave_score = mean(scores),
        maf = mean(af),
        assocs = paste(unique(assoc), collapse = ", "),
        lines = paste(unique(line), collapse = ", ")) %>%
      as.data.frame(.)
    
    SNPs_ann_ps$host <- ifelse(SNPs_ann_ps$lines == "A17, DZA" | SNPs_ann_ps$lines == "DZA, A17", "Both",
                               ifelse(SNPs_ann_ps$lines == "A17", "A17 only", "DZA only"))
    
    SNPs_ann_ps <- SNPs_ann_ps %>%
      select(-lines)
    
    ## save
    write.csv(SNPs_ann_ps, paste0("./Data_output/Var-lvl_summaries/SNPs_ann_ps_", traits, ".csv"), row.names = FALSE)
    return(SNPs_ann_ps)
  }
  
  else if (traits == "chloro1" || traits == "height1" || traits == "leaf1" || traits == "shoot" ||
      traits == "chloro.plast" || traits == "height.plast" || traits == "leaf.plast" || traits == "shoot.plast"){
  
  # Summarize at the ps-level
  SNPs_ann_ps <- df %>%
    filter(trait == traits) %>%
    droplevels(.) %>%
    group_by(region, rs, gene_name, RefSeq_ID, ncbi_func, start_pos, end_pos, effect) %>% ## summarized to rs-level
    summarize(
      no_effects = n(),
      ave_score = mean(scores),
      maf = mean(af),
      assocs = paste(unique(assoc), collapse = ", "),
      lines = paste(unique(line), collapse = ", ")) %>%
    as.data.frame(.)
  
  ## designate host (DZA only, A17 only, Both)
  SNPs_ann_ps$host <- ifelse(SNPs_ann_ps$lines == "A17, DZA" | SNPs_ann_ps$lines == "DZA, A17", "Both",
                             ifelse(SNPs_ann_ps$lines == "A17", "A17 only", "DZA only"))
  
  ## get rid of lines col
  SNPs_ann_ps <- SNPs_ann_ps %>%
    select(-lines)
  
  ## add in direction of associations (same sign = con, different signs = dis)
  SNPs_ann_ps$dir <- ifelse(SNPs_ann_ps$no_effects > 1 &
                              grepl(paste0(traits,"[-]"), 
                                    SNPs_ann_ps$assocs, fixed = FALSE) &
                              grepl(paste0(traits,"[+]"), 
                                    SNPs_ann_ps$assocs, fixed = FALSE), "dis",
                            ifelse(SNPs_ann_ps$no_effects == 1, "single","con"))
  
  ### categorize GxE into antagnoistic or concordant for both hosts
  SNPs_ann_ps$GxE <- ifelse(SNPs_ann_ps$dir == "dis" &
                              SNPs_ann_ps$host == "DZA only","antagonistic_DZA",
                            ifelse(SNPs_ann_ps$dir == "dis" &
                                     SNPs_ann_ps$host == "A17 only","antagonistic_A17",
                                   ifelse(SNPs_ann_ps$dir == "con" &
                                            SNPs_ann_ps$host == "DZA only","concordant_DZA", 
                                          ifelse(SNPs_ann_ps$dir == "con" &
                                                   SNPs_ann_ps$host == "A17 only","concordant_A17",
                                                 ## more than 2 assoc
                                                 ### DZA antagonistic
                                                 ifelse(SNPs_ann_ps$no_effects > 2 &
                                                          SNPs_ann_ps$dir == "dis" & 
                                                          grepl("DZA-1", SNPs_ann_ps$assocs, 
                                                                fixed = FALSE) &
                                                          grepl("DZA-3", SNPs_ann_ps$assocs, 
                                                                fixed = FALSE),
                                                        "antagonistic_DZA",
                                                        ### A17 antagonistic
                                                        ifelse(SNPs_ann_ps$no_effects > 2 &
                                                                 SNPs_ann_ps$dir == "dis" & 
                                                                 grepl("A17-2", SNPs_ann_ps$assocs, 
                                                                       fixed = FALSE) &
                                                                 grepl("A17-4", SNPs_ann_ps$assocs, 
                                                                       fixed = FALSE),
                                                               "antagonistic_A17",
                                                               ### DZA concordant
                                                               ifelse(SNPs_ann_ps$no_effects > 2 &
                                                                        SNPs_ann_ps$dir == "con" & 
                                                                        grepl("DZA-1", SNPs_ann_ps$assocs, 
                                                                              fixed = FALSE) &
                                                                        grepl("DZA-3", SNPs_ann_ps$assocs, 
                                                                              fixed = FALSE),
                                                                      "concordant_DZA",
                                                                      ### A17 antagonistic
                                                                      ifelse(SNPs_ann_ps$no_effects > 2 &
                                                                               SNPs_ann_ps$dir == "con" & 
                                                                               grepl("A17-2", SNPs_ann_ps$assocs, 
                                                                                     fixed = FALSE) &
                                                                               grepl("A17-4", SNPs_ann_ps$assocs, 
                                                                                     fixed = FALSE),
                                                                             "concordant_A17","other"))))))))
  
  
  ## create summary for counts in each cat
  SNPs_ann_ps.sum <- SNPs_ann_ps %>%
    group_by(GxE, host, dir) %>%
    summarize(count = n())
  
  ## save
  write.csv(SNPs_ann_ps, paste0("./Data_output/Var-lvl_summaries/SNPs_ann_ps_", traits, ".csv"), row.names = FALSE)
  
  return(list(SNPs_ann_ps.sum, SNPs_ann_ps))
  }
  
}