sig_snp_func <- function(combs, trait_y, trait_x, line, df){
  
  # print combinations
  print(paste0(combs)) # x vs y
  
  # format traits
  sig_y <- paste0("sig.", trait_y,".", line)
  scores_y <- paste0("beta.", trait_y,".", line)
  sig_x <- paste0("sig.", trait_x,".", line)
  scores_x <- paste0("beta.", trait_x,".", line)
  
  ## subset to significant SNPs (at least one trait)
  sig_snps <- df %>%
    filter(get(sig_y) == "TRUE" | get(sig_x) == "TRUE") %>%
    droplevels(.)
  
  ## add col for plant line
  sig_snps$line <- paste0(line)
  sig_snps$line <- factor(sig_snps$line, levels = c("DZA","A17"))
  
  # add in effect direction (major allele!!!)
  sig_snps$trait_y_dir <- ifelse(sig_snps[[scores_y]] > 0, "-","+")
  sig_snps$trait_x_dir <- ifelse(sig_snps[[scores_x]] > 0, "-","+")
  
  # add cols for pw trait combs
  sig_snps$trait_comb <- paste0(trait_y,"_vs_",trait_x)
  sig_snps$trait_comb_line <- paste0(combs)
  sig_snps$trait_comb_dir <- paste0(trait_y,sig_snps$trait_y_dir,"_vs_",
                                    trait_x,sig_snps$trait_x_dir,"_",line)
  
  ## significant SNPs, either trait
  sig_snps$sig_traits <- ifelse(sig_snps[[sig_y]] == "TRUE" &
                              sig_snps[[sig_x]] == "TRUE", "Both_traits",
                            ifelse(sig_snps[[sig_y]] == "TRUE" &
                                     sig_snps[[sig_x]] == "FALSE", paste0(trait_y, "_only"),
                                   paste0(trait_x, "_only")))
  ## order
  sig_snps$sig_traits <- factor(sig_snps$sig_traits, levels = c(paste0(trait_x, "_only"), 
                                                        paste0(trait_y, "_only"), 
                                                        "Both_traits"))
  
  ## specify whether variants are pleiotropic or not:
  sig_snps$pleio <- ifelse(sig_snps$sig_traits == "Both_traits","Y","N")
  
  # categorize effects
  sig_snps$effect <- ifelse(sig_snps[[scores_x]] > 0 & 
                              sig_snps[[scores_y]] > 0, "Top_right", 
                            ifelse(sig_snps[[scores_x]] < 0 & 
                                     sig_snps[[scores_y]] > 0, "Top_left", 
                                   ifelse(sig_snps[[scores_x]] > 0 &
                                            sig_snps[[scores_y]] < 0, "Bottom_right",
                                          "Bottom_left")))
  ## order
  sig_snps$effect <- factor(sig_snps$effect, 
                            levels = c("Top_right", "Bottom_right", "Bottom_left", "Top_left"))
  
  ## categorize into diagonal (D) and off-diagonal (OD)
  sig_snps$cat <- ifelse(sig_snps$effect == "Top_left" |
                           sig_snps$effect == "Bottom_right",
                         "OD","D")
  ## order
  sig_snps$cat <- factor(sig_snps$cat, levels = c("D","OD"))
  
  ## add scores to the df
  sig_snps$betas_y <- sig_snps[[scores_y]]
  sig_snps$betas_x <- sig_snps[[scores_x]]
  
  ## correct to put into two quadrants
  sig_snps$betas_y.c <- abs(sig_snps$betas_y)
                               
  sig_snps$betas_x.c <- ifelse(sig_snps$cat == "D", abs(sig_snps$betas_x),
                               (abs(sig_snps$betas_x))*-1)
  
  ## create line_cat
  sig_snps$line_cat <- paste0(sig_snps$line, "_", sig_snps$cat)
  ## order
  sig_snps$line_cat <- factor(sig_snps$line_cat, levels = c("DZA_D", "A17_D",
                                                            "DZA_OD","A17_OD"))
  
  ## specify type of pleiotropy (rhiz-only versus symbiotic)
  sig_snps$pleio_type <- ifelse(grepl("shoot", sig_snps$trait_comb) | 
                                  grepl("chloro", sig_snps$trait_comb),
                                "sym","rhiz")
  
  ## categorize pleiotropic effects (Al = alignment; C = conflict; An = antagonistic; P = positive)
  sig_snps$pleio_cat <- ifelse(sig_snps$pleio_type == "sym" & sig_snps$effect == "Top_right" |
                                 sig_snps$pleio_type == "sym" & sig_snps$effect == "Bottom_left", "Al",
                               ifelse(sig_snps$pleio_type == "sym" & sig_snps$effect == "Top_left" |
                                        sig_snps$pleio_type == "sym" & sig_snps$effect == "Bottom_right", "C",
                                      ifelse(sig_snps$pleio_type == "rhiz" & sig_snps$effect == "Top_right" |
                                               sig_snps$pleio_type == "rhiz" & sig_snps$effect == "Bottom_left",
                                             "P","An")))
  ## order
  sig_snps$pleio_cat <- factor(sig_snps$pleio_cat, levels = c("Al","C","P","An"))
  
  ## rm unwanted columns
  sig_snps.f <- sig_snps %>%
    filter(pleio == "Y") %>%
    droplevels(.) %>%
    select(-contains("beta."),-contains("sig."),-trait_y_dir,-trait_x_dir)
  
  ## write to csv
  write.csv(sig_snps.f, file = paste("./Data_output/sig_snps_files/sig_snps_", 
                                     combs, ".csv"), 
            row.names = FALSE)
  
  return(list(sig_snps))
  
}