popgenome_func <- function(group){
  
  # Load the data
  snp <- readData(paste0("./FASTAS",group))
  
  ## This is complex object, with several slots
  sum <- get.sum.data(snp)
  sum <- as.data.frame(sum)
  sum$files <- rownames(sum)
  sum <- sum %>% separate(files, c("gene_code", "type"), sep = "[.]") %>%
    select(-type)
  
  ## Diversities (by gene)
  snp <- diversity.stats(snp, pi = TRUE) # this does the calculations and 
  # adds the results to the appropriate slots
  
  ## get diversities
  div <- get.diversity(snp)[[1]]
  div <- as.data.frame(div)
  div$files <- rownames(div)
  div <- div %>% separate(files, c("gene_code", "type"), sep = "[.]") %>%
    select(-type, -contains("F_ST"))
  
  ## combine
  sum_div <- full_join(sum, div, by = "gene_code")
  
  ## neutrality stats
  snp <- neutrality.stats(snp)
  neut <- get.neutrality(snp)[[1]]
  neut <- as.data.frame(neut)
  neut$files <- rownames(neut)
  neut <- neut %>% separate(files, c("gene_code", "type"), sep = "[.]") %>%
    select(-type, -contains("Roza"),-contains("F_S"),-contains(".H"),
           -contains(".E"),-contains(".S"))
  
  ## combine
  real_gene_stats <- full_join(sum_div, neut, by = "gene_code")
  
  ## save to out
  return(list(real_gene_stats))
}