ttest_func <- function(combs, stats, sizes, cats, df){
  
  print(paste0(combs))
  
  ## subset real df to stat and cat
  ## df subsets 
  
  if (cats == "GWAS"){
  real <- df %>%
    filter(cat != "all") %>%
    droplevels(.)
  }
  
  if (cats == "sig"){
  real <- df %>%
    filter(cat == "sig" | cat == "pleio") %>%
    droplevels(.)
  }
  
  if (cats == "alignment"){
  real <- df %>%
    filter(sym == "alignment") %>%
    droplevels(.)
  }
  
  if (cats == "conflict"){
  real <- df %>%
    filter(sym == "conflict") %>%
    droplevels(.)
  }
  
  if (cats == "positive"){
  real <- df %>%
    filter(rhiz == "positive") %>%
    droplevels(.)
  }
  
  if (cats == "antagonistic"){
  real <- df %>%
    filter(rhiz == "antagonistic") %>%
    droplevels(.)
  }
  
  ## randomly sample gene sets according to sample size in focal gene cat
  ttests_out <- sapply(1:1000, function(x){
    sample <- sample(df[[stats]], size = sizes, replace = FALSE)
    ttest <- t.test(sample, real[[stats]])
    return(ttest$p.value)
  })
  
  ## format, calculate global p
  ttests <- as.data.frame(ttests_out)
  ttests$cutoff <- ifelse(ttests$ttests_out < 0.05, "sig","NS")
  ttests$cutoff <- factor(ttests$cutoff, levels = c("NS","sig"))
  ttests_sum <- ttests %>%
    group_by(cutoff, .drop = FALSE) %>%
    summarize(counts = n())

  g_pval <- ttests_sum[[1,2]]/1000
  return(g_pval)
}
