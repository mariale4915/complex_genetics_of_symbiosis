mean_func <- function(cats, df){
  
  print(paste0(cats))
  
  ## create long version
  df.l <- df %>%
    gather(stat, vals, Pi:Fu.Li.D)
  
  ## subset real df to stat and cat
  ## df subsets 
  
  if (cats == "GWAS"){
    real <- df.l %>%
      filter(cat != "all") %>%
      droplevels(.)
  }
  
  if (cats == "sig"){
    real <- df.l %>%
      filter(cat == "sig" | cat == "pleio") %>%
      droplevels(.)
  }
  
  if (cats == "alignment"){
    real <- df.l %>%
      filter(sym == "alignment") %>%
      droplevels(.)
  }
  
  if (cats == "conflict"){
    real <- df.l %>%
      filter(sym == "conflict") %>%
      droplevels(.)
  }
  
  if (cats == "positive"){
    real <- df.l %>%
      filter(rhiz == "positive") %>%
      droplevels(.)
  }
  
  if (cats == "antagonistic"){
    real <- df.l %>%
      filter(rhiz == "antagonistic") %>%
      droplevels(.)
  }
  
  ## calculate the mean values for each stat
  real.sum <- real %>%
    group_by(stat) %>%
    summarize(mean = mean(vals))
  
  real.sum$cat <- paste0(cats)
  
  return(real.sum)
}