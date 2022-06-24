qq_func <- function(traits){
  
  ## print traits
  print(paste0(traits))
                   
  ## qq plot function
  source('../Source_code/qqunif.plot.R')
  
  traits_comb <- comb_out[[traits]]
  
  my.pvalue.list<-list("Score"=traits_comb$p_score, 
                      "LRT"=traits_comb$p_lrt, 
                      "Wald"=traits_comb$p_wald)
  ## make qq-plots
  qq_plot <- qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)),
                        main = paste0(traits))
  
  return(qq_plot)
  
}