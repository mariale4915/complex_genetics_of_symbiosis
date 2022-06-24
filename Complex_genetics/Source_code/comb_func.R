comb_func <- function(traits_exps, traits, nums){
  
  ## print trait_exp being iterated
  print(paste0(traits_exps))
  
  ## import files for each trait (across experiments)
  file_list <- list.files(path="./Data_input/betasigs_29Mar2021/", 
                          pattern = paste0(nums, ".assoc"), full.names = TRUE)
  
  traits_out <- sapply(file_list, function(files){
    traits <- read.delim(files) # read in the file
    return(list(traits))
  })
  
  ## combine regions for each trait
  traits_comb <- bind_rows(traits_out)
  return(list(traits_comb))
}