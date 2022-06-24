comp.sigs_func <- function(traits_regions, traits, regions){
  
    ## print traits
    print(paste0(traits_regions))
                      
    ## load data
    mlmm <- read.delim(paste0("./Data_input/betasigs_29Mar2021/", regions, "_subset191_LD_mlmm_", traits, ".assoc.txt"))
    load(paste0("./Data_input/betasigs_29Mar2021/", regions, "_betasigs.Rdata")) ## realres.psFIsig
    
    ulmm <- realres.psFIsig %>%
      select(rs, contains(paste0(traits)), -contains("plast"))
    
    ## add in gene cats
    
    ### link up gene and var-level information
    var <- read_csv(paste0("./Data_output/Var-lvl_summaries/SNPs_ann_ps_", traits, ".csv"))
    gene <- read_csv(paste0("./Data_output/Gene-lvl_summaries/genes_", traits, ".csv"))
    
    ## match info
    ulmm <- left_join(ulmm, var[,c("rs","RefSeq_ID")], by = "rs")
    ulmm <- left_join(ulmm, gene[,c("RefSeq_ID","cat","host")], by = "RefSeq_ID")
    
    ## merge
    comb <- left_join(mlmm, ulmm, by = "rs")
    
    return(list(comb))
    
}