circos_func <- function(traits){
  
  print(paste0(traits))
  
  ## load genes (trait)
  load(file = paste0("./Data_output/Gene-lvl_summaries/genes_",traits,".Rdata")) ## loads SNPs_ann_gene
  trait <- SNPs_ann_gene
  trait$type <- "trait"
  trait <- trait %>%
    select(-host_exps)
  
  ## load genes (plast)
  load(file = paste0("./Data_output/Gene-lvl_summaries/genes_",traits,".plast.Rdata")) ## loads SNPs_ann_gene
  plast <- SNPs_ann_gene
  plast$cat <- plast$host
  plast$type <- "plast"
  
  df <- rbind(trait, plast)
  
  df$min_ps.t <- df$min_ps/1000000
  
  ### data subsets
  #### GxE
  gxe_1 <- df %>%
    filter(cat == "GxE" & exps == "1") %>%
    droplevels(.)
  gxe_2 <- df %>%
    filter(cat == "GxE" & exps == "2") %>%
    droplevels(.)
  gxe_3 <- df %>%
    filter(cat == "GxE" & exps == "3") %>%
    droplevels(.)
  gxe_4 <- df %>%
    filter(cat == "GxE" & exps == "4") %>%
    droplevels(.)
  #### GxG
  gxg_DZA <- df %>%
    filter(cat == "GxG" & host == "DZA only") %>%
    droplevels(.)
  gxg_A17 <- df %>%
    filter(cat == "GxG" & host == "A17 only") %>%
    droplevels(.)
  #### partial and uni
  pUni <- df %>%
    filter(cat == "Partially Universal") %>%
    droplevels(.)
  Uni <- df %>%
    filter(cat == "Universal") %>%
    droplevels(.)
  #### plasticity
  plast_DZA <- df %>%
    filter(host == "DZA only" & type == "plast") %>%
    droplevels(.)
  plast_A17 <- df %>%
    filter(host == "A17 only" & type == "plast") %>%
    droplevels(.)
  plast_both <- df %>%
    filter(host == "Both" & type == "plast") %>%
    droplevels(.)
  
  ## source functions
  source('../Source_code/circos_sector_func.R')
  source('../Source_code/circos_points_func.R')
  
  ### start save plot
  pdf(paste0("./Figures_tables/circos_plots/circos_",traits,".pdf"), 
      width = 5, height = 5, pointsize = 15)

  ### plotting code
  circos.par("track.height" = 0.15, circle.margin = c(0.15,0.15,0.5,0.35))
  circos.initialize(df$region, x = df$min_ps.t)
  ## add first track (GxE)
  circos.track(df$region, y = df$ave_abs_score,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter,
                             CELL_META$cell.ylim[2] + mm_y(18),
                             CELL_META$sector.index)
                 circos.axis(labels.cex = 0.7)
               })
  ## add points
  points_func(gxe_1, "#A4D466")
  points_func(gxe_3, "#4E820C")
  points_func(gxe_2, "#DE6BC4")
  points_func(gxe_4, "#AD208D")
  ## add sector highlights and text
  sector_func(0.778662, 0.813557,"Chromosome",1,"rsm-1")
  sector_func(0.908969, 0.948816,"Chromosome",1,"mot")
  sector_func(1.240720, 1.308367,"Chromosome",1,"rsm-2")
  sector_func(1.563939, 1.660057,"Chromosome",1,"16S-1")
  sector_func(2.117855, 2.180500,"Chromosome",1,"16S-2")
  sector_func(2.504874, 2.580316,"Chromosome",1,"16S-3")
  sector_func(0.014986, 0.070216,"pSymA",1,"fix")
  sector_func(0.227313, 0.273691,"pSymA",1,"nod/nif")
  sector_func(0.297595, 0.320934,"pSymA",1,"noe")
  sector_func(0.710907, 0.779780,"pSymA",1,"rhb")
  sector_func(0.000352, 0.054158,"pSymB",1,"dct/thi")
  sector_func(0.257079, 0.284461,"pSymB",1,"ccb/pqq")
  sector_func(0.500359, 0.553950,"pSymB",1,"nodU")
  sector_func(0.924092, 1.002841,"pSymB",1,"HK")
  sector_func(1.222343, 1.321863,"pSymB",1,"exo")
  sector_func(1.474514, 1.521610,"pSymB",1,"cyo/nad")
 
  ## add second track (GxG)
  circos.track(df$region, y = df$ave_abs_score)
  points_func(gxg_A17, "#C646A9")
  points_func(gxg_DZA, "#79AB39")
  ## add sector highlights
  sector_func(0.778662, 0.813557,"Chromosome",2,NULL)
  sector_func(0.908969, 0.948816,"Chromosome",2,NULL)
  sector_func(1.240720, 1.308367,"Chromosome",2,NULL)
  sector_func(1.563939, 1.660057,"Chromosome",2,NULL)
  sector_func(2.117855, 2.180500,"Chromosome",2,NULL)
  sector_func(2.504874, 2.580316,"Chromosome",2,NULL)
  sector_func(0.014986, 0.070216,"pSymA",2,NULL)
  sector_func(0.227313, 0.273691,"pSymA",2,NULL)
  sector_func(0.297595, 0.320934,"pSymA",2,NULL)
  sector_func(0.710907, 0.779780,"pSymA",2,NULL)
  sector_func(0.000352, 0.054158,"pSymB",2,NULL)
  sector_func(0.257079, 0.284461,"pSymB",2,NULL)
  sector_func(0.500359, 0.553950,"pSymB",2,NULL)
  sector_func(0.924092, 1.002841,"pSymB",2,NULL)
  sector_func(1.222343, 1.321863,"pSymB",2,NULL)
  sector_func(1.474514, 1.521610,"pSymB",2,NULL)
  
  ## add third track (pUni and uni)
  circos.track(df$region, y = df$ave_abs_score)
  points_func(pUni, "#9F7871")
  points_func(Uni, "#3B2C2A")
  ## add sector highlights
  sector_func(0.778662, 0.813557,"Chromosome",3,NULL)
  sector_func(0.908969, 0.948816,"Chromosome",3,NULL)
  sector_func(1.240720, 1.308367,"Chromosome",3,NULL)
  sector_func(1.563939, 1.660057,"Chromosome",3,NULL)
  sector_func(2.117855, 2.180500,"Chromosome",3,NULL)
  sector_func(2.504874, 2.580316,"Chromosome",3,NULL)
  sector_func(0.014986, 0.070216,"pSymA",3,NULL)
  sector_func(0.227313, 0.273691,"pSymA",3,NULL)
  sector_func(0.297595, 0.320934,"pSymA",3,NULL)
  sector_func(0.710907, 0.779780,"pSymA",3,NULL)
  sector_func(0.000352, 0.054158,"pSymB",3,NULL)
  sector_func(0.257079, 0.284461,"pSymB",3,NULL)
  sector_func(0.500359, 0.553950,"pSymB",3,NULL)
  sector_func(0.924092, 1.002841,"pSymB",3,NULL)
  sector_func(1.222343, 1.321863,"pSymB",3,NULL)
  sector_func(1.474514, 1.521610,"pSymB",3,NULL)
  
  ## add fourth track (plasticity)
  circos.track(df$region, y = df$ave_abs_score)
  points_func(plast_DZA, "#4E820C")
  points_func(plast_A17, "#AD208D")
  points_func(plast_both, "#9F7871")
  ## add sector highlights
  sector_func(0.778662, 0.813557,"Chromosome",4,NULL)
  sector_func(0.908969, 0.948816,"Chromosome",4,NULL)
  sector_func(1.240720, 1.308367,"Chromosome",4,NULL)
  sector_func(1.563939, 1.660057,"Chromosome",4,NULL)
  sector_func(2.117855, 2.180500,"Chromosome",4,NULL)
  sector_func(2.504874, 2.580316,"Chromosome",4,NULL)
  sector_func(0.014986, 0.070216,"pSymA",4,NULL)
  sector_func(0.227313, 0.273691,"pSymA",4,NULL)
  sector_func(0.297595, 0.320934,"pSymA",4,NULL)
  sector_func(0.710907, 0.779780,"pSymA",4,NULL)
  sector_func(0.000352, 0.054158,"pSymB",4,NULL)
  sector_func(0.257079, 0.284461,"pSymB",4,NULL)
  sector_func(0.500359, 0.553950,"pSymB",4,NULL)
  sector_func(0.924092, 1.002841,"pSymB",4,NULL)
  sector_func(1.222343, 1.321863,"pSymB",4,NULL)
  sector_func(1.474514, 1.521610,"pSymB",4,NULL)
  
  ## add histrogram
  circos.trackHist(df$region, df$min_ps.t, bin.size = 0.1, col = NA)
  
  ## clear
  circos.clear()
  
  ## end save plot
  dev.off()
  
}