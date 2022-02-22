## vcf table file for Rizwan

library(tidyverse)

## get strain names
GxG_data <- read_csv("../../Complex_genetics/Phenotypic_analyses/Raw_data/GxG_data.csv")
GxG_data$strain_ID <- as.factor(GxG_data$strain_ID)
strains <- levels(GxG_data$strain_ID) 
strains_sub <- as.data.frame(strains[!grepl("C", strains)])
### add MAG to strains
strains_sub$add <- "MAG"
strains_sub$add1 <- ".GT"
strains_sub$strains1 <- str_c(strains_sub$add, 
                              strains_sub$'strains[!grepl(\"C\", strains)]', sep = "", collapse = NULL)
strains_sub$strains <- str_c(strains_sub$strains1, 
                             strains_sub$add1, sep = "", collapse = NULL)
### convert to character
strains_inc <- as.character(strains_sub$strains)

## add in snpEff info:
VCF_ann <- read.delim("./Data_input/vcf_annotations_29Mar2021/subset191_LD.ann.TABLE")
## extract gene, impact, function from snpEff annotation
VCF_ann$ANN <- as.character(VCF_ann$ANN)
VCF_ann <- separate(data = VCF_ann, col = ANN, 
                    into = c("allele", "effect","impact","gene_name","feature_type"), 
                    sep = "\\|")
## create rs
VCF_ann$region <- ifelse(VCF_ann$CHROM == "NZ_CP021797.1", "chr", 
                         ifelse(VCF_ann$CHROM == "NZ_CP021798.1", 
                                "psyma","psymb"))
VCF_ann$add1 <- "freebayes"
VCF_ann$add2 <- "snp"
VCF_ann$rs <- str_c(VCF_ann$add1, VCF_ann$add2, 
                    VCF_ann$region, VCF_ann$POS,
                    VCF_ann$ALT,sep = "-", 
                    collapse = NULL)

## select colnames based on identifiers
VCF_ann.sub <- VCF_ann %>%
  select(rs, REF, ALT, matches(strains_inc))
### use strain names as rownames
row.names(VCF_ann.sub) <- VCF_ann.sub$rs

VCF_ann.sub <- VCF_ann.sub %>% mutate_if(is.factor,as.character)

## replace with 0's and 1's
VCF_recode <- sapply(strains_inc, df = VCF_ann.sub, function(cols, df){
  ifelse(df[[cols]] == df$REF, 0, 1)
})

VCF <- cbind(VCF_ann.sub %>% select(rs), as.data.frame(VCF_recode))

## add in trait info
SNPs_ann_ps_all <- read_csv("./Data_output/Var-lvl_summaries/SNPs_ann_ps_all.csv")
### left join 
SNP_file <- left_join(SNPs_ann_ps_all, VCF, by = "rs")

## save
write.csv(SNP_file, file = "./Data_output/SNPs_RNAseq.csv", row.names = FALSE)

## filter to only shoot DZA (example)

SNP_file$shoot_DZA <- ifelse(SNP_file$host == "DZA only" &
                               grepl("shoot+", SNP_file$assocs, fixed = TRUE), "YES",
                             ifelse(SNP_file$host == "DZA only" &
                                      grepl("shoot-", SNP_file$assocs, fixed = TRUE),
                                "YES", "NO"))