## harmonic mean test code:

library(harmonicmeanp)
setwd("C:/Users/IGB/Dropbox/complex_genetics_of_symbiosis/Complex_genetics/GWAS")

### load for DZA experiment 1
gwas_chr <- read.delim("./Data_input/betasigs_29Mar2021/chr_subset191_LD_ulmm_trait_13.assoc.txt")
gwas_psyma <- read.delim("./Data_input/betasigs_29Mar2021/psyma_subset191_LD_ulmm_trait_13.assoc.txt")
gwas_psymb <- read.delim("./Data_input/betasigs_29Mar2021/psymb_subset191_LD_ulmm_trait_13.assoc.txt")

### combine across regions
gwas <- rbind(gwas_chr, gwas_psyma, gwas_psymb)

### Harmonic mean p across whole genome
R <- 1:nrow(gwas)
L <- 7965
gwas$w <- 1/L
(HMP.R <- sum(gwas$w[R])/sum(gwas$w[R]/gwas$p_lrt[R])) ## 0.1061016
# Specify the false positive rate
alpha <-  0.05
# Compute the HMP significance threshold
(alpha.L  <-  qharmonicmeanp(alpha, L))
# Test whether the HMP for subset R is significance
w.R <- sum(gwas$w[R]) ## equals 1 (bc 7965 / 1/7965)
(alpha.L * w.R) ## same as the harmonic mean alpha
## Therefore after adjusting for multiple comparison we cannot reject the null hypothesis of no association across the genome at level α=0.05 because 0.03139138 (HMP.R) is not below 0.03139138 (HMP significance threshold for subset).

# Note that the p.hmp function has been redefined to take argument L. Omitting L will issue a warning.
(p.hmp(gwas$p_lrt[R],gwas$w[R],L))
## The combined p-value across the genome is useful because if the combined p-value is not significant, neither is any constituent p-value, after multiple testing correction, as always. 

### subset to the chr
R  <- which(gwas$chr=="chr")
(p.R <- p.hmp(gwas$p_lrt[R],gwas$w[R],L))
w.R <- sum(gwas$w[R])
alpha*w.R
(p.R.adjust.chr <- p.R/w.R) ## NS
## for positions on the chrom, the combined p-value was 0.08826114 which was above the significance threshold of 0.004413057 - thus, the test is non-significant
## subset to the psyma
R  <- which(gwas$chr=="psyma")
(p.R <- p.hmp(gwas$p_lrt[R],gwas$w[R],L))
w.R <- sum(gwas$w[R])
alpha*w.R
(p.R.adjust.psyma <- p.R/w.R) ## NS
## subset to the psymb
R  <- which(gwas$chr=="psymb")
(p.R <- p.hmp(gwas$p_lrt[R],gwas$w[R],L))
w.R <- sum(gwas$w[R])
alpha*w.R
(p.R.adjust.psymb <- p.R/w.R) ## NS