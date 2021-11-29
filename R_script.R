#step 1 : setting environment
library("ggplot2")
library("data.table")

#step 2 : loading data
genotypes.vcf <- fread("genotypes.vcf", header = F, data.table = F, na.string = ".")
phenotypes.txt <- fread("phenotypes.txt", header = T, data.table = F)
covariates.txt <- fread("covariates.txt", header = T, data.table = F)

sub1 <- genotypes.vcf[genotypes.vcf[2:nrow(genotypes.vcf),10:ncol(genotypes.vcf)] == "NA" ]

apply(z, 1, function(x) sum(is.na(x)))