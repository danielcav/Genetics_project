#step 1 : setting environment
library("ggplot2")
library("data.table")
#step 2 : loading data
genotypes.vcf <- fread("genotypes.vcf", header = F, data.table = F, na.string = ".")
phenotypes.txt <- fread("phenotypes.txt", header = T, data.table = F)
covariates.txt <- fread("covariates.txt", header = T, data.table = F)
#Clean the file in order to obtain only the wanted data.
genotypes.clean <- genotypes.vcf[2:nrow(genotypes.vcf),10:ncol(genotypes.vcf)]
#Calculate the call rate.
call.rate <- apply(genotypes.clean, 1, function(x) sum(!is.na(x))/ncol(genotypes.clean))
#We keep the variants with a 100% call rate.
variants.clean = call.rate[call.rate == 1.]
#We plot the histogram.
call.rate <- as.data.frame(call.rate)
call.rate$idx = 1:25000
ggplot(data = call.rate, mapping=aes(x=call.rate)) + geom_histogram()
#How many variants are removed ? 
print(paste("We've removed :", 25000-length(variants.clean),"variants"))
#3. SNP-level filtering: minor allele frequency.
