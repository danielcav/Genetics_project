#Step 1 : setting the environment
library("ggplot2")
library("data.table")
#Step 2 : loading data
genotypes.vcf <- fread("genotypes.vcf", header = F, data.table = F, na.string = ".")
phenotypes.txt <- fread("phenotypes.txt", header = T, data.table = F)
covariates.txt <- fread("covariates.txt", header = T, data.table = F)

#2 Data preprocessing.
#SNP-level filtering: call rate
#Clean the file in order to obtain only the wanted data.
genotypes.clean <- genotypes.vcf[2:nrow(genotypes.vcf),10:ncol(genotypes.vcf)]
#Calculate the call rate for every SNP.
call.rate <- apply(genotypes.clean, 1, function(x) sum(!is.na(x))/ncol(genotypes.clean))
#Keep the variants with a 100% call rate.
variants.clean = call.rate[call.rate == 1.]
genotypes.without.na <- genotypes.clean[rownames(genotypes.clean)%in%names(variants.clean),]
#Plot the histogram.
call.rate <- as.data.frame(call.rate)
call.rate$idx = 1:25000
ggplot(data = call.rate, mapping=aes(x=call.rate)) + geom_histogram() + ggtitle("Call rate distribution") + xlab("Call rate") + ylab("Number of SNPs")
#How many variants are removed ? 
print(paste("We've removed :", 25000-length(variants.clean),"variants"))

#SNP-level filtering: minor allele frequency.
#Function that computes the MAF values.
minors <- function (x){
  one.one = length(which(x == "1/1"))
  one.zero = length(which(x == "0/1")) + length(which(x == "1/0")) 
  zero.zero = length(which(x== "0/0"))
  nb.zero = zero.zero*2 + one.zero
  nb.one = one.one*2 + one.zero
  if(nb.zero == 0){return(0)}
  if (nb.zero<nb.one){
    return((nb.zero)/(2*length(x)))
  }else{return(nb.ones/(2*length(x)))}
}
#Calculate the MAF for every SNP.
#We also calculate the MAF for SNPs with a call rate < 1.
maf <- apply(genotypes.without.na, 1, function(x) (minors(x)))
maf.clean = maf[maf > 0.01]
#Filter the genotype data frame and plot the histogram.
genotypes.filtered <- genotypes.without.na[rownames(genotypes.without.na)%in%names(maf.clean),]
ggplot(data = as.data.frame(maf), mapping=aes(x=maf)) + geom_histogram(binwidth = 0.01) + ggtitle("Minor allele frequency distribution") + xlab("Minor allele frequency") + ylab("SNPs")
#How many variants are removed ?
print(paste("We've removed :", nrow(genotypes.without.na)-length(maf.clean),"variants"))

#3 Genome Wide Association Studies (GWAS).
#Set the gender as a factor.
covariates.txt$gender <- as.factor(covariates.txt$gender)
#Create the model
linmodel <- merge(phenotypes.txt, covariates.txt)
#Coefficient of determination.
r_squared <- summary(lm("Cholesterol ~ gender", data=linmodel))$r.squared
# Based on the statistics, do you expect the gender to impact the cholesterol level ?
print(paste0("Since the R-squared is close to 0 (",r_squared,") we can assume that there is no obvious relation between the gender and the cholesterol level"))

#Plot the boxplot.
ggplot(data = linmodel, mapping=aes(y=Cholesterol,x=gender)) + geom_boxplot() + xlab("Gender") + ylab("Cholesterol level") + ggtitle("Boxplot of the cholesterol level in different genders")
#Plot the distribution of cholesterol for different genders.
ggplot(linmodel, aes(Cholesterol)) + geom_density() + geom_density(aes(col = gender)) + geom_density() + ggtitle("Distribution of cholesterol for different genders") + xlab("Cholesterol level") + ylab("Density")
# Is gender a covariate for cholesterol level? Why?
print("No, both distributions are homogenous to a normal distribution, which is random.")

#Principal Components Analysis (PCA).
#Create the data for a PCA.
data.for.pca = as.data.frame(t(genotypes.filtered))
data.for.pca[data.for.pca == "0/0"] = 0
data.for.pca[data.for.pca == "0/1"] = 1
data.for.pca[data.for.pca == "1/0"] = 1
data.for.pca[data.for.pca == "1/1"] = 2
data.for.pca.clean <- as.matrix(data.for.pca)
storage.mode(data.for.pca.clean) <- "numeric"
#Calculate the principal components (PCs).
data.pca = prcomp(data.for.pca.clean)
#Plot the first and second principal components.
df <- data.pca$x
ggplot(as.data.frame(df),aes(x=PC1,y=PC2)) + geom_point() + ggtitle("Clusters after PCA") 
#How many clusters do you see ?
print("We can clearly observe 3 different clusters.")
#If more than one, what do the clusters represent ?
print("The clusters represent the different genotypes (3 genotypes ->3 clusters).")
#Should we correct for the population structure? Why?
print("No, we don't need to correct the population structure because the data correctly represents different groups of people.")

#GWAS.
#Function returning the coefficient of association and the p-value.
beta_p <- function(x){
  temp <- summary(lm(phenotypes.txt$Cholesterol ~ as.numeric(x)))
  return(c(temp$coefficients[2,1],-log10(temp$coefficients[2,4])))#(beta, p-value)
}
#Create the data for a GWAS.
data.for.gwas <- as.data.frame(t(data.for.pca))
data.gwas <- apply(data.for.gwas, 1, function(x) beta_p(x))
#Rank of SNPs.
data.gwas <- rbind(data.gwas,c(1:ncol(data.gwas)))
chr.pos <- genotypes.vcf[rownames(genotypes.vcf) %in% colnames(data.gwas),1]
data.gwas <- rbind(data.gwas,as.numeric(chr.pos))
data.gwas.t <- as.data.frame(t(data.gwas))
#Bonferroni corrected threshold with α = 0.05.
threshold <- -log10(0.05/nrow(data.gwas.t))
data.gwas.t$color[(data.gwas.t$V4 %% 2) == 0] = "odd"
data.gwas.t$color[(data.gwas.t$V4 %% 2) != 0] = "even"
data.gwas.t$color[data.gwas.t$V2 > threshold] = "significant" #Bonus

#Manhattan plot
ggplot(as.data.frame(data.gwas.t),aes(x=V3, y=V2, col = factor(color))) + geom_point() + geom_hline(yintercept=threshold, linetype="dashed", color = "red") + ggtitle("Manhattan plot without covariates") + xlab("Chromosome position") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "None") + ylab("-log10(p-value)")

#Repeat the GWAS and Manhattan plot considering top 10 principal components as covariates.
#Extracting the p-value
beta_p_cov <- function(x){
  temp <- summary(lm(phenotypes.txt$Cholesterol ~ as.numeric(x) + df[,1:10]))
  return(c(temp$coefficients[2,1],-log10(temp$coefficients[2,4])))
}
#Create the data for the GWAS.
gwas.cov <- apply(data.for.gwas, 1, function(x) beta_p_cov(x))
#Rank of SNPs
gwas.cov <- rbind(gwas.cov,c(1:ncol(gwas.cov)))
gwas.cov <- rbind(gwas.cov,as.numeric(chr.pos))
gwas.cov.t <- as.data.frame(t(gwas.cov))
gwas.cov.t$color[(gwas.cov.t$V4 %% 2) == 0] = "odd"
gwas.cov.t$color[(gwas.cov.t$V4 %% 2) != 0] = "even"
gwas.cov.t$color[gwas.cov.t$V2 > threshold] = "significant"
ggplot(as.data.frame(gwas.cov.t),aes(x=V3, y=V2, col = factor(color))) + geom_point() + geom_hline(yintercept=threshold, linetype="dashed", color = "red")  + ggtitle("Manhattan plot with covariates") + xlab("Chromosome position") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "None") + ylab("-log10(p-value)")

#Compare the results with and without covariates.
print("We see that the GWAS with covariates have less significant values (3 vs 4 without covariates).")

#Bonus: QQ plot.
#Generate the expected values .
generated.p.values <- sort(-log10(ppoints(nrow(gwas.cov.t))))
#Observed p-values with covariates.
observed.p.values <- sort(gwas.cov.t$V2)
qq.data <- cbind(generated.p.values, observed.p.values)
#QQ plot.
ggplot(as.data.frame(qq.data), aes(x=observed.p.values,y=generated.p.values)) + geom_point() + geom_abline() + ggtitle("QQ-plot") + xlab("Observed values (-log10 scale)") + ylab("Generated values (-log10 scale)")