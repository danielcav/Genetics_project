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
genotypes.without.na <- genotypes.clean[rownames(genotypes.clean)%in%names(variants.clean),]
#We plot the histogram.
call.rate <- as.data.frame(call.rate)
call.rate$idx = 1:25000
ggplot(data = call.rate, mapping=aes(x=call.rate)) + geom_histogram()
#How many variants are removed ? 
print(paste("We've removed :", 25000-length(variants.clean),"variants"))
#SNP-level filtering: minor allele frequency.
sum(!is.na(genotypes.clean[1,]))
minors <- function (x){
  one.one = length(which(x == "1/1"))
  one.zero = length(which(x == "0/1"))
  zero.zero = length(which(x== "0/0"))
  nb.zero = zero.zero*2 + one.zero
  nb.one = one.one*2 + one.zero
  if(nb.zero == 0){return(0)}
  if (nb.zero<nb.one){
    return((nb.zero)/(2*sum(!is.na(x))))
  }else{return(nb.ones/(2*sum(!is.na(x))))}
}
counted <- apply(genotypes.clean, 1, function(x) (minors(x)))
counted.clean = counted[counted > 0.01]
genotypes.filtered <- genotypes.without.na[rownames(genotypes.without.na)%in%names(counted.clean),]
print(paste("We've removed :", 25000-length(counted.clean),"variants"))
cou<-as.data.frame(counted)
ggplot(data = cou, mapping=aes(x=counted)) + geom_histogram(binwidth = 0.01)
#3 Genome Wide Association Studies.
covariates.txt$gender <- as.factor(covariates.txt$gender)
linmodel <- merge(phenotypes.txt, covariates.txt)
summary(lm("Cholesterol ~ gender", data=linmodel))
r_squared <- summary(lm("Cholesterol ~ gender", data=linmodel))$r.squared
print("Since the R-squared is close to 0, we can assume that there is no obvious relation between the gender and the cholesterol level")
#Boxplot
ggplot(data = linmodel, mapping=aes(y=Cholesterol,x=gender)) + geom_boxplot()
ggplot(linmodel, aes(Cholesterol)) + geom_density(aes(col = 'all')) + geom_density(aes(col = gender)) + geom_density()
print("No because the distribution is homogenous to a normal distribution, which is random.")
#1 = alt, 0 = ref
data.for.pca = as.data.frame(t(genotypes.filtered))
data.for.pca[data.for.pca == "0/0"] = 0
data.for.pca[data.for.pca == "0/1"] = 1
data.for.pca[data.for.pca == "1/0"] = 1
data.for.pca[data.for.pca == "1/1"] = 2
val <- as.matrix(data.for.pca)
storage.mode(val) <- "numeric"
data.pca = prcomp(val)
df <- data.pca$x
ggplot(as.data.frame(df),aes(x=PC1,y=PC2)) + geom_point()
#GWAS:
data.for.gwas <- as.data.frame(t(data.for.pca))
beta_p <- function(x){
  temp <- summary(lm(phenotypes.txt$Cholesterol ~ as.numeric(x)))
  return(c(temp$coefficients[2,1],-log10(temp$coefficients[2,4])))#beta,pvalue
}
data.gwas <- apply(data.for.gwas, 1, function(x) beta_p(x))
data.gwas <- rbind(data.gwas,c(1:ncol(data.gwas)))
chr.pos <- genotypes.vcf[rownames(genotypes.vcf) %in% colnames(data.gwas),1]
data.gwas <- rbind(data.gwas,as.numeric(chr.pos))
data.gwas.t <- as.data.frame(t(data.gwas))
treshold <- -log10(0.05/nrow(data.gwas.t))
data.gwas.t$color[(data.gwas.t$V4 %% 2) == 0] = "odd"
data.gwas.t$color[(data.gwas.t$V4 %% 2) != 0] = "even"
data.gwas.t$color[data.gwas.t$V2 > treshold] = "significant"
ggplot(as.data.frame(data.gwas.t),aes(x=V3, y=V2, col = factor(color))) + geom_point() + geom_hline(yintercept=treshold, linetype="dashed", color = "red")
#6
