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
#SNP-level filtering: minor allele frequency.
counted <- apply(genotypes.clean, 1, function(x) (min(table(x))/sum(!is.na(x))))
counted.clean = counted[counted > 0.01]
print(paste("We've removed :", 25000-length(counted.clean),"variants"))
cou<-as.data.frame(counted)
ggplot(data = cou, mapping=aes(x=counted)) + geom_histogram(binwidth = 0.011)
#3 Genome Wide Association Studies.
covariates.txt$gender <- as.factor(covariates.txt$gender)
linmodel <- merge(phenotypes.txt, covariates.txt)
summary(lm("Cholesterol ~ gender", data=linmodel))
r_squared <- summary(lm("Cholesterol ~ gender", data=linmodel))$r.squared
print("Since the R-squared is close to 0, we can assume that there is no obvious relation between the gender and the cholesterol level")
#Boxplot
ggplot(data = linmodel, mapping=aes(y=Cholesterol,x=gender)) + geom_boxplot()
ggplot(linmodel, aes(Cholesterol)) + geom_density(aes(col = 'all')) + geom_density(aes(col = gender)) + geom_density()
print("No because the distribution is homogenous to a normal distribution, which is random")
#1 = alt, 0 = ref
data.for.pca = t(genotypes.clean)
colnames(data.for.pca) = genotypes.vcf$V3[2:25001]
data.for.pca[data.for.pca == "0/0"] = 0
data.for.pca[data.for.pca == "0/1"] = 1
data.for.pca[data.for.pca == "1/1"] = 2
data.for.pca[is.na(data.for.pca)] <- 3
str(data.for.pca)
storage.mode(data.for.pca) <- "numeric"
str(data.for.pca)
data.pca = prcomp(data.for.pca, center = T)
#data.plot = data.pca$x
#plot(data.plot[,1], data.plot[,2])

df <- data.frame(data.pca$x)
ggplot(df,aes(x=PC1,y=PC2)) + geom_point()

