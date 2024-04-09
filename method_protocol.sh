########################################################################################################################
# Documents and scripts were written by: Ruirui Fu
# For manuscript: Ruirui Fu, et al. (2024). "Shared xerophytic genes and their re-use in local adaptation to aridity in a desert plant Gymnocarpos przewalskii" (under review). 
# email: furuirui@zju.edu.cn
# Jun Chen Lab. 
########################################################################################################################



########################################################################################################################
# Part 1 SNP identification in G. przewalskii
########################################################################################################################
## call snp 
bcftools mpileup --fasta-ref /data/frr/Gymnocarpos/genome/genome.review.assembly.FINAL.fasta  /data/frr/Gymnocarpos/mapping_result/*.bam -a FORMAT/AD,ADF,ADR,DP,SP,SCR -a INFO/AD -a INFO/ADF -a INFO/ADR -a INFO/SCR --threads 60 | bcftools call --threads 60 -mO z -o /data/frr/Gymnocarpos/mapping_result/bcf-callsnp-result/Gym178indcombine_bcftools_mpileup_call_allchr.vcf.gz
tabix /data/frr/Gymnocarpos/mapping_result/bcf-callsnp-result/Gym178indcombine_bcftools_mpileup_call_allchr.vcf.gz
## extract snp
bcftools view -v snps /data/frr/Gymnocarpos/mapping_result/bcf-callsnp-result/Gym178indcombine_bcftools_mpileup_call_allchr.vcf.gz -O z -o /data/frr/Gymnocarpos/mapping_result/all178_fltSNP/Gym_all178_fltSNP_allchr.vcf.gz
## filter by 'bcftools filter' or R
keep.snp<-snp[snp$qual>=30 & !is.na(snp$qual)& snp$DP>=500 & snp$DP<=5000 & snp$RPB>=0.05 & !is.na(snp$RPB) & snp$MQB>=0.03 & !is.na(snp$MQB) & snp$BQB>=0.1 & !is.na(snp$BQB) & snp$MQSB>=0.03 & !is.na(snp$MQSB) & snp$HOB<=0.1 & !is.na(snp$HOB) & snp$MQ>=25 & !is.na(snp$MQ) & snp$SOR<3 & !is.na(snp$SOR) & snp$mean_depth>=8 & !is.na(snp$mean_depth) & snp$QD>=2 & !is.na(snp$QD) & snp$GQ_mean>=20 & !is.na(snp$GQ_mean) ,]
## extract biallelic snp
vcftools --gzvcf /data/frr/Gymnocarpos/mapping_result/all178_fltSNP/Gym_all178_fltSNP_PASS_allchr.vcf.gz --remove-filtered-all  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c >/data/frr/Gymnocarpos/mapping_result/all178_flt2SNP_pass/Gym_all178_flt2SNP_PASS_allchr.vcf.gz
## fliter missing and minor allele frequency
vcftools --gzvcf Gym_all178_flt2SNP_PASS_allchr.vcf.gz --keep /data/frr/Gymnocarpos/mapping_result/Gym_177sample_id.txt --max-missing 0.9 --maf 0.01 --recode --stdout | bgzip -c >/data/frr/Gymnocarpos/mapping_result/all177_flt2SNP_miss09_maf001/Gym177_allchr_missing09maf001.vcf.gz

ps: Gym177_allchr_missing09maf001.vcf.gz is available for download from 'https://figshare.com/account/home#/data'

########################################################################################################################
# Part 2 Mantel and partial Mantel regression tests by R package 'vegan'
########################################################################################################################
### The correlation between genetic distance (FST/(1-Fst)) and environmental/geographic variables
library(fossil)
library(vegan)
library(dplyr)
library(reshape2)
#Geo. distance
locations <- read.table("D:/OneDrive - zju.edu.cn/Gym/environment/get_env_factor/gym_22pop_location.txt",header = TRUE) %>%tibble::column_to_rownames('pop') 
Geodistancematrix <- as.matrix( earth.dist(locations))
Fst_scalar<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/gym_22pop_fstlongitude_scalar.txt",header = TRUE) %>% as.matrix
#Env. distance
all104env<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/envfactor104_pop22_PCA_matrix.txt",header = TRUE)%>% as.matrix
longitude<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/longitude_pop22_matrix.txt",header = TRUE)%>% as.matrix
latitude<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/latitude_pop22_matrix.txt",header = TRUE)%>% as.matrix
pre<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_pre_PCA_matrix.txt",header = TRUE)%>% as.matrix
tem<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_tem_PCA_matrix.txt",header = TRUE)%>% as.matrix
srad<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_srad_PCA_matrix.txt",header = TRUE)%>% as.matrix
wind<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_wind_PCA_matrix.txt",header = TRUE)%>% as.matrix
vapr<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_vapr_PCA_matrix.txt",header = TRUE)%>% as.matrix
salt<- read.table("D:/OneDrive - zju.edu.cn/Gym/environment/get_env_factor/salt_matrix.txt",header = TRUE)%>% as.matrix
AI<- read.table("D:/OneDrive - zju.edu.cn/Gym/environment/get_env_factor/AI_matrix.txt",header = TRUE)%>% as.matrix
## Mantel test
mantel(Geodistancematrix, Fst_scalar)
mantel(Geodistancematrix, all104env) 
mantel(all104env, Fst_scalar)
mantel(longitude, Fst_scalar)
mantel(latitude, Fst_scalar)
mantel(pre, Fst_scalar)
mantel(tem, Fst_scalar)
mantel(srad, Fst_scalar)
mantel(vapr, Fst_scalar)
mantel(wind, Fst_scalar)
mantel(salt, Fst_scalar)
mantel(AI, Fst_scalar)
## partial Mantel test
mantel.partial(Geodistancematrix,Fst_scalar,all104env)
mantel.partial(all104env,Fst_scalar,Geodistancematrix)
mantel.partial(longitude, Fst_scalar,Geodistancematrix)
mantel.partial(latitude, Fst_scalar,Geodistancematrix)
mantel.partial(pre, Fst_scalar,Geodistancematrix)
mantel.partial(tem, Fst_scalar,Geodistancematrix)
mantel.partial(srad, Fst_scalar,Geodistancematrix)
mantel.partial(vapr, Fst_scalar,Geodistancematrix)
mantel.partial(wind, Fst_scalar,Geodistancematrix)
mantel.partial(salt, Fst_scalar,Geodistancematrix)
mantel.partial(AI, Fst_scalar,Geodistancematrix)

########################################################################################################################
# Part 3 The evolution of genes shared by xerophytic species
########################################################################################################################
### generalized linear model and linear mixed effect model
library(ggplot2)
library(MASS)
library(dplyr)
library(abc)
library(ape)
library(sde)
###read input file
overlap_family <- read.table("D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/pairspecies_overlap_genefamily_regroup.txt",header=T)
Xe <- overlap_family %>% filter(group1 == "Xero-Xero")
mean(Xe$overlap_genefamily)#11417.47
sd(Xe$overlap_genefamily)#575.122
noxe <-  overlap_family %>% filter(group1 == "NonX-NonX")
mean(noxe$overlap_genefamily)#10703.17
sd(noxe$overlap_genefamily)#470.0589
###Compare two generalized linear models
m1<-glm.nb(overlap_family$overlap_genefamily ~ overlap_family$branch)
m2<-glm.nb(overlap_family$overlap_genefamily  ~ overlap_family$branch + overlap_family$group1)
summary(m1)
# AIC: 1629.1
summary(m2)
# AIC: 1605.5
anova(m1,m2)
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: overlap_family$overlap_genefamily
# Model    theta Resid. df    2 x log-lik.   Test    df
# 1                         overlap_family$branch 426.9552       103       -1623.057             
# 2 overlap_family$branch + overlap_family$group1 561.5540       101       -1595.496 1 vs 2     2
# LR stat.      Pr(Chi)
# 1                      
# 2 27.56088 1.035692e-06
rsq(m1,adj=T) # variance explained by branch length, a pseduo R2
rsq(m2,adj=T) # by both branch length and selection 

#we can also test for random effect of branch length
m4=lme(sharedGenes~group1, random=~1|branch, data=overlap_family)
anova.lme(m4)
summary(m4) 

### Permutation analysis for a generalized linear model with randomly assigned xerophytic to the tips
overlap_family$group3=overlap_family$group1[sample(1:105,105,rep=F)]
m3<-glm.nb(overlap_family$overlap_genefamily  ~ overlap_family$branch + overlap_family$group3)
anova(m1,m3)
###10000 bootstrap
boot = function(x,size=105){
  index = sample(1:length(x),size,rep=F)
  return(x[index])
}
boot_type = replicate(10000,boot(overlap_family$group1))
boot_type <- data.frame(boot_type)
write.table(boot_type,file = "D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/pairspecies_overlap_genefamily_regroup_type_bootstrap.txt",quote = FALSE, row.name=F, sep = "\t")
###calculated pvalue
pvalue=NULL
for (i in 1:length(boot_type)){
m3<-glm.nb(overlap_family$overlap_genefamily  ~ overlap_family$branch + boot_type[,i])
com=anova(m1,m3)
p_value=com[2,8]
pvalue <- rbind(pvalue,data.frame(p_value))
}
write.table(pvalue,file = "D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/pairspecies_overlap_genefamily_regroup_type_bootstrap_pvalue.txt",quote = FALSE, row.name=F, sep = "\t")
summary(pvalue)

###The distribution for pvalue
png("D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/pvalue_distribution.png",width=6000,height=5500,res=600)
plot(density(pvalue$p_value),xlab="P value",main="Comparison of p values between different models",ylab="Density")
m2p=1.035692e-06
abline(v=m2p,lty=1,col='red',lwd=1.5)
dev.off()

### Ornstein-Uhlenbeck process see details in folder 'Ornstein-Uhlenbeck model'















