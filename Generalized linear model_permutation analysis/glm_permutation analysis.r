########################################################################################################################
# Documents and scripts were written by: Ruirui Fu 
# For manuscript: Ruirui Fu, et al. (2024). "Shared xerophytic genes and their re-use in local adaptation to aridity in a desert plant Gymnocarpos przewalskii" (under review). 
# email: furuirui@zju.edu.cn
# Jun Chen Lab. 
########################################################################################################################
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

###Permutation analysis for a generalized linear model with randomly assigned xerophytic to the tips
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
