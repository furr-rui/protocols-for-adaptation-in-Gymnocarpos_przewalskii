########################################################################################################################
# Documents and scripts were written by: Ruirui Fu 
# For manuscript: Ruirui Fu, et al. (2024). "Shared xerophytic genes and their re-use in local adaptation to aridity in a desert plant Gymnocarpos przewalskii" (under review). 
# email: furuirui@zju.edu.cn
# Jun Chen Lab. 
########################################################################################################################
library(ggplot2)
library(MASS)
library(dplyr)
library(rsq)
library(abc)
library(ape)
library(sde)
####################################Eleven xerophytes and four non-xerophytes####################################
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



####################################Four non-xerophytes and four xerophytes were randomly assigned####################################
###read input file
overlap_family <- read.table("D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/pairspecies_overlap_genefamily_regroup.txt",header=T)
xerophyte <- c("A.nanus","C.arietinum","G.przewalskii","L.ruth","P.vera","Pco" ,"Peu","Phacu","Pitaya","Rcu","jojoba")%>% as.data.frame()
non_xerophyte <- c("B.vulgaris", "D.caryophyllus","F.tataricum","Spov")%>% as.data.frame()
###Four xerophytes were randomly extracted
keep_xerophyte <- xerophyte$.[sample(1:11,4,rep=F)] %>% as.data.frame()
# > keep_xerophyte 
# [1] "G.przewalskii" "Phacu"         "P.vera"        "jojoba"
###get information for the number of shared genefamily and branch length 
overlap_family_fourxe_fournox1 <- NULL
overlap_family_fourxe_fournox3 <- NULL
 for(i in 1:3) {
   for (j in (i+1):4) {
     overlap_family_fourxe_fournox_tem1 <- overlap_family %>% filter((Species1==keep_xerophyte[i,] & Species2==keep_xerophyte[j,])|(Species2==keep_xerophyte[i,] & Species1==keep_xerophyte[j,]))
     overlap_family_fourxe_fournox1=rbind(overlap_family_fourxe_fournox1,overlap_family_fourxe_fournox_tem1)
     overlap_family_fourxe_fournox_tem3 <- overlap_family %>% filter((Species1==non_xerophyte[i,] & Species2==non_xerophyte[j,])|(Species2==non_xerophyte[i,] & Species1==non_xerophyte[j,]))
     overlap_family_fourxe_fournox3=rbind(overlap_family_fourxe_fournox3,overlap_family_fourxe_fournox_tem3)
     
   }
 }
overlap_family_fourxe_fournox2 <- NULL
for(i in 1:4) {
  for (j in 1:4) {
    overlap_family_fourxe_fournox_tem2 <- overlap_family %>% filter((Species1==keep_xerophyte[i,] & Species2==non_xerophyte[j,])|Species2==keep_xerophyte[i,] & Species1==non_xerophyte[j,])
    overlap_family_fourxe_fournox2=rbind(overlap_family_fourxe_fournox2,overlap_family_fourxe_fournox_tem2)
  }
}
overlap_family_fourxe_fournox <- rbind(overlap_family_fourxe_fournox1,overlap_family_fourxe_fournox2,overlap_family_fourxe_fournox3)
write.table(overlap_family_fourxe_fournox,file = "D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/noxe_xe_same4species/pairspecies_overlap_genefamily_regroup_type_same4species_extract.txt",quote = FALSE, row.name=F, sep = "\t")

###Compare two generalized linear models
m1<-glm.nb(overlap_family_fourxe_fournox$overlap_genefamily ~ overlap_family_fourxe_fournox$branch)
m2<-glm.nb(overlap_family_fourxe_fournox$overlap_genefamily  ~ overlap_family_fourxe_fournox$branch + overlap_family_fourxe_fournox$group1)
summary(m1)
# AIC: 436.32
summary(m2)
# AIC: 427.86
anova(m1,m2)
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: overlap_family_fourxe_fournox$overlap_genefamily
# Model    theta Resid. df    2 x log-lik.   Test
# 1                                        overlap_family_fourxe_fournox$branch 463.3162        26       -430.3224       
# 2 overlap_family_fourxe_fournox$branch + overlap_family_fourxe_fournox$group1 740.3854        24       -417.8596 1 vs 2
# df LR stat.     Pr(Chi)
# 1                           
# 2     2  12.4628 0.001966699
rsq(m1,adj=T) # variance explained by branch length, a pseduo R2
# 0.0102301
rsq(m2,adj=T) # by both branch length and selection 
# 0.3133765

#we can also test for random effect of branch length
m4=lme(sharedGenes~group1, random=~1|branch, data=overlap_family)
anova.lme(m4)
summary(m4) 

###Permutation analysis for a generalized linear model with randomly assigned xerophytic to the tips
overlap_family_fourxe_fournox$group3=overlap_family_fourxe_fournox$group1[sample(1:28,28,rep=F)]
m3<-glm.nb(overlap_family_fourxe_fournox$overlap_genefamily  ~ overlap_family_fourxe_fournox$branch + overlap_family_fourxe_fournox$group3)
anova(m1,m3)

####10000 bootstrap
boot = function(x,size=28){
  index = sample(1:length(x),size,rep=F)
  return(x[index])
}
boot_type = replicate(10000,boot(overlap_family_fourxe_fournox$group1))
boot_type <- data.frame(boot_type)
write.table(boot_type,file = "D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/noxe_xe_same4species/pairspecies_overlap_genefamily_regroup_type_same4species_bootstrap.txt",quote = FALSE, row.name=F, sep = "\t")
###calculated pvalue
pvalue=NULL
for (i in 1:length(boot_type)){
  m3<-glm.nb(overlap_family_fourxe_fournox$overlap_genefamily  ~ overlap_family_fourxe_fournox$branch + boot_type[,i])
  com=anova(m1,m3)
  p_value=com[2,8]
  pvalue <- rbind(pvalue,data.frame(p_value))
}
write.table(pvalue,file = "D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/noxe_xe_same4species/pairspecies_overlap_genefamily_regroup_type_same4species_bootstrap_pvalue.txt",quote = FALSE, row.name=F, sep = "\t")

###The distribution for pvalue
png("D:/OneDrive - zju.edu.cn/Gym/resequence/environment/Osdesert16sp/sp16desert_genefamily/species_genefamily_raw/anytwospecies_overlap_genefamily/boot/noxe_xe_same4species/pvalue_same4species_distribution.png",width=6000,height=5500,res=600)
plot(density(pvalue$p_value),xlab="P value",main="Comparison of p values between different models",ylab="Density")
m2p=0.001966699
abline(v=m2p,lty=1,col='red',lwd=1.5)
dev.off()

