##################
###Correcting for confounding factors
##################
##################
###Def Functions
##################
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
##################
###Read files
##################
setwd('/users/sinhas8/Project_Chromotrypsis/')
mat=read.csv('2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
age_sex=read.csv('2.Data/age_gender_info.csv')
ncimd=cbind(ncimd, age_sex[match(ncimd$acc,age_sex$accession),])
##################
###Process files
##################
mat=mat[mat$hist=='adeno' | mat$hist =='sq',]
mat$hist=as.character(mat$hist)
mat$hist=factor(mat$hist)
levels(mat$hist)=c('LUAD', 'LUSC')
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
mat=mat[order(mat$hist),]
mat$Normalized_gi=scaling_cancerType(mat$CNV_burden, mat$hist)
mat$Normalized_hrd.loh=scaling_cancerType(mat$HRD_by_LOH, mat$hist)
mat$Normalized_Chromothripsis_Presence=scaling_cancerType(mat$chtp_Quan, mat$hist)
prev_mat=mat
mat=cbind(mat, age_sex[match(mat$acc, age_sex$accession),])

#################
#Multivariate regession
#################
####For NCIMD
summary(lm(Normalized_gi~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs))), 
           data = mat[mat$hist=='LUAD',]))$coefficients
summary(lm(Normalized_gi~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs))), 
           data = mat[mat$hist=='LUAD',]))$coefficients

summary(lm(Normalized_hrd.loh~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs))), 
           data = mat[mat$hist=='LUAD',]))$coefficients

summary(lm(Normalized_gi~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs)))+purity, 
           data = mat[mat$hist=='LUSC',]))$coefficients
summary(lm(Normalized_hrd.loh~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs)))+purity,  
           data = mat[mat$hist=='LUSC',]))$coefficients
summary(lm(Normalized_gi~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs)))+purity,  
           data = mat[mat$hist=='LUSC',]))$coefficients

summary(lm(chtp_Quan~race+stage+
             age+GENDER+smoke+as.numeric(as.character(unlist(packyrs))), 
           data = mat[mat$hist=='LUSC',]))$coefficients



###For TCGA
load('/Users/sinhas8/ISLE-Basic/data/TCGA.RData')
demo_tcga=data.frame(sex=prob$sex[match(tcga$info_tcga.Sample_Name, prob$samples)],
           age=prob$age[match(tcga$info_tcga.Sample_Name, prob$samples)],
           #race=prob$race[match(tcga$info_tcga.Sample_Name, prob$samples)],
           stage=prob$stage[match(tcga$info_tcga.Sample_Name, prob$samples)]
           #types=prob$types[match(tcga$info_tcga.Sample_Name, prob$samples)]
           )
head(demo_tcga)

tcga=cbind(tcga, demo_tcga)
df2write1_tcga=summary(lm(Normalized_gi~race+stage+age+sex, 
                     data = tcga))$coefficients
df2write1_tcga=summary(lm(Normalized_hrd.loh ~ race+stage+age+sex, 
                          data = tcga))$coefficients
df2write1_tcga=summary(lm(Normalized_hrd.LST ~ race+stage+age+sex, 
                          data = tcga))$coefficients
summary(lm(Normalized_hrd.AIL ~ race+stage+age+sex, 
           data = tcga))$coefficients
summary(lm(Normalized_Chromothripsis_Presence ~ race+stage+age+sex, 
           data = tcga))$coefficients

wilcox.test(tcga$HRD_by_LST ~tcga$CHTP_Canonical_Definition_Presence, alternative='l')$p.value
wilcox.test(tcga$HRD_by_LOH ~tcga$CHTP_Canonical_Definition_Presence, alternative='l')$p.value
wilcox.test(tcga$HRD_by_AIL ~tcga$CHTP_Canonical_Definition_Presence, alternative='l')$p.value

##For AA
wilcox.test(tcga$HRD_by_LST[tcga$race=='AA'] ~ tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA'], alternative='l')$p.value
wilcox.test(tcga$HRD_by_LOH[tcga$race=='AA'] ~ tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA'], alternative='l')$p.value
wilcox.test(tcga$HRD_by_AIL[tcga$race=='AA'] ~ tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA'], alternative='l')$p.value

wilcox.test(tcga$HRD_by_LST[tcga$race=='EA'] ~tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA'], alternative='l')$p.value
wilcox.test(tcga$HRD_by_LOH[tcga$race=='EA'] ~tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA'], alternative='l')$p.value
wilcox.test(tcga$HRD_by_AIL[tcga$race=='EA'] ~tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA'], alternative='l')$p.value
