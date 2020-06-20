# WES TCGA matrix
##################
#TCGA:in-Pan cancer Procescessing
##################
require(ggplot2)
require('ggpubr')
require(grid)
require(gridExtra)
require(cowplot)

colorType_Set='Set1'
setwd('/Users/sinhas8/Project_Chromotrypsis/')
tcga=read.csv('/Users/sinhas8/Project_Chromotrypsis/Results_New/Nov_28/Supp_Table3.csv')

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
tcga$Normalized_gi=scaling_cancerType(tcga$CNV_Burden, tcga$info_tcga.hist)
tcga$Normalized_hrd.loh=scaling_cancerType(tcga$HRD_by_LOH, tcga$info_tcga.hist)
tcga$Normalized_hrd.LST=scaling_cancerType(tcga$HRD_by_LST, tcga$info_tcga.hist)
tcga$Normalized_hrd.AIL=scaling_cancerType(tcga$HRD_by_AIL, tcga$info_tcga.hist)
tcga$Normalized_Chromothripsis_Presence=scaling_cancerType(tcga$CHTP_Canonical_Definition_Presence, tcga$info_tcga.hist)
colnames(tcga)[6]='race'
levels(tcga$race)[2]=c('EA')
colnames(tcga)[7]='hist'
lung_can=tcga$hist=='LUSC'

##################
#TCGA:Mutation :: GEne counts
##################
prob_wdMut=readRDS('/Users/sinhas8/Downloads/TCGA_withMut.RDS')
MutLoad_NOG=apply(prob_wdMut$Mut, 2, sum)
tcga$mutLoad=MutLoad_NOG[match(tcga$info_tcga.Sample_Name, names(MutLoad_NOG))]
tcga=tcga[!is.na(tcga$mutLoad),]
tcga$scaled_mutLoad=scaling_cancerType(tcga$mutLoad, tcga$hist)
tcga$scaled_mutLoad=range01(tcga$mutLoad)

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Mut_Load_Feb7.tif', width = 1800, height = 600)
ggplot(tcga, aes(x=as.character(hist),y=scaled_mutLoad, fill=race))+
  geom_boxplot(data=tcga, aes(fill=race))+
  labs(title="MutLoad in 23 cancer types (TCGA)",x="Race", y = "Scaled MutLoad")+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     label.y=rep(c(0.8,1.0),  length(levels(tcga$hist))/2, 0.8),
                     size = 6)+
  # guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()
length(prob_wdMut$types)

Mut_LUSC=prob_wdMut$Mut[,prob_wdMut$types=='LUSC']
write.csv(Mut_LUSC, '/Users/sinhas8/data_For_Brid/LUSC/Mut_LUSC.csv')
metadata=data.frame(sample_ID=colnames(prob_wdMut$Mut)[prob_wdMut$types=='LUSC'], 
                    race=prob_wdMut$race[prob_wdMut$types=='LUSC'],
                    hist=prob_wdMut$types[prob_wdMut$types=='LUSC'],
                    stage=prob_wdMut$stage[prob_wdMut$types=='LUSC'],
                    age=prob_wdMut$age[prob_wdMut$types=='LUSC'],
                    sex=prob_wdMut$sex[prob_wdMut$types=='LUSC'])
write.csv(metadata, '/Users/sinhas8/data_For_Brid/LUSC/metadata_LUSC.csv')

##################
#TCGA:Mutation :: GEne counts
##################
