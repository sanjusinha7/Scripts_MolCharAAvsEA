####################################################################
# Get Files needed
####################################################################
setwd('/Users/sinhas8/Project_Chromotrypsis/2.Data/')
DDR=read.csv('/Users/sinhas8/Downloads/DDR_Genes.csv')
df=read.csv('PCA_pathVar_integrated_filtered_adjusted.tsv', sep='\t')
demo=read.csv('clinical_PANCAN_patient_with_followup.tsv', sep='\t')
####################################################################
# Preprocessing
####################################################################
colnames(DDR$Gene.Symbol); colnames(DDR)
DDR_Geneset_List=apply(DDR[,c(11:29)], 2, function(x) unlist(x[x!='']) )
HR_Genes=unlist(DDR[,15][DDR[,15]!=''])
core_HRGenes=unlist(DDR[,25][DDR[,25]!=''])
levels(demo[,c('race')])[c(6,8)]=c('AA','EA')
demo=demo[demo[,c('race')]=='AA' | demo[,c('race')]=='EA',]
df=df[!is.na(match(df$bcr_patient_barcode, demo$bcr_patient_barcode)),]
levels(df$race)[c(7, 9)]=c('AA', 'EA')
dim(df); dim(demo)
demo[1:5, 1:5]
####################################################################
# Functions defined
####################################################################
# Pan cancer analysis
PanCan_Enrichmentof_Germline_Deficiency_byRace<-function(Geneset_OI=HR_Genes){
  Geneset_OI=as.character(unlist(Geneset_OI))
  totalAA=table(demo$race)['AA']
  totalEA=table(demo$race)['EA']
  GermlineHRD_Total=length(na.omit(match(df$HUGO_Symbol,
                                         Geneset_OI)))
  GermlineHRD_inAA=length(na.omit(match(df$HUGO_Symbol[df$race=='AA'],
                                        Geneset_OI)))
  GermlineHRD_inEA=length(na.omit(match(df$HUGO_Symbol[df$race=='EA'],
                                        Geneset_OI)))
  pvalue=phyper(GermlineHRD_inAA, 
                totalAA, totalEA,
                GermlineHRD_Total, 
                lower.tail = F)
  if(pvalue==0){
    return(NA)
  } else{
    return(data.frame(pvalue, GermlineHRD_inAA/totalAA, GermlineHRD_inEA/totalEA ))
    }
  
}
Enrichmentof_Germline_Deficiency_byRace<-function(Geneset_OI=as.character(unlist(il6$Genes)),
                                                  cancer_type='LUAD'){
  Geneset_OI=as.character(unlist(Geneset_OI))
  totalAA=table(demo$race[demo$acronym==cancer_type])['AA']
  totalEA=table(demo$race[demo$acronym==cancer_type])['EA']
  GermlineHRD_Total=length(na.omit(match(df$HUGO_Symbol[df$cancer==cancer_type],
                                         Geneset_OI)))
  GermlineHRD_inAA=length(na.omit(match(df$HUGO_Symbol[df$race=='AA' &
                                                         df$cancer==cancer_type],
                                        Geneset_OI)))
  GermlineHRD_inEA=length(na.omit(match(df$HUGO_Symbol[df$race=='EA' & 
                                                         df$cancer==cancer_type],
                                        Geneset_OI)))
  
  pvalue=phyper(GermlineHRD_inAA, totalAA, totalEA, GermlineHRD_Total, lower.tail = F)
  if(pvalue==0 | GermlineHRD_inAA==0){
    return(data.frame(NA, NA, NA))
  } else{ return(data.frame(pvalue, GermlineHRD_inAA/totalAA, GermlineHRD_inEA/totalEA ))}
}
####################################################################
# calling functions
####################################################################
# Cancer type specific analysis
RaceCount_across_canType=aggregate(demo$race~ demo$acronym, 
                                   demo, 
                                   function(x) c(sum(x=='AA'),
                                                 sum(x=='EA')) )
sum(RaceCount_across_canType[,2][,2])
CanSpec_df2write=lapply(DDR_Geneset_List, 
                        function(x) sapply(levels(demo$acronym), function(y) 
                          Enrichmentof_Germline_Deficiency_byRace(Geneset_OI = x,
                                                                  cancer_type = y)))

# For core HRD
CanSpec_df2write$Homology.dependent.recombination..HDR.=CanSpec_df2write$Homology.dependent.recombination..HDR.[,unlist(!is.na(CanSpec_df2write$Homology.dependent.recombination..HDR.)[1,])]
CompHRD_pvalues=unlist(CanSpec_df2write$Homology.dependent.recombination..HDR.[1,])
HRD=CanSpec_df2write$Homology.dependent.recombination..HDR.[2:3,]
HRD=HRD[,apply(HRD, 2, function(x) sum(unlist(x))>0 ) ]
df_HRD=data.frame(Prop = unlist(c(HRD)),
                  Race=rep(c('AA', 'EA'), ncol(HRD)),
                  cancer_type=rep(colnames(HRD), each=2)  )
annotation=data.frame(pvalues=CompHRD_pvalues[colnames(HRD)],
                xstar=1.5,
                ystar=0.20,
                cancer_type=levels(df_HRD$cancer_type),
                Race='AA')
####################################################################
# Extended Figure 10
####################################################################
colorType_Set='Set1'
Plot1<-ggplot(df_HRD, aes(y=Prop, x=Race, fill=Race))+
  geom_bar(stat = 'identity')+
  labs(x="Race", y = "Germline HR-deficiency Proportion")+
  theme_bw(base_size = 25)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  facet_wrap(cancer_type~., ncol=4)+
  ylim(... = c(0,0.25))+
  #  geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=p_values), size=8)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank())
    
Plot1=Plot1+geom_text(data=annotation,
                  aes(x=xstar, y=ystar,
                      label = paste0('p= ', format.pval(pvalues, digits=1))), size=6)+
    geom_text(aes(label=round(as.numeric(Prop), 2)),
              position=position_dodge(width=0.9), vjust=-0.25, size=6)
  
  
setEPS()
postscript('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Ext_Figure10.eps',
     width = 12, height = 12)
Plot1
dev.off()
####################################################################
# For LUSC
####################################################################
df_HRD=data.frame(Prop=unlist(c(HRD)),
                  Race=rep(c('AA', 'EA'), ncol(HRD)),
                  cancer_type=rep(colnames(HRD), each=2),
                  Type= 'Germline HRD')
df_HRD=df_HRD[df_HRD$cancer_type=='LUSC',]
pvalues=CompHRD_pvalues[as.character(df_HRD$cancer_type[1]) ]
df_HRD$cancer_type='LUSC TCGA'
df_HRD$cancer_type=as.character(df_HRD$cancer_type)
annotation=data.frame(pvalues,
                      xstar=1.5,
                      ystar=max(df_HRD$Prop)+0.01,
                      cancer_type=levels(factor(df_HRD$cancer_type)),
                      Race='AA')
####################################################################
##Pancancer
####################################################################
df_HRD_PanCancer=data.frame(Prop=c(unlist(PanCan_Enrichmentof_Germline_Deficiency_byRace(HR_Genes)[-1])),
                           Race=c('AA','EA'), cancer_type=c('PanCancer TCGA'),
                  Type= 'Germline HRD')
df=rbind(df_HRD_PanCancer, df_HRD)
annotation_panCancer=data.frame(pvalues=0.02,
                      xstar=1.5,
                      ystar=max(df_HRD_PanCancer$Prop)+0.01,
                      cancer_type=levels(factor(df_HRD_PanCancer$cancer_type)),
                      Race='AA')
ann= rbind(annotation, annotation_panCancer)
colorType_Set='Set1'
table(df$Race[df$cancer_type==levels(df$cancer_type)[1]])
germHRD_plot=ggplot(df, aes(y=Prop, x=Race, fill=Race))+
  geom_bar(stat = 'identity')+
  labs(x="", y = "Germline HR-Deficiency")+
  theme_bw(base_size = 25)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  facet_wrap(.~cancer_type, scales = 'free')+
  #  geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=p_values), size=8)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )+
  geom_text(data=ann,
            aes(x=xstar, y=ystar,label = format.pval(pvalues, digits=1)), size=9)
##################
# Getting Legend
##################
# tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/germLing_PanandLUSC.tif',
#      width = 600, height = 400)
# plot_grid(leg, germHRD_plot, ncol=1, rel_heights = c(c(4/50, 46/50)))
# dev.off()
setEPS()
postscript('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure6_vecForm.eps',
           width = 8, height = 5.3)
plot_grid(leg, germHRD_plot, ncol=1, rel_heights = c(c(4/50, 46/50)))
dev.off()
##################
# For core HRD
##################
CompHRD_pvalues=unlist(CanSpec_df2write$Homology.dependent.recombination..HDR.[1,])
HRD=CanSpec_df2write$Homology.dependent.recombination..HDR.[2:3,]
HRD=HRD[,apply(HRD, 2, function(x) sum(unlist(x))>0 ) ]
df_HRD=data.frame(Prop = unlist(c(HRD)),
                  Race=rep(c('AA', 'EA'), ncol(HRD)),
                  cancer_type=rep(colnames(HRD), each=2)  )
df_HRD=na.omit(df_HRD)
annotation=data.frame(pvalues=CompHRD_pvalues[levels(df_HRD$cancer_type)],
                      xstar=1.5,
                      ystar=sapply(split(df_HRD$Prop,
                                         df_HRD$cancer_type), function(x) max(x)+0.01 ),
                      cancer_type=levels(df_HRD$cancer_type),
                      Race='AA')
####################################################################
# For all the cancer types
####################################################################
colorType_Set='Set1'
Plot1<-ggplot(df_HRD, aes(y=Prop, x=Race, fill=Race))+
  geom_bar(stat = 'identity')+
  labs(x="Race", y = "Germline HR-deficiency Proportion")+
  theme_bw(base_size = 25)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  facet_wrap(cancer_type~., ncol=4, scales = 'free')
#  geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=p_values), size=8)+
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
)

Plot1=Plot1+geom_text(data=annotation,
                      aes(x=xstar, y=ystar,
                          label = format.pval(pvalues, digits=1)), size=9)

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/germline_Figurev1_May8th.tif',
     width = 900, height = 900)
Plot1
dev.off()

####################################################################
##For LUSC
####################################################################
df_HRD=data.frame(Prop=unlist(c(HRD)),
                  Race=rep(c('AA', 'EA'), ncol(HRD)),
                  cancer_type=rep(colnames(HRD), each=2),
                  Type= 'Germline HRD')
df_HRD=df_HRD[df_HRD$cancer_type=='LUSC',]
pvalues=CompHRD_pvalues[as.character(df_HRD$cancer_type[1]) ]
df_HRD$cancer_type='LUSC TCGA'
df_HRD$cancer_type=as.character(df_HRD$cancer_type)
annotation=data.frame(pvalues,
                      xstar=1.5,
                      ystar=max(df_HRD$Prop)+0.01,
                      cancer_type=levels(factor(df_HRD$cancer_type)),
                      Race='AA')
####################################################################
##Pancancer
####################################################################
df_HRD_PanCancer=data.frame(Prop=c(unlist(PanCan_Enrichmentof_Germline_Deficiency_byRace(HR_Genes)[-1])),
                            Race=c('AA','EA'), cancer_type=c('PanCancer TCGA'),
                            Type= 'Germline HRD')
df=rbind(df_HRD_PanCancer, df_HRD)
annotation_panCancer=data.frame(pvalues=0.02,
                                xstar=1.5,
                                ystar=max(df_HRD_PanCancer$Prop)+0.01,
                                cancer_type=levels(factor(df_HRD_PanCancer$cancer_type)),
                                Race='AA')
ann= rbind(annotation, annotation_panCancer)
########
colorType_Set='Set1'
germHRD_plot=ggplot(df, aes(y=Prop, x=Race, fill=Race))+
  geom_bar(stat = 'identity')+
  labs(x="", y = "Germline HR-Deficiency")+
  theme_bw(base_size = 25)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  facet_wrap(.~cancer_type, scales = 'free')+
  #  geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=p_values), size=8)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )+
  geom_text(data=ann,
            aes(x=xstar, y=ystar,label = format.pval(pvalues, digits=1)), size=9)
##################
###Getting Legend
##################
ylim1=c(0,1)
leg=get_legend(ggplot(tcga, aes(x=race,y=Normalized_gi, fill=race))+
                 geom_violin(trim=FALSE)+
                 geom_boxplot(width=0.3, fill="white")+
                 labs(title="GI in PanCan (TCGA)",x="Race", y = "Scaled GI")+
                 theme_classic(base_size = 25)+
                 stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                    label.y=ylim1[2], size = 10)+
                 theme(legend.box = "horizontal", legend.position = "bottom",
                       legend.text=element_text(size=20), legend.title=element_text(size=25),
                       legend.key.size = unit(2,"line"))+
                 guides(fill=guide_legend(title="Event Type"))+
                 scale_fill_manual(values = c("blue","limegreen"), 
                                   labels=c('Del', 'Amp'),
                                   name="Event Type")
)

##################
###Getting Legend
##################
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/germLing_PanandLUSC.tif',
     width = 600, height = 400)
plot_grid(leg, germHRD_plot, ncol=1, rel_heights = c(c(4/50, 46/50)))
dev.off()


tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Events_Legend.tif',
     width = 300, height = 50)
plot(leg)
dev.off()
