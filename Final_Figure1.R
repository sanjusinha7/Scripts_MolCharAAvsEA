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
tcga=read.csv('Results_New/Nov_28/Supp_Table3.csv')

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


##############
######TCGA GI
##############
##Selecting colorSet
ylim1 = boxplot.stats(tcga$Normalized_gi)$stats[c(1, 5)]
A=ggplotGrob(
  ggplot(tcga, aes(x=hist,y=Normalized_gi, fill=race))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.3, fill="white")+
  labs(title="GI in PanCan (TCGA)",x="Race", y = "Scaled GI")+
  facet_wrap(.~info_tcga.Tissue_Type+info_tcga.CellofOrigin)+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5,  method.args = list(alternative = "greater"),
                     label.y=ylim1[2], size = 10)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  coord_cartesian(ylim = ylim1*1.15)
  )
# A$layout$l[A$layout$name == "title"] <- 1

##############
####TCGA CHTP
##############
###Approprite test for barplot
fisher.test(cbind(table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA']),table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA'])),
            alternative = 'l')
###
chtp_mat=aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), data = tcga)
C=ggplotGrob(ggplot(chtp_mat, aes(y=CHTP_Canonical_Definition_Presence, x=race, fill=race))+
  geom_bar(stat = 'identity')+
  labs(title="CHTP in PanCan (TCGA)",x="Race", y = "Proportion with Chromothripsis")+
  theme_classic(base_size = 25)+
  guides(fill=FALSE)+
  annotate("text", x = 1.5, y = max(chtp_mat[,2])*1.05, label = "*", size=10)+
  coord_cartesian(ylim = c(0, max(chtp_mat[,2])*1.15))+
  scale_fill_brewer(palette=colorType_Set))
##############
####TCGA HRD
##############
ylim1 = boxplot.stats(tcga$Normalized_hrd.loh)$stats[c(1, 5)]
B=ggplotGrob(ggplot(tcga, aes(x=race,y=Normalized_hrd.loh, fill=race))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.3, fill="white")+
  labs(title="HRD in PanCan (TCGA)",x="Race", y = "Scaled HRD")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     label.y=ylim1[2], size = 10, method.args = list(alternative = "greater"))+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  coord_cartesian(ylim = ylim1*1.15))
# B$layout$l[B$layout$name == "title"] <- 1

wilcox.test(tcga$Normalized_hrd.loh~tcga$race, alternative='g')

##################
#TCGA GI in lung cancer
##################
ylim1 = boxplot.stats(tcga$Normalized_gi[lung_can])$stats[c(1, 5)]
D1=ggplotGrob(
  ggplot(tcga[lung_can,], aes(x=race,y=Normalized_gi, fill=race))+
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.3, fill="white")+
                labs(title="GI in LUSC (TCGA)",x="Race", y = "Scaled GI")+
                theme_classic(base_size = 25)+
                stat_compare_means(method = "wilcox.test", 
                                   label = "p",
                                   label.x = 1.4,
                                   label.y=ylim1[2], method.args = list(alternative = "greater"),
                                   size = 10)+
                guides(fill=FALSE)+
                scale_fill_brewer(palette=colorType_Set)+
                coord_cartesian(ylim = ylim1*1.15)
)

##################
#TCGA HRD in lung cancer
##################
ylim1 = boxplot.stats(tcga$Normalized_hrd.loh[lung_can])$stats[c(1, 5)]
E1=ggplotGrob(ggplot(tcga[lung_can,], aes(x=race,y=Normalized_hrd.loh, fill=race))+
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.3, fill="white")+
                #  facet_grid(.~hist)+
                labs(title="HRD in LUSC (TCGA)",x="Race", y = "Scaled HRD")+
                theme_classic(base_size = 25)+
                stat_compare_means(method = "wilcox.test", label = "p",
                                   label.x = 1.5, method.args = list(alternative = "greater"),
                                   label.y=ylim1[2], size = 10)+
                guides(fill=FALSE)+
                scale_fill_brewer(palette=colorType_Set) +
                coord_cartesian(ylim = ylim1*1.15)
              )
##################
#TCGA CHTP in lung cancer
##################
fisher.test(cbind(table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA' & lung_can]),
                  table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA' & lung_can])),
            alternative = 'l')
chtp_mat=aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), 
                   data = tcga[lung_can,])
F1=ggplotGrob(ggplot(chtp_mat, aes(y=CHTP_Canonical_Definition_Presence, x=race, fill=race))+
                geom_bar(stat = 'identity')+
                labs(title="CHTP in LUSC (TCGA)",x="Race", y = "Proportion with Chromothripsis")+
                theme_classic(base_size = 25)+
                guides(fill=FALSE)+
                annotate("text", x = 1.5, y = max(chtp_mat[,2])*1.05, label = "p = 0.11", size=10)+
                coord_cartesian(ylim = c(0, max(chtp_mat[,2])*1.15))+
                scale_fill_brewer(palette=colorType_Set))
##################
#NCIMDProcessing
##################
mat=read.csv('/users/sinhas8/Project_Chromotrypsis/2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
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
##################
#NCIMD- GI in lung cancer
##################
ylim1=c(0,1)

D2=ggplotGrob(
  ggplot(mat[mat$hist=='LUSC',], aes(x=race, y=Normalized_gi, fill=race))+
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.3, fill="white")+
                labs(title="GI in LUSC (NCI-MD)",x="Race", y = "Scaled GI")+
                theme_classic(base_size = 25)+
                stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                   label.y=ylim1[2], size = 10)+
                guides(fill=FALSE)+
                scale_fill_brewer(palette=colorType_Set)+
                coord_cartesian(ylim = ylim1*1.15)
)
##################
#NCIMD- GI in lung cancer
##################
E2=ggplotGrob(ggplot(mat[mat$hist=='LUSC',], aes(x=race, y=Normalized_hrd.loh, fill=race))+
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.3, fill="white")+
                labs(title="HRD in LUSC (NCI-MD)",x="Race", y = "Scaled HRD")+
                theme_classic(base_size = 25)+
                stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                   label.y=ylim1[2], size = 10)+
                guides(fill=FALSE)+
                scale_fill_brewer(palette=colorType_Set)+
                coord_cartesian(ylim = ylim1*1.15))
##################
#NCIMD- CHTP in lung cancer
##################
mat$chtp_Quan=as.numeric(mat$chtp_Quan>0)
fisher.test(cbind(table(mat$chtp_Quan[mat$hist=='LUSC' & mat$race=='EA']),
                  table(mat$chtp_Quan[mat$hist=='LUSC' & mat$race=='AA'])),
            alternative = 'g')
chtp_mat=aggregate(chtp_Quan ~ race, function(x) sum(x)/length(x), data = mat[mat$hist=='LUSC',])
F2=ggplotGrob(ggplot(chtp_mat, aes(y=chtp_Quan, x=race, fill=race))+
                geom_bar(stat = 'identity')+
                labs(title="CHTP in LUSC (NCI-MD) ",x="Race", y = "Proportion with Chromothripsis")+
                theme_classic(base_size = 25)+
                guides(fill=FALSE)+
                annotate("text", x = 1.5, y = max(chtp_mat[,2])*1.05, label = "p = 0.12", size=10)+
                coord_cartesian(ylim = c(0, max(chtp_mat[,2])*1.15))+
                scale_fill_brewer(palette=colorType_Set) 
)
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
  guides(fill=guide_legend(title="Race"))+
  scale_fill_manual(values = c('#e41a1c','#377eb8'), 
                    labels=c('African Americans', 'European Americans'),
                    name="Race")
)

##################
###Previous Final Figure
##################
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/fig1ABC_new_Apr2V2_Set4.tif', width = 1440, height = 2000)
plot_grid(leg,
          plot_grid(NULL,A,B,NULL, align='h', nrow=1, rel_widths = c(1, 2, 2, 1),
                    labels = c('',LETTERS[1:2],'',''), label_size = 35 
                    # label_x = 0, label_y = 0,  hjust = -0.5, vjust = -0.5
                    ) , 
          plot_grid(NULL,D1, E1, NULL, align='h', nrow=1, rel_widths = c(1, 2, 2, 1),
                    labels = c('',LETTERS[3:4],'',''), label_size = 35 
                    # label_x = 0, label_y = 0,  hjust = -0.5, vjust = -0.5
                    ) , 
          plot_grid(NULL,D2, E2, NULL, align='h', nrow=1, rel_widths = c(1, 2, 2, 1),
                    labels = c('',LETTERS[5:6],'',''), label_size = 35 
                    # label_x = 0, label_y = 0,  hjust = -0.5, vjust = -0.5
                    ) , 
          plot_grid(C, F1, F2,align='h', nrow=1,
                    labels = c(LETTERS[seq(7:9)]), label_size = 35
                    # ,label_x = 0, label_y = 0,  hjust = -0.5, vjust = -0.5
                    ),
          rel_heights = c(4/50, 11.5/50, 11.5/50, 11.5/50, 11.5/50), nrow=5, ncol=1
)
dev.off()


##################
#Merging all together
##################
df1=data.frame(Value=tcga$Normalized_gi, Type='Genomic Instability', Cohort='PanCan TCGA', race=tcga$race)
df2=data.frame(Value=tcga$Normalized_hrd.loh, Type='HR-deficiency', Cohort='PanCan TCGA', race=tcga$race)
df3=data.frame(Value=tcga[lung_can,]$Normalized_gi, Type='Genomic Instability', Cohort='LUSC TCGA', race=tcga[lung_can,]$race)
df4=data.frame(Value=tcga[lung_can,]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUSC TCGA', race=tcga[lung_can,]$race)
df5=data.frame(Value=mat[mat$hist=='LUSC',]$Normalized_gi, Type='Genomic Instability', Cohort='LUSC NCIMD', 
               race=mat[mat$hist=='LUSC',]$race)
df6=data.frame(Value=mat[mat$hist=='LUSC',]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUSC NCIMD',
               race=mat[mat$hist=='LUSC',]$race)

df7=data.frame(aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), data = tcga), Cohort='PanCan TCGA')
df8=data.frame(aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), data = tcga[lung_can,]), Cohort='LUSC TCGA')
df9=data.frame(aggregate(chtp_Quan ~ race, function(x) sum(x)/length(x), data = mat[mat$hist=='LUSC',]), Cohort='LUSC NCIMD')
colnames(df9)[2]='CHTP_Canonical_Definition_Presence'
##
df_CHTP=cbind(rbind(df7, df8, df9),Type='Chromothripsis')
df_HRD=rbind(df2, df4, df6)
df_GI=rbind(df1, df3, df5)
##############
######GI
##############
GI_plot=ggplotGrob(
  ggplot(df_GI, aes(x=race,y=Value, fill=race))+
               geom_violin(trim=FALSE)+
               geom_boxplot(width=0.3, fill="white")+
               labs(x="", y = "Scaled GI")+
               theme_bw(base_size = 25)+
               stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                  #label.y=ylim1[2],
                                  size = 8)+
               guides(fill=FALSE)+
               scale_fill_brewer(palette=colorType_Set)+
               facet_grid(Cohort~Type, switch = 'y')+
  theme(#strip.background.y = element_blank(),
        #strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#               coord_cartesian(ylim = ylim1*1.15))
)
##############
######HRD
##############
HRD_plot=ggplotGrob(
  ggplot(df_HRD, aes(x=race,y=Value, fill=race))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.3, fill="white")+
    labs(x="Race", y = "Scaled HR-Deficiency")+
    theme_bw(base_size = 25)+
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                       #label.y=ylim1[2],
                       size = 8)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type)+
    theme(strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
)
##############
######CHTP
##############
###Approprite test for barplot
Pvalues=data.frame(p_values=c(P1=round(fisher.test(cbind(table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA']),table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA'])),
                                                   alternative = 'l')[1]$p.value, 2),
                              P2=round(fisher.test(cbind(table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA' & lung_can]),
                                                         table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA' & lung_can])),
                                                   alternative = 'l')[1]$p.value, 2),
                              P3=round(fisher.test(cbind(table(mat$chtp_Quan[mat$hist=='LUSC' & mat$race=='EA']),
                                                         table(mat$chtp_Quan[mat$hist=='LUSC' & mat$race=='AA'])),
                                                   alternative = 'g')[1]$p.value, 2) ),
                   xstar=c(1.5), ystar=c(0.2, 0.5, 0.5), 
                   Cohort=c('PanCan TCGA', 'LUSC TCGA', 'LUSC NCIMD'),
                   Type=rep('Chromothripsis',3),
                   race=rep('AA', 3)
)
Pvalues$p_values=ifelse(Pvalues$p_values<0.05, '*','ns')
###
CHTP_plot=ggplotGrob(
  ggplot(df_CHTP, aes(y=CHTP_Canonical_Definition_Presence, x=race, fill=race))+
               geom_bar(stat = 'identity')+
               labs(x="", y = "Chromothripsis Frequency")+
               theme_bw(base_size = 25)+
               guides(fill=FALSE)+
               scale_fill_brewer(palette=colorType_Set)+
               facet_grid(Cohort~Type, scales = 'free')+
    geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=p_values), size=8)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank())
)

##############
######Put all together
##############
tiff('prep_final_figures/Delthis_Figure1_May3.tif', width = 900, height = 900)
plot_grid(leg,
          plot_grid(GI_plot, HRD_plot, CHTP_plot, align='h', nrow=1,
                    labels = 'AUTO', label_size = 35, rel_widths = c(14/45, 13/45, 10/45) ), 
          nrow=2, ncol=1, rel_heights = c(3/50, 47/50))
dev.off()



##############
######Correlation Tests
##############

cor.test(tcga$Normalized_hrd.AIL+tcga$Normalized_hrd.loh+tcga$Normalized_hrd.LST, tcga$Normalized_gi)
cor.test(tcga$Normalized_hrd.AIL[tcga$race=='EA'], tcga$Normalized_gi[tcga$race=='EA'])
cor.test(tcga$Normalized_hrd.AIL[tcga$race=='AA'], tcga$Normalized_gi[tcga$race=='AA'])

cor.test(mat$Normalized_hrd.loh[mat$race=='AA' & mat$hist=='LUSC' ], mat$Normalized_gi[mat$race=='AA' & mat$hist=='LUSC'])
cor.test(mat$Normalized_hrd.loh[mat$race=='EA' & mat$hist=='LUSC' ], mat$Normalized_gi[mat$race=='EA' & mat$hist=='LUSC'])

cancer_type_specific_corr<-function(cancer_type='LUSC'){
  A=cor.test(tcga$Normalized_hrd.AIL[tcga$race=='EA' & tcga$hist==cancer_type], 
           tcga$Normalized_gi[tcga$race=='EA' & tcga$hist==cancer_type])
  B=cor.test(tcga$Normalized_hrd.AIL[tcga$race=='AA' & tcga$hist==cancer_type],
           tcga$Normalized_gi[tcga$race=='AA' & tcga$hist==cancer_type])
  C=cor.test(tcga$Normalized_hrd.AIL[tcga$hist==cancer_type],
             tcga$Normalized_gi[tcga$hist==cancer_type])
  cbind(data.frame(All_Samples_rho=C$estimate, All_Samples_Sig=C$p.value),
  data.frame(AA_Samples_rho=B$estimate, AA_Samples_Sig=B$p.value),
  data.frame(EA_Samples_rho=A$estimate, EA_Samples_Sig=A$p.value))
}
df2write=t(sapply(levels(tcga$hist), function(x) cancer_type_specific_corr(cancer_type=x)))
write.csv(df2write, '/Users/sinhas8/tableS19_AIL.csv')

df2write1=summary(lm(Normalized_gi~race+stage+
                       age+GENDER+smoke+as.numeric(as.character(unlist(packyrs))), 
                     data = tcga$X ))$coefficients
