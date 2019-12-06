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
##################
#Merging all together
##################
df1=data.frame(Value=tcga$Normalized_gi, Type='Genomic Instability', Cohort='PanCan TCGA', race=tcga$race)
df2=data.frame(Value=tcga$Normalized_hrd.loh, Type='HR-deficiency', Cohort='PanCan TCGA', race=tcga$race)
df3=data.frame(Value=tcga[lung_can,]$Normalized_gi, Type='Genomic Instability', Cohort='LUSC TCGA', race=tcga[lung_can,]$race)
df4=data.frame(Value=tcga[lung_can,]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUSC TCGA', race=tcga[lung_can,]$race)
# df5=data.frame(Value=mat[mat$hist=='LUSC',]$Normalized_gi, Type='Genomic Instability', Cohort='LUSC NCIMD', 
#                race=mat[mat$hist=='LUSC',]$race)
# df6=data.frame(Value=mat[mat$hist=='LUSC',]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUSC NCIMD',
#                race=mat[mat$hist=='LUSC',]$race)
# 
df7=data.frame(aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), data = tcga), Cohort='PanCan TCGA')
df8=data.frame(aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), data = tcga[lung_can,]), Cohort='LUSC TCGA')
# df9=data.frame(aggregate(chtp_Quan ~ race, function(x) sum(x)/length(x), data = mat[mat$hist=='LUSC',]), Cohort='LUSC NCIMD')
# colnames(df9)[2]='CHTP_Canonical_Definition_Presence'
##
df_CHTP=cbind(rbind(df7),Type='Chromothripsis')
df_HRD=rbind(df2)
df_GI=rbind(df1)
##############
######GI
##############
GI_plot=ggplotGrob(
  ggplot(df_GI, aes(x=race,y=Value, fill=race))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.3, fill="white")+
    labs(x="", y = "Scaled GI")+
    theme_bw(base_size = 25)+
    stat_compare_means(method = "wilcox.test", label = "p",
                       label.x = 1,
                       label.y=1.3,
                       size = 8)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type, switch = 'y')+
    theme(axis.title.y = element_text(size=18),
#strip.background.y = element_blank(),
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
    stat_compare_means(method = "wilcox.test", label = "p", 
                       label.x = 1.2,
                       label.y=1.2,
                       size = 8)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type)+
    theme(axis.title.y = element_text(size=18),
          strip.background.y = element_blank(),
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
                                                   alternative = 'l')[1]$p.value, 2)),
                   xstar=c(1.5), 
                   ystar=c(0.25), 
                   Cohort=c('PanCan TCGA'),
                   Type=rep('Chromothripsis'),
                   race=rep('AA')
)
#Pvalues$p_values=ifelse(Pvalues$p_values<0.05, '*','ns')
###

df_counts = cbind(rbind(data.frame(aggregate(X ~ race, 
                                       function(x) length(x), data = tcga),
                                   Cohort='PanCan TCGA')),
                  df_CHTP[,-c(1,3)], x_coor=c(1, 2))
CHTP_plot=ggplotGrob(
  ggplot(df_CHTP, aes(y=CHTP_Canonical_Definition_Presence, x=race, fill=race))+
    geom_bar(stat = 'identity')+
    labs(x="", y = "Chromothripsis Frequency")+
    theme_bw(base_size = 25)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type, scales = 'free')+
    geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=paste0('p=', p_values, sep='')), size=8)+
    geom_text(data=df_counts, aes(label=paste('n=',X, sep=''), x=x_coor), vjust=-0.3, size=5)+
    theme(axis.title.y = element_text(size=18),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank())
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
##############
######Put all together
##############
tiff('prep_final_figures/Delthis_Figure5_Nov3.tif', width = 900, height = 300)
plot_grid(leg,
          plot_grid(GI_plot, HRD_plot, CHTP_plot, align='h', nrow=1,
                    labels = 'AUTO', label_size = 35, rel_widths = c(14/45, 13/45, 10/45),
                    label_x =  c(0,-0.08,-0.08) ), 
          nrow=2, ncol=1, rel_heights = c(5/50, 47/50))
dev.off()
##############
######Create a vector version
##############
setEPS()
postscript('prep_final_figures/Fig5_vectorForm.eps', height=4, width = 12)
plot_grid(leg,
          plot_grid(GI_plot, HRD_plot, CHTP_plot, align='h', nrow=1,
                    labels = 'AUTO', label_size = 35, rel_widths = c(14/45, 13/45, 10/45),
                    label_x =  c(0,-0.08,-0.08) ), 
          nrow=2, ncol=1, rel_heights = c(5/50, 47/50))
dev.off()
